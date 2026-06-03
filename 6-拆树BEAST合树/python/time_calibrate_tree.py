#!/usr/bin/env python3
"""Scale an ultrametric tree from substitution units into years."""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

from phylo_merge_common import validate_branch_lengths_complete
from phylo_split_common import (
    PipelineError,
    assign_node_ids,
    compute_tip_hash,
    get_tip_names_from_clade,
    load_table,
    read_newick_tree,
    setup_logger,
    write_table,
    write_tree,
)


NODE_TIME_TABLE_COLUMNS = [
    "node_id",
    "node_type",
    "n_descendant_tips",
    "height_years",
    "depth_from_root_years",
    "hpd_low_years",
    "hpd_high_years",
    "source_height_subst_per_site",
    "hpd_source",
    "branch_length_years",
    "representative_tips",
]


TIME_CALIBRATION_SUMMARY_COLUMNS = [
    "input_tree",
    "output_tree",
    "method",
    "branch_length_unit",
    "substitution_rate_per_site_per_year",
    "sequence_length",
    "node_calibration_tip_name",
    "node_calibration_divergence_years",
    "node_calibration_input_root_to_tip",
    "scale_factor_years_per_input_branch_unit",
    "nonroot_edge_count",
    "input_branch_min",
    "input_branch_median",
    "input_branch_max",
    "output_branch_min_years",
    "output_branch_median_years",
    "output_branch_max_years",
]

TIME_CALIBRATION_EDGE_COLUMNS = [
    "node_id",
    "edge_role",
    "tip_name",
    "input_branch_length",
    "calibrated_branch_length_years",
]


def resolve_scale_factor(
    tree,
    method: str,
    substitution_rate_per_site_per_year: float,
    sequence_length: int,
    branch_length_unit: str,
    node_calibration_tip_name: str,
    node_calibration_divergence_years: float,
) -> tuple[float, str]:
    if method == "molecular_clock":
        if substitution_rate_per_site_per_year <= 0:
            raise PipelineError("substitution_rate_per_site_per_year must be positive.")
        if sequence_length <= 0:
            raise PipelineError("sequence_length must be a positive integer.")
        if branch_length_unit == "substitutions_per_site":
            return 1.0 / substitution_rate_per_site_per_year, ""
        if branch_length_unit == "substitutions_per_sequence":
            return 1.0 / (substitution_rate_per_site_per_year * float(sequence_length)), ""
        raise PipelineError(
            "Unsupported branch_length_unit. Expected substitutions_per_site or substitutions_per_sequence."
        )
    if method == "node_calibration":
        if node_calibration_divergence_years <= 0:
            raise PipelineError("node_calibration_divergence_years must be positive.")
        tip_lookup = {str(tip.name): tip for tip in tree.get_terminals()}
        if node_calibration_tip_name not in tip_lookup:
            raise PipelineError(
                f"Calibration tip {node_calibration_tip_name} is missing from the input tree."
            )
        input_root_to_tip = float(tree.distance(tree.root, tip_lookup[node_calibration_tip_name]))
        if input_root_to_tip <= 0:
            raise PipelineError(
                f"Calibration tip {node_calibration_tip_name} has a non-positive root-to-tip depth."
            )
        return float(node_calibration_divergence_years) / input_root_to_tip, f"{input_root_to_tip:.12g}"
    raise PipelineError("Unsupported method. Expected molecular_clock or node_calibration.")


def iter_nonroot_clades(tree):
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        yield clade


def summarize(values: list[float]) -> tuple[float, float, float]:
    if not values:
        raise PipelineError("Tree does not contain any non-root branches to calibrate.")
    return min(values), statistics.median(values), max(values)


def load_node_hpd_registry(hpd_paths):
    """读取一个或多个 ``beast_node_hpd.tsv`` 构建 {tip_hash: 注释} 注册表。

    多源同 clade 时的优先级：按 ``hpd_paths`` 给定顺序后写覆盖。stage5 调用方应把
    backbone 文件排在 target 子树文件之前，从而让 target 子树（局部采样更密）覆盖
    backbone 估计。注册表值为 (height_subst, hpd_low_subst, hpd_high_subst, source_id)。
    """
    registry = {}
    for path in hpd_paths:
        path = Path(path)
        if not path.exists():
            continue
        for row in load_table(path):
            tip_hash = row.get("tip_hash", "").strip()
            if not tip_hash:
                continue

            def _to_float(key):
                raw = (row.get(key) or "").strip()
                if raw == "":
                    return None
                try:
                    return float(raw)
                except ValueError:
                    return None

            registry[tip_hash] = {
                "height_subst": _to_float("height_subst"),
                "hpd_low_subst": _to_float("hpd_low_subst"),
                "hpd_high_subst": _to_float("hpd_high_subst"),
                "source_id": (row.get("source_id") or "").strip(),
            }
    return registry


def _format_hpd_years(value):
    return "" if value is None else f"{value:.12g}"


def build_figtree_nexus(tree, annotation_by_clade_id, tree_name="merged_time_tree"):
    """把校准后的年代树写成 FigTree 风格 NEXUS（每节点带 height/HPD 注释）。

    ``annotation_by_clade_id`` 以 ``id(clade)`` 为键，值含 ``height_years``、
    ``hpd_low_years``、``hpd_high_years``（后两者可为 None）。枝长单位为年。
    Bio.Phylo 无法写出 ``[&...]`` 注释，故此处手工拼装 newick。
    """
    terminals = tree.get_terminals()
    taxon_names = [str(tip.name) for tip in terminals]
    label_by_name = {name: str(index) for index, name in enumerate(taxon_names, start=1)}

    def _node_comment(clade):
        info = annotation_by_clade_id.get(id(clade))
        if info is None:
            return ""
        parts = [f"height={info['height_years']:.12g}"]
        if info["hpd_low_years"] is not None and info["hpd_high_years"] is not None:
            parts.append(
                "height_95%_HPD={" f"{info['hpd_low_years']:.12g},{info['hpd_high_years']:.12g}" "}"
            )
        return "[&" + ",".join(parts) + "]"

    def _emit(clade):
        comment = _node_comment(clade)
        branch = "" if clade.branch_length is None else f":{float(clade.branch_length):.12g}"
        if clade.is_terminal():
            return f"{label_by_name[str(clade.name)]}{comment}{branch}"
        inner = ",".join(_emit(child) for child in clade.clades)
        return f"({inner}){comment}{branch}"

    newick_body = _emit(tree.root)
    translate_lines = ",\n".join(f"\t\t{index} {name}" for index, name in enumerate(taxon_names, start=1))
    lines = [
        "#NEXUS",
        "",
        "Begin taxa;",
        f"\tDimensions ntax={len(taxon_names)};",
        "\tTaxlabels",
    ]
    lines.extend(f"\t\t{name}" for name in taxon_names)
    lines.append("\t\t;")
    lines.append("End;")
    lines.append("")
    lines.append("Begin trees;")
    lines.append("\tTranslate")
    lines.append(translate_lines)
    lines.append("\t\t;")
    lines.append(f"tree {tree_name} = [&R] {newick_body};")
    lines.append("End;")
    return "\n".join(lines) + "\n"


def run(
    input_tree,
    output_dir,
    output_tree_name,
    method,
    substitution_rate_per_site_per_year,
    sequence_length,
    branch_length_unit="substitutions_per_site",
    node_calibration_tip_name="RSRS",
    node_calibration_divergence_years=180000,
    beast_node_hpd_files=None,
    nexus_output_name="merged_time_tree.nexus",
    node_table_name="node_time_table.tsv",
    log_level="INFO",
):
    input_tree = Path(input_tree).resolve()
    output_dir = Path(output_dir).resolve()
    output_tree_name = str(output_tree_name)
    method = str(method)
    substitution_rate_per_site_per_year = float(substitution_rate_per_site_per_year)
    sequence_length = int(sequence_length)
    node_calibration_divergence_years = float(node_calibration_divergence_years)
    beast_node_hpd_files = list(beast_node_hpd_files or [])
    nexus_output_name = str(nexus_output_name)
    node_table_name = str(node_table_name)

    output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("time_calibrate_tree", output_dir / "time_calibration.log", log_level)
    logger.info("Starting time calibration for %s", input_tree)

    tree = read_newick_tree(input_tree)
    validate_branch_lengths_complete(tree)
    scale_factor, node_calibration_input_root_to_tip = resolve_scale_factor(
        tree=tree,
        method=method,
        substitution_rate_per_site_per_year=substitution_rate_per_site_per_year,
        sequence_length=sequence_length,
        branch_length_unit=branch_length_unit,
        node_calibration_tip_name=node_calibration_tip_name,
        node_calibration_divergence_years=node_calibration_divergence_years,
    )

    calibrated_tree = read_newick_tree(input_tree)
    input_node_id_map, _, _ = assign_node_ids(tree)
    calibrated_node_id_map, _, _ = assign_node_ids(calibrated_tree)
    calibrated_lookup = {node_id: clade for clade, node_id in calibrated_node_id_map.items()}

    edge_rows = []
    input_values = []
    output_values = []
    for clade, node_id in input_node_id_map.items():
        if clade is tree.root:
            continue
        input_branch_length = 0.0 if clade.branch_length is None else float(clade.branch_length)
        if input_branch_length < 0:
            raise PipelineError(f"Negative branch length detected at node {node_id}: {input_branch_length}")
        calibrated_branch_length = input_branch_length * scale_factor
        calibrated_lookup[node_id].branch_length = calibrated_branch_length
        input_values.append(input_branch_length)
        output_values.append(calibrated_branch_length)
        edge_rows.append(
            {
                "node_id": node_id,
                "edge_role": "terminal" if clade.is_terminal() else "internal",
                "tip_name": "" if not clade.is_terminal() else str(clade.name),
                "input_branch_length": f"{input_branch_length:.12g}",
                "calibrated_branch_length_years": f"{calibrated_branch_length:.12g}",
            }
        )

    input_min, input_median, input_max = summarize(input_values)
    output_min, output_median, output_max = summarize(output_values)

    output_tree = output_dir / output_tree_name
    write_tree(calibrated_tree, output_tree)
    write_table(edge_rows, TIME_CALIBRATION_EDGE_COLUMNS, output_dir / "time_calibration_edge_report.tsv")
    write_table(
        [
            {
                "input_tree": input_tree.as_posix(),
                "output_tree": output_tree.as_posix(),
                "method": method,
                "branch_length_unit": branch_length_unit,
                "substitution_rate_per_site_per_year": f"{substitution_rate_per_site_per_year:.12g}",
                "sequence_length": str(sequence_length),
                "node_calibration_tip_name": node_calibration_tip_name,
                "node_calibration_divergence_years": f"{node_calibration_divergence_years:.12g}",
                "node_calibration_input_root_to_tip": node_calibration_input_root_to_tip,
                "scale_factor_years_per_input_branch_unit": f"{scale_factor:.12g}",
                "nonroot_edge_count": str(len(edge_rows)),
                "input_branch_min": f"{input_min:.12g}",
                "input_branch_median": f"{input_median:.12g}",
                "input_branch_max": f"{input_max:.12g}",
                "output_branch_min_years": f"{output_min:.12g}",
                "output_branch_median_years": f"{output_median:.12g}",
                "output_branch_max_years": f"{output_max:.12g}",
            }
        ],
        TIME_CALIBRATION_SUMMARY_COLUMNS,
        output_dir / "time_calibration_summary.tsv",
    )

    # ---------------------------------------------------------------------
    # 生成带节点时间注释的 BEAST 风格 NEXUS + 全节点时间信息 TSV（含 HPD）。
    # HPD 传播：固定拓扑下最终合并树每个内部节点 = 一个确定 clade，按 tip_hash
    # 映射回源 BEAST MCC 的高度 95% HPD（substitutions/site），再按本节点的
    # 高度变换比例 f = height_years / height_subst 线性传播到年。无法匹配者记 NA。
    # ---------------------------------------------------------------------
    registry = load_node_hpd_registry(beast_node_hpd_files)

    depths = calibrated_tree.depths()  # clade -> 根到该 clade 的枝长之和（年）
    terminal_depths = [depths[tip] for tip in calibrated_tree.get_terminals()]
    total_height_years = max(terminal_depths) if terminal_depths else 0.0

    annotation_by_clade_id = {}
    node_time_rows = []
    matched_count = 0
    for clade, node_id in calibrated_node_id_map.items():
        depth_from_root = float(depths.get(clade, 0.0))
        height_years = max(0.0, total_height_years - depth_from_root)
        tip_names = get_tip_names_from_clade(clade)
        n_tips = len(tip_names)
        if clade.is_terminal():
            node_type = "tip"
        elif clade is calibrated_tree.root:
            node_type = "root"
        else:
            node_type = "internal"
        branch_length = 0.0 if clade.branch_length is None else float(clade.branch_length)

        hpd_low_years = None
        hpd_high_years = None
        source_height_subst = None
        hpd_source = "unmatched"
        if n_tips >= 2:
            entry = registry.get(compute_tip_hash(tip_names))
            if entry is not None:
                source_height_subst = entry["height_subst"]
                if (
                    entry["hpd_low_subst"] is not None
                    and entry["hpd_high_subst"] is not None
                    and source_height_subst not in (None, 0.0)
                ):
                    factor = height_years / source_height_subst
                    hpd_low_years = entry["hpd_low_subst"] * factor
                    hpd_high_years = entry["hpd_high_subst"] * factor
                    hpd_source = entry["source_id"] or "matched"
                    matched_count += 1
                else:
                    hpd_source = (entry["source_id"] or "matched") + ":no_hpd"

        # 内部节点（含根）写入 NEXUS 注释；叶节点 height 恒为 0，不注释 HPD。
        if not clade.is_terminal():
            annotation_by_clade_id[id(clade)] = {
                "height_years": height_years,
                "hpd_low_years": hpd_low_years,
                "hpd_high_years": hpd_high_years,
            }

        node_time_rows.append(
            {
                "node_id": node_id,
                "node_type": node_type,
                "n_descendant_tips": str(n_tips),
                "height_years": f"{height_years:.12g}",
                "depth_from_root_years": f"{depth_from_root:.12g}",
                "hpd_low_years": _format_hpd_years(hpd_low_years),
                "hpd_high_years": _format_hpd_years(hpd_high_years),
                "source_height_subst_per_site": "" if source_height_subst is None else f"{source_height_subst:.12g}",
                "hpd_source": hpd_source,
                "branch_length_years": f"{branch_length:.12g}",
                "representative_tips": ";".join(sorted(tip_names)[:3]),
            }
        )

    node_time_rows.sort(key=lambda item: item["node_id"])
    write_table(node_time_rows, NODE_TIME_TABLE_COLUMNS, output_dir / node_table_name)

    nexus_text = build_figtree_nexus(
        calibrated_tree, annotation_by_clade_id, tree_name=Path(nexus_output_name).stem,
    )
    (output_dir / nexus_output_name).write_text(nexus_text, encoding="utf-8")

    internal_total = sum(1 for clade in calibrated_tree.find_clades() if not clade.is_terminal())
    logger.info(
        "Wrote annotated NEXUS (%s) and node time table (%s); HPD matched %d/%d internal nodes from %d source file(s).",
        nexus_output_name,
        node_table_name,
        matched_count,
        internal_total,
        len(beast_node_hpd_files),
    )

    logger.info(
        "Time calibration finished: method=%s unit=%s rate=%.12g seq_len=%d calibration_tip=%s calibration_years=%.12g scale_factor=%.12g",
        method,
        branch_length_unit,
        substitution_rate_per_site_per_year,
        sequence_length,
        node_calibration_tip_name,
        node_calibration_divergence_years,
        scale_factor,
    )
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Scale a tree from substitution units into years.")
    parser.add_argument("--input-tree", required=True, help="Input ultrametric tree.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--output-tree-name", default="merged_ultrametric_tree_years.nwk", help="Output tree filename.")
    parser.add_argument("--method", default="molecular_clock", help="molecular_clock or node_calibration.")
    parser.add_argument(
        "--substitution-rate-per-site-per-year",
        required=True,
        type=float,
        help="Substitution rate in substitutions/site/year.",
    )
    parser.add_argument("--sequence-length", required=True, type=int, help="Sequence length.")
    parser.add_argument(
        "--branch-length-unit",
        default="substitutions_per_site",
        help="substitutions_per_site or substitutions_per_sequence.",
    )
    parser.add_argument(
        "--node-calibration-tip-name",
        default="RSRS",
        help="Tip used as the calibrated singleton-outgroup lineage in node_calibration mode.",
    )
    parser.add_argument(
        "--node-calibration-divergence-years",
        default=180000,
        type=float,
        help="Target divergence time in years for the calibration tip versus the remaining taxa.",
    )
    parser.add_argument(
        "--beast-node-hpd",
        action="append",
        default=None,
        help=(
            "Path to a beast_node_hpd.tsv (stage2). May be given multiple times; later files "
            "override earlier ones for shared clades. List the backbone file first so that "
            "target subtrees take precedence."
        ),
    )
    parser.add_argument(
        "--nexus-output-name",
        default="merged_time_tree.nexus",
        help="Filename for the FigTree-style annotated NEXUS time tree.",
    )
    parser.add_argument(
        "--node-table-name",
        default="node_time_table.tsv",
        help="Filename for the per-node time information table.",
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            input_tree=args.input_tree,
            output_dir=args.output_dir,
            output_tree_name=args.output_tree_name,
            method=args.method,
            substitution_rate_per_site_per_year=args.substitution_rate_per_site_per_year,
            sequence_length=args.sequence_length,
            branch_length_unit=args.branch_length_unit,
            node_calibration_tip_name=args.node_calibration_tip_name,
            node_calibration_divergence_years=args.node_calibration_divergence_years,
            beast_node_hpd_files=args.beast_node_hpd,
            nexus_output_name=args.nexus_output_name,
            node_table_name=args.node_table_name,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

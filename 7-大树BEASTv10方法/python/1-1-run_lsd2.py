#!/usr/bin/env python3
"""运行 IQ-TREE 内置 LSD2 并输出 BEAST startingTree。"""

from __future__ import annotations

import argparse
import logging
import math
import shutil
from pathlib import Path

from beast_pipeline_utils import (
    PipelineError,
    check_monophyly,
    convert_tree_to_newick,
    parse_calibration_nodes,
    read_fasta_ids,
    read_membership,
    read_tree,
    run_iqtree_lsd2,
    tree_tip_names,
    validate_calibration_nesting,
    validate_starting_tree,
    write_filtered_fasta,
    write_lsd2_date_file,
    write_validation_table,
)


log = logging.getLogger(__name__)


def write_report(path: Path, values: dict[str, object]) -> None:
    """写 LSD2 步骤报告。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# LSD2 定年树构建报告",
        "",
        f"- FASTA 序列数: {values['fasta_count']}",
        f"- ML tree tips: {values['tree_tip_count']}",
        f"- FASTA 中未进入树的序列数: {values['fasta_not_in_tree']}",
    ]
    for name, count in values["node_counts"].items():  # type: ignore[union-attr]
        lines.append(f"- {name} 校准 tips: {count}")
    lines += [
        f"- LSD2 threads: {values['threads']}",
        f"- startingTree root height: {values['starting_tree_root_height']:.6f}",
        f"- startingTree: `{values['starting_tree_path']}`",
        f"- LSD2 NEXUS: `{values['timetree_nexus_path']}`",
        "",
        "IQ-TREE 内置 LSD2 使用固定 ML 拓扑和 --date-root 根年龄约束生成时间尺度 startingTree。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(
    *,
    fasta: str,
    ml_tree: str,
    calibration_membership: str,
    calibration_nodes: str,
    output_table: str,
    output_report: str,
    temp_dir: str,
    iqtree_bin: str,
    iqtree_model: str,
    iqtree_threads: str,
    iqtree_outgroup: str,
    iqtree_date_options: str,
    overwrite: str = "false",
) -> int:
    """运行 LSD2 步骤并写出第二步所需的正式产物。"""
    fasta_path = Path(fasta)
    ml_tree_path = Path(ml_tree)
    membership_path = Path(calibration_membership)
    output_table_path = Path(output_table)
    output_report_path = Path(output_report)
    temp_path = Path(temp_dir)
    for required in [fasta_path, ml_tree_path, membership_path]:
        if not required.exists():
            raise PipelineError(f"缺少输入文件: {required}")
    if overwrite.lower() == "true" and temp_path.exists():
        shutil.rmtree(temp_path)
    output_table_path.mkdir(parents=True, exist_ok=True)
    output_report_path.mkdir(parents=True, exist_ok=True)
    temp_path.mkdir(parents=True, exist_ok=True)

    nodes = parse_calibration_nodes(calibration_nodes)
    non_root = [n for n in nodes if not n["is_root"]]
    root_node = next(n for n in nodes if n["is_root"])
    if root_node["age"] is None:
        raise PipelineError("root 节点缺少 age，无法用 --date-root 锁定根年龄")

    fasta_ids, fasta_lengths = read_fasta_ids(fasta_path)
    membership, hap_map = read_membership(membership_path)
    ml_tree_obj = read_tree(ml_tree_path)
    taxa = tree_tip_names(ml_tree_obj)
    taxa_set = set(taxa)
    fasta_set = set(fasta_ids)
    if missing := sorted(taxa_set - fasta_set):
        raise PipelineError(f"树中样本不在 FASTA: {missing[:10]}")

    # 成员表 ↔ yaml 一致性：成员表中出现的节点必须在 yaml 中定义。
    yaml_names = {n["name"] for n in nodes}
    for tsv_node in membership:
        if tsv_node not in yaml_names:
            raise PipelineError(f"成员表中的节点 {tsv_node} 未在 yaml calibration.nodes 定义，请保持一致")

    # 每个非 root 节点的成员（与树 tip 取交集）。
    calibration_taxa: dict[str, list[str]] = {}
    for node in non_root:
        name = node["name"]
        if node["age"] is None:
            raise PipelineError(f"校准节点 {name} 缺少 age")
        if name not in membership:
            raise PipelineError(f"yaml 定义的校准节点 {name} 在成员表中无记录，请保持一致")
        members = sorted(set(membership[name]) & taxa_set)
        if not members:
            raise PipelineError(f"校准节点 {name} 在树中没有成员（成员表与树不匹配）")
        calibration_taxa[name] = members

    # 逐节点单系检查。
    for name, members in calibration_taxa.items():
        result = check_monophyly(ml_tree_path, set(members))
        if not result["is_monophyletic"]:
            conflict_path = output_table_path / f"{name}_non_monophyletic_extra_taxa.tsv"
            with conflict_path.open("w", encoding="utf-8") as handle:
                handle.write("sample_id\thaplogroup\n")
                for sample_id in result["extra_taxa"]:
                    handle.write(f"{sample_id}\t{hap_map.get(sample_id, 'NA')}\n")
            raise PipelineError(f"{name} 在 ML 树中不是单系，冲突样本写入 {conflict_path}")

    # 嵌套关系与年龄单调性校验。
    validate_calibration_nesting(nodes, {name: set(members) for name, members in calibration_taxa.items()})

    filtered_fasta = temp_path / "lsd2_input" / "alignment.tree_tips.fasta"
    write_filtered_fasta(fasta_path, taxa_set, filtered_fasta)
    date_file = output_table_path / "lsd2_dates.txt"
    node_ages = {node["name"]: node["age"] for node in non_root}
    write_lsd2_date_file(
        taxa=taxa,
        calibration_taxa=calibration_taxa,
        node_ages=node_ages,
        output_path=date_file,
    )
    lsd2_outdir = temp_path / "lsd2"
    timetree_nexus, lsd2_log = run_iqtree_lsd2(
        iqtree_bin=Path(iqtree_bin),
        alignment=filtered_fasta,
        ml_tree=ml_tree_path,
        date_file=date_file,
        root_age=root_node["age"],
        outdir=lsd2_outdir,
        model=iqtree_model,
        threads=iqtree_threads,
        outgroup=iqtree_outgroup,
        date_options=iqtree_date_options,
    )
    (output_table_path / "iqtree_lsd2.stdout_stderr.log").write_text(lsd2_log, encoding="utf-8")
    timetree_nexus_out = output_table_path / "lsd.timetree.nex"
    shutil.copy2(timetree_nexus, timetree_nexus_out)

    starting_tree_path = output_table_path / "startingTree.lsd2.nwk"
    convert_tree_to_newick(timetree_nexus, starting_tree_path)
    tree_validation = validate_starting_tree(ml_tree_path, starting_tree_path, taxa_set)

    validation_rows = [
        ("fasta_count", str(len(fasta_ids)), "PASS"),
        ("tree_tip_count", str(len(taxa)), "PASS"),
        ("membership_rows", str(sum(len(v) for v in membership.values())), "PASS"),
        ("fasta_not_in_tree", str(len(fasta_set - taxa_set)), "WARN"),
        ("sequence_lengths", ",".join(f"{k}:{v}" for k, v in sorted({v: list(fasta_lengths.values()).count(v) for v in set(fasta_lengths.values())}.items())), "PASS"),
    ]
    for name, members in calibration_taxa.items():
        validation_rows.append((f"{name}_monophyletic_tip_count", str(len(members)), "PASS"))
    validation_rows += [
        ("starting_tree_root_height", f"{tree_validation['root_height']:.6f}", "PASS" if math.isfinite(tree_validation["root_height"]) else "FAIL"),
        ("starting_tree_topology", "compatible", "PASS"),
        ("lsd2_threads", str(iqtree_threads), "PASS"),
    ]
    write_validation_table(output_table_path / "input_validation.tsv", validation_rows)
    write_report(
        output_report_path / "lsd2_report.md",
        {
            "fasta_count": len(fasta_ids),
            "tree_tip_count": len(taxa),
            "fasta_not_in_tree": len(fasta_set - taxa_set),
            "node_counts": {name: len(members) for name, members in calibration_taxa.items()},
            "threads": iqtree_threads,
            "starting_tree_root_height": tree_validation["root_height"],
            "starting_tree_path": starting_tree_path,
            "timetree_nexus_path": timetree_nexus_out,
        },
    )
    log.info("LSD2 startingTree 已生成: %s", starting_tree_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description="运行 IQ-TREE 内置 LSD2 生成 BEAST startingTree")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--ml-tree", required=True)
    parser.add_argument("--calibration-membership", required=True)
    parser.add_argument("--calibration-nodes", required=True)
    parser.add_argument("--output-table", required=True)
    parser.add_argument("--output-report", required=True)
    parser.add_argument("--temp-dir", required=True)
    parser.add_argument("--iqtree-bin", required=True)
    parser.add_argument("--iqtree-model", required=True)
    parser.add_argument("--iqtree-threads", required=True)
    parser.add_argument("--iqtree-outgroup", default="")
    parser.add_argument("--iqtree-date-options", default="")
    parser.add_argument("--overwrite", default="false")
    return parser


def main(argv: list[str] | None = None) -> int:
    """命令行入口。"""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(**vars(args))
    except PipelineError as exc:
        log.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

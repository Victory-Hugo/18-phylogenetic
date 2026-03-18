#!/usr/bin/env python3
"""Scale an ultrametric tree from substitution units into years."""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

from phylo_merge_common import validate_branch_lengths_complete
from phylo_split_common import PipelineError, assign_node_ids, read_newick_tree, setup_logger, write_table, write_tree


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
    log_level="INFO",
):
    input_tree = Path(input_tree).resolve()
    output_dir = Path(output_dir).resolve()
    output_tree_name = str(output_tree_name)
    method = str(method)
    substitution_rate_per_site_per_year = float(substitution_rate_per_site_per_year)
    sequence_length = int(sequence_length)
    node_calibration_divergence_years = float(node_calibration_divergence_years)

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
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

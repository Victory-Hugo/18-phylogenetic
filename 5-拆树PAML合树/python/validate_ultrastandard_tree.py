#!/usr/bin/env python3
"""Validate ultrametric projection outputs derived from merged_ml_tree."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import validate_branch_lengths_complete, validate_rooted_tree_with_outgroup
from phylo_split_common import PipelineError, get_tip_names_from_tree, load_table, setup_logger, write_validation_report
from phylo_ultrastandard_common import infer_unique_outgroup, read_tree, validate_ultrametric


def run(
    split_output_dir,
    ultrastandard_output_dir,
    ultrametric_tolerance=1e-8,
    log_level="INFO",
):
    split_output_dir = Path(split_output_dir).resolve()
    ultrastandard_output_dir = Path(ultrastandard_output_dir).resolve()
    ultrametric_tolerance = float(ultrametric_tolerance)

    logger = setup_logger("validate_ultrastandard_tree", None, log_level)
    logger.info("Starting ultrastandard validation.")

    outgroup_tip_name = infer_unique_outgroup(split_output_dir / "paml_subtree_summary.tsv")
    master_tree = read_tree(split_output_dir / "intermediate" / "rooted.tree")
    projection_input_tree = read_tree(ultrastandard_output_dir / "projection_input_tree.nwk")
    merged_ultra_tree_path = ultrastandard_output_dir / "merged_ultrametric_tree.nwk"
    merged_tree = read_tree(merged_ultra_tree_path if merged_ultra_tree_path.exists() else (ultrastandard_output_dir / "merged_tree.nwk"))
    projection_rows = load_table(ultrastandard_output_dir / "ultrametric_projection_report.tsv")

    report_rows = []
    validate_rooted_tree_with_outgroup(master_tree, outgroup_tip_name)
    report_rows.append(("master_rooted_topology", "PASS", "Master tree remains rooted with singleton outgroup"))
    validate_rooted_tree_with_outgroup(projection_input_tree, outgroup_tip_name)
    report_rows.append(("projection_input_rooted_topology", "PASS", "Projection input tree rooted topology verified"))
    validate_rooted_tree_with_outgroup(merged_tree, outgroup_tip_name)
    report_rows.append(("merged_rooted_topology", "PASS", "Merged ultrametric tree rooted topology verified"))

    validate_ultrametric(merged_tree, ultrametric_tolerance, "merged_ultrametric_tree")
    report_rows.append(("merged_ultrametric", "PASS", "Final merged tree is ultrametric"))

    master_tip_set = set(get_tip_names_from_tree(master_tree))
    input_tip_set = set(get_tip_names_from_tree(projection_input_tree))
    merged_tip_set = set(get_tip_names_from_tree(merged_tree))
    if input_tip_set != master_tip_set:
        missing = sorted(master_tip_set - input_tip_set)
        extra = sorted(input_tip_set - master_tip_set)
        report_rows.append(("projection_input_tip_set", "FAIL", f"missing={missing[:5]} extra={extra[:5]}"))
    else:
        report_rows.append(("projection_input_tip_set", "PASS", f"{len(input_tip_set)} input tips verified"))
    if merged_tip_set != master_tip_set:
        missing = sorted(master_tip_set - merged_tip_set)
        extra = sorted(merged_tip_set - master_tip_set)
        report_rows.append(("merged_tip_set", "FAIL", f"missing={missing[:5]} extra={extra[:5]}"))
    else:
        report_rows.append(("merged_tip_set", "PASS", f"{len(merged_tip_set)} merged tips verified"))

    validate_branch_lengths_complete(projection_input_tree)
    report_rows.append(("projection_input_branch_lengths_complete", "PASS", "Projection input branches are complete"))
    validate_branch_lengths_complete(merged_tree)
    report_rows.append(("merged_branch_lengths_complete", "PASS", "Projected branches are complete"))

    if not projection_rows:
        report_rows.append(("projection_report", "FAIL", "ultrametric_projection_report.tsv is empty"))
    else:
        changed_edges = sum(1 for row in projection_rows if row.get("changed") == "true")
        report_rows.append(("projection_report", "PASS", f"{len(projection_rows)} projected edges recorded; changed={changed_edges}"))

    report_path = ultrastandard_output_dir / "ultrastandard_validation_report.tsv"
    write_validation_report(report_rows, report_path)
    failures = [row for row in report_rows if row[1] == "FAIL"]
    if failures:
        details = "; ".join(f"{name}: {detail}" for name, _, detail in failures[:10])
        raise PipelineError(f"Ultrastandard validation failed with {len(failures)} issue(s): {details}")
    logger.info("Ultrastandard validation passed.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Validate ultrametric projection outputs.")
    parser.add_argument("--split-output-dir", required=True, help="Split output directory.")
    parser.add_argument("--ultrastandard-output-dir", required=True, help="Ultrastandard output directory.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            split_output_dir=args.split_output_dir,
            ultrastandard_output_dir=args.ultrastandard_output_dir,
            ultrametric_tolerance=args.ultrametric_tolerance,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

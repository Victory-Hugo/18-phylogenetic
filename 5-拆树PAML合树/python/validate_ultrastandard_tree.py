#!/usr/bin/env python3
"""Validate the ultrastandard merged tree and its intermediate ultrametric products."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import (
    load_backbone_summary,
    load_paml_summary,
    load_target_summary,
    validate_branch_lengths_complete,
    validate_rooted_tree_with_outgroup,
)
from phylo_split_common import PipelineError, get_tip_names_from_tree, setup_logger, write_validation_report
from phylo_ultrastandard_common import infer_unique_outgroup, read_tree, validate_ultrametric


def run(
    split_output_dir,
    paml_output_dir,
    ultrastandard_output_dir,
    ultrametric_tolerance=1e-8,
    backbone_ultrametric_tree=None,
    log_level="INFO",
):
    split_output_dir = Path(split_output_dir).resolve()
    paml_output_dir = Path(paml_output_dir).resolve()
    ultrastandard_output_dir = Path(ultrastandard_output_dir).resolve()
    ultrametric_tolerance = float(ultrametric_tolerance)

    logger = setup_logger("validate_ultrastandard_tree", None, log_level)
    logger.info("Starting ultrastandard validation.")

    paml_summary_tsv = split_output_dir / "paml_subtree_summary.tsv"
    outgroup_tip_name = infer_unique_outgroup(paml_summary_tsv)

    master_tree = read_tree(split_output_dir / "intermediate" / "rooted.tree")
    assembly_scaffold = read_tree(ultrastandard_output_dir / "assembly_scaffold.nwk")
    backbone_assigned_scaffold = read_tree(ultrastandard_output_dir / "backbone_assigned_scaffold.nwk")
    merged_tree = read_tree(ultrastandard_output_dir / "merged_tree.nwk")

    backbone_records = load_backbone_summary(split_output_dir / "backbone_summary.tsv")
    target_records = load_target_summary(split_output_dir / "target_subtree_summary.tsv")
    load_paml_summary(paml_summary_tsv)
    backbone_tip_names = [record.backbone_tip_name for record in sorted(backbone_records, key=lambda item: item.selection_rank)]
    backbone_ultra_path = (
        Path(backbone_ultrametric_tree).resolve()
        if backbone_ultrametric_tree not in (None, "", "null")
        else (paml_output_dir / "backbone_analysis" / "backbone_ultrametric_tree.nwk")
    )
    backbone_tree = read_tree(backbone_ultra_path)

    report_rows = []
    validate_rooted_tree_with_outgroup(master_tree, outgroup_tip_name)
    report_rows.append(("master_rooted_topology", "PASS", "Master tree remains rooted with singleton outgroup"))
    validate_rooted_tree_with_outgroup(assembly_scaffold, outgroup_tip_name)
    report_rows.append(("assembly_scaffold_rooted_topology", "PASS", "Assembly scaffold rooted topology verified"))
    validate_rooted_tree_with_outgroup(backbone_assigned_scaffold, outgroup_tip_name)
    report_rows.append(("backbone_assigned_scaffold_rooted_topology", "PASS", "Backbone-assigned scaffold topology verified"))
    validate_rooted_tree_with_outgroup(merged_tree, outgroup_tip_name)
    report_rows.append(("merged_rooted_topology", "PASS", "Merged tree rooted topology verified"))

    validate_ultrametric(backbone_tree, ultrametric_tolerance, "backbone_ultrametric_tree")
    report_rows.append(("backbone_ultrametric", "PASS", "Backbone ultrametric tree verified"))
    validate_ultrametric(merged_tree, ultrametric_tolerance, "ultrastandard_merged_tree")
    report_rows.append(("merged_ultrametric", "PASS", "Final merged tree is ultrametric"))

    expected_backbone_tip_set = set(backbone_tip_names) | {outgroup_tip_name}
    actual_backbone_tip_set = set(get_tip_names_from_tree(backbone_tree))
    if actual_backbone_tip_set != expected_backbone_tip_set:
        missing = sorted(expected_backbone_tip_set - actual_backbone_tip_set)
        extra = sorted(actual_backbone_tip_set - expected_backbone_tip_set)
        report_rows.append(("backbone_tip_set", "FAIL", f"missing={missing[:5]} extra={extra[:5]}"))
    else:
        report_rows.append(("backbone_tip_set", "PASS", f"{len(actual_backbone_tip_set)} backbone tips verified"))

    master_tip_set = set(get_tip_names_from_tree(master_tree))
    merged_tip_set = set(get_tip_names_from_tree(merged_tree))
    if merged_tip_set != master_tip_set:
        missing = sorted(master_tip_set - merged_tip_set)
        extra = sorted(merged_tip_set - master_tip_set)
        report_rows.append(("merged_tip_set", "FAIL", f"missing={missing[:5]} extra={extra[:5]}"))
    else:
        report_rows.append(("merged_tip_set", "PASS", f"{len(merged_tip_set)} merged tips verified"))

    validate_branch_lengths_complete(merged_tree)
    report_rows.append(("merged_branch_lengths_complete", "PASS", "All merged non-root branches are non-negative"))

    relative_target_dir = ultrastandard_output_dir / "relative_target_trees"
    relative_failures = []
    for target_record in target_records:
        relative_path = relative_target_dir / f"{target_record.target_subtree_id}.relative_ultrametric.nwk"
        if not relative_path.exists():
            relative_failures.append(f"Missing {relative_path.name}")
            continue
        try:
            relative_tree = read_tree(relative_path)
            validate_ultrametric(relative_tree, ultrametric_tolerance, target_record.target_subtree_id)
        except PipelineError as exc:
            relative_failures.append(str(exc))
    if relative_failures:
        report_rows.append(("relative_target_trees_ultrametric", "FAIL", "; ".join(relative_failures[:5])))
    else:
        report_rows.append(("relative_target_trees_ultrametric", "PASS", f"{len(target_records)} relative target trees verified"))

    report_path = ultrastandard_output_dir / "ultrastandard_validation_report.tsv"
    write_validation_report(report_rows, report_path)
    failures = [row for row in report_rows if row[1] == "FAIL"]
    if failures:
        details = "; ".join(f"{name}: {detail}" for name, _, detail in failures[:10])
        raise PipelineError(f"Ultrastandard validation failed with {len(failures)} issue(s): {details}")
    logger.info("Ultrastandard validation passed.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Validate ultrastandard merged tree outputs.")
    parser.add_argument("--split-output-dir", required=True, help="Split output directory.")
    parser.add_argument("--paml-output-dir", required=True, help="PAML output directory.")
    parser.add_argument("--ultrastandard-output-dir", required=True, help="Ultrastandard output directory.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--backbone-ultrametric-tree", default=None, help="Optional external backbone ultrametric tree.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            split_output_dir=args.split_output_dir,
            paml_output_dir=args.paml_output_dir,
            ultrastandard_output_dir=args.ultrastandard_output_dir,
            ultrametric_tolerance=args.ultrametric_tolerance,
            backbone_ultrametric_tree=args.backbone_ultrametric_tree,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Project merged_ml_tree onto an ultrametric tree without changing topology."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import (
    validate_branch_lengths_complete,
    validate_rooted_tree_with_outgroup,
    validate_tree_against_expected_tip_set,
    write_rows,
)
from phylo_split_common import PipelineError, assign_node_ids, get_tip_names_from_tree, setup_logger
from phylo_ultrastandard_common import (
    build_root_to_tip_row,
    infer_unique_outgroup,
    project_tree_to_ultrametric,
    read_tree,
    validate_ultrametric,
    write_root_to_tip_report,
    write_simple_tree,
)


ULTRAMETRIC_PROJECTION_REPORT_COLUMNS = [
    "node_id",
    "edge_role",
    "tip_name",
    "old_branch_length",
    "projected_branch_length",
    "delta",
    "changed",
    "details",
]


def resolve_primary_tree_input(merge_output_dir: Path, configured_path: str | None) -> Path:
    if configured_path not in (None, "", "null"):
        path = Path(configured_path)
        return path if path.is_absolute() else path.resolve()
    default_path = merge_output_dir / "merged_ml_tree.nwk"
    if default_path.exists():
        return default_path
    legacy_path = merge_output_dir / "merged_tree.nwk"
    return legacy_path


def build_projection_rows(input_tree, projected_tree, projection_mode: str):
    input_node_id_map, _, _ = assign_node_ids(input_tree)
    _, projected_reverse_map, _ = assign_node_ids(projected_tree)
    rows = []
    for clade, node_id in input_node_id_map.items():
        if clade is input_tree.root:
            continue
        projected_clade = projected_reverse_map.get(node_id)
        if projected_clade is None:
            raise PipelineError(f"Projected tree is missing node id {node_id}")
        old_branch_length = 0.0 if clade.branch_length is None else float(clade.branch_length)
        new_branch_length = 0.0 if projected_clade.branch_length is None else float(projected_clade.branch_length)
        rows.append(
            {
                "node_id": node_id,
                "edge_role": "terminal" if clade.is_terminal() else "internal",
                "tip_name": "" if not clade.is_terminal() else str(clade.name),
                "old_branch_length": f"{old_branch_length:.12g}",
                "projected_branch_length": f"{new_branch_length:.12g}",
                "delta": f"{(new_branch_length - old_branch_length):.12g}",
                "changed": "true" if abs(new_branch_length - old_branch_length) > 0 else "false",
                "details": (
                    "Terminal-extension ultrametric projection."
                    if projection_mode == "extend_terminal_to_max_depth" and clade.is_terminal()
                    else (
                        "Internal branch preserved."
                        if projection_mode == "extend_terminal_to_max_depth"
                        else "Constrained-height ultrametric projection."
                    )
                ),
            }
        )
    return rows


def run(
    split_output_dir,
    merge_output_dir,
    ultrastandard_output_dir,
    min_branch_length=1e-8,
    ultrametric_tolerance=1e-8,
    projection_mode="extend_terminal_to_max_depth",
    primary_tree_input=None,
    log_level="INFO",
):
    if projection_mode not in {"extend_terminal_to_max_depth", "constrained_height_fit"}:
        raise PipelineError(f"Unsupported ultrastandard.projection_mode: {projection_mode}")

    split_output_dir = Path(split_output_dir).resolve()
    merge_output_dir = Path(merge_output_dir).resolve()
    ultrastandard_output_dir = Path(ultrastandard_output_dir).resolve()
    min_branch_length = float(min_branch_length)
    ultrametric_tolerance = float(ultrametric_tolerance)
    ultrastandard_output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logger("merge_ultrastandard_tree", ultrastandard_output_dir / "ultrastandard.log", log_level)
    logger.info("Starting ultrametric projection from merged_ml_tree.")

    paml_summary_tsv = split_output_dir / "beast_subtree_summary.tsv"
    outgroup_tip_name = infer_unique_outgroup(paml_summary_tsv)
    master_tree = read_tree(split_output_dir / "intermediate" / "rooted.tree")
    input_tree_path = resolve_primary_tree_input(merge_output_dir, primary_tree_input)
    if not input_tree_path.exists():
        raise PipelineError(f"Primary tree input not found: {input_tree_path}")
    merged_ml_tree = read_tree(input_tree_path)

    validate_rooted_tree_with_outgroup(master_tree, outgroup_tip_name)
    validate_rooted_tree_with_outgroup(merged_ml_tree, outgroup_tip_name)
    validate_branch_lengths_complete(merged_ml_tree)
    validate_tree_against_expected_tip_set(
        merged_ml_tree,
        get_tip_names_from_tree(master_tree),
        "merged_ml_tree",
        outgroup_tip_name,
    )

    write_simple_tree(ultrastandard_output_dir / "projection_input_tree.nwk", merged_ml_tree)
    root_to_tip_rows = [build_root_to_tip_row("merged_ml_tree", merged_ml_tree, ultrametric_tolerance)]

    projected_tree, target_height, normalization_rows = project_tree_to_ultrametric(
        merged_ml_tree,
        projection_mode=projection_mode,
        min_branch_length=min_branch_length,
        tolerance=ultrametric_tolerance,
    )
    validate_rooted_tree_with_outgroup(projected_tree, outgroup_tip_name)
    validate_branch_lengths_complete(projected_tree)
    validate_ultrametric(projected_tree, ultrametric_tolerance, "merged_ultrametric_tree")

    projection_rows = build_projection_rows(merged_ml_tree, projected_tree, projection_mode=projection_mode)
    write_rows(
        projection_rows,
        ULTRAMETRIC_PROJECTION_REPORT_COLUMNS,
        ultrastandard_output_dir / "ultrametric_projection_report.tsv",
    )

    root_to_tip_rows.append(normalization_rows["pre"])
    root_to_tip_rows.append(normalization_rows["post"])
    root_to_tip_rows.append(build_root_to_tip_row("merged_ultrametric_tree", projected_tree, ultrametric_tolerance))
    write_root_to_tip_report(root_to_tip_rows, ultrastandard_output_dir / "root_to_tip_report.tsv")

    write_simple_tree(ultrastandard_output_dir / "merged_ultrametric_tree.nwk", projected_tree)
    write_simple_tree(ultrastandard_output_dir / "merged_tree.nwk", projected_tree)
    logger.info(
        "Ultrametric projection finished: input=%s target_height=%.12g changed_edges=%d",
        input_tree_path,
        target_height,
        sum(1 for row in projection_rows if row["changed"] == "true"),
    )
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Project merged_ml_tree onto an ultrametric tree.")
    parser.add_argument("--split-output-dir", required=True, help="Output directory from split pipeline.")
    parser.add_argument("--merge-output-dir", required=True, help="Output directory from merge pipeline.")
    parser.add_argument("--ultrastandard-output-dir", required=True, help="Ultrastandard output directory.")
    parser.add_argument("--min-branch-length", default=1e-8, type=float, help="Minimum positive branch length.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--projection-mode", default="extend_terminal_to_max_depth", help="Ultrametric projection mode.")
    parser.add_argument("--primary-tree-input", default=None, help="Optional path to merged_ml_tree input.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            split_output_dir=args.split_output_dir,
            merge_output_dir=args.merge_output_dir,
            ultrastandard_output_dir=args.ultrastandard_output_dir,
            min_branch_length=args.min_branch_length,
            ultrametric_tolerance=args.ultrametric_tolerance,
            projection_mode=args.projection_mode,
            primary_tree_input=args.primary_tree_input,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

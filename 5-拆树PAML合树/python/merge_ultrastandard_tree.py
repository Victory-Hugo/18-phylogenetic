#!/usr/bin/env python3
"""Merge subtree PAML outputs under a backbone ultrametric time framework."""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

from Bio.Phylo.BaseTree import Tree

from phylo_merge_common import (
    build_analysis_tree_path,
    build_scaffold_tree,
    extract_monophyletic_target_clade,
    load_backbone_summary,
    load_paml_summary,
    load_target_summary,
    validate_branch_lengths_complete,
    validate_rooted_tree_with_outgroup,
    validate_tree_against_expected_tip_set,
)
from phylo_split_common import PipelineError, assign_node_ids, get_tip_names_from_tree, setup_logger
from phylo_ultrastandard_common import (
    BACKBONE_EDGE_ASSIGNMENT_COLUMNS,
    GRAFT_REPORT_COLUMNS,
    TARGET_SCALING_REPORT_COLUMNS,
    assign_backbone_lengths_to_scaffold,
    build_root_to_tip_row,
    collapse_tree,
    count_backbone_descendants,
    extend_terminal_branches_to_height,
    infer_unique_outgroup,
    normalize_relative_ultrametric_tree,
    replace_placeholder_with_scaled_tree,
    resolve_backbone_ultrametric_tree_path,
    scale_tree_nonroot,
    tree_height,
    validate_ultrametric,
    write_root_to_tip_report,
    write_simple_tree,
)
from phylo_merge_common import ensure_positive_branch_length, write_rows


def run(
    split_output_dir,
    paml_output_dir,
    ultrastandard_output_dir,
    min_branch_length=1e-8,
    ultrametric_tolerance=1e-8,
    relative_ultrametric_method="extend_terminal_to_max_depth",
    require_backbone_descendant_anchor=True,
    allow_multifurcation=True,
    backbone_ultrametric_tree=None,
    log_level="INFO",
):
    if relative_ultrametric_method != "extend_terminal_to_max_depth":
        raise PipelineError(f"Unsupported ultrastandard.relative_ultrametric_method: {relative_ultrametric_method}")
    del allow_multifurcation

    split_output_dir = Path(split_output_dir).resolve()
    paml_output_dir = Path(paml_output_dir).resolve()
    ultrastandard_output_dir = Path(ultrastandard_output_dir).resolve()
    min_branch_length = float(min_branch_length)
    ultrametric_tolerance = float(ultrametric_tolerance)

    logger = setup_logger("merge_ultrastandard_tree", ultrastandard_output_dir / "ultrastandard.log", log_level)
    logger.info("Starting ultrastandard backbone-anchored merge.")

    master_tree_path = split_output_dir / "intermediate" / "rooted.tree"
    paml_summary_tsv = split_output_dir / "paml_subtree_summary.tsv"
    outgroup_tip_name = infer_unique_outgroup(paml_summary_tsv)

    master_tree = read_tree(master_tree_path)
    validate_rooted_tree_with_outgroup(master_tree, outgroup_tip_name)

    backbone_records = load_backbone_summary(split_output_dir / "backbone_summary.tsv")
    target_records = load_target_summary(split_output_dir / "target_subtree_summary.tsv")
    paml_records = load_paml_summary(paml_summary_tsv)
    paml_by_target_id = {record.target_subtree_id: record for record in paml_records}
    backbone_tip_names = [record.backbone_tip_name for record in sorted(backbone_records, key=lambda item: item.selection_rank)]
    backbone_tip_set = set(backbone_tip_names)

    backbone_ultrametric_tree_path = resolve_backbone_ultrametric_tree_path(paml_output_dir, backbone_ultrametric_tree)
    backbone_tree = read_tree(backbone_ultrametric_tree_path)
    validate_rooted_tree_with_outgroup(backbone_tree, outgroup_tip_name)
    validate_tree_against_expected_tip_set(
        backbone_tree,
        backbone_tip_names + [outgroup_tip_name],
        "backbone_ultrametric_tree",
        outgroup_tip_name,
    )
    validate_ultrametric(backbone_tree, ultrametric_tolerance, "backbone_ultrametric_tree")

    for target_record in target_records:
        target_tree = read_tree(split_output_dir / target_record.target_tree_file)
        backbone_descendants = count_backbone_descendants(target_tree, backbone_tip_set)
        if require_backbone_descendant_anchor and backbone_descendants < 1:
            raise PipelineError(
                f"{target_record.target_subtree_id} has no backbone descendants in its target subtree and cannot be anchored."
            )

    scaffold_tree, placeholder_to_target_id = build_scaffold_tree(
        master_tree=master_tree,
        target_records=target_records,
        backbone_tip_names=backbone_tip_names,
        outgroup_tip_name=outgroup_tip_name,
    )
    ultrastandard_output_dir.mkdir(parents=True, exist_ok=True)
    relative_target_dir = ultrastandard_output_dir / "relative_target_trees"
    relative_target_dir.mkdir(parents=True, exist_ok=True)
    write_simple_tree(ultrastandard_output_dir / "assembly_scaffold.nwk", scaffold_tree)

    node_id_map, _, _ = assign_node_ids(scaffold_tree)
    backbone_assignment_rows = assign_backbone_lengths_to_scaffold(
        scaffold_tree=scaffold_tree,
        backbone_tree=backbone_tree,
        backbone_tip_set=backbone_tip_set,
        outgroup_tip_name=outgroup_tip_name,
        min_branch_length=min_branch_length,
        node_id_map=node_id_map,
    )
    write_rows(
        backbone_assignment_rows,
        BACKBONE_EDGE_ASSIGNMENT_COLUMNS,
        ultrastandard_output_dir / "backbone_edge_assignment_report.tsv",
    )

    backbone_terminals = [
        tip for tip in scaffold_tree.get_terminals() if str(tip.name) in backbone_tip_set or str(tip.name) == outgroup_tip_name
    ]
    H = max(float(scaffold_tree.distance(scaffold_tree.root, tip)) for tip in backbone_terminals)
    for tip in backbone_terminals:
        current_depth = float(scaffold_tree.distance(scaffold_tree.root, tip))
        deficit = H - current_depth
        tip.branch_length = ensure_positive_branch_length(tip.branch_length, min_branch_length) + max(0.0, deficit)

    write_simple_tree(ultrastandard_output_dir / "backbone_assigned_scaffold.nwk", scaffold_tree)

    H = max(
        float(scaffold_tree.distance(scaffold_tree.root, tip))
        for tip in scaffold_tree.get_terminals()
        if str(tip.name) in backbone_tip_set or str(tip.name) == outgroup_tip_name
    )
    scaling_rows = []
    graft_rows = []
    root_to_tip_rows = [build_root_to_tip_row("backbone_ultrametric_tree", backbone_tree, ultrametric_tolerance)]

    for placeholder_name, target_id in sorted(placeholder_to_target_id.items()):
        target_record = next(record for record in target_records if record.target_subtree_id == target_id)
        paml_record = paml_by_target_id[target_id]
        analysis_tree_path = build_analysis_tree_path(
            merge_output_dir=ultrastandard_output_dir,
            paml_subtree_id=paml_record.paml_subtree_id,
            analysis_tree_source="external",
            external_result_dir=paml_output_dir / "analysis_trees",
        )
        analysis_tree = read_tree(analysis_tree_path)
        validate_tree_against_expected_tip_set(
            analysis_tree,
            paml_record.total_tip_names,
            paml_record.paml_subtree_id,
            outgroup_tip_name,
        )

        extracted_clade = extract_monophyletic_target_clade(analysis_tree, target_record.target_nonbackbone_tip_names)
        target_tree = collapse_tree(Tree(root=copy.deepcopy(extracted_clade), rooted=True))
        relative_tree, relative_height, _ = normalize_relative_ultrametric_tree(
            target_tree,
            min_branch_length=min_branch_length,
            tolerance=ultrametric_tolerance,
        )
        validate_ultrametric(relative_tree, ultrametric_tolerance, f"{target_id}.relative_tree")
        relative_target_path = relative_target_dir / f"{target_record.target_subtree_id}.relative_ultrametric.nwk"
        write_simple_tree(relative_target_path, relative_tree)
        root_to_tip_rows.append(build_root_to_tip_row(f"{target_id}.relative_ultrametric", relative_tree, ultrametric_tolerance))

        parent_clade, _ = find_placeholder_parent(scaffold_tree, placeholder_name)
        anchor_height = float(scaffold_tree.distance(scaffold_tree.root, parent_clade))
        available_depth = H - anchor_height
        if available_depth <= 0 and abs(available_depth) <= ultrametric_tolerance:
            available_depth = min_branch_length
        if available_depth <= 0:
            raise PipelineError(
                f"{target_id} anchor height is not below the tree height: anchor={anchor_height:.12g}, H={H:.12g}"
            )
        scaled_tree = scale_tree_nonroot(relative_tree, available_depth, min_branch_length=min_branch_length)
        scaled_root = copy.deepcopy(scaled_tree.root)
        scaled_root.branch_length = 0.0
        replace_placeholder_with_scaled_tree(scaffold_tree, placeholder_name, scaled_root)

        scaling_rows.append(
            {
                "target_subtree_id": target_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "attachment_placeholder": placeholder_name,
                "anchor_node_id": target_record.target_root_node_id,
                "anchor_height": f"{anchor_height:.12g}",
                "tree_height": f"{H:.12g}",
                "available_depth": f"{available_depth:.12g}",
                "relative_target_height": f"{relative_height:.12g}",
                "scaling_factor": f"{available_depth:.12g}",
                "status": "scaled",
                "details": "Relative ultrametric target tree scaled to anchor depth.",
            }
        )
        graft_rows.append(
            {
                "target_subtree_id": target_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "attachment_placeholder": placeholder_name,
                "anchor_node_id": target_record.target_root_node_id,
                "graft_status": "success",
                "details": "Scaled relative ultrametric target grafted onto backbone scaffold.",
            }
        )
        root_to_tip_rows.append(build_root_to_tip_row(f"{target_id}.scaled_tree", scaled_tree, ultrametric_tolerance))
        logger.info(
            "Grafted %s at anchor=%s available_depth=%.12g",
            target_id,
            target_record.target_root_node_id,
            available_depth,
        )

    scaffold_tree = collapse_tree(scaffold_tree)
    validate_rooted_tree_with_outgroup(scaffold_tree, outgroup_tip_name)
    validate_branch_lengths_complete(scaffold_tree)

    master_tip_set = set(get_tip_names_from_tree(master_tree))
    merged_tip_set = set(get_tip_names_from_tree(scaffold_tree))
    if merged_tip_set != master_tip_set:
        missing = sorted(master_tip_set - merged_tip_set)
        extra = sorted(merged_tip_set - master_tip_set)
        raise PipelineError(f"Ultrastandard merged tip set mismatch. missing={missing[:10]} extra={extra[:10]}")

    scaffold_tree, _, _ = extend_terminal_branches_to_height(
        scaffold_tree,
        min_branch_length=min_branch_length,
        tolerance=ultrametric_tolerance,
    )
    validate_ultrametric(scaffold_tree, ultrametric_tolerance, "ultrastandard_merged_tree")
    write_rows(scaling_rows, TARGET_SCALING_REPORT_COLUMNS, ultrastandard_output_dir / "target_scaling_report.tsv")
    write_rows(graft_rows, GRAFT_REPORT_COLUMNS, ultrastandard_output_dir / "graft_report.tsv")
    root_to_tip_rows.append(build_root_to_tip_row("final_merged_tree", scaffold_tree, ultrametric_tolerance))
    write_root_to_tip_report(root_to_tip_rows, ultrastandard_output_dir / "root_to_tip_report.tsv")
    write_simple_tree(ultrastandard_output_dir / "merged_tree.nwk", scaffold_tree)
    logger.info("Ultrastandard merge finished with %d grafted target clades.", len(graft_rows))
    return 0


def read_tree(path: Path):
    from phylo_ultrastandard_common import read_tree as _read_tree

    return _read_tree(path)


def find_placeholder_parent(tree, placeholder_name):
    from phylo_ultrastandard_common import find_placeholder_parent as _find_placeholder_parent

    return _find_placeholder_parent(tree, placeholder_name)


def build_parser():
    parser = argparse.ArgumentParser(description="Merge subtree PAML outputs under a backbone ultrametric framework.")
    parser.add_argument("--split-output-dir", required=True, help="Output directory from split pipeline.")
    parser.add_argument("--paml-output-dir", required=True, help="Output directory from paml pipeline.")
    parser.add_argument("--ultrastandard-output-dir", required=True, help="Ultrastandard output directory.")
    parser.add_argument("--min-branch-length", default=1e-8, type=float, help="Minimum positive branch length.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--relative-ultrametric-method", default="extend_terminal_to_max_depth", help="Relative ultrametric method.")
    parser.add_argument("--require-backbone-descendant-anchor", default="true", help="Require at least one backbone descendant in each target subtree.")
    parser.add_argument("--allow-multifurcation", default="true", help="Allow multifurcations in the final tree.")
    parser.add_argument("--backbone-ultrametric-tree", default=None, help="Optional path to an existing backbone ultrametric tree.")
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
            min_branch_length=args.min_branch_length,
            ultrametric_tolerance=args.ultrametric_tolerance,
            relative_ultrametric_method=args.relative_ultrametric_method,
            require_backbone_descendant_anchor=str(args.require_backbone_descendant_anchor).lower() == "true",
            allow_multifurcation=str(args.allow_multifurcation).lower() == "true",
            backbone_ultrametric_tree=args.backbone_ultrametric_tree,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

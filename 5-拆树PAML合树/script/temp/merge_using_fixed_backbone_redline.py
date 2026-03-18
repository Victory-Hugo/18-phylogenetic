#!/usr/bin/env python3
"""Redline-constrained merge on a fixed ultrametric backbone."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path
from typing import Dict, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = PROJECT_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from phylo_merge_common import (  # noqa: E402
    build_analysis_tree_path,
    build_scaffold_tree,
    ensure_positive_branch_length,
    extract_monophyletic_target_clade,
    get_tip_names_from_clade,
    get_tip_names_from_tree,
    load_paml_summary,
    load_target_summary,
    replace_placeholder_with_clade,
    scale_tree_branch_lengths,
    validate_rooted_tree_with_outgroup,
    validate_tree_against_expected_tip_set,
    write_tree_file,
)
from phylo_split_common import PipelineError, clone_clade, read_newick_tree  # noqa: E402


def write_tsv(rows, fieldnames, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_root_to_tip(tree, excluded_tip_names=None):
    excluded_tip_names = set() if excluded_tip_names is None else {str(name) for name in excluded_tip_names}
    depths = [
        float(tree.distance(tree.root, tip))
        for tip in tree.get_terminals()
        if str(tip.name) not in excluded_tip_names
    ]
    return {
        "n_tips": len(depths),
        "min_depth": min(depths),
        "median_depth": statistics.median(depths),
        "max_depth": max(depths),
        "delta": max(depths) - min(depths),
    }


def build_tip_depth_map(tree) -> Dict[str, float]:
    return {str(tip.name): float(tree.distance(tree.root, tip)) for tip in tree.get_terminals()}


def estimate_scale_factor_from_shared_backbone_tips(reference_tree, analysis_tree, shared_tip_names):
    unique_tip_names = sorted({str(name) for name in shared_tip_names})
    if len(unique_tip_names) < 2:
        return {
            "scale_factor": 1.0,
            "shared_tip_n": len(unique_tip_names),
            "shared_pair_n": 0,
            "median_log_ratio": 0.0,
            "residual_mad": 0.0,
            "status": "insufficient_backbone_overlap",
        }

    log_ratios = []
    for index, left_tip in enumerate(unique_tip_names):
        for right_tip in unique_tip_names[index + 1 :]:
            ref_distance = float(reference_tree.distance(left_tip, right_tip))
            qry_distance = float(analysis_tree.distance(left_tip, right_tip))
            if ref_distance <= 1e-12 or qry_distance <= 1e-12:
                continue
            log_ratios.append(math.log(ref_distance / qry_distance))

    if not log_ratios:
        return {
            "scale_factor": 1.0,
            "shared_tip_n": len(unique_tip_names),
            "shared_pair_n": 0,
            "median_log_ratio": 0.0,
            "residual_mad": 0.0,
            "status": "no_usable_backbone_pairs",
        }

    median_log_ratio = float(statistics.median(log_ratios))
    residual_mad = float(statistics.median(abs(value - median_log_ratio) for value in log_ratios))
    return {
        "scale_factor": math.exp(median_log_ratio),
        "shared_tip_n": len(unique_tip_names),
        "shared_pair_n": len(log_ratios),
        "median_log_ratio": median_log_ratio,
        "residual_mad": residual_mad,
        "status": "scaled_by_shared_backbone_pairs",
    }


def get_descendant_backbone_tips(clade, backbone_tip_set):
    return [tip for tip in get_tip_names_from_clade(clade) if tip in backbone_tip_set]


def get_single_tip_parent_depth(backbone_tree, target_tip, min_branch_length):
    stack = [(backbone_tree.root, 0.0)]
    while stack:
        parent, parent_depth = stack.pop()
        for child in parent.clades:
            child_branch = ensure_positive_branch_length(child.branch_length, min_branch_length)
            if child.is_terminal() and str(child.name) == target_tip:
                return float(parent_depth), float(parent_depth + child_branch)
            if not child.is_terminal():
                stack.append((child, parent_depth + child_branch))
    raise PipelineError(f"Could not locate backbone tip {target_tip}")


def compute_attachment_depth_from_backbone(clade, fixed_backbone_tree, backbone_tip_set, min_branch_length) -> Optional[float]:
    descendant_backbone = get_descendant_backbone_tips(clade, backbone_tip_set)
    if not descendant_backbone:
        return None
    if len(descendant_backbone) == 1:
        parent_depth, _ = get_single_tip_parent_depth(
            fixed_backbone_tree,
            descendant_backbone[0],
            min_branch_length=min_branch_length,
        )
        return parent_depth
    mrca = fixed_backbone_tree.common_ancestor(descendant_backbone)
    return float(fixed_backbone_tree.distance(fixed_backbone_tree.root, mrca))


def compute_backbone_node_depth(clade, fixed_backbone_tree, backbone_tip_set, min_branch_length) -> Optional[float]:
    descendant_backbone = get_descendant_backbone_tips(clade, backbone_tip_set)
    if not descendant_backbone:
        return None
    if clade.is_terminal() and len(descendant_backbone) == 1:
        _, tip_depth = get_single_tip_parent_depth(
            fixed_backbone_tree,
            descendant_backbone[0],
            min_branch_length=min_branch_length,
        )
        return tip_depth
    if len(descendant_backbone) == 1:
        return compute_attachment_depth_from_backbone(
            clade,
            fixed_backbone_tree=fixed_backbone_tree,
            backbone_tip_set=backbone_tip_set,
            min_branch_length=min_branch_length,
        )
    mrca = fixed_backbone_tree.common_ancestor(descendant_backbone)
    return float(fixed_backbone_tree.distance(fixed_backbone_tree.root, mrca))


def rewrite_scaffold_backbone_depths(scaffold_tree, fixed_backbone_tree, backbone_tip_set, min_branch_length):
    desired_depths = {scaffold_tree.root: 0.0}
    stack = [scaffold_tree.root]
    while stack:
        parent = stack.pop()
        parent_depth = desired_depths[parent]
        for child in parent.clades:
            if child.is_terminal() and str(child.name).startswith("TARGET_"):
                child.branch_length = None
                continue

            child_depth = compute_backbone_node_depth(
                child,
                fixed_backbone_tree=fixed_backbone_tree,
                backbone_tip_set=backbone_tip_set,
                min_branch_length=min_branch_length,
            )
            if child_depth is None:
                child.branch_length = float(min_branch_length)
                desired_depths[child] = parent_depth + float(min_branch_length)
            else:
                branch_length = child_depth - parent_depth
                if branch_length < -1e-10:
                    raise PipelineError(
                        f"Backbone depth decreased along scaffold path: parent={parent_depth:.12g} child={child_depth:.12g}"
                    )
                child.branch_length = max(float(min_branch_length), branch_length)
                desired_depths[child] = child_depth
            if not child.is_terminal():
                stack.append(child)
    return desired_depths


def find_placeholder_parent(tree, placeholder_name, min_branch_length):
    stack = [(tree.root, 0.0)]
    while stack:
        parent, parent_depth = stack.pop()
        for child in parent.clades:
            if child.is_terminal() and str(child.name) == placeholder_name:
                return parent, parent_depth
            child_depth = parent_depth + ensure_positive_branch_length(child.branch_length, min_branch_length)
            if not child.is_terminal():
                stack.append((child, child_depth))
    raise PipelineError(f"Could not find placeholder {placeholder_name}")


def compute_clade_height(clade, min_branch_length):
    def recurse(node, prefix):
        current = prefix + ensure_positive_branch_length(node.branch_length, min_branch_length)
        if node.is_terminal():
            return [current]
        values = []
        for child in node.clades:
            values.extend(recurse(child, current))
        return values

    depths = recurse(clade, 0.0)
    return min(depths), statistics.median(depths), max(depths)


def adjust_ingroup_terminal_branches_to_redline(tree, redline_height, outgroup_tip_name, min_branch_length):
    for tip in tree.get_terminals():
        if str(tip.name) == outgroup_tip_name:
            continue
        current_depth = float(tree.distance(tree.root, tip))
        delta = float(redline_height) - current_depth
        current_branch = ensure_positive_branch_length(tip.branch_length, min_branch_length)
        new_branch = current_branch + delta
        tip.branch_length = max(float(min_branch_length), new_branch)


def run(split_output_dir, backbone_tree, analysis_tree_dir, output_dir, min_branch_length=1e-8, redline_tolerance=1e-7):
    split_output_dir = Path(split_output_dir).resolve()
    backbone_tree = Path(backbone_tree).resolve()
    analysis_tree_dir = Path(analysis_tree_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    master_tree = read_newick_tree(split_output_dir / "intermediate" / "rooted.tree")
    fixed_backbone_tree = read_newick_tree(backbone_tree)
    target_records = load_target_summary(split_output_dir / "target_subtree_summary.tsv")
    paml_records = load_paml_summary(split_output_dir / "paml_subtree_summary.tsv")

    unique_outgroups = sorted({record.outgroup_tip for record in paml_records})
    if len(unique_outgroups) != 1:
        raise PipelineError(f"Expected exactly one unique outgroup, got {unique_outgroups}")
    outgroup_tip_name = unique_outgroups[0]

    validate_rooted_tree_with_outgroup(master_tree, outgroup_tip_name)
    validate_rooted_tree_with_outgroup(fixed_backbone_tree, outgroup_tip_name)

    backbone_tip_names = sorted(set(get_tip_names_from_tree(fixed_backbone_tree)) - {outgroup_tip_name})
    scaffold_tree, placeholder_to_target_id = build_scaffold_tree(
        master_tree=master_tree,
        target_records=target_records,
        backbone_tip_names=backbone_tip_names,
        outgroup_tip_name=outgroup_tip_name,
    )
    desired_depths = rewrite_scaffold_backbone_depths(
        scaffold_tree,
        fixed_backbone_tree=fixed_backbone_tree,
        backbone_tip_set=set(backbone_tip_names),
        min_branch_length=min_branch_length,
    )

    target_to_placeholder = {target_id: placeholder for placeholder, target_id in placeholder_to_target_id.items()}
    target_by_id = {record.target_subtree_id: record for record in target_records}
    paml_by_target_id = {record.target_subtree_id: record for record in paml_records}
    redline_height = summarize_root_to_tip(fixed_backbone_tree, excluded_tip_names={outgroup_tip_name})["max_depth"]

    scale_rows = []
    graft_rows = []

    for target_id, paml_record in sorted(paml_by_target_id.items()):
        analysis_tree_path = build_analysis_tree_path(
            merge_output_dir=output_dir,
            paml_subtree_id=paml_record.paml_subtree_id,
            analysis_tree_source="external",
            external_result_dir=analysis_tree_dir,
        )
        analysis_tree = read_newick_tree(analysis_tree_path)
        validate_tree_against_expected_tip_set(
            analysis_tree,
            paml_record.total_tip_names,
            paml_record.paml_subtree_id,
            outgroup_tip_name,
        )
        validate_rooted_tree_with_outgroup(analysis_tree, outgroup_tip_name)

        initial_scale = estimate_scale_factor_from_shared_backbone_tips(
            reference_tree=fixed_backbone_tree,
            analysis_tree=analysis_tree,
            shared_tip_names=paml_record.global_backbone_tip_names,
        )
        pre_scaled_tree = scale_tree_branch_lengths(
            analysis_tree,
            factor=float(initial_scale["scale_factor"]),
            min_branch_length=min_branch_length,
        )

        target_record = target_by_id[target_id]
        target_clade = extract_monophyletic_target_clade(
            pre_scaled_tree,
            target_record.target_nonbackbone_tip_names,
        )
        target_min, target_median, target_max = compute_clade_height(target_clade, min_branch_length)
        if target_max <= 0:
            raise PipelineError(f"Target clade height is non-positive for {target_id}: {target_max}")

        placeholder = target_to_placeholder[target_id]
        parent_clade, parent_depth = find_placeholder_parent(
            scaffold_tree,
            placeholder_name=placeholder,
            min_branch_length=min_branch_length,
        )
        attachment_depth = desired_depths.get(parent_clade, parent_depth)
        available_depth = redline_height - attachment_depth
        if available_depth <= 0:
            raise PipelineError(
                f"No available redline depth for {target_id}: attachment={attachment_depth:.12g} redline={redline_height:.12g}"
            )

        redline_factor = available_depth / target_max
        scaled_analysis_tree = scale_tree_branch_lengths(
            analysis_tree,
            factor=float(initial_scale["scale_factor"]) * float(redline_factor),
            min_branch_length=min_branch_length,
        )
        replacement = extract_monophyletic_target_clade(
            scaled_analysis_tree,
            target_record.target_nonbackbone_tip_names,
        )
        replacement = clone_clade(replacement)
        replacement.branch_length = ensure_positive_branch_length(replacement.branch_length, min_branch_length)
        repl_min, repl_median, repl_max = compute_clade_height(replacement, min_branch_length)
        if not replace_placeholder_with_clade(scaffold_tree, placeholder, replacement):
            raise PipelineError(f"Failed to replace placeholder {placeholder} for {target_id}")

        scale_rows.append(
            {
                "target_subtree_id": target_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "initial_scale_factor": f"{float(initial_scale['scale_factor']):.12g}",
                "redline_factor": f"{float(redline_factor):.12g}",
                "final_scale_factor": f"{float(initial_scale['scale_factor']) * float(redline_factor):.12g}",
                "shared_tip_n": str(initial_scale["shared_tip_n"]),
                "shared_pair_n": str(initial_scale["shared_pair_n"]),
                "median_log_ratio": f"{float(initial_scale['median_log_ratio']):.12g}",
                "residual_mad": f"{float(initial_scale['residual_mad']):.12g}",
                "status": "scaled_to_redline",
            }
        )
        graft_rows.append(
            {
                "target_subtree_id": target_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "attachment_placeholder": placeholder,
                "attachment_depth": f"{attachment_depth:.12g}",
                "available_depth": f"{available_depth:.12g}",
                "target_tip_n": str(target_record.target_nonbackbone_n_tips),
                "target_height_pre_min": f"{target_min:.12g}",
                "target_height_pre_median": f"{target_median:.12g}",
                "target_height_pre_max": f"{target_max:.12g}",
                "target_height_post_min": f"{repl_min:.12g}",
                "target_height_post_median": f"{repl_median:.12g}",
                "target_height_post_max": f"{repl_max:.12g}",
            }
        )

    validate_rooted_tree_with_outgroup(scaffold_tree, outgroup_tip_name)
    merged_tip_set = set(get_tip_names_from_tree(scaffold_tree))
    master_tip_set = set(get_tip_names_from_tree(master_tree))
    if merged_tip_set != master_tip_set:
        missing = sorted(master_tip_set - merged_tip_set)
        extra = sorted(merged_tip_set - master_tip_set)
        raise PipelineError(
            f"Final merged tip set mismatch. missing={missing[:10]} extra={extra[:10]}"
        )

    adjust_ingroup_terminal_branches_to_redline(
        scaffold_tree,
        redline_height=redline_height,
        outgroup_tip_name=outgroup_tip_name,
        min_branch_length=min_branch_length,
    )
    summary = summarize_root_to_tip(scaffold_tree, excluded_tip_names={outgroup_tip_name})
    if abs(summary["delta"]) > float(redline_tolerance):
        tip_depth_map = build_tip_depth_map(scaffold_tree)
        tip_depth_map = {tip: depth for tip, depth in tip_depth_map.items() if tip != outgroup_tip_name}
        min_depth = min(tip_depth_map.values())
        max_depth = max(tip_depth_map.values())
        min_tips = sorted([tip for tip, depth in tip_depth_map.items() if abs(depth - min_depth) <= 1e-12])
        max_tips = sorted([tip for tip, depth in tip_depth_map.items() if abs(depth - max_depth) <= 1e-12])
        write_tsv(
            [
                {"tip_name": tip, "root_to_tip_depth": f"{depth:.12g}"}
                for tip, depth in sorted(tip_depth_map.items(), key=lambda item: item[1])
            ],
            ["tip_name", "root_to_tip_depth"],
            output_dir / "redline_failed_tip_depths.tsv",
        )
        raise PipelineError(
            "Redline merge failed to produce an ultrametric ingroup tree: "
            f"delta={summary['delta']:.12g} min_tips={min_tips[:5]} max_tips={max_tips[:5]}"
        )

    summary_rows = [
        {
            "tree_label": "fixed_backbone_redline_merge",
            "n_tips": str(len(scaffold_tree.get_terminals())),
            "ingroup_min_depth": f"{summary['min_depth']:.12g}",
            "ingroup_median_depth": f"{summary['median_depth']:.12g}",
            "ingroup_max_depth": f"{summary['max_depth']:.12g}",
            "ingroup_delta": f"{summary['delta']:.12g}",
            "backbone_tree": backbone_tree.as_posix(),
            "analysis_tree_dir": analysis_tree_dir.as_posix(),
        }
    ]

    write_tree_file(scaffold_tree, output_dir / "merged_ml_tree_fixed_backbone_redline.nwk")
    write_tree_file(scaffold_tree, output_dir / "merged_tree_fixed_backbone_redline.nwk")
    write_tsv(
        scale_rows,
        [
            "target_subtree_id",
            "paml_subtree_id",
            "initial_scale_factor",
            "redline_factor",
            "final_scale_factor",
            "shared_tip_n",
            "shared_pair_n",
            "median_log_ratio",
            "residual_mad",
            "status",
        ],
        output_dir / "fixed_backbone_redline_scale_report.tsv",
    )
    write_tsv(
        graft_rows,
        [
            "target_subtree_id",
            "paml_subtree_id",
            "attachment_placeholder",
            "attachment_depth",
            "available_depth",
            "target_tip_n",
            "target_height_pre_min",
            "target_height_pre_median",
            "target_height_pre_max",
            "target_height_post_min",
            "target_height_post_median",
            "target_height_post_max",
        ],
        output_dir / "fixed_backbone_redline_graft_report.tsv",
    )
    write_tsv(
        summary_rows,
        [
            "tree_label",
            "n_tips",
            "ingroup_min_depth",
            "ingroup_median_depth",
            "ingroup_max_depth",
            "ingroup_delta",
            "backbone_tree",
            "analysis_tree_dir",
        ],
        output_dir / "fixed_backbone_redline_summary.tsv",
    )


def build_parser():
    parser = argparse.ArgumentParser(description="Redline-constrained graft merge on a fixed backbone tree.")
    parser.add_argument("--split-output-dir", required=True)
    parser.add_argument("--backbone-tree", required=True)
    parser.add_argument("--analysis-tree-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-branch-length", default=1e-8, type=float)
    parser.add_argument("--redline-tolerance", default=1e-7, type=float)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run(
            split_output_dir=args.split_output_dir,
            backbone_tree=args.backbone_tree,
            analysis_tree_dir=args.analysis_tree_dir,
            output_dir=args.output_dir,
            min_branch_length=args.min_branch_length,
            redline_tolerance=args.redline_tolerance,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Redline-constrained merge on a fixed ultrametric backbone."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Optional

from Bio.Phylo.BaseTree import Tree

from phylo_merge_common import (
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
from phylo_split_common import PipelineError, clone_clade, read_newick_tree


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
    edge_rows = []
    while stack:
        parent = stack.pop()
        parent_depth = desired_depths[parent]
        for child in parent.clades:
            if child.is_terminal() and str(child.name).startswith("TARGET_"):
                child.branch_length = None
                continue

            old_length = ensure_positive_branch_length(child.branch_length, min_branch_length)
            child_depth = compute_backbone_node_depth(
                child,
                fixed_backbone_tree=fixed_backbone_tree,
                backbone_tip_set=backbone_tip_set,
                min_branch_length=min_branch_length,
            )
            if child_depth is None:
                new_length = float(min_branch_length)
                child.branch_length = new_length
                desired_depths[child] = parent_depth + new_length
            else:
                new_length = child_depth - parent_depth
                if new_length < -1e-10:
                    raise PipelineError(
                        f"Backbone depth decreased along scaffold path: parent={parent_depth:.12g} child={child_depth:.12g}"
                    )
                new_length = max(float(min_branch_length), new_length)
                child.branch_length = new_length
                desired_depths[child] = child_depth
            edge_rows.append(
                {
                    "signature_key": "",
                    "old_branch_length": f"{old_length:.12g}",
                    "new_branch_length": f"{float(new_length):.12g}",
                    "descendant_backbone_n_tips": str(len(get_descendant_backbone_tips(child, backbone_tip_set))),
                }
            )
            if not child.is_terminal():
                stack.append(child)
    return desired_depths, edge_rows


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
        tip.branch_length = max(float(min_branch_length), current_branch + delta)


def compute_subtree_replacement(task):
    analysis_tree_path = Path(task["analysis_tree_path"])
    fixed_backbone_tree = read_newick_tree(Path(task["backbone_tree_path"]))
    analysis_tree = read_newick_tree(analysis_tree_path)
    validate_tree_against_expected_tip_set(
        analysis_tree,
        task["total_tip_names"],
        task["beast_subtree_id"],
        task["outgroup_tip_name"],
    )
    validate_rooted_tree_with_outgroup(analysis_tree, task["outgroup_tip_name"])

    initial_scale = estimate_scale_factor_from_shared_backbone_tips(
        reference_tree=fixed_backbone_tree,
        analysis_tree=analysis_tree,
        shared_tip_names=task["global_backbone_tip_names"],
    )
    pre_scaled_tree = scale_tree_branch_lengths(
        analysis_tree,
        factor=float(initial_scale["scale_factor"]),
        min_branch_length=task["min_branch_length"],
    )
    target_clade = extract_monophyletic_target_clade(
        pre_scaled_tree,
        task["target_nonbackbone_tip_names"],
    )
    target_min, target_median, target_max = compute_clade_height(target_clade, task["min_branch_length"])
    if target_max <= 0:
        raise PipelineError(f"Target clade height is non-positive for {task['target_subtree_id']}: {target_max}")

    available_depth = float(task["redline_height"]) - float(task["attachment_depth"])
    if available_depth <= 0:
        raise PipelineError(
            f"No available redline depth for {task['target_subtree_id']}: "
            f"attachment={float(task['attachment_depth']):.12g} redline={float(task['redline_height']):.12g}"
        )

    redline_factor = available_depth / target_max
    scaled_analysis_tree = scale_tree_branch_lengths(
        analysis_tree,
        factor=float(initial_scale["scale_factor"]) * float(redline_factor),
        min_branch_length=task["min_branch_length"],
    )
    replacement = extract_monophyletic_target_clade(
        scaled_analysis_tree,
        task["target_nonbackbone_tip_names"],
    )
    replacement = clone_clade(replacement)
    replacement.branch_length = ensure_positive_branch_length(replacement.branch_length, task["min_branch_length"])
    repl_min, repl_median, repl_max = compute_clade_height(replacement, task["min_branch_length"])
    del repl_min, repl_median

    replacement_clade_path = Path(task["replacement_clade_path"])
    write_tree_file(Tree(root=clone_clade(replacement), rooted=True), replacement_clade_path)

    final_scale_factor = float(initial_scale["scale_factor"]) * float(redline_factor)
    return {
        "target_subtree_id": task["target_subtree_id"],
        "beast_subtree_id": task["beast_subtree_id"],
        "placeholder": task["placeholder"],
        "replacement_clade_path": replacement_clade_path.as_posix(),
        "subtree_scale_row": {
            "target_subtree_id": task["target_subtree_id"],
            "beast_subtree_id": task["beast_subtree_id"],
            "scale_factor": f"{final_scale_factor:.12g}",
            "shared_paths_n": str(initial_scale["shared_pair_n"]),
            "median_log_ratio": f"{float(initial_scale['median_log_ratio']):.12g}",
            "residual_mad": f"{float(initial_scale['residual_mad']):.12g}",
            "status": "scaled_to_redline",
            "details": (
                f"initial_scale={float(initial_scale['scale_factor']):.12g}; "
                f"redline_factor={float(redline_factor):.12g}; "
                f"attachment_depth={float(task['attachment_depth']):.12g}; "
                f"available_depth={available_depth:.12g}; "
                f"target_height_post={repl_max:.12g}"
            ),
        },
        "graft_row": {
            "target_subtree_id": task["target_subtree_id"],
            "target_root_node_id": task["target_root_node_id"],
            "beast_subtree_id": task["beast_subtree_id"],
            "analysis_tree_file": analysis_tree_path.as_posix(),
            "target_nonbackbone_n_tips": str(task["target_nonbackbone_n_tips"]),
            "graft_status": "success",
            "attachment_placeholder": task["placeholder"],
            "topology_match": "true",
            "details": (
                f"attachment_depth={float(task['attachment_depth']):.12g}; "
                f"available_depth={available_depth:.12g}; "
                f"target_pre_max={target_max:.12g}; "
                f"target_post_max={repl_max:.12g}"
            ),
        },
        "edge_update_row": {
            "edge_role": "target_graft",
            "target_subtree_id": task["target_subtree_id"],
            "beast_subtree_id": task["beast_subtree_id"],
            "node_id": f"GRAFT::{task['target_subtree_id']}",
            "descendant_tip_hash": task["target_tip_hash"],
            "descendant_n_tips": str(task["target_nonbackbone_n_tips"]),
            "old_branch_length": f"{float(task['min_branch_length']):.12g}",
            "new_branch_length": f"{ensure_positive_branch_length(replacement.branch_length, task['min_branch_length']):.12g}",
            "aggregation_source": task["placeholder"],
            "changed": "true",
        },
    }


def run(
    split_output_dir,
    merge_output_dir,
    analysis_tree_source,
    external_result_dir,
    backbone_tree,
    min_branch_length=1e-8,
    redline_tolerance=1e-6,
    parallel_jobs=8,
):
    if analysis_tree_source not in {"external", "simulated"}:
        raise PipelineError(f"Unsupported analysis_tree_source: {analysis_tree_source}")

    split_output_dir = Path(split_output_dir).resolve()
    merge_output_dir = Path(merge_output_dir).resolve()
    merge_output_dir.mkdir(parents=True, exist_ok=True)
    backbone_tree = Path(backbone_tree).resolve()
    external_result_dir = None if external_result_dir in (None, "", "null") else Path(external_result_dir).resolve()

    master_tree = read_newick_tree(split_output_dir / "intermediate" / "rooted.tree")
    fixed_backbone_tree = read_newick_tree(backbone_tree)
    target_records = load_target_summary(split_output_dir / "target_subtree_summary.tsv")
    paml_records = load_paml_summary(split_output_dir / "beast_subtree_summary.tsv")

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
    write_tree_file(scaffold_tree, merge_output_dir / "assembly_scaffold.nwk")
    desired_depths, backbone_edge_rows = rewrite_scaffold_backbone_depths(
        scaffold_tree,
        fixed_backbone_tree=fixed_backbone_tree,
        backbone_tip_set=set(backbone_tip_names),
        min_branch_length=min_branch_length,
    )

    target_to_placeholder = {target_id: placeholder for placeholder, target_id in placeholder_to_target_id.items()}
    target_by_id = {record.target_subtree_id: record for record in target_records}
    paml_by_target_id = {record.target_subtree_id: record for record in paml_records}
    redline_height = summarize_root_to_tip(fixed_backbone_tree, excluded_tip_names={outgroup_tip_name})["max_depth"]
    total_subtrees = len(paml_by_target_id)
    print(
        f"[INFO] Redline merge started: total_subtrees={total_subtrees} "
        f"backbone_tips={len(backbone_tip_names)} parallel_jobs={int(parallel_jobs)} "
        f"redline_tolerance={float(redline_tolerance):.12g}",
        flush=True,
    )

    subtree_scale_rows = []
    graft_rows = []
    edge_update_rows = []
    ordered_pairs = sorted(paml_by_target_id.items())
    tasks = []
    replacement_clade_dir = merge_output_dir / "intermediate" / "replacement_clades"
    for target_id, paml_record in ordered_pairs:
        target_record = target_by_id[target_id]
        placeholder = target_to_placeholder[target_id]
        parent_clade, parent_depth = find_placeholder_parent(
            scaffold_tree,
            placeholder_name=placeholder,
            min_branch_length=min_branch_length,
        )
        attachment_depth = desired_depths.get(parent_clade, parent_depth)
        tasks.append(
            {
                "target_subtree_id": target_id,
                "beast_subtree_id": paml_record.paml_subtree_id,
                "target_root_node_id": target_record.target_root_node_id,
                "target_tip_hash": target_record.target_tip_hash,
                "target_nonbackbone_n_tips": target_record.target_nonbackbone_n_tips,
                "target_nonbackbone_tip_names": list(target_record.target_nonbackbone_tip_names),
                "total_tip_names": list(paml_record.total_tip_names),
                "global_backbone_tip_names": list(paml_record.global_backbone_tip_names),
                "outgroup_tip_name": outgroup_tip_name,
                "analysis_tree_path": build_analysis_tree_path(
                    merge_output_dir=merge_output_dir,
                    paml_subtree_id=paml_record.paml_subtree_id,
                    analysis_tree_source=analysis_tree_source,
                    external_result_dir=external_result_dir,
                ).as_posix(),
                "backbone_tree_path": backbone_tree.as_posix(),
                "attachment_depth": float(attachment_depth),
                "redline_height": float(redline_height),
                "min_branch_length": float(min_branch_length),
                "placeholder": placeholder,
                "replacement_clade_path": (replacement_clade_dir / f"{target_id}.nwk").as_posix(),
            }
        )

    if int(parallel_jobs) <= 1:
        result_iter = map(compute_subtree_replacement, tasks)
        for index, result in enumerate(result_iter, start=1):
            if index == 1 or index % 5 == 0 or index == total_subtrees:
                print(
                    f"[INFO] Redline merge progress: subtree {index}/{total_subtrees} "
                    f"{result['target_subtree_id']} ({result['beast_subtree_id']})",
                    flush=True,
                )
            replacement_path = Path(result["replacement_clade_path"])
            replacement = read_newick_tree(replacement_path).root
            replacement_path.unlink(missing_ok=True)
            replacement.branch_length = ensure_positive_branch_length(replacement.branch_length, min_branch_length)
            if not replace_placeholder_with_clade(scaffold_tree, result["placeholder"], replacement):
                raise PipelineError(
                    f"Failed to replace placeholder {result['placeholder']} for {result['target_subtree_id']}"
                )
            subtree_scale_rows.append(result["subtree_scale_row"])
            graft_rows.append(result["graft_row"])
            edge_update_rows.append(result["edge_update_row"])
    else:
        with ProcessPoolExecutor(max_workers=int(parallel_jobs)) as executor:
            result_iter = executor.map(compute_subtree_replacement, tasks)
            for index, result in enumerate(result_iter, start=1):
                if index == 1 or index % 5 == 0 or index == total_subtrees:
                    print(
                        f"[INFO] Redline merge progress: subtree {index}/{total_subtrees} "
                        f"{result['target_subtree_id']} ({result['beast_subtree_id']})",
                        flush=True,
                    )
                replacement_path = Path(result["replacement_clade_path"])
                replacement = read_newick_tree(replacement_path).root
                replacement_path.unlink(missing_ok=True)
                replacement.branch_length = ensure_positive_branch_length(replacement.branch_length, min_branch_length)
                if not replace_placeholder_with_clade(scaffold_tree, result["placeholder"], replacement):
                    raise PipelineError(
                        f"Failed to replace placeholder {result['placeholder']} for {result['target_subtree_id']}"
                    )
                subtree_scale_rows.append(result["subtree_scale_row"])
                graft_rows.append(result["graft_row"])
                edge_update_rows.append(result["edge_update_row"])

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
            merge_output_dir / "redline_failed_tip_depths.tsv",
        )
        raise PipelineError(
            "Redline merge failed to produce an ultrametric ingroup tree: "
            f"delta={summary['delta']:.12g} min_tips={min_tips[:5]} max_tips={max_tips[:5]}"
        )

    write_tree_file(scaffold_tree, merge_output_dir / "merged_ml_tree.nwk")
    write_tree_file(scaffold_tree, merge_output_dir / "merged_tree.nwk")
    write_tsv(backbone_edge_rows, ["signature_key", "old_branch_length", "new_branch_length", "descendant_backbone_n_tips"], merge_output_dir / "backbone_edge_estimates.tsv")
    write_tsv(graft_rows, ["target_subtree_id", "target_root_node_id", "beast_subtree_id", "analysis_tree_file", "target_nonbackbone_n_tips", "graft_status", "attachment_placeholder", "topology_match", "details"], merge_output_dir / "graft_report.tsv")
    write_tsv(edge_update_rows, ["edge_role", "target_subtree_id", "beast_subtree_id", "node_id", "descendant_tip_hash", "descendant_n_tips", "old_branch_length", "new_branch_length", "aggregation_source", "changed"], merge_output_dir / "edge_update_report.tsv")
    write_tsv(subtree_scale_rows, ["target_subtree_id", "beast_subtree_id", "scale_factor", "shared_paths_n", "median_log_ratio", "residual_mad", "status", "details"], merge_output_dir / "subtree_scale_report.tsv")
    print(
        f"[INFO] Redline merge finished: merged_tips={len(scaffold_tree.get_terminals())} "
        f"ingroup_delta={float(summary['delta']):.12g}",
        flush=True,
    )
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Redline-constrained merge on a fixed ultrametric backbone.")
    parser.add_argument("--split-output-dir", required=True)
    parser.add_argument("--merge-output-dir", required=True)
    parser.add_argument("--analysis-tree-source", required=True)
    parser.add_argument("--external-result-dir", default=None)
    parser.add_argument("--backbone-tree", required=True)
    parser.add_argument("--min-branch-length", default=1e-8, type=float)
    parser.add_argument("--redline-tolerance", default=1e-6, type=float)
    parser.add_argument("--parallel-jobs", default=8, type=int)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            split_output_dir=args.split_output_dir,
            merge_output_dir=args.merge_output_dir,
            analysis_tree_source=args.analysis_tree_source,
            external_result_dir=args.external_result_dir,
            backbone_tree=args.backbone_tree,
            min_branch_length=args.min_branch_length,
            redline_tolerance=args.redline_tolerance,
            parallel_jobs=args.parallel_jobs,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

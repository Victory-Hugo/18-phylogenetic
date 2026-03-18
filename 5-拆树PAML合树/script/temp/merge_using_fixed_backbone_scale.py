#!/usr/bin/env python3
"""Simplified merge experiment with a fixed ultrametric backbone.

Rule set:
1. Keep the backbone tree unchanged.
2. For each subtree, estimate one global scale factor from shared backbone tips.
3. Extract the exclusive target clade from the scaled subtree.
4. Replace the corresponding placeholder node in the scaffold.
5. Require the final merged tree to contain exactly the same tip set as the
   original rooted master tree.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = PROJECT_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from phylo_merge_common import (  # noqa: E402
    build_analysis_tree_path,
    build_scaffold_tree,
    ensure_positive_branch_length,
    extract_monophyletic_target_clade,
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


def summarize_root_to_tip(tree):
    depths = [float(tree.distance(tree.root, tip)) for tip in tree.get_terminals()]
    return {
        "n_tips": len(depths),
        "min_depth": min(depths),
        "median_depth": statistics.median(depths),
        "max_depth": max(depths),
        "delta": max(depths) - min(depths),
    }


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


def run(split_output_dir, backbone_tree, analysis_tree_dir, output_dir, min_branch_length=1e-8):
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
    target_to_placeholder = {target_id: placeholder for placeholder, target_id in placeholder_to_target_id.items()}
    target_by_id = {record.target_subtree_id: record for record in target_records}
    paml_by_target_id = {record.target_subtree_id: record for record in paml_records}

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

        scale_summary = estimate_scale_factor_from_shared_backbone_tips(
            reference_tree=fixed_backbone_tree,
            analysis_tree=analysis_tree,
            shared_tip_names=paml_record.global_backbone_tip_names,
        )
        scaled_tree = scale_tree_branch_lengths(
            analysis_tree,
            factor=float(scale_summary["scale_factor"]),
            min_branch_length=min_branch_length,
        )

        target_record = target_by_id[target_id]
        target_clade = extract_monophyletic_target_clade(
            scaled_tree,
            target_record.target_nonbackbone_tip_names,
        )
        replacement = clone_clade(target_clade)
        replacement.branch_length = ensure_positive_branch_length(replacement.branch_length, min_branch_length)

        placeholder = target_to_placeholder[target_id]
        if not replace_placeholder_with_clade(scaffold_tree, placeholder, replacement):
            raise PipelineError(f"Failed to replace placeholder {placeholder} for {target_id}")

        scale_rows.append(
            {
                "target_subtree_id": target_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "scale_factor": f"{float(scale_summary['scale_factor']):.12g}",
                "shared_tip_n": str(scale_summary["shared_tip_n"]),
                "shared_pair_n": str(scale_summary["shared_pair_n"]),
                "median_log_ratio": f"{float(scale_summary['median_log_ratio']):.12g}",
                "residual_mad": f"{float(scale_summary['residual_mad']):.12g}",
                "status": str(scale_summary["status"]),
            }
        )
        graft_rows.append(
            {
                "target_subtree_id": target_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "attachment_placeholder": placeholder,
                "target_tip_n": str(target_record.target_nonbackbone_n_tips),
                "replacement_root_branch_length": f"{replacement.branch_length:.12g}",
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

    summary = summarize_root_to_tip(scaffold_tree)
    summary_rows = [
        {
            "tree_label": "fixed_backbone_simple_graft",
            "n_tips": str(summary["n_tips"]),
            "min_depth": f"{summary['min_depth']:.12g}",
            "median_depth": f"{summary['median_depth']:.12g}",
            "max_depth": f"{summary['max_depth']:.12g}",
            "delta": f"{summary['delta']:.12g}",
            "backbone_tree": backbone_tree.as_posix(),
            "analysis_tree_dir": analysis_tree_dir.as_posix(),
        }
    ]

    write_tree_file(scaffold_tree, output_dir / "merged_ml_tree_fixed_backbone_simple_graft.nwk")
    write_tree_file(scaffold_tree, output_dir / "merged_tree_fixed_backbone_simple_graft.nwk")
    write_tsv(
        scale_rows,
        ["target_subtree_id", "paml_subtree_id", "scale_factor", "shared_tip_n", "shared_pair_n", "median_log_ratio", "residual_mad", "status"],
        output_dir / "fixed_backbone_simple_graft_scale_report.tsv",
    )
    write_tsv(
        graft_rows,
        ["target_subtree_id", "paml_subtree_id", "attachment_placeholder", "target_tip_n", "replacement_root_branch_length"],
        output_dir / "fixed_backbone_simple_graft_report.tsv",
    )
    write_tsv(
        summary_rows,
        ["tree_label", "n_tips", "min_depth", "median_depth", "max_depth", "delta", "backbone_tree", "analysis_tree_dir"],
        output_dir / "fixed_backbone_simple_graft_summary.tsv",
    )


def build_parser():
    parser = argparse.ArgumentParser(description="Simplified graft merge on a fixed backbone tree.")
    parser.add_argument("--split-output-dir", required=True)
    parser.add_argument("--backbone-tree", required=True)
    parser.add_argument("--analysis-tree-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-branch-length", default=1e-8, type=float)
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
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

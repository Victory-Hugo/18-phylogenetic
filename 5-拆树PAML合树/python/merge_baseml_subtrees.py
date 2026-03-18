#!/usr/bin/env python3
"""Merge PAML subtree results by backbone aggregation and target grafting."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from phylo_merge_common import (
    BACKBONE_EDGE_ESTIMATE_COLUMNS,
    EDGE_UPDATE_COLUMNS,
    GRAFT_REPORT_COLUMNS,
    SUBTREE_SCALE_REPORT_COLUMNS,
    build_analysis_tree_path,
    build_backbone_signature_key,
    build_edge_update_row,
    build_scale_reference_tip_names,
    build_scaffold_tree,
    build_tip_distance_context,
    compute_branch_by_signature,
    estimate_scale_factor_from_reference_tips,
    ensure_positive_branch_length,
    extract_monophyletic_target_clade,
    get_tip_names_from_clade,
    get_tip_names_from_tree,
    load_backbone_summary,
    load_paml_summary,
    load_target_summary,
    read_newick_tree,
    replace_placeholder_with_clade,
    scale_tree_branch_lengths,
    validate_rooted_tree_with_outgroup,
    validate_rooted_binary_tree,
    validate_tree_against_expected_tip_set,
    weighted_geometric_mean,
    write_rows,
    write_tree_file,
)
from phylo_split_common import PipelineError, assign_node_ids, clone_clade, compute_tip_hash, setup_logger


def _resolve_external_result_dir(external_result_dir):
    if external_result_dir in (None, "", "null"):
        return None
    return Path(external_result_dir).resolve()


def run(
    split_output_dir,
    merge_output_dir,
    analysis_tree_source,
    external_result_dir,
    merge_mode="backbone_graft",
    subtree_scale_method="median_log_path_ratio",
    backbone_edge_aggregation="weighted_geometric_mean",
    min_branch_length=1e-8,
    log_level="INFO",
    outgroup_tip_name=None,
):
    split_output_dir = Path(split_output_dir).resolve()
    merge_output_dir = Path(merge_output_dir).resolve()
    external_result_dir = _resolve_external_result_dir(external_result_dir)
    min_branch_length = float(min_branch_length)

    if merge_mode != "backbone_graft":
        raise PipelineError(f"Unsupported merge.mode: {merge_mode}")
    if subtree_scale_method not in {"identity_graft", "median_log_path_ratio"}:
        raise PipelineError(f"Unsupported merge.subtree_scale_method: {subtree_scale_method}")
    if backbone_edge_aggregation != "weighted_geometric_mean":
        raise PipelineError(f"Unsupported merge.backbone_edge_aggregation: {backbone_edge_aggregation}")

    merge_output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("merge_baseml_subtrees", merge_output_dir / "merge.log", log_level)
    logger.info("Starting backbone-graft merge.")

    master_tree_path = split_output_dir / "intermediate" / "rooted.tree"
    if not master_tree_path.exists():
        raise PipelineError(f"Required rooted master tree not found: {master_tree_path}")
    master_tree = read_newick_tree(master_tree_path)

    backbone_records = load_backbone_summary(split_output_dir / "backbone_summary.tsv")
    target_records = load_target_summary(split_output_dir / "target_subtree_summary.tsv")
    paml_records = load_paml_summary(split_output_dir / "paml_subtree_summary.tsv")
    if outgroup_tip_name in (None, "", "null"):
        unique_outgroups = sorted({record.outgroup_tip for record in paml_records})
        if len(unique_outgroups) != 1:
            raise PipelineError(
                f"Could not infer a unique outgroup tip from {split_output_dir / 'paml_subtree_summary.tsv'}: {unique_outgroups}"
            )
        outgroup_tip_name = unique_outgroups[0]
    validate_rooted_binary_tree(master_tree, outgroup_tip_name)

    backbone_tips = [record.backbone_tip_name for record in sorted(backbone_records, key=lambda item: item.selection_rank)]
    backbone_tip_set = set(backbone_tips)
    master_distance_context = build_tip_distance_context(master_tree)
    target_by_id = {record.target_subtree_id: record for record in target_records}
    paml_by_target_id = {record.target_subtree_id: record for record in paml_records}
    if set(target_by_id) != set(paml_by_target_id):
        raise PipelineError("target_subtree_summary.tsv and paml_subtree_summary.tsv do not cover the same target ids.")

    scaffold_tree, placeholder_to_target_id = build_scaffold_tree(master_tree, target_records, backbone_tips, outgroup_tip_name)
    target_to_placeholder = {target_id: placeholder for placeholder, target_id in placeholder_to_target_id.items()}
    write_tree_file(scaffold_tree, merge_output_dir / "assembly_scaffold.nwk")

    backbone_estimates: dict[str, list[tuple[float, int, str]]] = {}
    graft_report_rows = []
    backbone_edge_rows = []
    edge_update_rows = []
    subtree_scale_rows = []
    extracted_grafts = {}

    for target_id, paml_record in sorted(paml_by_target_id.items()):
        analysis_tree_path = build_analysis_tree_path(
            merge_output_dir=merge_output_dir,
            paml_subtree_id=paml_record.paml_subtree_id,
            analysis_tree_source=analysis_tree_source,
            external_result_dir=external_result_dir,
        )
        if not analysis_tree_path.exists():
            raise PipelineError(f"Analysis tree file not found: {analysis_tree_path}")
        analysis_tree = read_newick_tree(analysis_tree_path)
        validate_tree_against_expected_tip_set(
            analysis_tree,
            paml_record.total_tip_names,
            paml_record.paml_subtree_id,
            outgroup_tip_name,
        )
        validate_rooted_binary_tree(analysis_tree, outgroup_tip_name)

        scale_summary = {
            "scale_factor": 1.0,
            "shared_paths_n": 0,
            "median_log_ratio": 0.0,
            "residual_mad": 0.0,
            "status": "identity_graft",
            "details": "Current merge stage preserves direct subtree branch lengths.",
        }
        scaled_analysis_tree = analysis_tree
        if subtree_scale_method == "median_log_path_ratio":
            reference_tip_names = build_scale_reference_tip_names(paml_record)
            scale_summary = estimate_scale_factor_from_reference_tips(
                reference_tree=master_tree,
                analysis_tree=analysis_tree,
                reference_tip_names=reference_tip_names,
                reference_context=master_distance_context,
            )
            if scale_summary["status"] == "scaled_by_reference_paths":
                scaled_analysis_tree = scale_tree_branch_lengths(
                    analysis_tree,
                    factor=float(scale_summary["scale_factor"]),
                    min_branch_length=min_branch_length,
                )
            else:
                scale_summary["details"] = (
                    f"{scale_summary['details']} Falling back to identity graft for {paml_record.paml_subtree_id}."
                )
        logger.info(
            "%s scale_factor=%.12g shared_paths=%s status=%s",
            paml_record.paml_subtree_id,
            float(scale_summary["scale_factor"]),
            scale_summary["shared_paths_n"],
            scale_summary["status"],
        )

        signature_map = compute_branch_by_signature(
            scaled_analysis_tree,
            backbone_tip_set=backbone_tip_set,
            outgroup_tip_name=outgroup_tip_name,
            min_branch_length=min_branch_length,
        )
        for signature_key, branch_length in signature_map.items():
            backbone_estimates.setdefault(signature_key, []).append(
                (branch_length, paml_record.target_nonbackbone_n_tips, paml_record.paml_subtree_id)
            )

        target_record = target_by_id[target_id]
        try:
            target_clade = extract_monophyletic_target_clade(scaled_analysis_tree, target_record.target_nonbackbone_tip_names)
            extracted_grafts[target_id] = clone_clade(target_clade)
            graft_report_rows.append(
                {
                    "target_subtree_id": target_record.target_subtree_id,
                    "target_root_node_id": target_record.target_root_node_id,
                    "paml_subtree_id": paml_record.paml_subtree_id,
                    "analysis_tree_file": analysis_tree_path.as_posix(),
                    "target_nonbackbone_n_tips": str(target_record.target_nonbackbone_n_tips),
                    "graft_status": "success",
                    "attachment_placeholder": target_to_placeholder[target_id],
                    "topology_match": "true",
                    "details": "strict monophyletic target clade extracted",
                }
            )
        except Exception as exc:
            graft_report_rows.append(
                {
                    "target_subtree_id": target_record.target_subtree_id,
                    "target_root_node_id": target_record.target_root_node_id,
                    "paml_subtree_id": paml_record.paml_subtree_id,
                    "analysis_tree_file": analysis_tree_path.as_posix(),
                    "target_nonbackbone_n_tips": str(target_record.target_nonbackbone_n_tips),
                    "graft_status": "failed",
                    "attachment_placeholder": target_to_placeholder[target_id],
                    "topology_match": "false",
                    "details": str(exc),
                }
            )
            raise
        subtree_scale_rows.append(
            {
                "target_subtree_id": target_record.target_subtree_id,
                "paml_subtree_id": paml_record.paml_subtree_id,
                "scale_factor": f"{float(scale_summary['scale_factor']):.12g}",
                "shared_paths_n": str(scale_summary["shared_paths_n"]),
                "median_log_ratio": (
                    f"{float(scale_summary['median_log_ratio']):.12g}"
                    if scale_summary["shared_paths_n"]
                    else ""
                ),
                "residual_mad": (
                    f"{float(scale_summary['residual_mad']):.12g}"
                    if scale_summary["shared_paths_n"]
                    else ""
                ),
                "status": str(scale_summary["status"]),
                "details": str(scale_summary["details"]),
            }
        )

    aggregated_backbone = {}
    for signature_key, values in sorted(backbone_estimates.items()):
        aggregated_backbone[signature_key] = weighted_geometric_mean([(value, weight) for value, weight, _ in values])
        support_ids = [paml_id for _, _, paml_id in values]
        descendant_n_tips = 0 if signature_key == "__OUTGROUP__" else len(values)
        backbone_edge_rows.append(
            {
                "signature_key": signature_key,
                "descendant_backbone_n_tips": str(descendant_n_tips),
                "supporting_paml_subtree_ids": json.dumps(support_ids, ensure_ascii=False),
                "n_estimates": str(len(values)),
                "aggregated_branch_length": f"{aggregated_backbone[signature_key]:.12g}",
                "details": "weighted_geometric_mean",
            }
        )

    scaffold_node_id_map, _, _ = assign_node_ids(scaffold_tree)
    for clade, node_id in scaffold_node_id_map.items():
        if clade is scaffold_tree.root:
            continue
        signature_key = build_backbone_signature_key(clade, backbone_tip_set, outgroup_tip_name)
        if signature_key is None or signature_key not in aggregated_backbone:
            continue
        old_branch_length = ensure_positive_branch_length(clade.branch_length, min_branch_length)
        new_branch_length = aggregated_backbone[signature_key]
        clade.branch_length = new_branch_length
        descendant_tips = [tip_name for tip_name in get_tip_names_from_clade(clade) if not str(tip_name).startswith("TARGET_")]
        edge_update_rows.append(
            build_edge_update_row(
                edge_role="backbone",
                target_subtree_id="",
                paml_subtree_id="",
                node_id=node_id,
                descendant_tip_hash=compute_tip_hash(descendant_tips or [signature_key]),
                descendant_n_tips=len(descendant_tips),
                old_branch_length=old_branch_length,
                new_branch_length=new_branch_length,
                aggregation_source=signature_key,
            )
        )

    for target_record in target_records:
        placeholder = target_to_placeholder[target_record.target_subtree_id]
        replacement = clone_clade(extracted_grafts[target_record.target_subtree_id])
        replacement.branch_length = ensure_positive_branch_length(replacement.branch_length, min_branch_length)
        if not replace_placeholder_with_clade(scaffold_tree, placeholder, replacement):
            raise PipelineError(f"Failed to replace placeholder {placeholder} in scaffold tree.")
        edge_update_rows.append(
            build_edge_update_row(
                edge_role="target_graft",
                target_subtree_id=target_record.target_subtree_id,
                paml_subtree_id=paml_by_target_id[target_record.target_subtree_id].paml_subtree_id,
                node_id=f"GRAFT::{target_record.target_subtree_id}",
                descendant_tip_hash=compute_tip_hash(target_record.target_nonbackbone_tip_names),
                descendant_n_tips=target_record.target_nonbackbone_n_tips,
                old_branch_length=min_branch_length,
                new_branch_length=ensure_positive_branch_length(replacement.branch_length, min_branch_length),
                aggregation_source=placeholder,
            )
        )

    validate_rooted_tree_with_outgroup(scaffold_tree, outgroup_tip_name)
    master_tip_set = set(get_tip_names_from_tree(master_tree))
    merged_tip_set = set(get_tip_names_from_tree(scaffold_tree))
    if merged_tip_set != master_tip_set:
        missing = sorted(master_tip_set - merged_tip_set)
        extra = sorted(merged_tip_set - master_tip_set)
        raise PipelineError(f"Merged tip set mismatch. missing={missing[:10]} extra={extra[:10]}")

    write_rows(backbone_edge_rows, BACKBONE_EDGE_ESTIMATE_COLUMNS, merge_output_dir / "backbone_edge_estimates.tsv")
    write_rows(graft_report_rows, GRAFT_REPORT_COLUMNS, merge_output_dir / "graft_report.tsv")
    write_rows(edge_update_rows, EDGE_UPDATE_COLUMNS, merge_output_dir / "edge_update_report.tsv")
    write_rows(subtree_scale_rows, SUBTREE_SCALE_REPORT_COLUMNS, merge_output_dir / "subtree_scale_report.tsv")
    write_tree_file(scaffold_tree, merge_output_dir / "merged_ml_tree.nwk")
    write_tree_file(scaffold_tree, merge_output_dir / "merged_tree.nwk")
    logger.info("Merged %d graft targets onto scaffold.", len(target_records))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Merge PAML subtree results by backbone grafting.")
    parser.add_argument("--split-output-dir", required=True, help="Output directory from the split pipeline.")
    parser.add_argument("--merge-output-dir", required=True, help="Merge output directory.")
    parser.add_argument("--analysis-tree-source", required=True, help="external or simulated")
    parser.add_argument("--external-result-dir", default=None, help="Directory containing parsed analysis trees.")
    parser.add_argument("--merge-mode", default="backbone_graft", help="Merge mode.")
    parser.add_argument("--subtree-scale-method", default="median_log_path_ratio", help="Subtree scale unification strategy.")
    parser.add_argument("--backbone-edge-aggregation", default="weighted_geometric_mean", help="Backbone edge aggregation strategy.")
    parser.add_argument("--min-branch-length", default=1e-8, type=float, help="Minimum positive branch length.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    parser.add_argument("--outgroup-tip-name", default=None, help="Singleton outgroup tip retained at the root.")
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
            merge_mode=args.merge_mode,
            subtree_scale_method=args.subtree_scale_method,
            backbone_edge_aggregation=args.backbone_edge_aggregation,
            min_branch_length=args.min_branch_length,
            log_level=args.log_level,
            outgroup_tip_name=args.outgroup_tip_name,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

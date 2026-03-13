#!/usr/bin/env python3
"""Merge simulated or external baseml subtree results back onto the master tree."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import (
    EDGE_UPDATE_COLUMNS,
    MERGE_REPORT_COLUMNS,
    BasemlSummaryRecord,
    CoreSummaryRecord,
    build_analysis_tree_path,
    build_clade_signature_maps,
    build_descendant_core_ids,
    build_edge_update_row,
    compute_tip_hash,
    compute_multiplier,
    edge_role_for_core,
    edge_role_for_scaffold,
    ensure_positive_branch_length,
    find_core_mrca,
    get_outgroup_child,
    load_baseml_summary,
    load_core_summary,
    read_newick_tree,
    validate_rooted_binary_tree,
    validate_tree_against_expected_tip_set,
    weighted_geometric_mean,
    write_rows,
    write_tree_file,
)
from phylo_split_common import GLOBAL_OUTGROUP_TIP, PipelineError, assign_node_ids, setup_logger


def _resolve_external_result_dir(external_result_dir):
    if external_result_dir in (None, "", "null"):
        return None
    return Path(external_result_dir).resolve()


def _build_core_lookup(core_records: list[CoreSummaryRecord]) -> dict[str, CoreSummaryRecord]:
    lookup = {}
    for record in core_records:
        lookup[record.core_subtree_id] = record
    return lookup


def _load_analysis_tree(record: BasemlSummaryRecord, merge_output_dir: Path, analysis_tree_source: str, external_result_dir: Path | None):
    analysis_tree_path = build_analysis_tree_path(
        merge_output_dir=merge_output_dir,
        baseml_subtree_id=record.baseml_subtree_id,
        analysis_tree_source=analysis_tree_source,
        external_result_dir=external_result_dir,
    )
    if not analysis_tree_path.exists():
        raise PipelineError(f"Analysis tree file not found: {analysis_tree_path}")
    analysis_tree = read_newick_tree(analysis_tree_path)
    validate_tree_against_expected_tip_set(analysis_tree, record.total_tip_names, record.baseml_subtree_id)
    validate_rooted_binary_tree(analysis_tree, GLOBAL_OUTGROUP_TIP)
    return analysis_tree_path, analysis_tree


def run(
    input_tree,
    outgroup_tip_file,
    split_output_dir,
    merge_output_dir,
    analysis_tree_source,
    external_result_dir,
    conda_env,
    gotree_bin,
    threads,
    preserve_master_topology=True,
    strict_core_topology_match=True,
    scaffold_scale_aggregation="weighted_geometric_mean",
    scaffold_scale_clamp=None,
    min_branch_length=1e-8,
    log_level="INFO",
):
    del input_tree, outgroup_tip_file, conda_env, gotree_bin, threads

    split_output_dir = Path(split_output_dir).resolve()
    merge_output_dir = Path(merge_output_dir).resolve()
    external_result_dir = _resolve_external_result_dir(external_result_dir)
    min_branch_length = float(min_branch_length)

    if not preserve_master_topology:
        raise PipelineError("This workflow requires preserve_master_topology=true.")
    if not strict_core_topology_match:
        raise PipelineError("This workflow requires strict_core_topology_match=true.")
    if scaffold_scale_aggregation != "weighted_geometric_mean":
        raise PipelineError(f"Unsupported scaffold_scale_aggregation: {scaffold_scale_aggregation}")
    if scaffold_scale_clamp not in (None, "", "null"):
        raise PipelineError("This workflow does not support scaffold_scale_clamp.")

    merge_output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("merge_baseml_subtrees", merge_output_dir / "merge.log", log_level)
    logger.info("Starting baseml subtree merge.")

    master_tree_path = split_output_dir / "intermediate" / "rooted.tree"
    if not master_tree_path.exists():
        raise PipelineError(f"Required rooted master tree not found: {master_tree_path}")

    master_tree = read_newick_tree(master_tree_path)
    validate_rooted_binary_tree(master_tree, GLOBAL_OUTGROUP_TIP)

    core_records = load_core_summary(split_output_dir / "core_subtree_summary.tsv")
    baseml_records = load_baseml_summary(split_output_dir / "baseml_subtree_summary.tsv")
    core_lookup = _build_core_lookup(core_records)

    if len(core_records) != len(baseml_records):
        raise PipelineError("core_subtree_summary.tsv and baseml_subtree_summary.tsv row counts do not match.")

    node_id_map, reverse_node_id_map, _ = assign_node_ids(master_tree)
    original_branch_by_node_id = {
        node_id: clade.branch_length
        for clade, node_id in node_id_map.items()
        if clade is not master_tree.root
    }

    updated_node_ids = set()
    core_root_multipliers = {}
    merge_report_rows = []
    edge_update_rows = []

    for baseml_record in baseml_records:
        if baseml_record.core_subtree_id not in core_lookup:
            raise PipelineError(f"Missing core summary for {baseml_record.core_subtree_id}")
        core_record = core_lookup[baseml_record.core_subtree_id]
        analysis_tree_path, analysis_tree = _load_analysis_tree(
            record=baseml_record,
            merge_output_dir=merge_output_dir,
            analysis_tree_source=analysis_tree_source,
            external_result_dir=external_result_dir,
        )

        master_core_clade = reverse_node_id_map.get(core_record.core_root_node_id)
        if master_core_clade is None:
            raise PipelineError(f"Master core root node id not found: {core_record.core_root_node_id}")

        analysis_core_mrca = find_core_mrca(analysis_tree, baseml_record.core_tip_names)
        master_signature_map, master_signature_counts = build_clade_signature_maps(master_core_clade)
        analysis_signature_map, _ = build_clade_signature_maps(analysis_core_mrca)

        if set(master_signature_map) != set(analysis_signature_map):
            raise PipelineError(f"Core topology mismatch for {core_record.core_subtree_id}")

        matched_internal_edges = 0
        matched_tip_edges = 0
        replaced_core_root_edge = False

        for signature, master_clade in master_signature_map.items():
            analysis_clade = analysis_signature_map[signature]
            node_id = node_id_map[master_clade]
            old_branch_length = ensure_positive_branch_length(original_branch_by_node_id[node_id], min_branch_length)
            new_branch_length = ensure_positive_branch_length(analysis_clade.branch_length, min_branch_length)
            scale_multiplier = new_branch_length / old_branch_length

            master_clade.branch_length = new_branch_length
            updated_node_ids.add(node_id)

            if master_clade.is_terminal():
                matched_tip_edges += 1
            elif node_id != core_record.core_root_node_id:
                matched_internal_edges += 1
            if node_id == core_record.core_root_node_id:
                replaced_core_root_edge = True

            edge_update_rows.append(
                build_edge_update_row(
                    edge_role=edge_role_for_core(node_id, master_clade, core_record.core_root_node_id),
                    core_subtree_id=core_record.core_subtree_id,
                    baseml_subtree_id=baseml_record.baseml_subtree_id,
                    node_id=node_id,
                    descendant_tip_hash=signature,
                    descendant_n_tips=master_signature_counts[signature],
                    old_branch_length=old_branch_length,
                    new_branch_length=new_branch_length,
                    scale_multiplier=scale_multiplier,
                    scale_source_core_ids=[],
                )
            )

        core_root_multipliers[core_record.core_subtree_id] = compute_multiplier(
            original_branch_by_node_id[core_record.core_root_node_id],
            master_core_clade.branch_length,
            min_branch_length,
        )
        merge_report_rows.append(
            {
                "core_subtree_id": core_record.core_subtree_id,
                "core_root_node_id": core_record.core_root_node_id,
                "baseml_subtree_id": baseml_record.baseml_subtree_id,
                "analysis_tree_file": analysis_tree_path.as_posix(),
                "core_n_tips": str(core_record.core_n_tips),
                "analysis_total_n_tips": str(baseml_record.total_n_tips),
                "matched_internal_edges": str(matched_internal_edges),
                "matched_tip_edges": str(matched_tip_edges),
                "replaced_core_root_edge": "true" if replaced_core_root_edge else "false",
                "topology_match": "true",
                "merge_status": "success",
                "details": "core edges replaced",
            }
        )
        logger.info(
            "Merged %s from %s: internal=%d tip=%d",
            core_record.core_subtree_id,
            baseml_record.baseml_subtree_id,
            matched_internal_edges,
            matched_tip_edges,
        )

    if len(core_root_multipliers) != len(core_records):
        raise PipelineError("Not all core subtree root multipliers were computed.")

    core_root_node_to_core_id = {record.core_root_node_id: record.core_subtree_id for record in core_records}
    descendant_core_ids = build_descendant_core_ids(master_tree, core_root_node_to_core_id, node_id_map)
    outgroup_child = get_outgroup_child(master_tree, GLOBAL_OUTGROUP_TIP)
    all_core_ids = [record.core_subtree_id for record in core_records]
    core_tip_weights = {record.core_subtree_id: int(record.core_n_tips) for record in core_records}

    for clade in master_tree.find_clades(order="preorder"):
        if clade is master_tree.root:
            continue
        node_id = node_id_map[clade]
        if node_id in updated_node_ids:
            continue
        if clade is outgroup_child:
            source_core_ids = list(all_core_ids)
        else:
            source_core_ids = list(descendant_core_ids[clade])
        if not source_core_ids:
            raise PipelineError(f"No descendant core roots found for scaffold edge {node_id}")
        multiplier = weighted_geometric_mean(
            [(core_root_multipliers[core_id], core_tip_weights[core_id]) for core_id in source_core_ids]
        )
        old_branch_length = ensure_positive_branch_length(original_branch_by_node_id[node_id], min_branch_length)
        new_branch_length = old_branch_length * multiplier
        clade.branch_length = new_branch_length
        tip_names = sorted(str(tip.name) for tip in clade.get_terminals())
        edge_update_rows.append(
            build_edge_update_row(
                edge_role=edge_role_for_scaffold(clade),
                core_subtree_id="",
                baseml_subtree_id="",
                node_id=node_id,
                descendant_tip_hash=compute_tip_hash(tip_names),
                descendant_n_tips=len(tip_names),
                old_branch_length=old_branch_length,
                new_branch_length=new_branch_length,
                scale_multiplier=multiplier,
                scale_source_core_ids=source_core_ids,
            )
        )

    expected_edge_count = sum(1 for _ in master_tree.find_clades()) - 1
    if len(edge_update_rows) != expected_edge_count:
        raise PipelineError(
            f"edge_update_report row count mismatch: expected {expected_edge_count}, got {len(edge_update_rows)}"
        )

    write_tree_file(master_tree, merge_output_dir / "merged_tree.nwk")
    write_rows(merge_report_rows, MERGE_REPORT_COLUMNS, merge_output_dir / "merge_report.tsv")
    write_rows(edge_update_rows, EDGE_UPDATE_COLUMNS, merge_output_dir / "edge_update_report.tsv")
    logger.info(
        "Merged tree written with %d core root multipliers and %d edge updates.",
        len(core_root_multipliers),
        len(edge_update_rows),
    )
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Merge baseml subtree outputs back onto the master tree.")
    parser.add_argument("--input-tree", required=True, help="Original input tree path.")
    parser.add_argument("--outgroup-tip-file", default=None, help="Outgroup tip file path.")
    parser.add_argument("--split-output-dir", required=True, help="Output directory from the split pipeline.")
    parser.add_argument("--merge-output-dir", required=True, help="Merge output directory.")
    parser.add_argument("--analysis-tree-source", required=True, help="simulated or external.")
    parser.add_argument("--external-result-dir", default=None, help="Directory of external baseml result trees.")
    parser.add_argument("--conda-env", required=True, help="Conda environment name.")
    parser.add_argument("--gotree-bin", required=True, help="gotree executable name.")
    parser.add_argument("--threads", required=True, type=int, help="Thread count for compatibility.")
    parser.add_argument("--preserve-master-topology", action="store_true", help="Preserve the master tree topology.")
    parser.add_argument("--strict-core-topology-match", action="store_true", help="Require exact core topology match.")
    parser.add_argument("--scaffold-scale-aggregation", default="weighted_geometric_mean", help="Scaffold scaling strategy.")
    parser.add_argument("--scaffold-scale-clamp", default=None, help="Optional clamp value. Unsupported when set.")
    parser.add_argument("--min-branch-length", type=float, default=1e-8, help="Minimum positive branch length.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            input_tree=args.input_tree,
            outgroup_tip_file=args.outgroup_tip_file,
            split_output_dir=args.split_output_dir,
            merge_output_dir=args.merge_output_dir,
            analysis_tree_source=args.analysis_tree_source,
            external_result_dir=args.external_result_dir,
            conda_env=args.conda_env,
            gotree_bin=args.gotree_bin,
            threads=args.threads,
            preserve_master_topology=args.preserve_master_topology,
            strict_core_topology_match=args.strict_core_topology_match,
            scaffold_scale_aggregation=args.scaffold_scale_aggregation,
            scaffold_scale_clamp=args.scaffold_scale_clamp,
            min_branch_length=args.min_branch_length,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

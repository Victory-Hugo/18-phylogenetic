#!/usr/bin/env python3
"""Validate the merged tree generated from baseml subtree outputs."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from phylo_merge_common import (
    build_descendant_core_ids,
    build_tree_signature_set,
    decode_json_list,
    get_outgroup_child,
    load_baseml_summary,
    load_core_summary,
    load_table,
    read_newick_tree,
    validate_branch_lengths_complete,
    validate_rooted_binary_tree,
)
from phylo_split_common import GLOBAL_OUTGROUP_TIP, PipelineError, assign_node_ids, setup_logger, write_validation_report


def run(
    input_tree,
    outgroup_tip_file,
    split_output_dir,
    merge_output_dir,
    conda_env,
    gotree_bin,
    threads,
    strict_core_topology_match=True,
    scaffold_scale_aggregation="weighted_geometric_mean",
    log_level="INFO",
):
    del input_tree, outgroup_tip_file, conda_env, gotree_bin, threads

    if not strict_core_topology_match:
        raise PipelineError("This workflow requires strict_core_topology_match=true.")
    if scaffold_scale_aggregation != "weighted_geometric_mean":
        raise PipelineError(f"Unsupported scaffold_scale_aggregation: {scaffold_scale_aggregation}")

    split_output_dir = Path(split_output_dir).resolve()
    merge_output_dir = Path(merge_output_dir).resolve()
    logger = setup_logger("validate_merged_tree", None, log_level)
    logger.info("Starting merged tree validation.")

    master_tree_path = split_output_dir / "intermediate" / "rooted.tree"
    merged_tree_path = merge_output_dir / "merged_tree.nwk"
    if not master_tree_path.exists():
        raise PipelineError(f"Required rooted master tree not found: {master_tree_path}")
    if not merged_tree_path.exists():
        raise PipelineError(f"Merged tree not found: {merged_tree_path}")

    master_tree = read_newick_tree(master_tree_path)
    merged_tree = read_newick_tree(merged_tree_path)
    validate_rooted_binary_tree(master_tree, GLOBAL_OUTGROUP_TIP)
    validate_rooted_binary_tree(merged_tree, GLOBAL_OUTGROUP_TIP)
    validate_branch_lengths_complete(merged_tree)

    core_records = load_core_summary(split_output_dir / "core_subtree_summary.tsv")
    baseml_records = load_baseml_summary(split_output_dir / "baseml_subtree_summary.tsv")
    merge_report_rows = load_table(merge_output_dir / "merge_report.tsv")
    edge_update_rows = load_table(merge_output_dir / "edge_update_report.tsv")

    master_tip_names = [str(tip.name) for tip in master_tree.get_terminals()]
    merged_tip_names = [str(tip.name) for tip in merged_tree.get_terminals()]
    report_rows = []

    if len(merged_tip_names) != len(master_tip_names):
        report_rows.append(("merged_tip_count", "FAIL", "Merged tree tip count differs from master tree"))
    else:
        report_rows.append(("merged_tip_count", "PASS", f"{len(merged_tip_names)} tips"))

    if len(set(merged_tip_names)) != len(merged_tip_names):
        report_rows.append(("merged_tip_duplicates", "FAIL", "Merged tree contains duplicate tip names"))
    else:
        report_rows.append(("merged_tip_duplicates", "PASS", "No duplicate tip names"))

    if set(merged_tip_names) != set(master_tip_names):
        missing = sorted(set(master_tip_names) - set(merged_tip_names))
        extra = sorted(set(merged_tip_names) - set(master_tip_names))
        report_rows.append(("merged_tip_set", "FAIL", f"missing={missing[:5]} extra={extra[:5]}"))
    else:
        report_rows.append(("merged_tip_set", "PASS", "Merged tree tip set matches master tree"))

    report_rows.append(("merged_rooted", "PASS", "Merged tree is rooted"))
    report_rows.append(("merged_outgroup_singleton", "PASS", "RSRS remains a singleton root child"))
    report_rows.append(("merged_binary", "PASS", "Merged tree is strictly bifurcating"))

    if build_tree_signature_set(master_tree) != build_tree_signature_set(merged_tree):
        report_rows.append(("merged_topology_exact", "FAIL", "Merged tree topology differs from master tree"))
    else:
        report_rows.append(("merged_topology_exact", "PASS", "Merged tree topology matches master tree"))

    report_rows.append(("merged_branch_lengths_complete", "PASS", "All non-root branch lengths are present and non-negative"))

    expected_core_ids = {record.core_subtree_id for record in core_records}
    actual_core_ids = {row["core_subtree_id"] for row in merge_report_rows}
    if actual_core_ids != expected_core_ids:
        report_rows.append(("core_mapping_complete", "FAIL", "merge_report.tsv core ids do not match expected core ids"))
    elif any(row["merge_status"] != "success" for row in merge_report_rows):
        report_rows.append(("core_mapping_complete", "FAIL", "merge_report.tsv contains failed merge rows"))
    else:
        report_rows.append(("core_mapping_complete", "PASS", f"{len(actual_core_ids)} core clades merged"))

    expected_baseml_ids = {record.baseml_subtree_id for record in baseml_records}
    actual_baseml_ids = {row["baseml_subtree_id"] for row in merge_report_rows}
    if actual_baseml_ids != expected_baseml_ids:
        report_rows.append(("analysis_tree_coverage", "FAIL", "merge_report.tsv baseml subtree ids do not match expected ids"))
    else:
        report_rows.append(("analysis_tree_coverage", "PASS", f"{len(actual_baseml_ids)} analysis trees covered"))

    master_node_id_map, _, _ = assign_node_ids(master_tree)
    merged_node_id_map, _, _ = assign_node_ids(merged_tree)
    master_node_ids = {node_id for clade, node_id in master_node_id_map.items() if clade is not master_tree.root}
    merged_node_ids = {node_id for clade, node_id in merged_node_id_map.items() if clade is not merged_tree.root}
    edge_rows_by_node_id = {row["node_id"]: row for row in edge_update_rows}
    if set(edge_rows_by_node_id) != master_node_ids or merged_node_ids != master_node_ids:
        report_rows.append(("edge_update_coverage", "FAIL", "edge_update_report.tsv node ids do not cover all non-root edges"))
    else:
        report_rows.append(("edge_update_coverage", "PASS", f"{len(master_node_ids)} non-root edges covered"))

    core_root_node_to_core_id = {record.core_root_node_id: record.core_subtree_id for record in core_records}
    descendant_core_ids = build_descendant_core_ids(master_tree, core_root_node_to_core_id, master_node_id_map)
    updated_core_node_ids = set()
    for record in core_records:
        clade = next(clade for clade, node_id in master_node_id_map.items() if node_id == record.core_root_node_id)
        for descendant in clade.find_clades(order="preorder"):
            updated_core_node_ids.add(master_node_id_map[descendant])
    outgroup_child = get_outgroup_child(master_tree, GLOBAL_OUTGROUP_TIP)

    scaffold_failures = []
    merged_branch_by_node_id = {
        node_id: clade.branch_length
        for clade, node_id in merged_node_id_map.items()
        if clade is not merged_tree.root
    }
    for clade, node_id in master_node_id_map.items():
        if clade is master_tree.root:
            continue
        row = edge_rows_by_node_id.get(node_id)
        if row is None:
            scaffold_failures.append(f"{node_id}: missing row")
            continue
        merged_value = merged_branch_by_node_id[node_id]
        reported_value = float(row["new_branch_length"])
        if abs(float(merged_value) - reported_value) > 1e-5:
            scaffold_failures.append(f"{node_id}: merged tree branch length does not match edge_update_report.tsv")
        if node_id in updated_core_node_ids:
            if not row["edge_role"].startswith("core_"):
                scaffold_failures.append(f"{node_id}: expected core edge role")
        else:
            if not row["edge_role"].startswith("scaffold_"):
                scaffold_failures.append(f"{node_id}: expected scaffold edge role")
            source_core_ids = decode_json_list(row["scale_source_core_ids"])
            if clade is outgroup_child:
                expected_sources = sorted(expected_core_ids)
            else:
                expected_sources = sorted(descendant_core_ids[clade])
            if sorted(source_core_ids) != expected_sources:
                scaffold_failures.append(f"{node_id}: scaffold source core ids mismatch")
    if scaffold_failures:
        report_rows.append(("scaffold_rescaled", "FAIL", "; ".join(scaffold_failures[:10])))
    else:
        report_rows.append(("scaffold_rescaled", "PASS", "All scaffold edges were rescaled and reported correctly"))

    n_master_edges_total = len(master_node_ids)
    n_core_edges_replaced = sum(1 for row in edge_update_rows if row["edge_role"].startswith("core_"))
    n_scaffold_edges_rescaled = sum(1 for row in edge_update_rows if row["edge_role"].startswith("scaffold_"))
    global_core_multiplier_values = [float(row["scale_multiplier"]) for row in edge_update_rows if row["edge_role"] == "core_root"]
    if global_core_multiplier_values:
        global_core_multiplier_geomean = sum(math.log(value) for value in global_core_multiplier_values) / len(global_core_multiplier_values)
        geomean_str = f"{math.exp(global_core_multiplier_geomean):.12g}"
    else:
        geomean_str = "nan"
    report_rows.append(("n_master_edges_total", "PASS", str(n_master_edges_total)))
    report_rows.append(("n_core_edges_replaced", "PASS", str(n_core_edges_replaced)))
    report_rows.append(("n_scaffold_edges_rescaled", "PASS", str(n_scaffold_edges_rescaled)))
    report_rows.append(("n_failed_core_subtrees", "PASS", str(sum(1 for row in merge_report_rows if row["merge_status"] != "success"))))
    report_rows.append(("global_core_multiplier_geomean", "PASS", geomean_str))

    write_validation_report(report_rows, merge_output_dir / "merge_validation_report.tsv")
    failed = [row for row in report_rows if row[1] == "FAIL"]
    if failed:
        details = "; ".join(f"{name}: {detail}" for name, _, detail in failed[:10])
        raise PipelineError(f"Merged tree validation failed with {len(failed)} issue(s): {details}")
    logger.info("Merged tree validation passed.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Validate the merged tree from baseml subtree results.")
    parser.add_argument("--input-tree", required=True, help="Original input tree path.")
    parser.add_argument("--outgroup-tip-file", default=None, help="Outgroup tip file path.")
    parser.add_argument("--split-output-dir", required=True, help="Output directory from the split pipeline.")
    parser.add_argument("--merge-output-dir", required=True, help="Merge output directory.")
    parser.add_argument("--conda-env", required=True, help="Conda environment name.")
    parser.add_argument("--gotree-bin", required=True, help="gotree executable name.")
    parser.add_argument("--threads", required=True, type=int, help="Thread count for compatibility.")
    parser.add_argument("--strict-core-topology-match", action="store_true", help="Require exact core topology match.")
    parser.add_argument("--scaffold-scale-aggregation", default="weighted_geometric_mean", help="Scaffold scaling strategy.")
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
            conda_env=args.conda_env,
            gotree_bin=args.gotree_bin,
            threads=args.threads,
            strict_core_topology_match=args.strict_core_topology_match,
            scaffold_scale_aggregation=args.scaffold_scale_aggregation,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

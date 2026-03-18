#!/usr/bin/env python3
"""Validate backbone-target split outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from phylo_split_common import (
    ANCHOR_MANIFEST_COLUMNS,
    BACKBONE_SUMMARY_COLUMNS,
    PAML_MANIFEST_COLUMNS,
    PAML_SUMMARY_COLUMNS,
    SUBTREE_DESIGN_COLUMNS,
    TARGET_MANIFEST_COLUMNS,
    TARGET_SUMMARY_COLUMNS,
    PipelineError,
    assign_node_ids,
    build_anchor_manifest_rows,
    build_paml_subtrees,
    build_subtree_design_rows,
    build_target_manifest_rows,
    build_target_partition,
    compute_target_partition_profiles,
    compute_tip_hash,
    detect_tree_rooted_status,
    get_tip_names_from_clade,
    get_tip_names_from_tree,
    is_binary_rooted_with_outgroup,
    load_json,
    load_table,
    parse_tip_file,
    prepare_rooted_tree,
    read_newick_tree,
    resolve_backbone_selection,
    resolve_outgroup_tip_name,
    resolve_root_distance_by_clade,
    setup_logger,
    validate_outgroup_tips,
    write_validation_report,
)


def eval_tip_list(value: str):
    parsed = json.loads(value)
    return [str(item) for item in parsed]


def _compare_rows(actual_rows, expected_rows, key_fields, report_rows, label):
    if isinstance(key_fields, str):
        key_fields = [key_fields]

    def normalize(row):
        return {key: str(value) for key, value in row.items()}

    def build_key(row):
        return tuple(str(row[field]) for field in key_fields)

    actual_map = {build_key(row): normalize(row) for row in actual_rows}
    expected_map = {build_key(row): normalize(row) for row in expected_rows}
    if set(actual_map) != set(expected_map):
        missing = sorted(set(expected_map) - set(actual_map))
        extra = sorted(set(actual_map) - set(expected_map))
        report_rows.append((label, "FAIL", f"Key mismatch. missing={missing[:5]} extra={extra[:5]}"))
        return
    report_rows.append((label, "PASS", f"{len(actual_map)} rows present as expected"))
    for key in sorted(expected_map):
        if actual_map[key] != expected_map[key]:
            report_rows.append((f"{label}:{key}", "FAIL", "Row contents do not match expected values"))


def _validate_tree_files(
    actual_rows,
    output_dir: Path,
    tip_field: str,
    tree_field: str,
    outgroup_tip_name: str,
    max_tips: int,
    report_rows,
    label_prefix: str,
    expected_tip_sets=None,
    allowed_extra_tip_set=None,
):
    for row in actual_rows:
        tree_path = output_dir / row[tree_field]
        subtree_id = row.get("paml_subtree_id") or row.get("target_subtree_id") or tree_path.stem
        if not tree_path.exists():
            report_rows.append((f"{label_prefix}_tree_exists", "FAIL", f"{subtree_id}: missing file"))
            continue
        tree = read_newick_tree(tree_path)
        tip_names = get_tip_names_from_tree(tree)
        if expected_tip_sets is None:
            expected_tip_names = set(eval_tip_list(row[tip_field]))
        else:
            expected_tip_names = expected_tip_sets[subtree_id]
        if allowed_extra_tip_set is None:
            comparable_tip_names = set(tip_names)
            if comparable_tip_names != expected_tip_names:
                report_rows.append((f"{label_prefix}_tree_tip_set", "FAIL", f"{subtree_id}: tip set mismatch"))
            else:
                report_rows.append((f"{label_prefix}_tree_tip_set", "PASS", f"{subtree_id}: tip set verified"))
        else:
            comparable_tip_names = {tip_name for tip_name in tip_names if tip_name not in allowed_extra_tip_set}
            unexpected_tip_names = sorted(
                tip_name
                for tip_name in tip_names
                if tip_name not in expected_tip_names and tip_name not in allowed_extra_tip_set
            )
            missing_tip_names = sorted(expected_tip_names - comparable_tip_names)
            if missing_tip_names or unexpected_tip_names:
                report_rows.append(
                    (
                        f"{label_prefix}_tree_tip_set",
                        "FAIL",
                        f"{subtree_id}: missing={missing_tip_names[:5]} unexpected={unexpected_tip_names[:5]}",
                    )
                )
            else:
                report_rows.append(
                    (
                        f"{label_prefix}_tree_tip_set",
                        "PASS",
                        f"{subtree_id}: target tips verified with allowed backbone descendants",
                    )
                )
        if compute_tip_hash(sorted(comparable_tip_names)) != compute_tip_hash(sorted(expected_tip_names)):
            report_rows.append((f"{label_prefix}_tree_hash", "FAIL", f"{subtree_id}: hash mismatch"))
        else:
            report_rows.append((f"{label_prefix}_tree_hash", "PASS", f"{subtree_id}: hash verified"))
        if tree_field == "paml_tree_file":
            ok, reason = is_binary_rooted_with_outgroup(tree, outgroup_tip_name)
            if not ok:
                report_rows.append((f"{label_prefix}_rooted_binary", "FAIL", f"{subtree_id}: {reason}"))
            elif len(tip_names) > max_tips:
                report_rows.append((f"{label_prefix}_max_tips", "FAIL", f"{subtree_id}: exceeds max_tips"))
            else:
                report_rows.append((f"{label_prefix}_rooted_binary", "PASS", f"{subtree_id}: rooted binary"))


def _run_fast_validation(output_dir: Path, rooted_tree, outgroup_tip_name: str, max_tips: int, report_rows):
    actual_backbone_rows = load_table(output_dir / "backbone_summary.tsv")
    actual_target_rows = load_table(output_dir / "target_subtree_summary.tsv")
    actual_target_manifest_rows = load_table(output_dir / "target_tree_manifest.tsv")
    actual_subtree_design_rows = load_table(output_dir / "subtree_design_summary.tsv")
    actual_anchor_manifest_rows = load_table(output_dir / "anchor_manifest.tsv")
    actual_paml_rows = load_table(output_dir / "paml_subtree_summary.tsv")
    actual_paml_manifest_rows = load_table(output_dir / "paml_tree_manifest.tsv")

    required_files = [
        output_dir / "backbone_tips.txt",
        output_dir / "backbone_tree.nwk",
        output_dir / "backbone_summary.tsv",
        output_dir / "target_subtree_summary.tsv",
        output_dir / "target_tree_manifest.tsv",
        output_dir / "subtree_design_summary.tsv",
        output_dir / "anchor_manifest.tsv",
        output_dir / "paml_subtree_summary.tsv",
        output_dir / "paml_tree_manifest.tsv",
    ]
    for path in required_files:
        if not path.exists():
            report_rows.append(("required_files", "FAIL", f"Missing {path.name}"))
    if any(name == "required_files" and status == "FAIL" for name, status, _ in report_rows):
        return

    backbone_tip_names = [row["backbone_tip_name"] for row in actual_backbone_rows]
    backbone_tip_file = [line.strip() for line in (output_dir / "backbone_tips.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
    if backbone_tip_file != backbone_tip_names:
        report_rows.append(("backbone_tip_file", "FAIL", "backbone_tips.txt contents differ from backbone_summary.tsv"))
    else:
        report_rows.append(("backbone_tip_file", "PASS", f"{len(backbone_tip_file)} backbone tips verified"))

    backbone_tree_obj = read_newick_tree(output_dir / "backbone_tree.nwk")
    ok, reason = is_binary_rooted_with_outgroup(backbone_tree_obj, outgroup_tip_name)
    expected_backbone_tip_set = set(backbone_tip_names) | {outgroup_tip_name}
    if set(get_tip_names_from_tree(backbone_tree_obj)) != expected_backbone_tip_set:
        report_rows.append(("backbone_tree_tip_set", "FAIL", "backbone_tree.nwk tip set mismatch"))
    elif not ok:
        report_rows.append(("backbone_tree_rooted_binary", "FAIL", reason))
    else:
        report_rows.append(("backbone_tree_rooted_binary", "PASS", "backbone_tree.nwk verified"))

    target_by_id = {row["target_subtree_id"]: row for row in actual_target_rows}
    paml_by_target_id = {row["target_subtree_id"]: row for row in actual_paml_rows}
    subtree_design_by_target_id = {row["target_subtree_id"]: row for row in actual_subtree_design_rows}
    anchor_manifest_by_target_id = {}
    for row in actual_anchor_manifest_rows:
        anchor_manifest_by_target_id.setdefault(row["target_subtree_id"], []).append(row)
    paml_manifest_by_id = {}
    for row in actual_paml_manifest_rows:
        paml_manifest_by_id.setdefault(row["paml_subtree_id"], []).append(row)

    if set(target_by_id) != set(paml_by_target_id):
        report_rows.append(("target_vs_paml_target_ids", "FAIL", "target_subtree_summary and paml_subtree_summary target ids differ"))
    else:
        report_rows.append(("target_vs_paml_target_ids", "PASS", f"{len(target_by_id)} target ids aligned"))

    expected_target_tip_sets = {
        row["target_subtree_id"]: set(eval_tip_list(row["target_nonbackbone_tip_names"]))
        for row in actual_target_rows
    }
    _validate_tree_files(
        actual_target_rows,
        output_dir,
        "target_nonbackbone_tip_names",
        "target_tree_file",
        outgroup_tip_name,
        max_tips,
        report_rows,
        "target",
        expected_tip_sets=expected_target_tip_sets,
        allowed_extra_tip_set=set(backbone_tip_names),
    )
    _validate_tree_files(actual_paml_rows, output_dir, "total_tip_names", "paml_tree_file", outgroup_tip_name, max_tips, report_rows, "paml")

    nonbackbone_union = set()
    for row in actual_target_rows:
        tip_names = eval_tip_list(row["target_nonbackbone_tip_names"])
        if int(row["target_nonbackbone_n_tips"]) != len(tip_names):
            report_rows.append(("target_summary_counts", "FAIL", f"{row['target_subtree_id']}: target_nonbackbone_n_tips mismatch"))
        elif row["target_tip_hash"] != compute_tip_hash(tip_names):
            report_rows.append(("target_summary_hash", "FAIL", f"{row['target_subtree_id']}: target_tip_hash mismatch"))
        else:
            report_rows.append(("target_summary_counts", "PASS", f"{row['target_subtree_id']}: counts/hash verified"))
        for tip_name in tip_names:
            if tip_name in nonbackbone_union:
                report_rows.append(("target_overlap", "FAIL", f"{tip_name} appears in more than one target subtree"))
                break
            nonbackbone_union.add(tip_name)
    if not any(name == "target_overlap" and status == "FAIL" for name, status, _ in report_rows):
        report_rows.append(("target_overlap", "PASS", f"{len(nonbackbone_union)} non-backbone tips assigned exactly once"))

    for target_id, paml_row in paml_by_target_id.items():
        total_tip_names = eval_tip_list(paml_row["total_tip_names"])
        backbone_tip_names_row = eval_tip_list(paml_row["backbone_tip_names"])
        target_tip_names = eval_tip_list(paml_row["target_nonbackbone_tip_names"])
        global_backbone_tip_names = eval_tip_list(paml_row["global_backbone_tip_names"])
        core_target_tip_names = eval_tip_list(paml_row["core_target_tip_names"])
        local_anchor_tip_names = eval_tip_list(paml_row["local_anchor_tip_names"])
        if int(paml_row["total_n_tips"]) != len(total_tip_names):
            report_rows.append(("paml_summary_counts", "FAIL", f"{paml_row['paml_subtree_id']}: total_n_tips mismatch"))
            continue
        if paml_row["tip_hash"] != compute_tip_hash(total_tip_names):
            report_rows.append(("paml_summary_hash", "FAIL", f"{paml_row['paml_subtree_id']}: tip_hash mismatch"))
            continue
        if int(paml_row["local_anchor_n_tips"]) != len(local_anchor_tip_names):
            report_rows.append(("paml_anchor_counts", "FAIL", f"{paml_row['paml_subtree_id']}: local_anchor_n_tips mismatch"))
            continue
        if int(paml_row["core_target_n_tips"]) != len(core_target_tip_names):
            report_rows.append(("paml_core_counts", "FAIL", f"{paml_row['paml_subtree_id']}: core_target_n_tips mismatch"))
            continue
        if int(paml_row["global_backbone_n_tips"]) != len(global_backbone_tip_names):
            report_rows.append(("paml_backbone_counts", "FAIL", f"{paml_row['paml_subtree_id']}: global_backbone_n_tips mismatch"))
            continue
        if paml_row["outgroup_tip"] not in total_tip_names:
            report_rows.append(("paml_outgroup_presence", "FAIL", f"{paml_row['paml_subtree_id']}: outgroup missing from total_tip_names"))
            continue
        if set(backbone_tip_names_row) != set(global_backbone_tip_names):
            report_rows.append(("paml_backbone_consistency", "FAIL", f"{paml_row['paml_subtree_id']}: backbone_tip_names vs global_backbone_tip_names mismatch"))
            continue
        if set(target_tip_names) != set(core_target_tip_names):
            report_rows.append(("paml_target_consistency", "FAIL", f"{paml_row['paml_subtree_id']}: target_nonbackbone_tip_names vs core_target_tip_names mismatch"))
            continue
        if len(set(total_tip_names)) != len(total_tip_names):
            report_rows.append(("paml_total_tip_uniqueness", "FAIL", f"{paml_row['paml_subtree_id']}: duplicate tips in total_tip_names"))
            continue
        report_rows.append(("paml_summary_counts", "PASS", f"{paml_row['paml_subtree_id']}: counts/hash verified"))

        subtree_design_row = subtree_design_by_target_id.get(target_id)
        if subtree_design_row is None:
            report_rows.append(("subtree_design_presence", "FAIL", f"{target_id}: missing subtree_design row"))
        elif (
            subtree_design_row["paml_tree_file"] != paml_row["paml_tree_file"]
            or subtree_design_row["tip_hash"] != paml_row["tip_hash"]
            or int(subtree_design_row["local_anchor_n_tips"]) != len(local_anchor_tip_names)
        ):
            report_rows.append(("subtree_design_consistency", "FAIL", f"{target_id}: subtree_design row inconsistent with paml_summary"))
        else:
            report_rows.append(("subtree_design_consistency", "PASS", f"{target_id}: subtree_design row verified"))

        manifest_rows = anchor_manifest_by_target_id.get(target_id, [])
        manifest_tip_names = sorted(row["tip_name"] for row in manifest_rows)
        if manifest_tip_names != sorted(local_anchor_tip_names):
            report_rows.append(("anchor_manifest_consistency", "FAIL", f"{target_id}: anchor_manifest tips mismatch"))
        else:
            report_rows.append(("anchor_manifest_consistency", "PASS", f"{target_id}: anchor_manifest verified"))

        manifest_rows = paml_manifest_by_id.get(paml_row["paml_subtree_id"], [])
        if len(manifest_rows) != len(total_tip_names):
            report_rows.append(("paml_manifest_counts", "FAIL", f"{paml_row['paml_subtree_id']}: manifest tip count mismatch"))
        else:
            manifest_tip_names = sorted(row["tip_name"] for row in manifest_rows)
            if manifest_tip_names != sorted(total_tip_names):
                report_rows.append(("paml_manifest_consistency", "FAIL", f"{paml_row['paml_subtree_id']}: manifest tips mismatch"))
            else:
                report_rows.append(("paml_manifest_consistency", "PASS", f"{paml_row['paml_subtree_id']}: manifest verified"))

    master_tip_set = set(get_tip_names_from_tree(rooted_tree)) - {outgroup_tip_name}
    if nonbackbone_union != (master_tip_set - set(backbone_tip_names)):
        missing = sorted((master_tip_set - set(backbone_tip_names)) - nonbackbone_union)
        extra = sorted(nonbackbone_union - (master_tip_set - set(backbone_tip_names)))
        report_rows.append(("target_cover_nonbackbone_tips", "FAIL", f"missing={missing[:5]} extra={extra[:5]}"))
    else:
        report_rows.append(("target_cover_nonbackbone_tips", "PASS", f"{len(nonbackbone_union)} non-backbone tips covered"))


def _run_strict_validation(
    input_tree,
    outgroup_tip_file,
    output_dir,
    max_tips,
    backbone_size,
    conda_env,
    gotree_bin,
    threads,
    backbone_tree,
    backbone_tip_id_file,
    backbone_sampling_strategy,
    target_partition_mode,
    local_anchor_count,
    anchor_selection_strategy,
    pre_binarize_rooted_tree,
    outgroup_tip_name,
    precomputed_context_json,
    logger,
    report_rows,
):
    input_tree = Path(input_tree).resolve()
    outgroup_tip_file = Path(outgroup_tip_file).resolve() if outgroup_tip_file else None
    backbone_tree = Path(backbone_tree).resolve() if backbone_tree else None
    backbone_tip_id_file = Path(backbone_tip_id_file).resolve() if backbone_tip_id_file else None
    precomputed_context = Path(precomputed_context_json).resolve() if precomputed_context_json else None

    original_tree = read_newick_tree(input_tree)
    original_tip_names = get_tip_names_from_tree(original_tree)
    if outgroup_tip_name not in original_tip_names:
        raise PipelineError(f"Required global outgroup tip {outgroup_tip_name} is missing from input tree.")

    rooted_status, _ = detect_tree_rooted_status(input_tree, conda_env, gotree_bin, threads, logger)
    if rooted_status == "unrooted":
        if not outgroup_tip_file:
            raise PipelineError("Input tree is unrooted and --outgroup-tip-file is required.")
        outgroup_tips = parse_tip_file(outgroup_tip_file)
        validate_outgroup_tips(outgroup_tips, original_tip_names)
        if outgroup_tip_name not in outgroup_tips:
            raise PipelineError(f"Outgroup tip file must contain {outgroup_tip_name}.")

    rooted_tree_path = output_dir / "intermediate" / "rooted.validation.tree"
    prepare_rooted_tree(
        input_tree=input_tree,
        rooted_tree=rooted_tree_path,
        conda_env=conda_env,
        gotree_bin=gotree_bin,
        threads=threads,
        logger=logger,
        outgroup_tip_file=outgroup_tip_file,
        rooted_status=rooted_status,
        outgroup_tip_name=outgroup_tip_name,
        pre_binarize_rooted_tree=pre_binarize_rooted_tree,
    )
    rooted_tree = read_newick_tree(rooted_tree_path)
    node_id_map, _, parent_map = assign_node_ids(rooted_tree)

    backbone_tips, backbone_records, _ = resolve_backbone_selection(
        rooted_tree=rooted_tree,
        outgroup_tip_name=outgroup_tip_name,
        backbone_tree=backbone_tree,
        backbone_tip_id_file=backbone_tip_id_file,
        backbone_size=int(backbone_size),
        logger=logger,
    )
    partition_profiles = compute_target_partition_profiles(
        rooted_tree,
        set(backbone_tips),
        outgroup_tip_name,
    )
    target_records = build_target_partition(
        tree=rooted_tree,
        node_id_map=node_id_map,
        parent_map=parent_map,
        backbone_tip_set=set(backbone_tips),
        outgroup_tip=outgroup_tip_name,
        target_capacity=int(max_tips) - len(backbone_tips) - 1 - int(local_anchor_count),
        logger=logger,
        partition_profiles=partition_profiles,
    )
    precomputed_context_payload = load_json(precomputed_context) if precomputed_context and precomputed_context.exists() else None
    root_distance_by_clade, tip_root_distance_by_name = resolve_root_distance_by_clade(
        rooted_tree,
        node_id_map,
        precomputed_context=precomputed_context_payload,
    )
    paml_records, paml_manifest_rows = build_paml_subtrees(
        tree=rooted_tree,
        backbone_tip_names=backbone_tips,
        target_records=target_records,
        parent_map=parent_map,
        ordered_nonbackbone_tips=partition_profiles.ordered_nonbackbone_tips,
        root_distance_by_clade=root_distance_by_clade,
        tip_root_distance_by_name=tip_root_distance_by_name,
        outgroup_tip=outgroup_tip_name,
        max_tips=int(max_tips),
        local_anchor_count=int(local_anchor_count),
        anchor_selection_strategy=anchor_selection_strategy,
        logger=logger,
    )
    expected_backbone_rows = [record.to_row() for record in backbone_records]
    expected_target_rows = [record.to_row() for record in target_records]
    expected_target_manifest_rows = build_target_manifest_rows(target_records)
    expected_subtree_design_rows = build_subtree_design_rows(paml_records, target_records)
    expected_anchor_manifest_rows = build_anchor_manifest_rows(paml_records)
    expected_paml_rows = [record.to_row() for record in paml_records]
    expected_target_tip_sets = {
        record.target_subtree_id: set(get_tip_names_from_clade(record.clade))
        for record in target_records
    }

    actual_backbone_rows = load_table(output_dir / "backbone_summary.tsv")
    actual_target_rows = load_table(output_dir / "target_subtree_summary.tsv")
    actual_target_manifest_rows = load_table(output_dir / "target_tree_manifest.tsv")
    actual_subtree_design_rows = load_table(output_dir / "subtree_design_summary.tsv")
    actual_anchor_manifest_rows = load_table(output_dir / "anchor_manifest.tsv")
    actual_paml_rows = load_table(output_dir / "paml_subtree_summary.tsv")
    actual_paml_manifest_rows = load_table(output_dir / "paml_tree_manifest.tsv")

    _compare_rows(actual_backbone_rows, expected_backbone_rows, "backbone_tip_name", report_rows, "backbone_summary")
    _compare_rows(actual_target_rows, expected_target_rows, "target_subtree_id", report_rows, "target_summary")
    _compare_rows(actual_target_manifest_rows, expected_target_manifest_rows, ["target_subtree_id", "tip_name"], report_rows, "target_manifest")
    _compare_rows(actual_subtree_design_rows, expected_subtree_design_rows, "target_subtree_id", report_rows, "subtree_design")
    _compare_rows(actual_anchor_manifest_rows, expected_anchor_manifest_rows, ["target_subtree_id", "tip_name"], report_rows, "anchor_manifest")
    _compare_rows(actual_paml_rows, expected_paml_rows, "paml_subtree_id", report_rows, "paml_summary")
    _compare_rows(actual_paml_manifest_rows, paml_manifest_rows, ["paml_subtree_id", "tip_name"], report_rows, "paml_manifest")

    _validate_tree_files(
        actual_target_rows,
        output_dir,
        "target_nonbackbone_tip_names",
        "target_tree_file",
        outgroup_tip_name,
        int(max_tips),
        report_rows,
        "target",
        expected_tip_sets=expected_target_tip_sets,
    )
    _validate_tree_files(actual_paml_rows, output_dir, "total_tip_names", "paml_tree_file", outgroup_tip_name, int(max_tips), report_rows, "paml")


def run(
    input_tree,
    outgroup_tip_file,
    output_dir,
    max_tips,
    backbone_size,
    conda_env,
    gotree_bin,
    threads,
    backbone_tree=None,
    backbone_tip_id_file=None,
    backbone_sampling_strategy="greedy_farthest_patristic",
    target_partition_mode="recursive_monophyletic",
    local_anchor_count=16,
    anchor_selection_strategy="nearest_boundary_patristic",
    benchmark_tree_tips=300,
    pre_binarize_rooted_tree=True,
    validation_mode="fast",
    precomputed_context_json=None,
    log_level="INFO",
    outgroup_tip_name=None,
):
    output_dir = Path(output_dir).resolve()
    outgroup_tip_file = Path(outgroup_tip_file).resolve() if outgroup_tip_file else None
    outgroup_tip_name = resolve_outgroup_tip_name(outgroup_tip_file, outgroup_tip_name)
    benchmark_tree_tips = int(benchmark_tree_tips)
    del benchmark_tree_tips
    validation_mode = str(validation_mode).strip().lower()
    if validation_mode not in {"fast", "strict"}:
        raise PipelineError(f"Unsupported validation_mode: {validation_mode}")

    logger = setup_logger("validate_phylo_split", None, log_level)
    rooted_tree = read_newick_tree(output_dir / "intermediate" / "rooted.tree")
    ok, reason = is_binary_rooted_with_outgroup(rooted_tree, outgroup_tip_name)
    report_rows = []
    if not ok:
        report_rows.append(("rooted_tree_rooted_binary", "FAIL", reason))
    else:
        report_rows.append(("rooted_tree_rooted_binary", "PASS", "rooted.tree verified"))

    if validation_mode == "fast":
        _run_fast_validation(output_dir, rooted_tree, outgroup_tip_name, int(max_tips), report_rows)
    else:
        _run_strict_validation(
            input_tree=input_tree,
            outgroup_tip_file=outgroup_tip_file,
            output_dir=output_dir,
            max_tips=max_tips,
            backbone_size=backbone_size,
            conda_env=conda_env,
            gotree_bin=gotree_bin,
            threads=threads,
            backbone_tree=backbone_tree,
            backbone_tip_id_file=backbone_tip_id_file,
            backbone_sampling_strategy=backbone_sampling_strategy,
            target_partition_mode=target_partition_mode,
            local_anchor_count=local_anchor_count,
            anchor_selection_strategy=anchor_selection_strategy,
            pre_binarize_rooted_tree=pre_binarize_rooted_tree,
            outgroup_tip_name=outgroup_tip_name,
            precomputed_context_json=precomputed_context_json,
            logger=logger,
            report_rows=report_rows,
        )

    write_validation_report(report_rows, output_dir / "split_validation_report.tsv")
    failures = [row for row in report_rows if row[1] == "FAIL"]
    if failures:
        details = "; ".join(f"{name}: {detail}" for name, _, detail in failures[:10])
        raise PipelineError(f"Split validation failed with {len(failures)} issue(s): {details}")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Validate backbone-target split outputs.")
    parser.add_argument("--input-tree", required=True, help="Path to the input Newick tree.")
    parser.add_argument("--outgroup-tip-file", default=None, help="Path to the outgroup tip file.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--max-tips", required=True, type=int, help="Maximum total tips per PAML subtree.")
    parser.add_argument("--backbone-size", required=True, type=int, help="Number of ingroup tips in the backbone.")
    parser.add_argument("--threads", required=True, type=int, help="Thread count passed to gotree.")
    parser.add_argument("--conda-env", required=True, help="Conda environment containing gotree.")
    parser.add_argument("--gotree-bin", required=True, help="gotree executable name.")
    parser.add_argument("--backbone-tree", default=None, help="Optional user-provided backbone tree.")
    parser.add_argument("--backbone-tip-id-file", default=None, help="Optional user-provided backbone tip text file.")
    parser.add_argument("--backbone-sampling-strategy", default="greedy_farthest_patristic", help="Backbone sampling strategy.")
    parser.add_argument("--target-partition-mode", default="recursive_monophyletic", help="Target partition mode.")
    parser.add_argument("--local-anchor-count", default=16, type=int, help="Reserved overlap anchors per subtree.")
    parser.add_argument("--anchor-selection-strategy", default="nearest_boundary_patristic", help="Local anchor selection strategy.")
    parser.add_argument("--benchmark-tree-tips", default=300, type=int, help="Reserved benchmark tree size metadata for future full-tree comparisons.")
    parser.add_argument("--validation-mode", default="fast", help="fast or strict validation mode.")
    parser.add_argument("--precomputed-context-json", default=None, help="Optional rooted tree precompute context JSON.")
    parser.add_argument(
        "--pre-binarize-rooted-tree",
        default="true",
        help="Normalize intermediate rooted.validation.tree into a rooted binary tree before validation.",
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    parser.add_argument("--outgroup-tip-name", default=None, help="Singleton outgroup tip retained at the root.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            input_tree=args.input_tree,
            outgroup_tip_file=args.outgroup_tip_file,
            output_dir=args.output_dir,
            max_tips=args.max_tips,
            backbone_size=args.backbone_size,
            conda_env=args.conda_env,
            gotree_bin=args.gotree_bin,
            threads=args.threads,
            backbone_tree=args.backbone_tree,
            backbone_tip_id_file=args.backbone_tip_id_file,
            backbone_sampling_strategy=args.backbone_sampling_strategy,
            target_partition_mode=args.target_partition_mode,
            local_anchor_count=args.local_anchor_count,
            anchor_selection_strategy=args.anchor_selection_strategy,
            benchmark_tree_tips=args.benchmark_tree_tips,
            pre_binarize_rooted_tree=str(args.pre_binarize_rooted_tree).lower() == "true",
            validation_mode=args.validation_mode,
            precomputed_context_json=args.precomputed_context_json,
            log_level=args.log_level,
            outgroup_tip_name=args.outgroup_tip_name,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

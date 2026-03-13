#!/usr/bin/env python3
"""Validate core partitions and baseml-ready subtree outputs.

Dependencies:
  pip install biopython
  conda run -n BigLin gotree --help
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import (
    BASEML_SUMMARY_COLUMNS,
    CORE_SUMMARY_COLUMNS,
    GLOBAL_OUTGROUP_TIP,
    OVERLAP_COLUMNS,
    PipelineError,
    assign_node_ids,
    build_baseml_subtrees,
    build_core_partition,
    build_overlap_rows,
    build_tip_owner_map,
    compute_tip_hash,
    compute_tip_counts,
    decode_json_list,
    detect_tree_rooted_status,
    get_tip_names_from_tree,
    is_binary_rooted_with_outgroup,
    load_table,
    load_root_list,
    load_backbone_tips,
    parse_tip_file,
    prepare_rooted_tree,
    read_newick_tree,
    setup_logger,
    summarize_context,
    validate_outgroup_tips,
    write_validation_report,
)


def compare_row_dicts(actual_rows, expected_rows, key_fields, report_rows, label):
    if isinstance(key_fields, str):
        key_fields = [key_fields]

    def normalize_row(row):
        return {key: str(value) for key, value in row.items()}

    def build_key(row):
        return tuple(row[field] for field in key_fields)

    actual_map = {build_key(row): normalize_row(row) for row in actual_rows}
    expected_map = {build_key(normalize_row(row)): normalize_row(row) for row in expected_rows}
    actual_keys = set(actual_map)
    expected_keys = set(expected_map)
    if actual_keys != expected_keys:
        missing = sorted(expected_keys - actual_keys)
        extra = sorted(actual_keys - expected_keys)
        report_rows.append((label, "FAIL", f"Key mismatch. missing={missing[:5]} extra={extra[:5]}"))
        return
    report_rows.append((label, "PASS", f"{len(actual_keys)} rows present as expected"))
    for key in sorted(expected_keys):
        if actual_map[key] != expected_map[key]:
            report_rows.append((f"{label}:{key}", "FAIL", "Row contents do not match expected values"))
        else:
            report_rows.append((f"{label}:{key}", "PASS", "Row matches expected values"))


def validate_core_tree_files(actual_rows, output_dir: Path, expected_tip_map, effective_core_max_tips: int, report_rows):
    for row in actual_rows:
        subtree_path = output_dir / row["core_tree_file"]
        if not subtree_path.exists():
            report_rows.append(("core_tree_exists", "FAIL", f"{row['core_subtree_id']}: missing file"))
            continue
        tree = read_newick_tree(subtree_path)
        tips = get_tip_names_from_tree(tree)
        if len(tips) != int(row["core_n_tips"]):
            report_rows.append(("core_tree_count", "FAIL", f"{row['core_subtree_id']}: tip count mismatch"))
        else:
            report_rows.append(("core_tree_count", "PASS", f"{row['core_subtree_id']}: tip count verified"))
        if len(tips) > effective_core_max_tips:
            report_rows.append(("core_tree_limit", "FAIL", f"{row['core_subtree_id']}: exceeds effective core max"))
        else:
            report_rows.append(("core_tree_limit", "PASS", f"{row['core_subtree_id']}: within effective core max"))
        if set(tips) != expected_tip_map[row["core_subtree_id"]]:
            report_rows.append(("core_tree_tip_set", "FAIL", f"{row['core_subtree_id']}: tip set mismatch"))
        else:
            report_rows.append(("core_tree_tip_set", "PASS", f"{row['core_subtree_id']}: tip set verified"))
        if compute_tip_hash(tips) != row["tip_hash"]:
            report_rows.append(("core_tree_hash", "FAIL", f"{row['core_subtree_id']}: hash mismatch"))
        else:
            report_rows.append(("core_tree_hash", "PASS", f"{row['core_subtree_id']}: hash verified"))


def validate_baseml_tree_files(actual_rows, output_dir: Path, min_baseml_tips: int, max_tips: int, total_input_tips: int, report_rows):
    for row in actual_rows:
        tree_path = output_dir / row["baseml_tree_file"]
        if not tree_path.exists():
            report_rows.append(("baseml_tree_exists", "FAIL", f"{row['baseml_subtree_id']}: missing file"))
            continue
        tree = read_newick_tree(tree_path)
        tips = get_tip_names_from_tree(tree)
        ok, reason = is_binary_rooted_with_outgroup(tree, GLOBAL_OUTGROUP_TIP)
        if ok:
            report_rows.append(("baseml_rooted_binary", "PASS", f"{row['baseml_subtree_id']}: rooted binary"))
        else:
            report_rows.append(("baseml_rooted_binary", "FAIL", f"{row['baseml_subtree_id']}: {reason}"))
        if tips.count(GLOBAL_OUTGROUP_TIP) != 1:
            report_rows.append(("baseml_outgroup_once", "FAIL", f"{row['baseml_subtree_id']}: RSRS count != 1"))
        else:
            report_rows.append(("baseml_outgroup_once", "PASS", f"{row['baseml_subtree_id']}: RSRS present exactly once"))
        if total_input_tips >= min_baseml_tips and len(tips) < min_baseml_tips:
            report_rows.append(("baseml_min_tips", "FAIL", f"{row['baseml_subtree_id']}: below min_baseml_tips"))
        else:
            report_rows.append(("baseml_min_tips", "PASS", f"{row['baseml_subtree_id']}: min size acceptable"))
        if len(tips) > max_tips:
            report_rows.append(("baseml_max_tips", "FAIL", f"{row['baseml_subtree_id']}: exceeds max_tips"))
        else:
            report_rows.append(("baseml_max_tips", "PASS", f"{row['baseml_subtree_id']}: within max_tips"))
        if compute_tip_hash(tips) != row["tip_hash"]:
            report_rows.append(("baseml_hash", "FAIL", f"{row['baseml_subtree_id']}: hash mismatch"))
        else:
            report_rows.append(("baseml_hash", "PASS", f"{row['baseml_subtree_id']}: hash verified"))


def run(
    input_tree,
    outgroup_tip_file,
    output_dir,
    max_tips,
    min_baseml_tips,
    conda_env,
    gotree_bin,
    threads,
    backbone_tree=None,
    backbone_mode="reference_only",
    log_level="INFO",
    construct_baseml_subtrees=True,
    always_include_outgroup=True,
    anchor_strategy="rsrs_plus_local",
    reserve_slots_for_outgroup=1,
    enable_merge_small_blocks=False,
):
    input_tree = Path(input_tree).resolve()
    outgroup_tip_file = Path(outgroup_tip_file).resolve() if outgroup_tip_file else None
    output_dir = Path(output_dir).resolve()
    backbone_tree = Path(backbone_tree).resolve() if backbone_tree else None

    if backbone_mode != "reference_only":
        raise PipelineError(f"Unsupported backbone mode: {backbone_mode}")
    if not always_include_outgroup:
        raise PipelineError("This workflow requires always_include_outgroup=true.")
    if anchor_strategy != "rsrs_plus_local":
        raise PipelineError(f"Unsupported anchor strategy: {anchor_strategy}")

    logger = setup_logger("validate_phylo_split", None, log_level)
    logger.info("Starting validation for baseml-oriented subtree outputs.")

    original_tree = read_newick_tree(input_tree)
    original_tip_names = get_tip_names_from_tree(original_tree)
    if GLOBAL_OUTGROUP_TIP not in original_tip_names:
        raise PipelineError(f"Required global outgroup tip {GLOBAL_OUTGROUP_TIP} is missing from input tree.")
    if backbone_tree:
        _ = load_backbone_tips(backbone_tree, original_tip_names)

    rooted_status, _ = detect_tree_rooted_status(input_tree, conda_env, gotree_bin, threads, logger)
    if rooted_status == "unrooted":
        if not outgroup_tip_file:
            raise PipelineError("Input tree is unrooted and --outgroup-tip-file is required.")
        outgroup_tips = parse_tip_file(outgroup_tip_file)
        validate_outgroup_tips(outgroup_tips, original_tip_names)
        if GLOBAL_OUTGROUP_TIP not in outgroup_tips:
            raise PipelineError(f"Outgroup tip file must contain {GLOBAL_OUTGROUP_TIP}.")

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
    )
    rooted_tree = read_newick_tree(rooted_tree_path)
    context = summarize_context(rooted_tree, GLOBAL_OUTGROUP_TIP, int(max_tips), int(reserve_slots_for_outgroup))
    node_id_map, _, parent_map = assign_node_ids(rooted_tree)
    tip_counts = compute_tip_counts(rooted_tree)

    expected_core_records = build_core_partition(
        tree=rooted_tree,
        node_id_map=node_id_map,
        parent_map=parent_map,
        tip_counts=tip_counts,
        effective_core_max_tips=int(context["effective_core_max_tips"]),
        outgroup_tip=GLOBAL_OUTGROUP_TIP,
        enable_merge_small_blocks=enable_merge_small_blocks,
        logger=logger,
    )
    expected_baseml_records, expected_manifest_rows = build_baseml_subtrees(
        tree=rooted_tree,
        core_records=expected_core_records,
        parent_map=parent_map,
        tip_counts=tip_counts,
        node_id_map=node_id_map,
        outgroup_tip=GLOBAL_OUTGROUP_TIP,
        min_baseml_tips=int(min_baseml_tips),
        max_tips=int(max_tips),
        logger=logger,
    )
    expected_overlap_rows = build_overlap_rows(expected_manifest_rows, build_tip_owner_map(expected_core_records))

    actual_core_rows = load_table(output_dir / "core_subtree_summary.tsv")
    actual_baseml_rows = load_table(output_dir / "baseml_subtree_summary.tsv")
    actual_manifest_rows = load_table(output_dir / "baseml_tree_manifest.tsv")
    actual_overlap_rows = load_table(output_dir / "overlap_report.tsv")
    actual_root_ids = load_root_list(output_dir / "core_subtree_roots.txt")

    report_rows = []

    expected_core_rows = [record.to_row() for record in expected_core_records]
    expected_baseml_rows = [record.to_row() for record in expected_baseml_records]

    compare_row_dicts(actual_core_rows, expected_core_rows, "core_subtree_id", report_rows, "core_summary")
    compare_row_dicts(actual_baseml_rows, expected_baseml_rows, "baseml_subtree_id", report_rows, "baseml_summary")
    compare_row_dicts(actual_manifest_rows, expected_manifest_rows, ["baseml_subtree_id", "tip_name"], report_rows, "manifest")
    compare_row_dicts(actual_overlap_rows, expected_overlap_rows, "tip_name", report_rows, "overlap")

    expected_root_ids = [record.core_root_node_id for record in expected_core_records]
    if actual_root_ids != expected_root_ids:
        report_rows.append(("core_root_order", "FAIL", "core_subtree_roots.txt order does not match expected values"))
    else:
        report_rows.append(("core_root_order", "PASS", "core_subtree_roots.txt order verified"))

    expected_core_tip_map = {
        record.core_subtree_id: set(record.core_tip_names)
        for record in expected_core_records
    }
    validate_core_tree_files(
        actual_rows=actual_core_rows,
        output_dir=output_dir,
        expected_tip_map=expected_core_tip_map,
        effective_core_max_tips=int(context["effective_core_max_tips"]),
        report_rows=report_rows,
    )
    validate_baseml_tree_files(
        actual_rows=actual_baseml_rows,
        output_dir=output_dir,
        min_baseml_tips=int(min_baseml_tips),
        max_tips=int(max_tips),
        total_input_tips=len(original_tip_names),
        report_rows=report_rows,
    )

    for row in actual_baseml_rows:
        core_tip_names = set(decode_json_list(row["core_tip_names"]))
        anchor_tip_names = set(decode_json_list(row["anchor_tip_names"]))
        total_tip_names = set(decode_json_list(row["total_tip_names"]))
        if GLOBAL_OUTGROUP_TIP in core_tip_names:
            report_rows.append(("baseml_core_outgroup", "FAIL", f"{row['baseml_subtree_id']}: RSRS appears in core tips"))
        else:
            report_rows.append(("baseml_core_outgroup", "PASS", f"{row['baseml_subtree_id']}: RSRS not in core tips"))
        if core_tip_names & anchor_tip_names:
            report_rows.append(("baseml_core_anchor_overlap", "FAIL", f"{row['baseml_subtree_id']}: core/anchor overlap exists"))
        else:
            report_rows.append(("baseml_core_anchor_overlap", "PASS", f"{row['baseml_subtree_id']}: core/anchor disjoint"))
        if not core_tip_names.issubset(total_tip_names):
            report_rows.append(("baseml_core_subset", "FAIL", f"{row['baseml_subtree_id']}: core tips not subset of total"))
        else:
            report_rows.append(("baseml_core_subset", "PASS", f"{row['baseml_subtree_id']}: core tips subset of total"))

    owner_map = build_tip_owner_map(expected_core_records)
    manifest_by_tip = {}
    for row in actual_manifest_rows:
        manifest_by_tip.setdefault(row["tip_name"], []).append(row)
    for tip_name, entries in manifest_by_tip.items():
        if tip_name == GLOBAL_OUTGROUP_TIP:
            if any(entry["tip_role"] != "outgroup" for entry in entries):
                report_rows.append(("overlap_roles", "FAIL", f"{tip_name}: RSRS appears with non-outgroup role"))
            else:
                report_rows.append(("overlap_roles", "PASS", f"{tip_name}: RSRS roles verified"))
            continue
        owner = owner_map.get(tip_name)
        for entry in entries:
            role = entry["tip_role"]
            source_core_subtree_id = entry["source_core_subtree_id"]
            if role == "core" and source_core_subtree_id != owner:
                report_rows.append(("overlap_roles", "FAIL", f"{tip_name}: core role owner mismatch"))
            elif role != "core" and source_core_subtree_id == owner:
                report_rows.append(("overlap_roles", "PASS", f"{tip_name}: owner preserved"))
        non_owner_roles = [entry["tip_role"] for entry in entries if entry["source_core_subtree_id"] != owner]
        if any(role == "core" for role in non_owner_roles):
            report_rows.append(("overlap_roles", "FAIL", f"{tip_name}: appears as core outside owner subtree"))
        else:
            report_rows.append(("overlap_roles", "PASS", f"{tip_name}: overlap roles valid"))

    write_validation_report(report_rows, output_dir / "baseml_validation_report.tsv")
    failed = [row for row in report_rows if row[1] == "FAIL"]
    if failed:
        details = "; ".join(f"{name}: {detail}" for name, _, detail in failed[:10])
        raise PipelineError(f"Validation failed with {len(failed)} issue(s): {details}")

    logger.info("Validation completed successfully for %d core and %d baseml subtrees.", len(actual_core_rows), len(actual_baseml_rows))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Validate core partitions and baseml-ready subtree outputs.")
    parser.add_argument("--input-tree", required=True, help="Path to the input Newick tree.")
    parser.add_argument("--outgroup-tip-file", default=None, help="Path to the outgroup tip file.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--max-tips", required=True, type=int, help="Maximum allowed tips per baseml subtree.")
    parser.add_argument("--min-baseml-tips", required=True, type=int, help="Minimum desired tips per baseml subtree.")
    parser.add_argument("--threads", required=True, type=int, help="Thread count passed to gotree.")
    parser.add_argument("--conda-env", required=True, help="Conda environment containing gotree.")
    parser.add_argument("--gotree-bin", required=True, help="gotree executable name.")
    parser.add_argument("--backbone-tree", default=None, help="Optional backbone tree path.")
    parser.add_argument("--backbone-mode", default="reference_only", help="Backbone handling mode.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    parser.add_argument("--construct-baseml-subtrees", action="store_true", help="Construct overlapping baseml analysis trees.")
    parser.add_argument("--always-include-outgroup", action="store_true", help="Always include RSRS in each baseml subtree.")
    parser.add_argument("--anchor-strategy", default="rsrs_plus_local", help="Anchor construction strategy.")
    parser.add_argument("--reserve-slots-for-outgroup", type=int, default=1, help="Tips reserved from core max to hold the outgroup.")
    parser.add_argument("--enable-merge-small-blocks", action="store_true", help="Enable deterministic core subtree merging.")
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
            min_baseml_tips=args.min_baseml_tips,
            conda_env=args.conda_env,
            gotree_bin=args.gotree_bin,
            threads=args.threads,
            backbone_tree=args.backbone_tree,
            backbone_mode=args.backbone_mode,
            log_level=args.log_level,
            construct_baseml_subtrees=args.construct_baseml_subtrees,
            always_include_outgroup=args.always_include_outgroup,
            anchor_strategy=args.anchor_strategy,
            reserve_slots_for_outgroup=args.reserve_slots_for_outgroup,
            enable_merge_small_blocks=args.enable_merge_small_blocks,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

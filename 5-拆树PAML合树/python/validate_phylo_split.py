#!/usr/bin/env python3
"""Validate backbone-target split outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import (
    BACKBONE_SUMMARY_COLUMNS,
    PAML_MANIFEST_COLUMNS,
    PAML_SUMMARY_COLUMNS,
    TARGET_MANIFEST_COLUMNS,
    TARGET_SUMMARY_COLUMNS,
    PipelineError,
    assign_node_ids,
    build_paml_subtrees,
    build_target_manifest_rows,
    build_target_partition,
    compute_tip_hash,
    detect_tree_rooted_status,
    get_tip_names_from_clade,
    get_tip_names_from_tree,
    is_binary_rooted_with_outgroup,
    load_table,
    parse_tip_file,
    prepare_rooted_tree,
    read_newick_tree,
    resolve_backbone_selection,
    resolve_outgroup_tip_name,
    setup_logger,
    validate_outgroup_tips,
    write_validation_report,
)


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


def _validate_tree_files(actual_rows, output_dir: Path, tip_field: str, tree_field: str, outgroup_tip_name: str, max_tips: int, report_rows, label_prefix: str, expected_tip_sets=None):
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
        if set(tip_names) != expected_tip_names:
            report_rows.append((f"{label_prefix}_tree_tip_set", "FAIL", f"{subtree_id}: tip set mismatch"))
        else:
            report_rows.append((f"{label_prefix}_tree_tip_set", "PASS", f"{subtree_id}: tip set verified"))
        if compute_tip_hash(tip_names) != compute_tip_hash(sorted(expected_tip_names)):
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


def eval_tip_list(value: str):
    import json

    parsed = json.loads(value)
    return [str(item) for item in parsed]


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
    log_level="INFO",
    outgroup_tip_name=None,
):
    input_tree = Path(input_tree).resolve()
    outgroup_tip_file = Path(outgroup_tip_file).resolve() if outgroup_tip_file else None
    output_dir = Path(output_dir).resolve()
    backbone_tree = Path(backbone_tree).resolve() if backbone_tree else None
    backbone_tip_id_file = Path(backbone_tip_id_file).resolve() if backbone_tip_id_file else None
    outgroup_tip_name = resolve_outgroup_tip_name(outgroup_tip_file, outgroup_tip_name)

    logger = setup_logger("validate_phylo_split", None, log_level)
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
    target_records = build_target_partition(
        tree=rooted_tree,
        node_id_map=node_id_map,
        parent_map=parent_map,
        backbone_tip_set=set(backbone_tips),
        outgroup_tip=outgroup_tip_name,
        target_capacity=int(max_tips) - len(backbone_tips) - 1,
        logger=logger,
    )
    paml_records, paml_manifest_rows = build_paml_subtrees(
        tree=rooted_tree,
        backbone_tip_names=backbone_tips,
        target_records=target_records,
        outgroup_tip=outgroup_tip_name,
        max_tips=int(max_tips),
        logger=logger,
    )
    expected_backbone_rows = [record.to_row() for record in backbone_records]
    expected_target_rows = [record.to_row() for record in target_records]
    expected_target_manifest_rows = build_target_manifest_rows(target_records)
    expected_paml_rows = [record.to_row() for record in paml_records]
    expected_target_tip_sets = {
        record.target_subtree_id: set(get_tip_names_from_clade(record.clade))
        for record in target_records
    }

    actual_backbone_rows = load_table(output_dir / "backbone_summary.tsv")
    actual_target_rows = load_table(output_dir / "target_subtree_summary.tsv")
    actual_target_manifest_rows = load_table(output_dir / "target_tree_manifest.tsv")
    actual_paml_rows = load_table(output_dir / "paml_subtree_summary.tsv")
    actual_paml_manifest_rows = load_table(output_dir / "paml_tree_manifest.tsv")

    report_rows = []
    _compare_rows(actual_backbone_rows, expected_backbone_rows, "backbone_tip_name", report_rows, "backbone_summary")
    _compare_rows(actual_target_rows, expected_target_rows, "target_subtree_id", report_rows, "target_summary")
    _compare_rows(actual_target_manifest_rows, expected_target_manifest_rows, ["target_subtree_id", "tip_name"], report_rows, "target_manifest")
    _compare_rows(actual_paml_rows, expected_paml_rows, "paml_subtree_id", report_rows, "paml_summary")
    _compare_rows(actual_paml_manifest_rows, paml_manifest_rows, ["paml_subtree_id", "tip_name"], report_rows, "paml_manifest")

    backbone_tip_file = output_dir / "backbone_tips.txt"
    if not backbone_tip_file.exists():
        report_rows.append(("backbone_tip_file", "FAIL", "backbone_tips.txt is missing"))
    else:
        actual_backbone_tip_file = [line.strip() for line in backbone_tip_file.read_text(encoding="utf-8").splitlines() if line.strip()]
        if actual_backbone_tip_file != backbone_tips:
            report_rows.append(("backbone_tip_file", "FAIL", "backbone_tips.txt contents differ from expected backbone"))
        else:
            report_rows.append(("backbone_tip_file", "PASS", f"{len(actual_backbone_tip_file)} backbone tips verified"))

    backbone_tree_path = output_dir / "backbone_tree.nwk"
    if not backbone_tree_path.exists():
        report_rows.append(("backbone_tree_exists", "FAIL", "backbone_tree.nwk is missing"))
    else:
        backbone_tree_obj = read_newick_tree(backbone_tree_path)
        ok, reason = is_binary_rooted_with_outgroup(backbone_tree_obj, outgroup_tip_name)
        expected_tip_set = set(backbone_tips) | {outgroup_tip_name}
        if set(get_tip_names_from_tree(backbone_tree_obj)) != expected_tip_set:
            report_rows.append(("backbone_tree_tip_set", "FAIL", "backbone_tree.nwk tip set mismatch"))
        elif not ok:
            report_rows.append(("backbone_tree_rooted_binary", "FAIL", reason))
        else:
            report_rows.append(("backbone_tree_rooted_binary", "PASS", "backbone_tree.nwk verified"))

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

    nonbackbone_union = set()
    for row in actual_target_rows:
        for tip_name in eval_tip_list(row["target_nonbackbone_tip_names"]):
            if tip_name in nonbackbone_union:
                report_rows.append(("target_overlap", "FAIL", f"{tip_name} appears in more than one target subtree"))
                break
            nonbackbone_union.add(tip_name)
    if not any(name == "target_overlap" and status == "FAIL" for name, status, _ in report_rows):
        report_rows.append(("target_overlap", "PASS", f"{len(nonbackbone_union)} non-backbone tips assigned exactly once"))

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
            log_level=args.log_level,
            outgroup_tip_name=args.outgroup_tip_name,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

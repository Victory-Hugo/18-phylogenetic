#!/usr/bin/env python3
"""Build backbone-driven target partitions and PAML-ready subtrees."""

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
    build_induced_tree_file,
    build_paml_subtrees,
    build_target_manifest_rows,
    build_target_partition,
    clean_generated_outputs,
    detect_tree_rooted_status,
    get_tip_names_from_tree,
    parse_tip_file,
    prepare_rooted_tree,
    read_newick_tree,
    resolve_backbone_selection,
    resolve_outgroup_tip_name,
    setup_logger,
    validate_outgroup_tips,
    write_table,
    write_tip_list,
    write_tree,
    write_clade_tree,
)


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
    clean_split_output_dir=True,
    log_level="INFO",
    outgroup_tip_name=None,
):
    input_tree = Path(input_tree).resolve()
    outgroup_tip_file = Path(outgroup_tip_file).resolve() if outgroup_tip_file else None
    output_dir = Path(output_dir).resolve()
    backbone_tree = Path(backbone_tree).resolve() if backbone_tree else None
    backbone_tip_id_file = Path(backbone_tip_id_file).resolve() if backbone_tip_id_file else None
    outgroup_tip_name = resolve_outgroup_tip_name(outgroup_tip_file, outgroup_tip_name)
    max_tips = int(max_tips)
    backbone_size = int(backbone_size)

    if backbone_sampling_strategy != "greedy_farthest_patristic":
        raise PipelineError(f"Unsupported runtime.backbone_sampling_strategy: {backbone_sampling_strategy}")
    if target_partition_mode != "recursive_monophyletic":
        raise PipelineError(f"Unsupported runtime.target_partition_mode: {target_partition_mode}")

    if clean_split_output_dir:
        clean_generated_outputs(output_dir)

    for relative in ["target_subtrees", "paml_subtrees", "intermediate"]:
        (output_dir / relative).mkdir(parents=True, exist_ok=True)

    logger = setup_logger("split_phylo_tree", output_dir / "split_tree.log", log_level)
    logger.info("Starting backbone-target split workflow.")
    logger.info("Input tree: %s", input_tree)
    logger.info("max_tips=%d backbone_size=%d", max_tips, backbone_size)

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

    rooted_tree_path = output_dir / "intermediate" / "rooted.tree"
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
    master_tip_order = get_tip_names_from_tree(rooted_tree)
    backbone_tips, backbone_records, backbone_source = resolve_backbone_selection(
        rooted_tree=rooted_tree,
        outgroup_tip_name=outgroup_tip_name,
        backbone_tree=backbone_tree,
        backbone_tip_id_file=backbone_tip_id_file,
        backbone_size=backbone_size,
        logger=logger,
    )
    logger.info("Backbone source: %s", backbone_source)
    logger.info("Backbone tips selected: %d", len(backbone_tips))

    if len(set(backbone_tips)) != len(backbone_tips):
        raise PipelineError("Backbone tip list contains duplicates.")
    if outgroup_tip_name in backbone_tips:
        raise PipelineError("Backbone tips must not contain the singleton outgroup tip.")

    target_capacity = max_tips - len(backbone_tips) - 1
    logger.info("Target capacity per subtree: %d non-backbone tips", target_capacity)
    target_records = build_target_partition(
        tree=rooted_tree,
        node_id_map=node_id_map,
        parent_map=parent_map,
        backbone_tip_set=set(backbone_tips),
        outgroup_tip=outgroup_tip_name,
        target_capacity=target_capacity,
        logger=logger,
    )
    target_manifest_rows = build_target_manifest_rows(target_records)
    paml_records, paml_manifest_rows = build_paml_subtrees(
        tree=rooted_tree,
        backbone_tip_names=backbone_tips,
        target_records=target_records,
        outgroup_tip=outgroup_tip_name,
        max_tips=max_tips,
        logger=logger,
    )

    write_tip_list(backbone_tips, output_dir / "backbone_tips.txt")
    write_table([record.to_row() for record in backbone_records], BACKBONE_SUMMARY_COLUMNS, output_dir / "backbone_summary.tsv")
    write_table([record.to_row() for record in target_records], TARGET_SUMMARY_COLUMNS, output_dir / "target_subtree_summary.tsv")
    write_table(target_manifest_rows, TARGET_MANIFEST_COLUMNS, output_dir / "target_tree_manifest.tsv")
    write_table([record.to_row() for record in paml_records], PAML_SUMMARY_COLUMNS, output_dir / "paml_subtree_summary.tsv")
    write_table(paml_manifest_rows, PAML_MANIFEST_COLUMNS, output_dir / "paml_tree_manifest.tsv")

    backbone_tree_path = output_dir / "backbone_tree.nwk"
    build_induced_tree_file(
        rooted_master_tree=rooted_tree_path,
        selected_tip_names=[tip_name for tip_name in master_tip_order if tip_name in set(backbone_tips) or tip_name == outgroup_tip_name],
        outgroup_tip=outgroup_tip_name,
        destination=backbone_tree_path,
        intermediate_dir=output_dir / "intermediate",
        conda_env=conda_env,
        gotree_bin=gotree_bin,
        threads=int(threads),
        logger=logger,
    )

    for record in target_records:
        write_clade_tree(record.clade, output_dir / record.target_tree_file)
    for record in paml_records:
        build_induced_tree_file(
            rooted_master_tree=rooted_tree_path,
            selected_tip_names=record.total_tip_names,
            outgroup_tip=outgroup_tip_name,
            destination=output_dir / record.paml_tree_file,
            intermediate_dir=output_dir / "intermediate",
            conda_env=conda_env,
            gotree_bin=gotree_bin,
            threads=int(threads),
            logger=logger,
        )

    logger.info("Constructed %d target subtrees and %d PAML subtrees.", len(target_records), len(paml_records))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Build backbone-driven target partitions and PAML-ready subtrees.")
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
    parser.add_argument("--clean-split-output-dir", action="store_true", help="Clean generated outputs before running.")
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
            clean_split_output_dir=args.clean_split_output_dir,
            log_level=args.log_level,
            outgroup_tip_name=args.outgroup_tip_name,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build core partitions and analysis-ready baseml subtrees from a large rooted tree.

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
    MANIFEST_COLUMNS,
    OVERLAP_COLUMNS,
    PipelineError,
    assign_node_ids,
    build_baseml_subtrees,
    build_core_partition,
    build_overlap_rows,
    build_tip_owner_map,
    build_induced_tree_file,
    clean_generated_outputs,
    compute_tip_counts,
    detect_tree_rooted_status,
    get_tip_names_from_tree,
    load_backbone_tips,
    parse_tip_file,
    prepare_rooted_tree,
    read_newick_tree,
    setup_logger,
    summarize_context,
    validate_outgroup_tips,
    write_clade_tree,
    write_root_list,
    write_table,
)


def run(
    input_tree,
    outgroup_tip_file,
    output_dir,
    max_tips,
    min_baseml_tips,
    conda_env,
    gotree_bin,
    threads,
    enable_merge_small_blocks=False,
    backbone_tree=None,
    backbone_mode="reference_only",
    clean_subtree_dir=True,
    log_level="INFO",
    construct_baseml_subtrees=True,
    always_include_outgroup=True,
    anchor_strategy="rsrs_plus_local",
    reserve_slots_for_outgroup=1,
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

    if clean_subtree_dir:
        clean_generated_outputs(output_dir)

    for relative in ["core_subtrees", "baseml_subtrees", "intermediate"]:
        (output_dir / relative).mkdir(parents=True, exist_ok=True)

    logger = setup_logger("split_phylo_tree", output_dir / "split_tree.log", log_level)
    logger.info("Starting baseml-oriented phylogenetic split.")
    logger.info("Input tree: %s", input_tree)
    logger.info("Output dir: %s", output_dir)
    logger.info("max_tips=%s min_baseml_tips=%s reserve_slots_for_outgroup=%s", max_tips, min_baseml_tips, reserve_slots_for_outgroup)

    original_tree = read_newick_tree(input_tree)
    original_tip_names = get_tip_names_from_tree(original_tree)
    logger.info("Original tree tips: %d", len(original_tip_names))
    logger.info("Original tree internal nodes: %d", len(original_tree.get_nonterminals()))

    if GLOBAL_OUTGROUP_TIP not in original_tip_names:
        raise PipelineError(f"Required global outgroup tip {GLOBAL_OUTGROUP_TIP} is missing from input tree.")

    if backbone_tree:
        backbone_tips = load_backbone_tips(backbone_tree, original_tip_names)
        logger.info("Loaded backbone tips: %d", len(backbone_tips))

    rooted_status, _ = detect_tree_rooted_status(input_tree, conda_env, gotree_bin, threads, logger)
    if rooted_status == "unrooted":
        if not outgroup_tip_file:
            raise PipelineError("Input tree is unrooted and --outgroup-tip-file is required.")
        outgroup_tips = parse_tip_file(outgroup_tip_file)
        validate_outgroup_tips(outgroup_tips, original_tip_names)
        if GLOBAL_OUTGROUP_TIP not in outgroup_tips:
            raise PipelineError(f"Outgroup tip file must contain {GLOBAL_OUTGROUP_TIP}.")
        logger.info("Validated outgroup tips for rerooting: %d", len(outgroup_tips))
    else:
        logger.info("Input tree already rooted; using existing rooting.")

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
    context = summarize_context(rooted_tree, GLOBAL_OUTGROUP_TIP, int(max_tips), int(reserve_slots_for_outgroup))
    logger.info("Effective core max tips: %d", context["effective_core_max_tips"])

    node_id_map, _, parent_map = assign_node_ids(rooted_tree)
    tip_counts = compute_tip_counts(rooted_tree)

    core_records = build_core_partition(
        tree=rooted_tree,
        node_id_map=node_id_map,
        parent_map=parent_map,
        tip_counts=tip_counts,
        effective_core_max_tips=int(context["effective_core_max_tips"]),
        outgroup_tip=GLOBAL_OUTGROUP_TIP,
        enable_merge_small_blocks=enable_merge_small_blocks,
        logger=logger,
    )

    for record in core_records:
        write_clade_tree(record.clade, output_dir / record.core_tree_file)
        logger.info("Core %s %s %d", record.core_subtree_id, record.core_root_node_id, record.core_n_tips)

    write_table([record.to_row() for record in core_records], CORE_SUMMARY_COLUMNS, output_dir / "core_subtree_summary.tsv")
    write_root_list([record.core_root_node_id for record in core_records], output_dir / "core_subtree_roots.txt")

    if not construct_baseml_subtrees:
        logger.info("construct_baseml_subtrees=false; stopping after core partition.")
        return 0

    baseml_records, manifest_rows = build_baseml_subtrees(
        tree=rooted_tree,
        core_records=core_records,
        parent_map=parent_map,
        tip_counts=tip_counts,
        node_id_map=node_id_map,
        outgroup_tip=GLOBAL_OUTGROUP_TIP,
        min_baseml_tips=int(min_baseml_tips),
        max_tips=int(max_tips),
        logger=logger,
    )

    for record in baseml_records:
        build_induced_tree_file(
            rooted_master_tree=rooted_tree_path,
            selected_tip_names=record.total_tip_names,
            outgroup_tip=GLOBAL_OUTGROUP_TIP,
            destination=output_dir / record.baseml_tree_file,
            intermediate_dir=output_dir / "intermediate",
            conda_env=conda_env,
            gotree_bin=gotree_bin,
            threads=int(threads),
            logger=logger,
        )

    tip_owner_map = build_tip_owner_map(core_records)
    overlap_rows = build_overlap_rows(manifest_rows, tip_owner_map)

    write_table([record.to_row() for record in baseml_records], BASEML_SUMMARY_COLUMNS, output_dir / "baseml_subtree_summary.tsv")
    write_table(manifest_rows, MANIFEST_COLUMNS, output_dir / "baseml_tree_manifest.tsv")
    write_table(overlap_rows, OVERLAP_COLUMNS, output_dir / "overlap_report.tsv")

    logger.info("Constructed %d core subtrees and %d baseml subtrees.", len(core_records), len(baseml_records))
    logger.info("Split pipeline completed successfully.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Build core partitions and baseml-ready subtrees.")
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
    parser.add_argument("--enable-merge-small-blocks", action="store_true", help="Enable deterministic core subtree merging.")
    parser.add_argument("--clean-subtree-dir", action="store_true", help="Clean generated outputs before running.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    parser.add_argument("--construct-baseml-subtrees", action="store_true", help="Construct overlapping baseml analysis trees.")
    parser.add_argument("--always-include-outgroup", action="store_true", help="Always include RSRS in each baseml subtree.")
    parser.add_argument("--anchor-strategy", default="rsrs_plus_local", help="Anchor construction strategy.")
    parser.add_argument("--reserve-slots-for-outgroup", type=int, default=1, help="Tips reserved from core max to hold the outgroup.")
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
            enable_merge_small_blocks=args.enable_merge_small_blocks,
            backbone_tree=args.backbone_tree,
            backbone_mode=args.backbone_mode,
            clean_subtree_dir=args.clean_subtree_dir,
            log_level=args.log_level,
            construct_baseml_subtrees=args.construct_baseml_subtrees,
            always_include_outgroup=args.always_include_outgroup,
            anchor_strategy=args.anchor_strategy,
            reserve_slots_for_outgroup=args.reserve_slots_for_outgroup,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

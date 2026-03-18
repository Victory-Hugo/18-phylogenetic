#!/usr/bin/env python3
"""Precompute rooted tree context for split-stage acceleration."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import (
    PipelineError,
    assign_node_ids,
    build_split_precompute_payload,
    detect_tree_rooted_status,
    get_tip_names_from_tree,
    parse_tip_file,
    prepare_rooted_tree,
    read_newick_tree,
    resolve_outgroup_tip_name,
    setup_logger,
    validate_outgroup_tips,
    write_json,
)


def run(
    input_tree,
    outgroup_tip_file,
    output_dir,
    conda_env,
    gotree_bin,
    threads,
    pre_binarize_rooted_tree=True,
    log_level="INFO",
    outgroup_tip_name=None,
):
    input_tree = Path(input_tree).resolve()
    outgroup_tip_file = Path(outgroup_tip_file).resolve() if outgroup_tip_file else None
    output_dir = Path(output_dir).resolve()
    outgroup_tip_name = resolve_outgroup_tip_name(outgroup_tip_file, outgroup_tip_name)

    logger = setup_logger("precompute_split_context", output_dir / "split_precompute.log", log_level)
    logger.info("Starting split precompute stage.")

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
        outgroup_tip_name=outgroup_tip_name,
        pre_binarize_rooted_tree=pre_binarize_rooted_tree,
    )

    rooted_tree = read_newick_tree(rooted_tree_path)
    node_id_map, _, _ = assign_node_ids(rooted_tree)
    payload = build_split_precompute_payload(
        input_tree=input_tree,
        rooted_tree=rooted_tree_path,
        tree=rooted_tree,
        node_id_map=node_id_map,
        outgroup_tip_name=outgroup_tip_name,
    )
    context_path = output_dir / "intermediate" / "split_context.json"
    write_json(payload, context_path)
    logger.info("Wrote split precompute context: %s", context_path)
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Precompute rooted tree context for split-stage acceleration.")
    parser.add_argument("--input-tree", required=True, help="Path to the input Newick tree.")
    parser.add_argument("--outgroup-tip-file", default=None, help="Path to the outgroup tip file.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--conda-env", required=True, help="Conda environment containing gotree.")
    parser.add_argument("--gotree-bin", required=True, help="gotree executable name.")
    parser.add_argument("--threads", required=True, type=int, help="Thread count passed to gotree.")
    parser.add_argument(
        "--pre-binarize-rooted-tree",
        default="true",
        help="Normalize intermediate rooted.tree into a rooted binary tree before split logic.",
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
            conda_env=args.conda_env,
            gotree_bin=args.gotree_bin,
            threads=args.threads,
            pre_binarize_rooted_tree=str(args.pre_binarize_rooted_tree).lower() == "true",
            log_level=args.log_level,
            outgroup_tip_name=args.outgroup_tip_name,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

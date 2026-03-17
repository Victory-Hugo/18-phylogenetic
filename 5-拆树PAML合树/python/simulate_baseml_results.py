#!/usr/bin/env python3
"""Simulate PAML baseml subtree outputs by randomizing branch lengths."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import (
    SIMULATION_MANIFEST_COLUMNS,
    build_analysis_tree_path,
    load_baseml_summary,
    randomize_tree_branch_lengths,
    read_newick_tree,
    stable_hash_int,
    validate_rooted_binary_tree,
    validate_tree_against_expected_tip_set,
    write_rows,
    write_tree_file,
)
from phylo_split_common import PipelineError, clone_tree, setup_logger


def run(
    paml_summary_tsv,
    paml_tree_dir,
    output_dir,
    randomization_model,
    randomization_sigma,
    randomization_seed,
    min_branch_length,
    log_level="INFO",
    outgroup_tip_name=None,
):
    paml_summary_tsv = Path(paml_summary_tsv).resolve()
    paml_tree_dir = Path(paml_tree_dir).resolve()
    output_dir = Path(output_dir).resolve()
    simulated_tree_dir = output_dir / "simulated_baseml_subtrees"
    simulated_tree_dir.mkdir(parents=True, exist_ok=True)

    if randomization_model != "lognormal_multiplicative":
        raise PipelineError(f"Unsupported randomization_model: {randomization_model}")

    logger = setup_logger("simulate_baseml_results", output_dir / "merge.log", log_level)
    logger.info("Starting simulated baseml subtree generation.")

    baseml_records = load_baseml_summary(paml_summary_tsv)
    if outgroup_tip_name in (None, "", "null"):
        unique_outgroups = sorted({record.outgroup_tip for record in baseml_records})
        if len(unique_outgroups) != 1:
            raise PipelineError(
                f"Could not infer a unique outgroup tip from {paml_summary_tsv}: {unique_outgroups}"
            )
        outgroup_tip_name = unique_outgroups[0]
    manifest_rows = []

    for record in baseml_records:
        source_tree_path = paml_tree_dir / Path(record.baseml_tree_file).name
        if not source_tree_path.exists():
            raise PipelineError(f"Baseml subtree file not found: {source_tree_path}")
        source_tree = read_newick_tree(source_tree_path)
        validate_tree_against_expected_tip_set(source_tree, record.total_tip_names, record.baseml_subtree_id, outgroup_tip_name)
        validate_rooted_binary_tree(source_tree, outgroup_tip_name)

        tree_seed = int(randomization_seed) + (stable_hash_int(record.baseml_subtree_id) % (2 ** 31))
        simulated_tree = clone_tree(source_tree)
        n_edges, min_multiplier, max_multiplier, mean_multiplier = randomize_tree_branch_lengths(
            simulated_tree,
            sigma=float(randomization_sigma),
            seed=tree_seed,
            min_branch_length=float(min_branch_length),
        )
        destination = build_analysis_tree_path(output_dir, record.baseml_subtree_id, "simulated", None)
        write_tree_file(simulated_tree, destination)
        logger.info(
            "%s randomized %d edges seed=%d multiplier_range=[%.6f, %.6f]",
            record.baseml_subtree_id,
            n_edges,
            tree_seed,
            min_multiplier,
            max_multiplier,
        )
        manifest_rows.append(
            {
                "baseml_subtree_id": record.baseml_subtree_id,
                "source_tree_file": source_tree_path.as_posix(),
                "simulated_tree_file": destination.as_posix(),
                "tree_seed": str(tree_seed),
                "n_edges_randomized": str(n_edges),
                "min_multiplier": f"{min_multiplier:.12g}",
                "max_multiplier": f"{max_multiplier:.12g}",
                "mean_multiplier": f"{mean_multiplier:.12g}",
            }
        )

    write_rows(manifest_rows, SIMULATION_MANIFEST_COLUMNS, output_dir / "simulation_manifest.tsv")
    logger.info("Simulated %d baseml subtree result files.", len(manifest_rows))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Simulate baseml subtree outputs by randomizing branch lengths.")
    parser.add_argument("--paml-summary-tsv", required=True, help="Path to paml_subtree_summary.tsv.")
    parser.add_argument("--paml-tree-dir", required=True, help="Directory containing PAML subtree tree files.")
    parser.add_argument("--output-dir", required=True, help="Merge output directory.")
    parser.add_argument("--randomization-model", required=True, help="Randomization model name.")
    parser.add_argument("--randomization-sigma", required=True, type=float, help="Sigma for lognormal perturbation.")
    parser.add_argument("--randomization-seed", required=True, type=int, help="Base random seed.")
    parser.add_argument("--min-branch-length", required=True, type=float, help="Minimum branch length floor.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    parser.add_argument("--outgroup-tip-name", default=None, help="Singleton outgroup tip retained at the root.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            paml_summary_tsv=args.paml_summary_tsv,
            paml_tree_dir=args.paml_tree_dir,
            output_dir=args.output_dir,
            randomization_model=args.randomization_model,
            randomization_sigma=args.randomization_sigma,
            randomization_seed=args.randomization_seed,
            min_branch_length=args.min_branch_length,
            log_level=args.log_level,
            outgroup_tip_name=args.outgroup_tip_name,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

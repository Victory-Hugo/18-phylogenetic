#!/usr/bin/env python3
"""Compare rooted tree topologies using non-trivial clade signatures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import build_clade_signature_maps
from phylo_split_common import PipelineError, get_tip_names_from_tree, setup_logger
from phylo_ultrastandard_common import read_tree


def build_nontrivial_signature_set(tree) -> set[str]:
    total_tip_count = len(get_tip_names_from_tree(tree))
    signature_map, signature_to_count = build_clade_signature_maps(tree.root)
    return {
        signature
        for signature in signature_map
        if 1 < signature_to_count[signature] < total_tip_count
    }


def run(reference_tree, query_tree, log_level="INFO") -> int:
    reference_tree = Path(reference_tree).resolve()
    query_tree = Path(query_tree).resolve()
    logger = setup_logger("check_tree_topology_match", None, log_level)

    if not reference_tree.exists():
        raise PipelineError(f"Reference tree not found: {reference_tree}")
    if not query_tree.exists():
        raise PipelineError(f"Query tree not found: {query_tree}")

    reference = read_tree(reference_tree)
    query = read_tree(query_tree)

    reference_tip_set = set(get_tip_names_from_tree(reference))
    query_tip_set = set(get_tip_names_from_tree(query))
    if reference_tip_set != query_tip_set:
        missing = sorted(reference_tip_set - query_tip_set)
        extra = sorted(query_tip_set - reference_tip_set)
        logger.warning("Topology comparison skipped because tip sets differ.")
        print(
            f"[WARN] Tree tip sets differ. missing={missing[:5]} extra={extra[:5]}",
            file=sys.stderr,
        )
        return 2

    reference_signatures = build_nontrivial_signature_set(reference)
    query_signatures = build_nontrivial_signature_set(query)
    if reference_signatures == query_signatures:
        logger.info(
            "Tree topology match verified with %d shared non-trivial clades.",
            len(reference_signatures),
        )
        print(
            f"[OK] Tree topology matches: {len(reference_signatures)} non-trivial clades verified."
        )
        return 0

    missing = sorted(reference_signatures - query_signatures)
    extra = sorted(query_signatures - reference_signatures)
    logger.warning(
        "Tree topology mismatch detected. missing_clades=%d extra_clades=%d",
        len(missing),
        len(extra),
    )
    print(
        "[WARN] Tree topology mismatch detected: "
        f"missing_clades={len(missing)} extra_clades={len(extra)}.",
        file=sys.stderr,
    )
    return 2


def build_parser():
    parser = argparse.ArgumentParser(description="Compare rooted tree topologies.")
    parser.add_argument("--reference-tree", required=True, help="Reference rooted tree.")
    parser.add_argument("--query-tree", required=True, help="Tree to compare against the reference.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            reference_tree=args.reference_tree,
            query_tree=args.query_tree,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Prepare a single backbone-only baseml job."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import PipelineError, get_tip_names_from_tree, read_newick_tree, setup_logger
from phylo_ultrastandard_common import BACKBONE_ANALYSIS_MANIFEST_COLUMNS, infer_unique_outgroup
from paml_common import (
    map_tip_to_source_id,
    read_fasta_record_map,
    render_ctl_template,
    resolve_seq_id_strategy,
    write_paml_treefile,
    write_rows,
    write_subtree_fasta,
)


def run(
    backbone_tree,
    backbone_summary_tsv,
    paml_summary_tsv,
    input_fasta,
    ctl_template,
    job_dir,
    analysis_dir,
    seq_id_strategy,
    log_level="INFO",
):
    backbone_tree = Path(backbone_tree).resolve()
    backbone_summary_tsv = Path(backbone_summary_tsv).resolve()
    paml_summary_tsv = Path(paml_summary_tsv).resolve()
    input_fasta = Path(input_fasta).resolve()
    ctl_template = Path(ctl_template).resolve()
    job_dir = Path(job_dir).resolve()
    analysis_dir = Path(analysis_dir).resolve()

    logger = setup_logger("prepare_backbone_paml_input", analysis_dir / "prepare_backbone_paml_input.log", log_level)
    strategy = resolve_seq_id_strategy(seq_id_strategy)
    outgroup_tip = infer_unique_outgroup(paml_summary_tsv)

    if not backbone_tree.exists():
        raise PipelineError(f"Backbone tree not found: {backbone_tree}")
    if not backbone_summary_tsv.exists():
        raise PipelineError(f"Backbone summary not found: {backbone_summary_tsv}")

    tree = read_newick_tree(backbone_tree)
    tree_tip_names = get_tip_names_from_tree(tree)
    if outgroup_tip not in tree_tip_names:
        raise PipelineError(f"Backbone tree does not contain inferred outgroup tip: {outgroup_tip}")

    fasta_record_map = read_fasta_record_map(input_fasta)
    tip_to_source_id = map_tip_to_source_id(tree_tip_names, fasta_record_map, strategy)

    job_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    treefile = job_dir / "subtree.treefile"
    seqfile = job_dir / "subtree.fasta"
    ctlfile = job_dir / "baseml.ctl"
    outfile = job_dir / "mlb"

    write_paml_treefile(backbone_tree, len(tree_tip_names), treefile)
    write_subtree_fasta(tree_tip_names, fasta_record_map, tip_to_source_id, seqfile)
    render_ctl_template(
        template_path=ctl_template,
        seqfile_name=seqfile.name,
        treefile_name=treefile.name,
        outfile_name=outfile.name,
        destination=ctlfile,
    )

    manifest_rows = []
    for tip_name in tree_tip_names:
        manifest_rows.append(
            {
                "backbone_job_id": "backbone_only",
                "tip_name": tip_name,
                "tip_role": "outgroup" if tip_name == outgroup_tip else "backbone",
                "source_id": tip_to_source_id[tip_name],
            }
        )
    write_rows(manifest_rows, BACKBONE_ANALYSIS_MANIFEST_COLUMNS, analysis_dir / "backbone_analysis_manifest.tsv")
    logger.info("Prepared backbone-only baseml job with %d tips.", len(tree_tip_names))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Prepare a single backbone-only baseml job.")
    parser.add_argument("--backbone-tree", required=True, help="Path to split backbone_tree.nwk.")
    parser.add_argument("--backbone-summary-tsv", required=True, help="Path to backbone_summary.tsv.")
    parser.add_argument("--paml-summary-tsv", required=True, help="Path to paml_subtree_summary.tsv.")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA containing all sequences.")
    parser.add_argument("--ctl-template", required=True, help="Template baseml ctl file.")
    parser.add_argument("--job-dir", required=True, help="Directory for the backbone baseml job.")
    parser.add_argument("--analysis-dir", required=True, help="Directory for backbone analysis outputs.")
    parser.add_argument("--seq-id-strategy", required=True, help="exact or prefix_before_first_underscore.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            backbone_tree=args.backbone_tree,
            backbone_summary_tsv=args.backbone_summary_tsv,
            paml_summary_tsv=args.paml_summary_tsv,
            input_fasta=args.input_fasta,
            ctl_template=args.ctl_template,
            job_dir=args.job_dir,
            analysis_dir=args.analysis_dir,
            seq_id_strategy=args.seq_id_strategy,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

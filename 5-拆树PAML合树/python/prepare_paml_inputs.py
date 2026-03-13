#!/usr/bin/env python3
"""Prepare baseml job directories from split subtree outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from paml_common import (
    PAML_JOB_MANIFEST_COLUMNS,
    load_baseml_records,
    map_tip_to_source_id,
    read_fasta_record_map,
    render_ctl_template,
    resolve_seq_id_strategy,
    write_paml_treefile,
    write_rows,
    write_subtree_fasta,
)


def run(
    paml_summary_tsv,
    paml_tree_dir,
    input_fasta,
    ctl_template,
    output_dir,
    seq_id_strategy,
    log_level="INFO",
):
    paml_summary_tsv = Path(paml_summary_tsv).resolve()
    paml_tree_dir = Path(paml_tree_dir).resolve()
    input_fasta = Path(input_fasta).resolve()
    ctl_template = Path(ctl_template).resolve()
    output_dir = Path(output_dir).resolve()
    jobs_dir = output_dir / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)

    strategy = resolve_seq_id_strategy(seq_id_strategy)
    logger = setup_logger("prepare_paml_inputs", output_dir / "prepare_paml_inputs.log", log_level)
    logger.info("Preparing baseml inputs from %s", paml_summary_tsv)

    baseml_records = load_baseml_records(paml_summary_tsv)
    fasta_record_map = read_fasta_record_map(input_fasta)

    manifest_rows = []
    for record in baseml_records:
        job_dir = jobs_dir / record.baseml_subtree_id
        job_dir.mkdir(parents=True, exist_ok=True)
        treefile = job_dir / "subtree.treefile"
        seqfile = job_dir / "subtree.fasta"
        ctlfile = job_dir / "baseml.ctl"
        outfile = job_dir / "mlb"

        source_tree_path = paml_tree_dir / Path(record.baseml_tree_file).name
        if not source_tree_path.exists():
            raise PipelineError(f"Baseml subtree file not found: {source_tree_path}")

        tip_to_source_id = map_tip_to_source_id(record.total_tip_names, fasta_record_map, strategy)
        write_paml_treefile(source_tree_path, record.total_n_tips, treefile)
        write_subtree_fasta(record.total_tip_names, fasta_record_map, tip_to_source_id, seqfile)
        render_ctl_template(
            template_path=ctl_template,
            seqfile_name=seqfile.name,
            treefile_name=treefile.name,
            outfile_name=outfile.name,
            destination=ctlfile,
        )

        manifest_rows.append(
            {
                "baseml_subtree_id": record.baseml_subtree_id,
                "job_dir": job_dir.as_posix(),
                "treefile": treefile.as_posix(),
                "seqfile": seqfile.as_posix(),
                "ctlfile": ctlfile.as_posix(),
                "outfile": outfile.as_posix(),
                "expected_tip_count": str(record.total_n_tips),
                "tip_hash": record.tip_hash,
                "seq_id_strategy": strategy,
            }
        )
        logger.info(
            "Prepared %s with %d tips using strategy=%s",
            record.baseml_subtree_id,
            record.total_n_tips,
            strategy,
        )

    manifest_path = output_dir / "paml_job_manifest.tsv"
    write_rows(manifest_rows, PAML_JOB_MANIFEST_COLUMNS, manifest_path)
    logger.info("Prepared %d baseml jobs.", len(manifest_rows))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Prepare baseml job directories.")
    parser.add_argument("--paml-summary-tsv", required=True, help="Path to paml_subtree_summary.tsv.")
    parser.add_argument("--paml-tree-dir", required=True, help="Directory containing PAML subtree trees.")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA containing all sequences.")
    parser.add_argument("--ctl-template", required=True, help="Template baseml ctl file.")
    parser.add_argument("--output-dir", required=True, help="PAML output directory.")
    parser.add_argument("--seq-id-strategy", required=True, help="exact or prefix_before_first_underscore.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            paml_summary_tsv=args.paml_summary_tsv,
            paml_tree_dir=args.paml_tree_dir,
            input_fasta=args.input_fasta,
            ctl_template=args.ctl_template,
            output_dir=args.output_dir,
            seq_id_strategy=args.seq_id_strategy,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

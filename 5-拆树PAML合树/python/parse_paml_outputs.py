#!/usr/bin/env python3
"""Parse baseml outputs into analysis trees and tabular summaries."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from phylo_merge_common import load_paml_summary, validate_tip_hash, validate_tree_against_expected_tip_set
from phylo_split_common import PipelineError, read_newick_tree, setup_logger
from paml_common import (
    PAML_PARSE_STATUS_COLUMNS,
    build_parameter_rows,
    extract_tip_length_rows,
    parse_baseml_output,
    write_analysis_tree,
    write_rows,
)


def write_csv_rows(rows, fieldnames, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def dump_parse_artifacts(parsed_output, destination_dir: Path) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)
    (destination_dir / "branch_tags.txt").write_text(parsed_output.branch_tags_line + "\n", encoding="utf-8")
    (destination_dir / "branch_params.txt").write_text(parsed_output.branch_params_line + "\n", encoding="utf-8")
    (destination_dir / "se.txt").write_text(parsed_output.se_line + "\n", encoding="utf-8")
    (destination_dir / "named_tree.nwk").write_text(parsed_output.named_tree_line + "\n", encoding="utf-8")
    (destination_dir / "indexed_tree.nwk").write_text(parsed_output.indexed_tree_line + "\n", encoding="utf-8")


def run(
    paml_summary_tsv,
    manifest_path,
    output_dir,
    tmp_dir,
    outgroup_tip_name=None,
    log_level="INFO",
):
    paml_summary_tsv = Path(paml_summary_tsv).resolve()
    manifest_path = Path(manifest_path).resolve()
    output_dir = Path(output_dir).resolve()
    tmp_dir = Path(tmp_dir).resolve()

    logger = setup_logger("parse_paml_outputs", output_dir / "parse_paml_outputs.log", log_level)
    logger.info("Parsing baseml outputs from manifest: %s", manifest_path)

    baseml_summary_records = load_paml_summary(paml_summary_tsv)
    if outgroup_tip_name in (None, "", "null"):
        unique_outgroups = sorted({record.outgroup_tip for record in baseml_summary_records})
        if len(unique_outgroups) != 1:
            raise PipelineError(
                f"Could not infer a unique outgroup tip from {paml_summary_tsv}: {unique_outgroups}"
            )
        outgroup_tip_name = unique_outgroups[0]
    baseml_records = {
        record.baseml_subtree_id: record
        for record in baseml_summary_records
    }
    manifest_rows = list(csv.DictReader(manifest_path.open("r", encoding="utf-8"), delimiter="\t"))
    if not manifest_rows:
        raise PipelineError(f"PAML job manifest is empty: {manifest_path}")

    analysis_tree_dir = output_dir / "analysis_trees"
    tables_dir = output_dir / "tables"
    parse_tmp_root = tmp_dir / "paml_parse"
    analysis_tree_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    parse_tmp_root.mkdir(parents=True, exist_ok=True)

    status_rows = []
    failures = []

    for row in manifest_rows:
        baseml_subtree_id = row["baseml_subtree_id"]
        record = baseml_records.get(baseml_subtree_id)
        if record is None:
            raise PipelineError(f"Manifest subtree id not present in baseml summary: {baseml_subtree_id}")
        if record.outgroup_tip != outgroup_tip_name:
            raise PipelineError(
                f"Outgroup mismatch for {baseml_subtree_id}: summary={record.outgroup_tip}, expected={outgroup_tip_name}"
            )

        result_file = Path(row["outfile"]).resolve()
        analysis_tree_file = analysis_tree_dir / f"{baseml_subtree_id}.nwk"
        parameters_csv = tables_dir / baseml_subtree_id / "Parameters.csv"
        tip_length_csv = tables_dir / baseml_subtree_id / "TIP_Length.csv"
        parse_tmp_dir = parse_tmp_root / baseml_subtree_id

        try:
            if not result_file.exists() or result_file.stat().st_size == 0:
                raise PipelineError(f"Missing or empty baseml result file: {result_file}")

            parsed_output = parse_baseml_output(result_file)
            dump_parse_artifacts(parsed_output, parse_tmp_dir)

            parameter_rows = build_parameter_rows(parsed_output)
            tip_length_rows = extract_tip_length_rows(parsed_output.named_tree_line)
            write_analysis_tree(parsed_output.named_tree_line, analysis_tree_file)
            write_csv_rows(
                parameter_rows,
                ["Branch", "Param", "SE", "Parent", "Child", "Parent_ID", "Child_ID"],
                parameters_csv,
            )
            write_csv_rows(tip_length_rows, ["ID", "Length"], tip_length_csv)

            analysis_tree = read_newick_tree(analysis_tree_file)
            validate_tree_against_expected_tip_set(
                analysis_tree,
                record.total_tip_names,
                f"analysis tree {baseml_subtree_id}",
                outgroup_tip_name,
            )
            validate_tip_hash(record.total_tip_names, record.tip_hash, f"analysis tree {baseml_subtree_id}")

            actual_tip_count = len(analysis_tree.get_terminals())
            if actual_tip_count != record.total_n_tips:
                raise PipelineError(
                    f"{baseml_subtree_id}: expected {record.total_n_tips} tips, got {actual_tip_count}"
                )

            status_rows.append(
                {
                    "baseml_subtree_id": baseml_subtree_id,
                    "result_file": result_file.as_posix(),
                    "analysis_tree_file": analysis_tree_file.as_posix(),
                    "parameters_csv": parameters_csv.as_posix(),
                    "tip_length_csv": tip_length_csv.as_posix(),
                    "status": "success",
                    "detail": "Parsed successfully",
                }
            )
            logger.info("Parsed %s successfully", baseml_subtree_id)
        except Exception as exc:
            detail = str(exc)
            failures.append((baseml_subtree_id, detail))
            status_rows.append(
                {
                    "baseml_subtree_id": baseml_subtree_id,
                    "result_file": result_file.as_posix(),
                    "analysis_tree_file": analysis_tree_file.as_posix(),
                    "parameters_csv": parameters_csv.as_posix(),
                    "tip_length_csv": tip_length_csv.as_posix(),
                    "status": "failed",
                    "detail": detail,
                }
            )
            logger.error("Failed to parse %s: %s", baseml_subtree_id, detail)

    status_rows.sort(key=lambda item: item["baseml_subtree_id"])
    write_rows(status_rows, PAML_PARSE_STATUS_COLUMNS, output_dir / "paml_parse_status.tsv")

    if failures:
        preview = "; ".join(f"{subtree_id}: {detail}" for subtree_id, detail in failures[:5])
        raise PipelineError(f"{len(failures)} baseml output(s) failed to parse: {preview}")

    logger.info("Parsed %d baseml outputs successfully.", len(status_rows))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Parse baseml outputs into merge-ready trees.")
    parser.add_argument("--paml-summary-tsv", required=True, help="Path to paml_subtree_summary.tsv.")
    parser.add_argument("--manifest", required=True, help="Path to paml_job_manifest.tsv.")
    parser.add_argument("--output-dir", required=True, help="PAML output directory.")
    parser.add_argument("--tmp-dir", required=True, help="Pipeline temporary directory.")
    parser.add_argument("--outgroup-tip-name", default=None, help="Global outgroup tip name.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            paml_summary_tsv=args.paml_summary_tsv,
            manifest_path=args.manifest,
            output_dir=args.output_dir,
            tmp_dir=args.tmp_dir,
            outgroup_tip_name=args.outgroup_tip_name,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

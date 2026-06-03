#!/usr/bin/env python3
"""根据 stage1 拆树结果为每棵子树准备 BEAST 作业目录。

对应原 PAML 版本的 ``prepare_paml_inputs.py``：逐子树写出固定拓扑的超度量起始树、
子树 FASTA，并渲染 BEAST XML（注入 taxa / alignment / 起始树 / 文件名 / 链长），
最后生成作业清单 ``beast_job_manifest.tsv``。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from beast_common import (
    BEAST_JOB_MANIFEST_COLUMNS,
    build_alignment_block,
    build_calibrated_starting_newick,
    build_calibration_block,
    build_taxa_block,
    build_ultrametric_starting_newick,
    classify_m_n_tips,
    load_beast_records,
    map_tip_to_source_id,
    read_fasta_record_map,
    render_beast_xml,
    render_calibrated_beast_xml,
    resolve_seq_id_strategy,
    write_status_rows,
    write_subtree_fasta,
)


def run(
    beast_summary_tsv,
    beast_tree_dir,
    input_fasta,
    xml_template,
    output_dir,
    seq_id_strategy,
    chain_length,
    log_every,
    clock_rate=1.0,
    calibration_enabled=False,
    calibrated_xml_template=None,
    calibration_id_hap_tsv=None,
    calibration_phylotree_json=None,
    calibration_warm_start=True,
    calibration_popsize_start=100000.0,
    calibration_root_age=180000.0,
    calibration_root_stdev=20000.0,
    calibration_m_age=50000.0,
    calibration_m_stdev=5000.0,
    calibration_n_age=57000.0,
    calibration_n_stdev=5000.0,
    min_branch_length=1e-8,
    ultrametric_tolerance=1e-8,
    log_level="INFO",
):
    beast_summary_tsv = Path(beast_summary_tsv).resolve()
    beast_tree_dir = Path(beast_tree_dir).resolve()
    input_fasta = Path(input_fasta).resolve()
    xml_template = Path(xml_template).resolve()
    output_dir = Path(output_dir).resolve()
    calibration_enabled = _as_bool(calibration_enabled)
    calibrated_xml_template = Path(calibrated_xml_template).resolve() if calibrated_xml_template else None
    calibration_id_hap_tsv = Path(calibration_id_hap_tsv).resolve() if calibration_id_hap_tsv else None
    calibration_phylotree_json = Path(calibration_phylotree_json).resolve() if calibration_phylotree_json else None
    jobs_dir = output_dir / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)

    strategy = resolve_seq_id_strategy(seq_id_strategy)
    logger = setup_logger("prepare_beast_inputs", output_dir / "prepare_beast_inputs.log", log_level)
    logger.info("Preparing BEAST inputs from %s", beast_summary_tsv)

    records = load_beast_records(beast_summary_tsv)
    fasta_record_map = read_fasta_record_map(input_fasta)
    if calibration_enabled:
        for required_path, label in (
            (calibrated_xml_template, "calibrated XML template"),
            (calibration_id_hap_tsv, "calibration id_hap_tsv"),
            (calibration_phylotree_json, "calibration phylotree_json"),
        ):
            if required_path is None or not required_path.exists():
                raise PipelineError(f"Missing {label}: {required_path}")

    manifest_rows = []
    for record in records:
        subtree_id = record.baseml_subtree_id
        job_dir = jobs_dir / subtree_id
        job_dir.mkdir(parents=True, exist_ok=True)

        file_prefix = subtree_id
        starting_tree_file = job_dir / "starting_tree.nwk"
        seqfile = job_dir / "subtree.fasta"
        xmlfile = job_dir / "beast.xml"
        trees_file = job_dir / f"{file_prefix}.trees"
        mcc_file = job_dir / f"{file_prefix}.mcc.tree"

        source_tree_path = beast_tree_dir / Path(record.baseml_tree_file).name
        if not source_tree_path.exists():
            raise PipelineError(f"Subtree file not found: {source_tree_path}")

        tip_to_source_id = map_tip_to_source_id(record.total_tip_names, fasta_record_map, strategy)

        # 1) 子树 FASTA
        write_subtree_fasta(record.total_tip_names, fasta_record_map, tip_to_source_id, seqfile)

        # 2) 固定拓扑起始树 + XML。校准模式使用年尺度 UCLD 热启动；否则使用 strict clock substitutions/site。
        taxa_block = build_taxa_block(record.total_tip_names)
        alignment_block = build_alignment_block(record.total_tip_names, fasta_record_map, tip_to_source_id)
        if calibration_enabled:
            m_tips, n_tips, missing_tips = classify_m_n_tips(
                record.total_tip_names,
                id_hap_tsv=calibration_id_hap_tsv,
                phylotree_json=calibration_phylotree_json,
            )
            if missing_tips:
                logger.warning("%s has %d tips without M/N lineage assignment.", subtree_id, len(missing_tips))
            starting_newick, ucld_mean, root_subst = build_calibrated_starting_newick(
                source_tree_path=source_tree_path,
                expected_tip_count=record.total_n_tips,
                min_branch_length=min_branch_length,
                ultrametric_tolerance=ultrametric_tolerance,
                m_tips=m_tips,
                n_tips=n_tips,
                root_age=float(calibration_root_age),
                m_age=float(calibration_m_age),
                n_age=float(calibration_n_age),
                warm_start=_as_bool(calibration_warm_start),
                source_id=subtree_id,
            )
            starting_tree_file.write_text(starting_newick + "\n", encoding="utf-8")
            render_calibrated_beast_xml(
                template_path=calibrated_xml_template,
                taxa_block=taxa_block,
                alignment_block=alignment_block,
                starting_tree_newick=starting_newick,
                calibration_block=build_calibration_block(m_tips, n_tips),
                file_prefix=file_prefix,
                chain_length=chain_length,
                log_every=log_every,
                ucld_mean=ucld_mean,
                popsize=float(calibration_popsize_start),
                root_age=float(calibration_root_age),
                root_stdev=float(calibration_root_stdev),
                m_age=float(calibration_m_age),
                m_stdev=float(calibration_m_stdev),
                n_age=float(calibration_n_age),
                n_stdev=float(calibration_n_stdev),
                destination=xmlfile,
            )
            logger.info(
                "Prepared calibrated XML for %s: M=%d N=%d root_subst=%.6g ucld_mean=%.6g",
                subtree_id, len(m_tips), len(n_tips), root_subst, ucld_mean,
            )
        else:
            starting_newick = build_ultrametric_starting_newick(
                source_tree_path,
                record.total_n_tips,
                min_branch_length=min_branch_length,
                ultrametric_tolerance=ultrametric_tolerance,
            )
            starting_tree_file.write_text(starting_newick + "\n", encoding="utf-8")
            render_beast_xml(
                template_path=xml_template,
                taxa_block=taxa_block,
                alignment_block=alignment_block,
                starting_tree_newick=starting_newick,
                file_prefix=file_prefix,
                chain_length=chain_length,
                log_every=log_every,
                clock_rate=clock_rate,
                destination=xmlfile,
            )

        manifest_rows.append(
            {
                "beast_subtree_id": subtree_id,
                "job_dir": job_dir.as_posix(),
                "xmlfile": xmlfile.as_posix(),
                "starting_tree_file": starting_tree_file.as_posix(),
                "seqfile": seqfile.as_posix(),
                "trees_file": trees_file.as_posix(),
                "mcc_file": mcc_file.as_posix(),
                "file_prefix": file_prefix,
                "expected_tip_count": str(record.total_n_tips),
                "tip_hash": record.tip_hash,
                "seq_id_strategy": strategy,
            }
        )
        logger.info("Prepared %s with %d tips (strategy=%s)", subtree_id, record.total_n_tips, strategy)

    manifest_path = output_dir / "beast_job_manifest.tsv"
    write_status_rows(manifest_rows, BEAST_JOB_MANIFEST_COLUMNS, manifest_path)
    logger.info("Prepared %d BEAST jobs.", len(manifest_rows))
    return 0


def _as_bool(value) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def build_parser():
    parser = argparse.ArgumentParser(description="Prepare BEAST job directories from split subtrees.")
    parser.add_argument("--beast-summary-tsv", required=True, help="Path to beast_subtree_summary.tsv (stage1).")
    parser.add_argument("--beast-tree-dir", required=True, help="Directory containing split subtree trees.")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA with all sequences.")
    parser.add_argument("--xml-template", required=True, help="BEAST XML template path.")
    parser.add_argument("--output-dir", required=True, help="BEAST output directory.")
    parser.add_argument("--seq-id-strategy", required=True, help="exact or prefix_before_first_underscore.")
    parser.add_argument("--chain-length", required=True, type=int, help="MCMC chain length.")
    parser.add_argument("--log-every", required=True, type=int, help="Sampling interval.")
    parser.add_argument("--clock-rate", default=1.0, type=float, help="Fixed strict clock rate (default 1.0 => substitutions/site).")
    parser.add_argument("--calibration-enabled", default="false", help="Whether to render calibrated/UCLD XML.")
    parser.add_argument("--calibrated-xml-template", help="Calibrated/UCLD BEAST XML template path.")
    parser.add_argument("--calibration-id-hap-tsv", help="Sample ID to haplogroup TSV for M/N calibration.")
    parser.add_argument("--calibration-phylotree-json", help="Phylotree lineage JSON for M/N calibration.")
    parser.add_argument("--calibration-warm-start", default="true", help="Whether to warm-start calibrated M/N/root ages.")
    parser.add_argument("--calibration-popsize-start", default=100000.0, type=float, help="Initial constant population size.")
    parser.add_argument("--calibration-root-age", default=180000.0, type=float, help="Root calibration age.")
    parser.add_argument("--calibration-root-stdev", default=20000.0, type=float, help="Root calibration normal-prior stdev.")
    parser.add_argument("--calibration-m-age", default=50000.0, type=float, help="M clade calibration age.")
    parser.add_argument("--calibration-m-stdev", default=5000.0, type=float, help="M clade calibration normal-prior stdev.")
    parser.add_argument("--calibration-n-age", default=57000.0, type=float, help="N clade calibration age.")
    parser.add_argument("--calibration-n-stdev", default=5000.0, type=float, help="N clade calibration normal-prior stdev.")
    parser.add_argument("--min-branch-length", default=1e-8, type=float, help="Minimum positive branch length.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            beast_summary_tsv=args.beast_summary_tsv,
            beast_tree_dir=args.beast_tree_dir,
            input_fasta=args.input_fasta,
            xml_template=args.xml_template,
            output_dir=args.output_dir,
            seq_id_strategy=args.seq_id_strategy,
            chain_length=args.chain_length,
            log_every=args.log_every,
            clock_rate=args.clock_rate,
            calibration_enabled=args.calibration_enabled,
            calibrated_xml_template=args.calibrated_xml_template,
            calibration_id_hap_tsv=args.calibration_id_hap_tsv,
            calibration_phylotree_json=args.calibration_phylotree_json,
            calibration_warm_start=args.calibration_warm_start,
            calibration_popsize_start=args.calibration_popsize_start,
            calibration_root_age=args.calibration_root_age,
            calibration_root_stdev=args.calibration_root_stdev,
            calibration_m_age=args.calibration_m_age,
            calibration_m_stdev=args.calibration_m_stdev,
            calibration_n_age=args.calibration_n_age,
            calibration_n_stdev=args.calibration_n_stdev,
            min_branch_length=args.min_branch_length,
            ultrametric_tolerance=args.ultrametric_tolerance,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

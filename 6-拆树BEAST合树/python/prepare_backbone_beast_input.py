#!/usr/bin/env python3
"""为 backbone（骨架）树准备单个 BEAST 作业。

对应原 PAML 版本的 ``prepare_backbone_paml_input.py``：把骨架树超度量化后作为固定拓扑
起始树，连同骨架序列与渲染好的 BEAST XML 写入作业目录，并记录 backbone 分析清单。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import PipelineError, get_tip_names_from_tree, read_newick_tree, setup_logger
from phylo_ultrastandard_common import BACKBONE_ANALYSIS_MANIFEST_COLUMNS, infer_unique_outgroup
from beast_common import (
    build_alignment_block,
    build_calibrated_starting_newick,
    build_calibration_block,
    build_taxa_block,
    build_ultrametric_starting_newick,
    classify_m_n_tips,
    map_tip_to_source_id,
    read_fasta_record_map,
    render_beast_xml,
    render_calibrated_beast_xml,
    resolve_seq_id_strategy,
    write_status_rows,
    write_subtree_fasta,
)


def run(
    backbone_tree,
    backbone_summary_tsv,
    beast_summary_tsv,
    input_fasta,
    xml_template,
    job_dir,
    analysis_dir,
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
    backbone_tree = Path(backbone_tree).resolve()
    backbone_summary_tsv = Path(backbone_summary_tsv).resolve()
    beast_summary_tsv = Path(beast_summary_tsv).resolve()
    input_fasta = Path(input_fasta).resolve()
    xml_template = Path(xml_template).resolve()
    job_dir = Path(job_dir).resolve()
    analysis_dir = Path(analysis_dir).resolve()
    calibration_enabled = _as_bool(calibration_enabled)
    calibrated_xml_template = Path(calibrated_xml_template).resolve() if calibrated_xml_template else None
    calibration_id_hap_tsv = Path(calibration_id_hap_tsv).resolve() if calibration_id_hap_tsv else None
    calibration_phylotree_json = Path(calibration_phylotree_json).resolve() if calibration_phylotree_json else None

    logger = setup_logger("prepare_backbone_beast_input", analysis_dir / "prepare_backbone_beast_input.log", log_level)
    strategy = resolve_seq_id_strategy(seq_id_strategy)
    outgroup_tip = infer_unique_outgroup(beast_summary_tsv)

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
    if calibration_enabled:
        for required_path, label in (
            (calibrated_xml_template, "calibrated XML template"),
            (calibration_id_hap_tsv, "calibration id_hap_tsv"),
            (calibration_phylotree_json, "calibration phylotree_json"),
        ):
            if required_path is None or not required_path.exists():
                raise PipelineError(f"Missing {label}: {required_path}")

    job_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    file_prefix = "backbone_only"
    starting_tree_file = job_dir / "starting_tree.nwk"
    seqfile = job_dir / "subtree.fasta"
    xmlfile = job_dir / "beast.xml"

    write_subtree_fasta(tree_tip_names, fasta_record_map, tip_to_source_id, seqfile)
    taxa_block = build_taxa_block(tree_tip_names)
    alignment_block = build_alignment_block(tree_tip_names, fasta_record_map, tip_to_source_id)
    if calibration_enabled:
        m_tips, n_tips, missing_tips = classify_m_n_tips(
            tree_tip_names,
            id_hap_tsv=calibration_id_hap_tsv,
            phylotree_json=calibration_phylotree_json,
        )
        if missing_tips:
            logger.warning("backbone has %d tips without M/N lineage assignment.", len(missing_tips))
        starting_newick, ucld_mean, root_subst = build_calibrated_starting_newick(
            source_tree_path=backbone_tree,
            expected_tip_count=len(tree_tip_names),
            min_branch_length=min_branch_length,
            ultrametric_tolerance=ultrametric_tolerance,
            m_tips=m_tips,
            n_tips=n_tips,
            root_age=float(calibration_root_age),
            m_age=float(calibration_m_age),
            n_age=float(calibration_n_age),
            warm_start=_as_bool(calibration_warm_start),
            source_id="backbone",
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
            "Prepared calibrated backbone XML: M=%d N=%d root_subst=%.6g ucld_mean=%.6g",
            len(m_tips), len(n_tips), root_subst, ucld_mean,
        )
    else:
        starting_newick = build_ultrametric_starting_newick(
            backbone_tree, len(tree_tip_names),
            min_branch_length=min_branch_length, ultrametric_tolerance=ultrametric_tolerance,
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

    manifest_rows = []
    for tip_name in tree_tip_names:
        manifest_rows.append({
            "backbone_job_id": file_prefix,
            "tip_name": tip_name,
            "tip_role": "outgroup" if tip_name == outgroup_tip else "backbone",
            "source_id": tip_to_source_id[tip_name],
        })
    write_status_rows(manifest_rows, BACKBONE_ANALYSIS_MANIFEST_COLUMNS, analysis_dir / "backbone_analysis_manifest.tsv")
    logger.info("Prepared backbone-only BEAST job with %d tips.", len(tree_tip_names))
    return 0


def _as_bool(value) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def build_parser():
    parser = argparse.ArgumentParser(description="Prepare a single backbone-only BEAST job.")
    parser.add_argument("--backbone-tree", required=True, help="Path to split backbone_tree.nwk.")
    parser.add_argument("--backbone-summary-tsv", required=True, help="Path to backbone_summary.tsv.")
    parser.add_argument("--beast-summary-tsv", required=True, help="Path to beast_subtree_summary.tsv (stage1).")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA with all sequences.")
    parser.add_argument("--xml-template", required=True, help="BEAST XML template path.")
    parser.add_argument("--job-dir", required=True, help="Directory for the backbone BEAST job.")
    parser.add_argument("--analysis-dir", required=True, help="Directory for backbone analysis outputs.")
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
            backbone_tree=args.backbone_tree,
            backbone_summary_tsv=args.backbone_summary_tsv,
            beast_summary_tsv=args.beast_summary_tsv,
            input_fasta=args.input_fasta,
            xml_template=args.xml_template,
            job_dir=args.job_dir,
            analysis_dir=args.analysis_dir,
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

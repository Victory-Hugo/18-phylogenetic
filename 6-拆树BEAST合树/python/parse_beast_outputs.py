#!/usr/bin/env python3
"""把 BEAST 的 MCC 树解析为合并阶段可用的超度量分析树。

对应原 PAML 版本的 ``parse_paml_outputs.py``：逐子树读取 TreeAnnotator 的 MCC NEXUS，
还原 tip 名并剥离注释得到 Newick，再做超度量收敛与一致性校验，写出
``analysis_trees/{subtree_id}.nwk``（单位 substitutions/site，供 stage3 合并）。
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

from Bio import Phylo

from phylo_merge_common import (
    load_paml_summary,
    validate_tip_hash,
    validate_tree_against_expected_tip_set,
    write_tree_file,
)
from phylo_split_common import PipelineError, setup_logger
from phylo_ultrastandard_common import extend_terminal_branches_to_height, validate_ultrametric
from beast_common import (
    BEAST_NODE_HPD_COLUMNS,
    BEAST_PARSE_STATUS_COLUMNS,
    build_node_hpd_rows,
    load_rows,
    parse_mcc_node_annotations,
    parse_mcc_to_newick,
    write_status_rows,
)


def run(
    beast_summary_tsv,
    manifest_path,
    output_dir,
    min_branch_length=1e-8,
    ultrametric_tolerance=1e-8,
    outgroup_tip_name=None,
    log_level="INFO",
):
    beast_summary_tsv = Path(beast_summary_tsv).resolve()
    manifest_path = Path(manifest_path).resolve()
    output_dir = Path(output_dir).resolve()
    min_branch_length = float(min_branch_length)
    ultrametric_tolerance = float(ultrametric_tolerance)

    logger = setup_logger("parse_beast_outputs", output_dir / "parse_beast_outputs.log", log_level)
    logger.info("Parsing BEAST MCC outputs from manifest: %s", manifest_path)

    summary_records = load_paml_summary(beast_summary_tsv)
    if outgroup_tip_name in (None, "", "null"):
        unique_outgroups = sorted({record.outgroup_tip for record in summary_records})
        if len(unique_outgroups) != 1:
            raise PipelineError(
                f"Could not infer a unique outgroup tip from {beast_summary_tsv}: {unique_outgroups}"
            )
        outgroup_tip_name = unique_outgroups[0]
    records = {record.baseml_subtree_id: record for record in summary_records}

    manifest_rows = load_rows(manifest_path)
    if not manifest_rows:
        raise PipelineError(f"BEAST job manifest is empty: {manifest_path}")

    analysis_tree_dir = output_dir / "analysis_trees"
    analysis_tree_dir.mkdir(parents=True, exist_ok=True)

    status_rows = []
    failures = []
    node_hpd_rows = []

    for row in manifest_rows:
        subtree_id = row["beast_subtree_id"]
        record = records.get(subtree_id)
        if record is None:
            raise PipelineError(f"Manifest subtree id not present in summary: {subtree_id}")

        mcc_file = Path(row["mcc_file"]).resolve()
        analysis_tree_file = analysis_tree_dir / f"{subtree_id}.nwk"

        try:
            newick = parse_mcc_to_newick(mcc_file)
            mcc_tree = Phylo.read(io.StringIO(newick + "\n"), "newick")

            # 超度量收敛（BEAST MCC 已近似超度量，这里抹平数值容差）
            ultrametric_tree, _, _ = extend_terminal_branches_to_height(
                mcc_tree, min_branch_length=min_branch_length, tolerance=ultrametric_tolerance,
            )
            validate_ultrametric(ultrametric_tree, ultrametric_tolerance, f"analysis tree {subtree_id}")
            validate_tree_against_expected_tip_set(
                ultrametric_tree, record.total_tip_names, f"analysis tree {subtree_id}", outgroup_tip_name,
            )
            validate_tip_hash(record.total_tip_names, record.tip_hash, f"analysis tree {subtree_id}")

            actual_tip_count = len(ultrametric_tree.get_terminals())
            if actual_tip_count != record.total_n_tips:
                raise PipelineError(
                    f"{subtree_id}: expected {record.total_n_tips} tips, got {actual_tip_count}"
                )

            write_tree_file(ultrametric_tree, analysis_tree_file)
            # 抽取该子树 MCC 各内部节点的高度与 95% HPD（substitutions/site），
            # 供 stage5 按 clade 把 HPD 线性传播到最终年代树。
            node_hpd_rows.extend(build_node_hpd_rows(subtree_id, parse_mcc_node_annotations(mcc_file)))
            status_rows.append({
                "beast_subtree_id": subtree_id,
                "mcc_file": mcc_file.as_posix(),
                "analysis_tree_file": analysis_tree_file.as_posix(),
                "status": "success",
                "detail": "Parsed MCC into ultrametric analysis tree",
            })
            logger.info("Parsed %s successfully", subtree_id)
        except Exception as exc:
            detail = str(exc)
            failures.append((subtree_id, detail))
            status_rows.append({
                "beast_subtree_id": subtree_id,
                "mcc_file": mcc_file.as_posix(),
                "analysis_tree_file": analysis_tree_file.as_posix(),
                "status": "failed",
                "detail": detail,
            })
            logger.error("Failed to parse %s: %s", subtree_id, detail)

    status_rows.sort(key=lambda item: item["beast_subtree_id"])
    write_status_rows(status_rows, BEAST_PARSE_STATUS_COLUMNS, output_dir / "beast_parse_status.tsv")
    write_status_rows(node_hpd_rows, BEAST_NODE_HPD_COLUMNS, output_dir / "beast_node_hpd.tsv")

    if failures:
        preview = "; ".join(f"{sid}: {detail}" for sid, detail in failures[:5])
        raise PipelineError(f"{len(failures)} BEAST output(s) failed to parse: {preview}")

    logger.info("Parsed %d BEAST outputs successfully.", len(status_rows))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Parse BEAST MCC outputs into merge-ready ultrametric trees.")
    parser.add_argument("--beast-summary-tsv", required=True, help="Path to beast_subtree_summary.tsv (stage1).")
    parser.add_argument("--manifest", required=True, help="Path to beast_job_manifest.tsv.")
    parser.add_argument("--output-dir", required=True, help="BEAST output directory.")
    parser.add_argument("--min-branch-length", default=1e-8, type=float, help="Minimum positive branch length.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--outgroup-tip-name", default=None, help="Global outgroup tip name.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            beast_summary_tsv=args.beast_summary_tsv,
            manifest_path=args.manifest,
            output_dir=args.output_dir,
            min_branch_length=args.min_branch_length,
            ultrametric_tolerance=args.ultrametric_tolerance,
            outgroup_tip_name=args.outgroup_tip_name,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

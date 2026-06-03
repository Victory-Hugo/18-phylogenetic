#!/usr/bin/env python3
"""把 backbone BEAST 的 MCC 树解析为规范化的超度量骨架树。

对应原 PAML 版本的 ``parse_backbone_paml_output.py``：产出
``backbone_raw_analysis_tree.nwk``、``backbone_ultrametric_tree.nwk`` 与超度量报表，
其中 ``backbone_ultrametric_tree.nwk`` 供 stage3 合并作为骨架。
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

from Bio import Phylo

from phylo_merge_common import validate_tree_against_expected_tip_set
from phylo_split_common import PipelineError, get_tip_names_from_tree, read_newick_tree, setup_logger
from phylo_ultrastandard_common import (
    BACKBONE_ULTRAMETRIC_REPORT_COLUMNS,
    build_root_to_tip_row,
    extend_terminal_branches_to_height,
    infer_unique_outgroup,
    validate_ultrametric,
    write_simple_tree,
)
from beast_common import (
    BEAST_NODE_HPD_COLUMNS,
    build_node_hpd_rows,
    parse_mcc_node_annotations,
    parse_mcc_to_newick,
    write_status_rows,
)


def run(
    backbone_tree,
    beast_summary_tsv,
    job_dir,
    analysis_dir,
    min_branch_length=1e-8,
    ultrametric_tolerance=1e-8,
    log_level="INFO",
):
    backbone_tree = Path(backbone_tree).resolve()
    beast_summary_tsv = Path(beast_summary_tsv).resolve()
    job_dir = Path(job_dir).resolve()
    analysis_dir = Path(analysis_dir).resolve()
    min_branch_length = float(min_branch_length)
    ultrametric_tolerance = float(ultrametric_tolerance)

    logger = setup_logger("parse_backbone_beast_output", analysis_dir / "parse_backbone_beast_output.log", log_level)
    mcc_file = job_dir / "backbone_only.mcc.tree"
    if not mcc_file.exists():
        raise PipelineError(f"Backbone MCC result file not found: {mcc_file}")

    analysis_dir.mkdir(parents=True, exist_ok=True)
    outgroup_tip = infer_unique_outgroup(beast_summary_tsv)
    expected_tip_names = get_tip_names_from_tree(read_newick_tree(backbone_tree))

    newick = parse_mcc_to_newick(mcc_file)
    raw_tree = Phylo.read(io.StringIO(newick + "\n"), "newick")
    raw_tree_path = analysis_dir / "backbone_raw_analysis_tree.nwk"
    write_simple_tree(raw_tree_path, raw_tree)

    raw_tree = read_newick_tree(raw_tree_path)
    validate_tree_against_expected_tip_set(raw_tree, expected_tip_names, "backbone_raw_analysis_tree", outgroup_tip)

    ultrametric_tree, _, normalization_rows = extend_terminal_branches_to_height(
        raw_tree, min_branch_length=min_branch_length, tolerance=ultrametric_tolerance,
    )
    validate_tree_against_expected_tip_set(
        ultrametric_tree, expected_tip_names, "backbone_ultrametric_tree", outgroup_tip,
    )

    ultrametric_tree_path = analysis_dir / "backbone_ultrametric_tree.nwk"
    validate_ultrametric(ultrametric_tree, ultrametric_tolerance, "backbone_ultrametric_tree")
    write_simple_tree(ultrametric_tree_path, ultrametric_tree)

    report_rows = [
        normalization_rows["pre"],
        normalization_rows["post"],
        build_root_to_tip_row("backbone_ultrametric_tree", ultrametric_tree, ultrametric_tolerance),
    ]
    write_status_rows(report_rows, BACKBONE_ULTRAMETRIC_REPORT_COLUMNS, analysis_dir / "backbone_ultrametric_report.tsv")

    # 抽取 backbone MCC 各内部节点的高度与 95% HPD（substitutions/site），单列一份，
    # 供 stage5 与各 target 子树的 beast_node_hpd.tsv 合并构建 HPD 注册表。
    backbone_node_hpd_rows = build_node_hpd_rows("backbone", parse_mcc_node_annotations(mcc_file))
    write_status_rows(backbone_node_hpd_rows, BEAST_NODE_HPD_COLUMNS, analysis_dir / "backbone_node_hpd.tsv")

    logger.info("Parsed backbone BEAST output and wrote ultrametric backbone tree.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Parse backbone BEAST output into an ultrametric backbone tree.")
    parser.add_argument("--backbone-tree", required=True, help="Path to split backbone_tree.nwk.")
    parser.add_argument("--beast-summary-tsv", required=True, help="Path to beast_subtree_summary.tsv (stage1).")
    parser.add_argument("--job-dir", required=True, help="Backbone BEAST job directory.")
    parser.add_argument("--analysis-dir", required=True, help="Backbone analysis output directory.")
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
            beast_summary_tsv=args.beast_summary_tsv,
            job_dir=args.job_dir,
            analysis_dir=args.analysis_dir,
            min_branch_length=args.min_branch_length,
            ultrametric_tolerance=args.ultrametric_tolerance,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

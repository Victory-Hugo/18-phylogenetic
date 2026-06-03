#!/usr/bin/env python3
"""生成可供合并阶段使用的伪 BEAST 输出（保持拓扑、超度量枝长）。

对应原 PAML 版本的 ``fake_paml_outputs.py``：在不真正运行 BEAST MCMC 的情况下，
直接由 stage1 子树拓扑随机化枝长并超度量化，快速产出 ``analysis_trees/{id}.nwk`` 与
backbone 超度量树，用于联调 stage1→3→4→5 全链路。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_merge_common import (
    load_paml_summary,
    randomize_tree_branch_lengths,
    read_newick_tree,
    stable_hash_int,
    validate_rooted_binary_tree,
    validate_tree_against_expected_tip_set,
    write_rows,
    write_tree_file,
)
from phylo_split_common import PipelineError, clone_tree, setup_logger
from phylo_ultrastandard_common import (
    BACKBONE_ULTRAMETRIC_REPORT_COLUMNS,
    build_root_to_tip_row,
    extend_terminal_branches_to_height,
    validate_ultrametric,
)
from beast_common import (
    BEAST_JOB_MANIFEST_COLUMNS,
    BEAST_PARSE_STATUS_COLUMNS,
    BEAST_RUN_STATUS_COLUMNS,
)


FAKE_BEAST_MANIFEST_COLUMNS = [
    "beast_subtree_id",
    "source_tree_file",
    "analysis_tree_file",
    "tree_seed",
    "n_edges_randomized",
    "min_multiplier",
    "max_multiplier",
    "mean_multiplier",
    "pre_ultrametric_delta",
    "post_ultrametric_delta",
    "target_height",
]


def _infer_outgroup(records, summary_tsv: Path) -> str:
    unique_outgroups = sorted({record.outgroup_tip for record in records})
    if len(unique_outgroups) != 1:
        raise PipelineError(f"Could not infer a unique outgroup tip from {summary_tsv}: {unique_outgroups}")
    return unique_outgroups[0]


def _materialize_fake_ultrametric_tree(source_tree, tree_seed, branch_length_sigma, min_branch_length, ultrametric_tolerance):
    randomized_tree = clone_tree(source_tree)
    n_edges, min_multiplier, max_multiplier, mean_multiplier = randomize_tree_branch_lengths(
        randomized_tree, sigma=branch_length_sigma, seed=tree_seed, min_branch_length=min_branch_length,
    )
    pre_row = build_root_to_tip_row("pre_ultrametric", randomized_tree, ultrametric_tolerance)
    ultrametric_tree, target_height, normalization_rows = extend_terminal_branches_to_height(
        randomized_tree, min_branch_length=min_branch_length, tolerance=ultrametric_tolerance,
    )
    post_row = build_root_to_tip_row("post_ultrametric", ultrametric_tree, ultrametric_tolerance)
    validate_ultrametric(ultrametric_tree, ultrametric_tolerance, "fake_ultrametric_tree")
    return {
        "tree": ultrametric_tree,
        "tree_seed": tree_seed,
        "n_edges_randomized": n_edges,
        "min_multiplier": min_multiplier,
        "max_multiplier": max_multiplier,
        "mean_multiplier": mean_multiplier,
        "target_height": target_height,
        "pre_row": normalization_rows.get("pre", pre_row),
        "post_row": normalization_rows.get("post", post_row),
    }


def run(
    beast_summary_tsv,
    beast_tree_dir,
    backbone_tree,
    output_dir,
    backbone_analysis_dir,
    branch_length_model,
    branch_length_sigma,
    random_seed,
    min_branch_length,
    ultrametric_tolerance,
    backbone_only_enabled=True,
    log_level="INFO",
):
    if branch_length_model != "lognormal_then_extend_terminal":
        raise PipelineError(f"Unsupported beast.fake_branch_length_model: {branch_length_model}")

    beast_summary_tsv = Path(beast_summary_tsv).resolve()
    beast_tree_dir = Path(beast_tree_dir).resolve()
    backbone_tree = Path(backbone_tree).resolve()
    output_dir = Path(output_dir).resolve()
    backbone_analysis_dir = Path(backbone_analysis_dir).resolve()
    branch_length_sigma = float(branch_length_sigma)
    random_seed = int(random_seed)
    min_branch_length = float(min_branch_length)
    ultrametric_tolerance = float(ultrametric_tolerance)

    logger = setup_logger("fake_beast_outputs", output_dir / "fake_beast_outputs.log", log_level)
    analysis_tree_dir = output_dir / "analysis_trees"
    analysis_tree_dir.mkdir(parents=True, exist_ok=True)
    backbone_analysis_dir.mkdir(parents=True, exist_ok=True)

    records = load_paml_summary(beast_summary_tsv)
    outgroup_tip_name = _infer_outgroup(records, beast_summary_tsv)

    manifest_rows = []
    run_status_rows = []
    parse_status_rows = []
    fake_rows = []

    for record in records:
        subtree_id = record.baseml_subtree_id
        source_tree_path = beast_tree_dir / Path(record.baseml_tree_file).name
        if not source_tree_path.exists():
            raise PipelineError(f"Subtree file not found: {source_tree_path}")
        source_tree = read_newick_tree(source_tree_path)
        validate_tree_against_expected_tip_set(source_tree, record.total_tip_names, subtree_id, outgroup_tip_name)
        validate_rooted_binary_tree(source_tree, outgroup_tip_name)

        tree_seed = random_seed + (stable_hash_int(subtree_id) % (2**31))
        fake_result = _materialize_fake_ultrametric_tree(
            source_tree=source_tree, tree_seed=tree_seed, branch_length_sigma=branch_length_sigma,
            min_branch_length=min_branch_length, ultrametric_tolerance=ultrametric_tolerance,
        )
        analysis_tree_file = analysis_tree_dir / f"{subtree_id}.nwk"
        write_tree_file(fake_result["tree"], analysis_tree_file)

        manifest_rows.append({
            "beast_subtree_id": subtree_id, "job_dir": "", "xmlfile": "",
            "starting_tree_file": source_tree_path.as_posix(), "seqfile": "",
            "trees_file": "", "mcc_file": "", "file_prefix": subtree_id,
            "expected_tip_count": str(record.total_n_tips), "tip_hash": record.tip_hash, "seq_id_strategy": "fake",
        })
        run_status_rows.append({
            "beast_subtree_id": subtree_id, "job_dir": "", "trees_file": "", "mcc_file": "",
            "return_code": "0", "status": "success", "detail": "Generated by fake_beast_outputs.py",
        })
        parse_status_rows.append({
            "beast_subtree_id": subtree_id, "mcc_file": "", "analysis_tree_file": analysis_tree_file.as_posix(),
            "status": "success", "detail": "Ultrametric fake analysis tree generated directly from split subtree topology.",
        })
        fake_rows.append({
            "beast_subtree_id": subtree_id,
            "source_tree_file": source_tree_path.as_posix(),
            "analysis_tree_file": analysis_tree_file.as_posix(),
            "tree_seed": str(fake_result["tree_seed"]),
            "n_edges_randomized": str(fake_result["n_edges_randomized"]),
            "min_multiplier": f"{fake_result['min_multiplier']:.12g}",
            "max_multiplier": f"{fake_result['max_multiplier']:.12g}",
            "mean_multiplier": f"{fake_result['mean_multiplier']:.12g}",
            "pre_ultrametric_delta": fake_result["pre_row"]["delta"],
            "post_ultrametric_delta": fake_result["post_row"]["delta"],
            "target_height": f"{float(fake_result['target_height']):.12g}",
        })
        logger.info("Generated %s seed=%d n_edges=%d", subtree_id, fake_result["tree_seed"], fake_result["n_edges_randomized"])

    write_rows(manifest_rows, BEAST_JOB_MANIFEST_COLUMNS, output_dir / "beast_job_manifest.tsv")
    write_rows(run_status_rows, BEAST_RUN_STATUS_COLUMNS, output_dir / "beast_run_status.tsv")
    write_rows(parse_status_rows, BEAST_PARSE_STATUS_COLUMNS, output_dir / "beast_parse_status.tsv")
    write_rows(fake_rows, FAKE_BEAST_MANIFEST_COLUMNS, output_dir / "fake_beast_manifest.tsv")

    if backbone_only_enabled:
        source_backbone_tree = read_newick_tree(backbone_tree)
        validate_rooted_binary_tree(source_backbone_tree, outgroup_tip_name)
        backbone_seed = random_seed + 1_000_000_007
        fake_backbone = _materialize_fake_ultrametric_tree(
            source_tree=source_backbone_tree, tree_seed=backbone_seed, branch_length_sigma=branch_length_sigma,
            min_branch_length=min_branch_length, ultrametric_tolerance=ultrametric_tolerance,
        )
        raw_backbone_path = backbone_analysis_dir / "backbone_raw_analysis_tree.nwk"
        ultrametric_backbone_path = backbone_analysis_dir / "backbone_ultrametric_tree.nwk"
        write_tree_file(fake_backbone["tree"], raw_backbone_path)
        write_tree_file(fake_backbone["tree"], ultrametric_backbone_path)
        backbone_report_rows = [
            fake_backbone["pre_row"],
            fake_backbone["post_row"],
            build_root_to_tip_row("backbone_ultrametric_tree", fake_backbone["tree"], ultrametric_tolerance),
        ]
        write_rows(backbone_report_rows, BACKBONE_ULTRAMETRIC_REPORT_COLUMNS, backbone_analysis_dir / "backbone_ultrametric_report.tsv")
        logger.info("Generated fake backbone ultrametric tree seed=%d", backbone_seed)

    logger.info("Generated %d fake ultrametric analysis tree(s).", len(records))
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Generate fake ultrametric BEAST outputs for merge-stage development.")
    parser.add_argument("--beast-summary-tsv", required=True, help="Path to beast_subtree_summary.tsv (stage1).")
    parser.add_argument("--beast-tree-dir", required=True, help="Directory containing split subtree trees.")
    parser.add_argument("--backbone-tree", required=True, help="Path to split backbone_tree.nwk.")
    parser.add_argument("--output-dir", required=True, help="BEAST output directory.")
    parser.add_argument("--backbone-analysis-dir", required=True, help="Backbone analysis output directory.")
    parser.add_argument("--branch-length-model", required=True, help="Fake branch length model.")
    parser.add_argument("--branch-length-sigma", default=0.25, type=float, help="Sigma for branch-length perturbation.")
    parser.add_argument("--random-seed", default=42, type=int, help="Base random seed.")
    parser.add_argument("--min-branch-length", default=1e-8, type=float, help="Minimum positive branch length.")
    parser.add_argument("--ultrametric-tolerance", default=1e-8, type=float, help="Ultrametric tolerance.")
    parser.add_argument("--backbone-only-enabled", default="true", help="Whether to generate backbone-only fake output.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            beast_summary_tsv=args.beast_summary_tsv,
            beast_tree_dir=args.beast_tree_dir,
            backbone_tree=args.backbone_tree,
            output_dir=args.output_dir,
            backbone_analysis_dir=args.backbone_analysis_dir,
            branch_length_model=args.branch_length_model,
            branch_length_sigma=args.branch_length_sigma,
            random_seed=args.random_seed,
            min_branch_length=args.min_branch_length,
            ultrametric_tolerance=args.ultrametric_tolerance,
            backbone_only_enabled=str(args.backbone_only_enabled).lower() == "true",
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

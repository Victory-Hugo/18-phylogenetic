#!/usr/bin/env python3
"""根据 LSD2 输出和 XML 模板生成 BEAST v10 Thorney XML。"""

from __future__ import annotations

import argparse
import logging
import shutil
from pathlib import Path
from typing import Any

from beast_pipeline_utils import (
    PipelineError,
    parse_calibration_nodes,
    read_membership,
    read_tree,
    read_validation_table,
    render_beast_xml,
    tree_tip_names,
    tree_to_newick_text,
    validate_calibration_nesting,
    validate_starting_tree,
)


log = logging.getLogger(__name__)


def write_report(path: Path, values: dict[str, object]) -> None:
    """写 XML 构建报告。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# BEAST XML 构建报告",
        "",
        f"- ML tree tips: {values['tree_tip_count']}",
        f"- 用于 XML 的 taxa: {values['xml_taxa_count']}",
        f"- FASTA 中未进入树的序列数: {values['fasta_not_in_tree']}",
    ]
    for name, count in values["node_counts"].items():  # type: ignore[union-attr]
        lines.append(f"- {name} 校准 tips: {count}")
    lines += [
        f"- startingTree root height: {values['starting_tree_root_height']:.6f}",
        f"- XML 模板: `{values['template_path']}`",
        f"- XML: `{values['xml_path']}`",
        "",
        "XML 由 conf 中的模板渲染；每个校准节点用 normalPrior 软约束其 TMRCA。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_xml_text(
    *,
    template_path: Path,
    ml_tree_path: Path,
    starting_tree_path: Path,
    membership_path: Path,
    nodes: list[dict[str, Any]],
    scale: float,
    chain_length: int,
    log_every: int,
    tree_log_every: int,
    clock_rate: float,
    skygrid_points: int,
    skygrid_init_logpop: float,
    output_prefix: str,
) -> tuple[str, dict[str, object]]:
    """准备模板变量并渲染 XML。"""
    ml_tree_obj = read_tree(ml_tree_path)
    taxa = tree_tip_names(ml_tree_obj)
    taxa_set = set(taxa)
    tree_validation = validate_starting_tree(ml_tree_path, starting_tree_path, taxa_set)
    membership, _hap_map = read_membership(membership_path)
    data_tree = tree_to_newick_text(ml_tree_path)
    starting_tree = tree_to_newick_text(starting_tree_path)

    non_root = [n for n in nodes if not n["is_root"]]
    root_node = next(n for n in nodes if n["is_root"])
    for node in nodes:
        if node["age"] is None or node["stdev"] is None:
            raise PipelineError(f"校准节点 {node['name']} 缺少 age 或 stdev（XML 步骤两者都需要）")

    # 成员表 ↔ yaml 一致性。
    yaml_names = {n["name"] for n in nodes}
    for tsv_node in membership:
        if tsv_node not in yaml_names:
            raise PipelineError(f"成员表中的节点 {tsv_node} 未在 yaml calibration.nodes 定义，请保持一致")

    calibration_taxa: dict[str, list[str]] = {}
    for node in non_root:
        name = node["name"]
        if name not in membership:
            raise PipelineError(f"yaml 定义的校准节点 {name} 在成员表中无记录，请保持一致")
        members = sorted(set(membership[name]) & taxa_set)
        if not members:
            raise PipelineError(f"校准节点 {name} 在树中没有成员（成员表与树不匹配）")
        calibration_taxa[name] = members
    validate_calibration_nesting(nodes, {name: set(members) for name, members in calibration_taxa.items()})

    xml_text = render_beast_xml(
        template_path=template_path,
        taxa=taxa,
        starting_tree=starting_tree,
        data_tree=data_tree,
        calibration_taxa=calibration_taxa,
        nodes=nodes,
        scale=scale,
        clock_rate=clock_rate,
        chain_length=chain_length,
        log_every=log_every,
        tree_log_every=tree_log_every,
        output_prefix=output_prefix,
        # cutOff 略大于 root 年龄，避免最老 coalescent interval 恰好压在网格边界报错。
        skygrid_cutoff=root_node["age"] * 1.1,
        skygrid_init_logpop=skygrid_init_logpop,
        skygrid_points=skygrid_points,
    )
    metadata = {
        "tree_tip_count": len(taxa),
        "xml_taxa_count": len(taxa),
        "node_counts": {name: len(members) for name, members in calibration_taxa.items()},
        "starting_tree_root_height": tree_validation["root_height"],
    }
    return xml_text, metadata


def run(
    *,
    ml_tree: str,
    starting_tree: str,
    calibration_membership: str,
    calibration_nodes: str,
    input_validation: str,
    template_xml: str,
    output_table: str,
    output_report: str,
    temp_dir: str,
    output_prefix: str,
    scale: float,
    chain_length: int,
    log_every: int,
    tree_log_every: int,
    clock_rate: float,
    skygrid_points: int,
    skygrid_init_logpop: float,
    overwrite: str = "false",
) -> int:
    """根据 Step 1 输出生成正式 XML 和短 chain 校验 XML。"""
    ml_tree_path = Path(ml_tree)
    starting_tree_path = Path(starting_tree)
    membership_path = Path(calibration_membership)
    input_validation_path = Path(input_validation)
    template_path = Path(template_xml)
    output_table_path = Path(output_table)
    output_report_path = Path(output_report)
    temp_path = Path(temp_dir)
    for required in [ml_tree_path, starting_tree_path, membership_path, input_validation_path, template_path]:
        if not required.exists():
            raise PipelineError(f"缺少输入文件: {required}")
    if overwrite.lower() == "true" and temp_path.exists():
        shutil.rmtree(temp_path)
    output_table_path.mkdir(parents=True, exist_ok=True)
    output_report_path.mkdir(parents=True, exist_ok=True)
    temp_path.mkdir(parents=True, exist_ok=True)

    nodes = parse_calibration_nodes(calibration_nodes)

    xml_text, metadata = build_xml_text(
        template_path=template_path,
        ml_tree_path=ml_tree_path,
        starting_tree_path=starting_tree_path,
        membership_path=membership_path,
        nodes=nodes,
        scale=scale,
        chain_length=chain_length,
        log_every=log_every,
        tree_log_every=tree_log_every,
        clock_rate=clock_rate,
        skygrid_points=skygrid_points,
        skygrid_init_logpop=skygrid_init_logpop,
        output_prefix=output_prefix,
    )
    xml_path = output_table_path / f"{output_prefix}.xml"
    xml_path.write_text(xml_text, encoding="utf-8")

    validate_xml_text, _ = build_xml_text(
        template_path=template_path,
        ml_tree_path=ml_tree_path,
        starting_tree_path=starting_tree_path,
        membership_path=membership_path,
        nodes=nodes,
        scale=scale,
        chain_length=200,
        log_every=100,
        tree_log_every=200,
        clock_rate=clock_rate,
        skygrid_points=skygrid_points,
        skygrid_init_logpop=skygrid_init_logpop,
        output_prefix=f"{output_prefix}.validate",
    )
    validate_xml_path = temp_path / "beast_validate" / f"{output_prefix}.validate.xml"
    validate_xml_path.parent.mkdir(parents=True, exist_ok=True)
    validate_xml_path.write_text(validate_xml_text, encoding="utf-8")

    validation_values = read_validation_table(input_validation_path)
    write_report(
        output_report_path / "build_beast_xml_report.md",
        {
            **metadata,
            "fasta_not_in_tree": validation_values.get("fasta_not_in_tree", "NA"),
            "template_path": template_path,
            "xml_path": xml_path,
        },
    )
    log.info("XML 已生成: %s", xml_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description="根据 LSD2 输出生成 BEAST v10 Thorney XML")
    parser.add_argument("--ml-tree", required=True)
    parser.add_argument("--starting-tree", required=True)
    parser.add_argument("--calibration-membership", required=True)
    parser.add_argument("--calibration-nodes", required=True)
    parser.add_argument("--input-validation", required=True)
    parser.add_argument("--template-xml", required=True)
    parser.add_argument("--output-table", required=True)
    parser.add_argument("--output-report", required=True)
    parser.add_argument("--temp-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--scale", type=float, required=True)
    parser.add_argument("--chain-length", type=int, required=True)
    parser.add_argument("--log-every", type=int, required=True)
    parser.add_argument("--tree-log-every", type=int, required=True)
    parser.add_argument("--clock-rate", type=float, required=True)
    parser.add_argument("--skygrid-points", type=int, required=True)
    parser.add_argument("--skygrid-init-logpop", type=float, required=True)
    parser.add_argument("--overwrite", default="false")
    return parser


def main(argv: list[str] | None = None) -> int:
    """命令行入口。"""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(**vars(args))
    except PipelineError as exc:
        log.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

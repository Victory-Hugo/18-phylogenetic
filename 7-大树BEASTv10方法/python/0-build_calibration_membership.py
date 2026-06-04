#!/usr/bin/env python3
"""从 ID-Hap 映射与 phylotree JSON 生成校准成员表（多遗传标记通用预处理）。

把「哪些样本属于哪个校准 clade」一次性物化成 long 格式 tsv，供下游
1-lsd2 / 2-beast_xml 使用。此后主流程不再读 phylotree JSON，从而与具体
遗传标记解耦：其它标记的用户只需用自己的方式产出同格式的成员表即可。
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
from pathlib import Path

from beast_pipeline_utils import (
    PipelineError,
    haplogroup_in_node,
    parse_calibration_nodes,
    read_id_hap,
    write_membership,
)


log = logging.getLogger(__name__)


def write_report(path: Path, values: dict[str, object]) -> None:
    """写预处理步骤报告。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# 校准成员表构建报告",
        "",
        f"- ID-Hap 样本数: {values['id_hap_count']}",
        f"- 校准节点数（不含 root）: {values['node_count']}",
        f"- 成员表总行数: {values['row_count']}",
        "",
        "## 各节点成员数",
        "",
    ]
    for node, count in values["counts"].items():  # type: ignore[union-attr]
        lines.append(f"- {node}: {count}")
    lines += [
        "",
        "成员关系基于 phylotree lineage 判定（haplogroup 的 lineage 是否包含目标节点）。",
        "嵌套节点的样本会在多个节点下各出现一行；下游会再与树 tip 取交集。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(
    *,
    id_hap: str,
    phylotree_json: str,
    calibration_nodes: str,
    output_table: str,
    output_report: str,
    temp_dir: str,
    overwrite: str = "false",
) -> int:
    """生成校准成员表 calibration_membership.tsv。"""
    id_hap_path = Path(id_hap)
    phylotree_path = Path(phylotree_json)
    output_table_path = Path(output_table)
    output_report_path = Path(output_report)
    temp_path = Path(temp_dir)
    for required in [id_hap_path, phylotree_path]:
        if not required.exists():
            raise PipelineError(f"缺少输入文件: {required}")
    if overwrite.lower() == "true" and temp_path.exists():
        shutil.rmtree(temp_path)
    output_table_path.mkdir(parents=True, exist_ok=True)
    output_report_path.mkdir(parents=True, exist_ok=True)
    temp_path.mkdir(parents=True, exist_ok=True)

    nodes = parse_calibration_nodes(calibration_nodes)
    target_nodes = [n["name"] for n in nodes if not n["is_root"]]
    if not target_nodes:
        raise PipelineError("除 root 外没有可生成成员的校准节点")

    id_hap_map = read_id_hap(id_hap_path)
    phylotree = json.loads(phylotree_path.read_text(encoding="utf-8"))
    haplogroups = phylotree.get("haplogroups", {})
    if not haplogroups:
        raise PipelineError(f"phylotree JSON 缺少 haplogroups: {phylotree_path}")

    rows: list[tuple[str, str, str]] = []
    counts: dict[str, int] = {}
    for node in target_nodes:
        members = sorted(
            sample_id
            for sample_id, hap in id_hap_map.items()
            if haplogroup_in_node(hap, node, haplogroups)
        )
        if not members:
            raise PipelineError(f"校准节点 {node} 在 ID-Hap 中没有任何成员，请检查节点名是否为有效 haplogroup")
        counts[node] = len(members)
        for sample_id in members:
            rows.append((node, sample_id, id_hap_map[sample_id]))

    membership_path = output_table_path / "calibration_membership.tsv"
    write_membership(membership_path, rows)
    write_report(
        output_report_path / "calibration_membership_report.md",
        {
            "id_hap_count": len(id_hap_map),
            "node_count": len(target_nodes),
            "row_count": len(rows),
            "counts": counts,
        },
    )
    log.info("成员表已生成: %s (%d 行, %d 个节点)", membership_path, len(rows), len(target_nodes))
    return 0


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description="从 ID-Hap + phylotree JSON 生成校准成员表")
    parser.add_argument("--id-hap", required=True)
    parser.add_argument("--phylotree-json", required=True)
    parser.add_argument("--calibration-nodes", required=True)
    parser.add_argument("--output-table", required=True)
    parser.add_argument("--output-report", required=True)
    parser.add_argument("--temp-dir", required=True)
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

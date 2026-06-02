"""
s08_prepare_ml_inputs.py
准备 IQ-TREE 3 ML 分析的输入文件：净化名称，提取叶节点序列。

IQ-TREE -te 模式要求 FASTA 序列与树的叶节点一一对应。
本脚本从 haplogroup_graph.json 中识别叶节点（无子节点），
仅输出这些叶节点的净化序列 FASTA。

输出两个文件：
1. ancestors_sanitized_leaves.fasta  —— 仅叶节点（~3019 条），供 IQ-TREE 使用
2. name_sanitize_map.tsv             —— 全部节点净化映射（5439 行），供 s10 反映射

用法（CLI）：
    python s08_prepare_ml_inputs.py \
        --fasta      data/ancestral_sequences.fasta \
        --hap-graph  output/stage_1/intermediate/haplogroup_graph.json \
        --output-fasta  output/stage_2/intermediate/ml/ancestors_sanitized_leaves.fasta \
        --output-map    output/stage_2/intermediate/ml/name_sanitize_map.tsv

用法（import）：
    from s08_prepare_ml_inputs import run
    run(fasta_in=..., hap_graph=..., out_fasta=..., out_map=...)
"""

import argparse
import json
import logging
import re
from pathlib import Path

from Bio import SeqIO

log = logging.getLogger(__name__)

_SPECIAL = re.compile(r"[()[\]!',\s]")
_MULTI_UNDERSCORE = re.compile(r"_+")


def sanitize_name(name: str) -> str:
    """将单倍群名中的 Newick 禁用字符替换为 '_'，去除首尾和连续下划线。"""
    s = _SPECIAL.sub("_", name)
    s = _MULTI_UNDERSCORE.sub("_", s)
    return s.strip("_")


def run(fasta_in: str, hap_graph: str, out_fasta: str, out_map: str) -> None:
    fasta_in = Path(fasta_in)
    hap_graph = Path(hap_graph)
    out_fasta = Path(out_fasta)
    out_map = Path(out_map)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    # 读取图，识别叶节点（无子节点的单倍群）
    with open(hap_graph) as f:
        graph = json.load(f)
    leaf_set = {k for k, v in graph.items() if not v["children"]}
    log.info("图节点：%d（叶节点 %d，内部节点 %d）", len(graph), len(leaf_set), len(graph) - len(leaf_set))

    records = list(SeqIO.parse(fasta_in, "fasta"))
    log.info("读入序列：%d 条", len(records))

    # 生成净化映射（全部节点，用于 s10 反映射），检查碰撞
    sanitize_map: dict[str, str] = {}
    seen: dict[str, str] = {}
    for r in records:
        s = sanitize_name(r.id)
        if s in seen:
            raise ValueError(f"净化后名称冲突：'{r.id}' 与 '{seen[s]}' 均映射到 '{s}'")
        seen[s] = r.id
        sanitize_map[r.id] = s

    changed = sum(1 for orig, san in sanitize_map.items() if orig != san)
    log.info("净化名称：%d 条需要修改", changed)

    # 写全部节点映射表（s10 反映射用）
    with open(out_map, "w") as f:
        f.write("original\tsanitized\n")
        for orig, san in sanitize_map.items():
            f.write(f"{orig}\t{san}\n")
    log.info("名称映射表写出：%s（%d 行）", out_map, len(sanitize_map))

    # 只输出叶节点序列（供 IQ-TREE 使用）
    leaf_records = []
    skipped = 0
    for r in records:
        if r.id in leaf_set:
            r.id = sanitize_map[r.id]
            r.description = ""
            leaf_records.append(r)
        else:
            skipped += 1

    SeqIO.write(leaf_records, out_fasta, "fasta")
    log.info(
        "叶节点 FASTA 写出：%s（%d 条，跳过内部节点 %d 条）",
        out_fasta, len(leaf_records), skipped,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta",        required=True, help="输入 FASTA（ancestral_sequences.fasta）")
    p.add_argument("--hap-graph",    required=True, help="haplogroup_graph.json")
    p.add_argument("--output-fasta", required=True, help="输出叶节点净化 FASTA")
    p.add_argument("--output-map",   required=True, help="输出全节点名称映射 TSV")
    return p


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args()
    run(
        fasta_in=args.fasta,
        hap_graph=args.hap_graph,
        out_fasta=args.output_fasta,
        out_map=args.output_map,
    )


if __name__ == "__main__":
    main()

"""
s03_build_sample_mapping.py
建立现代样本 → 单倍群 → 后代集合 的映射关系。

用法（CLI）：
    python s03_build_sample_mapping.py \
        --haplogrep META/2-单倍群信息.tsv \
        --haplogroup-graph output/stage_1/intermediate/haplogroup_graph.json \
        --vcf-sample-list output/stage_1/intermediate/all_samples.txt \
        --output-modern output/stage_1/intermediate/modern_samples.txt \
        --output-ancestor output/stage_1/intermediate/ancestor_samples.txt \
        --output-sample-hap output/stage_1/intermediate/sample_to_haplogroup.tsv \
        --output-hap-desc output/stage_1/intermediate/haplogroup_to_descendants.json

用法（import）：
    from s03_build_sample_mapping import run
    run(...)
"""

import argparse
import json
import logging
import subprocess
from collections import deque
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

# 祖先样本ID前缀
ANCESTOR_PREFIX = "Haplogroup_"


def _split_samples(all_samples: list[str]) -> tuple[list[str], list[str]]:
    """按前缀分离祖先节点样本和现代样本。"""
    ancestors = [s for s in all_samples if s.startswith(ANCESTOR_PREFIX)]
    modern    = [s for s in all_samples if not s.startswith(ANCESTOR_PREFIX)]
    log.info(f"  祖先节点样本: {len(ancestors)}")
    log.info(f"  现代样本:     {len(modern)}")
    return ancestors, modern


def _load_haplogrep(tsv: str, modern_set: set) -> dict[str, str]:
    """
    读取Haplogrep3分型结果，返回 sample_id -> haplogroup 映射。
    若同一样本有多行，取Quality最高的Rank=1行。
    只保留在现代样本集合中的样本。
    """
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    # 标准化列名（Haplogrep3输出列：SampleID, Haplogroup, Rank, Quality, ...）
    col_map = {}
    for col in df.columns:
        lower = col.lower().strip()
        if "sampleid" in lower or lower == "sampleid":
            col_map[col] = "SampleID"
        elif "haplogroup" in lower:
            col_map[col] = "Haplogroup"
        elif "quality" in lower:
            col_map[col] = "Quality"
        elif "rank" in lower:
            col_map[col] = "Rank"
    df = df.rename(columns=col_map)

    # 过滤出在VCF中存在的现代样本
    df = df[df["SampleID"].isin(modern_set)].copy()

    # 只取Rank=1的行，若有多行取Quality最高者
    if "Rank" in df.columns:
        df = df[df["Rank"].astype(str).str.strip() == "1"]
    if "Quality" in df.columns:
        df["Quality"] = pd.to_numeric(df["Quality"], errors="coerce")
        df = df.sort_values("Quality", ascending=False)
    # 去重，保留每个样本的第一行（Quality最高）
    df = df.drop_duplicates(subset=["SampleID"], keep="first")

    sample_to_hap = dict(zip(df["SampleID"], df["Haplogroup"].str.strip()))
    log.info(f"  Haplogrep3映射: {len(sample_to_hap)} 个现代样本有分型结果")
    return sample_to_hap


def _bfs_descendants(graph: dict, haplogroup: str) -> set:
    """BFS获取某单倍群的所有后代节点名称（含自身）。"""
    result = set()
    queue = deque([haplogroup])
    while queue:
        node = queue.popleft()
        if node in result:
            continue
        result.add(node)
        for child in graph.get(node, {}).get("children", []):
            queue.append(child)
    return result


def _build_hap_to_descendants(
    graph: dict,
    sample_to_hap: dict[str, str],
) -> dict[str, list[str]]:
    """
    对每个PhyloTree节点H，找出所有后代节点，
    再从sample_to_hap中筛选属于这些后代节点的现代样本。

    返回 {haplogroup: [sample_id, ...]}
    """
    # 预计算：haplogroup -> 属于该haplogroup的现代样本列表（直接分型）
    hap_to_direct_samples: dict[str, list[str]] = {}
    for sample_id, hap in sample_to_hap.items():
        hap_to_direct_samples.setdefault(hap, []).append(sample_id)

    result = {}
    n_skipped = 0

    for hap_name in graph:
        desc_nodes = _bfs_descendants(graph, hap_name)
        # 收集所有属于后代节点的现代样本
        samples = []
        for desc_node in desc_nodes:
            samples.extend(hap_to_direct_samples.get(desc_node, []))
        # 去重
        samples = list(set(samples))
        result[hap_name] = samples

    log.info(f"  共建立 {len(result)} 个单倍群的后代样本映射")
    n_with_samples = sum(1 for v in result.values() if v)
    log.info(f"  有现代样本的单倍群: {n_with_samples}")
    return result


def run(
    haplogrep: str,
    haplogroup_graph: str,
    vcf_sample_list: str,
    output_modern: str,
    output_ancestor: str,
    output_sample_hap: str,
    output_hap_desc: str,
) -> int:
    """
    建立样本-单倍群-后代映射。

    参数：
        haplogrep:         Haplogrep3 TSV路径
        haplogroup_graph:  s02输出的单倍群图JSON路径
        vcf_sample_list:   VCF全样本列表文件路径（每行一个样本ID）
        output_modern:     输出现代样本列表文件
        output_ancestor:   输出祖先节点样本列表文件
        output_sample_hap: 输出 sample_id<tab>haplogroup TSV
        output_hap_desc:   输出 haplogroup_to_descendants JSON

    返回：
        0 表示成功
    """
    for path in [output_modern, output_ancestor, output_sample_hap, output_hap_desc]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    # 读取VCF样本列表
    log.info(f"读取VCF样本列表: {vcf_sample_list}")
    with open(vcf_sample_list) as f:
        all_samples = [line.strip() for line in f if line.strip()]
    log.info(f"  VCF总样本数: {len(all_samples)}")

    # 分离祖先和现代样本
    ancestors, modern = _split_samples(all_samples)

    # 写出样本列表
    Path(output_modern).write_text("\n".join(modern) + "\n")
    Path(output_ancestor).write_text("\n".join(ancestors) + "\n")
    log.info(f"现代样本列表 -> {output_modern}")
    log.info(f"祖先样本列表 -> {output_ancestor}")

    # 读取单倍群图
    log.info(f"读取单倍群图: {haplogroup_graph}")
    with open(haplogroup_graph) as f:
        graph = json.load(f)

    # 读取Haplogrep3结果
    log.info(f"读取Haplogrep3分型结果: {haplogrep}")
    modern_set = set(modern)
    sample_to_hap = _load_haplogrep(haplogrep, modern_set)

    # 写出sample_to_haplogroup.tsv
    with open(output_sample_hap, "w") as f:
        f.write("sample_id\thaplogroup\n")
        for sid, hap in sorted(sample_to_hap.items()):
            f.write(f"{sid}\t{hap}\n")
    log.info(f"样本-单倍群映射 -> {output_sample_hap}")

    # 构建 haplogroup_to_descendants
    hap_desc = _build_hap_to_descendants(graph, sample_to_hap)

    with open(output_hap_desc, "w") as f:
        json.dump(hap_desc, f, indent=2, ensure_ascii=False)
    log.info(f"单倍群后代映射 -> {output_hap_desc}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="建立样本-单倍群-后代映射（s03）")
    p.add_argument("--haplogrep",         required=True, help="Haplogrep3 TSV路径")
    p.add_argument("--haplogroup-graph",  required=True, help="单倍群图JSON路径")
    p.add_argument("--vcf-sample-list",   required=True, help="VCF全样本列表文件")
    p.add_argument("--output-modern",     required=True, help="输出现代样本列表")
    p.add_argument("--output-ancestor",   required=True, help="输出祖先样本列表")
    p.add_argument("--output-sample-hap", required=True, help="输出sample→haplogroup TSV")
    p.add_argument("--output-hap-desc",   required=True, help="输出haplogroup→descendants JSON")
    return p


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        haplogrep=args.haplogrep,
        haplogroup_graph=args.haplogroup_graph,
        vcf_sample_list=args.vcf_sample_list,
        output_modern=args.output_modern,
        output_ancestor=args.output_ancestor,
        output_sample_hap=args.output_sample_hap,
        output_hap_desc=args.output_hap_desc,
    )


if __name__ == "__main__":
    raise SystemExit(main())

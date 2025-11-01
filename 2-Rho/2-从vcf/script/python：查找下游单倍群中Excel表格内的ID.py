#!/usr/bin/env python3
r"""
脚本：extract_haplogroup_ids.py
功能：从 phylotree 文件中提取给定单倍群及其下游（或直接子）单倍群对应的样本 ID
用法：
    python extract_haplogroup_ids.py \
        --haplo-file /path/to/phylotree.txt \
        --input-file /path/to/Chenjing_ID_Hap.tsv \
        --queries D4 B5a1 B4 \  # 多个单倍群
        [--direct-only] \
        [--output-file /path/to/output.tsv]
"""
import os
import argparse


def extract_haplogroup_ids(
    haplogroup_file_path: str,
    input_txt_path: str,
    query_haplogroups: list,
    direct_only: bool = False,
    output_file_path: str = None,
    encoding: str = 'utf-8'
) -> str:
    """
    从 phylotree 文件中找出给定单倍群的下游单倍群，
    并从输入 ID–单倍群 列表中提取对应 ID，
    最终输出 ID, QueryHaplogroup, MatchedDownstreamHaplogroup TSV。

    参数:
        haplogroup_file_path: 单倍群树文件路径（每行一条，前导\t表示层级）
        input_txt_path:      输入文件路径(ID<TAB>Haplogroup)
        query_haplogroups:   查询的单倍群列表
        direct_only:         True: 仅直接子单倍群；False: 包括所有后代及自身
        output_file_path:    输出文件路径（不指定则写到桌面）
        encoding:            文件编码

    返回:
        写入的 output_file_path
    """
    # 1. 确定输出文件路径
    if output_file_path is None:
        desktop = os.path.join(os.path.expanduser("~"), "Desktop")
        queries_str = "_".join(query_haplogroups)
        suffix = "direct" if direct_only else "all"
        fname = f"{queries_str}_{suffix}_downstream_IDs.txt"
        output_file_path = os.path.join(desktop, fname)

    # 2. 解析 phylotree，构建 (level, name) 列表
    haplo_tree = []
    with open(haplogroup_file_path, 'r', encoding=encoding) as f:
        for raw in f:
            name = raw.strip()
            if not name:
                continue
            level = raw.count('\t')
            haplo_tree.append((level, name))

    # 3. 找到每个 query 的下游单倍群（包含自身）
    query_to_downstreams = {}
    for query in query_haplogroups:
        found = False
        base_level = None
        downstream = []
        for level, name in haplo_tree:
            if not found:
                if name == query:
                    found = True
                    base_level = level
                    # 始终包含自身
                    downstream.append(name)
                continue
            # 已找到 query 之后
            if direct_only:
                if level == base_level + 1:
                    downstream.append(name)
                elif level <= base_level:
                    break
            else:
                if level > base_level:
                    downstream.append(name)
                else:
                    break
        query_to_downstreams[query] = downstream

    # 4. 读取 ID–Hap 对照表并匹配
    results = []
    with open(input_txt_path, 'r', encoding=encoding) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 2:
                continue
            sample_id, haplo = parts
            for query, downstream in query_to_downstreams.items():
                if haplo in downstream:
                    results.append((sample_id, query, haplo))

    # 5. 写入输出文件 (TSV)
    with open(output_file_path, 'w', encoding=encoding) as fout:
        fout.write("ID\tQueryHaplogroup\tMatchedDownstreamHaplogroup\n")
        for sid, q, m in results:
            fout.write(f"{sid}\t{q}\t{m}\n")

    print(f"✅ Results written to: {output_file_path}")
    return output_file_path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Extract haplogroup downstream IDs from phylotree and sample list.'
    )
    parser.add_argument(
        '--haplo-file', '-hf', required=True,
        help='Path to haplogroup tree text file'
    )
    parser.add_argument(
        '--input-file', '-i', required=True,
        help='Path to sample ID-Haplogroup TSV'
    )
    parser.add_argument(
        '--queries', '-q', required=True, nargs='+',
        help='One or more query haplogroups'
    )
    parser.add_argument(
        '--direct-only', '-d', action='store_true',
        help='If set, only extract direct child haplogroups'
    )
    parser.add_argument(
        '--output-file', '-o', default=None,
        help='Output file path (default: Desktop)'
    )
    parser.add_argument(
        '--encoding', default='utf-8',
        help='File encoding (default: utf-8)'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    extract_haplogroup_ids(
        haplogroup_file_path=args.haplo_file,
        input_txt_path=args.input_file,
        query_haplogroups=args.queries,
        direct_only=args.direct_only,
        output_file_path=args.output_file,
        encoding=args.encoding
    )

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本：自动提取单倍群列表.py
功能：对输入的 CSV（ID,Haplogroup）和 Phylotree 树文件：
    1. 去重获取所有目标单倍群；
    2. 分别提取每个目标单倍群的直接子节点和所有后代节点对应样本，
       输出到指定目录下的 direct_children.csv 与 all_descendants.csv。
用法：
    python3 自动提取单倍群列表.py \
        --haplo-file /path/to/phylotree.txt \
        --input-csv  /path/to/samples.csv \
        --out-dir     /path/to/output_directory \
        [--encoding utf-8]
"""
import os
import csv
import argparse


def parse_tree(haplo_file, encoding='utf-8'):
    """解析 Phylotree 树文件，返回 (level, name) 列表"""
    tree = []
    with open(haplo_file, 'r', encoding=encoding) as f:
        for raw in f:
            if not raw.strip():
                continue
            level = raw.count('\t')
            name = raw.strip()
            tree.append((level, name))
    return tree


def collect_downstreams(tree, query, direct_only=False):
    """
    收集 query 自身及下游节点名称列表
    direct_only=True 仅直接子节点；False 则所有后代（含自身）
    """
    found = False
    base_level = None
    downstream = []
    for level, name in tree:
        if not found:
            if name == query:
                found = True
                base_level = level
                downstream.append(name)  # 包含自身
            continue
        # 已经找到 query
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
    return downstream


def main():
    parser = argparse.ArgumentParser(
        description="自动提取直接子单倍群与所有后代对应样本"
    )
    parser.add_argument('--haplo-file', '-hf', required=True,
                        help="Phylotree 树文件路径")
    parser.add_argument('--input-csv', '-i', required=True,
                        help="样本 CSV，第一列 ID，第二列 Haplogroup")
    parser.add_argument('--out-dir', '-o', required=True,
                        help="输出目录，文件名自动生成")
    parser.add_argument('--encoding', default="utf-8",
                        help="文件编码（默认 utf-8）")
    args = parser.parse_args()

    # 1. 确保输出目录存在
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    # 2. 解析树文件
    tree = parse_tree(args.haplo_file, encoding=args.encoding)

    # 3. 读取原始 CSV 并去重获取目标单倍群列表
    samples = []
    targets = set()
    with open(args.input_csv, 'r', encoding=args.encoding, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2:
                continue
            sid, hap = row[0].strip(), row[1].strip()
            samples.append((sid, hap))
            targets.add(hap)
    targets = sorted(targets)

    # 4. 自动生成两个输出文件路径
    direct_path = os.path.join(out_dir, "直接后代.csv")
    all_path    = os.path.join(out_dir, "所有后代.csv")

    # 5. 打开输出文件并写入表头
    header = ['Target_Haplogroup', 'ID', 'Haplogroup']
    with open(direct_path, 'w', encoding=args.encoding, newline='') as fd, \
         open(all_path,    'w', encoding=args.encoding, newline='') as fa:

        writer_direct = csv.writer(fd)
        writer_all    = csv.writer(fa)
        writer_direct.writerow(header)
        writer_all.writerow(header)

        # 6. 遍历每个目标单倍群，匹配并写入
        for tgt in targets:
            downstream_direct = collect_downstreams(tree, tgt, direct_only=True)
            downstream_all    = collect_downstreams(tree, tgt, direct_only=False)

            for sid, hap in samples:
                if hap in downstream_direct:
                    writer_direct.writerow([tgt, sid, hap])
                if hap in downstream_all:
                    writer_all.writerow([tgt, sid, hap])

    print(f"✅ 直接子单倍群结果已写入：{os.path.abspath(direct_path)}")
    print(f"✅ 所有后代结果已写入：{os.path.abspath(all_path)}")


if __name__ == '__main__':
    main()

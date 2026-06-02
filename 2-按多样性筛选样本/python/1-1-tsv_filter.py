#!/usr/bin/env python3
"""
1-1-tsv_filter.py — 第一层：TSV 硬过滤（双模块）
==================================================
从 meta/总表.tsv 中筛选同时满足以下三项的样本：
  1. QC_Other  == "YES"
  2. QC_Haplogrep >= qc_min（默认 0.9）
  3. Sequence_method == "WHOLE"（排除 CHIP/CONTROL）

输出：过滤后的样本列表 TSV，保留全部原始列，并追加 major_haplogroup 列。
"""

import argparse
from pathlib import Path

import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# 核心逻辑
# ─────────────────────────────────────────────────────────────────────────────

def load_and_filter(tsv_path: str, qc_min: float) -> pd.DataFrame:
    """读取总表，应用三项硬过滤，返回过滤后的 DataFrame。"""
    df = pd.read_csv(tsv_path, sep='\t', low_memory=False)

    # QC_Haplogrep 列转数值（非数值条目置为 NaN，直接被 >= 过滤掉）
    df['QC_Haplogrep'] = pd.to_numeric(df['QC_Haplogrep'], errors='coerce')

    mask = (
        (df['QC_Other'] == 'YES') &
        (df['QC_Haplogrep'] >= qc_min) &
        (df['Sequence_method'] == 'WHOLE')
    )
    filtered = df[mask].copy()

    # 追加主单倍群标签（Haplogroup_17.2 首字母大写）
    filtered['major_haplogroup'] = (
        filtered['Haplogroup_17.2']
        .astype(str)
        .str[0]
        .str.upper()
        .fillna('Unknown')
    )

    return filtered


def run(tsv_path: str, output_path: str, qc_min: float) -> int:
    """执行过滤并写出结果，返回 0 表示成功。"""
    print(f"[1-1] 读取 TSV: {tsv_path}")
    filtered = load_and_filter(tsv_path, qc_min)

    n_total = sum(1 for _ in open(tsv_path)) - 1   # 减去表头
    n_pass  = len(filtered)
    print(f"[1-1] 原始样本: {n_total:,}  →  硬过滤后: {n_pass:,}"
          f"  (QC_Haplogrep≥{qc_min}, WHOLE, QC_Other=YES)")

    # 按主单倍群打印分布（前 15 个）
    top = filtered['major_haplogroup'].value_counts().head(15)
    print("[1-1] 主单倍群分布（前 15）:")
    for hap, cnt in top.items():
        print(f"       {hap}: {cnt:,}")

    # 只保留下游需要的列
    keep = ['ID', 'Haplogroup_17.2', 'QC_Haplogrep', 'major_haplogroup']
    out = filtered[[c for c in keep if c in filtered.columns]]

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep='\t', index=False)
    print(f"[1-1] 结果写入: {output_path}")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# CLI 入口
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description='Step 1-1: TSV 硬过滤',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--tsv',     required=True, help='总表.tsv 路径')
    p.add_argument('--output',  required=True, help='输出 TSV 路径')
    p.add_argument('--qc-min',  type=float, default=0.9,
                   help='QC_Haplogrep 最低阈值（含）')
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return run(
        tsv_path=args.tsv,
        output_path=args.output,
        qc_min=args.qc_min,
    )


if __name__ == '__main__':
    raise SystemExit(main())

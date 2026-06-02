#!/usr/bin/env python3
"""
1-2-vcf_fps_select.py — 第二层：VCF 成对差异 + FPS 降采样（双模块）
====================================================================
输入：第一层硬过滤后的样本列表（含 major_haplogroup 列）。
流程（按主单倍群分块）：
  1. 用 bcftools 从 VCF 提取块内样本基因型
  2. 计算块内成对差异矩阵（忽略缺失位点，WHOLE 样本缺失率 ≈ 0%）
  3. 并查集折叠差异=0 的相同单倍型，每组保留一个代表
  4. 自适应 FPS：直到所有未选样本与已选集合最近距离 <= min_dist 时停止

输出 selected_samples.tsv 列说明：
  ID                  ← 样本 ID
  major_haplogroup    ← 主单倍群（首字母）
  haplogroup          ← 原始 Haplogroup_17.2
  QC_Haplogrep        ← 质量分
  cluster_size        ← 该样本所代表的相同单倍型数量（1 = 唯一型）
  fps_rank            ← FPS 选取顺序（1 = 起始种子点）
  fps_dist_to_nearest ← 被选中时与最近已选样本的成对差异（种子点为 -1）
  reason              ← fps_seed | fps_diverse
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# VCF 基因型提取
# ─────────────────────────────────────────────────────────────────────────────

def extract_gt_matrix(vcf_path: str, sample_ids: list,
                      bcftools: str) -> tuple:
    """
    用 bcftools 提取指定样本的基因型矩阵。

    返回：
      G           : ndarray, shape (n_sites, n_samples), dtype int8
                    0=参考态, k=alt等位序号, -1=缺失
      vcf_samples : 与 G 列对应的样本 ID 列表（只含 VCF 中实际存在的样本）
    """
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as f:
        f.write('\n'.join(sample_ids) + '\n')
        tmp = f.name

    try:
        # 获取 VCF 中实际存在的样本顺序
        vcf_samples = subprocess.run(
            [bcftools, 'query', '-S', tmp, '-l', vcf_path],
            capture_output=True, text=True, check=True
        ).stdout.strip().split('\n')
        vcf_samples = [s for s in vcf_samples if s]

        if not vcf_samples:
            return np.empty((0, 0), dtype=np.int8), []

        # 提取基因型（每行一个位点，tab 分隔 GT 字段）
        gt_text = subprocess.run(
            [bcftools, 'query', '-S', tmp, '-f', '[%GT\t]\n', vcf_path],
            capture_output=True, text=True, check=True
        ).stdout

        n = len(vcf_samples)
        rows = []
        for line in gt_text.splitlines():
            parts = line.rstrip('\t').split('\t')
            row = np.empty(n, dtype=np.int8)
            for i, gt in enumerate(parts[:n]):
                a = gt.split('/')[0]
                row[i] = -1 if a == '.' else int(a)
            rows.append(row)

        if not rows:
            return np.empty((0, n), dtype=np.int8), vcf_samples

        return np.array(rows, dtype=np.int8), vcf_samples

    finally:
        os.unlink(tmp)


def get_segregating_mask(G: np.ndarray) -> np.ndarray:
    """返回在当前块内有变异的位点布尔掩码（行索引）。"""
    return np.array(
        [len(np.unique(G[i][G[i] >= 0])) > 1 for i in range(G.shape[0])]
    )


# ─────────────────────────────────────────────────────────────────────────────
# 成对差异矩阵
# ─────────────────────────────────────────────────────────────────────────────

def pairwise_diff(Gt: np.ndarray, chunk: int = 2000) -> np.ndarray:
    """
    计算成对差异矩阵。Gt shape: (n_samples, n_var_sites), dtype int8, -1=缺失。
    返回 D: shape (n, n), dtype int16，对称。
    分块计算以控制内存峰值。
    """
    n = Gt.shape[0]
    D = np.zeros((n, n), dtype=np.int16)

    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        Gi = Gt[start:end]          # (chunk, S)
        valid_i = Gi >= 0           # (chunk, S)
        for j in range(n):
            gj = Gt[j]
            both = valid_i & (gj >= 0)
            D[start:end, j] = ((Gi != gj) & both).sum(axis=1)

    return np.minimum(D, D.T)      # 对称化（消除缺失引起的微小不对称）


# ─────────────────────────────────────────────────────────────────────────────
# 相同单倍型折叠（并查集）
# ─────────────────────────────────────────────────────────────────────────────

def collapse_identical(D: np.ndarray) -> tuple:
    """
    将成对差异=0 的样本归为同一单倍型组，返回：
      reps     : 各组代表索引 ndarray
      clusters : {代表索引: [所有成员索引]}
    """
    n = D.shape[0]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    zi, zj = np.where(np.triu(D == 0, k=1))
    for a, b in zip(zi, zj):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[ra] = rb

    clusters: dict = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(i)

    reps = np.array(list(clusters.keys()))
    return reps, clusters


# ─────────────────────────────────────────────────────────────────────────────
# 自适应 FPS
# ─────────────────────────────────────────────────────────────────────────────

def fps_adaptive(D_reps: np.ndarray, min_dist: int,
                 max_tips: int) -> list:
    """
    自适应最远点采样，返回 (sel_indices, dist_when_selected) 两个列表。

    sel_indices         : 选中点在 D_reps 中的行索引（按选入顺序）
    dist_when_selected  : 各点被选中时与已选集合的最近距离
                          （种子点记为 -1）

    停止条件（满足任一）：
      A. 所有未选样本与已选集合最近距离 <= min_dist
      B. 已选数量达到 max_tips（0 = 不限）
    """
    m = D_reps.shape[0]
    if m == 0:
        return [], []
    if m == 1:
        return [0], [-1]

    sel_indices = [0]
    dist_record = [-1]                       # 种子点距离记为 -1

    mind = D_reps[0].astype(np.int32).copy()
    mind[0] = 0

    while True:
        best_dist = int(mind.max())
        if best_dist <= min_dist:
            break
        if max_tips > 0 and len(sel_indices) >= max_tips:
            break

        nxt = int(np.argmax(mind))
        dist_record.append(best_dist)
        sel_indices.append(nxt)
        mind = np.minimum(mind, D_reps[nxt].astype(np.int32))
        mind[nxt] = 0

    return sel_indices, dist_record


# ─────────────────────────────────────────────────────────────────────────────
# 单块处理
# ─────────────────────────────────────────────────────────────────────────────

def process_block(block_name: str, block_df: pd.DataFrame,
                  vcf_path: str, bcftools: str,
                  min_dist: int, max_tips: int,
                  temp_dir: str, overwrite: bool) -> pd.DataFrame:
    """
    处理单个主单倍群块，返回筛选结果 DataFrame（含 reason 等列）。
    若 overwrite=False 且已有缓存，直接读取缓存文件。
    """
    cache = Path(temp_dir) / f"block_{block_name}_selected.tsv"
    if (not overwrite) and cache.exists():
        print(f"  [块 {block_name}] 使用缓存: {cache}")
        return pd.read_csv(cache, sep='\t')

    sample_ids = block_df['ID'].tolist()
    n_input = len(sample_ids)
    print(f"\n  [块 {block_name}] 输入样本: {n_input:,}")

    # ── VCF 基因型提取 ──────────────────────────────────────────────────────
    G, vcf_samples = extract_gt_matrix(vcf_path, sample_ids, bcftools)
    if not vcf_samples:
        print(f"  [块 {block_name}] 警告：VCF 中无匹配样本，跳过。")
        return pd.DataFrame()

    n_vcf = len(vcf_samples)
    if n_vcf < n_input:
        print(f"  [块 {block_name}] {n_input - n_vcf} 个样本不在 VCF 中，已排除。")

    # 将 DataFrame 对齐到 VCF 实际样本顺序
    vcf_set = set(vcf_samples)
    block_df = (block_df[block_df['ID'].isin(vcf_set)]
                .set_index('ID')
                .loc[vcf_samples]
                .reset_index())

    # ── 块内分离位点 ────────────────────────────────────────────────────────
    seg_mask = get_segregating_mask(G)
    n_seg = int(seg_mask.sum())
    print(f"  [块 {block_name}] 块内分离位点: {n_seg:,}")

    if n_seg == 0:
        # 所有样本完全相同，取第 1 个作为唯一代表
        result = block_df.iloc[[0]].copy()
        result['cluster_size']        = n_vcf
        result['fps_rank']            = 1
        result['fps_dist_to_nearest'] = -1
        result['reason']              = 'fps_seed'
        result.to_csv(cache, sep='\t', index=False)
        return result

    Gt = G[seg_mask].T.copy()   # (n_samples, n_var_sites)

    # ── 成对差异矩阵 ────────────────────────────────────────────────────────
    print(f"  [块 {block_name}] 计算成对差异矩阵 ({n_vcf} × {n_vcf})...")
    D = pairwise_diff(Gt)

    # ── 折叠相同单倍型 ──────────────────────────────────────────────────────
    reps, clusters = collapse_identical(D)
    n_unique = len(reps)
    # 每个代表对应的簇大小（相同单倍型数）
    cluster_size_map = {rep: len(members) for rep, members in clusters.items()}
    print(f"  [块 {block_name}] 唯一单倍型: {n_unique:,}  "
          f"（折叠 {n_vcf - n_unique:,} 个重复）")

    Dr = D[np.ix_(reps, reps)]

    # ── 自适应 FPS ──────────────────────────────────────────────────────────
    sel_local, dist_record = fps_adaptive(Dr, min_dist=min_dist, max_tips=max_tips)
    sel_global = reps[sel_local]   # 在 vcf_samples 数组中的原始索引
    n_sel = len(sel_global)
    print(f"  [块 {block_name}] FPS 选中: {n_sel:,}  (min_dist={min_dist})")

    # ── 组装结果 ────────────────────────────────────────────────────────────
    rows = block_df.iloc[sel_global].copy()

    reasons = []
    for rank, (loc_i, dist) in enumerate(zip(sel_local, dist_record), start=1):
        csize = cluster_size_map[reps[loc_i]]
        if rank == 1:
            reasons.append(f'fps_seed;block={block_name};cluster={csize}')
        else:
            reasons.append(
                f'fps_diverse;block={block_name};rank={rank};dist={dist};cluster={csize}'
            )
    rows['Reason'] = reasons

    rows.to_csv(cache, sep='\t', index=False)
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────────────────────

def run(filtered_list: str, vcf_path: str, output_path: str,
        summary_path: str, temp_dir: str, bcftools: str,
        min_dist: int, max_tips: int, overwrite: bool) -> int:
    """完整第二层流程，返回 0 表示成功。"""
    print(f"[1-2] 读取硬过滤结果: {filtered_list}")
    df = pd.read_csv(filtered_list, sep='\t', low_memory=False)

    if 'major_haplogroup' not in df.columns:
        df['major_haplogroup'] = (
            df['Haplogroup_17.2'].astype(str).str[0].str.upper().fillna('Unknown')
        )

    Path(temp_dir).mkdir(parents=True, exist_ok=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(summary_path).parent.mkdir(parents=True, exist_ok=True)

    blocks = sorted(df['major_haplogroup'].unique())
    print(f"[1-2] 主单倍群块: {blocks}")

    all_results = []
    summary_rows = []

    for block_name in blocks:
        block_df = df[df['major_haplogroup'] == block_name].copy()
        n_input = len(block_df)

        result = process_block(
            block_name=block_name,
            block_df=block_df,
            vcf_path=vcf_path,
            bcftools=bcftools,
            min_dist=min_dist,
            max_tips=max_tips,
            temp_dir=temp_dir,
            overwrite=overwrite,
        )
        if len(result) == 0:
            continue

        all_results.append(result)

        # 汇总行（取第一行的统计列，它们对整块相同）
        first = result.iloc[0]
        summary_rows.append({
            'major_haplogroup': block_name,
            'n_layer1_input':   n_input,
            'n_fps_selected':   len(result),
        })

    if not all_results:
        print("[ERROR] 所有块均为空，无输出。")
        return 1

    # ── 输出合并结果 ────────────────────────────────────────────────────────
    combined = pd.concat(all_results, ignore_index=True)

    # 最终输出仅保留 4 列
    out_cols = ['ID', 'Haplogroup_17.2', 'QC_Haplogrep', 'Reason']
    combined[out_cols].to_csv(output_path, sep='\t', index=False)
    print(f"[1-2] 最终选中样本: {len(combined):,}  →  {output_path}")

    # ── 汇总报告 ─────────────────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"[1-2] 汇总报告: {summary_path}")

    # 终端打印摘要
    print("\n─── 各块筛选摘要 ────────────────────────────────────────────────")
    print(f"  {'块':>4}  {'Layer1输入':>10}  {'FPS选中':>8}")
    for _, row in summary_df.iterrows():
        print(f"  {row['major_haplogroup']:>4}  "
              f"{int(row['n_layer1_input']):>10,}  "
              f"{int(row['n_fps_selected']):>8,}")
    print(f"\n  合计: {summary_df['n_fps_selected'].sum():,} 个样本")
    return 0


# ─────────────────────────────────────────────────────────────────────────────
# CLI 入口
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description='Step 1-2: VCF 成对差异 + FPS 降采样',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--filtered-list', required=True,
                   help='第一层过滤结果 TSV（含 major_haplogroup 列）')
    p.add_argument('--vcf',           required=True,
                   help='merged_clean.vcf.gz 路径')
    p.add_argument('--output',        required=True,
                   help='最终输出 TSV 路径（selected_samples.tsv）')
    p.add_argument('--summary',       required=True,
                   help='汇总报告 TSV 路径（selection_summary.tsv）')
    p.add_argument('--temp',          required=True,
                   help='临时目录（各块缓存存放于此）')
    p.add_argument('--bcftools',      required=True,
                   help='bcftools 可执行路径')
    p.add_argument('--min-dist',      type=int, default=1,
                   help='FPS 停止阈值（未选样本最近距离 <= 此值时停止）')
    p.add_argument('--max-tips',      type=int, default=0,
                   help='每块最大 tip 数（0 = 不限）')
    p.add_argument('--overwrite',     default='false',
                   help='是否覆盖块缓存（true/false）')
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    overwrite = args.overwrite.lower() in ('true', '1', 'yes')
    return run(
        filtered_list=args.filtered_list,
        vcf_path=args.vcf,
        output_path=args.output,
        summary_path=args.summary,
        temp_dir=args.temp,
        bcftools=args.bcftools,
        min_dist=args.min_dist,
        max_tips=args.max_tips,
        overwrite=overwrite,
    )


if __name__ == '__main__':
    raise SystemExit(main())

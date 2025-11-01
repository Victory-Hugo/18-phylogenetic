
# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# 按“单倍群 × 基因区段”并行计算 mtDNA ρ (rho) 与 σ。
# 输入：对齐的多序列 FASTA
# 输出：与旧 VCF 版完全一致，支持断点续跑 (*.ckpt)
# """

# import os, sys, csv, math, argparse, logging
# import numpy as np
# from multiprocessing import Pool

# # =========================================
# # ❶ 五个区段及 FASTA 文件名
# # =========================================
# regions = [
#     ('CS',          'RHO_Complete_sequence.fasta'),
#     ('HVSI',        'RHO_HVS-I.fasta'),
#     ('HVSI_TRANSI', 'RHO_HVS-I_transitions.fasta'),
#     ('HVSII',       'RHO_HVS-II.fasta'),
#     ('CR',          'RHO_Control_Region.fasta'),
# ]

# # =========================================
# # ❷ FASTA → 整型矩阵
# # =========================================
# def _read_fasta(path):
#     """读取多序列 FASTA，返回 {seq_id: sequence}；要求所有序列等长"""
#     seqs, cur_id, buf = {}, None, []
#     with open(path) as f:
#         for ln in f:
#             if ln.startswith('>'):
#                 if cur_id:
#                     seqs[cur_id] = ''.join(buf).upper()
#                 cur_id = ln[1:].split()[0]
#                 buf = []
#             else:
#                 buf.append(ln.strip())
#         if cur_id:
#             seqs[cur_id] = ''.join(buf).upper()

#     lens = {len(s) for s in seqs.values()}
#     if len(lens) != 1:
#         raise ValueError(f'FASTA {path} 含不同长度序列，无法对齐')
#     return seqs

# def build_int_matrix_fasta(fa_file, samples):
#     """返回 shape=(多态位点数, 样本数) 的 np.int8 矩阵"""
#     seqs = _read_fasta(fa_file)
#     present = [s for s in samples if s in seqs]
#     if len(present) < 2:
#         logging.warning(f'{fa_file}: 可用样本 {len(present)} < 2，跳过')
#         return np.zeros((0, len(samples)), dtype=np.int8)

#     base2int = {b:i for i, b in enumerate('ACGT')}
#     def enc(ch): return base2int.get(ch, -1)

#     L = len(next(iter(seqs.values())))
#     rows = []
#     for pos in range(L):
#         col = np.fromiter((enc(seqs[s][pos]) for s in present),
#                           dtype=np.int8, count=len(present))
#         if len({x for x in col if x >= 0}) >= 2:
#             rows.append(col)
#     if not rows:
#         logging.warning(f'{fa_file}: 无多态位点')
#         return np.zeros((0, len(present)), dtype=np.int8)

#     logging.debug(f'{fa_file}: 构建矩阵 {len(rows)}×{len(present)} 完成')
#     return np.vstack(rows)

# # =========================================
# # ❸ rho / sigma 计算
# # =========================================
# def rho_sigma_fast(A: np.ndarray, block_size: int = 20):
#     if A.shape[1] < 2 or A.shape[0] == 0:
#         return float('nan'), float('nan')
#     M = np.where(A < 0, 127, A).astype(np.int8)
#     n_sites, n_samples = M.shape
#     D = np.zeros((n_samples, n_samples), dtype=np.int32)
#     for st in range(0, n_sites, block_size):
#         Mb = M[st:st+block_size]
#         diff = Mb[:, :, None] != Mb[:, None, :]
#         D += diff.sum(axis=0).astype(np.int32)
#     tri = D[np.triu_indices(n_samples, 1)]
#     rho   = tri.mean()
#     sigma = math.sqrt((tri.astype(np.int64)**2).sum() / (tri.size**2))
#     return rho, sigma

# def compute_one(args):
#     hap, samples, fa_dir, label, fname = args
#     path = os.path.join(fa_dir, fname)
#     if not os.path.isfile(path):
#         logging.error(f'缺失文件 {path}')
#         return hap, label, float('nan'), float('nan')
#     logging.info(f'[{hap}] {label}: 处理 {fname}')
#     G = build_int_matrix_fasta(path, samples)
#     rho, sigma = rho_sigma_fast(G)
#     logging.info(f'[{hap}] {label}: rho={rho:.4f}, sigma={sigma:.4f}')
#     return hap, label, rho, sigma

# # =========================================
# # ❹ 主程序
# # =========================================
# def main():
#     ap = argparse.ArgumentParser(
#         description='基于 FASTA 的 mtDNA ρ/σ 计算（支持断点续跑）')
#     ap.add_argument('hap_csv',   help='两列 CSV：ID,Haplogroup')
#     ap.add_argument('fasta_dir', help='存放 RHO_*.fasta 的目录')
#     ap.add_argument('output_tsv',help='最终 TSV 输出路径')
#     ap.add_argument('-t','--threads', type=int, default=1, help='并行进程数')
#     ap.add_argument('-v','--verbose', action='store_true', help='调试模式')
#     args = ap.parse_args()

#     logging.basicConfig(
#         format='[%(asctime)s] %(levelname)s: %(message)s',
#         level=logging.DEBUG if args.verbose else logging.INFO,
#         datefmt='%Y-%m-%d %H:%M:%S')

#     # ---------- 读取 CSV ----------
#     hap2samples = {}
#     with open(args.hap_csv) as f:
#         rdr = csv.DictReader(f)
#         if {'ID','Haplogroup'} - set(rdr.fieldnames):
#             logging.error('CSV 必须含列 ID 和 Haplogroup')
#             sys.exit(1)
#         for r in rdr:
#             hap2samples.setdefault(r['Haplogroup'], []).append(r['ID'])

#     # ---------- 过滤样本不足 2 ----------
#     skipped = []
#     for h, s in list(hap2samples.items()):
#         if len(s) < 2:
#             skipped.append(h)
#             logging.warning(f'跳过 {h}（样本数 {len(s)} < 2）')
#             hap2samples.pop(h)
#     if not hap2samples:
#         logging.info('所有单倍群样本数不足 2 或已完成，将直接进入汇总')

#     # ---------- 断点文件 ----------
#     ckpt = args.output_tsv + '.ckpt'
#     done = set()
#     if os.path.isfile(ckpt):
#         with open(ckpt) as f:
#             for r in csv.DictReader(f, delimiter='\t'):
#                 done.add((r['Haplogroup'], r['Region']))
#     else:
#         with open(ckpt,'w') as f:
#             csv.writer(f, delimiter='\t').writerow(
#                 ['Haplogroup','Region','RHO','RHO_SE'])

#     # ---------- 待执行任务 ----------
#     tasks = [(hap, samples, args.fasta_dir, lbl, fname)
#              for hap, samples in hap2samples.items()
#              for lbl, fname in regions
#              if (hap, lbl) not in done]
#     logging.info(f'待执行任务 {len(tasks)}')

#     if tasks:
#         with Pool(args.threads) as pool, open(ckpt,'a') as f:
#             wr = csv.writer(f, delimiter='\t')
#             for hap, lbl, rho, sig in pool.imap_unordered(compute_one, tasks):
#                 wr.writerow([hap, lbl, f'{rho:.6g}', f'{sig:.6g}'])
#                 f.flush()

#     # ---------- 汇总 ----------
#     summary = {}
#     with open(ckpt) as f:
#         for r in csv.DictReader(f, delimiter='\t'):
#             summary.setdefault(r['Haplogroup'], {})[r['Region']] = (
#                 float(r['RHO']), float(r['RHO_SE']))

#     with open(args.output_tsv,'w') as out:
#         header = ['Haplogroup'] + [
#             f'RHO_{l}\tRHO_{l}_SE' for l, _ in regions]
#         out.write('\t'.join(header) + '\n')
#         for hap in sorted(summary):
#             row = [hap]
#             for lbl, _ in regions:
#                 rho, sig = summary[hap].get(lbl, (float('nan'), float('nan')))
#                 row += [str(rho), str(sig)]
#             out.write('\t'.join(row) + '\n')
#     logging.info('全部完成 ✔️')

# if __name__ == '__main__':
#     main()

#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
基于方案 A：全局预加载 FASTA + 复用子进程
按“单倍群 × 基因区段”并行计算 mtDNA ρ (rho) 与 σ。
输入：对齐的多序列 FASTA（只读一次）
输出：与旧 VCF 版完全一致，支持断点续跑 (*.ckpt)
"""

import os
import sys
import csv
import math
import argparse
import logging
import numpy as np
import multiprocessing as mp
from functools import lru_cache
from multiprocessing import Pool

# =========================================
# ❶ 区段列表 & FASTA 文件名
# =========================================
regions = [
    ('CS',          'RHO_Complete_sequence.fasta'),
    ('HVSI',        'RHO_HVS-I.fasta'),
    ('HVSI_TRANSI', 'RHO_HVS-I_transitions.fasta'),
    ('HVSII',       'RHO_HVS-II.fasta'),
    ('CR',          'RHO_Control_Region.fasta'),
]

# =========================================
# ❷ 全局加载 FASTA（只读一次）
# =========================================
@lru_cache(maxsize=None)
def load_region(fasta_path):
    """
    读取多序列 FASTA，返回：
      ids: [seq_id, ...]
      arr: np.ndarray(uint8) shape=(n_ids, seq_len)
    其中 A/C/G/T 分别编码为 0/1/2/3，其他（缺失）为 255
    """
    seqs = {}
    cur_id = None
    buf = []
    with open(fasta_path) as fh:
        for ln in fh:
            if ln.startswith('>'):
                if cur_id:
                    seqs[cur_id] = ''.join(buf).upper()
                cur_id = ln[1:].split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if cur_id:
            seqs[cur_id] = ''.join(buf).upper()

    # 检查对齐
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(f'FASTA {fasta_path} 含不同长度序列，无法对齐')
    seq_len = lengths.pop()
    ids = list(seqs.keys())

    # 构建查表
    look = np.full(256, 255, dtype=np.uint8)
    look[ord('A')] = 0
    look[ord('C')] = 1
    look[ord('G')] = 2
    look[ord('T')] = 3

    # 按 ids 顺序编码成 ndarray
    arr = np.zeros((len(ids), seq_len), dtype=np.uint8)
    for i, sid in enumerate(ids):
        data = seqs[sid].encode('ascii')
        arr[i, :] = look[np.frombuffer(data, dtype=np.uint8)]

    return ids, arr

# =========================================
# ❸ 构建每个区段的子矩阵
# =========================================
def build_int_matrix_cached(fasta_dir, fname, samples):
    """
    基于 load_region 缓存，针对所给 samples：
    返回 shape=(n_polymorphic_sites, n_present_samples) 的 np.int8 矩阵
    """
    path = os.path.join(fasta_dir, fname)
    ids, full = load_region(path)  # 共享内存
    # 找到需要的样本行索引
    idx = [ids.index(s) for s in samples if s in ids]
    n = len(idx)
    if n < 2:
        logging.warning(f'{fname}：可用样本 {n} < 2，跳过')
        return np.zeros((0, n), dtype=np.int8)

    sub = full[idx, :]  # shape = (n, seq_len)
    seq_len = sub.shape[1]
    rows = []
    for pos in range(seq_len):
        col = sub[:, pos]
        # 255 为缺失
        valid = col != 255
        alleles = col[valid]
        if alleles.size < 2:
            continue
        # 判断多态
        if np.unique(alleles).size >= 2:
            rows.append(alleles.astype(np.int8))
    if not rows:
        logging.warning(f'{fname}：无多态位点')
        return np.zeros((0, n), dtype=np.int8)

    mat = np.vstack(rows)  # shape = (#sites, n)
    logging.debug(f'{fname}：构建矩阵 {mat.shape[0]}×{mat.shape[1]} 完成')
    return mat

# =========================================
# ❹ rho / sigma 计算
# =========================================
def rho_sigma_fast(A: np.ndarray, block_size: int = 20):
    if A.shape[1] < 2 or A.shape[0] == 0:
        return float('nan'), float('nan')
    M = np.where(A < 0, 127, A).astype(np.int8)
    n_sites, n_samples = M.shape
    D = np.zeros((n_samples, n_samples), dtype=np.int32)
    for st in range(0, n_sites, block_size):
        Mb = M[st:st+block_size]
        diff = Mb[:, :, None] != Mb[:, None, :]
        D += diff.sum(axis=0).astype(np.int32)
    iu = np.triu_indices(n_samples, 1)
    tri = D[iu]
    rho   = tri.mean()
    sigma = math.sqrt((tri.astype(np.int64)**2).sum() / (tri.size**2))
    return rho, sigma

def compute_one(args):
    hap, samples, fasta_dir, label, fname = args
    logging.info(f'[{hap}] 区域 {label} 开始：{fname}')
    G = build_int_matrix_cached(fasta_dir, fname, samples)
    rho, sigma = rho_sigma_fast(G)
    logging.info(f'[{hap}] 区域 {label} 完成：rho={rho:.4f}, sigma={sigma:.4f}')
    return hap, label, rho, sigma

# =========================================
# ❺ 主程序
# =========================================
def main():
    # 切换为 forkserver，以便预加载数组在子进程共享
    mp.set_start_method('forkserver', force=True)

    p = argparse.ArgumentParser(
        description='基于预加载 FASTA 的 mtDNA ρ/σ 计算（方案 A，支持断点续跑）')
    p.add_argument('hap_csv',    help='两列 CSV：ID,Haplogroup')
    p.add_argument('fasta_dir',  help='存放 RHO_*.fasta 的目录')
    p.add_argument('output_tsv', help='最终 TSV 输出路径')
    p.add_argument('-t','--threads', type=int, default=1,
                   help='并行进程数')
    p.add_argument('-v','--verbose', action='store_true',
                   help='调试日志')
    args = p.parse_args()

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)s: %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    # --- 读取单倍群与样本列表 ---
    hap2samples = {}
    with open(args.hap_csv) as f:
        rdr = csv.DictReader(f)
        if {'ID','Haplogroup'} - set(rdr.fieldnames):
            logging.error('CSV 必须包含列 ID 和 Haplogroup')
            sys.exit(1)
        for r in rdr:
            hap2samples.setdefault(r['Haplogroup'], []).append(r['ID'])

    # --- 过滤样本不足 2 的单倍群 ---
    for h in list(hap2samples):
        if len(hap2samples[h]) < 2:
            logging.warning(f'跳过 {h}（样本数 {len(hap2samples[h])} < 2）')
            hap2samples.pop(h)

    # --- 准备断点文件 ---
    ckpt = args.output_tsv + '.ckpt'
    done = set()
    if os.path.isfile(ckpt):
        with open(ckpt) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                done.add((r['Haplogroup'], r['Region']))
    else:
        with open(ckpt,'w') as f:
            csv.writer(f, delimiter='\t').writerow(
                ['Haplogroup','Region','RHO','RHO_SE'])

    # --- 构建任务列表 ---
    tasks = []
    for hap, samples in hap2samples.items():
        for label, fname in regions:
            if (hap, label) in done:
                logging.debug(f'跳过 {(hap,label)}（已完成）')
                continue
            tasks.append((hap, samples, args.fasta_dir, label, fname))

    logging.info(f'待执行任务数：{len(tasks)}')
    if tasks:
        with Pool(args.threads) as pool, open(ckpt,'a') as f:
            writer = csv.writer(f, delimiter='\t')
            for hap, label, rho, sigma in pool.imap_unordered(compute_one, tasks):
                writer.writerow([hap, label, f'{rho:.6g}', f'{sigma:.6g}'])
                f.flush()

    # --- 汇总到最终 TSV ---
    summary = {}
    with open(ckpt) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            summary.setdefault(r['Haplogroup'], {})[r['Region']] = (
                float(r['RHO']), float(r['RHO_SE'])
            )

    with open(args.output_tsv, 'w') as out:
        header = ['Haplogroup']
        for label, _ in regions:
            header += [f'RHO_{label}', f'RHO_{label}_SE']
        out.write('\t'.join(header) + '\n')
        for hap in sorted(summary):
            row = [hap]
            for label, _ in regions:
                rho, sigma = summary[hap].get(label, (float('nan'), float('nan')))
                row += [str(rho), str(sigma)]
            out.write('\t'.join(row) + '\n')

    logging.info('全部完成 ✔️')

if __name__ == '__main__':
    main()

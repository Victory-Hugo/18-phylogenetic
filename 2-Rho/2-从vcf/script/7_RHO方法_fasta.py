#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按"单倍群 × 基因区段"并行计算 mtDNA ρ (rho) 与 σ。
输入：对齐的多序列 FASTA
输出：支持断点续跑 (*.ckpt)
"""

import os
import sys
import csv
import math
import argparse
import logging
import numpy as np
from functools import lru_cache
from multiprocessing import Pool

regions = [
    ('CS',          'RHO_Complete_sequence.fasta'),
    ('HVSI',        'RHO_HVS-I.fasta'),
    ('HVSI_TRANSI', 'RHO_HVS-I_transitions.fasta'),
    ('HVSII',       'RHO_HVS-II.fasta'),
    ('CR',          'RHO_Control_Region.fasta'),
]


def validate_csv_file(csv_path):
    """验证CSV文件格式"""
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"CSV文件不存在: {csv_path}")
    
    logging.info(f"验证CSV文件: {csv_path}")
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} 是空文件")
        
        if {'ID', 'Haplogroup'} - set(reader.fieldnames):
            raise ValueError(f"{csv_path} 必须包含列 'ID' 和 'Haplogroup'")
        
        row_num = 1
        for row in reader:
            row_num += 1
            sample_id = row.get('ID', '').strip()
            haplogroup = row.get('Haplogroup', '').strip()
            
            if not sample_id:
                raise ValueError(f"{csv_path} 第{row_num}行: ID为空")
            if not haplogroup:
                raise ValueError(f"{csv_path} 第{row_num}行: Haplogroup为空")
            
            for char in ['-', '.', '\n', '\r', '\t']:
                if char in sample_id:
                    raise ValueError(
                        f"{csv_path} 第{row_num}行: ID '{sample_id}' 包含禁止字符 '{repr(char)}'\n"
                        f"请清理CSV文件并重新运行"
                    )
                if char in haplogroup:
                    raise ValueError(
                        f"{csv_path} 第{row_num}行: Haplogroup '{haplogroup}' 包含禁止字符 '{repr(char)}'\n"
                        f"请清理CSV文件并重新运行"
                    )
    
    logging.info(f"✓ CSV文件格式正确（{row_num-1}条数据）")


def validate_fasta_files(fasta_dir):
    """验证FASTA文件格式"""
    logging.info(f"验证FASTA文件目录: {fasta_dir}")
    
    if not os.path.isdir(fasta_dir):
        raise NotADirectoryError(f"FASTA目录不存在: {fasta_dir}")
    
    for label, fname in regions:
        fasta_path = os.path.join(fasta_dir, fname)
        
        if not os.path.isfile(fasta_path):
            raise FileNotFoundError(f"FASTA文件不存在: {fasta_path}")
        
        logging.info(f"  验证 {fname} ...")
        
        seqs = {}
        cur_id = None
        buf = []
        line_num = 0
        
        with open(fasta_path) as f:
            for line in f:
                line_num += 1
                if line.startswith('>'):
                    if cur_id:
                        seq_str = ''.join(buf).upper()
                        for char in seq_str:
                            if char not in 'ACGT':
                                raise ValueError(
                                    f"{fasta_path} 序列 '{cur_id}' 中包含无效字符 '{char}'（ASCII: {ord(char)}）\n"
                                    f"序列只允许包含 A/C/G/T\n"
                                    f"请检查是否包含缺失标记 -, ., N 或其他字符\n"
                                    f"请格式化FASTA文件后重新运行"
                                )
                        seqs[cur_id] = seq_str
                    
                    cur_id = line[1:].split()[0]
                    if not cur_id:
                        raise ValueError(f"{fasta_path} 第{line_num}行: 序列名为空")
                    buf = []
                else:
                    buf.append(line.rstrip('\n\r'))
            
            if cur_id:
                seq_str = ''.join(buf).upper()
                for char in seq_str:
                    if char not in 'ACGT':
                        raise ValueError(
                            f"{fasta_path} 序列 '{cur_id}' 中包含无效字符 '{char}'（ASCII: {ord(char)}）\n"
                            f"序列只允许包含 A/C/G/T\n"
                            f"请检查是否包含缺失标记 -, ., N 或其他字符\n"
                            f"请格式化FASTA文件后重新运行"
                        )
                seqs[cur_id] = seq_str
        
        if not seqs:
            raise ValueError(f"{fasta_path} 文件为空或不包含有效序列")
        
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            length_dist = {}
            for sid, seq in seqs.items():
                length_dist.setdefault(len(seq), []).append(sid)
            
            detail = "\n".join(
                f"  {length} bp: {len(samples)} sequences"
                for length, samples in sorted(length_dist.items())
            )
            raise ValueError(
                f"{fasta_path} 包含不同长度的序列，无法对齐：\n{detail}\n"
                f"请对齐所有序列后重新运行"
            )
        
        logging.info(f"    ✓ {len(seqs)} sequences, 均为 {list(lengths)[0]} bp")


@lru_cache(maxsize=None)
def load_region(fasta_path):
    """读取FASTA并编码为整型矩阵"""
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

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(f'FASTA {fasta_path} 含不同长度序列，无法对齐')
    seq_len = lengths.pop()
    ids = list(seqs.keys())

    look = np.full(256, 255, dtype=np.uint8)
    look[ord('A')] = 0
    look[ord('C')] = 1
    look[ord('G')] = 2
    look[ord('T')] = 3

    arr = np.zeros((len(ids), seq_len), dtype=np.uint8)
    for i, sid in enumerate(ids):
        data = seqs[sid].encode('ascii')
        arr[i, :] = look[np.frombuffer(data, dtype=np.uint8)]

    return ids, arr


def build_int_matrix_cached(fasta_dir, fname, samples):
    """构建多态位点矩阵"""
    path = os.path.join(fasta_dir, fname)
    ids, full = load_region(path)
    idx = [ids.index(s) for s in samples if s in ids]
    n = len(idx)
    if n < 2:
        logging.warning(f'{fname}：可用样本 {n} < 2，跳过')
        return np.zeros((0, n), dtype=np.int8)

    sub = full[idx, :]
    seq_len = sub.shape[1]
    rows = []
    for pos in range(seq_len):
        col = sub[:, pos]
        valid = col != 255
        alleles = col[valid]
        if alleles.size < 2:
            continue
        if np.unique(alleles).size >= 2:
            rows.append(alleles.astype(np.int8))
    if not rows:
        logging.warning(f'{fname}：无多态位点')
        return np.zeros((0, n), dtype=np.int8)

    mat = np.vstack(rows)
    logging.debug(f'{fname}：构建矩阵 {mat.shape[0]}×{mat.shape[1]} 完成')
    return mat


def rho_sigma_fast(A: np.ndarray, block_size: int = 20):
    """计算ρ和σ"""
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
    """处理单个任务"""
    hap, samples, fasta_dir, label, fname = args
    logging.info(f'[{hap}] 区域 {label} 开始：{fname}')
    G = build_int_matrix_cached(fasta_dir, fname, samples)
    rho, sigma = rho_sigma_fast(G)
    logging.info(f'[{hap}] 区域 {label} 完成：rho={rho:.4f}, sigma={sigma:.4f}')
    return hap, label, rho, sigma


def main():
    p = argparse.ArgumentParser(
        description='基于预加载 FASTA 的 mtDNA ρ/σ 计算（支持断点续跑）')
    p.add_argument('hap_csv',    help='两列 CSV：ID,Haplogroup')
    p.add_argument('fasta_dir',  help='存放 RHO_*.fasta 的目录')
    p.add_argument('output_tsv', help='最终 TSV 输出路径')
    p.add_argument('-t','--threads', type=int, default=1, help='并行进程数')
    p.add_argument('-v','--verbose', action='store_true', help='调试日志')
    args = p.parse_args()

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)s: %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    try:
        logging.info("=" * 60)
        logging.info("开始验证输入文件格式")
        logging.info("=" * 60)
        
        validate_csv_file(args.hap_csv)
        validate_fasta_files(args.fasta_dir)
        
        logging.info("=" * 60)
        logging.info("✓ 所有输入文件验证通过")
        logging.info("=" * 60)
        
    except (FileNotFoundError, NotADirectoryError, ValueError) as e:
        logging.error("\n" + "=" * 60)
        logging.error("❌ 输入文件验证失败")
        logging.error("=" * 60)
        logging.error(str(e))
        logging.error("=" * 60)
        sys.exit(1)

    hap2samples = {}
    with open(args.hap_csv) as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            hap2samples.setdefault(r['Haplogroup'], []).append(r['ID'])

    for h in list(hap2samples):
        if len(hap2samples[h]) < 2:
            logging.warning(f'跳过 {h}（样本数 {len(hap2samples[h])} < 2）')
            hap2samples.pop(h)

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

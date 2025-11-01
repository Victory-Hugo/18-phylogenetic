#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import csv
import argparse
import logging
import numpy as np
from cyvcf2 import VCF
from multiprocessing import Pool

# 预定义的五个区域及其对应 VCF 文件名
regions = [
    ('CS',         'RHO_Complete_sequence.vcf.gz'),
    ('HVSI',       'RHO_HVS-I.vcf.gz'),
    ('HVSI_TRANSI','RHO_HVS-I_transitions.vcf.gz'),
    ('HVSII',      'RHO_HVS-II.vcf.gz'),
    ('CR',         'RHO_Control_Region.vcf.gz')
]

def build_int_matrix(vcf_file, samples):
    """
    使用 cyvcf2 读取 VCF，构建 (sites x samples) 的 int8 矩阵：
      0 => REF, 1+ => ALT, -1 => 缺失或杂合
    samples 可以是一个字符串路径（每行一个样本名），也可以是一个样本名列表
    """
    logging.debug(f"加载样本列表，来自: {samples}")
    if isinstance(samples, str):
        with open(samples) as f:
            samples = [s.strip() for s in f if s.strip()]
    logging.debug(f"打开 VCF 文件: {vcf_file}")
    vcf = VCF(vcf_file)
    idx = [vcf.samples.index(s) for s in samples if s in vcf.samples]
    if len(idx) < 2:
        raise RuntimeError(f"样本不足 2 个: {vcf_file}")
    rows = []
    for var in vcf:
        gt = var.genotypes
        row = np.fromiter(
            (gt[i][0] if gt[i][0] == gt[i][1] and gt[i][0] >= 0 else -1
             for i in idx),
            dtype=np.int8, count=len(idx)
        )
        rows.append(row)
    vcf.close()
    logging.debug(f"构建矩阵完成: {len(rows)} 个位点 × {len(idx)} 个样本")
    return np.vstack(rows)

# def rho_sigma_fast(A: np.ndarray):
#     """
#     计算平均核苷酸差异 ρ 和 σ：
#       ρ = mean(Hamming distances)
#       σ = sqrt(sum(sq)/n^2)
#     """
#     logging.debug("开始计算 rho 和 sigma")
#     M = np.where(A < 0, 127, A)
#     D = (M[:, :, None] != M[:, None, :]).sum(axis=0)
#     tri = D[np.triu_indices_from(D, 1)]
#     rho   = tri.mean()
#     sigma = math.sqrt((tri.astype(np.int64)**2).sum() / (tri.size**2))
#     logging.debug(f"计算结果 -- rho: {rho}, sigma: {sigma}")
#     return rho, sigma
def rho_sigma_fast(A: np.ndarray, block_size: int = 20):
    """
    分块计算样本间海明距离矩阵 D，再从 D 中抽出上三角，
    最终计算 rho 和 sigma。

    参数
    ----
    A : np.ndarray, shape=(n_sites, n_samples)
      原始整数矩阵，-1（缺失）会被替换为 127。
    block_size : int
      每次处理多少行以控制内存。

    返回
    ----
    rho : float
    sigma : float
    """
    logging.debug(f"分块大小：{block_size} 行")
    # 将缺失值 A<0 替换成一个不会相等的常数
    M = np.where(A < 0, 127, A).astype(np.int8)

    n_sites, n_samples = M.shape
    # 用 int32 足够存放 site 数之和
    D = np.zeros((n_samples, n_samples), dtype=np.int32)

    # 按行分块累加每个子块的样本对差异数
    for start in range(0, n_sites, block_size):
        end = min(start + block_size, n_sites)
        Mb = M[start:end]                         # shape = (b, n)
        # 生成布尔差异：shape = (b, n, n)
        diff = Mb[:, :, None] != Mb[:, None, :]
        # 每对样本在这一块上的差异总和，累加到 D 中
        D += diff.sum(axis=0).astype(np.int32)
        # 手动释放中间变量
        del diff, Mb

    # 取上三角 (i<j) 部分
    iu = np.triu_indices(n_samples, k=1)
    tri = D[iu]

    # ρ = mean distance
    rho = tri.mean()
    # σ = sqrt( sum(d^2) / n_pairs^2 )
    sigma = math.sqrt((tri.astype(np.int64)**2).sum() / (tri.size**2))

    logging.debug(f"分块后计算结果 -- rho: {rho}, sigma: {sigma}")
    return rho, sigma

def compute_one(args):
    """
    对单个 (haplogroup, region) 任务执行计算，返回 (hap, label, rho, sigma)
    args = (hap, samples, vcf_dir, label, fname)
    """
    hap, samples, vcf_dir, label, fname = args
    path = os.path.join(vcf_dir, fname)
    if not os.path.isfile(path):
        logging.error(f"未找到 VCF 文件 {path}")
        sys.exit(1)
    logging.info(f"[{hap}] 区域 {label}: 读取 {fname}")
    G = build_int_matrix(path, samples)
    rho, sigma = rho_sigma_fast(G)
    logging.info(f"[{hap}] 区域 {label} 完成: rho={rho:.4f}, sigma={sigma:.4f}")
    return hap, label, rho, sigma

def main():
    parser = argparse.ArgumentParser(
        description="按单倍群×区域并行计算 ρ (rho) 和 σ，输出汇总 TSV。"
    )
    parser.add_argument('hap_csv',
        help='CSV 文件，第一列 ID，第二列 Haplogroup，无多余列')
    parser.add_argument('vcf_dir',
        help='存放 RHO_*.vcf.gz 文件的目录')
    parser.add_argument('output_tsv',
        help='最终输出的 TSV 文件路径')
    parser.add_argument('--threads', '-t', type=int, default=1,
        help='并行进程数（默认 1）')
    parser.add_argument('--verbose', '-v', action='store_true',
        help='显示调试日志')
    args = parser.parse_args()

    # 配置日志
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format='[%(asctime)s] %(levelname)s: %(message)s',
        level=level,
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info("脚本启动")

    # 读取 CSV 并按 Haplogroup 分组
    logging.info(f"读取 CSV: {args.hap_csv}")
    hap2samples = {}
    with open(args.hap_csv) as f:
        reader = csv.DictReader(f)
        if 'ID' not in reader.fieldnames or 'Haplogroup' not in reader.fieldnames:
            logging.error("CSV 必须包含 'ID' 和 'Haplogroup' 两列")
            sys.exit(1)
        for row in reader:
            hid = row['Haplogroup']
            sid = row['ID']
            hap2samples.setdefault(hid, []).append(sid)

    # 构建扁平化任务列表
    tasks = [
        (hap, samples, args.vcf_dir, label, fname)
        for hap, samples in hap2samples.items()
        for label, fname in regions
    ]

    logging.info(f"共 {len(tasks)} 个 (hap,region) 任务，使用 {args.threads} 个进程并行")
    with Pool(args.threads) as pool:
        all_results = pool.map(compute_one, tasks)

    # 汇总结果：hap -> { region: (rho,sigma) }
    summary = {}
    for hap, label, rho, sigma in all_results:
        summary.setdefault(hap, {})[label] = (rho, sigma)

    # 写入 TSV
    header = ['Haplogroup']
    for label, _ in regions:
        header += [f"RHO_{label}", f"RHO_{label}_SE"]
    logging.info(f"写入输出文件: {args.output_tsv}")
    with open(args.output_tsv, 'w') as out:
        out.write("\t".join(header) + "\n")
        for hap in sorted(summary):
            row = [hap]
            for label, _ in regions:
                rho, sigma = summary[hap][label]
                row += [str(rho), str(sigma)]
            out.write("\t".join(row) + "\n")

    logging.info("全部完成 ✔️")

if __name__ == '__main__':
    main()


# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# import os
# import sys
# import math
# import csv
# import argparse
# import logging
# import numpy as np
# from cyvcf2 import VCF
# from multiprocessing import Pool

# # 预定义的五个区域及其对应 VCF 文件名
# regions = [
#     ('CS',         'RHO_Complete_sequence.vcf.gz'),
#     ('HVSI',       'RHO_HVS-I.vcf.gz'),
#     ('HVSI_TRANSI','RHO_HVS-I_transitions.vcf.gz'),
#     ('HVSII',      'RHO_HVS-II.vcf.gz'),
#     ('CR',         'RHO_Control_Region.vcf.gz')
# ]

# def build_int_matrix(vcf_file, samples):
#     logging.debug(f"加载样本列表，来自: {samples}")
#     if isinstance(samples, str):
#         with open(samples) as f:
#             samples = [s.strip() for s in f if s.strip()]
#     vcf = VCF(vcf_file)
#     idx = [vcf.samples.index(s) for s in samples if s in vcf.samples]
#     n = len(idx)
#     if n < 2:
#         logging.warning(f"{vcf_file} 样本数 {n} < 2，跳过构建矩阵")
#         vcf.close()
#         return np.zeros((0, n), dtype=np.int8)
#     rows = []
#     for var in vcf:
#         gt = var.genotypes
#         row = np.fromiter(
#             (gt[i][0] if (gt[i][0] == gt[i][1] and gt[i][0] >= 0) else -1
#              for i in idx),
#             dtype=np.int8, count=n
#         )
#         rows.append(row)
#     vcf.close()
#     if not rows:
#         logging.warning(f"{vcf_file} 无变异位点，返回 0×{n} 矩阵")
#         return np.zeros((0, n), dtype=np.int8)
#     logging.debug(f"构建矩阵完成: {len(rows)} 个位点 × {n} 个样本")
#     return np.vstack(rows)

# def rho_sigma_fast(A: np.ndarray, block_size: int = 20):
#     if A.shape[1] < 2 or A.shape[0] < 1:
#         logging.warning(f"A.shape={A.shape}，样本或位点不足，返回 NaN")
#         return float('nan'), float('nan')
#     M = np.where(A < 0, 127, A).astype(np.int8)
#     n_sites, n_samples = M.shape
#     D = np.zeros((n_samples, n_samples), dtype=np.int32)
#     for start in range(0, n_sites, block_size):
#         end = min(start + block_size, n_sites)
#         Mb = M[start:end]
#         diff = Mb[:, :, None] != Mb[:, None, :]
#         D += diff.sum(axis=0).astype(np.int32)
#         del diff, Mb
#     iu = np.triu_indices(n_samples, k=1)
#     tri = D[iu]
#     rho   = tri.mean()
#     sigma = math.sqrt((tri.astype(np.int64)**2).sum() / (tri.size**2))
#     logging.debug(f"计算结果 -- rho: {rho:.4f}, sigma: {sigma:.4f}")
#     return rho, sigma

# def compute_one(args):
#     hap, samples, vcf_dir, label, fname = args
#     path = os.path.join(vcf_dir, fname)
#     if not os.path.isfile(path):
#         logging.error(f"未找到 VCF 文件 {path}")
#         return hap, label, float('nan'), float('nan')
#     logging.info(f"[{hap}] 区域 {label}: 读取 {fname}")
#     G = build_int_matrix(path, samples)
#     rho, sigma = rho_sigma_fast(G)
#     logging.info(f"[{hap}] 区域 {label} 完成: rho={rho:.4f}, sigma={sigma:.4f}")
#     return hap, label, rho, sigma

# def main():
#     parser = argparse.ArgumentParser(
#         description="按单倍群×区域并行计算 ρ (rho) 和 σ，支持断点续跑，输出汇总 TSV。"
#     )
#     parser.add_argument('hap_csv',
#                         help='CSV 文件，第一列 ID，第二列 Haplogroup，无多余列')
#     parser.add_argument('vcf_dir',
#                         help='存放 RHO_*.vcf.gz 文件的目录')
#     parser.add_argument('output_tsv',
#                         help='最终输出的 TSV 文件路径')
#     parser.add_argument('--threads', '-t', type=int, default=1,
#                         help='并行进程数（默认 1）')
#     parser.add_argument('--verbose', '-v', action='store_true',
#                         help='显示调试日志')
#     args = parser.parse_args()

#     # 日志配置
#     level = logging.DEBUG if args.verbose else logging.INFO
#     logging.basicConfig(
#         format='[%(asctime)s] %(levelname)s: %(message)s',
#         level=level, datefmt='%Y-%m-%d %H:%M:%S'
#     )
#     logging.info("脚本启动")

#     # 读取 CSV 并按 Haplogroup 分组
#     logging.info(f"读取 CSV: {args.hap_csv}")
#     hap2samples = {}
#     with open(args.hap_csv) as f:
#         reader = csv.DictReader(f)
#         if 'ID' not in reader.fieldnames or 'Haplogroup' not in reader.fieldnames:
#             logging.error("CSV 必须包含 'ID' 和 'Haplogroup' 两列")
#             sys.exit(1)
#         for row in reader:
#             hap2samples.setdefault(row['Haplogroup'], []).append(row['ID'])

#     # 过滤样本数不足 2 的单倍群
#     skipped = [hap for hap, smp in hap2samples.items() if len(smp) < 2]
#     for hap in skipped:
#         logging.warning(f"跳过 {hap}：样本数 {len(hap2samples[hap])} < 2")
#     hap2samples = {hap: smp for hap, smp in hap2samples.items() if len(smp) >= 2}

#     # 准备断点续跑
#     checkpoint = args.output_tsv + '.ckpt'
#     done = set()
#     if os.path.isfile(checkpoint):
#         logging.info(f"加载断点文件: {checkpoint}")
#         with open(checkpoint) as ck:
#             rdr = csv.DictReader(ck, delimiter='\t')
#             for r in rdr:
#                 done.add((r['Haplogroup'], r['Region']))
#     else:
#         with open(checkpoint, 'w') as ck:
#             writer = csv.writer(ck, delimiter='\t')
#             writer.writerow(['Haplogroup', 'Region', 'RHO', 'RHO_SE'])

#     # 构建待计算任务，跳过已完成
#     tasks = []
#     for hap, samples in hap2samples.items():
#         for label, fname in regions:
#             if (hap, label) in done:
#                 logging.debug(f"已完成 {(hap,label)}，跳过")
#                 continue
#             tasks.append((hap, samples, args.vcf_dir, label, fname))
#     logging.info(f"待执行任务数: {len(tasks)}")

#     # 并行计算并实时写入断点
#     if tasks:
#         with Pool(args.threads) as pool, open(checkpoint, 'a') as ck:
#             writer = csv.writer(ck, delimiter='\t')
#             for hap, label, rho, sigma in pool.imap_unordered(compute_one, tasks):
#                 writer.writerow([hap, label, f"{rho:.6g}", f"{sigma:.6g}"])
#                 ck.flush()
#                 done.add((hap, label))

#     # 从断点文件汇总到最终 TSV
#     summary = {}
#     with open(checkpoint) as ck:
#         rdr = csv.DictReader(ck, delimiter='\t')
#         for r in rdr:
#             summary.setdefault(r['Haplogroup'], {})[r['Region']] = (
#                 float(r['RHO']), float(r['RHO_SE'])
#             )

#     logging.info(f"写入汇总文件: {args.output_tsv}")
#     header = ['Haplogroup']
#     for label, _ in regions:
#         header += [f"RHO_{label}", f"RHO_{label}_SE"]
#     with open(args.output_tsv, 'w') as out:
#         out.write('\t'.join(header) + '\n')
#         for hap in sorted(summary):
#             row = [hap]
#             for label, _ in regions:
#                 rho, sigma = summary.get(hap, {}).get(label, (float('nan'), float('nan')))
#                 row += [str(rho), str(sigma)]
#             out.write('\t'.join(row) + '\n')

#     logging.info("全部完成 ✔️")

# if __name__ == '__main__':
#     main()

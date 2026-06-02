"""
s05_calc_rho_dating.py
对所有单倍群节点计算ρ统计量，并套用Soares 2009区域公式输出TMRCA估计。

方法：
    ρ = Haplogroup_H 到所有后代现代样本的平均有效突变距离
    SE = Saillard-style 标准误
    age = Soares 2009 区域校正公式

用法（CLI）：
    python s05_calc_rho_dating.py \
        --geno-matrix output/stage_1/intermediate/geno_matrix.npz \
        --geno-positions output/stage_1/intermediate/geno_positions.npy \
        --geno-samples output/stage_1/intermediate/geno_samples.txt \
        --site-masks output/stage_1/intermediate/site_masks.json \
        --hap-desc output/stage_1/intermediate/haplogroup_to_descendants.json \
        --output output/stage_1/results/rho_dating.tsv

用法（import）：
    from s05_calc_rho_dating import run
    run(...)
"""

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

log = logging.getLogger(__name__)

ANCESTOR_PREFIX = "Haplogroup_"

# 所有区域名及对应Soares线性换算系数（complete用非线性公式）
REGIONS = ["complete", "synonymous", "hvsi_full", "hvsi_trans", "hvsii", "control"]
LINEAR_RATES = {
    "synonymous": 7872,
    "hvsi_full":  16677,
    "hvsi_trans": 18845,
    "hvsii":      22388,
    "control":    9058,
}


# ── Soares 2009 年龄公式 ─────────────────────────────────────────────────

def _soares_complete(x: float) -> float:
    """Complete sequence 非线性校正公式（Soares 2009）。"""
    return 3624.0 * x * math.exp(-math.exp((x + 40.2789) * -0.0263))


def _soares_age(rho: float, se: float, region: str) -> dict:
    """
    计算给定区域的年龄和95%CI。
    CI计算：对ρ±1.96×SE代入公式（不截断负值）。
    """
    if region == "complete":
        F = _soares_complete
    else:
        rate = LINEAR_RATES[region]
        F = lambda x: x * rate

    lo = rho - 1.96 * se
    hi = rho + 1.96 * se
    return {
        "age_years":       F(rho),
        "ci95_lower_years": F(lo),
        "ci95_upper_years": F(hi),
        "age_kya":         F(rho) / 1000.0,
        "ci95_lower_kya":  F(lo) / 1000.0,
        "ci95_upper_kya":  F(hi) / 1000.0,
    }


# ── 置信度和异常标记 ─────────────────────────────────────────────────────

def _confidence(n: int) -> str:
    if n < 2:
        return "no_ci"
    if n < 3:
        return "very_low"
    if n < 10:
        return "low"
    if n < 50:
        return "medium"
    return "high"


def _flags(n: int, age_kya: float, ci_lo_kya: float, has_ancestor: bool) -> list[str]:
    flags = []
    if not has_ancestor:
        flags.append("no_ancestor")
    if ci_lo_kya < 0:
        flags.append("negative_lower_ci")
    if age_kya > 300:
        flags.append("exceeds_300kya")
    return flags


# ── 核心距离计算（向量化）────────────────────────────────────────────────

def _calc_distances_vectorized(
    anc_vec: np.ndarray,          # shape=(n_sites,)，祖先基因型
    desc_matrix: np.ndarray,      # shape=(n_desc, n_sites)，后代基因型
) -> np.ndarray:
    """
    批量计算祖先到每个后代的有效突变距离。
    有效位点：祖先和后代均非 -1（缺失）。
    距离 = 有效位点中基因型不同的位点数。

    返回 shape=(n_desc,) 的距离数组。
    """
    # 有效掩码：双方均非缺失（-1）
    valid_mask = (anc_vec != -1) & (desc_matrix != -1)  # shape=(n_desc, n_sites)
    # 不同碱基的位点
    diff_mask  = (anc_vec != desc_matrix) & valid_mask   # shape=(n_desc, n_sites)
    distances  = diff_mask.sum(axis=1)                   # shape=(n_desc,)
    return distances.astype(np.float64)


def run(
    geno_matrix: str,
    geno_positions: str,
    geno_samples: str,
    site_masks: str,
    hap_desc: str,
    output: str,
) -> int:
    """
    对所有单倍群计算ρ dating并输出TSV。

    参数：
        geno_matrix:    s04输出的.npz矩阵路径
        geno_positions: s04输出的POS数组路径
        geno_samples:   s04输出的样本列表路径
        site_masks:     s01输出的site_masks.json路径
        hap_desc:       s03输出的haplogroup_to_descendants.json路径
        output:         输出TSV路径

    返回：
        0 表示成功
    """
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    # ── 加载数据 ──────────────────────────────────────────────────────────
    log.info("加载基因型矩阵...")
    data = np.load(geno_matrix)
    matrix = data["matrix"]                     # shape=(n_samples, n_sites)
    positions = np.load(geno_positions)         # shape=(n_sites,)，1-based POS

    with open(geno_samples) as f:
        samples = [l.strip() for l in f if l.strip()]
    sample_index = {s: i for i, s in enumerate(samples)}
    log.info(f"  矩阵: {matrix.shape[0]} 样本 × {matrix.shape[1]} SNP位点")

    with open(site_masks) as f:
        masks = json.load(f)
    # 构建每个区域的列索引（在矩阵中的列号）
    pos_to_col = {int(p): idx for idx, p in enumerate(positions)}
    region_cols = {}
    for region in REGIONS:
        cols = [pos_to_col[p] for p in masks[region] if p in pos_to_col]
        region_cols[region] = np.array(cols, dtype=np.int32)
        log.info(f"  {region}: {len(cols)} 个有效列")

    with open(hap_desc) as f:
        hap_desc_map = json.load(f)
    log.info(f"  单倍群后代映射: {len(hap_desc_map)} 个节点")

    # hvsi_trans需要transition过滤：位点已在s01筛选，矩阵列已是trans位点，无需再处理

    # ── 主计算循环 ────────────────────────────────────────────────────────
    rows = []
    haplogroups = sorted(hap_desc_map.keys())

    for hap in tqdm(haplogroups, desc="计算ρ dating"):
        anc_id = f"{ANCESTOR_PREFIX}{hap}"
        has_ancestor = anc_id in sample_index

        if not has_ancestor:
            # 没有祖先节点：写NA行
            for region in REGIONS:
                rows.append(_na_row(hap, region, 0, "no_ancestor"))
            continue

        desc_sample_ids = hap_desc_map[hap]
        # 过滤出在VCF中实际存在的后代样本
        desc_indices = [sample_index[s] for s in desc_sample_ids if s in sample_index]
        n = len(desc_indices)

        if n == 0:
            for region in REGIONS:
                rows.append(_na_row(hap, region, 0, "no_descendants"))
            continue

        anc_row = matrix[sample_index[anc_id]]   # shape=(n_sites,)

        for region in REGIONS:
            cols = region_cols[region]
            if len(cols) == 0:
                rows.append(_na_row(hap, region, n, "no_region_sites"))
                continue

            anc_vec = anc_row[cols]               # shape=(n_region_sites,)
            desc_mat = matrix[np.array(desc_indices)][:, cols]  # shape=(n, n_region_sites)

            distances = _calc_distances_vectorized(anc_vec, desc_mat)

            rho = float(np.mean(distances))
            conf = _confidence(n)

            if n >= 2:
                se = float(np.sqrt(
                    np.sum((distances - rho) ** 2) / (n * (n - 1))
                ))
                ages = _soares_age(rho, se, region)
            else:
                # n=1：计算年龄但无SE/CI
                se = float("nan")
                ages = {
                    "age_years":        _soares_complete(rho) if region == "complete" else rho * LINEAR_RATES.get(region, 0),
                    "ci95_lower_years": float("nan"),
                    "ci95_upper_years": float("nan"),
                    "age_kya":          (_soares_complete(rho) if region == "complete" else rho * LINEAR_RATES.get(region, 0)) / 1000.0,
                    "ci95_lower_kya":   float("nan"),
                    "ci95_upper_kya":   float("nan"),
                }

            # 异常标记
            flag_list = _flags(
                n=n,
                age_kya=ages["age_kya"],
                ci_lo_kya=ages["ci95_lower_kya"] if not math.isnan(ages["ci95_lower_kya"]) else 0.0,
                has_ancestor=True,
            )

            rows.append({
                "haplogroup":       hap,
                "region":           region,
                "n_samples":        n,
                "rho":              round(rho, 6),
                "se":               round(se, 6) if not math.isnan(se) else "NA",
                "age_years":        round(ages["age_years"], 1),
                "ci95_lower_years": round(ages["ci95_lower_years"], 1) if not math.isnan(ages["ci95_lower_years"]) else "NA",
                "ci95_upper_years": round(ages["ci95_upper_years"], 1) if not math.isnan(ages["ci95_upper_years"]) else "NA",
                "age_kya":          round(ages["age_kya"], 4),
                "ci95_lower_kya":   round(ages["ci95_lower_kya"], 4) if not math.isnan(ages["ci95_lower_kya"]) else "NA",
                "ci95_upper_kya":   round(ages["ci95_upper_kya"], 4) if not math.isnan(ages["ci95_upper_kya"]) else "NA",
                "confidence":       conf,
                "flag":             ";".join(flag_list) if flag_list else "",
            })

    # ── 写出结果 ──────────────────────────────────────────────────────────
    df = pd.DataFrame(rows)
    df.to_csv(output, sep="\t", index=False)
    log.info(f"ρ dating结果 -> {output} ({len(df)} 行)")
    return 0


def _na_row(hap: str, region: str, n: int, flag: str) -> dict:
    """生成全NA结果行。"""
    return {
        "haplogroup": hap, "region": region, "n_samples": n,
        "rho": "NA", "se": "NA",
        "age_years": "NA", "ci95_lower_years": "NA", "ci95_upper_years": "NA",
        "age_kya": "NA", "ci95_lower_kya": "NA", "ci95_upper_kya": "NA",
        "confidence": "no_ci", "flag": flag,
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="计算ρ dating（s05）")
    p.add_argument("--geno-matrix",    required=True, help="基因型矩阵.npz路径")
    p.add_argument("--geno-positions", required=True, help="POS数组.npy路径")
    p.add_argument("--geno-samples",   required=True, help="样本列表.txt路径")
    p.add_argument("--site-masks",     required=True, help="site_masks.json路径")
    p.add_argument("--hap-desc",       required=True, help="haplogroup_to_descendants.json路径")
    p.add_argument("--output",         required=True, help="输出TSV路径")
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
        geno_matrix=args.geno_matrix,
        geno_positions=args.geno_positions,
        geno_samples=args.geno_samples,
        site_masks=args.site_masks,
        hap_desc=args.hap_desc,
        output=args.output,
    )


if __name__ == "__main__":
    raise SystemExit(main())

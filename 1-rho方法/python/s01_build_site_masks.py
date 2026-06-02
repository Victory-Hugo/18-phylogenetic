"""
s01_build_site_masks.py
构建各区域位点掩码（complete/synonymous/hvsi_full/hvsi_trans/hvsii/control）。

用法（CLI）：
    python s01_build_site_masks.py \
        --vcf data/merged_clean.vcf.gz \
        --back-mutation META/1-回复突变.tsv \
        --output output/stage_1/intermediate/site_masks.json \
        --bcftools /path/to/bcftools

用法（import）：
    from s01_build_site_masks import run
    run(vcf="...", back_mutation="...", output="...", bcftools="...")
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

# ── 区域坐标常量（rCRS 1-based）─────────────────────────────────────────
HVSI_FULL_START    = 16051
HVSI_FULL_END      = 16400
HVSI_TRANS_START   = 16090
HVSI_TRANS_END     = 16365
HVSII_START        = 68
HVSII_END          = 263
CONTROL_RANGES     = [(16024, 16569), (1, 576)]  # 控制区两段（环状）

# 需要在complete/hvsi中排除的位点（整个位点保守排除）
EXCLUDE_COMPLETE   = {16182, 16183, 16194, 16519}
EXCLUDE_HVSI       = {16182, 16183, 16194}

# Transition碱基对：只要REF/ALT属于以下任一对，则为transition
TRANSITION_PAIRS   = {frozenset("AG"), frozenset("CT")}

# 编码区基因名（Location字段值），用于筛选synonymous
CODING_GENES = {
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "CO1", "CO2", "CO3",
    "CYTB",
    "ATP6", "ATP8",
    # ATP8/6重叠区保守处理：只计入纯S，不计S/NS或NS/S
}


def _run_bcftools_query(vcf: str, bcftools: str) -> list[dict]:
    """用bcftools query提取所有位点的POS、TYPE、REF、ALT信息。"""
    cmd = [
        bcftools, "query",
        "-f", "%POS\t%TYPE\t%REF\t%ALT\n",
        vcf,
    ]
    log.info("运行bcftools query提取位点信息...")
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    records = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        pos, typ, ref, alt_str = int(parts[0]), parts[1], parts[2], parts[3]
        alts = alt_str.split(",")  # 多等位位点有多个ALT
        records.append({"pos": pos, "type": typ, "ref": ref, "alts": alts})
    log.info(f"  共读取 {len(records)} 个位点")
    return records


def _is_transition(ref: str, alt: str) -> bool:
    """判断一个碱基替换是否为transition（A↔G 或 C↔T）。"""
    if len(ref) != 1 or len(alt) != 1:
        return False
    return frozenset([ref.upper(), alt.upper()]) in TRANSITION_PAIRS


def _load_back_mutation(back_mutation_tsv: str) -> dict[int, dict]:
    """
    读取回复突变表，返回 pos -> {location, aa_change} 映射。
    只保留Location为编码区基因的行，用于synonymous筛选。
    """
    df = pd.read_csv(back_mutation_tsv, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()
    # 列名：Location, Aamino_acid_change, Position, Number_of_occurrences
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df = df.dropna(subset=["Position"])
    df["Position"] = df["Position"].astype(int)

    pos_map = {}
    for _, row in df.iterrows():
        loc = str(row.get("Location", "")).strip()
        aac = str(row.get("Aamino_acid_change", "")).strip()
        pos = int(row["Position"])
        # 同一位点可能多行，取第一行（主要信息）
        if pos not in pos_map:
            pos_map[pos] = {"location": loc, "aa_change": aac}
    log.info(f"  回复突变表共 {len(pos_map)} 个位点注释")
    return pos_map


def _classify_sites(records: list[dict], pos_map: dict[int, dict]) -> dict:
    """
    对所有VCF位点按区域分类，返回各区域的POS列表。

    分类规则：
    - all_snp_positions: TYPE=SNP 的所有位点
    - indel_positions: TYPE=INDEL 的所有位点
    - complete: SNP，且不在EXCLUDE_COMPLETE中
    - synonymous: SNP，Location为编码基因，aa_change == 'S'（严格同义，排除S/NS等重叠）
    - hvsi_full: SNP，POS in [16051,16400]，不在EXCLUDE_HVSI中
    - hvsi_trans: SNP，POS in [16090,16365]，且所有ALT均为transition（多等位：有任一ALT为transition则纳入）
    - hvsii: SNP，POS in [68,263]
    - control: SNP，POS in [16024,16569] 或 [1,576]
    """
    all_snp, indels = [], []
    complete, synonymous = [], []
    hvsi_full, hvsi_trans = [], []
    hvsii, control = [], []

    for rec in records:
        pos, typ, ref, alts = rec["pos"], rec["type"], rec["ref"], rec["alts"]

        if typ == "INDEL":
            indels.append(pos)
            continue
        if typ != "SNP":
            # 其他类型（如MNP）跳过
            continue

        all_snp.append(pos)

        # ── complete ──────────────────────────────────────────────────
        if pos not in EXCLUDE_COMPLETE:
            complete.append(pos)

        # ── synonymous ────────────────────────────────────────────────
        ann = pos_map.get(pos)
        if ann:
            loc = ann["location"]
            aac = ann["aa_change"]
            # 只有Location是编码基因且aa_change严格为"S"时才计入
            if loc in CODING_GENES and aac == "S":
                synonymous.append(pos)

        # ── hvsi_full ─────────────────────────────────────────────────
        if HVSI_FULL_START <= pos <= HVSI_FULL_END and pos not in EXCLUDE_HVSI:
            hvsi_full.append(pos)

        # ── hvsi_trans ────────────────────────────────────────────────
        if HVSI_TRANS_START <= pos <= HVSI_TRANS_END:
            # 多等位位点：只要有至少一个ALT是transition就纳入
            if any(_is_transition(ref, alt) for alt in alts if len(alt) == 1):
                hvsi_trans.append(pos)

        # ── hvsii ─────────────────────────────────────────────────────
        if HVSII_START <= pos <= HVSII_END:
            hvsii.append(pos)

        # ── control region ────────────────────────────────────────────
        in_control = any(
            lo <= pos <= hi for lo, hi in CONTROL_RANGES
        )
        if in_control:
            control.append(pos)

    return {
        "all_snp_positions":  sorted(set(all_snp)),
        "indel_positions":    sorted(set(indels)),
        "complete":           sorted(set(complete)),
        "synonymous":         sorted(set(synonymous)),
        "hvsi_full":          sorted(set(hvsi_full)),
        "hvsi_trans":         sorted(set(hvsi_trans)),
        "hvsii":              sorted(set(hvsii)),
        "control":            sorted(set(control)),
    }


def run(
    vcf: str,
    back_mutation: str,
    output: str,
    bcftools: str = "bcftools",
) -> int:
    """
    构建各区域位点掩码并写出JSON。

    参数：
        vcf:          VCF文件路径
        back_mutation: 回复突变TSV路径
        output:       输出JSON路径
        bcftools:     bcftools可执行文件路径

    返回：
        0 表示成功
    """
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    # 提取VCF位点信息
    records = _run_bcftools_query(vcf, bcftools)

    # 加载位点注释
    pos_map = _load_back_mutation(back_mutation)

    # 分类
    masks = _classify_sites(records, pos_map)

    # 输出统计
    for region, positions in masks.items():
        log.info(f"  {region}: {len(positions)} 个位点")

    # 写出JSON
    with open(output, "w") as f:
        json.dump(masks, f, indent=2)
    log.info(f"位点掩码已写出至: {output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="构建各区域位点掩码（s01）"
    )
    p.add_argument("--vcf",           required=True, help="输入VCF路径")
    p.add_argument("--back-mutation", required=True, help="回复突变TSV路径")
    p.add_argument("--output",        required=True, help="输出JSON路径")
    p.add_argument("--bcftools",      default="bcftools", help="bcftools路径")
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
        vcf=args.vcf,
        back_mutation=args.back_mutation,
        output=args.output,
        bcftools=args.bcftools,
    )


if __name__ == "__main__":
    raise SystemExit(main())

"""
s07_merge_output.py
将ρ dating长格式结果转为宽格式（每行一个单倍群），同时添加汇总统计。

用法（CLI）：
    python s07_merge_output.py \
        --rho-dating output/stage_1/results/rho_dating.tsv \
        --output-wide output/stage_1/results/rho_dating_wide.tsv \
        --output-summary output/stage_1/results/rho_dating_summary.tsv

用法（import）：
    from s07_merge_output import run
    run(...)
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

REGIONS = ["complete", "synonymous", "hvsi_full", "hvsi_trans", "hvsii", "control"]


def run(
    rho_dating: str,
    output_wide: str,
    output_summary: str,
) -> int:
    """
    转换ρ dating结果为宽格式，并生成汇总统计。

    参数：
        rho_dating:     s05输出的rho_dating.tsv路径
        output_wide:    输出宽格式TSV路径
        output_summary: 输出汇总统计TSV路径

    返回：
        0 表示成功
    """
    for path in [output_wide, output_summary]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log.info(f"读取ρ dating结果: {rho_dating}")
    df = pd.read_csv(rho_dating, sep="\t", dtype=str)
    log.info(f"  共 {len(df)} 行（{df['haplogroup'].nunique()} 个单倍群 × {df['region'].nunique()} 个区域）")

    # ── 转为宽格式 ────────────────────────────────────────────────────────
    # 每行一个单倍群，每个区域的age_kya/ci95_lower_kya/ci95_upper_kya作为独立列
    numeric_cols = ["rho", "se", "age_years", "ci95_lower_years", "ci95_upper_years",
                    "age_kya", "ci95_lower_kya", "ci95_upper_kya"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["n_samples"] = pd.to_numeric(df["n_samples"], errors="coerce")

    wide_rows = []
    for hap, grp in df.groupby("haplogroup", sort=False):
        row = {"haplogroup": hap}
        grp_dict = {r: g.iloc[0] for r, g in grp.groupby("region")}

        # 取n_samples（以complete区域为准）
        complete_row = grp_dict.get("complete")
        row["n_samples"]  = int(complete_row["n_samples"]) if complete_row is not None and not pd.isna(complete_row["n_samples"]) else None
        row["confidence"] = complete_row["confidence"] if complete_row is not None else "NA"
        row["flag"]       = complete_row["flag"] if complete_row is not None else "NA"

        for region in REGIONS:
            r_row = grp_dict.get(region)
            pfx = region  # 列名前缀

            if r_row is None:
                row[f"{pfx}_rho"]        = None
                row[f"{pfx}_se"]         = None
                row[f"{pfx}_age_kya"]    = None
                row[f"{pfx}_ci_lo_kya"]  = None
                row[f"{pfx}_ci_hi_kya"]  = None
            else:
                row[f"{pfx}_rho"]       = r_row["rho"]
                row[f"{pfx}_se"]        = r_row["se"]
                row[f"{pfx}_age_kya"]   = r_row["age_kya"]
                row[f"{pfx}_ci_lo_kya"] = r_row["ci95_lower_kya"]
                row[f"{pfx}_ci_hi_kya"] = r_row["ci95_upper_kya"]

        wide_rows.append(row)

    wide_df = pd.DataFrame(wide_rows)
    # 按complete_age_kya排序（大→小，NA最后）
    wide_df = wide_df.sort_values("complete_age_kya", ascending=False, na_position="last")
    wide_df.to_csv(output_wide, sep="\t", index=False, float_format="%.4f")
    log.info(f"宽格式结果 -> {output_wide} ({len(wide_df)} 行)")

    # ── 汇总统计 ──────────────────────────────────────────────────────────
    summary_df = df[df["region"] == "complete"].copy()
    summary_df["n_samples"] = pd.to_numeric(summary_df["n_samples"], errors="coerce")
    summary_df["age_kya"]   = pd.to_numeric(summary_df["age_kya"],   errors="coerce")

    summary_stats = {
        "total_haplogroups":         summary_df["haplogroup"].nunique(),
        "with_ancestor":             int((summary_df["flag"].fillna("").str.contains("no_ancestor") == False).sum()),
        "with_descendants_n_ge_2":   int((summary_df["n_samples"] >= 2).sum()),
        "confidence_high":           int((summary_df["confidence"] == "high").sum()),
        "confidence_medium":         int((summary_df["confidence"] == "medium").sum()),
        "confidence_low":            int((summary_df["confidence"] == "low").sum()),
        "confidence_very_low":       int((summary_df["confidence"] == "very_low").sum()),
        "confidence_no_ci":          int((summary_df["confidence"] == "no_ci").sum()),
        "flag_no_ancestor":          int(summary_df["flag"].fillna("").str.contains("no_ancestor").sum()),
        "flag_negative_lower_ci":    int(summary_df["flag"].fillna("").str.contains("negative_lower_ci").sum()),
        "flag_exceeds_300kya":       int(summary_df["flag"].fillna("").str.contains("exceeds_300kya").sum()),
        "median_age_kya":            round(float(summary_df["age_kya"].median()), 2),
        "mean_age_kya":              round(float(summary_df["age_kya"].mean()), 2),
    }

    summary_rows = [{"metric": k, "value": v} for k, v in summary_stats.items()]
    pd.DataFrame(summary_rows).to_csv(output_summary, sep="\t", index=False)
    log.info(f"汇总统计 -> {output_summary}")

    # 控制台打印关键统计
    for k, v in summary_stats.items():
        log.info(f"  {k}: {v}")

    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="合并输出结果为宽格式（s07）")
    p.add_argument("--rho-dating",     required=True, help="rho_dating.tsv路径")
    p.add_argument("--output-wide",    required=True, help="宽格式输出TSV路径")
    p.add_argument("--output-summary", required=True, help="汇总统计TSV路径")
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
        rho_dating=args.rho_dating,
        output_wide=args.output_wide,
        output_summary=args.output_summary,
    )


if __name__ == "__main__":
    raise SystemExit(main())

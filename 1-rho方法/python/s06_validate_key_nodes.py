"""
s06_validate_key_nodes.py
验证关键单倍群节点的年龄估计是否符合文献预期范围。

用法（CLI）：
    python s06_validate_key_nodes.py \
        --rho-dating output/stage_1/results/rho_dating.tsv \
        --validate-nodes L3,M,N,R,H,D4,U \
        --output output/stage_1/results/validation_key_nodes.tsv

用法（import）：
    from s06_validate_key_nodes import run
    run(...)
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

# 关键节点预期年龄范围（kya），来自文献（主要Soares 2009及后续研究）
EXPECTED_RANGES_KYA = {
    "L3":  (60,  90),
    "M":   (55,  75),
    "N":   (55,  75),
    "R":   (50,  70),
    "H":   (12,  25),
    "D4":  (35,  55),
    "U":   (45,  65),
}


def run(
    rho_dating: str,
    validate_nodes: str,
    output: str,
) -> int:
    """
    验证关键节点年龄。

    参数：
        rho_dating:      s05输出的rho_dating.tsv路径
        validate_nodes:  逗号分隔的待验证节点名称（如"L3,M,N,R,H,D4,U"）
        output:          输出验证报告TSV路径

    返回：
        0 表示通过验证，1 表示有节点超出预期范围
    """
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    nodes = [n.strip() for n in validate_nodes.split(",") if n.strip()]
    log.info(f"验证节点: {nodes}")

    df = pd.read_csv(rho_dating, sep="\t", dtype=str)
    # 只看complete区域结果
    df_complete = df[df["region"] == "complete"].copy()
    df_complete["age_kya"] = pd.to_numeric(df_complete["age_kya"], errors="coerce")
    df_complete["ci95_lower_kya"] = pd.to_numeric(df_complete["ci95_lower_kya"], errors="coerce")
    df_complete["ci95_upper_kya"] = pd.to_numeric(df_complete["ci95_upper_kya"], errors="coerce")
    df_complete["n_samples"] = pd.to_numeric(df_complete["n_samples"], errors="coerce")

    rows = []
    all_pass = True

    for node in nodes:
        row_df = df_complete[df_complete["haplogroup"] == node]
        expected = EXPECTED_RANGES_KYA.get(node, None)

        if row_df.empty:
            status = "MISSING"
            all_pass = False
            rows.append({
                "haplogroup": node,
                "age_kya": "NA",
                "ci95_lower_kya": "NA",
                "ci95_upper_kya": "NA",
                "n_samples": "NA",
                "expected_min_kya": expected[0] if expected else "NA",
                "expected_max_kya": expected[1] if expected else "NA",
                "status": status,
                "note": "节点在结果中缺失",
            })
            log.warning(f"  {node}: 结果中缺失")
            continue

        row = row_df.iloc[0]
        age = row["age_kya"]
        lo  = row["ci95_lower_kya"]
        hi  = row["ci95_upper_kya"]
        n   = row["n_samples"]
        flag = str(row.get("flag", ""))
        conf = str(row.get("confidence", ""))

        if pd.isna(age):
            status = "NA_AGE"
            all_pass = False
            note = "age_kya为NA（可能无祖先节点或无后代样本）"
        elif expected:
            exp_lo, exp_hi = expected
            if exp_lo <= age <= exp_hi:
                status = "PASS"
                note = f"在预期范围[{exp_lo},{exp_hi}]内"
            else:
                status = "FAIL"
                all_pass = False
                note = f"超出预期范围[{exp_lo},{exp_hi}]，实际={age:.1f}kya"
        else:
            status = "NO_EXPECTED"
            note = "无预设预期范围"

        # 打印详细信息
        log.info(
            f"  {node}: age={age:.1f}kya CI=[{lo},{hi}] n={n} "
            f"conf={conf} flag={flag} → {status}"
        )

        rows.append({
            "haplogroup":       node,
            "age_kya":          round(age, 2) if not pd.isna(age) else "NA",
            "ci95_lower_kya":   round(lo, 2)  if not pd.isna(lo)  else "NA",
            "ci95_upper_kya":   round(hi, 2)  if not pd.isna(hi)  else "NA",
            "n_samples":        int(n) if not pd.isna(n) else "NA",
            "expected_min_kya": expected[0] if expected else "NA",
            "expected_max_kya": expected[1] if expected else "NA",
            "status":           status,
            "note":             note,
        })

    # 写出验证报告
    out_df = pd.DataFrame(rows)
    out_df.to_csv(output, sep="\t", index=False)
    log.info(f"验证报告 -> {output}")

    # 汇总
    pass_count = sum(1 for r in rows if r["status"] == "PASS")
    fail_count = sum(1 for r in rows if r["status"] == "FAIL")
    log.info(f"验证结果：PASS={pass_count}，FAIL={fail_count}，其他={len(rows)-pass_count-fail_count}")

    return 0 if all_pass else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="验证关键节点年龄（s06）")
    p.add_argument("--rho-dating",      required=True, help="rho_dating.tsv路径")
    p.add_argument("--validate-nodes",  required=True, help="逗号分隔的节点名称")
    p.add_argument("--output",          required=True, help="输出验证TSV路径")
    return p


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    parser = build_parser()
    args = parser.parse_args(argv)
    rc = run(
        rho_dating=args.rho_dating,
        validate_nodes=args.validate_nodes,
        output=args.output,
    )
    if rc != 0:
        log.warning("部分关键节点未通过验证，请检查结果！")
    return 0   # 验证失败不阻断pipeline，返回0继续后续步骤


if __name__ == "__main__":
    raise SystemExit(main())

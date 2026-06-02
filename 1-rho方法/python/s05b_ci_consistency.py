"""
s05b_ci_consistency.py
祖先-后代 CI 一致性检查：
  若某下游单倍群（基于全基因组 complete 区域）的 CI 上限超过任何直接/间接
  祖先单倍群的 CI 上限，或 CI 下限超过任何祖先的 CI 下限，则将该单倍群
  所有区域行的 confidence 降级为 low，并在 flag 中追加 ci_exceeds_ancestor。

支持两种输入格式：
  --mode rho  : rho_dating.tsv（长格式，多区域行），使用 complete 区域的
                ci95_upper_kya / ci95_lower_kya 做检查
  --mode ml   : ml_dating.tsv（宽格式，每行一个单倍群），使用
                ml_ci95_upper_kya / ml_ci95_lower_kya 做检查；
                若无 confidence/flag 列则自动添加

用法：
    python s05b_ci_consistency.py --mode rho \
        --input  output/stage_1/results/rho_dating.tsv \
        --hap-graph output/stage_1/intermediate/haplogroup_graph.json \
        --output output/stage_1/results/rho_dating.tsv

    python s05b_ci_consistency.py --mode ml \
        --input  output/stage_2/results/ml_dating.tsv \
        --hap-graph output/stage_1/intermediate/haplogroup_graph.json \
        --output output/stage_2/results/ml_dating.tsv
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

CONF_ORDER = ["no_ci", "very_low", "low", "medium", "high"]


# ── 工具函数 ─────────────────────────────────────────────────────────────────

def _load_graph(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def _ancestors(graph: dict, hap: str) -> list[str]:
    """返回 hap 的所有上游节点（从直接父节点到根），不含自身。"""
    result = []
    node = graph.get(hap, {})
    parent = node.get("parent")
    while parent:
        result.append(parent)
        node = graph.get(parent, {})
        parent = node.get("parent")
    return result


def _flagged_haps(hap_ci: dict[str, tuple[float, float]], graph: dict) -> set[str]:
    """
    找出 CI 与祖先不一致的单倍群集合。
    hap_ci: {haplogroup: (ci_upper_kya, ci_lower_kya)}
    规则：若 upper > 任意祖先 upper，或 lower > 任意祖先 lower → 标记。
    """
    flagged: set[str] = set()
    for hap, (upper, lower) in hap_ci.items():
        for anc in _ancestors(graph, hap):
            if anc not in hap_ci:
                continue
            anc_upper, anc_lower = hap_ci[anc]
            if upper > anc_upper or lower > anc_lower:
                flagged.add(hap)
                break
    return flagged


def _downgrade_conf(current: str) -> str:
    """将 confidence 降级至 low（若当前比 low 高）。"""
    if current not in CONF_ORDER:
        return "low"
    idx = CONF_ORDER.index(current)
    low_idx = CONF_ORDER.index("low")
    return "low" if idx > low_idx else current


def _append_flag(existing: str, new_flag: str) -> str:
    if pd.isna(existing) or str(existing).strip() == "" or str(existing) == "nan":
        return new_flag
    if new_flag in str(existing):
        return existing
    return f"{existing};{new_flag}"


# ── rho 模式 ─────────────────────────────────────────────────────────────────

def _run_rho(input_tsv: str, hap_graph: str, output: str) -> int:
    graph = _load_graph(hap_graph)

    df = pd.read_csv(input_tsv, sep="\t", dtype=str)
    for col in ["ci95_upper_kya", "ci95_lower_kya"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # 仅用 complete 区域 CI 做判断
    complete = df[df["region"] == "complete"].dropna(subset=["ci95_upper_kya", "ci95_lower_kya"])
    hap_ci: dict[str, tuple[float, float]] = {
        row["haplogroup"]: (row["ci95_upper_kya"], row["ci95_lower_kya"])
        for _, row in complete.iterrows()
    }

    flagged = _flagged_haps(hap_ci, graph)
    log.info(f"[rho 模式] 共 {len(hap_ci)} 个单倍群有 complete CI，"
             f"其中 {len(flagged)} 个 CI 超过祖先 → 降级为 low")

    def fix_row(row):
        if row["haplogroup"] not in flagged:
            return row
        row = row.copy()
        row["confidence"] = _downgrade_conf(str(row["confidence"]))
        row["flag"] = _append_flag(row["flag"], "ci_exceeds_ancestor")
        return row

    df = df.apply(fix_row, axis=1)
    # 恢复数值列格式
    for col in ["ci95_upper_kya", "ci95_lower_kya"]:
        df[col] = df[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "NA")

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)
    log.info(f"更新后 rho_dating → {output}")
    return 0


# ── ml 模式 ──────────────────────────────────────────────────────────────────

def _run_ml(input_tsv: str, hap_graph: str, output: str) -> int:
    graph = _load_graph(hap_graph)

    df = pd.read_csv(input_tsv, sep="\t", dtype=str)
    for col in ["ml_ci95_upper_kya", "ml_ci95_lower_kya"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    hap_ci: dict[str, tuple[float, float]] = {
        row["haplogroup"]: (row["ml_ci95_upper_kya"], row["ml_ci95_lower_kya"])
        for _, row in df.iterrows()
        if pd.notna(row["ml_ci95_upper_kya"]) and pd.notna(row["ml_ci95_lower_kya"])
    }

    flagged = _flagged_haps(hap_ci, graph)
    log.info(f"[ml 模式] 共 {len(hap_ci)} 个单倍群有 ML CI，"
             f"其中 {len(flagged)} 个 CI 超过祖先 → 降级为 low")

    # 若无 confidence/flag 列则初始化
    if "confidence" not in df.columns:
        df["confidence"] = "NA"
    if "flag" not in df.columns:
        df["flag"] = ""

    def fix_row(row):
        if row["haplogroup"] not in flagged:
            return row
        row = row.copy()
        row["confidence"] = _downgrade_conf(str(row["confidence"]))
        row["flag"] = _append_flag(row["flag"], "ci_exceeds_ancestor")
        return row

    df = df.apply(fix_row, axis=1)
    for col in ["ml_ci95_upper_kya", "ml_ci95_lower_kya"]:
        df[col] = df[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "NA")

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)
    log.info(f"更新后 ml_dating → {output}")
    return 0


# ── 入口 ─────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="祖先-后代 CI 一致性检查（s05b）")
    p.add_argument("--mode",      required=True, choices=["rho", "ml"],
                   help="rho：rho_dating 长格式；ml：ml_dating 宽格式")
    p.add_argument("--input",     required=True, help="输入 TSV 路径")
    p.add_argument("--hap-graph", required=True, help="haplogroup_graph.json 路径")
    p.add_argument("--output",    required=True, help="输出 TSV 路径（可与 --input 相同）")
    return p


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    args = build_parser().parse_args(argv)
    if args.mode == "rho":
        return _run_rho(args.input, args.hap_graph, args.output)
    return _run_ml(args.input, args.hap_graph, args.output)


if __name__ == "__main__":
    raise SystemExit(main())

"""
s10_ml_dating.py
解析 IQ-TREE 3 treefile，用「局部路径插值」进行分子钟校准，
输出全部 5,439 个单倍群的 ML 年龄及 CI，并与 ρ 结果比较。

校准方法（局部路径插值 Local Clock Interpolation）：
    对 ML 树中每个节点 V，找其最近高可信 ρ 祖先 A（n ≥ calib_min_n）
    和 V 子树中所有高可信 ρ 后代 D_i，用 ML 枝长做比例插值：

        f_i  = D_{A→V}^{ML} / D_{A→D_i}^{ML}
        age_V = median_i [ age_A − f_i × (age_A − age_D_i) ]

    CI 同比例线性传播。
    ML 枝长只用于确定相对位置，绝对时间由 ρ 锚点定标。

已知 IQ-TREE 行为处理：
    · IQ-TREE 删去根节点标签：重新以 L0 为外群定根
    · IQ-TREE 合并 1,059 个单子节点：用父子 ρ/ML 年龄中点插值

不使用 bootstrap（-b/-bb）原因：
    年龄 CI 主要来源是 ρ 的 CI，而非 ML 枝长不确定性；
    bootstrap 仅提供拓扑支持值，对固定拓扑无意义。

用法（CLI）：
    python s10_ml_dating.py \
        --treefile        output/stage_2/intermediate/ml/anc_tree.treefile \
        --map             output/stage_2/intermediate/ml/name_sanitize_map.tsv \
        --rho             output/stage_1/results/rho_dating.tsv \
        --hap-graph       output/stage_1/intermediate/haplogroup_graph.json \
        --root-name       "mt-MRCA(RSRS)" \
        --root-age        163927 \
        --calib-min-n     10 \
        --outgroup        "L0" \
        --output-dating     output/stage_2/results/ml_dating.tsv \
        --output-comparison output/stage_2/results/ml_vs_rho_comparison.tsv

用法（import）：
    from s10_ml_dating import run
    run(treefile=..., map_tsv=..., rho_tsv=..., hap_graph_json=...,
        root_name=..., root_age=..., calib_min_n=..., outgroup=...,
        out_dating=..., out_comparison=...)
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from ete3 import Tree

log = logging.getLogger(__name__)

ROOT_NAME_DEFAULT  = "mt-MRCA(RSRS)"
ROOT_AGE_DEFAULT   = 163927.0
CALIB_MIN_N_DEFAULT = 10
OUTGROUP_DEFAULT   = "L0"


# ── 辅助函数 ─────────────────────────────────────────────────────────────────

def _load_rho(rho_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(rho_tsv, sep="\t")
    return df[df["region"] == "complete"].copy()


def _build_calib_map(rho_complete: pd.DataFrame, calib_min_n: int,
                     root_name: str, root_age: float) -> dict:
    """构建 {haplogroup → {age_years, ci95_lower_years, ci95_upper_years}} 校准锚点表。"""
    cmap: dict[str, dict] = {}
    for _, row in rho_complete.iterrows():
        if row["n_samples"] >= calib_min_n:
            cmap[row["haplogroup"]] = {
                "age_years":        float(row["age_years"]),
                "ci95_lower_years": max(0.0, float(row["ci95_lower_years"])),
                "ci95_upper_years": float(row["ci95_upper_years"]),
                "n_samples":        int(row["n_samples"]),
            }
    # 根节点作为固定锚点
    cmap[root_name] = {
        "age_years":        root_age,
        "ci95_lower_years": root_age,
        "ci95_upper_years": root_age,
        "n_samples":        999999,
    }
    log.info("校准锚点：%d 个（n ≥ %d + 根节点）", len(cmap) - 1, calib_min_n)
    return cmap


def _reroot(t: Tree, outgroup_name: str, root_label: str) -> tuple[Tree, object]:
    """以 outgroup_name 为外群重新定根，给新根打上 root_label 标签。"""
    og_nodes = t.search_nodes(name=outgroup_name)
    if not og_nodes:
        log.warning("外群节点 '%s' 未找到，使用当前根", outgroup_name)
        root = t.get_tree_root()
        root.name = root_label
        return t, root
    old_root = t.get_tree_root()
    old_root_name = old_root.name          # e.g., 'A' (IQ-TREE artifact at root)
    old_root_child_names = {ch.name for ch in old_root.children}
    t.set_outgroup(og_nodes[0])
    root = t.get_tree_root()
    # ete3 reuses the old root object as the new root and creates a new unnamed
    # internal node to hold the old root's haplogroup position.  Find that unnamed
    # node (its children overlap with old root's original children) and restore label.
    if old_root_name:
        for node in t.traverse():
            if node.name == '' and not node.is_leaf():
                if {ch.name for ch in node.children} & old_root_child_names:
                    node.name = old_root_name
                    break
    root.name = root_label
    log.info("树已重新定根（外群=%s，根=%s）", outgroup_name, root_label)
    return t, root


def _interpolate(D_AV: float, D_AD: float,
                 age_A: float, age_D: float,
                 lo_A: float, hi_A: float,
                 lo_D: float, hi_D: float,
                 root_age: float) -> tuple[float, float, float]:
    """路径比例插值，返回 (age_V, ci_lower, ci_upper)，结果夹在 [0, root_age]。"""
    if D_AD <= 0 or D_AV < 0:
        return age_A, lo_A, hi_A
    f = min(1.0, D_AV / D_AD)
    age_V = age_A - f * (age_A - age_D)
    lo_V  = lo_A  - f * (lo_A  - lo_D)
    hi_V  = hi_A  - f * (hi_A  - hi_D)
    age_V = float(np.clip(age_V, 0.0, root_age))
    lo_V  = float(np.clip(lo_V,  0.0, root_age))
    hi_V  = float(np.clip(hi_V,  0.0, root_age))
    return age_V, lo_V, hi_V


# ── 主函数 ───────────────────────────────────────────────────────────────────

def run(
    treefile: str,
    map_tsv: str,
    rho_tsv: str,
    hap_graph_json: str,
    root_name: str,
    root_age: float,
    calib_min_n: int,
    outgroup: str,
    out_dating: str,
    out_comparison: str,
) -> None:
    treefile       = Path(treefile)
    map_tsv        = Path(map_tsv)
    rho_tsv        = Path(rho_tsv)
    hap_graph_json = Path(hap_graph_json)
    out_dating     = Path(out_dating)
    out_comparison = Path(out_comparison)
    out_dating.parent.mkdir(parents=True, exist_ok=True)
    out_comparison.parent.mkdir(parents=True, exist_ok=True)

    # ── 名称映射 ─────────────────────────────────────────────────────────────
    san_to_orig: dict[str, str] = {}
    orig_to_san: dict[str, str] = {}
    with open(map_tsv) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                san_to_orig[parts[1]] = parts[0]
                orig_to_san[parts[0]] = parts[1]

    # ── 单倍群图（单子节点处理用）────────────────────────────────────────────
    with open(hap_graph_json) as f:
        graph: dict = json.load(f)
    single_child = {k: (v["parent"], v["children"][0])
                    for k, v in graph.items() if len(v["children"]) == 1}
    log.info("单子节点总数：%d", len(single_child))

    # ── ρ 数据 ───────────────────────────────────────────────────────────────
    rho_complete = _load_rho(rho_tsv)
    rho_idx: dict[str, dict] = {
        row["haplogroup"]: row.to_dict()
        for _, row in rho_complete.iterrows()
    }
    calib_map = _build_calib_map(rho_complete, calib_min_n, root_name, root_age)

    # ── 加载 ML 树并重新定根 ─────────────────────────────────────────────────
    t, root_node = _reroot(
        Tree(str(treefile), format=1),
        outgroup_name=outgroup,
        root_label=orig_to_san.get(root_name, root_name),
    )
    root_san = root_node.name
    log.info("树节点总数（重定根后）：%d（叶 %d，内部 %d）",
             sum(1 for _ in t.traverse()),
             sum(1 for n in t.iter_leaves()),
             sum(1 for n in t.traverse() if not n.is_leaf()))

    # ── 预计算根到所有节点 ML 距离 ───────────────────────────────────────────
    dist_from_root: dict[int, float] = {}
    for node in t.traverse():
        dist_from_root[id(node)] = node.get_distance(root_node)

    # ── 后序遍历：为每个节点积累子树内校准叶节点（O(n)）────────────────────
    # 每个节点存 list of (D_root_to_L, calib_entry)
    calib_below: dict[int, list] = {}
    for node in t.traverse("postorder"):
        orig = san_to_orig.get(node.name, node.name)
        own = []
        if orig in calib_map and node != root_node:
            own = [(dist_from_root[id(node)], calib_map[orig])]
        if node.is_leaf():
            calib_below[id(node)] = own
        else:
            merged = own[:]
            for child in node.children:
                merged.extend(calib_below.get(id(child), []))
            calib_below[id(node)] = merged

    # ── 预序遍历：为每个节点记录最近校准祖先 ───────────────────────────────
    nearest_anc: dict[int, object] = {}  # node_id → nearest calibrated ancestor node
    for node in t.traverse("preorder"):
        orig = san_to_orig.get(node.name, node.name)
        if orig in calib_map or node == root_node:
            nearest_anc[id(node)] = node
        else:
            nearest_anc[id(node)] = nearest_anc.get(id(node.up)) if node.up else root_node

    # ── 主循环：对 ML 树中每个节点计算年龄 ──────────────────────────────────
    ml_result: dict[str, tuple] = {}   # orig_name → (age, ci_lo, ci_hi, method)

    for node in t.traverse():
        san = node.name or ""
        orig = san_to_orig.get(san, san)
        if not orig:
            continue

        # 根节点
        if node == root_node:
            ml_result[root_name] = (root_age, root_age, root_age, "root_calibration")
            continue

        # 高可信锚点节点：直接用 ρ
        if orig in calib_map:
            e = calib_map[orig]
            ml_result[orig] = (
                e["age_years"], e["ci95_lower_years"], e["ci95_upper_years"],
                "rho_direct"
            )
            continue

        # 找最近校准祖先 A
        A_node = nearest_anc.get(id(node), root_node)
        orig_A = san_to_orig.get(A_node.name, A_node.name)
        if orig_A == root_san or A_node == root_node:
            orig_A = root_name
        eA = calib_map.get(orig_A, calib_map[root_name])
        age_A = eA["age_years"]
        lo_A  = eA["ci95_lower_years"]
        hi_A  = eA["ci95_upper_years"]
        D_root_to_A = dist_from_root[id(A_node)]
        D_root_to_V = dist_from_root[id(node)]
        D_AV = D_root_to_V - D_root_to_A

        # 收集子树中所有校准后代
        interp_ages, interp_lo, interp_hi = [], [], []
        for D_root_to_L, eD in calib_below.get(id(node), []):
            age_D = eD["age_years"]
            if age_D >= age_A:        # 后代必须比祖先年轻
                continue
            D_AD = D_root_to_L - D_root_to_A
            if D_AD <= D_AV:          # 后代必须在 V 之下
                continue
            av, lo, hi = _interpolate(
                D_AV, D_AD,
                age_A, age_D,
                lo_A, hi_A,
                eD["ci95_lower_years"], eD["ci95_upper_years"],
                root_age,
            )
            interp_ages.append(av)
            interp_lo.append(lo)
            interp_hi.append(hi)

        if interp_ages:
            age_V = float(np.median(interp_ages))
            lo_V  = float(np.median(interp_lo))
            hi_V  = float(np.median(interp_hi))
            ml_result[orig] = (age_V, lo_V, hi_V, "local_interpolation")
        elif orig in rho_idx:
            r = rho_idx[orig]
            ml_result[orig] = (
                float(r["age_years"]),
                max(0.0, float(r["ci95_lower_years"])),
                float(r["ci95_upper_years"]),
                "rho_fallback_no_desc",
            )
        else:
            ml_result[orig] = (age_A, lo_A, hi_A, "ancestor_only")

    log.info("ML 树节点年龄计算完成：%d 个", len(ml_result))

    # ── 处理 1,059 个单子节点（IQ-TREE 合并的）──────────────────────────────
    def _get_age(name: str):
        if name in ml_result:
            return ml_result[name][:3]
        if name in rho_idx:
            r = rho_idx[name]
            return (float(r["age_years"]),
                    max(0.0, float(r["ci95_lower_years"])),
                    float(r["ci95_upper_years"]))
        return None, None, None

    n_single_done = 0
    for haplo, (parent_name, child_name) in single_child.items():
        if haplo in ml_result:
            continue
        aP, lP, hP = _get_age(parent_name)
        aC, lC, hC = _get_age(child_name)
        if aP is not None and aC is not None:
            ml_result[haplo] = (
                (aP + aC) / 2.0,
                (lP + lC) / 2.0,
                (hP + hC) / 2.0,
                "single_child_interpolation",
            )
            n_single_done += 1
        elif haplo in rho_idx:
            r = rho_idx[haplo]
            ml_result[haplo] = (
                float(r["age_years"]),
                max(0.0, float(r["ci95_lower_years"])),
                float(r["ci95_upper_years"]),
                "rho_fallback_single_child",
            )
            n_single_done += 1
    log.info("单子节点处理：%d / %d", n_single_done, len(single_child))

    # ── 构建输出 DataFrame ────────────────────────────────────────────────────
    rows = []
    for haplo, (age, lo, hi, method) in ml_result.items():
        n_samples = int(rho_idx[haplo]["n_samples"]) if haplo in rho_idx else 0
        rows.append({
            "haplogroup":            haplo,
            "ml_age_years":          round(age, 1),
            "ml_ci95_lower_years":   round(lo, 1),
            "ml_ci95_upper_years":   round(hi, 1),
            "ml_age_kya":            round(age / 1000.0, 4),
            "ml_ci95_lower_kya":     round(lo / 1000.0, 4),
            "ml_ci95_upper_kya":     round(hi / 1000.0, 4),
            "n_samples":             n_samples,
            "calibration_method":    method,
        })
    df = pd.DataFrame(rows).sort_values("haplogroup").reset_index(drop=True)
    df.to_csv(out_dating, sep="\t", index=False)
    log.info("ML dating 写出：%s（%d 行）", out_dating, len(df))
    log.info("校准方法分布：\n%s", df["calibration_method"].value_counts().to_string())

    # ── 与 ρ 结果合并 ─────────────────────────────────────────────────────────
    rho_out = rho_complete[
        ["haplogroup", "n_samples", "age_kya", "ci95_lower_kya", "ci95_upper_kya", "confidence"]
    ].rename(columns={
        "age_kya":        "rho_age_kya",
        "ci95_lower_kya": "rho_ci95_lower_kya",
        "ci95_upper_kya": "rho_ci95_upper_kya",
    })
    merged = rho_out.merge(
        df[["haplogroup", "ml_age_kya", "ml_ci95_lower_kya", "ml_ci95_upper_kya", "calibration_method"]],
        on="haplogroup", how="outer"
    )
    merged["diff_kya"] = merged["ml_age_kya"] - merged["rho_age_kya"]
    denom = merged["rho_age_kya"].abs()
    merged["rel_diff_pct"] = np.where(denom > 0.001,
                                       merged["diff_kya"].abs() / denom * 100,
                                       float("nan"))
    merged.sort_values("haplogroup", inplace=True)
    merged.reset_index(drop=True, inplace=True)
    merged.to_csv(out_comparison, sep="\t", index=False)
    log.info("比较表写出：%s（%d 行）", out_comparison, len(merged))

    # ── 关键节点验证摘要 ─────────────────────────────────────────────────────
    log.info("── 关键节点验证 ──")
    for hap in ["mt-MRCA(RSRS)", "L3", "M", "N", "R", "H", "D4", "U"]:
        row = merged[merged["haplogroup"] == hap]
        if row.empty:
            log.info("  %-25s  未在比较表中", hap)
            continue
        r = row.iloc[0]
        log.info(
            "  %-25s  rho=%6.1f kya  ML=%6.1f kya  diff=%+.1f kya  n=%-7s  方法=%s",
            hap,
            r.get("rho_age_kya", float("nan")),
            r.get("ml_age_kya",  float("nan")),
            r.get("diff_kya",    float("nan")),
            int(r["n_samples"]) if pd.notna(r.get("n_samples")) else "NA",
            r.get("calibration_method", "?"),
        )


# ── CLI ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--treefile",          required=True)
    p.add_argument("--map",               required=True)
    p.add_argument("--rho",               required=True)
    p.add_argument("--hap-graph",         required=True)
    p.add_argument("--root-name",         default=ROOT_NAME_DEFAULT)
    p.add_argument("--root-age",          type=float, default=ROOT_AGE_DEFAULT)
    p.add_argument("--calib-min-n",       type=int,   default=CALIB_MIN_N_DEFAULT)
    p.add_argument("--outgroup",          default=OUTGROUP_DEFAULT)
    p.add_argument("--output-dating",     required=True)
    p.add_argument("--output-comparison", required=True)
    return p


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args()
    run(
        treefile=args.treefile,
        map_tsv=args.map,
        rho_tsv=args.rho,
        hap_graph_json=args.hap_graph,
        root_name=args.root_name,
        root_age=args.root_age,
        calib_min_n=args.calib_min_n,
        outgroup=args.outgroup,
        out_dating=args.output_dating,
        out_comparison=args.output_comparison,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2-cluster.py — 逐时间窗聚类、跨时间对齐、后验重抽同簇频率与降维坐标。
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ne_core import (FEATURE_DISPLAY, FEATURE_KEYS,  # noqa: E402
                     align_cluster_labels_across_bins, cluster_features,
                     feature_distance, scale_features, segment_features, window_mask)
from pipeline_utils import load_config, resolve_paths, setup_logging, write_tsv  # noqa: E402

log = logging.getLogger(__name__)

HEADER = ["Quantity", "Window Type", "Time Window", "Window Start (Years)",
          "Category", "Second Category", "Group", "Cluster", "Component",
          "Score Component", "Value", "Reliable Window Coverage", "Usable"]


def build_chord_edges(co: np.ndarray, names: list[str], threshold: float
                      ) -> list[tuple[str, str, float]]:
    """同簇频率矩阵 → 弦图边表，只保留不低于阈值的类别对。

    阈值以下的连边在图上只会变成噪声：几十条几乎透明的弦盖住真正稳定的关系。
    """
    edges = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            v = float(co[i, j])
            if np.isfinite(v) and v >= threshold:
                edges.append((names[i], names[j], v))
    return edges


# 两个分数各自的特征来源与定向特征（定向特征的载荷强制为正，否则同一批数据在不同次
# 运行里可能整体反号，颜色含义会翻转）
SCORE_COMPONENTS = {
    "Level Component": (("level", "amplitude"), "level"),
    "Change Component": (("delta", "rate", "max_expansion_rate",
                          "max_contraction_rate"), "rate"),
}


def composite_score(cells: np.ndarray, keys: tuple[str, ...], anchor: str
                    ) -> tuple[np.ndarray, np.ndarray, float]:
    """把一组特征压成一个分数，返回 (逐格分数, 权重, 解释率)。

    cells (n_cells, 8) 是全部"类别 × 时间 bin"的特征，先按 keys 取子集再逐特征标准化
    ——特征量纲不同，不标准化时数值量级大的那个会独占分数——最后取第一主成分得分。
    """
    from sklearn.decomposition import PCA

    idx = [FEATURE_KEYS.index(k) for k in keys]
    z = scale_features(cells[:, idx], "zscore")
    pca = PCA(n_components=1)
    score = pca.fit_transform(z)[:, 0]
    load = pca.components_[0]
    if load[keys.index(anchor)] < 0:
        score, load = -score, -load
    return score, load, float(pca.explained_variance_ratio_[0])


def standardize_per_period(scores: np.ndarray, windows: list[str]) -> np.ndarray:
    """逐时间段做 z-score：同一个时间 bin 内的类别之间比较。

    跨时间的整体漂移会淹没同一时段的类别差异，热图上就成了一片同色。逐段标准化后颜色
    只在列内可比，读法是"这个时段谁偏高、谁偏低"，不能跨列读绝对高低。
    """
    out = np.full(len(scores), np.nan)
    win = np.array(windows)
    for w in np.unique(win):
        m = win == w
        x = scores[m]
        sd = x.std()
        out[m] = (x - x.mean()) / (sd if sd > 1e-12 else 1.0)
    return out


def scatter_coordinates(mat: np.ndarray, method: str = "pca"
                        ) -> tuple[np.ndarray, np.ndarray]:
    """特征矩阵 → 二维坐标与解释率（UMAP 无解释率，返回 NaN）。"""
    if method == "umap":
        try:
            from umap import UMAP
        except ImportError as exc:                      # noqa: BLE001
            raise ImportError("scatter_method = umap 需要 umap-learn："
                              "conda install -c conda-forge umap-learn") from exc
        n_neighbors = max(2, min(15, mat.shape[0] - 1))
        emb = UMAP(n_components=2, n_neighbors=n_neighbors, random_state=0).fit_transform(mat)
        return np.asarray(emb), np.array([np.nan, np.nan])
    from sklearn.decomposition import PCA
    pca = PCA(n_components=min(2, mat.shape[0] - 1, mat.shape[1]))
    emb = pca.fit_transform(mat)
    if emb.shape[1] == 1:                               # 类别极少时补一列 0，图仍可画
        emb = np.column_stack([emb, np.zeros(len(emb))])
    ratio = np.full(2, np.nan)
    ratio[:len(pca.explained_variance_ratio_)] = pca.explained_variance_ratio_
    return emb, ratio


def coclustering_frequency(curves: dict[str, np.ndarray], rates: dict[str, np.ndarray],
                           grid: np.ndarray, mask: np.ndarray, names: list[str],
                           k: int, metric: str, method: str, scaling: str,
                           n_resample: int, rng: np.random.Generator) -> np.ndarray:
    """后验重抽同簇频率：每次每个类别各抽 1 条曲线，重算特征后重新聚类。"""
    co = np.zeros((len(names), len(names)))
    for _ in range(n_resample):
        feats = []
        for c in names:
            pick = rng.integers(curves[c].shape[0])
            f = segment_features(curves[c][pick:pick + 1], rates[c][pick:pick + 1],
                                 grid, mask)
            feats.append([float(f[key][0]) for key in FEATURE_KEYS])
        lab, _, _ = cluster_features(scale_features(np.array(feats), scaling),
                                     metric, method, n_clusters=k)
        co += (lab[:, None] == lab[None, :]).astype(float)
    return co / max(n_resample, 1)


def run(config: str, project_root: str = ".", log_file: str | None = None) -> int:
    cfg = resolve_paths(load_config(config), project_root)
    setup_logging(cfg["runtime"]["log_level"], log_file)
    io, po, se, cl, fg = (cfg["io"], cfg["posterior"], cfg["select"], cfg["cluster"],
                          cfg["figure"])
    out_dir = Path(io["output_dir"])
    rng = np.random.default_rng(int(po["seed"]))

    with np.load(out_dir / ".cache" / "curves.npz", allow_pickle=True) as z:
        grid = z["grid"]
        cats = [str(x) for x in z["categories"]]
        groups = {c: str(g) for c, g in zip(cats, z["group"])}
        has_post = bool(z["has_posterior"])
        wmeta = z["window_meta"]
        fmat = z["feature_matrix"]                     # (n_windows, n_cats, 8)
        coverage = z["coverage"]
        curves = {c: z[f"curve__{c}"].astype(float) for c in cats}
        rates = {c: z[f"rate__{c}"].astype(float) for c in cats}
    if not has_post:
        log.warning("输入不含真后验，同簇频率由伪样本得到，只能作相似性强度读，"
                    "不能称作后验概率")

    rows: list[list] = []
    prev_by_type: dict[str, dict[str, int]] = {}
    debug_rows: list[list] = []
    score_cells: list[np.ndarray] = []
    score_index: list[tuple[str, str, float, str, float]] = []
    for wi, (wtype, wname, t0, t1) in enumerate(wmeta):
        t0, t1 = float(t0), float(t1)
        cov = coverage[wi]
        usable = [i for i, c in enumerate(cats) if cov[i] >= float(se["min_coverage"])]
        # 合并分数按可用的格子逐个收集，不可靠的格子留空（热图上显示为灰）
        if wtype == "Time Bin":
            for i in usable:
                score_cells.append(fmat[wi][i])
                score_index.append((wtype, wname, t0, cats[i], float(cov[i])))
        if len(usable) < 3:
            log.warning("%s 只有 %d 个类别达到覆盖度 %.2f，该窗口不聚类",
                        wname, len(usable), float(se["min_coverage"]))
            for i, c in enumerate(cats):
                rows.append(["Cluster Membership", wtype, wname, f"{t0:g}", c, "",
                             groups[c], "Not Assessed", "", "", f"{cov[i]:.3f}", "No"])
            continue
        names = [cats[i] for i in usable]
        scaled = scale_features(fmat[wi][usable], cl["feature_scaling"])
        lab, k_used, scores = cluster_features(
            scaled, cl["distance_metric"], cl["linkage_method"],
            int(cl["n_clusters"]), int(cl["k_max"]))
        cur = {n: int(v) for n, v in zip(names, lab)}
        if prev_by_type.get(wtype):
            cur = align_cluster_labels_across_bins(prev_by_type[wtype], cur)
        prev_by_type[wtype] = cur
        for i, c in enumerate(cats):
            rows.append(["Cluster Membership", wtype, wname, f"{t0:g}", c, "", groups[c],
                         f"Cluster {cur[c]}" if c in cur else "Not Assessed", "", "", "",
                         f"{cov[i]:.3f}", "Yes" if c in cur else "No"])
        if cfg["runtime"]["debug_mode"]:
            for kk, sc in scores.items():
                debug_rows.append(["Silhouette Coefficient", wtype, wname, kk, "", f"{sc:.4f}"])
            d = feature_distance(scaled, cl["distance_metric"])
            for a in range(len(names)):
                for b in range(a + 1, len(names)):
                    debug_rows.append(["Feature Distance", wtype, wname, "", names[a],
                                       f"{names[b]}\t{d[a, b]:.5f}"])

        # 时间段另外给弦图与散点图所需的量：同簇频率与降维坐标
        if wtype != "Time Period":
            continue
        mask = window_mask(grid, t0, t1)
        # 两个配对量都写进长表：弦图按 chord_metric 二选一，配对矩阵热图要把两者并排对照
        pair_mats = {}
        pair_mats["Co-clustering Frequency"] = coclustering_frequency(
            curves, rates, grid, mask, names, k_used, cl["distance_metric"],
            cl["linkage_method"], cl["feature_scaling"], int(po["n_resample"]), rng)
        d = feature_distance(scaled, cl["distance_metric"])
        pair_mats["Similarity"] = 1.0 - d / (d.max() or 1.0)
        co = pair_mats["Co-clustering Frequency" if fg["chord_metric"] == "coclustering"
                       else "Similarity"]
        for quantity, mat in pair_mats.items():
            for a in range(len(names)):
                for b in range(a + 1, len(names)):
                    rows.append([quantity, wtype, wname, f"{t0:g}", names[a], names[b],
                                 groups[names[a]], "", "", "", f"{mat[a, b]:.4f}", "",
                                 "Yes"])
        emb, ratio = scatter_coordinates(scaled, fg["scatter_method"])
        for a, c in enumerate(names):
            for comp in (0, 1):
                rows.append(["Scatter Coordinate", wtype, wname, f"{t0:g}", c, "",
                             groups[c], f"Cluster {cur[c]}", f"Component {comp + 1}", "",
                             f"{emb[a, comp]:.5f}", f"{cov[usable[a]]:.3f}", "Yes"])
        for comp in (0, 1):
            rows.append(["Explained Variance", wtype, wname, f"{t0:g}", "", "", "", "",
                         f"Component {comp + 1}", "",
                         "" if not np.isfinite(ratio[comp]) else f"{ratio[comp]:.5f}",
                         "", ""])
        # 节点大小 = 该时间段的平均 log Ne，弦图用它表示类别的规模
        for a, c in enumerate(names):
            rows.append(["Node Size", wtype, wname, f"{t0:g}", c, "", groups[c],
                         f"Cluster {cur[c]}", "", "",
                         f"{fmat[wi][usable[a]][FEATURE_KEYS.index('level')]:.5f}",
                         "", "Yes"])
        edges = build_chord_edges(co, names, float(fg["chord_threshold"]))
        log.info("%s：%d 个可用类别，K=%d，弦图保留 %d 条连边（阈值 %.2f）",
                 wname, len(names), k_used, len(edges), float(fg["chord_threshold"]))

    # 全时间轴热图的颜色：水平与变化两个分数，各自逐时间段标准化
    if score_cells:
        cells = np.array(score_cells)
        win_names = [x[1] for x in score_index]
        for comp_name, (keys, anchor) in SCORE_COMPONENTS.items():
            raw, weights, ratio = composite_score(cells, keys, anchor)
            std = standardize_per_period(raw, win_names)
            for (wtype, wname, t0, cat, cov_v), v, z in zip(score_index, raw, std):
                rows.append(["Composite Trajectory Score", wtype, wname, f"{t0:g}", cat,
                             "", groups[cat], "", "", comp_name, f"{v:.5f}",
                             f"{cov_v:.3f}", "Yes"])
                rows.append(["Standardized Composite Score", wtype, wname, f"{t0:g}", cat,
                             "", groups[cat], "", "", comp_name, f"{z:.5f}",
                             f"{cov_v:.3f}", "Yes"])
            for key, w in zip(keys, weights):
                rows.append(["Score Weight", "", "", "", "", "", "", "",
                             FEATURE_DISPLAY[key], comp_name, f"{w:.5f}", "", ""])
            rows.append(["Score Explained Variance", "", "", "", "", "", "", "", "",
                         comp_name, f"{ratio:.5f}", "", ""])
            log.info("%s：%d 个可用格子，第一主成分解释 %.1f%%（已逐时间段标准化）",
                     comp_name, len(score_cells), 100 * ratio)

    n = write_tsv(str(out_dir / "⭐2-Cluster-And-Similarity.tsv"), HEADER, rows)
    log.info("聚类与相似性长表 -> %s（%d 行）",
             out_dir / "⭐2-Cluster-And-Similarity.tsv", n)
    if cfg["runtime"]["debug_mode"] and debug_rows:
        write_tsv(str(out_dir / "debug" / "cluster-diagnostics.tsv"),
                  ["Quantity", "Window Type", "Time Window", "Number of Clusters",
                   "Category", "Value"], debug_rows)
        log.info("诊断表 -> %s", out_dir / "debug" / "cluster-diagnostics.tsv")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="逐时间窗聚类、跨时间对齐与同簇频率")
    p.add_argument("--config", required=True)
    p.add_argument("--project-root", default=".")
    p.add_argument("--log-file", default=None)
    return p


def main(argv=None) -> int:
    return run(**vars(build_parser().parse_args(argv)))


if __name__ == "__main__":
    raise SystemExit(main())

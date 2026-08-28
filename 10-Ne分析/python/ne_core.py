#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ne_core.py — Ne 轨迹分析的通用算法库（不含 main，只被 import）。
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


# ============================================================ 1. Skygrid 展开


def skygrid_boundaries(n_popsize: int, cutoff: float) -> np.ndarray:
    m = n_popsize - 1
    if m < 1:
        raise ValueError("logPopSize 至少要有 2 个维度才谈得上分段")
    return np.array([cutoff * i / m for i in range(1, m + 1)], dtype=float)


def skygrid_log_popsize(values: np.ndarray, cutoff: float,
                        grid: np.ndarray) -> np.ndarray:
    k = values.shape[1]
    bounds = np.append(skygrid_boundaries(k, cutoff), np.inf)
    idx = np.searchsorted(bounds, grid, side="right")
    return values[:, np.clip(idx, 0, k - 1)]


_XML_DIM = re.compile(r'id="skygrid\.logPopSize"\s+dimension="(\d+)"')
_XML_NGP = re.compile(r'id="skygrid\.numGridPoints"\s+value="([\d.eE+-]+)"')
_XML_CUT = re.compile(r'id="skygrid\.cutOff"\s+value="([\d.eE+-]+)"')


def read_xml_skygrid(xml_path: str | Path) -> dict[str, float] | None:
    p = Path(xml_path)
    if not p.is_file():
        return None
    text = p.read_text(encoding="utf-8", errors="ignore")
    dim, ngp, cut = _XML_DIM.search(text), _XML_NGP.search(text), _XML_CUT.search(text)
    if not cut:
        return None
    out = {"cutoff": float(cut.group(1))}
    out["dimension"] = float(dim.group(1)) if dim else float("nan")
    out["n_grid_points"] = float(ngp.group(1)) if ngp else float("nan")
    if dim and ngp and int(out["dimension"]) != int(out["n_grid_points"]) + 1:
        log.warning("%s 的 XML 自身不自洽：dimension=%d, numGridPoints=%d",
                    p.name, int(out["dimension"]), int(out["n_grid_points"]))
    return out


def check_grid_against_xml(name: str, n_popsize: int, xml_info: dict | None,
                           analysis_max: float) -> float:
    if xml_info is None:
        raise FileNotFoundError(
            f"{name}: 找不到可解析的 XML（需要 skygrid.cutOff）。请把与 .log 同名的 "
            f".xml 放在一起，或在配置里显式给出 time.cutoff_override")
    dim = xml_info.get("dimension", float("nan"))
    if np.isfinite(dim) and int(dim) != n_popsize:
        raise ValueError(f"{name}: XML dimension={int(dim)} 与日志的 logPopSize 列数 "
                         f"{n_popsize} 不一致")
    cutoff = float(xml_info["cutoff"])
    if analysis_max > cutoff + 1e-9:
        raise ValueError(f"{name}: 分析上界 {analysis_max:g} 超出 Skygrid cutOff "
                         f"{cutoff:g}；cutOff 之外是延续最后一段的常数，不是人口历史")
    return cutoff


def thin_indices(n: int, keep: int) -> np.ndarray:
    if keep <= 0 or n <= keep:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, keep).astype(int))


def display_name(raw: str) -> str:
    txt = re.sub(r"[_\-]+", " ", str(raw)).strip()
    return " ".join(w if w.isupper() else w.capitalize() for w in txt.split()) or str(raw)


# ============================================================ 2. 平滑与求导


def loess_operators(x: np.ndarray, bandwidth: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(x)
    span_raw = float(x[-1] - x[0])
    span = bandwidth / span_raw if span_raw > 0 else 1.0
    span = min(1.0, max(span, 8.0 / n))
    n_neighbors = max(int(round(span * n)), 3)

    h0 = np.zeros((n, n))
    h1 = np.zeros((n, n))
    for j in range(n):
        dist = np.abs(x - x[j])
        idx = np.argpartition(dist, n_neighbors - 1)[:n_neighbors]
        d = dist[idx]
        dmax = d.max() or 1e-12
        w = np.clip(1.0 - (d / dmax) ** 3, 0.0, None) ** 3
        dx = x[idx] - x[j]
        design = np.column_stack([np.ones_like(dx), dx, dx ** 2])
        wd = design * w[:, None]
        try:
            coef = np.linalg.solve(design.T @ wd, wd.T)
        except np.linalg.LinAlgError:
            coef = np.linalg.pinv(design.T @ wd) @ wd.T
        h0[j, idx] = coef[0, :]
        h1[j, idx] = coef[1, :]
    return h0, h1


# ============================================================ 3. 后验汇总


def hpd_interval(x: np.ndarray, cred: float = 0.95) -> tuple[float, float]:
    """最高后验密度区间（单峰假设下的最短覆盖区间）。"""
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0:
        return (np.nan, np.nan)
    m = max(int(np.floor(cred * n)), 1)
    if m >= n:
        return (float(x[0]), float(x[-1]))
    widths = x[m:] - x[:n - m]
    i = int(np.argmin(widths))
    return (float(x[i]), float(x[i + m]))


def hpd_width_curve(draws: np.ndarray, cred: float = 0.95) -> np.ndarray:
    """逐时间点的 HPD 宽度，用于判定曲线在哪一段还"说得清"。"""
    return np.array([np.subtract(*reversed(hpd_interval(draws[:, k], cred)))
                     for k in range(draws.shape[1])])


def reliable_mask(draws: np.ndarray, max_width: float,
                  cred: float = 0.95) -> np.ndarray:
    if draws.shape[0] < 2:
        return np.ones(draws.shape[1], dtype=bool)
    ok = hpd_width_curve(draws, cred) <= max_width
    out = np.zeros_like(ok)
    for i in range(len(ok)):
        if not ok[i]:
            break
        out[i] = True
    return out


# ============================================================ 4. 分段特征


def window_mask(grid: np.ndarray, t_start: float, t_end: float) -> np.ndarray:
    """闭区间 [t_start, t_end] 在网格上的布尔掩码。"""
    return (grid >= t_start - 1e-9) & (grid <= t_end + 1e-9)


def segment_features(y: np.ndarray, r: np.ndarray, grid: np.ndarray,
                     mask: np.ndarray) -> dict[str, np.ndarray]:
    g = grid[mask]
    yy = y[:, mask]
    rr = r[:, mask]
    span = float(g[-1] - g[0])

    level = np.trapezoid(yy, g, axis=1) / span if span > 0 else yy[:, 0]
    delta = yy[:, 0] - yy[:, -1]                        # 由老到新的净变化
    rate = delta / span if span > 0 else np.zeros(len(yy))
    amplitude = yy.max(axis=1) - yy.min(axis=1)

    i_exp = np.argmax(rr, axis=1)
    i_con = np.argmin(rr, axis=1)
    rows = np.arange(len(rr))
    return {
        "level": level,
        "delta": delta,
        "rate": rate,
        "amplitude": amplitude,
        "max_expansion_time": g[i_exp],
        "max_expansion_rate": rr[rows, i_exp],
        "max_contraction_time": g[i_con],
        "max_contraction_rate": rr[rows, i_con],
    }


FEATURE_KEYS = ("level", "delta", "rate", "amplitude", "max_expansion_time",
                "max_expansion_rate", "max_contraction_time", "max_contraction_rate")

FEATURE_DISPLAY = {
    "level": "Mean Log Effective Population Size",
    "delta": "Net Change in Log Effective Population Size",
    "rate": "Mean Rate of Change",
    "amplitude": "Fluctuation Amplitude",
    "max_expansion_time": "Time of Maximum Expansion",
    "max_expansion_rate": "Maximum Expansion Rate",
    "max_contraction_time": "Time of Maximum Contraction",
    "max_contraction_rate": "Maximum Contraction Rate",
}

TIME_FEATURES = frozenset({"max_expansion_time", "max_contraction_time"})

# 方向性后验概率怎么定义：(比较方向, 用来判方向的特征, 自然语言说明)。
# 事件时间本身没有正负，方向支持度取自对应的极值速率；水平与幅度恒为正，不给概率。
PROBABILITY_SPEC: dict[str, tuple[str | None, str, str]] = {
    "level": (None, "", ""),
    "delta": ("greater", "delta", "Net increase over the window"),
    "rate": ("greater", "rate", "Positive mean rate of change"),
    "amplitude": (None, "", ""),
    "max_expansion_time": ("greater", "max_expansion_rate",
                           "An expansion occurs within the window"),
    "max_expansion_rate": ("greater", "max_expansion_rate",
                           "An expansion occurs within the window"),
    "max_contraction_time": ("less", "max_contraction_rate",
                             "A contraction occurs within the window"),
    "max_contraction_rate": ("less", "max_contraction_rate",
                             "A contraction occurs within the window"),
}


def directional_probability(key: str, feats: dict[str, np.ndarray]
                            ) -> tuple[float, str]:
    direction, source, statement = PROBABILITY_SPEC.get(key, (None, "", ""))
    if direction is None:
        return float("nan"), ""
    x = np.asarray(feats[source], dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), statement
    prob = float(np.mean(x > 0)) if direction == "greater" else float(np.mean(x < 0))
    return prob, statement


def at_window_boundary(key: str, median_time: float, grid_window: np.ndarray) -> str:
    if key not in TIME_FEATURES or grid_window.size == 0:
        return ""
    edge = min(abs(median_time - float(grid_window[0])),
               abs(median_time - float(grid_window[-1])))
    return "Yes" if edge < 1e-6 else "No"

# ============================================================ 5. 规模校正


def size_adjustment_fit(level: np.ndarray, log_size: np.ndarray) -> tuple[float, float]:
    x = np.asarray(log_size, dtype=float)
    y = np.asarray(level, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.allclose(x[ok], x[ok][0]):
        return 0.0, 0.0
    slope, intercept = np.polyfit(x[ok], y[ok], 1)
    return float(slope), float(intercept)


def apply_size_adjustment(curves: np.ndarray, log_size: np.ndarray,
                          slope: float, intercept: float) -> np.ndarray:
    shift = slope * np.asarray(log_size, dtype=float) + intercept
    return curves - shift[:, None]


# ============================================================ 6. 距离与聚类


def scale_features(mat: np.ndarray, method: str = "zscore") -> np.ndarray:
    x = np.asarray(mat, dtype=float)
    if method in ("none", "", None):
        return x
    if method == "zscore":
        sd = x.std(axis=0, ddof=0)
        sd[sd < 1e-12] = 1.0
        return (x - x.mean(axis=0)) / sd
    if method == "minmax":
        lo, hi = x.min(axis=0), x.max(axis=0)
        rng = np.where(hi - lo < 1e-12, 1.0, hi - lo)
        return (x - lo) / rng
    raise ValueError(f"未知的 feature_scaling：{method}")


# 只有这些度量是（加权）欧氏距离，Ward 的方差递推公式仅在欧氏距离下成立
EUCLIDEAN_METRICS = frozenset({"euclidean", "seuclidean"})


def linkage_is_valid(metric: str, method: str) -> bool:
    """Ward 只能配欧氏距离；其余连接方式对任意距离矩阵都有定义。"""
    return method != "ward" or metric in EUCLIDEAN_METRICS


def feature_distance(mat: np.ndarray, metric: str = "euclidean") -> np.ndarray:
    """特征矩阵 → 方阵距离矩阵，直接调 scipy 的 pdist。"""
    from scipy.spatial.distance import pdist, squareform
    return squareform(pdist(np.asarray(mat, dtype=float), metric=metric))


def condensed(mat: np.ndarray) -> np.ndarray:
    m = np.array(mat, dtype=float, copy=True)
    finite = m[np.isfinite(m)]
    fill = float(finite.max()) if finite.size else 0.0
    m[~np.isfinite(m)] = fill
    np.fill_diagonal(m, 0.0)
    m = (m + m.T) / 2.0
    iu = np.triu_indices_from(m, k=1)
    return m[iu]


def cluster_features(mat: np.ndarray, metric: str = "euclidean",
                     method: str = "average", n_clusters: int = 0,
                     k_max: int = 8) -> tuple[np.ndarray, int, dict[int, float]]:
    from scipy.cluster.hierarchy import fcluster, linkage
    from sklearn.metrics import silhouette_score

    if not linkage_is_valid(metric, method):
        raise ValueError(f"连接方式 {method} 要求欧氏距离，与 {metric} 不兼容")
    n = mat.shape[0]
    dist = feature_distance(mat, metric)
    z = linkage(condensed(dist), method=method)
    upper = min(k_max, n - 1)
    scores: dict[int, float] = {}
    for k in range(2, upper + 1):
        lab = fcluster(z, k, criterion="maxclust")
        if len(set(lab)) > 1:
            scores[k] = float(silhouette_score(dist, lab, metric="precomputed"))
    k_use = int(n_clusters) if n_clusters and n_clusters >= 2 else (
        max(scores, key=scores.get) if scores else 2)
    k_use = max(2, min(k_use, upper)) if upper >= 2 else 1
    return fcluster(z, k_use, criterion="maxclust"), k_use, scores


def align_cluster_labels_across_bins(prev: dict[str, int], cur: dict[str, int]
                                     ) -> dict[str, int]:
    from scipy.optimize import linear_sum_assignment

    shared = sorted(set(prev) & set(cur))
    pl = sorted({prev[c] for c in shared})
    cl = sorted(set(cur.values()))
    mapping: dict[int, int] = {}
    if shared and pl:
        gain = np.zeros((len(cl), len(pl)))
        for a, c in enumerate(cl):
            for b, q in enumerate(pl):
                gain[a, b] = sum(1 for s in shared if cur[s] == c and prev[s] == q)
        rows, cols = linear_sum_assignment(-gain)
        mapping = {cl[a]: pl[b] for a, b in zip(rows, cols) if gain[a, b] > 0}
    free = [x for x in range(1, len(cl) + len(pl) + 1) if x not in mapping.values()]
    for c in cl:
        if c not in mapping:
            mapping[c] = free.pop(0)
    return {s: mapping[c] for s, c in cur.items()}

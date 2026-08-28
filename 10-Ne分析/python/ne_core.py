#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ne_core.py — Ne 轨迹分析的通用算法库（不含 main，只被 import）。

本文件只放**与具体课题无关**的算法：Skygrid 后验展开、loess 平滑与求导、分段特征、
事件方向判据、可靠时间窗、规模校正、距离与聚类的适用性检查。任何分组名、时间范围、
文件命名都不出现在这里，一律由调用方从配置传入。

聚类、轮廓系数、主成分、匈牙利匹配一律调用 scipy / scikit-learn 的成熟实现。
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


# ============================================================ 1. Skygrid 展开


def skygrid_boundaries(n_popsize: int, cutoff: float) -> np.ndarray:
    """Skygrid 的内部分段边界（不含 +inf 尾段），单位与 cutOff 一致。

    BEAST 的 gmrfSkyGridLikelihood 里 logPopSize 的 dimension = numGridPoints + 1，
    边界落在 cutOff * i / numGridPoints（i = 1..numGridPoints）。
    """
    m = n_popsize - 1
    if m < 1:
        raise ValueError("logPopSize 至少要有 2 个维度才谈得上分段")
    return np.array([cutoff * i / m for i in range(1, m + 1)], dtype=float)


def skygrid_log_popsize(values: np.ndarray, cutoff: float,
                        grid: np.ndarray) -> np.ndarray:
    """把 logPopSize 后验矩阵展开到任意时间网格上。

    Skygrid 是分段常数：Ne(t) = theta_k 当 x_{k-1} <= t < x_k。时间恰好等于内部边界
    时应归入**后**一段，因此 searchsorted 必须用 side="right"——side="left" 会把边界
    点错误地留在前一段，而格宽与分析网格步长常有整除关系，边界重合并非边角情况。

    values (n_draws, n_popsize)；返回 (n_draws, n_grid)。
    """
    k = values.shape[1]
    bounds = np.append(skygrid_boundaries(k, cutoff), np.inf)
    idx = np.searchsorted(bounds, grid, side="right")
    return values[:, np.clip(idx, 0, k - 1)]


_XML_DIM = re.compile(r'id="skygrid\.logPopSize"\s+dimension="(\d+)"')
_XML_NGP = re.compile(r'id="skygrid\.numGridPoints"\s+value="([\d.eE+-]+)"')
_XML_CUT = re.compile(r'id="skygrid\.cutOff"\s+value="([\d.eE+-]+)"')


def read_xml_skygrid(xml_path: str | Path) -> dict[str, float] | None:
    """从 BEAST XML 读出 Skygrid 的时间轴定义，读不到返回 None。

    XML 是 cutOff 与网格数的**权威来源**：.log 文件里只有 logPopSize 的取值，没有任何
    时间信息，没有 cutOff 就无从确定每一段对应哪个时间区间。

    返回 {"dimension": …, "n_grid_points": …, "cutoff": …}。dimension 与
    n_grid_points 自身不自洽（应相差 1）时只记警告，由调用方结合 .log 的列数裁决。
    """
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
    """核对 .log 的 logPopSize 列数与 XML 的维度，返回该文件应使用的 cutOff。

    两处致命错误，直接抛异常而不是降级：
      1. XML 的 dimension 与 .log 的列数对不上——两者不是同一次分析；
      2. 分析上界超出 cutOff——cutOff 之后 Skygrid 不再有新参数，只是延续最后一段的
         常数值，把那一段当"人口长期平稳"解读是纯粹的假象。
    """
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
    """等间隔抽稀到 keep 条；n <= keep 或 keep <= 0 时原样返回。

    等间隔而非随机：MCMC 样本自相关，等间隔正是标准 thinning。
    """
    if keep <= 0 or n <= keep:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, keep).astype(int))


def display_name(raw: str) -> str:
    """把机器命名归一成可以直接进图表的展示名：下划线/连字符换空格并首字母大写。

    这只是**兜底**： 里给出的 display_name 永远优先。机械去下划线得到的
    名字往往仍不像自然语言，凡是要放进论文或汇报的类别名都应在清单里显式写好。
    """
    txt = re.sub(r"[_\-]+", " ", str(raw)).strip()
    return " ".join(w if w.isupper() else w.capitalize() for w in txt.split()) or str(raw)


# ============================================================ 2. 平滑与求导


def loess_operators(x: np.ndarray, bandwidth: float) -> tuple[np.ndarray, np.ndarray]:
    """构造 loess（局部二次、tricube 核）的"平滑"与"求导"两个线性算子。

    权重只由 x 决定、与 y 无关，因此算子只需算一次，对任意多条后验曲线都只是一次矩阵
    乘法。平滑值取局部二次拟合的常数项、导数取一次项，两者来自同一次拟合——否则峰位置
    与导数过零点会对不上。

    返回 (H0, H1)：smoothed = y @ H0.T，derivative = y @ H1.T。
    """
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
    """可靠时间窗：HPD 宽度不超过阈值、且**含最近端**的最长连续区间。

    不取全局最长区间：曲线可信度从现在往过去单调衰减，中间冒出的孤立合格区间通常是
    HPD 宽度的偶然回落，不构成可以引用的时段。单条曲线（无后验）时全部视为可靠。
    """
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
    """在给定时间窗内提取 8 个分段特征，每个特征返回长度 n_draws 的向量。

    y (n_draws, n_grid) 对数尺度 Ne；r 同形状的变化速率（现在方向为正）；网格按年代
    递增（0 = 现在），因此窗内"起点"是较老的一端、"终点"是较新的一端，净变化与速率
    都按"由老到新"定义。

    事件时间的陷阱：argmax / argmin 在任何窗口里都必然返回一个位置，但"窗内速率最大的
    时刻"不等于"发生了扩张"——若整窗速率都是负的，argmax 找到的只是**收缩最弱**的时刻。
    方向是否成立由 max_expansion_rate > 0 / max_contraction_rate < 0 逐条后验样本判定，
    汇总成方向支持度（见 directional_probability）。
    """
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
    """特征的方向性后验概率，返回 (概率, 自然语言说明)；无定义时概率为 NaN。

    收缩类取 P(速率 < 0) 而不是 P(> 0)：后者说的是"连最低速率都还是正的"，与"这里发生
    了收缩"恰好相反，拿来当收缩支持度会把结论读反。
    """
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
    """极值时间是否贴在窗口最边缘的一格上。

    贴边意味着真正的极值多半落在窗外、只是被窗口截断，不能当作"事件发生在这个年代"。
    """
    if key not in TIME_FEATURES or grid_window.size == 0:
        return ""
    edge = min(abs(median_time - float(grid_window[0])),
               abs(median_time - float(grid_window[-1])))
    return "Yes" if edge < 1e-6 else "No"


# ============================================================ 5. 规模校正


def size_adjustment_fit(level: np.ndarray, log_size: np.ndarray) -> tuple[float, float]:
    """拟合"窗内平均 log Ne ~ log 样本量"的一元回归，返回 (斜率, 截距)。

    动机：合并率 Ne 对群体结构敏感，越宽泛、内部越异质的类别估出的 Ne 越高（Wahlund
    类效应），而样本量往往正是"覆盖广度"的代理变量，不校正时聚类可能读到的是规模梯度
    而不是人口史模式。

    代价必须写明：回归**不能**证明由 log n 解释的那部分一定是估计偏差而不是真实的 Ne
    差异，两者在同一批数据里不可识别。因此本项默认关闭，开启时应与未校正结果并列报告。

    样本量少于 3 个或全部相同时无法识别，返回 (0, 0) 即不做校正。
    """
    x = np.asarray(log_size, dtype=float)
    y = np.asarray(level, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.allclose(x[ok], x[ok][0]):
        return 0.0, 0.0
    slope, intercept = np.polyfit(x[ok], y[ok], 1)
    return float(slope), float(intercept)


def apply_size_adjustment(curves: np.ndarray, log_size: np.ndarray,
                          slope: float, intercept: float) -> np.ndarray:
    """逐条曲线减去 slope * log n + intercept。

    每条曲线只被平移一个常数，形状、幅度、速率完全不受影响，只有整体水平被重设。
    参数在**全分析窗上只拟合一次**，所有时间段共用——逐窗重拟合会让校正函数本身随时间
    漂移，窗间差异就无法归因于人口史了。
    """
    shift = slope * np.asarray(log_size, dtype=float) + intercept
    return curves - shift[:, None]


# ============================================================ 6. 距离与聚类


def scale_features(mat: np.ndarray, method: str = "zscore") -> np.ndarray:
    """特征标准化。8 个特征量纲不同（对数尺度、每千年速率、年代），不标准化时
    年代类特征会凭数值量级独占距离。"""
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
    """对称距离方阵 → scipy 层次聚类要求的 condensed 向量；缺失记为最不相似。"""
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
    """特征矩阵 → 簇标签。n_clusters <= 0 时按轮廓系数在 [2, k_max] 内自选。

    返回 (标签, 实际簇数, 逐 K 的轮廓系数)。
    """
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
    """把当前时间 bin 的簇编号重排成与前一个 bin 最匹配的一套编号。

    层次聚类给出的簇编号是任意的，不对齐时"同一条带"会在图上反复换色。匹配只在**两个
    bin 共有的类别**上做（某类别在前一个 bin 不可用时不该参与匹配），以共有成员数为收益
    做一次匈牙利算法最大匹配；没匹配上的簇顺延到未占用的编号。
    """
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

#!/usr/bin/env python3
"""3-1 山脊密度：把节点年龄转成组间可比的脊线高度。

平滑尺度的权衡
  节点年龄极度右偏，线性尺度上的单一带宽要么把近期结构抹平，要么让深层尾部碎成锯齿；纯对数
  尺度上的固定带宽又会在年龄趋近 0 时把分辨率放到几十年，把密集采样的家族内合并放大成一根针。
  这里用“相对带宽 + 绝对下限”的自适应高斯核：

      h(a) = max(alpha * a, floor)

  alpha 由该维度全部类别合并后的对数年龄按 Silverman 规则求得（乘 ln10 换算成线性尺度上的
  相对宽度），floor 封住分辨率下限。

可分辨下界
  当 alpha * a < floor 时带宽被下限钉死，那一段曲线的平滑尺度不是数据设定的，而是下限本身；
  落在那里的节点会被同一个宽度的核摊成一坨，再被边界反射折回来加倍，使 0 附近的高度变成
  “有多少节点低于下限”的读数，而非近期事件强度的读数。因此密度只在

      cut = floor / alpha

  以上的节点参与构核（核仍在年龄 0 处镜像反射，积分恒为 1）。由于最近的参与节点距离现在还有
  好几个核宽，曲线在接近现在的一端自然衰减到零，脊线仍连续画到 0 kya。比 cut 更年轻的节点改
  按计数导出比例，标在纵轴的类别名上——只裁剪显示区间是不够的，那些节点的质量本来就外溢到
  cut 以外。

  同一个维度的全部类别、全部时间窗口共用同一个带宽：密度只算一次（覆盖全部年龄），时间窗口
  只是同一条曲线的不同裁剪，因此同一组在不同窗口的图里峰形完全一致。

定标
  Standardized density per 1,000 samples
      height(a) = f(a) * n_nodes / n_samples * 1000
      n_nodes 只计可分辨范围内的节点，故曲线积分等于该组落在可分辨范围内的每千人事件数。
      它比较的是**合并事件在时间上的集中程度**，不是事件总量，峰高也不能直接解读为人口扩张
      强度或有效群体大小的变化。

      抽稀时展示的是各重复曲线的逐点中位数。中位数与积分不可交换，中位曲线因此不再自动满足
      上面的定标，实测偏差可达 3%。既然纵轴的单位就是这个积分，展示曲线必须守住它：取完
      逐点中位数后，整条曲线按一个标量重标定回 median(定标基准)，区间带用同一个系数。
      常数倍缩放不改变峰形与相对高度，只恢复面积的物理含义。

计算实现
  抽稀有几百次重复，逐次做核密度太慢。这里把年龄先落到固定分箱上，预先算好一张
  “分箱中心 → 网格”的核矩阵，每次重复的密度就退化成一次矩阵乘法。分箱宽度远小于带宽下限。

产物
  output/result/3-ridge_density/⭐山脊密度长表.tsv

CLI
  python python/3-1-ridge_density.py --conf ... --node-ages ... --output-result ...
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm

log = logging.getLogger(__name__)

OUTPUT_NAME = "⭐山脊密度长表.tsv"
LN10 = math.log(10.0)

WITHIN_ALL = "Within-group subtree (all samples)"
WITHIN_RAREFIED = "Within-group subtree (rarefied)"
HEIGHT_UNIT = "Coalescent events per 1,000 sampled individuals per kyr"

# 区间的含义必须随数据一起走，不能只写在图注里：本流程读的是一棵固定的时间树，
# 区间来自抽稀重复之间的差异，**不含**后验定年不确定性。
UNCERTAINTY_RAREFIED = ("Rarefaction replicates on a single fixed time tree; "
                        "posterior dating uncertainty not propagated")
UNCERTAINTY_NONE = ("Single fixed time tree without rarefaction; "
                    "no interval estimated, ridge heights not comparable across groups")


def _silverman(values: np.ndarray) -> float:
    sd = float(np.std(values))
    low, high = np.percentile(values, [25, 75])
    iqr = high - low
    spread = min(sd, iqr / 1.349) if iqr > 0 else sd
    return 0.9 * spread * values.size ** (-0.2)


def _resolution_cut(spec: str, alpha: float, floor: float) -> float:
    """可分辨下界：带宽由下限而非数据决定的那段年龄不参与密度计算。"""
    return floor / alpha if str(spec).strip().lower() == "auto" else float(spec)


def _build_grid(windows, root_age: float, grid_points: int, alpha: float,
                floor: float) -> np.ndarray:
    """一条覆盖全部年龄的评估网格：每个时间窗口内部加密，窗口之外只需支撑积分自检。"""
    pieces = [np.linspace(0.0, high, grid_points) for _, _, high in windows]
    tail_start = max(high for _, _, high in windows)
    tail_stop = root_age * (1 + 6 * alpha) + 6 * floor
    pieces.append(np.linspace(tail_start, tail_stop, 240))
    return np.unique(np.concatenate(pieces))


def _build_bins(root_age: float, cut: float) -> np.ndarray:
    """分箱中心：近期用 0.02 kya，20 kya 之后用 0.5 kya（那里带宽已远大于分箱宽度）。"""
    fine = np.arange(cut, 20.0, 0.02) + 0.01
    coarse = np.arange(20.0, root_age * 1.05, 0.5) + 0.25
    return np.concatenate([fine, coarse])


def _kernel_matrix(grid: np.ndarray, bins: np.ndarray, alpha: float, floor: float) -> np.ndarray:
    """核矩阵 K[i, j]：分箱中心 j 的一个单位质量在网格点 i 上贡献的密度（含 0 处镜像反射）。

    每个核在 [0, inf) 上的积分恒为 1，与 h 无关，因此不同年龄的节点权重完全相同，
    密度也不会漏到负年龄。
    """
    bandwidth = np.maximum(alpha * bins, floor)
    matrix = np.empty((grid.size, bins.size))
    for start in range(0, bins.size, 512):
        block = slice(start, start + 512)
        h = bandwidth[None, block]
        matrix[:, block] = (norm.pdf((grid[:, None] - bins[None, block]) / h)
                            + norm.pdf((grid[:, None] + bins[None, block]) / h)) / h
    return matrix


def _bin_counts(ages: np.ndarray, bins: np.ndarray) -> np.ndarray:
    edges = np.concatenate([[bins[0] - 0.01], (bins[:-1] + bins[1:]) / 2, [bins[-1] * 2]])
    counts, _ = np.histogram(ages, bins=edges)
    return counts


def _parse_windows(spec: str):
    windows = []
    for item in spec.split(","):
        name, _, bounds = item.strip().partition(":")
        low, _, high = bounds.partition("-")
        windows.append((name, float(low), float(high)))
    return windows


def run(conf: str, node_ages: str, rarefied_node_ages: str, sample_counts: str,
        group_design: str, output_result: str) -> int:
    settings = yaml.safe_load(Path(conf).read_text(encoding="utf-8"))
    analysis = settings["analysis"]
    result_path = Path(output_result)
    result_path.mkdir(parents=True, exist_ok=True)

    windows = _parse_windows(analysis["time_windows"])
    window_label = {name: f"Last {high:g} kyr" for name, _, high in windows}
    grid_points = int(analysis["grid_points"])
    floor = float(analysis["bandwidth_floor_kya"])
    factor = float(analysis["bandwidth_factor"])
    low_pct, high_pct = [float(x) for x in str(analysis["interval_percentiles"]).split(",")]
    recent_marks = [float(x) for x in str(analysis["recent_windows_kya"]).split(",")]
    recent_columns = [f"Share of events within {x:g} kyr" for x in recent_marks]

    df1 = pd.read_csv(node_ages, sep="\t")
    df2 = pd.read_csv(sample_counts, sep="\t")
    df3 = pd.read_csv(group_design, sep="\t")
    denominator = {(r.Grouping, r.Group): r.Samples for r in df2.itertuples(index=False)}
    order = df3.sort_values("Order").groupby("Grouping")["Group"].apply(list).to_dict()

    rarefied_path = Path(rarefied_node_ages)
    use_rarefaction = rarefied_path.exists()
    df4 = pd.read_csv(rarefied_path, sep="\t") if use_rarefaction else None
    root_age = float(df1["Node age (kya)"].max())
    log.info("全样本节点 %d 行，%s，最深节点 %.1f kya", len(df1),
             f"抽稀节点 {len(df4)} 行" if use_rarefaction else "抽稀已关闭", root_age)

    records, renorm_log = [], []
    for grouping, df_group in df1.groupby("Grouping", sort=False):
        levels = order[grouping]

        # 该维度的带宽只算一次：由全样本诱导子树的全部节点年龄决定，
        # 之后该维度的每个类别、每个时间窗口都用它。
        alpha = factor * _silverman(np.log10(df_group["Node age (kya)"].to_numpy())) * LN10
        cut = _resolution_cut(analysis["resolution_cut"], alpha, floor)
        grid = _build_grid(windows, root_age, grid_points, alpha, floor)
        bins = _build_bins(root_age, cut)
        kernel = _kernel_matrix(grid, bins, alpha, floor)
        log.info("  %-30s 相对带宽=%.4f，下限=%.2f kya，可分辨下界=%.2f kya，网格 %d 点，分箱 %d 个",
                 grouping, alpha, floor, cut, grid.size, bins.size)

        for category in levels:
            n_samples = int(denominator[(grouping, category)])
            all_ages = df_group.loc[df_group["Group"] == category, "Node age (kya)"].to_numpy()
            if all_ages.size == 0:
                continue

            # 某个维度若没有抽稀数据（单类别的全样本总览，或全局关闭抽稀），
            # 就用全部节点直接出曲线
            subset = (df4[(df4["Grouping"] == grouping) & (df4["Group"] == category)]
                      if use_rarefaction else None)
            if subset is not None and not subset.empty:
                n_replicates = int(subset["Replicate"].max()) + 1
                counts = np.zeros((bins.size, n_replicates))
                per_replicate_nodes = np.zeros(n_replicates)
                resolved_nodes = np.zeros(n_replicates)
                share_window = {name: np.zeros(n_replicates) for name, _, _ in windows}
                share_recent = {col: np.zeros(n_replicates) for col in recent_columns}
                share_below = np.zeros(n_replicates)
                for replicate, block in subset.groupby("Replicate", sort=False):
                    ages = block["Node age (kya)"].to_numpy()
                    keep = ages[ages >= cut]
                    counts[:, replicate] = _bin_counts(keep, bins)
                    per_replicate_nodes[replicate] = ages.size
                    resolved_nodes[replicate] = keep.size
                    share_below[replicate] = 1.0 - keep.size / ages.size
                    for name, _, high in windows:
                        share_window[name][replicate] = (ages <= high).mean()
                    for col, mark in zip(recent_columns, recent_marks):
                        share_recent[col][replicate] = (ages <= mark).mean()
                if resolved_nodes.min() == 0:
                    continue
                target = int(per_replicate_nodes[0]) + 1

                # 每次重复各自按自己的可分辨节点数定标，被排除的比例因此体现为曲线下面积的缺口
                density = kernel @ counts / resolved_nodes[None, :]
                height_matrix = density * resolved_nodes[None, :] / target * 1000.0

                # 每一次重复的曲线积分都严格等于它自己的定标基准，这是核归一化的直接推论
                per_replicate_basis = resolved_nodes / target * 1000.0
                worst = float(np.max(np.abs(np.trapezoid(height_matrix, grid, axis=0)
                                            - per_replicate_basis) / per_replicate_basis))
                if worst > 0.01:
                    raise AssertionError(
                        f"{grouping}/{category}: 单次重复的全网格积分最大偏离定标基准 "
                        f"{worst:.1%}，超过 1%")

                height = np.median(height_matrix, axis=1)
                lower = np.percentile(height_matrix, low_pct, axis=1)
                upper = np.percentile(height_matrix, high_pct, axis=1)

                # 逐网格点取中位数与求积分不可交换：每个重复的曲线都严格归一，但它们的逐点
                # 中位数不再是。实测偏差可达 3%，而纵轴的单位是「每千人每千年的事件数」，
                # 面积必须守恒，否则这个单位就不成立。因此把整条中位曲线按一个标量重标定到
                # median(定标基准)。缩放是常数倍，峰形与相对高度不受影响；区间带用同一个
                # 系数，使带宽相对中位线的比例保持不变。
                basis = float(np.median(per_replicate_basis))
                total = float(np.trapezoid(height, grid))
                renorm = basis / total
                height, lower, upper = height * renorm, lower * renorm, upper * renorm
                renorm_log.append((grouping, category, abs(renorm - 1.0)))
                if abs(renorm - 1.0) > 0.10:
                    raise AssertionError(
                        f"{grouping}/{category}: 中位曲线积分 {total:.3f} 与定标基准 "
                        f"{basis:.3f} 相差 {abs(renorm - 1) :.1%}，超过 10%，"
                        f"重标定会掩盖真实问题，请先检查抽稀重复之间的离散程度")
                shared = {
                    "Node definition": WITHIN_RAREFIED,
                    "Rarefied sample size": target,
                    "Rarefaction replicates": n_replicates,
                    "Coalescent events": target - 1,
                    "Median node age (kya)": float(np.median(subset["Node age (kya)"].to_numpy())),
                    "Uncertainty source": UNCERTAINTY_RAREFIED,
                    "Share of events below the resolution limit": float(np.median(share_below)),
                    "Share of events below the resolution limit, lower":
                        float(np.percentile(share_below, low_pct)),
                    "Share of events below the resolution limit, upper":
                        float(np.percentile(share_below, high_pct)),
                    **{col: float(np.median(v)) for col, v in share_recent.items()},
                }
                window_share = {name: (float(np.median(share_window[name])),
                                       float(np.percentile(share_window[name], low_pct)),
                                       float(np.percentile(share_window[name], high_pct)))
                                for name, _, _ in windows}
            else:
                keep = all_ages[all_ages >= cut]
                if keep.size == 0:
                    continue
                n_nodes = float(keep.size)
                density = kernel @ _bin_counts(keep, bins) / n_nodes
                height = density * n_nodes / n_samples * 1000.0
                lower = upper = np.full(grid.size, np.nan)
                basis = n_nodes / n_samples * 1000.0
                total = float(np.trapezoid(height, grid))
                if abs(total - basis) / basis > 0.01:
                    raise AssertionError(
                        f"{grouping}/{category}: 全网格积分 {total:.3f} 偏离定标基准 "
                        f"{basis:.3f} 超过 1%")
                shared = {
                    "Node definition": WITHIN_ALL,
                    "Rarefied sample size": pd.NA,
                    "Rarefaction replicates": pd.NA,
                    "Coalescent events": n_nodes,
                    "Median node age (kya)": float(np.median(all_ages)),
                    "Uncertainty source": UNCERTAINTY_NONE,
                    "Share of events below the resolution limit": 1.0 - keep.size / all_ages.size,
                    "Share of events below the resolution limit, lower": np.nan,
                    "Share of events below the resolution limit, upper": np.nan,
                    **{col: float((all_ages <= mark).mean())
                       for col, mark in zip(recent_columns, recent_marks)},
                }
                window_share = {name: (float((all_ages <= high).mean()), np.nan, np.nan)
                                for name, _, high in windows}

            for name, _, high in windows:
                inside = grid <= high
                share, share_low, share_high = window_share[name]
                records.append(pd.DataFrame({
                    "Time window": window_label[name],
                    "Grouping": grouping,
                    "Group": category,
                    "Node age (kya)": grid[inside],
                    "Ridge height": height[inside],
                    "Ridge height lower": lower[inside],
                    "Ridge height upper": upper[inside],
                    "Ridge height unit": HEIGHT_UNIT,
                    "Samples": n_samples,
                    "Relative bandwidth": alpha,
                    "Bandwidth floor (kya)": floor,
                    "Resolution limit (kya)": cut,
                    "Share of events within window": share,
                    "Share of events within window, lower": share_low,
                    "Share of events within window, upper": share_high,
                    **shared,
                }))

    if renorm_log:
        worst = max(renorm_log, key=lambda x: x[2])
        log.info("中位曲线重标定：最大修正 %.2f%%（%s / %s），中位修正 %.2f%%",
                 worst[2] * 100, worst[0], worst[1],
                 float(np.median([x[2] for x in renorm_log])) * 100)

    df5 = pd.concat(records, ignore_index=True)
    output_file = result_path / OUTPUT_NAME
    df5.to_csv(output_file, sep="\t", index=False, float_format="%.8g")
    log.info("积分自检通过：每个抽稀重复的曲线、以及重标定后的中位曲线，"
             "积分均等于其定标基准")
    log.info("山脊密度长表 %d 行 → %s", len(df5), output_file)

    column = "Share of events below the resolution limit"
    excluded = df5.drop_duplicates(["Grouping", "Group"]).nlargest(3, column)
    for _, row in excluded.iterrows():
        log.info("  低于可分辨下界 %.2f kya 的事件比例最高：%s / %s = %.1f%%",
                 row["Resolution limit (kya)"], row["Grouping"], row["Group"], row[column] * 100)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="计算组间可比的山脊密度高度")
    parser.add_argument("--conf", required=True)
    parser.add_argument("--node-ages", required=True)
    parser.add_argument("--rarefied-node-ages", required=True)
    parser.add_argument("--sample-counts", required=True)
    parser.add_argument("--group-design", required=True)
    parser.add_argument("--output-result", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    return run(**vars(build_parser().parse_args(argv)))


if __name__ == "__main__":
    raise SystemExit(main())

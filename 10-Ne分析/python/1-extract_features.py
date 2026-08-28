#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1-extract_features.py — 读入 Skygrid 后验，按时间段与时间 bin 提取 8 个轨迹特征。
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from loaders import load_data  # noqa: E402
from ne_core import (FEATURE_DISPLAY, FEATURE_KEYS, apply_size_adjustment,  # noqa: E402
                     at_window_boundary, directional_probability, hpd_interval,
                     segment_features, size_adjustment_fit, window_mask)
from pipeline_utils import load_config, resolve_paths, setup_logging, write_tsv  # noqa: E402

log = logging.getLogger(__name__)

HEADER = ["Category", "Group", "Window Type", "Time Window", "Window Start (Years)",
          "Window End (Years)", "Feature", "Median", "Lower 95% HPD", "Upper 95% HPD",
          "Directional Posterior Probability", "Probability Statement",
          "Reliable Window Coverage", "Usable", "At Window Boundary"]


def build_windows(periods: list, period_width: float, period_max: float,
                  bin_width: float, analysis_max: float
                  ) -> list[tuple[str, str, float, float]]:
    """时间窗清单：时间段 + 覆盖整条时间轴的等宽 bin。

    时间段默认按 period_width 等宽切分到 period_max（配置里显式给了 time_periods 就用
    那一套，不再自动切）；bin 用于全时间轴的合并分数热图。两类窗口的边界可以完全相同，
    因此下游一律用 (窗口类型, 标签) 作键。
    """
    out: list[tuple[str, str, float, float]] = []
    if periods:
        out += [("Time Period", f"{t0:g}-{t1:g} Years", float(t0), float(t1))
                for t0, t1 in periods]
    else:
        t = 0.0
        while t + period_width <= period_max + 1e-9:
            out.append(("Time Period", f"{t:g}-{t + period_width:g} Years",
                        t, t + period_width))
            t += period_width
    t = 0.0
    while t + bin_width <= analysis_max + 1e-9:
        out.append(("Time Bin", f"{t:g}-{t + bin_width:g} Years", t, t + bin_width))
        t += bin_width
    return out


def extract_features_by_period(curves, windows: list, min_coverage: float,
                               cred_level: float) -> tuple[list[list], dict]:
    """逐类别、逐时间窗提取特征，返回 (长表行, 中位数特征矩阵字典)。

    特征矩阵按 {(窗口类型, 窗口标签): {类别: 长度 8 的中位数向量}} 组织，供第 2 步直接
    取用。键里必须带窗口类型：时间段与时间 bin 的边界可以完全相同，只用标签当键会让
    两者互相覆盖，特征矩阵与窗口清单随即错位。
    """
    rows: list[list] = []
    matrices: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for wtype, wname, t0, t1 in windows:
        mask = window_mask(curves.grid, t0, t1)
        if mask.sum() < 3:
            log.warning("%s 窗内不足 3 个网格点，跳过", wname)
            continue
        wkey = (wtype, wname)
        matrices[wkey] = {}
        for cat in curves.categories:
            feats = segment_features(curves.curves[cat], curves.rates[cat],
                                     curves.grid, mask)
            coverage = curves.coverage(mask, cat)
            usable = "Yes" if coverage >= min_coverage else "No"
            matrices[wkey][cat] = np.array([np.median(feats[k]) for k in FEATURE_KEYS])
            for key in FEATURE_KEYS:
                draws = feats[key]
                lo, hi = hpd_interval(draws, cred_level)
                med = float(np.median(draws))
                prob, statement = directional_probability(key, feats)
                rows.append([
                    cat, curves.group[cat], wtype, wname, f"{t0:g}", f"{t1:g}",
                    FEATURE_DISPLAY[key], f"{med:.5f}", f"{lo:.5f}", f"{hi:.5f}",
                    "" if not np.isfinite(prob) else f"{prob:.3f}", statement,
                    f"{coverage:.3f}", usable,
                    at_window_boundary(key, med, curves.grid[mask])])
    return rows, matrices


def run(config: str, project_root: str = ".", log_file: str | None = None) -> int:
    cfg = resolve_paths(load_config(config), project_root)
    setup_logging(cfg["runtime"]["log_level"], log_file)
    io, tm, po, se = cfg["io"], cfg["time"], cfg["posterior"], cfg["select"]
    unit = float(tm["years_per_unit"]) if str(tm["time_unit"]) != "years" else 1.0

    curves = load_data(
        input_dir=io["input_dir"], grid_step=float(tm["grid_step"]),
        analysis_max=float(tm["analysis_max"]), smooth_bandwidth=float(tm["smooth_bandwidth"]),
        input_format=io["input_format"], samples_table=io["samples_table"],
        years_per_unit=unit, burnin_frac=float(po["burnin_frac"]),
        max_draws=int(po["max_draws"]), hpd_width_max=float(se["hpd_width_max"]),
        cred_level=float(po["cred_level"]),
        selected_categories=str(se["selected_categories"]),
        cutoff_override=float(tm["cutoff_override"]), seed=int(po["seed"]))

    # 规模校正：参数在整条分析窗上只拟合一次，全部时间窗共用同一套
    if se["size_adjustment"]:
        sizes = np.array([curves.n_tips.get(c, np.nan) for c in curves.categories])
        if np.isfinite(sizes).sum() < 3:
            log.warning("开启了规模校正但可用样本量不足 3 个，已跳过")
        else:
            full = window_mask(curves.grid, 0.0, float(curves.grid[-1]))
            med = np.array([np.median(curves.curves[c][:, full], axis=0).mean()
                            for c in curves.categories])
            slope, intercept = size_adjustment_fit(med, np.log(sizes))
            for i, c in enumerate(curves.categories):
                curves.curves[c] = apply_size_adjustment(
                    curves.curves[c], np.full(curves.curves[c].shape[0], np.log(sizes[i])),
                    slope, intercept)
            log.info("规模校正（全窗一次拟合）：平均 log Ne = %.3f × log n + %.3f",
                     slope, intercept)

    # 时间段与 bin 都不能越过实际读到的时间轴：输入只到 20 kya 时，再切到 50 kya
    # 只会得到由端点值外推出来的假曲线
    axis_max = float(curves.grid[-1])
    windows = build_windows(tm["time_periods"], float(tm["period_width"]),
                            min(float(tm["period_max"]), axis_max),
                            float(tm["bin_width"]), axis_max)
    n_period = sum(1 for w in windows if w[0] == "Time Period")
    log.info("时间窗：%d 个时间段（宽 %g 年）+ %d 个时间 bin（宽 %g 年）",
             n_period, float(tm["period_width"]), len(windows) - n_period,
             float(tm["bin_width"]))
    rows, matrices = extract_features_by_period(curves, windows,
                                                float(se["min_coverage"]),
                                                float(po["cred_level"]))
    out_dir = Path(io["output_dir"])
    n = write_tsv(str(out_dir / "⭐1-Feature-Summary.tsv"), HEADER, rows)
    log.info("特征长表 -> %s（%d 行，%d 个时间窗）",
             out_dir / "⭐1-Feature-Summary.tsv", n, len(matrices))

    kept = [w for w in windows if (w[0], w[1]) in matrices]
    cache = out_dir / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache / "curves.npz",
        grid=curves.grid, categories=np.array(curves.categories, dtype=object),
        group=np.array([curves.group[c] for c in curves.categories], dtype=object),
        has_posterior=curves.has_posterior,
        window_meta=np.array([[w[0], w[1], str(w[2]), str(w[3])] for w in kept],
                             dtype=object),
        feature_matrix=np.array([[matrices[(w[0], w[1])][c] for c in curves.categories]
                                 for w in kept]),
        coverage=np.array([[curves.coverage(window_mask(curves.grid, float(w[2]),
                                                        float(w[3])), c)
                            for c in curves.categories] for w in kept]),
        **{f"curve__{c}": curves.curves[c].astype(np.float32) for c in curves.categories},
        **{f"rate__{c}": curves.rates[c].astype(np.float32) for c in curves.categories})
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="按时间段与时间 bin 提取 Ne 轨迹特征")
    p.add_argument("--config", required=True)
    p.add_argument("--project-root", default=".")
    p.add_argument("--log-file", default=None)
    return p


def main(argv=None) -> int:
    return run(**vars(build_parser().parse_args(argv)))


if __name__ == "__main__":
    raise SystemExit(main())

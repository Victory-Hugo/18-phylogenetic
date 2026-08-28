#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""loaders.py — 四种输入格式的统一入口（不含 main，只被 import）。
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ne_core import (check_grid_against_xml, display_name, loess_operators,
                     read_xml_skygrid, reliable_mask, skygrid_log_popsize,
                     thin_indices)

log = logging.getLogger(__name__)

RATE_UNIT_YEARS = 1000.0     # 速率统一表达为"每千年"，与年为单位的时间轴解耦


@dataclass
class CurveSet:
    """全部类别的统一表示：同一条时间网格、同样的平滑口径。"""

    grid: np.ndarray                                  # 年，升序，0 = 现在
    categories: list[str]
    curves: dict[str, np.ndarray]                     # 平滑后的 log Ne
    rates: dict[str, np.ndarray]                      # d log Ne / d(千年)，现在方向为正
    valid: dict[str, np.ndarray]                      # 可靠时间窗掩码
    group: dict[str, str] = field(default_factory=dict)
    n_tips: dict[str, float] = field(default_factory=dict)
    has_posterior: bool = True
    source_format: str = ""

    def coverage(self, mask: np.ndarray, category: str) -> float:
        """某个时间窗有多大比例落在该类别的可靠窗内。"""
        m = self.valid[category][mask]
        return float(m.mean()) if m.size else 0.0


def tsv_kind(path: Path) -> str:
    """按**表头**判断一张 TSV 是什么：分组清单 / 逐条后验 / 中位数汇总。

    不按文件名判断：清单文件叫什么都行（`samples.tsv`、`0-分组信息-可选.tsv` 都可以），
    靠文件名识别会在用户改名后悄悄把清单当成曲线数据读进来。
    """
    cols = {c.strip().lower() for c in pd.read_csv(path, sep="\t", nrows=1).columns}
    if "time" not in cols:
        return "samples"
    return "posterior_table" if "draw" in cols else "summary_table"


def split_tsv_inputs(input_dir: Path) -> tuple[list[Path], Path | None]:
    """把目录里的 TSV 分成"曲线数据"与"分组清单"两类。"""
    data, meta = [], None
    for f in sorted(input_dir.glob("*.tsv")):
        if tsv_kind(f) == "samples":
            meta = meta or f
        else:
            data.append(f)
    return data, meta


def detect_format(input_dir: Path) -> str:
    """按目录内容猜输入格式，顺序反映"信息量从多到少"。"""
    if list(input_dir.glob("*.log")):
        return "beast_log"
    if list(input_dir.glob("*.npz")):
        return "npz"
    data, _ = split_tsv_inputs(input_dir)
    if data:
        return tsv_kind(data[0])
    raise FileNotFoundError(f"{input_dir} 里没有 .log / .npz / .tsv 输入")


def read_samples_table(path: Path | None) -> pd.DataFrame | None:
    """可选清单：category / group / n_tips，只有需要分组或规模校正时才必要。"""
    if path is None or not Path(path).is_file():
        return None
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]
    if "category" not in df.columns:
        raise ValueError(f"{path} 缺少 category 列")
    return df.set_index("category")


# ------------------------------------------------------------------ 各格式读取


def _read_beast_log(path: Path, burnin_frac: float, max_draws: int) -> np.ndarray:
    """读一个 .log，去 burn-in、等间隔抽稀，返回 logPopSize 矩阵 (n_draws, dim)。"""
    df = pd.read_csv(path, sep="\t", comment="#")
    if df.empty:
        raise ValueError(f"{path.name} 为空")
    if "state" in df.columns and 0 < burnin_frac < 1:
        cut = int(len(df) * burnin_frac)
        df = df[df["state"] >= int(df["state"].iloc[cut])]
    cols = sorted((c for c in df.columns if c.startswith("skygrid.logPopSize")),
                  key=lambda c: int(c.replace("skygrid.logPopSize", "") or 0))
    if not cols:
        raise KeyError(f"{path.name} 里找不到 skygrid.logPopSize 列")
    df = df.iloc[thin_indices(len(df), max_draws)]
    return df[cols].to_numpy(dtype=float)


def _from_beast_logs(input_dir: Path, grid: np.ndarray, burnin_frac: float,
                     max_draws: int, years_per_unit: float,
                     cutoff_override: float) -> tuple[dict[str, np.ndarray], bool, float]:
    """input/*.log + 同名 *.xml → 每个类别的 log Ne 后验矩阵。"""
    out = {}
    axis_max = np.inf
    for f in sorted(input_dir.glob("*.log")):
        values = _read_beast_log(f, burnin_frac, max_draws)
        xml_info = read_xml_skygrid(f.with_suffix(".xml"))
        if xml_info is None and cutoff_override > 0:
            xml_info = {"cutoff": cutoff_override / years_per_unit,
                        "dimension": float(values.shape[1]),
                        "n_grid_points": float(values.shape[1] - 1)}
            log.warning("%s 没有可用 XML，改用配置里的 cutoff_override", f.name)
        # cutOff 与网格都在"树的时间单位"里，先换算成年再和分析网格比较
        cutoff_years = check_grid_against_xml(
            f.name, values.shape[1],
            None if xml_info is None else {**xml_info,
                                           "cutoff": xml_info["cutoff"] * years_per_unit},
            float(grid[-1]))
        out[f.stem] = skygrid_log_popsize(values, cutoff_years, grid)
        axis_max = min(axis_max, cutoff_years)
        log.info("%s：%d 条后验样本，cutOff %.0f 年", f.stem, values.shape[0], cutoff_years)
    return out, True, float(axis_max)


def _from_npz(input_dir: Path, grid: np.ndarray, max_draws: int
              ) -> tuple[dict[str, np.ndarray], dict[str, dict], bool, float]:
    """上游后验缓存：log_ne (n_draws, n_grid) + grid_kya + 可选 meta。"""
    curves, meta_all = {}, {}
    axis_max = np.inf
    for f in sorted(input_dir.glob("*.npz")):
        with np.load(f, allow_pickle=True) as z:
            raw = z["log_ne" if "log_ne" in z.files else "log_ne_sm"].astype(float)
            src_grid = z["grid_kya"].astype(float) * 1000.0 if "grid_kya" in z.files \
                else z["grid_years"].astype(float)
            meta = dict(z["meta"].item()) if "meta" in z.files else {}
        raw = raw[thin_indices(raw.shape[0], max_draws)]
        name = str(meta.get("stratum_id") or meta.get("category") or f.stem)
        curves[name] = np.column_stack([np.interp(grid, src_grid, row) for row in raw]).T
        meta_all[name] = meta
        axis_max = min(axis_max, float(src_grid.max()))
        log.info("%s：%d 条后验样本（缓存）", name, raw.shape[0])
    return curves, meta_all, True, float(axis_max)


def _from_posterior_table(files: list[Path], grid: np.ndarray, max_draws: int
                          ) -> tuple[dict[str, np.ndarray], bool, float]:
    """逐条后验长表：draw / time / log_ne，category 可以来自列，也可以来自文件名。

    一类别一文件是更常见的组织方式（文件名即类别），单文件多类别也照样支持。
    """
    curves = {}
    axis_max = np.inf
    for path in files:
        df = pd.read_csv(path, sep="\t")
        df.columns = [c.strip().lower() for c in df.columns]
        if "category" not in df.columns:
            df["category"] = path.stem
        for cat, g in df.groupby("category"):
            wide = g.pivot_table(index="draw", columns="time", values="log_ne")
            keep = thin_indices(len(wide), max_draws)
            src_grid = wide.columns.to_numpy(dtype=float)
            mat = wide.to_numpy(dtype=float)[keep]
            curves[str(cat)] = np.column_stack(
                [np.interp(grid, src_grid, r) for r in mat]).T
            axis_max = min(axis_max, float(src_grid.max()))
            log.info("%s：%d 条后验样本（长表）", cat, mat.shape[0])
    return curves, True, float(axis_max)


def _from_summary_table(files: list[Path], grid: np.ndarray, n_pseudo: int, seed: int
                        ) -> tuple[dict[str, np.ndarray], bool, float]:
    """中位数 + 区间汇总表 → 合成伪后验（整条曲线共用一个位移）。

    没有真正的后验样本时，逐时间点独立抽样会造出锯齿状的假曲线，速率与事件时间会被
    彻底污染。这里改为每条伪样本只抽一个标准正态位移、乘以逐点标准差后整体平移，
    形状严格保留，只有水平带上不确定性。这是近似，不是后验。
    """
    rng = np.random.default_rng(seed)
    curves = {}
    axis_max = np.inf
    for path in files:
        df = pd.read_csv(path, sep="\t")
        df.columns = [c.strip().lower() for c in df.columns]
        if "category" not in df.columns:
            df["category"] = path.stem
        lo_col = next((c for c in df.columns if "lower" in c), None)
        hi_col = next((c for c in df.columns if "upper" in c), None)
        for cat, g in df.groupby("category"):
            g = g.sort_values("time")
            src = g["time"].to_numpy(dtype=float)
            med = np.interp(grid, src, g["median"].to_numpy(dtype=float))
            if lo_col and hi_col:
                sd = np.interp(grid, src, (g[hi_col].to_numpy(dtype=float)
                                           - g[lo_col].to_numpy(dtype=float)) / 3.92)
            else:
                sd = np.zeros_like(med)
            z = rng.standard_normal(n_pseudo)[:, None]
            curves[str(cat)] = med[None, :] + z * sd[None, :]
            axis_max = min(axis_max, float(src.max()))
    return curves, False, float(axis_max)


# ------------------------------------------------------------------ 统一入口


def load_data(input_dir: str, grid_step: float, analysis_max: float,
              smooth_bandwidth: float, input_format: str = "auto",
              samples_table: str = "", years_per_unit: float = 1.0,
              burnin_frac: float = 0.1, max_draws: int = 6000,
              hpd_width_max: float = 2.0, cred_level: float = 0.95,
              selected_categories: str = "", cutoff_override: float = 0.0,
              seed: int = 20260828) -> CurveSet:
    """把任意一种输入读成同一套曲线：统一网格、统一平滑、统一速率单位。

    时间一律换算成年（0 = 现在）。平滑与求导用同一次 loess 局部二次拟合，速率的单位是
    "每千年"，与网格步长无关。
    """
    src = Path(input_dir)
    fmt = input_format if input_format not in ("", "auto") else detect_format(src)
    grid = np.arange(0.0, analysis_max + 1e-9, grid_step)
    meta_all: dict[str, dict] = {}

    if fmt == "beast_log":
        curves, has_post, axis_max = _from_beast_logs(src, grid, burnin_frac, max_draws,
                                                      years_per_unit, cutoff_override)
    elif fmt == "npz":
        curves, meta_all, has_post, axis_max = _from_npz(src, grid, max_draws)
    elif fmt in ("posterior_table", "summary_table"):
        files, _ = split_tsv_inputs(src)
        if not files:
            raise FileNotFoundError(f"{src} 里没有可读的长表")
        if fmt == "posterior_table":
            curves, has_post, axis_max = _from_posterior_table(files, grid, max_draws)
        else:
            curves, has_post, axis_max = _from_summary_table(files, grid, 200, seed)
    else:
        raise ValueError(f"未知的 input_format：{fmt}")
    if not curves:
        raise RuntimeError(f"{src} 没有读出任何类别")

    # 输入本身没覆盖到的年代必须截掉：np.interp 会把末端值一路平推出去，
    # 那段"平稳的 Ne"完全是插值的产物，不是数据
    if np.isfinite(axis_max) and axis_max < grid[-1] - 1e-9:
        keep_t = grid <= axis_max + 1e-9
        log.warning("输入只覆盖到 %g 年，分析上界由 %g 年截到该处", axis_max, grid[-1])
        grid = grid[keep_t]
        curves = {k: v[:, keep_t] for k, v in curves.items()}

    wanted = {c.strip() for c in selected_categories.split(",") if c.strip()}
    if wanted:
        missing = wanted - set(curves)
        if missing:
            raise KeyError(f"selected_categories 里这些类别不存在：{sorted(missing)}")
        curves = {k: v for k, v in curves.items() if k in wanted}

    # 分组清单：配置里指定优先，否则在输入目录里按表头自动认领
    meta_path = Path(samples_table) if samples_table else split_tsv_inputs(src)[1]
    samples = read_samples_table(meta_path)
    if samples is not None:
        log.info("分组清单：%s（%d 条）", Path(meta_path).name, len(samples))
    h0, h1 = loess_operators(grid, smooth_bandwidth)
    smoothed, rates, valid, group, n_tips = {}, {}, {}, {}, {}
    fallback_named = []
    for cat, mat in curves.items():
        meta = meta_all.get(cat, {})
        g = str(meta.get("stratum_class", "") or "")
        n = float(meta.get("n_tips", float("nan")) or float("nan"))
        name = ""
        if samples is not None and cat in samples.index:
            row = samples.loc[cat]
            name = str(row.get("display_name", "") or "")
            g = str(row.get("group", g) or g)
            n = float(row.get("n_tips", n) or n)
        if not name:
            name = display_name(cat)
            fallback_named.append(cat)
        if name in smoothed:
            raise ValueError(f"展示名 {name} 对应了多个类别，请在 samples.tsv 里区分开")
        smoothed[name] = mat @ h0.T
        # 网格按年代递增，取负号后现在方向为正；单位统一为"每千年"
        rates[name] = -(mat @ h1.T) * RATE_UNIT_YEARS
        valid[name] = reliable_mask(smoothed[name], hpd_width_max, cred_level)
        group[name] = g or "All Categories"
        n_tips[name] = n
    if fallback_named:
        log.info("%d 个类别没有 display_name，已按下划线换空格自动命名：%s",
                 len(fallback_named), "、".join(fallback_named[:5]))

    cats = sorted(smoothed)
    if not has_post:
        log.warning("输入不含后验样本，已按整体位移合成伪样本：区间宽度可用，"
                    "但同簇频率不再是后验概率，报告中须写明")
    log.info("输入格式 %s：%d 个类别，网格 0–%g 年、步长 %g 年，平滑带宽 %g 年",
             fmt, len(cats), grid[-1], grid_step, smooth_bandwidth)
    return CurveSet(grid=grid, categories=cats, curves=smoothed, rates=rates,
                    valid=valid, group=group, n_tips=n_tips,
                    has_posterior=has_post, source_format=fmt)

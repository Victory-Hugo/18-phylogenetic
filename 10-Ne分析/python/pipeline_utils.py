#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pipeline_utils.py — 配置读取、日志与写表（不含 main，只被 import）。"""
from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path
from typing import Any

import yaml

DEFAULTS: dict[str, dict[str, Any]] = {
    "io": {"input_dir": "input", "output_dir": "output", "input_format": "auto",
           "samples_table": ""},
    "time": {"time_unit": "years", "years_per_unit": 1.0, "time_periods": [],
             "period_width": 500.0, "period_max": 50000.0, "bin_width": 500.0,
             "analysis_max": 50000.0, "grid_step": 200.0,
             "smooth_bandwidth": 3000.0, "cutoff_override": 0.0},
    "posterior": {"burnin_frac": 0.1, "max_draws": 4000, "n_resample": 500,
                  "cred_level": 0.95, "seed": 20260828},
    "select": {"group_var": "", "selected_categories": "", "min_coverage": 0.8,
               "hpd_width_max": 2.0, "size_adjustment": False},
    "cluster": {"distance_metric": "euclidean", "linkage_method": "average",
                "n_clusters": 0, "k_max": 8, "feature_scaling": "zscore"},
    "figure": {"scatter_method": "pca", "chord_metric": "coclustering",
               "chord_threshold": 0.5, "heatmap_palette": "diverging",
               "score_palette": "diverging", "scatter_ellipse": True,
               "panels_per_row": 3},
    "runtime": {"debug_mode": False, "log_level": "INFO"},
    # 留空 = 自动认领（优先当前 conda 环境）；跨机器跑时可以在这里写死绝对路径
    "tools": {"python_bin": "", "rscript_bin": ""},
}


def setup_logging(level: str = "INFO", log_file: str | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(level=getattr(logging, str(level).upper(), logging.INFO),
                        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                        handlers=handlers, force=True)


def load_config(path: str) -> dict[str, dict[str, Any]]:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    cfg = {sec: dict(vals) for sec, vals in DEFAULTS.items()}
    for sec, vals in raw.items():
        if sec not in cfg:
            raise KeyError(f"配置里有未知的段落 {sec}，可用的是 {sorted(cfg)}")
        unknown = set(vals) - set(cfg[sec])
        if unknown:
            raise KeyError(f"{sec} 段里有未知的键 {sorted(unknown)}")
        cfg[sec].update(vals)
    if not cfg["time"]["time_periods"] and float(cfg["time"]["period_width"]) <= 0:
        raise ValueError("要么显式给出 time.time_periods，要么给一个正的 time.period_width")
    return cfg


def resolve_paths(cfg: dict, base: str) -> dict:
    """把配置里的相对路径按配置文件所在项目根目录展开。"""
    root = Path(base).resolve()
    for key in ("input_dir", "output_dir", "samples_table"):
        v = cfg["io"][key]
        if v and not Path(v).is_absolute():
            cfg["io"][key] = str(root / v)
    return cfg


def write_tsv(path: str, header: list[str], rows: list[list]) -> int:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return len(rows)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""selfcheck.py — 关键实现的最小测试，不依赖任何真实数据。

六项，全部用人工构造的最小样例：

  Test 1  Skygrid 分段边界归属。时间恰好落在内部边界上时必须归入**后**一段。
  Test 2  XML 与日志的一致性守卫。维度对不上、分析上界超出 cutOff 都必须报错。
  Test 3  事件方向语义。整窗只收缩时扩张支持度必须为 0、收缩为 1；水平不给概率。
  Test 4  跨时间 bin 的标签对齐。只在两个 bin 共有的类别上匹配，缺席者不参与。
  Test 5  Ward 的适用范围。非欧氏度量必须被挡在 Ward 之外。
  Test 6  端到端。合成三个类别的 mini .log + .xml，跑通特征提取并检查关键列。

    python python/selfcheck.py
返回码 0 表示全部通过。
"""
from __future__ import annotations

import argparse
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ne_core import (align_cluster_labels_across_bins, check_grid_against_xml,  # noqa: E402
                     directional_probability, linkage_is_valid, read_xml_skygrid,
                     segment_features, skygrid_log_popsize)
from pipeline_utils import setup_logging  # noqa: E402

log = logging.getLogger(__name__)

XML_TEMPLATE = """<beast>
  <parameter id="skygrid.logPopSize" dimension="{dim}" value="1.0"/>
  <parameter id="skygrid.numGridPoints" value="{ngp}"/>
  <parameter id="skygrid.cutOff" value="{cutoff}"/>
</beast>
"""


def _write_fake_run(dirpath: Path, name: str, level: float, trend: float,
                    dim: int = 6, cutoff: float = 24000.0, n_draws: int = 60) -> None:
    """造一个 mini BEAST 输出：.log 里是 logPopSize 的后验，.xml 里是时间轴定义。"""
    rng = np.random.default_rng(abs(hash(name)) % 2**32)
    # 列序号 1 = 最近的一段；trend > 0 表示越靠近现在 Ne 越高
    base = level + trend * np.linspace(1.0, 0.0, dim)
    draws = base[None, :] + rng.normal(0, 0.05, size=(n_draws, dim))
    df = pd.DataFrame(draws, columns=[f"skygrid.logPopSize{i + 1}" for i in range(dim)])
    df.insert(0, "state", np.arange(n_draws) * 1000)
    df.to_csv(dirpath / f"{name}.log", sep="\t", index=False)
    (dirpath / f"{name}.xml").write_text(
        XML_TEMPLATE.format(dim=dim, ngp=dim - 1, cutoff=cutoff), encoding="utf-8")


def test_skygrid_boundary() -> tuple[bool, str]:
    """dimension=3 → numGridPoints=2，cutOff=10000，边界在 5000 与 10000。"""
    values = np.array([[1.0, 2.0, 3.0]])
    t = np.array([0, 4999, 5000, 5001, 9999, 10000, 10001], dtype=float)
    got = skygrid_log_popsize(values, 10000.0, t)[0]
    want = np.array([1, 1, 2, 2, 2, 3, 3], dtype=float)
    return np.array_equal(got, want), f"边界归属 {got.astype(int).tolist()}"


def test_xml_guard() -> tuple[bool, str]:
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "a.xml"
        p.write_text(XML_TEMPLATE.format(dim=6, ngp=5, cutoff=24000), encoding="utf-8")
        info = read_xml_skygrid(p)
        ok_read = info is not None and abs(info["cutoff"] - 24000) < 1e-9
        ok_pass = abs(check_grid_against_xml("a", 6, info, 20000.0) - 24000.0) < 1e-9
        try:
            check_grid_against_xml("a", 5, info, 20000.0)
            ok_dim = False
        except ValueError:
            ok_dim = True
        try:
            check_grid_against_xml("a", 6, info, 30000.0)
            ok_cut = False
        except ValueError:
            ok_cut = True
        try:
            check_grid_against_xml("a", 6, None, 1000.0)
            ok_missing = False
        except FileNotFoundError:
            ok_missing = True
    ok = ok_read and ok_pass and ok_dim and ok_cut and ok_missing
    return ok, f"解析 {ok_read}、放行 {ok_pass}、维度不符报错 {ok_dim}、超 cutOff 报错 {ok_cut}、缺 XML 报错 {ok_missing}"


def test_event_direction() -> tuple[bool, str]:
    """整窗速率恒负：argmax 仍给得出时间，但扩张支持度必须为 0。"""
    grid = np.linspace(0, 4000, 9)
    rng = np.random.default_rng(7)
    r = -rng.uniform(0.2, 0.8, size=(30, len(grid)))
    y = np.cumsum(-r, axis=1) * 0.5
    feats = segment_features(y, r, grid, np.ones(len(grid), bool))
    p_exp, _ = directional_probability("max_expansion_time", feats)
    p_con, _ = directional_probability("max_contraction_rate", feats)
    p_lvl, _ = directional_probability("level", feats)
    ok = p_exp == 0.0 and p_con == 1.0 and not np.isfinite(p_lvl)
    return ok, (f"扩张支持度 {p_exp:.3f}、收缩支持度 {p_con:.3f}、"
                f"水平不给概率 {not np.isfinite(p_lvl)}")


def test_label_alignment() -> tuple[bool, str]:
    """前一个 bin 缺席的类别不参与匹配，共有类别的编号必须保持一致。"""
    prev = {"A": 1, "B": 1, "C": 2}
    cur = {"A": 2, "B": 2, "C": 1, "D": 1}      # 编号整体互换，且多出一个 D
    out = align_cluster_labels_across_bins(prev, cur)
    ok = out["A"] == out["B"] == 1 and out["C"] == 2 and out["D"] == 2
    return ok, f"对齐结果 {out}"


def test_ward_guard() -> tuple[bool, str]:
    cases = {("euclidean", "ward"): True, ("correlation", "ward"): False,
             ("correlation", "average"): True, ("cityblock", "complete"): True}
    bad = [k for k, want in cases.items() if linkage_is_valid(*k) != want]
    return not bad, f"不合预期的组合 {bad}" if bad else "全部组合判定正确"


def test_end_to_end() -> tuple[bool, str]:
    """合成三个类别的 mini log+xml，跑通特征提取。"""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "extract_features", Path(__file__).resolve().parent / "1-extract_features.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "input").mkdir()
        (root / "output").mkdir()
        for name, level, trend in [("Alpha", 9.0, 1.2), ("Beta", 9.1, 1.1),
                                   ("Gamma", 10.5, -0.9)]:
            _write_fake_run(root / "input", name, level, trend)
        (root / "config.yaml").write_text(
            "io:\n  input_dir: input\n  output_dir: output\n"
            "time:\n  time_periods:\n    - [0, 4000]\n    - [4000, 8000]\n"
            "  bin_width: 4000\n  analysis_max: 8000\n  grid_step: 400\n"
            "  smooth_bandwidth: 2000\n"
            "posterior:\n  max_draws: 40\n  n_resample: 20\n"
            "select:\n  min_coverage: 0.5\n", encoding="utf-8")
        rc = mod.run(str(root / "config.yaml"), str(root))
        tab = pd.read_csv(root / "output" / "⭐1-Feature-Summary.tsv", sep="\t")
    n_cat = tab["Category"].nunique()
    n_feat = tab["Feature"].nunique()
    has_stmt = tab.loc[tab["Feature"] == "Time of Maximum Contraction",
                       "Probability Statement"].notna().all()
    ok = rc == 0 and n_cat == 3 and n_feat == 8 and has_stmt
    return ok, f"退出码 {rc}，{n_cat} 个类别 × {n_feat} 个特征，方向说明齐全 {has_stmt}"


TESTS = [
    ("Test 1 Skygrid 分段边界归属", test_skygrid_boundary),
    ("Test 2 XML 与日志一致性守卫", test_xml_guard),
    ("Test 3 事件方向的后验支持度语义", test_event_direction),
    ("Test 4 跨时间 bin 的标签对齐", test_label_alignment),
    ("Test 5 Ward 连接要求欧氏距离", test_ward_guard),
    ("Test 6 mini BEAST 输入端到端", test_end_to_end),
]


def run(log_level: str = "INFO", log_file: str | None = None) -> int:
    setup_logging(log_level, log_file)
    failed = 0
    for name, fn in TESTS:
        try:
            ok, detail = fn()
        except Exception as exc:                        # noqa: BLE001
            ok, detail = False, f"抛出异常 {type(exc).__name__}: {exc}"
        setup_logging(log_level, log_file)              # 被测模块可能改过日志级别
        log.info("%s：%s —— %s", name, "PASS" if ok else "FAIL", detail)
        failed += 0 if ok else 1
    log.info("自检完成：%d/%d 通过", len(TESTS) - failed, len(TESTS))
    return 0 if failed == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="通用 Ne 轨迹工具的最小自检")
    p.add_argument("--log-level", default="INFO")
    p.add_argument("--log-file", default=None)
    return p


def main(argv=None) -> int:
    return run(**vars(build_parser().parse_args(argv)))


if __name__ == "__main__":
    raise SystemExit(main())

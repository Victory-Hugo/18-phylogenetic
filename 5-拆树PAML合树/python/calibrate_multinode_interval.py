#!/usr/bin/env python3
"""多节点时间校准 + 基于子树离散度的 Monte-Carlo 年代区间。

功能概述
========
1. 多节点校准（不重算 PAML）：
   - 读入超度量树（subst/site 单位），按给定单系（如 N、M）与 root 的目标
     TMRCA，对树做"分支系独立缩放"，使每个校准节点精确命中目标年龄，同时
     保持超度量与正枝长。
   - 之所以不用全局单调的"深度->年代"曲线：当前树中 M 节点比 N 节点更深，
     而目标年龄要求 N 老于 M，顺序相反，全局单调曲线会产生负枝长。N、M 为
     互不嵌套的单系，分别按各自 crown 目标年龄缩放 clade 内部、骨架按 root
     目标年龄缩放，N/M 茎枝自动吸收差值即可。

2. 年代区间（不重算 PAML，只用已有 output/）：
   - 每个 PAML 子树在合并阶段都留有相对骨架的尺度离散度
     （output/merge/subtree_scale_report.tsv 的 residual_mad / median_log_ratio）。
   - 以 sigma_s = residual_mad_s * 1.4826 作为该子树枝长的对数空间乘性散布，
     做 Monte-Carlo 自助：每个重复对枝长抽乘性扰动 -> 重算节点高度 ->
     重做同一套分区校准（N/M/root 仍固定）-> 记录每个非锚点节点年龄，
     最终取百分位得到置信区间。

本模块为"双模块"：既可作为库被 import 调用 run(...)，也可作为命令行脚本运行。
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# 复用既有公共工具，避免重复造轮子
from phylo_merge_common import validate_branch_lengths_complete
from phylo_split_common import (
    PipelineError,
    assign_node_ids,
    clone_tree,
    decode_json_list,
    load_json,
    load_table,
    read_newick_tree,
    setup_logger,
    write_table,
    write_tree,
)

# MAD -> 稳健标准差 的换算常数（正态分布下 SD = 1.4826 * MAD）
MAD_TO_SD = 1.4826

NODE_AGE_CI_COLUMNS = [
    "node_id",
    "clade_label",
    "zone",
    "is_anchor",
    "n_tips",
    "point_age_years",
    "age_median_years",
    "age_low_years",
    "age_high_years",
    "ci_width_years",
]

SUMMARY_COLUMNS = [
    "item",
    "value",
]


# ---------------------------------------------------------------------------
# 数据载入与单系定位
# ---------------------------------------------------------------------------
def load_id_haplogroup(id_hap_tsv: Path) -> Dict[str, str]:
    """读取 ID -> 单倍群 映射表（两列、制表符分隔、无表头）。"""
    mapping: Dict[str, str] = {}
    with id_hap_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0]:
                mapping[parts[0]] = parts[1]
    if not mapping:
        raise PipelineError(f"ID-单倍群映射为空: {id_hap_tsv}")
    return mapping


def resolve_lineage(haplogroup: str, haplogroups: Dict[str, dict]) -> List[str]:
    """返回某单倍群从 mt-MRCA 到自身的谱系列表；解析不到则返回空列表。"""
    record = haplogroups.get(haplogroup)
    if record is None:
        # 去掉变异后缀（如 "HV_A73G!" -> "HV"）再尝试
        record = haplogroups.get(haplogroup.split("_")[0])
    if record is None:
        return []
    return list(record.get("lineage", []))


def collect_clade_tip_names(
    tip_to_lineage: Dict[str, List[str]],
    target_haplogroup: str,
) -> List[str]:
    """收集谱系中包含 target_haplogroup 的全部 tip 名（即该单系成员）。"""
    return [
        tip_name
        for tip_name, lineage in tip_to_lineage.items()
        if target_haplogroup in lineage
    ]


# ---------------------------------------------------------------------------
# 树结构索引（用整数下标表达，便于 numpy 向量化）
# ---------------------------------------------------------------------------
class TreeIndex:
    """把 Biopython 树展开成下标数组，供向量化高度/年龄计算使用。"""

    def __init__(self, tree) -> None:
        clades = list(tree.find_clades(order="preorder"))
        self.clades = clades
        self.n_nodes = len(clades)
        self.clade_to_idx = {clade: i for i, clade in enumerate(clades)}

        self.parent_idx = np.full(self.n_nodes, -1, dtype=np.int64)
        self.children: List[List[int]] = [[] for _ in range(self.n_nodes)]
        self.is_tip = np.zeros(self.n_nodes, dtype=bool)
        self.tip_name: List[Optional[str]] = [None] * self.n_nodes
        self.base_branch = np.zeros(self.n_nodes, dtype=np.float64)

        for clade in clades:
            idx = self.clade_to_idx[clade]
            bl = 0.0 if clade.branch_length is None else float(clade.branch_length)
            self.base_branch[idx] = bl
            if clade.is_terminal():
                self.is_tip[idx] = True
                self.tip_name[idx] = str(clade.name)
            for child in clade.clades:
                cidx = self.clade_to_idx[child]
                self.parent_idx[cidx] = idx
                self.children[idx].append(cidx)

        # 后序遍历顺序（叶在前、根在后），用于自底向上聚合
        self.postorder = self._build_postorder()
        # 每个节点的后代 tip 数（带权求均路长用）
        self.n_tips = self._count_tips()

    def _build_postorder(self) -> np.ndarray:
        order: List[int] = []
        # 迭代式后序，避免深树递归溢出
        stack: List[Tuple[int, bool]] = [(self._root_idx(), False)]
        while stack:
            idx, processed = stack.pop()
            if processed:
                order.append(idx)
                continue
            stack.append((idx, True))
            for cidx in self.children[idx]:
                stack.append((cidx, False))
        return np.asarray(order, dtype=np.int64)

    def _root_idx(self) -> int:
        roots = np.where(self.parent_idx < 0)[0]
        if len(roots) != 1:
            raise PipelineError(f"树根数量异常: {len(roots)}")
        return int(roots[0])

    @property
    def root_idx(self) -> int:
        return int(np.where(self.parent_idx < 0)[0][0])

    def _count_tips(self) -> np.ndarray:
        n_tips = np.zeros(self.n_nodes, dtype=np.float64)
        for idx in self.postorder:
            if self.is_tip[idx]:
                n_tips[idx] = 1.0
            else:
                n_tips[idx] = sum(n_tips[c] for c in self.children[idx])
        return n_tips

    def descendants(self, idx: int) -> List[int]:
        """返回包含自身在内的所有后代下标。"""
        out: List[int] = []
        stack = [idx]
        while stack:
            cur = stack.pop()
            out.append(cur)
            stack.extend(self.children[cur])
        return out


def compute_heights(tree_index: TreeIndex, branch: np.ndarray) -> np.ndarray:
    """按"均路长"计算每个节点到其后代 tip 的平均距离（支持二维 [节点, 重复]）。

    branch 形状可为 (n_nodes,) 或 (n_nodes, R)。返回同形状的高度数组。
    对超度量树（无扰动）该均值即精确高度；扰动后取均值以得到唯一高度。
    """
    two_d = branch.ndim == 2
    width = branch.shape[1] if two_d else 1
    shape = (tree_index.n_nodes, width)
    path_sum = np.zeros(shape, dtype=np.float64)  # 后代 tip 路径长度之和
    br = branch.reshape(tree_index.n_nodes, width)
    n_tips = tree_index.n_tips
    for idx in tree_index.postorder:
        if tree_index.is_tip[idx]:
            continue
        acc = np.zeros(width, dtype=np.float64)
        for cidx in tree_index.children[idx]:
            acc += path_sum[cidx] + n_tips[cidx] * br[cidx]
        path_sum[idx] = acc
    heights = path_sum / n_tips.reshape(-1, 1)
    return heights if two_d else heights.reshape(-1)


# ---------------------------------------------------------------------------
# 分区与年龄映射
# ---------------------------------------------------------------------------
def build_zone_array(
    tree_index: TreeIndex,
    zones: Sequence[Tuple[int, float, str]],
) -> np.ndarray:
    """按"最近校准祖先"为每个节点分区，返回 zone_of（指向 zones 列表的下标）。

    zones 末项约定为 root 区（覆盖全树）。其余项需按 clade 由小到大（即由深到浅）
    排在前面：分配时"已分配则跳过"，于是每个节点归入包含它的最小校准 clade，
    即其最近校准祖先 —— 这天然支持嵌套约束（如 root ⊃ L3 ⊃ {N, M}）。
    """
    zone_of = np.full(tree_index.n_nodes, -1, dtype=np.int16)
    for zi, (crown_idx, _age, _label) in enumerate(zones):
        for idx in tree_index.descendants(crown_idx):
            if zone_of[idx] == -1:
                zone_of[idx] = zi
    # root 区 crown 即根，其 descendants 覆盖全树，理论上不会再有 -1
    if (zone_of < 0).any():
        zone_of[zone_of < 0] = len(zones) - 1
    return zone_of


def heights_to_ages(
    heights: np.ndarray,
    zone_of: np.ndarray,
    zones: Sequence[Tuple[int, float, str]],
) -> np.ndarray:
    """把节点高度按所属分区线性缩放成年龄（支持一维/二维）。

    每个区的缩放因子 = 该区校准 crown 的目标年龄 / 其当前高度；crown 自身经本区
    缩放后精确命中目标年龄。
    """
    two_d = heights.ndim == 2
    h = heights.reshape(heights.shape[0], -1)
    factor = np.empty_like(h)
    for zi, (crown_idx, age, _label) in enumerate(zones):
        f = age / h[crown_idx]            # 形状 (R,)
        factor[zone_of == zi] = f         # 按行广播到该区所有节点
    ages = h * factor
    return ages if two_d else ages.reshape(-1)


def ages_to_branch_tree(tree, tree_index: TreeIndex, ages: np.ndarray, min_branch: float):
    """用节点年龄重算枝长（bl = age(parent) - age(child)）并返回克隆树。"""
    calibrated = clone_tree(tree)
    calib_clades = list(calibrated.find_clades(order="preorder"))
    # 克隆树与原树 preorder 顺序一致，可直接对位
    n_negative = 0
    for i, clade in enumerate(calib_clades):
        pidx = tree_index.parent_idx[i]
        if pidx < 0:
            clade.branch_length = None  # 根无枝长
            continue
        bl = float(ages[pidx] - ages[i])
        if bl < 0:
            n_negative += 1
            bl = min_branch
        clade.branch_length = bl
    return calibrated, n_negative


# ---------------------------------------------------------------------------
# 子树离散度（sigma）归属
# ---------------------------------------------------------------------------
def load_subtree_sigma(
    subtree_scale_report: Path,
    use_systematic_offset: bool,
) -> Tuple[Dict[str, float], Dict[str, float], float]:
    """读取每个 target 子树的 residual_mad/median_log_ratio，换算 sigma/offset。

    返回 (sigma_map, offset_map, sigma_backbone_default)。
    """
    rows = load_table(subtree_scale_report)
    sigma_map: Dict[str, float] = {}
    offset_map: Dict[str, float] = {}
    for row in rows:
        sid = row.get("target_subtree_id", "")
        if not sid:
            continue
        try:
            residual_mad = float(row.get("residual_mad", "0") or 0.0)
        except ValueError:
            residual_mad = 0.0
        sigma_map[sid] = residual_mad * MAD_TO_SD
        if use_systematic_offset:
            try:
                offset_map[sid] = float(row.get("median_log_ratio", "0") or 0.0)
            except ValueError:
                offset_map[sid] = 0.0
        else:
            offset_map[sid] = 0.0
    if not sigma_map:
        raise PipelineError(f"未从 {subtree_scale_report} 解析到任何子树离散度")
    positive_sigmas = [s for s in sigma_map.values() if s > 0]
    sigma_backbone = float(np.median(positive_sigmas)) if positive_sigmas else 0.0
    return sigma_map, offset_map, sigma_backbone


def build_tip_to_subtree(target_subtree_summary: Path) -> Dict[str, str]:
    """从 target_subtree_summary.tsv 建立 tip 名 -> target 子树 id 的映射。"""
    rows = load_table(target_subtree_summary)
    tip_to_subtree: Dict[str, str] = {}
    for row in rows:
        sid = row.get("target_subtree_id", "")
        names_json = row.get("target_nonbackbone_tip_names", "")
        if not sid or not names_json:
            continue
        for tip_name in decode_json_list(names_json):
            tip_to_subtree[tip_name] = sid
    return tip_to_subtree


def assign_edge_sigma(
    tree_index: TreeIndex,
    tip_to_subtree: Dict[str, str],
    sigma_map: Dict[str, float],
    offset_map: Dict[str, float],
    sigma_backbone: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """为每条边赋 (sigma, offset)。

    节点归属：若其全部后代 tip 同属一个子树则归该子树，否则视为骨架。
    边的离散度采用该节点所属子树的 sigma；骨架边用 sigma_backbone。
    """
    MIXED = "::MIXED::"
    assignment: List[Optional[str]] = [None] * tree_index.n_nodes
    for idx in tree_index.postorder:
        if tree_index.is_tip[idx]:
            assignment[idx] = tip_to_subtree.get(tree_index.tip_name[idx], MIXED)
            continue
        child_assigns = {assignment[c] for c in tree_index.children[idx]}
        if len(child_assigns) == 1 and MIXED not in child_assigns:
            assignment[idx] = next(iter(child_assigns))
        else:
            assignment[idx] = MIXED

    sigma_edge = np.full(tree_index.n_nodes, sigma_backbone, dtype=np.float64)
    offset_edge = np.zeros(tree_index.n_nodes, dtype=np.float64)
    for idx in range(tree_index.n_nodes):
        sid = assignment[idx]
        if sid is not None and sid != MIXED and sid in sigma_map:
            sigma_edge[idx] = sigma_map[sid]
            offset_edge[idx] = offset_map.get(sid, 0.0)
    return sigma_edge, offset_edge


# ---------------------------------------------------------------------------
# 单系标签（取后代 tip 谱系的最长公共前缀末项）
# ---------------------------------------------------------------------------
def compute_clade_labels(
    tree_index: TreeIndex,
    tip_to_lineage: Dict[str, List[str]],
) -> List[str]:
    """每个节点的单倍群标签 = 其后代 tip 谱系的最长公共前缀的末项。"""
    prefixes: List[Optional[List[str]]] = [None] * tree_index.n_nodes
    for idx in tree_index.postorder:
        if tree_index.is_tip[idx]:
            prefixes[idx] = list(tip_to_lineage.get(tree_index.tip_name[idx], []))
            continue
        common: Optional[List[str]] = None
        for cidx in tree_index.children[idx]:
            cp = prefixes[cidx]
            if cp is None:
                continue
            if common is None:
                common = list(cp)
                continue
            limit = min(len(common), len(cp))
            k = 0
            while k < limit and common[k] == cp[k]:
                k += 1
            common = common[:k]
        prefixes[idx] = common if common is not None else []
    labels: List[str] = []
    for idx in range(tree_index.n_nodes):
        pre = prefixes[idx] or []
        labels.append(pre[-1] if pre else "")
    return labels


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def run(
    input_tree,
    id_hap_tsv,
    phylotree_json,
    merge_output_dir,
    split_output_dir,
    output_dir,
    calibrations,
    root_age,
    replicates=1000,
    seed=42,
    use_systematic_offset=False,
    ci_lower=2.5,
    ci_upper=97.5,
    output_tree_name="merged_multinode_tree_years.nwk",
    chunk_size=250,
    log_level="INFO",
):
    input_tree = Path(input_tree).resolve()
    id_hap_tsv = Path(id_hap_tsv).resolve()
    phylotree_json = Path(phylotree_json).resolve()
    merge_output_dir = Path(merge_output_dir).resolve()
    split_output_dir = Path(split_output_dir).resolve()
    output_dir = Path(output_dir).resolve()
    root_age = float(root_age)
    replicates = int(replicates)
    chunk_size = max(1, int(chunk_size))

    output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger("calibrate_multinode_interval", output_dir / "multinode_calibration.log", log_level)
    logger.info("开始多节点校准: 输入树=%s", input_tree)

    # --- 校准约束解析 ---
    if not calibrations:
        raise PipelineError("至少需要一个单系校准约束。")
    calib_by_hap = {str(c["haplogroup"]): float(c["age"]) for c in calibrations}
    logger.info("校准目标: %s, root=%.6g 年",
                ", ".join(f"{k}={v:.6g}" for k, v in calib_by_hap.items()), root_age)

    # --- 载入树与注释 ---
    tree = read_newick_tree(input_tree)
    validate_branch_lengths_complete(tree)
    id_to_hap = load_id_haplogroup(id_hap_tsv)
    haplogroups = load_json(phylotree_json)["haplogroups"]

    tip_clades = tree.get_terminals()
    tip_names_in_tree = [str(t.name) for t in tip_clades]
    tip_to_lineage: Dict[str, List[str]] = {}
    n_unresolved = 0
    for tip_name in tip_names_in_tree:
        hap = id_to_hap.get(tip_name)
        if hap is None:
            n_unresolved += 1
            tip_to_lineage[tip_name] = []
            continue
        lineage = resolve_lineage(hap, haplogroups)
        if not lineage:
            n_unresolved += 1
        tip_to_lineage[tip_name] = lineage
    logger.info("tip 总数=%d, 无法解析谱系=%d", len(tip_names_in_tree), n_unresolved)

    # --- 定位各单系 crown 并校验单系性 ---
    crown_clades = {}
    crown_member_counts = {}
    for hap in calib_by_hap:
        members = collect_clade_tip_names(tip_to_lineage, hap)
        if len(members) < 2:
            raise PipelineError(f"单系 {hap} 成员过少 (n={len(members)})，无法定位 crown。")
        member_set = set(members)
        member_clades = [t for t in tip_clades if str(t.name) in member_set]
        crown = tree.common_ancestor(member_clades)
        covered = len(crown.get_terminals())
        if covered != len(members):
            raise PipelineError(
                f"单系 {hap} 非单系: crown 覆盖 {covered} tip 但成员 {len(members)} 个。"
            )
        crown_clades[hap] = crown
        crown_member_counts[hap] = len(members)
        logger.info("单系 %s: 成员=%d, crown 覆盖=%d (单系校验通过)", hap, len(members), covered)

    # --- 建索引 ---
    tree_index = TreeIndex(tree)
    node_id_map, _, _ = assign_node_ids(tree)  # NODE_xxxxxx / TIP::name 标签
    node_ids = [node_id_map[clade] for clade in tree_index.clades]
    root_idx = tree_index.root_idx

    # --- 构造分区列表（按 clade 由小到大 = 由深到浅，root 区置于末尾）---
    # zones 元素: (crown_idx, 目标年龄, 标签)；支持嵌套约束 root ⊃ ... ⊃ 深层单系。
    hap_zones = []
    for hap, age in calib_by_hap.items():
        crown_idx = tree_index.clade_to_idx[crown_clades[hap]]
        hap_zones.append((crown_idx, age, hap))
    hap_zones.sort(key=lambda z: crown_member_counts[z[2]])  # 成员少者(更深)在前
    zones = hap_zones + [(root_idx, root_age, "root")]
    zone_of = build_zone_array(tree_index, zones)
    zone_label_of = {zi: z[2] for zi, z in enumerate(zones)}
    for zi, z in enumerate(zones):
        logger.info("分区 %s: 节点数=%d", z[2], int((zone_of == zi).sum()))

    # --- 点估计 ---
    base_heights = compute_heights(tree_index, tree_index.base_branch)
    point_ages = heights_to_ages(base_heights, zone_of, zones)
    calibrated_tree, n_neg_point = ages_to_branch_tree(tree, tree_index, point_ages, min_branch=0.0)
    if n_neg_point > 0:
        raise PipelineError(f"点估计出现 {n_neg_point} 条负枝长，校准约束相互冲突。")
    output_tree_path = output_dir / output_tree_name
    write_tree(calibrated_tree, output_tree_path)
    logger.info("已写出点估计年化树: %s", output_tree_path)

    # 校验点估计锚点命中
    anchor_idx = {tree_index.clade_to_idx[crown_clades[h]]: h for h in calib_by_hap}
    anchor_idx[root_idx] = "root"
    logger.info("锚点命中(点估计): %s",
                ", ".join(f"{lbl}={point_ages[i]:.2f}" for i, lbl in anchor_idx.items()))

    # --- 区间：MC 自助 ---
    sigma_map, offset_map, sigma_backbone = load_subtree_sigma(
        merge_output_dir / "subtree_scale_report.tsv", use_systematic_offset
    )
    tip_to_subtree = build_tip_to_subtree(split_output_dir / "target_subtree_summary.tsv")
    sigma_edge, offset_edge = assign_edge_sigma(
        tree_index, tip_to_subtree, sigma_map, offset_map, sigma_backbone
    )
    logger.info("子树离散度: 子树数=%d, sigma 中位=%.3g, sigma_backbone=%.3g",
                len(sigma_map), float(np.median(list(sigma_map.values()))), sigma_backbone)

    rng = np.random.default_rng(seed)
    ages_samples = np.empty((tree_index.n_nodes, replicates), dtype=np.float64)
    n_neg_mc = 0
    done = 0
    while done < replicates:
        c = min(chunk_size, replicates - done)
        # 乘性扰动：bl' = bl * exp(offset + sigma * Z)
        z = rng.standard_normal((tree_index.n_nodes, c))
        log_factor = offset_edge.reshape(-1, 1) + sigma_edge.reshape(-1, 1) * z
        branch = tree_index.base_branch.reshape(-1, 1) * np.exp(log_factor)
        heights = compute_heights(tree_index, branch)
        ages = heights_to_ages(heights, zone_of, zones)
        ages_samples[:, done:done + c] = ages
        done += c
    logger.info("MC 完成: 重复=%d", replicates)

    age_low = np.percentile(ages_samples, ci_lower, axis=1)
    age_median = np.percentile(ages_samples, 50.0, axis=1)
    age_high = np.percentile(ages_samples, ci_upper, axis=1)

    # 负枝长统计（基于中位年龄重算，作为健壮性诊断）
    for idx in range(tree_index.n_nodes):
        pidx = tree_index.parent_idx[idx]
        if pidx >= 0 and (age_median[pidx] - age_median[idx]) < 0:
            n_neg_mc += 1

    # --- 输出 node_age_ci.tsv ---
    clade_labels = compute_clade_labels(tree_index, tip_to_lineage)
    rows: List[Dict[str, object]] = []
    for idx in range(tree_index.n_nodes):
        if tree_index.is_tip[idx]:
            continue  # 仅报告内部节点（tip 年龄=0）
        is_anchor = idx in anchor_idx
        if is_anchor:
            low = med = high = ""
            width = ""
        else:
            low = f"{age_low[idx]:.6g}"
            med = f"{age_median[idx]:.6g}"
            high = f"{age_high[idx]:.6g}"
            width = f"{age_high[idx] - age_low[idx]:.6g}"
        rows.append({
            "node_id": node_ids[idx],
            "clade_label": anchor_idx[idx] if is_anchor else clade_labels[idx],
            "zone": zone_label_of[int(zone_of[idx])],
            "is_anchor": "true" if is_anchor else "false",
            "n_tips": str(int(tree_index.n_tips[idx])),
            "point_age_years": f"{point_ages[idx]:.6g}",
            "age_median_years": med,
            "age_low_years": low,
            "age_high_years": high,
            "ci_width_years": width,
        })
    write_table(rows, NODE_AGE_CI_COLUMNS, output_dir / "node_age_ci.tsv")
    logger.info("已写出节点年龄区间表: %s (%d 内部节点)", output_dir / "node_age_ci.tsv", len(rows))

    # --- 输出 summary ---
    summary_rows = [
        ("input_tree", input_tree.as_posix()),
        ("output_tree", output_tree_path.as_posix()),
        ("age_root_years", f"{root_age:.6g}"),
    ]
    for hap, age in calib_by_hap.items():
        cidx = tree_index.clade_to_idx[crown_clades[hap]]
        summary_rows.extend([
            (f"age_{hap}_years", f"{age:.6g}"),
            (f"{hap}_member_tips", str(crown_member_counts[hap])),
            (f"{hap}_input_height_subst", f"{base_heights[cidx]:.12g}"),
            (f"scale_factor_{hap}_years_per_subst", f"{age / base_heights[cidx]:.12g}"),
            (f"point_anchor_hit_{hap}", f"{point_ages[cidx]:.6g}"),
        ])
    summary_rows.extend([
        ("root_input_height_subst", f"{base_heights[root_idx]:.12g}"),
        ("scale_factor_backbone_years_per_subst", f"{root_age / base_heights[root_idx]:.12g}"),
        ("point_anchor_hit_root", f"{point_ages[root_idx]:.6g}"),
        ("point_negative_branches", str(n_neg_point)),
        ("mc_replicates", str(replicates)),
        ("mc_seed", str(seed)),
        ("mc_use_systematic_offset", str(use_systematic_offset)),
        ("mc_ci_lower_pct", f"{ci_lower:.6g}"),
        ("mc_ci_upper_pct", f"{ci_upper:.6g}"),
        ("mc_median_age_negative_branches", str(n_neg_mc)),
        ("subtree_sigma_median", f"{float(np.median(list(sigma_map.values()))):.6g}"),
        ("sigma_backbone", f"{sigma_backbone:.6g}"),
        ("tips_unresolved_haplogroup", str(n_unresolved)),
    ])
    write_table(
        [{"item": k, "value": v} for k, v in summary_rows],
        SUMMARY_COLUMNS,
        output_dir / "multinode_calibration_summary.tsv",
    )
    logger.info("已写出校准摘要: %s", output_dir / "multinode_calibration_summary.tsv")
    logger.info("多节点校准与区间估计完成。")
    return 0


# ---------------------------------------------------------------------------
# 命令行入口（双模块）
# ---------------------------------------------------------------------------
def parse_calibrations(value: str) -> List[dict]:
    """解析形如 'N=59000,M=49000' 的校准串。"""
    out: List[dict] = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise argparse.ArgumentTypeError(f"校准项格式应为 HAP=AGE: {item}")
        hap, age = item.split("=", 1)
        out.append({"haplogroup": hap.strip(), "age": float(age)})
    return out


def build_parser():
    parser = argparse.ArgumentParser(
        description="多节点时间校准 + 基于子树离散度的 Monte-Carlo 年代区间（不重算 PAML）。"
    )
    parser.add_argument("--input-tree", required=True, help="输入超度量树 (subst/site)。")
    parser.add_argument("--id-hap-tsv", required=True, help="ID->单倍群 映射表。")
    parser.add_argument("--phylotree-json", required=True, help="phylotree 结构 JSON。")
    parser.add_argument("--merge-output-dir", required=True, help="合并输出目录 (含 subtree_scale_report.tsv)。")
    parser.add_argument("--split-output-dir", required=True, help="拆树输出目录 (含 target_subtree_summary.tsv)。")
    parser.add_argument("--output-dir", required=True, help="输出目录。")
    parser.add_argument("--calibrations", required=True, type=parse_calibrations,
                        help="单系校准，如 'N=59000,M=49000'。")
    parser.add_argument("--root-age", required=True, type=float, help="root 目标年龄。")
    parser.add_argument("--replicates", default=1000, type=int, help="MC 重复数。")
    parser.add_argument("--seed", default=42, type=int, help="随机种子。")
    parser.add_argument("--use-systematic-offset", action="store_true",
                        help="是否叠加每子树系统尺度偏移 (median_log_ratio)。")
    parser.add_argument("--ci-lower", default=2.5, type=float, help="区间下分位 (%)。")
    parser.add_argument("--ci-upper", default=97.5, type=float, help="区间上分位 (%)。")
    parser.add_argument("--output-tree-name", default="merged_multinode_tree_years.nwk",
                        help="输出年化树文件名。")
    parser.add_argument("--chunk-size", default=250, type=int, help="MC 分块大小 (控制内存)。")
    parser.add_argument("--log-level", default="INFO", help="日志级别。")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            input_tree=args.input_tree,
            id_hap_tsv=args.id_hap_tsv,
            phylotree_json=args.phylotree_json,
            merge_output_dir=args.merge_output_dir,
            split_output_dir=args.split_output_dir,
            output_dir=args.output_dir,
            calibrations=args.calibrations,
            root_age=args.root_age,
            replicates=args.replicates,
            seed=args.seed,
            use_systematic_offset=args.use_systematic_offset,
            ci_lower=args.ci_lower,
            ci_upper=args.ci_upper,
            output_tree_name=args.output_tree_name,
            chunk_size=args.chunk_size,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

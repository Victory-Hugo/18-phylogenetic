#!/usr/bin/env python3
"""2-1 节点年龄：把时间树的内部节点摊到各分组类别上。

节点归属定义
  Within-group subtree (all samples)
      把树剪枝到某个类别自己的采样个体，其内部节点即该类别的合并事件。剪枝不改变节点高度，
      故无需真的构造剪枝树：某内部节点属于该类别的诱导子树 <=> 它至少有两个子节点的子树各含
      >=1 个该类别的 tip。每个类别恰好得到 n-1 个节点。

  Within-group subtree (rarefied)
      同上，但先把同一维度内的所有类别随机下采样到相同的样本量再取诱导子树，重复若干次。
      这一步是组间峰形可比的前提：诱导子树的形状本身依赖样本量，样本多的类别天然包含更多
      极近期分叉，不做等样本量抽稀就直接比较峰形并不公平。

产物
  output/result/2-node_ages/⭐节点年龄长表.tsv.gz     全样本诱导子树
  temp/2-node_ages/rarefied-node-ages.tsv.gz         抽稀重复的节点年龄（计算成本高，保留）

CLI
  python python/2-1-node_ages.py --conf ... --sample-groups ... --output-result ... --temp ...
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from collections import defaultdict
from pathlib import Path

import dendropy
import numpy as np
import yaml

log = logging.getLogger(__name__)

OUTPUT_NAME = "⭐节点年龄长表.tsv.gz"
RAREFIED_NAME = "rarefied-node-ages.tsv.gz"
FIELDNAMES = ["Grouping", "Group", "Node identifier", "Node age (kya)"]
RAREFIED_FIELDNAMES = ["Grouping", "Group", "Replicate", "Node age (kya)"]


class TreeArrays:
    """把时间树摊平成数组，后续所有遍历都在纯 Python/NumPy 上做，不再触碰 dendropy。"""

    def __init__(self, tree_path: Path):
        # dendropy.Tree.get 只取文件里的第一棵树，后面的会被静默丢掉。本流程按“一棵固定时间树”
        # 设计（下游区间只反映抽稀，不含后验定年不确定性），故 posterior tree set 必须显式报错。
        n_trees = tree_path.read_text().count(";")
        if n_trees > 1:
            raise ValueError(
                f"{tree_path} 含 {n_trees} 棵树；本流程只支持单棵时间树。若要把后验定年"
                f"不确定性纳入区间，需要先扩展本步骤按树循环，并相应改写 3-1 的区间定义。")
        tree = dendropy.Tree.get(path=str(tree_path), schema="newick",
                                 preserve_underscores=True, rooting="force-rooted")
        tree.calc_node_ages(ultrametricity_precision=1e-3)
        self.children: list[list[int]] = []
        self.age_kya: list[float] = []
        self.tip_label: dict[int, str] = {}
        index: dict[int, int] = {}
        for node in tree.postorder_node_iter():
            node_id = len(self.children)
            index[id(node)] = node_id
            self.children.append([index[id(child)] for child in node.child_node_iter()])
            self.age_kya.append(node.age / 1000.0)
            if node.is_leaf():
                self.tip_label[node_id] = node.taxon.label
        self.root_age_kya = self.age_kya[-1]
        self.n_nodes = len(self.children)
        self.internal = [i for i in range(self.n_nodes) if self.children[i]]


def _tip_category_map(groups_path: Path,
                      tip_labels: dict[int, str]) -> dict[str, dict[int, str]]:
    """读 sample-groups.tsv → {维度: {tip 节点编号: 类别}}。"""
    sample_to_tips: dict[str, list[int]] = defaultdict(list)
    for node_id, label in tip_labels.items():
        sample_to_tips[label].append(node_id)

    mapping: dict[str, dict[int, str]] = defaultdict(dict)
    with open(groups_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for node_id in sample_to_tips.get(row["Sample"], ()):
                mapping[row["Grouping"]][node_id] = row["Group"]
    return mapping


def _within_group_nodes(tree: TreeArrays, members: set[int]) -> list[float]:
    """某个 tip 集合的诱导子树内部节点年龄；至少两个子分支含成员才算一次合并事件。"""
    present = [False] * tree.n_nodes
    ages = []
    for node_id in range(tree.n_nodes):
        kids = tree.children[node_id]
        if not kids:
            present[node_id] = node_id in members
            continue
        hits = sum(1 for child in kids if present[child])
        present[node_id] = hits > 0
        if hits >= 2:
            ages.append(tree.age_kya[node_id])
    return ages


def _within_group_nodes_replicated(tree: TreeArrays, member_matrix: np.ndarray,
                                   members_index: np.ndarray) -> list[np.ndarray]:
    """一次遍历同时处理 B 个抽稀重复。member_matrix 形状 (n_tips, B)，行序对应 members_index。"""
    n_replicates = member_matrix.shape[1]
    present = np.zeros((tree.n_nodes, n_replicates), dtype=bool)
    present[members_index] = member_matrix
    collected: list[list[float]] = [[] for _ in range(n_replicates)]
    for node_id in tree.internal:
        kids = tree.children[node_id]
        hits = present[kids].sum(axis=0)
        present[node_id] = hits > 0
        emitting = np.flatnonzero(hits >= 2)
        if emitting.size:
            age = tree.age_kya[node_id]
            for replicate in emitting:
                collected[replicate].append(age)
    return [np.asarray(x) for x in collected]


def run(conf: str, sample_groups: str, output_result: str, temp: str) -> int:
    settings = yaml.safe_load(Path(conf).read_text(encoding="utf-8"))
    replicates = int(settings["analysis"]["rarefaction_replicates"])
    fraction = float(settings["analysis"].get("rarefaction_fraction", 1.0))
    fixed_size = int(settings["analysis"].get("rarefaction_size", 0) or 0)
    seed = int(settings["analysis"]["rarefaction_seed"])

    result_path, temp_path = Path(output_result), Path(temp)
    for directory in (result_path, temp_path):
        directory.mkdir(parents=True, exist_ok=True)

    arrays = TreeArrays(Path(settings["paths"]["tree"]))
    log.info("树读入完成：%d tips，%d 个内部节点，根深 %.1f kya",
             len(arrays.tip_label), len(arrays.internal), arrays.root_age_kya)

    mapping = _tip_category_map(Path(sample_groups), arrays.tip_label)
    rng = np.random.default_rng(seed)
    rarefied_file = temp_path / RAREFIED_NAME
    total_rows = 0

    with gzip.open(result_path / OUTPUT_NAME, "wt", newline="") as fh, \
         gzip.open(rarefied_file, "wt", newline="") as fh_rare:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(FIELDNAMES)
        writer_rare = csv.writer(fh_rare, delimiter="\t")
        writer_rare.writerow(RAREFIED_FIELDNAMES)

        for grouping, tip_category in mapping.items():
            categories = sorted(set(tip_category.values()))
            members = {c: {n for n, v in tip_category.items() if v == c} for c in categories}

            for category in categories:
                ages = _within_group_nodes(arrays, members[category])
                expected = len(members[category]) - 1
                if len(ages) != expected:
                    raise AssertionError(
                        f"{grouping}/{category}: 诱导子树节点数 {len(ages)} != n-1 = {expected}")
                for node_index, age in enumerate(ages):
                    writer.writerow([grouping, category, node_index, f"{age:.6f}"])
                total_rows += len(ages)
            log.info("  %-30s 全样本诱导子树完成（各类别节点数均等于 n-1）", grouping)

            if replicates <= 0:
                continue
            # 目标样本量若正好等于最小类别的样本量，该类别的每次重复都是同一批人：
            # 区间带宽恒为 0，曲线也完全没有被跨重复的中位数平滑，近期尖峰会以全高保留，
            # 而其余类别的尖峰被平均掉。乘一个小于 1 的系数，使每个类别都真正被抽稀。
            smallest = min(len(members[c]) for c in categories)
            if fixed_size > 0:
                # 绝对目标：各维度共用同一个样本量，维度之间也因此可比
                target = min(fixed_size, smallest)
                if fixed_size > smallest:
                    log.warning("  %-30s 目标 %d 超过最小类别 %d 例，已降到 %d",
                                grouping, fixed_size, smallest, target)
            else:
                target = max(2, int(fraction * smallest))
            for category in categories:
                pool = np.array(sorted(members[category]))
                matrix = np.zeros((pool.size, replicates), dtype=bool)
                for replicate in range(replicates):
                    matrix[rng.choice(pool.size, size=target, replace=False), replicate] = True
                for replicate, ages in enumerate(
                        _within_group_nodes_replicated(arrays, matrix, pool)):
                    if ages.size != target - 1:
                        raise AssertionError(
                            f"{grouping}/{category} 重复 {replicate}: 抽稀诱导子树节点数 "
                            f"{ages.size} != {target - 1}")
                    for age in ages:
                        writer_rare.writerow([grouping, category, replicate, f"{age:.6f}"])
            log.info("  %-30s 抽稀完成：每类降到 %d 例 × %d 次重复"
                     "（最小类别 %d 例，抽样比 %.0f%%）",
                     grouping, target, replicates, smallest, target / smallest * 100)

    if replicates <= 0:
        rarefied_file.unlink()
        log.info("抽稀已关闭，仅输出全样本诱导子树；组间峰形不可直接比较")
    else:
        log.info("抽稀节点年龄 → %s", rarefied_file)
    log.info("节点年龄长表 %d 行 → %s", total_rows, result_path / OUTPUT_NAME)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="按分组类别计算诱导子树的节点年龄")
    parser.add_argument("--conf", required=True)
    parser.add_argument("--sample-groups", required=True)
    parser.add_argument("--output-result", required=True)
    parser.add_argument("--temp", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    return run(**vars(build_parser().parse_args(argv)))


if __name__ == "__main__":
    raise SystemExit(main())

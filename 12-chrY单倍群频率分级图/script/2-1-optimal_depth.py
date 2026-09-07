#!/usr/bin/env python3
"""用交叉验证为 Y 染色体单倍群选择最优系统发育分辨率。

同时给出两种答案：

全局统一深度
    所有主干使用同一个 Level。对 Level 1..max_level 做重复分层交叉验证，
    以单倍群类别预测群体，比较 predictive log-loss，再用 one-standard-error
    规则在与最优无显著差异的候选中取最浅的 Level。

各分支自适应深度
    把 CART 的 cost-complexity pruning 搬到单倍群系统树上：对每个内部节点
    比较「停在此处」与「继续拆分」的目标函数

        J(T) = Σ_i log P(Class_i | g_T(H_i)) − λ·K_T

    用自底向上动态规划求给定 λ 的最优 tree cut，λ 同样由重复分层交叉验证
    加 one-standard-error 规则选出。因此 C 可以停在 Level 1，而 O2b1a 可以
    一直拆到 Level 5。

两种方案都遵循项目既有的 ISOGG 解析规则：节点被拆分时，恰好终止于该节点、
无法进入任何子节点的样本单独归入带星号的类别（如 O2*），不丢弃任何样本。

用法::

    python 2-1-optimal_depth.py --input-path input/example.tsv \
        --output-result output/result/2-optimal_depth \
        --output-figure output/figure/2-optimal_depth

也可作为模块导入后调用 :func:`run`，或用 --self-test 运行内置合成数据检验。
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import os
from typing import Dict, List, Sequence, Tuple

import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

sys.dont_write_bytecode = True

log = logging.getLogger(__name__)

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPEC = importlib.util.spec_from_file_location(
    "haplogroup_level_frequency",
    os.path.join(_HERE, "1-1-haplogroup_level_frequency.py"))
_BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASE)

tokenize = _BASE.tokenize
assign_level = _BASE.assign_level
UNRESOLVED_MARK = _BASE.UNRESOLVED_MARK
_category_sort_key = _BASE._category_sort_key

Key = Tuple[str, ...]


# ---------------------------------------------------------------- 概率与打分

def smoothed_log_prob(counts: np.ndarray, alpha: float) -> np.ndarray:
    """一个类别内群体分布的平滑对数概率（Dirichlet/Laplace 平滑）。"""
    total = counts.sum()
    return np.log((counts + alpha) / (total + alpha * counts.size))


def category_log_likelihood(counts: np.ndarray, alpha: float) -> float:
    """把一组样本当作单一类别时的训练对数似然。"""
    if counts.sum() == 0:
        return 0.0
    return float((counts * smoothed_log_prob(counts, alpha)).sum())


# ---------------------------------------------------------------- 系统树构建

class HaplogroupTree:
    """由 ISOGG 名称本身构建的单倍群前缀树，不使用任何群体标签。"""

    def __init__(self, haplogroups: Sequence[str]):
        self.paths: Dict[str, Key] = {h: tokenize(h) for h in set(haplogroups)}
        self.children: Dict[Key, List[Key]] = {}
        self.parent: Dict[Key, Key] = {}
        nodes = set()
        for tokens in self.paths.values():
            for depth in range(1, len(tokens) + 1):
                key = tokens[:depth]
                nodes.add(key)
                if depth > 1:
                    self.parent[key] = tokens[:depth - 1]
        self.nodes: List[Key] = sorted(nodes, key=lambda k: (len(k), k))
        for key in self.nodes:
            self.children.setdefault(key, [])
        for key, par in self.parent.items():
            self.children[par].append(key)
        for key in self.nodes:
            self.children[key].sort(key=_category_sort_key)
        self.roots: List[Key] = [k for k in self.nodes if len(k) == 1]
        # 自底向上的处理顺序
        self.post_order: List[Key] = sorted(
            self.nodes, key=lambda k: -len(k))

    @staticmethod
    def name(key: Key) -> str:
        return "".join(key)

    def count_matrix(self, tokens_per_sample: Sequence[Key],
                     class_index: np.ndarray, n_classes: int
                     ) -> Tuple[Dict[Key, np.ndarray], Dict[Key, np.ndarray]]:
        """返回每个节点的 terminal 计数与 subtree 计数矩阵。"""
        terminal = {key: np.zeros(n_classes) for key in self.nodes}
        for tokens, cls in zip(tokens_per_sample, class_index):
            terminal[tokens][cls] += 1
        subtree = {key: terminal[key].copy() for key in self.nodes}
        for key in self.post_order:
            par = self.parent.get(key)
            if par is not None:
                subtree[par] += subtree[key]
        return terminal, subtree


# ------------------------------------------------------- 自适应深度动态规划

class TreeCut:
    """一次 tree cut 的结果：每个节点停止还是继续拆分。"""

    def __init__(self, tree: HaplogroupTree, decision: Dict[Key, str]):
        self.tree = tree
        self.decision = decision

    def categories(self) -> List[Tuple[str, Key, bool]]:
        """返回 (类别名, 所在节点, 是否为祖先残留类别) 的 frontier 列表。"""
        out: List[Tuple[str, Key, bool]] = []
        stack = list(reversed(self.tree.roots))
        while stack:
            key = stack.pop()
            if self.decision[key] == "stop":
                out.append((HaplogroupTree.name(key), key, False))
            else:
                if self.decision.get(("*",) + key) == "leaf":
                    out.append((HaplogroupTree.name(key) + UNRESOLVED_MARK,
                                key, True))
                stack.extend(reversed(self.tree.children[key]))
        return out

    def assign(self, tokens: Key) -> str:
        """把一个样本的单倍群路径映射到 frontier 上的类别名。"""
        key = tokens[:1]
        if key not in self.decision:
            return ""
        while True:
            if self.decision[key] == "stop":
                return HaplogroupTree.name(key)
            depth = len(key)
            star = HaplogroupTree.name(key) + UNRESOLVED_MARK
            if len(tokens) == depth:
                return star
            child = tokens[:depth + 1]
            if child not in self.decision:
                # 训练时没见过的更深分支，退回该节点的祖先残留类别
                return star if self.decision.get(("*",) + key) == "leaf" else ""
            key = child


def fit_tree_cut(tree: HaplogroupTree, terminal: Dict[Key, np.ndarray],
                 subtree: Dict[Key, np.ndarray], lam: float,
                 alpha: float) -> Tuple[TreeCut, float]:
    """给定 λ，自底向上动态规划求最优 tree cut。"""
    stop_score: Dict[Key, float] = {}
    split_score: Dict[Key, float] = {}
    best: Dict[Key, float] = {}
    decision: Dict[Key, str] = {}

    for key in tree.post_order:
        s_stop = category_log_likelihood(subtree[key], alpha) - lam
        stop_score[key] = s_stop
        kids = tree.children[key]
        if not kids:
            best[key] = s_stop
            decision[key] = "stop"
            split_score[key] = float("-inf")
            continue
        s_split = sum(best[c] for c in kids)
        if terminal[key].sum() > 0:
            s_split += category_log_likelihood(terminal[key], alpha) - lam
        split_score[key] = s_split
        if s_split > s_stop:
            best[key] = s_split
            decision[key] = "split"
            if terminal[key].sum() > 0:
                decision[("*",) + key] = "leaf"
        else:
            best[key] = s_stop
            decision[key] = "stop"

    total = sum(best[r] for r in tree.roots)
    cut = TreeCut(tree, decision)
    cut.stop_score = stop_score
    cut.split_score = split_score
    return cut, total


# ------------------------------------------------------------ 交叉验证框架

def _stratify_labels(classes: np.ndarray, n_splits: int) -> np.ndarray:
    """样本量不足以分层的群体合并为一个稀有层，避免丢弃样本。"""
    values, counts = np.unique(classes, return_counts=True)
    rare = set(values[counts < n_splits])
    if rare:
        log.warning("%d 个群体样本量 < %d 折，交叉验证时合并为一个稀有层",
                    len(rare), n_splits)
    return np.array(["__rare__" if c in rare else c for c in classes])


def _resolve_folds(classes: np.ndarray, cv_folds: int) -> int:
    """必要时下调折数，保证分层交叉验证可执行。"""
    labels = _stratify_labels(classes, cv_folds)
    smallest = np.unique(labels, return_counts=True)[1].min()
    folds = int(min(cv_folds, smallest))
    if folds < cv_folds:
        log.warning("最小分层样本量为 %d，折数由 %d 下调为 %d",
                    smallest, cv_folds, folds)
    if folds < 2:
        raise ValueError("样本量不足以完成交叉验证")
    return folds


def _log_loss(train_counts: Dict[str, np.ndarray], prior: np.ndarray,
              val_categories: Sequence[str], val_class: np.ndarray,
              alpha: float) -> float:
    """验证集上的 multiclass log-loss，未见类别退回训练集群体先验。"""
    cache = {c: smoothed_log_prob(v, alpha) for c, v in train_counts.items()}
    prior_lp = smoothed_log_prob(prior, alpha)
    total = 0.0
    for cat, cls in zip(val_categories, val_class):
        lp = cache.get(cat)
        total -= (prior_lp if lp is None else lp)[cls]
    return total / len(val_class)


# ------------------------------------------------------------ 全局统一深度

def evaluate_global_levels(tokens_per_sample: Sequence[Key],
                           class_index: np.ndarray, n_classes: int,
                           max_level: int, alpha: float, cv_folds: int,
                           cv_repeats: int, random_seed: int) -> pd.DataFrame:
    """对每个候选 Level 做重复分层交叉验证，返回 log-loss 曲线。"""
    folds = _resolve_folds(class_index, cv_folds)
    strat = _stratify_labels(class_index, folds)
    level_cats = {
        lv: [assign_level(t, lv) for t in tokens_per_sample]
        for lv in range(1, max_level + 1)
    }
    records = []
    for lv in range(1, max_level + 1):
        cats = [c for c, _ in level_cats[lv]]
        resolved = np.array([r for _, r in level_cats[lv]])
        losses = []
        for rep in range(cv_repeats):
            splitter = StratifiedKFold(n_splits=folds, shuffle=True,
                                       random_state=random_seed + rep)
            for train_idx, val_idx in splitter.split(strat, strat):
                counts: Dict[str, np.ndarray] = {}
                prior = np.zeros(n_classes)
                for i in train_idx:
                    counts.setdefault(cats[i], np.zeros(n_classes))
                    counts[cats[i]][class_index[i]] += 1
                    prior[class_index[i]] += 1
                losses.append(_log_loss(counts, prior,
                                        [cats[i] for i in val_idx],
                                        class_index[val_idx], alpha))
        losses = np.asarray(losses)
        records.append({
            "Method": "Global uniform depth",
            "Parameter": "Level",
            "Parameter Value": float(lv),
            "Label": f"Level {lv}",
            "Number of Categories": len(set(cats)),
            "Mean Depth": float(lv),
            "Resolved Sample Count": int(resolved.sum()),
            "Unresolved Sample Count": int((~resolved).sum()),
            "Resolved Rate": float(resolved.mean()),
            "CV Mean Log Loss": float(losses.mean()),
            "CV SD": float(losses.std(ddof=1)),
            "CV SE": float(losses.std(ddof=1) / np.sqrt(losses.size)),
        })
    return pd.DataFrame(records)


# -------------------------------------------------------- 自适应深度交叉验证

def evaluate_adaptive_lambdas(tree: HaplogroupTree,
                              tokens_per_sample: Sequence[Key],
                              class_index: np.ndarray, n_classes: int,
                              lambdas: Sequence[float], alpha: float,
                              cv_folds: int, cv_repeats: int,
                              random_seed: int) -> pd.DataFrame:
    """对每个候选 λ 做重复分层交叉验证，返回 log-loss 曲线。"""
    folds = _resolve_folds(class_index, cv_folds)
    strat = _stratify_labels(class_index, folds)
    losses: Dict[float, List[float]] = {lam: [] for lam in lambdas}
    sizes: Dict[float, List[float]] = {lam: [] for lam in lambdas}
    depths: Dict[float, List[float]] = {lam: [] for lam in lambdas}

    for rep in range(cv_repeats):
        splitter = StratifiedKFold(n_splits=folds, shuffle=True,
                                   random_state=random_seed + rep)
        for train_idx, val_idx in splitter.split(strat, strat):
            train_tokens = [tokens_per_sample[i] for i in train_idx]
            terminal, subtree = tree.count_matrix(
                train_tokens, class_index[train_idx], n_classes)
            prior = subtree_root_prior(tree, subtree, n_classes)
            val_tokens = [tokens_per_sample[i] for i in val_idx]
            for lam in lambdas:
                cut, _ = fit_tree_cut(tree, terminal, subtree, lam, alpha)
                counts: Dict[str, np.ndarray] = {}
                for tokens, cls in zip(train_tokens, class_index[train_idx]):
                    cat = cut.assign(tokens)
                    counts.setdefault(cat, np.zeros(n_classes))
                    counts[cat][cls] += 1
                losses[lam].append(_log_loss(
                    counts, prior, [cut.assign(t) for t in val_tokens],
                    class_index[val_idx], alpha))
                frontier = cut.categories()
                sizes[lam].append(len(frontier))
                depths[lam].append(float(np.mean([len(k) for _, k, _ in
                                                  frontier])))

    records = []
    for lam in lambdas:
        arr = np.asarray(losses[lam])
        records.append({
            "Method": "Branch-adaptive depth",
            "Parameter": "Lambda",
            "Parameter Value": float(lam),
            "Label": f"Lambda = {lam:g}",
            "Number of Categories": float(np.mean(sizes[lam])),
            "Mean Depth": float(np.mean(depths[lam])),
            "CV Mean Log Loss": float(arr.mean()),
            "CV SD": float(arr.std(ddof=1)),
            "CV SE": float(arr.std(ddof=1) / np.sqrt(arr.size)),
        })
    return pd.DataFrame(records)


def subtree_root_prior(tree: HaplogroupTree, subtree: Dict[Key, np.ndarray],
                       n_classes: int) -> np.ndarray:
    """训练集整体的群体先验。"""
    prior = np.zeros(n_classes)
    for root in tree.roots:
        prior += subtree[root]
    return prior


# ------------------------------------------------------------ 1-SE 规则选择

def apply_one_se_rule(table: pd.DataFrame, simpler_is: str) -> pd.DataFrame:
    """标注 CV 最优与 one-standard-error 规则下最简单的候选。"""
    out = table.copy()
    best = out["CV Mean Log Loss"].idxmin()
    threshold = (out.loc[best, "CV Mean Log Loss"] + out.loc[best, "CV SE"])
    out["Delta From Best"] = (out["CV Mean Log Loss"]
                              - out.loc[best, "CV Mean Log Loss"])
    out["Within One SE"] = out["CV Mean Log Loss"] <= threshold
    within = out[out["Within One SE"]]
    if simpler_is == "smallest":
        pick = within["Parameter Value"].idxmin()
    else:
        pick = within["Parameter Value"].idxmax()
    out["Best Predictive"] = out.index == best
    out["Selected"] = out.index == pick
    return out


# ------------------------------------------------------------------- 主流程

def run(
    input_path: str,
    output_result: str,
    output_figure: str,
    id_column: str = "ID",
    haplogroup_column: str = "Haplogroup",
    population_column: str = "Class",
    label_column: str = "Label",
    max_level: int = 10,
    min_sample_size: int = 20,
    smoothing_alpha: float = 0.5,
    cv_folds: int = 5,
    cv_repeats: int = 10,
    random_seed: int = 42,
    lambda_min: float = 0.5,
    lambda_max: float = 2000.0,
    lambda_steps: int = 24,
    min_resolution_rate: float = 0.0,
    cluster_level: str = "Level 2",
) -> Dict[str, str]:
    """选出全局最优 Level 与各分支自适应深度，并写出结果表。"""
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"输入文件不存在: {input_path}")
    os.makedirs(output_result, exist_ok=True)
    os.makedirs(output_figure, exist_ok=True)

    samples = pd.read_csv(input_path, sep="\t", dtype=str)
    for col in (id_column, haplogroup_column, population_column):
        if col not in samples.columns:
            raise KeyError(f"输入文件缺少列: {col}")
    samples = samples.dropna(subset=[haplogroup_column, population_column])
    sizes = samples[population_column].value_counts()
    samples = samples[samples[population_column].isin(
        sizes[sizes >= min_sample_size].index)].reset_index(drop=True)
    if samples.empty:
        raise ValueError("按最小样本量筛选后无群体剩余")

    display_names = (
        samples[population_column].str.replace("_", " ", regex=False)
        .str.replace(r"(?<=[a-z])(?=[A-Z])", " ", regex=True)
    )
    classes = samples[population_column].to_numpy()
    class_levels = sorted(set(classes))
    class_lookup = {c: i for i, c in enumerate(class_levels)}
    class_index = np.array([class_lookup[c] for c in classes])
    n_classes = len(class_levels)
    tokens_per_sample = [tokenize(h) for h in samples[haplogroup_column]]
    log.info("样本 %d 条，群体 %d 个，原始单倍群 %d 种",
             len(samples), n_classes, samples[haplogroup_column].nunique())

    #* 全局统一深度
    global_cv = evaluate_global_levels(
        tokens_per_sample, class_index, n_classes, max_level,
        smoothing_alpha, cv_folds, cv_repeats, random_seed)
    if min_resolution_rate > 0:
        keep = global_cv["Resolved Rate"] >= min_resolution_rate
        if not keep.any():
            raise ValueError("没有 Level 满足 min_resolution_rate")
        global_cv = global_cv[keep].reset_index(drop=True)
    global_cv = apply_one_se_rule(global_cv, simpler_is="smallest")
    global_level = int(global_cv.loc[global_cv["Selected"],
                                     "Parameter Value"].iloc[0])
    log.info("全局最优 Level（1-SE 规则）: Level %d，CV log-loss %.4f",
             global_level,
             global_cv.loc[global_cv["Selected"],
                           "CV Mean Log Loss"].iloc[0])

    #* 各分支自适应深度
    tree = HaplogroupTree(samples[haplogroup_column])
    lambdas = np.unique(np.round(np.geomspace(
        lambda_min, lambda_max, lambda_steps), 4)).tolist()
    adaptive_cv = evaluate_adaptive_lambdas(
        tree, tokens_per_sample, class_index, n_classes, lambdas,
        smoothing_alpha, cv_folds, cv_repeats, random_seed)
    adaptive_cv = apply_one_se_rule(adaptive_cv, simpler_is="largest")
    best_lambda = float(adaptive_cv.loc[adaptive_cv["Selected"],
                                        "Parameter Value"].iloc[0])
    log.info("自适应深度最优 λ（1-SE 规则）: %g，CV log-loss %.4f",
             best_lambda,
             adaptive_cv.loc[adaptive_cv["Selected"],
                             "CV Mean Log Loss"].iloc[0])

    terminal, subtree = tree.count_matrix(
        tokens_per_sample, class_index, n_classes)
    final_cut, objective = fit_tree_cut(
        tree, terminal, subtree, best_lambda, smoothing_alpha)
    frontier = final_cut.categories()
    log.info("自适应 frontier 类别数 %d，训练目标函数 %.1f",
             len(frontier), objective)

    #* 样本归属
    global_cats = [assign_level(t, global_level) for t in tokens_per_sample]
    adaptive_cats = [final_cut.assign(t) for t in tokens_per_sample]
    if any(c == "" for c in adaptive_cats):
        raise ValueError("存在无法映射到 frontier 的样本")

    assignments = pd.DataFrame({
        "Sample": samples[id_column].to_numpy(),
        "Population": display_names.to_numpy(),
        "Source Haplogroup": samples[haplogroup_column].to_numpy(),
        "Global Haplogroup": [c for c, _ in global_cats],
        "Global Level": global_level,
        "Global Resolution": ["Resolved" if r else "Unresolved"
                              for _, r in global_cats],
        "Adaptive Haplogroup": adaptive_cats,
        "Adaptive Depth": [len(t) - t.count(UNRESOLVED_MARK)
                           for t in adaptive_cats],
        "Adaptive Resolution": ["Unresolved" if c.endswith(UNRESOLVED_MARK)
                                else "Resolved" for c in adaptive_cats],
    })
    depth_lookup = {name: len(key) for name, key, _ in frontier}
    assignments["Adaptive Depth"] = assignments["Adaptive Haplogroup"].map(
        depth_lookup)
    if len(assignments) != len(samples):
        raise ValueError("样本归属表行数与输入不一致")

    #* frontier 明细
    cat_counts = assignments["Adaptive Haplogroup"].value_counts()
    cat_pops = assignments.groupby("Adaptive Haplogroup")[
        "Population"].nunique()
    frontier_tbl = pd.DataFrame([
        {
            "Haplogroup": name,
            "Parent": HaplogroupTree.name(tree.parent[key])
            if key in tree.parent else "",
            "Major Haplogroup": key[0],
            "Depth": len(key),
            "Resolution": "Unresolved" if star else "Resolved",
            "Sample Count": int(cat_counts.get(name, 0)),
            "Number of Populations": int(cat_pops.get(name, 0)),
            "Stop Score": round(final_cut.stop_score[key], 4),
            "Split Score": (round(final_cut.split_score[key], 4)
                            if np.isfinite(final_cut.split_score[key])
                            else float("nan")),
        }
        for name, key, star in frontier
    ]).sort_values(["Major Haplogroup", "Depth", "Haplogroup"])

    #* 选择过程汇总
    summary = pd.concat([global_cv, adaptive_cv], ignore_index=True)
    summary = summary.loc[:, [
        "Method", "Parameter", "Parameter Value", "Label",
        "Number of Categories", "Mean Depth", "Resolved Sample Count",
        "Unresolved Sample Count", "Resolved Rate", "CV Mean Log Loss",
        "CV SD", "CV SE", "Delta From Best", "Within One SE",
        "Best Predictive", "Selected"]]

    paths = {
        "summary": os.path.join(
            output_result, "⭐2-1-Optimal-Depth-Selection.tsv"),
        "assignments": os.path.join(
            output_result, "⭐2-2-Optimal-Depth-Sample-Assignments.tsv"),
        "frontier": os.path.join(
            output_result, "2-3-Adaptive-Depth-Frontier.tsv"),
        "frequency": os.path.join(
            output_figure, "⭐2-4-Optimal-Depth-Frequency.tsv"),
    }
    summary.round(6).to_csv(paths["summary"], sep="\t", index=False)
    assignments.to_csv(paths["assignments"], sep="\t", index=False)
    frontier_tbl.to_csv(paths["frontier"], sep="\t", index=False)

    #* 绘图用频率长表
    freq = build_frequency_table(assignments, tree, cluster_level,
                                 tokens_per_sample, classes)
    if label_column in samples.columns:
        pop_label = pd.Series(samples[label_column].to_numpy(),
                              index=display_names.to_numpy())
        freq["Group Label"] = freq["Population"].map(
            pop_label.groupby(level=0).first())
    else:
        freq["Group Label"] = pd.NA
    freq.round(6).to_csv(paths["frequency"], sep="\t", index=False)

    resolved_rate = (assignments["Adaptive Resolution"] == "Resolved").mean()
    log.info("自适应深度：最浅 %d，最深 %d，平均 %.2f，resolved 比例 %.3f",
             frontier_tbl["Depth"].min(), frontier_tbl["Depth"].max(),
             float(assignments["Adaptive Depth"].mean()), resolved_rate)
    for key, path in paths.items():
        log.info("结果写出 [%s]: %s", key, path)
    return paths


def build_frequency_table(assignments: pd.DataFrame, tree: HaplogroupTree,
                          cluster_level: str,
                          tokens_per_sample: Sequence[Key],
                          classes: np.ndarray) -> pd.DataFrame:
    """把两种方案的样本归属汇总为绘图用的 tidy long 频率表。"""
    frames = []
    schemes = {
        "Global uniform depth": ("Global Haplogroup", "Global Resolution"),
        "Branch-adaptive depth": ("Adaptive Haplogroup",
                                  "Adaptive Resolution"),
    }
    for scheme, (col, res_col) in schemes.items():
        grouped = (
            assignments.groupby(["Population", col, res_col], as_index=False)
            .size().rename(columns={col: "Haplogroup", res_col: "Resolution",
                                    "size": "Count"})
        )
        totals = grouped.groupby("Population")["Count"].transform("sum")
        grouped["Sample Size"] = totals
        grouped["Frequency"] = grouped["Count"] / totals * 100.0
        grouped["Scheme"] = scheme
        frames.append(grouped)
    freq = pd.concat(frames, ignore_index=True)
    freq["Major Haplogroup"] = freq["Haplogroup"].str.extract(r"^([A-Z]+)")
    order = {name: i + 1 for i, name in enumerate(sorted(
        freq["Haplogroup"].unique(),
        key=lambda h: _category_sort_key(
            tokenize(h.rstrip(UNRESOLVED_MARK)))))}
    freq["Category Order"] = freq["Haplogroup"].map(order)

    pop_order = _BASE.cluster_population_order(
        freq[freq["Scheme"] == "Global uniform depth"].assign(
            **{"Population Display": lambda d: d["Population"],
               "Level": cluster_level}),
        cluster_level)
    rank = {p: i + 1 for i, p in enumerate(pop_order)}
    freq["Population Order"] = freq["Population"].map(rank)
    return freq.sort_values(["Scheme", "Population Order", "Category Order"])


# --------------------------------------------------------------- 合成数据自检

def self_test(alpha: float = 0.5) -> None:
    """用小规模合成数据检验算法在几种典型情形下的行为。"""
    rng = np.random.default_rng(0)

    def build(haplos: Sequence[str], pops: Sequence[str]):
        tree = HaplogroupTree(haplos)
        levels = sorted(set(pops))
        idx = np.array([levels.index(p) for p in pops])
        toks = [tokenize(h) for h in haplos]
        term, sub = tree.count_matrix(toks, idx, len(levels))
        return tree, term, sub, toks, idx, len(levels)

    # 1 所有分支都应停在浅层：群体与深层分支无关
    haplos = ["C1a"] * 40 + ["C1b"] * 40
    pops = (["P1", "P2"] * 40)
    tree, term, sub, toks, idx, nc = build(haplos, pops)
    cut, _ = fit_tree_cut(tree, term, sub, 20.0, alpha)
    assert [n for n, _, _ in cut.categories()] == ["C"], "浅层停止失败"

    # 2 只有一个分支值得继续拆
    haplos = ["C1a"] * 40 + ["C1b"] * 40 + ["O1a"] * 40 + ["O1b"] * 40
    pops = ["P1", "P2"] * 40 + ["P3"] * 40 + ["P4"] * 40
    tree, term, sub, toks, idx, nc = build(haplos, pops)
    cut, _ = fit_tree_cut(tree, term, sub, 5.0, alpha)
    names = {n for n, _, _ in cut.categories()}
    assert "C" in names and "O1a" in names and "O1b" in names, f"混合深度失败 {names}"

    # 3 祖先残留：O2 与 O2a 同时存在且 O2 被拆分
    haplos = ["O2"] * 30 + ["O2a"] * 30 + ["O2b"] * 30
    pops = ["P1"] * 30 + ["P2"] * 30 + ["P3"] * 30
    tree, term, sub, toks, idx, nc = build(haplos, pops)
    cut, _ = fit_tree_cut(tree, term, sub, 1.0, alpha)
    names = {n for n, _, _ in cut.categories()}
    assert "O2*" in names, f"祖先残留类别缺失 {names}"
    assert cut.assign(tokenize("O2")) == "O2*", "终止样本未归入星号类别"

    # 4 每个样本恰好归入一个类别，且不丢样本
    assigned = [cut.assign(t) for t in toks]
    assert all(a for a in assigned) and len(assigned) == len(toks), "样本丢失"

    # 5 验证集出现训练未见的更深分支时能回退
    assert cut.assign(tokenize("O2a1b")) in names, "未见分支回退失败"

    # 6 全局 Level 选择：群体完全由 Level 1 决定时应选择 Level 1
    haplos = ["C1a1"] * 60 + ["O1a1"] * 60
    pops = ["P1"] * 60 + ["P2"] * 60
    toks = [tokenize(h) for h in haplos]
    levels = sorted(set(pops))
    idx = np.array([levels.index(p) for p in pops])
    cv = evaluate_global_levels(toks, idx, len(levels), 4, alpha, 5, 2, 42)
    cv = apply_one_se_rule(cv, simpler_is="smallest")
    assert int(cv.loc[cv["Selected"], "Parameter Value"].iloc[0]) == 1, \
        "全局深度未选择 Level 1"

    # 7 群体需要 Level 2 才能区分时不应停在 Level 1
    haplos = ["O1a"] * 60 + ["O2a"] * 60
    pops = ["P1"] * 60 + ["P2"] * 60
    toks = [tokenize(h) for h in haplos]
    idx = np.array([sorted(set(pops)).index(p) for p in pops])
    cv = evaluate_global_levels(toks, idx, 2, 4, alpha, 5, 2, 42)
    cv = apply_one_se_rule(cv, simpler_is="smallest")
    assert int(cv.loc[cv["Selected"], "Parameter Value"].iloc[0]) >= 2, \
        "全局深度未识别 Level 2 的信息量"

    log.info("合成数据自检全部通过")


# ----------------------------------------------------------------------- CLI

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input-path", dest="input_path")
    p.add_argument("--output-result", dest="output_result")
    p.add_argument("--output-figure", dest="output_figure")
    p.add_argument("--id-column", dest="id_column", default="ID")
    p.add_argument("--haplogroup-column", dest="haplogroup_column",
                   default="Haplogroup")
    p.add_argument("--population-column", dest="population_column",
                   default="Class")
    p.add_argument("--label-column", dest="label_column", default="Label")
    p.add_argument("--max-level", dest="max_level", type=int, default=10)
    p.add_argument("--min-sample-size", dest="min_sample_size", type=int,
                   default=20)
    p.add_argument("--smoothing-alpha", dest="smoothing_alpha", type=float,
                   default=0.5)
    p.add_argument("--cv-folds", dest="cv_folds", type=int, default=5)
    p.add_argument("--cv-repeats", dest="cv_repeats", type=int, default=10)
    p.add_argument("--random-seed", dest="random_seed", type=int, default=42)
    p.add_argument("--lambda-min", dest="lambda_min", type=float, default=0.5)
    p.add_argument("--lambda-max", dest="lambda_max", type=float,
                   default=2000.0)
    p.add_argument("--lambda-steps", dest="lambda_steps", type=int, default=24)
    p.add_argument("--min-resolution-rate", dest="min_resolution_rate",
                   type=float, default=0.0)
    p.add_argument("--cluster-level", dest="cluster_level", default="Level 2")
    p.add_argument("--self-test", dest="self_test", action="store_true",
                   help="只运行合成数据自检")
    return p


def main(argv: List[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    args = vars(build_parser().parse_args(argv))
    if args.pop("self_test"):
        self_test(args["smoothing_alpha"])
        return 0
    for required in ("input_path", "output_result", "output_figure"):
        if not args.get(required):
            raise SystemExit(f"缺少必需参数: --{required.replace('_', '-')}")
    run(**args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

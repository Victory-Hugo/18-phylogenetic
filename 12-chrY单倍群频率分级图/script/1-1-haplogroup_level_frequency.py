#!/usr/bin/env python3
"""按 ISOGG 命名规则把 Y 染色体单倍群拆分到各系统发育层级并统计群体频率。

命名解析遵循 ISOGG 交替记号规则：首个大写字母块为主干单倍群，其后
数字块与小写字母块交替出现，每个记号即树上的一层分支。

    C1b1a1a  ->  ('C', '1', 'b', '1', 'a', '1', 'a')
    Level 1 = C, Level 2 = C1, Level 3 = C1b, Level 4 = C1b1

群体顺序可按名称字母序，或按指定层级的单倍群组成做层次聚类，
使组成相似的群体在图上相邻。

样本记号数不足当前层级时保留其祖先归属并标记星号（如 O2*），
而非并入全局 Other，因此每个群体在每个层级的频率之和恒为 100%。

用法::

    python 1-1-haplogroup_level_frequency.py \
        --input input/example.tsv --output-result output/result/1-xxx

也可作为模块导入后调用 :func:`run`。
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from typing import Dict, List, Sequence, Tuple

import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

log = logging.getLogger(__name__)

TOKEN_RE = re.compile(r"[0-9]+|[a-z]+")
ROOT_RE = re.compile(r"^[A-Z]+")

UNRESOLVED_MARK = "*"


class HaplogroupFormatError(ValueError):
    """单倍群名称不符合 ISOGG 命名规则。"""


def tokenize(haplogroup: str) -> Tuple[str, ...]:
    """把单倍群名称拆成系统发育记号序列。"""
    name = haplogroup.strip()
    root = ROOT_RE.match(name)
    if not root:
        raise HaplogroupFormatError(f"缺少主干单倍群字母: {haplogroup!r}")
    tokens = [root.group(0)]
    rest = name[root.end():]
    pos = 0
    for m in TOKEN_RE.finditer(rest):
        if m.start() != pos:
            raise HaplogroupFormatError(f"无法解析的字符: {haplogroup!r}")
        tokens.append(m.group(0))
        pos = m.end()
    if pos != len(rest):
        raise HaplogroupFormatError(f"无法解析的字符: {haplogroup!r}")
    return tuple(tokens)


def assign_level(tokens: Sequence[str], level: int) -> Tuple[str, bool]:
    """返回该记号序列在指定层级的类别名及其是否被完全解析。"""
    if len(tokens) >= level:
        return "".join(tokens[:level]), True
    return "".join(tokens) + UNRESOLVED_MARK, False


def _token_sort_key(token: str) -> Tuple[int, int, str]:
    """数字记号先于字母记号，数字按数值大小排序。"""
    if token.isdigit():
        return (0, int(token), "")
    return (1, 0, token)


def _category_sort_key(tokens: Sequence[str]) -> Tuple:
    """祖先节点排在其后代之前，同层按记号自然顺序。"""
    return (tokens[0],) + tuple(_token_sort_key(t) for t in tokens[1:])


def build_level_table(
    samples: pd.DataFrame,
    haplogroup_column: str,
    population_column: str,
    max_level: int,
) -> pd.DataFrame:
    """展开每个样本在 Level 1..max_level 上的类别归属。"""
    token_map = {h: tokenize(h) for h in samples[haplogroup_column].unique()}
    records = []
    for level in range(1, max_level + 1):
        for hg, tokens in token_map.items():
            name, resolved = assign_level(tokens, level)
            records.append(
                {
                    "Source Haplogroup": hg,
                    "Level": f"Level {level}",
                    "Haplogroup": name,
                    "Major Haplogroup": tokens[0],
                    "Resolution": "Resolved" if resolved else "Unresolved",
                    "Sort Key": _category_sort_key(tokens[:level]),
                }
            )
    expanded = pd.DataFrame(records)
    renamed = samples.rename(columns={haplogroup_column: "Source Haplogroup"})
    return renamed.merge(expanded, on="Source Haplogroup", how="left")


def summarise_frequency(
    long_table: pd.DataFrame,
    population_column: str,
) -> pd.DataFrame:
    """统计每个群体在每个层级上各类别的样本数与百分比频率。"""
    grouped = (
        long_table.groupby(
            [population_column, "Level", "Haplogroup", "Major Haplogroup",
             "Resolution", "Sort Key"],
            as_index=False,
        )
        .size()
        .rename(columns={"size": "Count", population_column: "Population"})
    )
    totals = grouped.groupby(["Population", "Level"])["Count"].transform("sum")
    grouped["Sample Size"] = totals
    grouped["Frequency"] = grouped["Count"] / totals * 100.0
    return grouped


def _global_category_order(table: pd.DataFrame) -> Dict[str, int]:
    """跨层级共享的类别顺序，保证配色与图例含义一致。"""
    unique = table[["Haplogroup", "Sort Key"]].drop_duplicates("Haplogroup")
    ordered = sorted(unique.itertuples(index=False), key=lambda r: r[1])
    return {row.Haplogroup: i + 1 for i, row in enumerate(ordered)}


def cluster_population_order(
    freq: pd.DataFrame,
    cluster_level: str,
) -> List[str]:
    """按指定层级的单倍群组成对群体做层次聚类，返回最优叶序。"""
    sub = freq[freq["Level"] == cluster_level]
    if sub.empty:
        raise ValueError(f"用于聚类的层级不存在: {cluster_level}")
    matrix = (
        sub.pivot_table(index="Population Display", columns="Haplogroup",
                        values="Frequency", fill_value=0.0)
        .sort_index()
    )
    if len(matrix) < 3:
        return list(matrix.index)
    dist = pdist(matrix.to_numpy(), metric="braycurtis")
    tree = linkage(dist, method="average", optimal_ordering=True)
    return list(matrix.index[leaves_list(tree)])


def check_totals(table: pd.DataFrame, tolerance: float = 1e-6) -> None:
    """校验每个群体每个层级的频率之和为 100%。"""
    sums = table.groupby(["Population", "Level"])["Frequency"].sum()
    bad = sums[(sums - 100.0).abs() > tolerance]
    if not bad.empty:
        raise ValueError(f"频率之和不等于 100%: \n{bad}")
    log.info("频率校验通过: %d 个群体 × 层级组合均为 100%%", len(sums))


def run(
    input_path: str,
    output_result: str,
    id_column: str = "ID",
    haplogroup_column: str = "Haplogroup",
    population_column: str = "Class",
    label_column: str = "Label",
    max_level: int = 4,
    min_sample_size: int = 20,
    population_order: str = "alphabetical",
    cluster_level: str = "Level 2",
) -> str:
    """主流程：读入样本表，输出层级频率长表，返回主结果文件路径。"""
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"输入文件不存在: {input_path}")
    if max_level < 1:
        raise ValueError("max_level 必须 >= 1")

    os.makedirs(output_result, exist_ok=True)
    samples = pd.read_csv(input_path, sep="\t", dtype=str)
    for col in (id_column, haplogroup_column, population_column):
        if col not in samples.columns:
            raise KeyError(f"输入文件缺少列: {col}")
    use_label = bool(label_column) and label_column in samples.columns
    if label_column and not use_label:
        log.warning("输入文件没有 %s 列，群体标签留空", label_column)
    samples = samples.dropna(subset=[haplogroup_column, population_column])
    log.info("读入样本 %d 条，群体 %d 个",
             len(samples), samples[population_column].nunique())

    sizes = samples[population_column].value_counts()
    keep = sizes[sizes >= min_sample_size].index
    dropped = sorted(set(sizes.index) - set(keep))
    if dropped:
        log.info("样本量 < %d 而剔除的群体 %d 个: %s",
                 min_sample_size, len(dropped), ", ".join(dropped[:10])
                 + (" ..." if len(dropped) > 10 else ""))
    samples = samples[samples[population_column].isin(keep)]
    if samples.empty:
        raise ValueError("按最小样本量筛选后无群体剩余")

    long_table = build_level_table(
        samples, haplogroup_column, population_column, max_level)
    freq = summarise_frequency(long_table, population_column)
    check_totals(freq)

    order_map = _global_category_order(freq)
    freq["Category Order"] = freq["Haplogroup"].map(order_map)
    freq["Population Display"] = (
        freq["Population"].str.replace("_", " ", regex=False)
        .str.replace(r"(?<=[a-z])(?=[A-Z])", " ", regex=True)
    )
    if population_order == "alphabetical":
        ordered_pops = sorted(freq["Population Display"].unique())
    elif population_order == "cluster":
        ordered_pops = cluster_population_order(freq, cluster_level)
        log.info("按 %s 的单倍群组成聚类排序群体（Bray-Curtis 距离）",
                 cluster_level)
    else:
        raise ValueError(f"暂不支持的群体排序方式: {population_order}")
    pop_rank = {p: i + 1 for i, p in enumerate(ordered_pops)}
    freq["Population Order"] = freq["Population Display"].map(pop_rank)
    if use_label:
        pop_label = samples.groupby(population_column)[label_column].first()
        freq["Group Label"] = freq["Population"].map(pop_label)
    else:
        freq["Group Label"] = pd.NA

    freq["Level Depth"] = (
        freq["Level"].str.replace("Level ", "", regex=False).astype(int))
    out = (
        freq.drop(columns=["Sort Key", "Population"])
        .rename(columns={"Population Display": "Population"})
        .loc[:, ["Population", "Population Order", "Group Label",
                 "Sample Size", "Level", "Level Depth",
                 "Haplogroup", "Major Haplogroup", "Resolution",
                 "Category Order", "Count", "Frequency"]]
        .sort_values(["Level Depth", "Population Order", "Category Order"])
    )
    out["Frequency"] = out["Frequency"].round(6)

    path = os.path.join(
        output_result, "⭐1-1-Haplogroup-Frequency-By-Level.tsv")
    out.to_csv(path, sep="\t", index=False)
    for level, sub in out.sort_values("Level Depth").groupby(
            "Level", sort=False):
        log.info("%s: %d 个类别, %d 个群体",
                 level, sub["Haplogroup"].nunique(), sub["Population"].nunique())
    log.info("结果写出: %s", path)
    return path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input-path", dest="input_path", required=True)
    p.add_argument("--output-result", dest="output_result", required=True)
    p.add_argument("--id-column", dest="id_column", default="ID")
    p.add_argument("--haplogroup-column", dest="haplogroup_column",
                   default="Haplogroup")
    p.add_argument("--population-column", dest="population_column",
                   default="Class")
    p.add_argument("--label-column", dest="label_column", default="Label",
                   help="群体分组标签列，缺失时忽略")
    p.add_argument("--max-level", dest="max_level", type=int, default=4)
    p.add_argument("--min-sample-size", dest="min_sample_size", type=int,
                   default=20)
    p.add_argument("--population-order", dest="population_order",
                   default="alphabetical", choices=["alphabetical", "cluster"])
    p.add_argument("--cluster-level", dest="cluster_level", default="Level 2",
                   help="population-order=cluster 时用于聚类的系统发育层级")
    return p


def main(argv: List[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    args = build_parser().parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

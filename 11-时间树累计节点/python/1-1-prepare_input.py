#!/usr/bin/env python3
"""1-1 准备输入：校验树与 metadata，生成样本分组长表、分母与配色。

分组维度来自 conf 的 groupings 段：一个维度就是 metadata 的一个列名。该列的类别自动从数据
里取全部取值，按样本量降序排列，样本量最大的画在最上面。

产物
  output/result/1-prepare_input/sample-groups.tsv   tidy long：样本 × 维度 × 类别
  output/result/1-prepare_input/sample-counts.tsv   各类别样本量（标准化定标的分母）
  output/result/1-prepare_input/group-design.tsv    脊线顺序、填充色与描边色

CLI
  python python/1-1-prepare_input.py --conf conf/1-conf.yaml --output-result ...
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import Counter
from pathlib import Path

import yaml

log = logging.getLogger(__name__)

TIP_PATTERN = re.compile(r"[(,]\s*([^(),:;]+):")
GROUP_FIELDS = ["Sample", "Grouping", "Group"]
ALL_SAMPLES = "All samples"
DESIGN_FIELDS = ["Grouping", "Group", "Order", "Colour", "Outline"]


def _read_meta(meta_path: Path) -> tuple[list[str], list[dict]]:
    with open(meta_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _to_rgb(value: str) -> tuple[int, int, int]:
    h = value.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def _interpolate(ramp: list[str], n: int) -> list[str]:
    """在颜色条上等距取 n 个色：相邻类别的颜色由同一条色带线性插值而来，因此不会跳色。"""
    anchors = [_to_rgb(c) for c in ramp]
    if n == 1:
        return [ramp[len(ramp) // 2]]
    picked = []
    for i in range(n):
        position = i * (len(anchors) - 1) / (n - 1)
        low = min(int(position), len(anchors) - 2)
        frac = position - low
        rgb = [round(a + (b - a) * frac) for a, b in zip(anchors[low], anchors[low + 1])]
        picked.append("#{:02X}{:02X}{:02X}".format(*rgb))
    return picked


def _darken(colour: str, factor: float) -> str:
    """描边色：同色相压暗一档，使脊线轮廓不至于用一个与填充无关的深灰去切割配色。"""
    return "#{:02X}{:02X}{:02X}".format(*(round(c * factor) for c in _to_rgb(colour)))


def run(conf: str, output_result: str) -> int:
    settings = yaml.safe_load(Path(conf).read_text(encoding="utf-8"))
    paths, cohort_conf = settings["paths"], settings["cohort"]
    groupings = list(settings["groupings"])
    plot_conf = settings["plot"]

    result_path = Path(output_result)
    result_path.mkdir(parents=True, exist_ok=True)

    tree_path, meta_path = Path(paths["tree"]), Path(paths["meta"])
    header, records = _read_meta(meta_path)

    # 源列存在性校验：登记的列名若在 metadata 中不存在，该维度会静默取到空值而整个消失，
    # 下游少出一组图却不报错。这类失败必须在第一步就炸掉。
    id_column = cohort_conf["id_column"]
    absent = sorted((set(groupings) | {id_column}) - set(header))
    if absent:
        raise ValueError(f"metadata 缺少所需的列：{absent}；实际列名为 {header}")

    # 实际分析的样本就是树的 tip 与 metadata 样本编号的交集
    tip_ids = set(TIP_PATTERN.findall(tree_path.read_text()))
    meta_ids = {str(r[id_column]).strip() for r in records}
    shared = tip_ids & meta_ids
    if not shared:
        raise ValueError(
            f"树的 tip 名与 metadata 的 {id_column} 列没有任何交集；"
            f"tip 示例 {sorted(tip_ids)[:3]}，{id_column} 示例 {sorted(meta_ids)[:3]}")
    log.info("树 tip %d 个，metadata %d 行，取交集后 %d 个样本进入分析"
             "（树上未匹配 %d，metadata 中未上树 %d）",
             len(tip_ids), len(records), len(shared),
             len(tip_ids - meta_ids), len(meta_ids - tip_ids))

    # 类别自动从数据里取：先数一遍各列取值的样本量，再按样本量降序定出脊线顺序
    counts: Counter[tuple[str, str]] = Counter()
    on_tree = [r for r in records if str(r[id_column]).strip() in shared]
    for record in on_tree:
        for column in groupings:
            value = str(record.get(column, "")).strip()
            if value:
                counts[(column, value)] += 1

    # 少于 2 例的类别在诱导子树里一个合并事件都构不成；更要紧的是抽稀的目标样本量取该维度
    # 最小类别的样本量，留着长尾里的小类别会把整个维度抽稀到无法出图。
    min_group = max(2, int(settings["analysis"].get("min_group_size", 2)))
    levels: dict[str, list[str]] = {}
    for column in groupings:
        kept = sorted(((v, n) for (c, v), n in counts.items() if c == column and n >= min_group),
                      key=lambda x: (-x[1], x[0]))
        dropped = [(v, n) for (c, v), n in counts.items() if c == column and n < min_group]
        if dropped:
            log.info("  %-28s 剔除 %d 个不足 %d 例的类别：%s", column, len(dropped), min_group,
                     ", ".join(f"{v}({n})" for v, n in
                               sorted(dropped, key=lambda x: -x[1])))
        levels[column] = [v for v, _ in kept]
        if kept:
            log.info("  %-28s %d 个类别，抽稀目标样本量受限于 %s（%d 例）",
                     column, len(kept), kept[-1][0], kept[-1][1])

    group_rows = []
    for record in on_tree:
        sample_id = str(record[id_column]).strip()
        for column in groupings:
            value = str(record.get(column, "")).strip()
            if value in levels[column]:
                group_rows.append({"Sample": sample_id, "Grouping": column, "Group": value})

    # 全样本总览：单类别维度，涵盖树上的每一个样本。单类别没有要均衡的对象，
    # 下游会跳过抽稀直接用全部个体。
    if settings.get("include_all_samples", False):
        groupings.insert(0, ALL_SAMPLES)
        levels[ALL_SAMPLES] = [ALL_SAMPLES]
        counts[(ALL_SAMPLES, ALL_SAMPLES)] = len(on_tree)
        group_rows.extend({"Sample": str(r[id_column]).strip(),
                           "Grouping": ALL_SAMPLES, "Group": ALL_SAMPLES} for r in on_tree)
        log.info("  %-28s 全样本总览，%d 个个体（单类别，不抽稀）", ALL_SAMPLES, len(on_tree))

    with open(result_path / "sample-groups.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=GROUP_FIELDS)
        writer.writeheader()
        writer.writerows(group_rows)

    # 维度级空值校验：源列存在但取值全部落在 levels 白名单之外（取值改了拼写、或队列筛选后
    # 一个都不剩），同样会让该维度静默消失。
    empty = [column for column in groupings if not levels[column]]
    if empty:
        raise ValueError(
            f"以下维度没有任何类别达到 min_group_size = {min_group}：{empty}")

    with open(result_path / "sample-counts.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Grouping", "Group", "Samples"])
        for column in groupings:
            for value in levels[column]:
                writer.writerow([column, value, counts[(column, value)]])

    # 脊线顺序与配色：category 模式下每条脊线一个纯色（沿色带插值），描边取同色压暗一档；
    # age 模式的填充色由节点年龄决定，故两列留空。
    by_category = plot_conf.get("fill_mode", "category") == "category"
    darken = float(plot_conf.get("outline_darken", 0.72))
    with open(result_path / "group-design.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(DESIGN_FIELDS)
        for column in groupings:
            values = levels[column]
            colours = _interpolate(plot_conf["colour_ramp"], len(values)) if by_category \
                else [""] * len(values)
            for order, (value, colour) in enumerate(zip(values, colours), start=1):
                writer.writerow([column, value, order, colour,
                                 _darken(colour, darken) if colour else ""])

    log.info("样本分组长表 %d 行，%d 个维度 → %s",
             len(group_rows), len(groupings), result_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="准备节点山脊图的输入数据与样本分组")
    parser.add_argument("--conf", required=True)
    parser.add_argument("--output-result", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    return run(**vars(build_parser().parse_args(argv)))


if __name__ == "__main__":
    raise SystemExit(main())

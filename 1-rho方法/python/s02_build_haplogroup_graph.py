"""
s02_build_haplogroup_graph.py
解析 PhyloTree Build 17 JSON，构建单倍群祖先-后代图。

用法（CLI）：
    python s02_build_haplogroup_graph.py \
        --phylotree data/phylotree_index_withacc.json \
        --output output/stage_1/intermediate/haplogroup_graph.json

用法（import）：
    from s02_build_haplogroup_graph import run
    run(phylotree="...", output="...")
"""

import argparse
import json
import logging
from collections import defaultdict, deque
from pathlib import Path

log = logging.getLogger(__name__)


def _build_graph(raw: dict) -> dict:
    """
    从原始PhyloTree JSON构建包含 parent/children/lineage 的完整图。

    原始JSON结构（每个节点）：
        {
            "parent": "H1",
            "lineage": ["mt-MRCA(RSRS)", ..., "H1", "H1a"],
            "level": 9,
            "mutations": "...",
            ...
        }

    返回格式：
        {
            "H1a": {
                "parent": "H1",
                "lineage": [...],
                "children": ["H1a1", "H1a2"]
            },
            ...
        }
    """
    graph = {}

    # 第一遍：初始化所有节点，记录parent和lineage
    for name, info in raw.items():
        graph[name] = {
            "parent":  info.get("parent", None),
            "lineage": info.get("lineage", []),
            "children": [],
        }

    # 第二遍：根据parent关系填充children列表
    for name, node in graph.items():
        parent = node["parent"]
        if parent and parent in graph:
            graph[parent]["children"].append(name)

    log.info(f"  共解析 {len(graph)} 个单倍群节点")
    return graph


def get_all_descendants(graph: dict, haplogroup: str) -> set:
    """
    BFS 获取 haplogroup 及其所有后代节点名称（含自身）。
    使用BFS避免递归栈溢出。
    """
    result = set()
    queue = deque([haplogroup])
    while queue:
        node = queue.popleft()
        if node in result:
            continue
        result.add(node)
        for child in graph.get(node, {}).get("children", []):
            queue.append(child)
    return result


def run(
    phylotree: str,
    output: str,
) -> int:
    """
    解析PhyloTree JSON，构建单倍群图，写出JSON。

    参数：
        phylotree: PhyloTree Build 17 JSON路径
        output:    输出JSON路径

    返回：
        0 表示成功
    """
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    log.info(f"读取PhyloTree JSON: {phylotree}")
    with open(phylotree) as f:
        root = json.load(f)

    # PhyloTree JSON 顶层结构为 {"haplogroups": {...}, ...}
    if isinstance(root, dict) and "haplogroups" in root:
        raw = root["haplogroups"]
        log.info(f"  使用顶层 'haplogroups' 键（含 {len(raw)} 个节点）")
    else:
        raw = root

    graph = _build_graph(raw)

    # 统计有children的节点数
    n_internal = sum(1 for v in graph.values() if v["children"])
    log.info(f"  内部节点（有子节点）: {n_internal}")
    log.info(f"  叶节点（无子节点）: {len(graph) - n_internal}")

    with open(output, "w") as f:
        json.dump(graph, f, indent=2, ensure_ascii=False)
    log.info(f"单倍群图已写出至: {output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="构建单倍群祖先-后代图（s02）")
    p.add_argument("--phylotree", required=True, help="PhyloTree JSON路径")
    p.add_argument("--output",    required=True, help="输出JSON路径")
    return p


def main(argv=None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(
        phylotree=args.phylotree,
        output=args.output,
    )


if __name__ == "__main__":
    raise SystemExit(main())

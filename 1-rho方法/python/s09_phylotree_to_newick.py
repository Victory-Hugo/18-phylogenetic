"""
s09_phylotree_to_newick.py
将 haplogroup_graph.json（父子关系图）转为 IQ-TREE 固定拓扑 Newick 文件。

节点名称使用与 s08_prepare_ml_inputs.py 相同的净化规则，保证与 FASTA 头一致。

用法（CLI）：
    python s09_phylotree_to_newick.py \
        --graph  output/stage_1/intermediate/haplogroup_graph.json \
        --map    output/stage_2/intermediate/ml/name_sanitize_map.tsv \
        --output output/stage_2/intermediate/ml/phylotree_fixed.nwk \
        --root   "mt-MRCA(RSRS)"

用法（import）：
    from s09_phylotree_to_newick import run
    run(graph_json=..., map_tsv=..., output_nwk=..., root=...)
"""

import argparse
import json
import logging
from pathlib import Path

log = logging.getLogger(__name__)

ROOT_DEFAULT = "mt-MRCA(RSRS)"


def _build_newick_iterative(graph: dict, root: str, sanitize_map: dict) -> str:
    """迭代后序遍历生成 Newick 字符串，避免 5000+ 层递归栈溢出。"""
    result: dict[str, str] = {}
    # 栈元素：(node, processed)
    stack = [(root, False)]

    while stack:
        node, processed = stack.pop()
        if processed:
            children = graph[node]["children"]
            tip = sanitize_map[node]
            if not children:
                result[node] = tip
            else:
                inner = ",".join(result[c] for c in children)
                result[node] = f"({inner}){tip}"
        else:
            stack.append((node, True))
            # 子节点逆序入栈，保证正序处理
            for child in reversed(graph[node]["children"]):
                stack.append((child, False))

    return result[root] + ";"


def run(
    graph_json: str,
    map_tsv: str,
    output_nwk: str,
    root: str = ROOT_DEFAULT,
) -> None:
    graph_json = Path(graph_json)
    map_tsv = Path(map_tsv)
    output_nwk = Path(output_nwk)
    output_nwk.parent.mkdir(parents=True, exist_ok=True)

    with open(graph_json) as f:
        graph = json.loads(f.read())
    log.info("读入单倍群图：%d 节点", len(graph))

    if root not in graph:
        raise ValueError(f"根节点 '{root}' 不在图中，可用节点样例：{list(graph.keys())[:5]}")

    sanitize_map: dict[str, str] = {}
    with open(map_tsv) as f:
        next(f)  # 跳过 header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                sanitize_map[parts[0]] = parts[1]
    log.info("读入名称映射：%d 条", len(sanitize_map))

    # 检查图中所有节点都在映射表内
    missing = [k for k in graph if k not in sanitize_map]
    if missing:
        raise ValueError(f"{len(missing)} 个图节点未在映射表中找到，样例：{missing[:5]}")

    newick = _build_newick_iterative(graph, root, sanitize_map)
    output_nwk.write_text(newick)
    log.info(
        "Newick 写出：%s（%d 字节，尖端数 = 叶节点数）",
        output_nwk,
        len(newick),
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--graph",  required=True, help="haplogroup_graph.json")
    p.add_argument("--map",    required=True, help="name_sanitize_map.tsv")
    p.add_argument("--output", required=True, help="输出 Newick 文件")
    p.add_argument("--root",   default=ROOT_DEFAULT, help=f"根节点名（默认: {ROOT_DEFAULT}）")
    return p


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args()
    run(
        graph_json=args.graph,
        map_tsv=args.map,
        output_nwk=args.output,
        root=args.root,
    )


if __name__ == "__main__":
    main()

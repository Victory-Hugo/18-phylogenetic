#!/usr/bin/env python3
"""BEAST v10 大树流程的共享函数。"""

from __future__ import annotations

import csv
import logging
import math
import re
import subprocess
from io import StringIO
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape, quoteattr

from Bio import Phylo, SeqIO


log = logging.getLogger(__name__)


class PipelineError(RuntimeError):
    """流程可预期错误。"""


def normalize_haplogroup(value: str) -> str:
    """将带突变后缀的 mtDNA haplogroup 还原为 phylotree 主节点名。"""
    hap = value.strip().replace(" ", "")
    hap = re.sub(r"!+$", "", hap)
    if "_" in hap:
        hap = hap.split("_", 1)[0]
        hap = re.sub(r"!+$", "", hap)
    return hap


def haplogroup_in_node(haplogroup: str, node: str, haplogroups: dict[str, Any]) -> bool:
    """判断一个 haplogroup 是否属于 phylotree JSON 中的目标节点。"""
    normalized = normalize_haplogroup(haplogroup)
    if normalized == node:
        return True
    info = haplogroups.get(normalized)
    if not info:
        return False
    return node in info.get("lineage", [])


def read_fasta_ids(path: Path) -> tuple[list[str], dict[str, int]]:
    """读取 FASTA ID 和序列长度。"""
    ids: list[str] = []
    lengths: dict[str, int] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        sample_id = record.id
        ids.append(sample_id)
        lengths[sample_id] = len(record.seq)
    if not ids:
        raise PipelineError(f"FASTA 为空: {path}")
    if len(ids) != len(set(ids)):
        raise PipelineError(f"FASTA 存在重复 ID: {path}")
    return ids, lengths


def read_id_hap(path: Path) -> dict[str, str]:
    """读取两列 ID-Hap TSV。"""
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PipelineError(f"{path}:{line_no} 不是两列 TSV")
            sample_id, hap = parts
            if sample_id in mapping:
                raise PipelineError(f"{path}:{line_no} 存在重复 ID: {sample_id}")
            mapping[sample_id] = hap
    if not mapping:
        raise PipelineError(f"ID-Hap TSV 为空: {path}")
    return mapping


def read_tree(path: Path):
    """读取 Newick 树。"""
    try:
        return Phylo.read(str(path), "newick")
    except Exception as exc:  # noqa: BLE001
        raise PipelineError(f"无法读取 Newick 树: {path}: {exc}") from exc


def tree_tip_names(tree) -> list[str]:
    """返回树 tip 名称并检查重复。"""
    names = [terminal.name for terminal in tree.get_terminals()]
    if any(name is None for name in names):
        raise PipelineError("树中存在无名称 tip")
    if len(names) != len(set(names)):
        raise PipelineError("树中存在重复 tip")
    return names


def check_monophyly(tree_path: Path, membership: set[str]) -> dict[str, Any]:
    """检查给定样本集合在树中是否为单系。"""
    tree = read_tree(tree_path)
    tips = {terminal.name: terminal for terminal in tree.get_terminals()}
    missing = sorted(membership - set(tips))
    if missing:
        raise PipelineError(f"单系检查成员不在树中: {missing[:10]}")
    if not membership:
        raise PipelineError("单系检查成员为空")
    mrca = tree.common_ancestor([tips[name] for name in sorted(membership)])
    descendant_taxa = {terminal.name for terminal in mrca.get_terminals()}
    extra = sorted(descendant_taxa - membership)
    return {
        "member_count": len(membership),
        "mrca_tip_count": len(descendant_taxa),
        "extra_taxa": extra,
        "is_monophyletic": len(extra) == 0,
    }


def write_filtered_fasta(input_fasta: Path, keep_ids: set[str], output_fasta: Path) -> None:
    """按树中 tip 过滤 FASTA。"""
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with output_fasta.open("w", encoding="utf-8") as handle:
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            if record.id in keep_ids:
                SeqIO.write(record, handle, "fasta")
                kept += 1
    if kept != len(keep_ids):
        raise PipelineError(f"过滤 FASTA 后数量不一致: expected={len(keep_ids)} observed={kept}")


def max_root_to_tip_distance(tree) -> float:
    """计算时间树最大 root-to-tip 距离。"""
    depths = tree.depths()
    terminal_depths = [depths[terminal] for terminal in tree.get_terminals()]
    max_depth = max(terminal_depths)
    if max_depth <= 0:
        raise PipelineError("定年树 root-to-tip 距离非正，无法度量树高")
    return float(max_depth)


def write_lsd2_date_file(
    *,
    taxa: list[str],
    calibration_taxa: dict[str, list[str]],
    node_ages: dict[str, float],
    output_path: Path,
) -> None:
    """写 IQ-TREE/LSD2 的 date 文件。

    约定 present=0、过去为负年（BP）。每行 ``taxon<TAB>date``；
    祖先节点用「逗号连接的 taxa 列表 + 日期」约束其 MRCA。
    root 不写入本文件，改用 ``--date-root`` 直接锁定。
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{taxon}\t0" for taxon in taxa]
    for node, members in calibration_taxa.items():
        if not members:
            continue
        members_token = ",".join(members)
        lines.append(f"{members_token}\t-{node_ages[node]:.0f}")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_iqtree_lsd2(
    *,
    iqtree_bin: Path,
    alignment: Path,
    ml_tree: Path,
    date_file: Path,
    root_age: float,
    outdir: Path,
    model: str,
    threads: str,
    outgroup: str,
    date_options: str,
) -> tuple[Path, str]:
    """用 IQ-TREE 内置 LSD2 对固定拓扑的 ML 树做最小二乘定年。"""
    if not iqtree_bin.exists():
        raise PipelineError(f"找不到 IQ-TREE: {iqtree_bin}")
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / "lsd"
    command = [
        str(iqtree_bin),
        "-s",
        str(alignment),
        "-te",
        str(ml_tree),
        "--tree-fix",
        "--date",
        str(date_file),
        "--date-tip",
        "0",
        "--date-root",
        f"-{root_age:.0f}",
        "-m",
        model,
        "-T",
        str(threads),
        "--prefix",
        str(prefix),
        "-redo",
    ]
    if outgroup.strip():
        # 以外群固定生根，避免 LSD 全分支搜根重排内群拓扑。
        command += ["-o", outgroup.strip()]
    if date_options.strip():
        command += ["--date-options", date_options.strip()]
    log.info("运行 IQ-TREE LSD2: %s", " ".join(command))
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    log_text = result.stdout + "\n" + result.stderr
    if result.returncode != 0:
        raise PipelineError(f"IQ-TREE LSD2 运行失败，退出码 {result.returncode}\n{log_text[-4000:]}")
    timetree_nexus = prefix.with_name(prefix.name + ".timetree.nex")
    if not timetree_nexus.exists():
        matches = sorted(outdir.glob("*timetree*"))
        raise PipelineError(f"IQ-TREE 未产生时间树 NEXUS，目录内容: {[p.name for p in matches]}")
    return timetree_nexus, log_text


def convert_tree_to_newick(tree_path: Path, output_newick: Path) -> str:
    """把 NEXUS/Newick 时间树转换成 BEAST 可读 Newick。"""
    fmt = "nexus" if tree_path.suffix.lower() in {".nex", ".nexus"} else "newick"
    tree = Phylo.read(str(tree_path), fmt)
    for clade in tree.find_clades():
        # IQ-TREE 时间树 NEXUS 把 [&date=...] 写进 comment / confidence，
        # 必须清除，否则导出的 newick 会带注释，BEAST 无法解析。
        clade.comment = None
        if clade.confidence is not None and not isinstance(clade.confidence, (int, float)):
            clade.confidence = None
    output_newick.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(output_newick), "newick")
    text = output_newick.read_text(encoding="utf-8").strip()
    if not text.endswith(";"):
        text += ";"
    return text


def clade_sets(tree) -> set[frozenset[str]]:
    """返回树中所有内部 clade 的 tip 集合。"""
    sets: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        names = frozenset(terminal.name for terminal in clade.get_terminals())
        if names:
            sets.add(names)
    return sets


def validate_starting_tree(original_tree_path: Path, timetree_newick: Path, expected_tips: set[str]) -> dict[str, Any]:
    """验证 startingTree 与 ML dataTree 的 tip 和拓扑兼容性。"""
    original = read_tree(original_tree_path)
    timed = read_tree(timetree_newick)
    original_tips = set(tree_tip_names(original))
    timed_tips = set(tree_tip_names(timed))
    if timed_tips != expected_tips or original_tips != timed_tips:
        raise PipelineError("定年树 tip 集合与输入 ML 树不一致")
    negative = [
        clade.name or "<internal>"
        for clade in timed.find_clades()
        if clade.branch_length is not None and clade.branch_length < -1e-9
    ]
    if negative:
        raise PipelineError(f"定年树存在负分支长度: {negative[:10]}")
    original_clades = clade_sets(original)
    timed_clades = clade_sets(timed)
    missing_clades = original_clades - timed_clades
    if missing_clades:
        sample = [sorted(c)[:5] for c in list(missing_clades)[:5]]
        raise PipelineError(f"定年树缺失 ML 树的 {len(missing_clades)} 个 clade，拓扑不兼容，示例: {sample}")
    root_height = max_root_to_tip_distance(timed)
    return {
        "root_height": root_height,
        "tip_count": len(timed_tips),
        "internal_clade_count": len(timed_clades),
        "extra_resolved_clades": len(timed_clades - original_clades),
    }


def escape_newick_for_xml(newick: str) -> str:
    """转义写入 XML 的 Newick 文本。"""
    return escape(newick.strip())


def taxa_block(taxa: list[str], block_id: str, use_idref: bool = False) -> str:
    """生成 BEAST taxa XML block。"""
    lines = [f'    <taxa id="{block_id}">']
    attr = "idref" if use_idref else "id"
    for taxon in taxa:
        lines.append(f"        <taxon {attr}={quoteattr(taxon)}/>")
    lines.append("    </taxa>")
    return "\n".join(lines)


def render_beast_xml(
    *,
    template_path: Path,
    taxa: list[str],
    starting_tree: str,
    data_tree: str,
    calibration_taxa: dict[str, list[str]],
    nodes: list[dict[str, Any]],
    scale: float,
    clock_rate: float,
    chain_length: int,
    log_every: int,
    tree_log_every: int,
    output_prefix: str,
    skygrid_points: int,
    skygrid_cutoff: float,
    skygrid_init_logpop: float,
) -> str:
    """用 conf 中的 XML 模板渲染 BEAST v10 Thorney XML。"""
    if not template_path.exists():
        raise PipelineError(f"缺少 XML 模板: {template_path}")
    skygrid_dim = max(2, int(skygrid_points))
    grid_intervals = skygrid_dim - 1
    skyline_dim = skygrid_dim
    skyline_init_popsize = math.exp(float(skygrid_init_logpop))
    skyline_upper_popsize = skyline_init_popsize * 1000.0
    taxa_xml = taxa_block(taxa, "taxa")
    # 按 yaml 中 calibration.nodes 的顺序遍历，循环生成任意数量（可嵌套）的
    # 校准 taxa block / tmrcaStatistic / normalPrior。root 节点用全树高度，
    # 无成员 taxa block；其余节点用 <mrca><taxa idref=节点名/></mrca>。
    cal_taxa_blocks: list[str] = []
    tmrca_blocks: list[str] = []
    prior_blocks: list[str] = []
    for node in nodes:
        name = node["name"]
        if node.get("is_root"):
            tmrca_blocks.append(
                f'    <tmrcaStatistic id="age({name})" absolute="true">\n'
                f'        <treeModel idref="treeModel"/>\n'
                f"    </tmrcaStatistic>"
            )
        else:
            cal_taxa_blocks.append(taxa_block(calibration_taxa[name], name, use_idref=True))
            tmrca_blocks.append(
                f'    <tmrcaStatistic id="age({name})" absolute="true">\n'
                f"        <mrca>\n"
                f'            <taxa idref="{name}"/>\n'
                f"        </mrca>\n"
                f'        <treeModel idref="treeModel"/>\n'
                f"    </tmrcaStatistic>"
            )
        prior_blocks.append(
            f'                <normalPrior id="{name}.prior" mean="{node["age"]}" stdev="{node["stdev"]}">\n'
            f'                    <tmrcaStatistic idref="age({name})"/>\n'
            f"                </normalPrior>"
        )
    calibration_taxa_xml = "\n".join(cal_taxa_blocks)
    tmrca_stats = "\n".join(tmrca_blocks).lstrip()
    normal_priors = "\n".join(prior_blocks)
    template = template_path.read_text(encoding="utf-8")
    return template.format(
        taxa_xml=taxa_xml,
        calibration_taxa_xml=calibration_taxa_xml,
        starting_tree=escape_newick_for_xml(starting_tree),
        data_tree=escape_newick_for_xml(data_tree),
        tmrca_stats=tmrca_stats,
        normal_priors=normal_priors,
        scale=float(scale),
        clock_rate=f"{clock_rate:.12g}",
        chain_length=int(chain_length),
        log_every=int(log_every),
        tree_log_every=int(tree_log_every),
        output_prefix=output_prefix,
        skygrid_dim=skygrid_dim,
        grid_intervals=grid_intervals,
        skygrid_cutoff=float(skygrid_cutoff),
        skygrid_init_logpop=float(skygrid_init_logpop),
        skyline_dim=skyline_dim,
        skyline_init_popsize=skyline_init_popsize,
        skyline_upper_popsize=skyline_upper_popsize,
    )


def write_validation_table(path: Path, rows: list[tuple[str, str, str]]) -> None:
    """写流程验证表。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["item", "value", "status"])
        writer.writerows(rows)


def read_validation_table(path: Path) -> dict[str, str]:
    """读取验证表为 item 到 value 的映射。"""
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            values[row["item"]] = row["value"]
    return values


def write_membership(path: Path, rows: list[tuple[str, str, str]]) -> None:
    """写校准成员表（long 格式：node, sample_id, haplogroup）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["node", "sample_id", "haplogroup"])
        writer.writerows(rows)


def read_membership(path: Path) -> tuple[dict[str, list[str]], dict[str, str]]:
    """读取校准成员表，返回 (node -> 样本列表, sample_id -> haplogroup)。

    取代旧的 read_calibration_taxa：不再假设固定的 M/N 节点，支持任意数量
    与嵌套节点（同一样本可在多个嵌套节点下出现多行）。
    """
    membership: dict[str, list[str]] = {}
    hap_map: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "node" not in reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise PipelineError(f"成员表缺少 node/sample_id 列: {path}")
        for row in reader:
            node = row["node"]
            sample_id = row["sample_id"]
            membership.setdefault(node, []).append(sample_id)
            hap_map[sample_id] = row.get("haplogroup", "NA")
    if not membership:
        raise PipelineError(f"成员表为空: {path}")
    for node in membership:
        membership[node] = sorted(membership[node])
    return membership, hap_map


def parse_calibration_nodes(spec: str) -> list[dict[str, Any]]:
    """解析 load_config.sh 传入的紧凑校准节点串。

    格式：``name:age:stdev:is_root``，以 ``;`` 分隔多个节点；age/stdev 可为空
    （某些步骤不需要 stdev）。返回保持原顺序的节点 dict 列表。
    """
    nodes: list[dict[str, Any]] = []
    for token in spec.split(";"):
        token = token.strip()
        if not token:
            continue
        parts = token.split(":")
        if len(parts) != 4:
            raise PipelineError(f"校准节点编码格式错误: {token!r}，应为 name:age:stdev:is_root")
        name, age_s, stdev_s, is_root_s = (p.strip() for p in parts)
        if not name:
            raise PipelineError(f"校准节点缺少名称: {token!r}")
        nodes.append({
            "name": name,
            "age": float(age_s) if age_s else None,
            "stdev": float(stdev_s) if stdev_s else None,
            "is_root": is_root_s == "1",
        })
    if not nodes:
        raise PipelineError("未解析到任何校准节点，请检查 yaml 的 calibration.nodes")
    roots = [n for n in nodes if n["is_root"]]
    if len(roots) != 1:
        raise PipelineError(f"calibration.nodes 必须且只能有一个 is_root 节点，当前 {len(roots)} 个")
    names = [n["name"] for n in nodes]
    if len(names) != len(set(names)):
        raise PipelineError(f"calibration.nodes 节点名重复: {names}")
    return nodes


def validate_calibration_nesting(nodes: list[dict[str, Any]], calibration_taxa: dict[str, set[str]]) -> None:
    """校验校准节点的嵌套关系与年龄单调性。

    - 任意两个非 root 节点的成员集合必须「互斥」或「真嵌套」，禁止部分重叠。
    - 若 A 真包含 B（A 是 B 的祖先 clade），要求 age(A) > age(B)。
    - root 年龄必须大于所有非 root 节点年龄。
    """
    age = {n["name"]: n["age"] for n in nodes}
    root = next(n for n in nodes if n["is_root"])
    names = list(calibration_taxa)
    for i, a in enumerate(names):
        sa = calibration_taxa[a]
        for b in names[i + 1:]:
            sb = calibration_taxa[b]
            if not (sa & sb):
                continue
            if sb < sa:
                if not age[a] > age[b]:
                    raise PipelineError(f"嵌套校准年龄非单调：{a}(age={age[a]:.0f}) 包含 {b}(age={age[b]:.0f})，祖先必须更老")
            elif sa < sb:
                if not age[b] > age[a]:
                    raise PipelineError(f"嵌套校准年龄非单调：{b}(age={age[b]:.0f}) 包含 {a}(age={age[a]:.0f})，祖先必须更老")
            else:
                raise PipelineError(f"校准节点 {a} 与 {b} 成员部分重叠且互不包含，违反 clade 嵌套/互斥约束")
    for name in calibration_taxa:
        if not age[root["name"]] > age[name]:
            raise PipelineError(f"root 年龄({age[root['name']]:.0f}) 必须大于节点 {name} 年龄({age[name]:.0f})")


def tree_to_newick_text(tree_path: Path) -> str:
    """读取树文件并返回单行 Newick 文本。"""
    text = tree_path.read_text(encoding="utf-8").strip()
    if not text.endswith(";"):
        text += ";"
    return text


def newick_text_to_tree(newick: str):
    """从 Newick 字符串读取树对象。"""
    return Phylo.read(StringIO(newick), "newick")

#!/usr/bin/env python3
"""BEAST(1.x / v10.5.0) stage2 公共工具模块。

本模块为「拆树 BEAST 合树」流水线第 2 阶段提供准备输入、运行 BEAST/TreeAnnotator、
解析 MCC 树所需的公共函数与常量，对应原 PAML 版本的 ``paml_common.py``。

设计要点：
- BEAST 严格分子钟要求**起始树超度量**（所有 tip 同时代）。因此本模块在准备阶段会先把
  stage1 给出的固定拓扑子树转换为超度量树，再写入 BEAST XML 的 ``<newick>`` 起始树。
- 设 strict clock 速率固定为 1.0（不给该参数加算子），BEAST 时间树的枝长单位即
  「期望替换数/位点」，与 baseml ``clock=1`` 的命名树等价，可直接供下游 merge 阶段使用。
- XML 算子集合只保留节点高度/根高度/模型参数算子，剔除一切改变拓扑的算子，从而**固定拓扑**。

该文件为库模块（被各 stage2 脚本 import），不单独作为命令行入口运行。
"""

from __future__ import annotations

import csv
import datetime as dt
import fcntl
import io
import json
import os
import re
import subprocess
import threading
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio import Phylo
from Bio.SeqRecord import SeqRecord

from phylo_merge_common import load_baseml_summary, write_rows  # noqa: F401  (write_rows 供外部复用)
from phylo_split_common import PipelineError, compute_tip_hash, read_newick_tree, write_table
from phylo_ultrastandard_common import extend_terminal_branches_to_height


# ---------------------------------------------------------------------------
# 列定义（对应 PAML 版本的同名常量，字段名保持中性的 subtree_id 以便下游复用）
# ---------------------------------------------------------------------------
BEAST_JOB_MANIFEST_COLUMNS = [
    "beast_subtree_id",
    "job_dir",
    "xmlfile",
    "starting_tree_file",
    "seqfile",
    "trees_file",
    "mcc_file",
    "file_prefix",
    "expected_tip_count",
    "tip_hash",
    "seq_id_strategy",
]

BEAST_RUN_STATUS_COLUMNS = [
    "beast_subtree_id",
    "job_dir",
    "trees_file",
    "mcc_file",
    "return_code",
    "status",
    "detail",
]

BEAST_PARSE_STATUS_COLUMNS = [
    "beast_subtree_id",
    "mcc_file",
    "analysis_tree_file",
    "status",
    "detail",
]

# 各子树/backbone MCC 节点的高度与 95% HPD 汇总表列定义。
# height 单位为 substitutions/site（clock 速率固定 1.0），供 stage5 把 HPD
# 按 clade（tip 集合 tip_hash）映射并线性传播到最终年代树。
BEAST_NODE_HPD_COLUMNS = [
    "source_id",        # 来源：beast_subtree_id 或 "backbone"
    "tip_hash",         # 该节点 clade 的规范化 tip 集合哈希（order-independent）
    "n_tips",           # clade 包含的 tip 数
    "height_subst",     # 节点高度（substitutions/site）
    "hpd_low_subst",    # 高度 95% HPD 下界（substitutions/site），缺失为空
    "hpd_high_subst",   # 高度 95% HPD 上界（substitutions/site），缺失为空
    "tip_names_json",   # 该 clade 的 tip 名（JSON 列表，便于审计）
]


# ---------------------------------------------------------------------------
# 序列与命名映射（与 paml_common 等价，确保 tip 名与 FASTA 记录一一对应）
# ---------------------------------------------------------------------------
def resolve_seq_id_strategy(strategy: str) -> str:
    """校验并归一化序列 ID 映射策略。"""
    normalized = str(strategy).strip()
    if normalized not in {"exact", "prefix_before_first_underscore"}:
        raise PipelineError(f"Unsupported beast.seq_id_strategy: {strategy}")
    return normalized


def read_fasta_record_map(fasta_path: Path) -> Dict[str, SeqRecord]:
    """读取 FASTA 为 {record_id: SeqRecord} 字典，禁止重复 ID。"""
    records: Dict[str, SeqRecord] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        if record.id in records:
            raise PipelineError(f"Duplicate FASTA record id: {record.id}")
        records[record.id] = record
    if not records:
        raise PipelineError(f"Input FASTA is empty: {fasta_path}")
    return records


def map_tip_to_source_id(
    tip_names: Sequence[str],
    fasta_ids: Iterable[str],
    seq_id_strategy: str,
) -> Dict[str, str]:
    """把树 tip 名映射到 FASTA 记录 ID。

    - ``exact``：tip 名必须与 FASTA ID 完全一致。
    - ``prefix_before_first_underscore``：优先精确匹配，否则取第一个下划线前的前缀。
    """
    fasta_id_set = set(fasta_ids)
    mapping: Dict[str, str] = {}
    used_source_ids: Dict[str, str] = {}

    for tip_name in tip_names:
        if seq_id_strategy == "exact":
            source_id = tip_name
        else:
            source_id = tip_name if tip_name in fasta_id_set else tip_name.split("_", 1)[0]
        if source_id not in fasta_id_set:
            raise PipelineError(f"Could not map tree tip to FASTA id: {tip_name} -> {source_id}")
        owner = used_source_ids.get(source_id)
        if owner is not None and owner != tip_name:
            raise PipelineError(
                f"FASTA id {source_id} maps to multiple tree tips ({owner}, {tip_name}); mapping is not unique."
            )
        used_source_ids[source_id] = tip_name
        mapping[tip_name] = source_id
    return mapping


def write_subtree_fasta(
    tip_names: Sequence[str],
    fasta_record_map: Dict[str, SeqRecord],
    tip_to_source_id: Dict[str, str],
    destination: Path,
) -> None:
    """按 tip 顺序写出子树 FASTA（序列头使用 tip 名，便于与 BEAST taxa 对应）。"""
    destination.parent.mkdir(parents=True, exist_ok=True)
    records: List[SeqRecord] = []
    for tip_name in tip_names:
        source_id = tip_to_source_id[tip_name]
        source_record = fasta_record_map[source_id]
        records.append(
            SeqRecord(seq=source_record.seq, id=tip_name, name=tip_name, description="")
        )
    SeqIO.write(records, str(destination), "fasta")


def normalize_newick_text(raw_text: str) -> str:
    """把多行 Newick 文本压成单行并保证以分号结尾。"""
    normalized = "".join(line.strip() for line in raw_text.splitlines() if line.strip())
    if not normalized.endswith(";"):
        normalized = f"{normalized};"
    return normalized


# ---------------------------------------------------------------------------
# 起始树（固定拓扑 + 超度量化），供 BEAST XML 注入
# ---------------------------------------------------------------------------
def build_ultrametric_starting_newick(
    source_tree_path: Path,
    expected_tip_count: int,
    min_branch_length: float,
    ultrametric_tolerance: float,
) -> str:
    """读取 stage1 固定拓扑子树，转换为超度量树并返回单行 Newick 字符串。

    BEAST 严格时钟要求所有现存 tip 同处时间 0，即起始树必须超度量。这里复用
    ``extend_terminal_branches_to_height`` 把终端枝延长到最大根到叶深度实现超度量化，
    拓扑保持不变。
    """
    tree = read_newick_tree(source_tree_path)
    actual_tip_count = len(tree.get_terminals())
    if actual_tip_count != int(expected_tip_count):
        raise PipelineError(
            f"Tree tip count mismatch for {source_tree_path}: "
            f"expected {expected_tip_count}, got {actual_tip_count}"
        )
    ultrametric_tree, _, _ = extend_terminal_branches_to_height(
        tree,
        min_branch_length=min_branch_length,
        tolerance=ultrametric_tolerance,
    )
    return _tree_to_single_line_newick(ultrametric_tree)


def _tree_to_single_line_newick(tree) -> str:
    """把 Bio.Phylo Tree 序列化为单行 Newick（不含内部节点名/置信度）。"""
    import io

    from Bio import Phylo

    buffer = io.StringIO()
    Phylo.write(tree, buffer, "newick", format_branch_length="%.12f")
    return normalize_newick_text(buffer.getvalue())


# ---------------------------------------------------------------------------
# BEAST XML 渲染
# ---------------------------------------------------------------------------
def build_taxa_block(tip_names: Sequence[str]) -> str:
    """生成 BEAST ``<taxa>`` 内部的 ``<taxon>`` 行集合。"""
    return "\n".join(f'\t\t<taxon id="{name}"/>' for name in tip_names)


def build_alignment_block(
    tip_names: Sequence[str],
    fasta_record_map: Dict[str, SeqRecord],
    tip_to_source_id: Dict[str, str],
) -> str:
    """生成 BEAST ``<alignment>`` 内部的 ``<sequence>`` 行集合。"""
    chunks: List[str] = []
    for tip_name in tip_names:
        source_record = fasta_record_map[tip_to_source_id[tip_name]]
        sequence_text = str(source_record.seq)
        chunks.append(
            f'\t\t<sequence>\n\t\t\t<taxon idref="{tip_name}"/>\n\t\t\t{sequence_text}\n\t\t</sequence>'
        )
    return "\n".join(chunks)


def render_beast_xml(
    template_path: Path,
    taxa_block: str,
    alignment_block: str,
    starting_tree_newick: str,
    file_prefix: str,
    chain_length: int,
    log_every: int,
    destination: Path,
    clock_rate: float = 1.0,
) -> None:
    """把模板中的占位符替换为具体内容，写出 BEAST 作业 XML。

    模板占位符（双花括号）：
    ``{{TAXA}}``、``{{ALIGNMENT}}``、``{{STARTING_TREE}}``、``{{FILE_PREFIX}}``、
    ``{{CHAIN_LENGTH}}``、``{{LOG_EVERY}}``、``{{CLOCK_RATE}}``。

    ``clock_rate`` 为 strict clock 速率（固定，不加算子）。默认 1.0 时枝长单位为
    「期望替换数/位点」，与下游合并约定一致；设为人类 mtDNA 速率（如 2.53e-8）则 BEAST
    以「年」为单位输出，但最终年代仍由 stage5 node_calibration 决定（见 DEC-009）。
    """
    text = template_path.read_text(encoding="utf-8")
    replacements = {
        "{{TAXA}}": taxa_block,
        "{{ALIGNMENT}}": alignment_block,
        "{{STARTING_TREE}}": starting_tree_newick,
        "{{FILE_PREFIX}}": file_prefix,
        "{{CHAIN_LENGTH}}": str(int(chain_length)),
        "{{LOG_EVERY}}": str(int(log_every)),
        "{{CLOCK_RATE}}": ("%g" % float(clock_rate)),
    }
    for key, value in replacements.items():
        text = text.replace(key, value)
    remaining = re.findall(r"\{\{[A-Z_]+\}\}", text)
    if remaining:
        raise PipelineError(f"Unresolved placeholders in rendered BEAST XML: {sorted(set(remaining))}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# 校准/UCLD XML 渲染与热启动
# ---------------------------------------------------------------------------
def classify_m_n_tips(
    tip_names: Sequence[str],
    id_hap_tsv: Path,
    phylotree_json: Path,
) -> Tuple[List[str], List[str], List[str]]:
    """根据样本单倍群与 phylotree lineage 把 tip 分为 M/N。

    返回 ``(M_tips, N_tips, missing_tips)``。RSRS 作为外群不参与 M/N 分类。
    """
    id_to_hap: Dict[str, str] = {}
    for line in Path(id_hap_tsv).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        id_to_hap[parts[0]] = parts[1]

    phylotree = json.loads(Path(phylotree_json).read_text(encoding="utf-8"))["haplogroups"]

    def resolve_lineage(haplogroup: Optional[str]) -> Optional[Sequence[str]]:
        if not haplogroup:
            return None
        if haplogroup in phylotree:
            return phylotree[haplogroup]["lineage"]
        for cut in range(len(haplogroup), 0, -1):
            prefix = haplogroup[:cut]
            if prefix in phylotree:
                return phylotree[prefix]["lineage"]
        return None

    m_tips: List[str] = []
    n_tips: List[str] = []
    missing: List[str] = []
    for tip_name in tip_names:
        if tip_name == "RSRS":
            continue
        lineage = resolve_lineage(id_to_hap.get(tip_name))
        if lineage is None:
            missing.append(tip_name)
        elif "M" in lineage:
            m_tips.append(tip_name)
        elif "N" in lineage:
            n_tips.append(tip_name)
        else:
            missing.append(tip_name)
    return m_tips, n_tips, missing


def build_calibration_block(m_tips: Sequence[str], n_tips: Sequence[str]) -> str:
    """生成 M/N 校准 taxa 与 tmrcaStatistic XML 块。"""
    _validate_calibration_tip_set(m_tips, "M")
    _validate_calibration_tip_set(n_tips, "N")

    def taxa_block(name: str, tips: Sequence[str]) -> str:
        rows = "\n".join(f'\t\t<taxon idref="{tip}"/>' for tip in tips)
        return f'\t<taxa id="{name}">\n{rows}\n\t</taxa>'

    def tmrca_block(stat_id: str, taxa_id: str) -> str:
        return (
            f'\t<tmrcaStatistic id="{stat_id}" absolute="true">\n'
            f'\t\t<mrca><taxa idref="{taxa_id}"/></mrca>\n'
            f'\t\t<treeModel idref="treeModel"/>\n'
            f'\t</tmrcaStatistic>'
        )

    return "\n".join(
        [
            taxa_block("M_clade", m_tips),
            taxa_block("N_clade", n_tips),
            tmrca_block("tmrca_M", "M_clade"),
            tmrca_block("tmrca_N", "N_clade"),
        ]
    )


def validate_calibration_clades(tree, m_tips: Sequence[str], n_tips: Sequence[str], source_id: str) -> None:
    """确认 M/N 校准集合在当前固定拓扑中可定义且单系。"""
    _validate_calibration_tip_set(m_tips, "M")
    _validate_calibration_tip_set(n_tips, "N")
    tree_tips = {tip.name for tip in tree.get_terminals()}
    for label, tips in (("M", m_tips), ("N", n_tips)):
        missing = sorted(set(tips) - tree_tips)
        if missing:
            raise PipelineError(f"{source_id}: {label} calibration tips missing from tree: {missing[:8]}")
        ancestor = tree.common_ancestor(list(tips))
        descendants = {tip.name for tip in ancestor.get_terminals()}
        extra = sorted(descendants - set(tips) - {"RSRS"})
        if extra:
            raise PipelineError(
                f"{source_id}: {label} calibration clade is not monophyletic; extra descendants: {extra[:8]}"
            )


def build_calibrated_starting_newick(
    source_tree_path: Path,
    expected_tip_count: int,
    min_branch_length: float,
    ultrametric_tolerance: float,
    m_tips: Sequence[str],
    n_tips: Sequence[str],
    root_age: float,
    m_age: float,
    n_age: float,
    warm_start: bool = True,
    source_id: str = "beast_job",
) -> Tuple[str, float, float]:
    """构建校准/UCLD XML 使用的年尺度热启动树。

    返回 ``(starting_newick_years, ucld_mean_start, root_height_subst)``。
    """
    subst_newick = build_ultrametric_starting_newick(
        source_tree_path=source_tree_path,
        expected_tip_count=expected_tip_count,
        min_branch_length=min_branch_length,
        ultrametric_tolerance=ultrametric_tolerance,
    )
    tree = Phylo.read(io.StringIO(subst_newick), "newick")
    validate_calibration_clades(tree, m_tips, n_tips, source_id)

    _, root_subst = _node_depths_and_root_height(tree)
    if root_subst <= 0:
        raise PipelineError(f"{source_id}: calibrated starting tree has non-positive root height.")
    if root_age <= 0:
        raise PipelineError(f"{source_id}: root calibration age must be positive.")

    scale = float(root_age) / root_subst
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= scale

    if warm_start:
        _rescale_calibrated_clade(tree, m_tips, float(m_age), "M", source_id)
        _rescale_calibrated_clade(tree, n_tips, float(n_age), "N", source_id)

    buffer = io.StringIO()
    Phylo.write(tree, buffer, "newick", format_branch_length="%.6f")
    return normalize_newick_text(buffer.getvalue()), root_subst / float(root_age), root_subst


def render_calibrated_beast_xml(
    template_path: Path,
    taxa_block: str,
    alignment_block: str,
    starting_tree_newick: str,
    calibration_block: str,
    file_prefix: str,
    chain_length: int,
    log_every: int,
    ucld_mean: float,
    popsize: float,
    root_age: float,
    root_stdev: float,
    m_age: float,
    m_stdev: float,
    n_age: float,
    n_stdev: float,
    destination: Path,
) -> None:
    """渲染校准/UCLD BEAST XML。"""
    text = Path(template_path).read_text(encoding="utf-8")
    replacements = {
        "{{TAXA}}": taxa_block,
        "{{ALIGNMENT}}": alignment_block,
        "{{STARTING_TREE}}": starting_tree_newick,
        "{{CALIBRATION_BLOCK}}": calibration_block,
        "{{UCLD_MEAN}}": _format_float(ucld_mean),
        "{{POPSIZE}}": _format_float(popsize),
        "{{ROOT_AGE}}": _format_float(root_age),
        "{{ROOT_SD}}": _format_float(root_stdev),
        "{{M_AGE}}": _format_float(m_age),
        "{{M_SD}}": _format_float(m_stdev),
        "{{N_AGE}}": _format_float(n_age),
        "{{N_SD}}": _format_float(n_stdev),
        "{{FILE_PREFIX}}": file_prefix,
        "{{CHAIN_LENGTH}}": str(int(chain_length)),
        "{{LOG_EVERY}}": str(int(log_every)),
    }
    for key, value in replacements.items():
        text = text.replace(key, value)
    remaining = re.findall(r"\{\{[A-Z_]+\}\}", text)
    if remaining:
        raise PipelineError(f"Unresolved placeholders in rendered calibrated BEAST XML: {sorted(set(remaining))}")
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    Path(destination).write_text(text, encoding="utf-8")


def _validate_calibration_tip_set(tips: Sequence[str], label: str) -> None:
    if len(set(tips)) != len(list(tips)):
        raise PipelineError(f"{label} calibration tips contain duplicates.")
    if len(tips) < 2:
        raise PipelineError(f"{label} calibration requires at least two tips; got {len(tips)}.")


def _node_depths_and_root_height(tree) -> Tuple[Dict[object, float], float]:
    depths: Dict[object, float] = {}

    def visit(clade, depth: float) -> None:
        depths[clade] = depth
        for child in clade.clades:
            visit(child, depth + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    root_height = max(depths[tip] for tip in tree.get_terminals())
    return depths, root_height


def _rescale_calibrated_clade(tree, tips: Sequence[str], target_age: float, label: str, source_id: str) -> None:
    if target_age <= 0:
        raise PipelineError(f"{source_id}: {label} calibration age must be positive.")
    ancestor = tree.common_ancestor(list(tips))
    depths, root_height = _node_depths_and_root_height(tree)
    current_age = root_height - depths[ancestor]
    if current_age <= 0:
        raise PipelineError(f"{source_id}: {label} current clade age is non-positive.")
    if target_age >= root_height:
        raise PipelineError(f"{source_id}: {label} calibration age must be younger than root age.")
    factor = float(target_age) / current_age
    for clade in ancestor.find_clades():
        if clade is ancestor:
            continue
        if clade.branch_length is not None:
            clade.branch_length *= factor
    ancestor.branch_length = float(ancestor.branch_length or 0.0) + (current_age - float(target_age))


def _format_float(value: float) -> str:
    return f"{float(value):.12g}"


# ---------------------------------------------------------------------------
# 运行 BEAST 与 TreeAnnotator
# ---------------------------------------------------------------------------
def run_beast_job(
    beast_bin: str,
    xml_path: Path,
    job_dir: Path,
    seed: Optional[int] = None,
    cpu_ids: Optional[Sequence[int]] = None,
    beagle_instances: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """在 ``job_dir`` 中运行 BEAST，``-overwrite`` 覆盖旧结果。

    通过 CPU 亲和把作业绑定到指定核。``beagle_instances`` 指定 BEAGLE 实例数（=每作业线程数）：
    单分区分析靠 ``-beagle_instances N`` 把位点切成 N 份并行计算，从而用满分配的多核；
    若为 None 则按绑定核数 ``len(cpu_ids)`` 自动设定（>1 时才启用，单核作业不加该参数）。
    实测大子树（~1936 tips）16 实例约 3.4x 加速；少作业多线程可避开多作业并行时的内存带宽竞争。
    """
    child_affinity = tuple(int(cpu_id) for cpu_id in cpu_ids) if cpu_ids else None

    def _apply_affinity():
        if child_affinity is not None:
            os.sched_setaffinity(0, set(child_affinity))

    if beagle_instances is None and child_affinity is not None:
        beagle_instances = len(child_affinity)

    command = [beast_bin, "-overwrite"]
    if seed is not None:
        command += ["-seed", str(int(seed))]
    if beagle_instances is not None and int(beagle_instances) > 1:
        command += ["-beagle_instances", str(int(beagle_instances))]
    command += [xml_path.name]

    return subprocess.run(
        command,
        cwd=str(job_dir),
        check=False,
        text=True,
        capture_output=True,
        preexec_fn=_apply_affinity if child_affinity is not None else None,
    )


def run_treeannotator(
    treeannotator_bin: str,
    trees_file: Path,
    mcc_file: Path,
    job_dir: Path,
    burnin_percent: float = 10.0,
    heights: str = "mean",
    target_tree: Optional[Path] = None,
) -> subprocess.CompletedProcess:
    """运行 TreeAnnotator 把后验树集汇总为 MCC 树（NEXUS）。

    本流程固定拓扑：若作业目录存在 ``starting_tree.nwk``（或显式给定 ``target_tree``），
    则用 ``-target`` 指定该固定拓扑作为标注目标，直接在已知拓扑上标注后验均值高度。
    这样既符合「固定拓扑」语义，又绕开 TreeAnnotator v10.5 默认 HIPSTR 目标树重建时
    偶发的 "Error annotating tree: null" 问题。
    """
    if target_tree is None:
        candidate = job_dir / "starting_tree.nwk"
        target_tree = candidate if candidate.exists() else None

    command = [
        treeannotator_bin,
        "-burninTrees",
        str(_burnin_trees_count(trees_file, burnin_percent)),
        "-heights",
        str(heights),
    ]
    if target_tree is not None:
        command += ["-target", Path(target_tree).name]
    command += [trees_file.name, mcc_file.name]

    return subprocess.run(
        command,
        cwd=str(job_dir),
        check=False,
        text=True,
        capture_output=True,
    )


def _burnin_trees_count(trees_file: Path, burnin_percent: float) -> int:
    """根据 .trees 文件中的 ``tree`` 行数与 burnin 百分比计算要丢弃的树数。"""
    if not trees_file.exists():
        return 0
    tree_count = 0
    with trees_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.lstrip()
            if stripped.startswith("tree ") or stripped.startswith("tree\t"):
                tree_count += 1
    return max(0, int(tree_count * float(burnin_percent) / 100.0))


# ---------------------------------------------------------------------------
# 解析 MCC（NEXUS）→ Newick
# ---------------------------------------------------------------------------
def _read_mcc_translate_and_tree(mcc_file: Path) -> Tuple[Dict[str, str], str]:
    """从 MCC NEXUS 中提取 ``Translate`` 映射表与首个 ``tree`` 行（含注释）。"""
    if not mcc_file.exists() or mcc_file.stat().st_size == 0:
        raise PipelineError(f"MCC tree file missing or empty: {mcc_file}")

    lines = mcc_file.read_text(encoding="utf-8", errors="replace").splitlines()
    translate: Dict[str, str] = {}
    tree_line: Optional[str] = None
    in_translate = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower().startswith("translate"):
            in_translate = True
            continue
        if in_translate:
            if stripped == ";" or stripped.startswith(")"):
                in_translate = False
                continue
            token = stripped.rstrip(",").strip()
            if not token:
                continue
            if token == ";":
                in_translate = False
                continue
            parts = token.split(None, 1)
            if len(parts) == 2:
                key, name = parts
                translate[key] = name.strip().strip("'\"")
            if stripped.endswith(";"):
                in_translate = False
            continue
        if stripped.lower().startswith("tree "):
            tree_line = stripped

    if tree_line is None:
        raise PipelineError(f"Could not locate a 'tree' line in MCC file: {mcc_file}")
    return translate, tree_line


def parse_mcc_to_newick(mcc_file: Path) -> str:
    """把 TreeAnnotator 输出的 MCC NEXUS 树解析为单行 Newick 字符串。

    处理要点：
    - 解析 ``Translate`` 表，把树串中的数字 tip 标签还原为真实 taxon 名。
    - 剥离 FigTree 风格的 ``[&...]`` 注释（节点高度、HPD 区间等）。
    - 仅保留拓扑与枝长（枝长即 substitutions/site，因 clock 速率固定为 1.0）。
    """
    translate, tree_line = _read_mcc_translate_and_tree(mcc_file)

    # 取 '=' 之后的 newick 部分
    newick = tree_line.split("=", 1)[1].strip()
    # 去掉根状态标记 [&R] / [&U]
    newick = re.sub(r"\[&[RU]\]\s*", "", newick)
    # 去掉所有 [&...] 注释（含可能的嵌套大括号 HPD）
    newick = _strip_square_bracket_comments(newick)
    newick = normalize_newick_text(newick)

    # 把数字 tip 标签还原为 taxon 名（仅替换出现在 '(' 或 ',' 之后、':' 之前的整数标签）
    if translate:
        def _replace(match: "re.Match") -> str:
            key = match.group(2)
            name = translate.get(key)
            return f"{match.group(1)}{name}" if name is not None else match.group(0)

        newick = re.sub(r"([(,])(\d+)(?=:)", _replace, newick)

    return newick


def _strip_square_bracket_comments(text: str) -> str:
    """删除 ``[...]`` 注释；正确处理注释内部可能出现的 ``{...}``（HPD 区间）。"""
    result: List[str] = []
    depth = 0
    for char in text:
        if char == "[":
            depth += 1
            continue
        if char == "]":
            if depth > 0:
                depth -= 1
            continue
        if depth == 0:
            result.append(char)
    return "".join(result)


# ---------------------------------------------------------------------------
# 解析 MCC 节点注释（高度 + 95% HPD），供 stage5 把 HPD 传播到最终年代树
# ---------------------------------------------------------------------------
class _MccNode:
    """MCC newick 解析的轻量节点：保留 FigTree 风格 [&...] 注释。"""

    __slots__ = ("children", "name", "comment", "tip_names")

    def __init__(self) -> None:
        self.children: List["_MccNode"] = []
        self.name: Optional[str] = None
        self.comment: Optional[str] = None
        self.tip_names: List[str] = []


def _read_bracket_comment(text: str, i: int) -> Tuple[str, int]:
    """从 ``text[i] == '['`` 处读取一段 ``[...]`` 注释，返回 (注释体, 结束后下标)。

    正确处理可能出现的嵌套方括号；HPD 的 ``{...}`` 原样保留在注释体内。
    """
    assert text[i] == "["
    depth = 0
    start = i
    while i < len(text):
        char = text[i]
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                return text[start + 1 : i], i + 1
        i += 1
    raise PipelineError("Unbalanced '[' in MCC newick comment.")


def _read_node_label(text: str, i: int) -> Tuple[str, int]:
    """读取叶/内部节点标签（支持单引号包裹），停在 ``[:,()`` 或末尾。"""
    if i < len(text) and text[i] == "'":
        j = i + 1
        chars: List[str] = []
        while j < len(text):
            if text[j] == "'":
                if j + 1 < len(text) and text[j + 1] == "'":  # 转义的单引号
                    chars.append("'")
                    j += 2
                    continue
                return "".join(chars), j + 1
            chars.append(text[j])
            j += 1
        raise PipelineError("Unbalanced quote in MCC newick label.")
    j = i
    while j < len(text) and text[j] not in "[:,();":
        j += 1
    return text[i:j], j


def _parse_mcc_clade(text: str, i: int, translate: Dict[str, str]) -> Tuple[_MccNode, int]:
    """递归下降解析带注释的 newick 子树。"""
    node = _MccNode()
    if text[i] == "(":
        i += 1
        while True:
            child, i = _parse_mcc_clade(text, i, translate)
            node.children.append(child)
            if i < len(text) and text[i] == ",":
                i += 1
                continue
            if i < len(text) and text[i] == ")":
                i += 1
                break
            raise PipelineError("Malformed MCC newick (expected ',' or ')').")
        # 内部节点可能有标签（BEAST 通常没有），随后是节点注释
        label, i = _read_node_label(text, i)
        if label.strip():
            node.name = label.strip().strip("'\"")
        if i < len(text) and text[i] == "[":
            node.comment, i = _read_bracket_comment(text, i)
    else:
        label, i = _read_node_label(text, i)
        raw = label.strip().strip("'\"")
        node.name = translate.get(raw, raw)
        if i < len(text) and text[i] == "[":
            node.comment, i = _read_bracket_comment(text, i)

    # 枝长 ``:数值``，其后可能还有一段（速率）注释，跳过
    if i < len(text) and text[i] == ":":
        i += 1
        j = i
        while j < len(text) and text[j] not in "[,();":
            j += 1
        i = j
        if i < len(text) and text[i] == "[":
            _, i = _read_bracket_comment(text, i)
    return node, i


def _assign_tip_names(node: _MccNode) -> List[str]:
    """后序遍历，为每个节点填充其后代 tip 名集合。"""
    if not node.children:
        node.tip_names = [node.name] if node.name else []
        return node.tip_names
    collected: List[str] = []
    for child in node.children:
        collected.extend(_assign_tip_names(child))
    node.tip_names = collected
    return collected


# TreeAnnotator 把节点高度写成 height_mean / height_median（取决于 -heights 选项），
# 并无裸 ``height=`` 键；按 mean→median→裸 height 的优先级提取。HPD 为 height_95%_HPD。
_HEIGHT_RES = [
    re.compile(r"height_mean=([-0-9.eEnN]+)"),
    re.compile(r"height_median=([-0-9.eEnN]+)"),
    re.compile(r"[&,]height=([-0-9.eEnN]+)"),
]
_HPD_RE = re.compile(r"height_95%_HPD=\{([^}]*)\}")


def _extract_height_hpd(comment: Optional[str]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """从节点注释中提取 (height, hpd_low, hpd_high)，缺失返回 None。"""
    if not comment:
        return None, None, None
    height: Optional[float] = None
    for pattern in _HEIGHT_RES:
        match = pattern.search(comment)
        if match:
            try:
                height = float(match.group(1))
            except ValueError:
                height = None
            break
    hpd_low: Optional[float] = None
    hpd_high: Optional[float] = None
    hpd_match = _HPD_RE.search(comment)
    if hpd_match:
        parts = [p.strip() for p in hpd_match.group(1).split(",") if p.strip()]
        if len(parts) == 2:
            try:
                hpd_low, hpd_high = float(parts[0]), float(parts[1])
            except ValueError:
                hpd_low, hpd_high = None, None
    return height, hpd_low, hpd_high


def parse_mcc_node_annotations(mcc_file: Path) -> List[Dict[str, object]]:
    """解析 MCC NEXUS，提取每个内部节点的高度与 95% HPD（substitutions/site）。

    返回的每条记录含规范化的 ``tip_hash``（order-independent，复用
    ``compute_tip_hash``）、tip 数、tip 名列表、节点高度及 HPD 上下界。仅返回
    含 ≥2 个 tip 的内部节点（叶节点高度恒为 0，无需记录）。
    """
    translate, tree_line = _read_mcc_translate_and_tree(mcc_file)
    newick = tree_line.split("=", 1)[1].strip()
    newick = re.sub(r"^\[&[RU]\]\s*", "", newick)

    root, _ = _parse_mcc_clade(newick, 0, translate)
    _assign_tip_names(root)

    annotations: List[Dict[str, object]] = []
    stack = [root]
    while stack:
        node = stack.pop()
        stack.extend(node.children)
        if len(node.tip_names) < 2:
            continue
        height, hpd_low, hpd_high = _extract_height_hpd(node.comment)
        annotations.append(
            {
                "tip_hash": compute_tip_hash(node.tip_names),
                "n_tips": len(node.tip_names),
                "height": height,
                "hpd_low": hpd_low,
                "hpd_high": hpd_high,
                "tip_names": list(node.tip_names),
            }
        )
    return annotations


def build_node_hpd_rows(source_id: str, annotations: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    """把 :func:`parse_mcc_node_annotations` 的结果转换为 ``BEAST_NODE_HPD_COLUMNS`` 行。"""
    rows: List[Dict[str, object]] = []
    for item in annotations:
        rows.append(
            {
                "source_id": source_id,
                "tip_hash": item["tip_hash"],
                "n_tips": item["n_tips"],
                "height_subst": "" if item["height"] is None else f"{float(item['height']):.12g}",
                "hpd_low_subst": "" if item["hpd_low"] is None else f"{float(item['hpd_low']):.12g}",
                "hpd_high_subst": "" if item["hpd_high"] is None else f"{float(item['hpd_high']):.12g}",
                "tip_names_json": json.dumps(item["tip_names"], ensure_ascii=False),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# 表 IO
# ---------------------------------------------------------------------------
def write_status_rows(rows: Sequence[Dict[str, object]], fieldnames: Sequence[str], destination: Path) -> None:
    """写出 TSV 状态/清单表。"""
    write_table(rows, fieldnames, destination)


def load_rows(path: Path) -> List[Dict[str, str]]:
    """读取 TSV 为字典列表。"""
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_beast_records(summary_path: Path):
    """读取 stage1 的子树汇总表（沿用 PAML 版 schema）。"""
    return load_baseml_summary(summary_path)


# ---------------------------------------------------------------------------
# CPU 亲和（与 paml_common.plan_cpu_slots 等价，用于并行作业绑核）
# ---------------------------------------------------------------------------
def _read_cpu_topology_key(cpu_id: int, key_name: str) -> Optional[int]:
    path = Path(f"/sys/devices/system/cpu/cpu{cpu_id}/topology/{key_name}")
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except (FileNotFoundError, ValueError):
        return None


def plan_cpu_slots(threads_per_job: int) -> List[Tuple[int, ...]]:
    """把可用 CPU 划分为若干「每作业 threads_per_job 核」的槽位，优先跨物理核。"""
    if int(threads_per_job) <= 0:
        raise PipelineError(f"threads_per_job must be positive: {threads_per_job}")

    available_cpus = sorted(os.sched_getaffinity(0))
    if len(available_cpus) < int(threads_per_job):
        raise PipelineError(
            f"Not enough CPUs for threads_per_job={threads_per_job}; available={len(available_cpus)}"
        )

    core_groups: Dict[Tuple[int, int], List[int]] = {}
    for cpu_id in available_cpus:
        package_id = _read_cpu_topology_key(cpu_id, "physical_package_id")
        core_id = _read_cpu_topology_key(cpu_id, "core_id")
        if package_id is None:
            package_id = 0
        if core_id is None:
            core_id = cpu_id
        core_groups.setdefault((package_id, core_id), []).append(cpu_id)

    ordered_groups = [sorted(group) for _, group in sorted(core_groups.items())]
    max_siblings = max(len(group) for group in ordered_groups)
    cpu_order: List[int] = []
    for sibling_index in range(max_siblings):
        for group in ordered_groups:
            if sibling_index < len(group):
                cpu_order.append(group[sibling_index])

    slots: List[Tuple[int, ...]] = []
    for index in range(0, len(cpu_order) - int(threads_per_job) + 1, int(threads_per_job)):
        slot = tuple(cpu_order[index : index + int(threads_per_job)])
        if len(slot) == int(threads_per_job):
            slots.append(slot)
    if not slots:
        raise PipelineError("Failed to derive any CPU slots from the current CPU affinity.")
    return slots


# ---------------------------------------------------------------------------
# 断点续跑追踪器（与 run_paml_jobs.ResumeTracker 等价：只按日志判完成）
# ---------------------------------------------------------------------------
class ResumeTracker:
    """基于 success.log / fail.log 的断点续跑状态追踪器（并行写入加文件锁）。

    完成判定**只**依据 success.log，绝不通过扫描输出目录或输出文件来推断完成状态。
    """

    def __init__(self, resume_log_dir: Path):
        self.log_dir = Path(resume_log_dir).resolve()
        self.lock_file = self.log_dir / ".resume.lock"
        self.success_log = self.log_dir / "success.log"
        self.fail_log = self.log_dir / "fail.log"
        self._state_lock = threading.Lock()
        self.completed_jobs: set = set()

    def initialize(self) -> None:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        for path in (self.lock_file, self.success_log, self.fail_log):
            path.touch(exist_ok=True)

        completed_jobs = set()
        with self.success_log.open("r", encoding="utf-8") as handle:
            for line in handle:
                job_id = line.split("\t", 1)[0].strip()
                if job_id:
                    completed_jobs.add(job_id)
        with self._state_lock:
            self.completed_jobs = completed_jobs

    def is_completed(self, job_id: str) -> bool:
        with self._state_lock:
            return job_id in self.completed_jobs

    def record_success(self, job_id: str) -> None:
        timestamp = dt.datetime.now().astimezone().isoformat(timespec="seconds")
        with self._state_lock:
            self._append_line(self.success_log, f"{job_id}\t{timestamp}\n")
            self.completed_jobs.add(job_id)

    def record_fail(self, job_id: str) -> None:
        with self._state_lock:
            self._append_line(self.fail_log, f"{job_id}\n")

    def _append_line(self, destination: Path, line: str) -> None:
        with self.lock_file.open("a", encoding="utf-8") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            try:
                with destination.open("a", encoding="utf-8") as destination_handle:
                    destination_handle.write(line)
                    destination_handle.flush()
            finally:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)

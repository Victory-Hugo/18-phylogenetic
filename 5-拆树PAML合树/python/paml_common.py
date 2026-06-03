#!/usr/bin/env python3
"""Utilities for preparing, running, and parsing baseml jobs."""

from __future__ import annotations

import csv
import io
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord

from phylo_merge_common import load_baseml_summary
from phylo_split_common import PipelineError, read_newick_tree, write_table, write_tree

PAML_JOB_MANIFEST_COLUMNS = [
    "baseml_subtree_id",
    "job_dir",
    "treefile",
    "seqfile",
    "ctlfile",
    "outfile",
    "expected_tip_count",
    "tip_hash",
    "seq_id_strategy",
]

PAML_RUN_STATUS_COLUMNS = [
    "baseml_subtree_id",
    "job_dir",
    "outfile",
    "return_code",
    "status",
    "detail",
]

PAML_PARSE_STATUS_COLUMNS = [
    "baseml_subtree_id",
    "result_file",
    "analysis_tree_file",
    "parameters_csv",
    "tip_length_csv",
    "status",
    "detail",
]


@dataclass
class ParsedBasemlOutput:
    branch_tags_line: str
    branch_params_line: str
    se_line: str
    named_tree_line: str
    indexed_tree_line: str


def resolve_seq_id_strategy(strategy: str) -> str:
    normalized = str(strategy).strip()
    if normalized not in {"exact", "prefix_before_first_underscore"}:
        raise PipelineError(f"Unsupported paml.seq_id_strategy: {strategy}")
    return normalized


def normalize_newick_text(raw_text: str) -> str:
    normalized = "".join(line.strip() for line in raw_text.splitlines() if line.strip())
    if not normalized.endswith(";"):
        normalized = f"{normalized};"
    return normalized


def read_fasta_record_map(fasta_path: Path) -> Dict[str, SeqRecord]:
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
    fasta_id_set = set(fasta_ids)
    mapping: Dict[str, str] = {}
    used_source_ids: Dict[str, str] = {}

    for tip_name in tip_names:
        if seq_id_strategy == "exact":
            source_id = tip_name
        else:
            # Prefer an exact FASTA id match when available. This keeps the
            # prefix strategy compatible with datasets whose tree tips and FASTA
            # headers already match in full, while still supporting prefix-only
            # FASTA ids such as JX289104_M33a3a -> JX289104.
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
    destination.parent.mkdir(parents=True, exist_ok=True)
    records: List[SeqRecord] = []
    for tip_name in tip_names:
        source_id = tip_to_source_id[tip_name]
        source_record = fasta_record_map[source_id]
        records.append(
            SeqRecord(
                seq=source_record.seq,
                id=tip_name,
                name=tip_name,
                description="",
            )
        )
    SeqIO.write(records, str(destination), "fasta")


def write_paml_treefile(
    source_tree_path: Path,
    expected_tip_count: int,
    destination: Path,
) -> None:
    tree = read_newick_tree(source_tree_path)
    actual_tip_count = len(tree.get_terminals())
    if actual_tip_count != int(expected_tip_count):
        raise PipelineError(
            f"Tree tip count mismatch for {source_tree_path}: expected {expected_tip_count}, got {actual_tip_count}"
        )
    newick_text = normalize_newick_text(source_tree_path.read_text(encoding="utf-8"))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(f"  {actual_tip_count} 1\n{newick_text}\n", encoding="utf-8")


def render_ctl_template(
    template_path: Path,
    seqfile_name: str,
    treefile_name: str,
    outfile_name: str,
    destination: Path,
) -> None:
    replacements = {
        "seqfile": seqfile_name,
        "treefile": treefile_name,
        "outfile": outfile_name,
    }
    lines = template_path.read_text(encoding="utf-8").splitlines()
    seen = set()
    rendered: List[str] = []
    for line in lines:
        stripped = line.lstrip()
        matched = False
        for key, value in replacements.items():
            if stripped.startswith(f"{key}"):
                prefix = line[: len(line) - len(stripped)]
                rendered.append(f"{prefix}{key} = {value}")
                seen.add(key)
                matched = True
                break
        if not matched:
            rendered.append(line)
    for key, value in replacements.items():
        if key not in seen:
            rendered.append(f"{key} = {value}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def write_rows(rows: Sequence[Dict[str, object]], fieldnames: Sequence[str], destination: Path) -> None:
    write_table(rows, fieldnames, destination)


def load_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def next_nonempty_index(lines: Sequence[str], start_index: int) -> int:
    for idx in range(start_index, len(lines)):
        if lines[idx].strip():
            return idx
    raise PipelineError("Unexpected end of baseml output while searching for the next non-empty line.")


def collect_nonempty_lines(
    lines: Sequence[str],
    start_index: int,
    end_index: Optional[int] = None,
) -> List[str]:
    stop_index = len(lines) if end_index is None else min(len(lines), end_index)
    return [line.strip() for line in lines[start_index:stop_index] if line.strip()]


def _parse_calibrated_rate(lines: Sequence[str]) -> Optional[float]:
    for idx, line in enumerate(lines):
        if "Substitution rate is per time unit" not in line:
            continue
        for next_line in lines[idx + 1 :]:
            stripped = next_line.strip()
            if not stripped:
                continue
            try:
                return float(stripped.split()[0])
            except (ValueError, IndexError):
                return None
    return None


def _parse_calibrated_node_times(text: str) -> Dict[int, float]:
    return {
        int(node_id): float(time_value)
        for node_id, time_value in re.findall(r"Node\s+(\d+)\s+Time\s+([0-9eE.+-]+)", text)
    }


def _clade_node_id(clade) -> Optional[int]:
    if clade.is_terminal() and clade.name:
        match = re.match(r"^(\d+)_", str(clade.name))
        if match:
            return int(match.group(1))
    if clade.confidence is not None:
        return int(clade.confidence)
    if clade.name and str(clade.name).isdigit():
        return int(str(clade.name))
    return None


def _strip_indexed_tip_name(name: str) -> str:
    match = re.match(r"^\d+_(.+)$", str(name))
    return match.group(1) if match else str(name)


def reconstruct_calibrated_named_tree(
    indexed_tree_line: str,
    node_times: Dict[int, float],
    substitution_rate: float,
    result_file: Path,
) -> str:
    """将 calibrated baseml 的节点年龄树重建为 substitutions/site 分支长度树。"""
    if substitution_rate <= 0:
        raise PipelineError(f"Invalid calibrated baseml substitution rate in {result_file}: {substitution_rate}")
    if not node_times:
        raise PipelineError(f"Could not parse calibrated baseml node times from {result_file}")

    tree = Phylo.read(io.StringIO(indexed_tree_line + "\n"), "newick")

    def assign_lengths(clade, parent_id: Optional[int], parent_time: Optional[float]) -> None:
        clade_id = _clade_node_id(clade)
        if clade.is_terminal():
            clade_time = 0.0
            clade.name = _strip_indexed_tip_name(clade.name)
        else:
            if clade_id not in node_times:
                raise PipelineError(f"Calibrated baseml node {clade_id} lacks time in {result_file}")
            clade_time = node_times[clade_id]
            clade.name = None
            clade.confidence = None

        if parent_id is None or parent_time is None:
            clade.branch_length = 0.0
        else:
            delta_time = parent_time - clade_time
            if delta_time < -1e-10:
                raise PipelineError(
                    f"Calibrated baseml has child time greater than parent time in {result_file}: "
                    f"parent={parent_id} child={clade_id}"
                )
            clade.branch_length = max(0.0, delta_time) * substitution_rate

        for child in clade.clades:
            assign_lengths(child, clade_id, clade_time)

    assign_lengths(tree.root, None, None)
    buffer = io.StringIO()
    Phylo.write(tree, buffer, "newick", format_branch_length="%1.12g")
    return normalize_newick_text(buffer.getvalue())


def parse_baseml_output(result_file: Path) -> ParsedBasemlOutput:
    lines = result_file.read_text(encoding="utf-8").splitlines()

    lnl_candidates = [
        idx
        for idx, line in enumerate(lines)
        if re.search(r"^\s*lnL(?:\(|\s*=)", line)
    ]
    lnl_index = lnl_candidates[-1] if lnl_candidates else None
    if lnl_index is None:
        raise PipelineError(f"Could not locate lnL section in baseml result: {result_file}")

    tree_anchor_index = next((idx for idx, line in enumerate(lines) if "tree length" in line), None)
    if tree_anchor_index is None:
        raise PipelineError(f"Could not locate tree length section in baseml result: {result_file}")
    if tree_anchor_index <= lnl_index:
        raise PipelineError(f"Invalid baseml output ordering around lnL/tree length: {result_file}")

    se_anchor_index = next(
        (
            idx
            for idx, line in enumerate(lines)
            if lnl_index < idx < tree_anchor_index and "SEs for parameters" in line
        ),
        None,
    )

    parameter_lines = collect_nonempty_lines(
        lines,
        lnl_index + 1,
        se_anchor_index if se_anchor_index is not None else tree_anchor_index,
    )
    if not parameter_lines:
        raise PipelineError(f"Could not locate branch parameter section in baseml result: {result_file}")

    first_parameter_tokens = parameter_lines[0].split()
    has_branch_tags = bool(first_parameter_tokens) and all(".." in token for token in first_parameter_tokens)
    if has_branch_tags:
        branch_tags_line = parameter_lines[0]
        if len(parameter_lines) < 2:
            raise PipelineError(f"Could not locate branch parameter values in baseml result: {result_file}")
        branch_params_line = parameter_lines[1]
    else:
        branch_tags_line = ""
        branch_params_line = parameter_lines[0]

    se_line = ""
    if se_anchor_index is not None:
        se_lines = collect_nonempty_lines(lines, se_anchor_index + 1, tree_anchor_index)
        if not se_lines:
            raise PipelineError(f"Could not locate SE values after SE header in baseml result: {result_file}")
        se_line = se_lines[0]

    named_tree_line = None
    indexed_tree_line = None
    found_named = False
    for idx in range(tree_anchor_index + 1, len(lines)):
        line = lines[idx].strip()
        if not line:
            continue
        if not found_named and ":" in line and line.endswith(";"):
            named_tree_line = line
            found_named = True
            continue
        if found_named and re.search(r"\b\d+_[^,\)\s:;]+", line) and line.endswith(";"):
            indexed_tree_line = line
            break

    if named_tree_line is None or indexed_tree_line is None:
        raise PipelineError(f"Could not extract named/indexed trees from baseml result: {result_file}")

    if any("Tree has ages" in line for line in lines):
        substitution_rate = _parse_calibrated_rate(lines)
        if substitution_rate is None:
            raise PipelineError(f"Could not parse calibrated baseml substitution rate from {result_file}")
        node_times = _parse_calibrated_node_times("\n".join(lines))
        named_tree_line = reconstruct_calibrated_named_tree(
            indexed_tree_line=indexed_tree_line,
            node_times=node_times,
            substitution_rate=substitution_rate,
            result_file=result_file,
        )

    return ParsedBasemlOutput(
        branch_tags_line=branch_tags_line,
        branch_params_line=branch_params_line,
        se_line=se_line,
        named_tree_line=named_tree_line,
        indexed_tree_line=indexed_tree_line,
    )


def baseml_output_looks_complete(result_file: Path) -> bool:
    if not result_file.exists() or result_file.stat().st_size == 0:
        return False
    try:
        parse_baseml_output(result_file)
    except PipelineError:
        return False
    return True


def parse_indexed_tip_mapping(indexed_tree_line: str) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for node_str, tip_name in re.findall(r"(\d+)_([^,\)\s:;]+)", indexed_tree_line):
        node_id = int(node_str)
        if node_id in mapping and mapping[node_id] != tip_name:
            raise PipelineError(f"Conflicting indexed tip mapping for node {node_id}")
        mapping[node_id] = tip_name
    if not mapping:
        raise PipelineError("No indexed tip mapping could be parsed from baseml output.")
    return mapping


def write_analysis_tree(named_tree_line: str, destination: Path) -> None:
    tree = Phylo.read(io.StringIO(named_tree_line + "\n"), "newick")
    write_tree(tree, destination)


def extract_tip_length_rows(named_tree_line: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for tip_name, length in re.findall(r"([^,:()\s]+):\s*([0-9eE.\-]+)", named_tree_line):
        rows.append({"ID": tip_name, "Length": length})
    if not rows:
        raise PipelineError("No tip lengths could be extracted from baseml named tree.")
    return rows


def build_parameter_rows(parsed_output: ParsedBasemlOutput) -> List[Dict[str, object]]:
    branch_tags = parsed_output.branch_tags_line.split()
    branch_params = parsed_output.branch_params_line.split()
    se_values = parsed_output.se_line.split()
    if not branch_tags or any(".." not in branch_tag for branch_tag in branch_tags):
        n_rows = len(branch_params)
        if n_rows == 0:
            raise PipelineError("No parameter rows could be parsed from baseml output.")
        rows: List[Dict[str, object]] = []
        for idx, branch_param in enumerate(branch_params[:n_rows], start=1):
            rows.append(
                {
                    "Branch": f"PARAM_{idx}",
                    "Param": branch_param,
                    "SE": se_values[idx - 1] if idx - 1 < len(se_values) else "",
                    "Parent": "",
                    "Child": "",
                    "Parent_ID": "",
                    "Child_ID": "",
                }
            )
        return rows

    indexed_tip_map = parse_indexed_tip_mapping(parsed_output.indexed_tree_line)
    n_rows = min(len(branch_tags), len(branch_params))
    if n_rows == 0:
        raise PipelineError("No branch parameter rows could be parsed from baseml output.")

    rows: List[Dict[str, object]] = []
    for idx, (branch_tag, branch_param) in enumerate(zip(branch_tags[:n_rows], branch_params[:n_rows])):
        parent_str, child_str = branch_tag.split("..", 1)
        parent_id = int(parent_str)
        child_id = int(child_str)
        rows.append(
            {
                "Branch": branch_tag,
                "Param": branch_param,
                "SE": se_values[idx] if idx < len(se_values) else "",
                "Parent": parent_id,
                "Child": child_id,
                "Parent_ID": indexed_tip_map.get(parent_id, ""),
                "Child_ID": indexed_tip_map.get(child_id, ""),
            }
        )
    return rows

def parse_cpu_list(cpu_list: str) -> Tuple[int, ...]:
    values: List[int] = []
    for token in str(cpu_list).split(","):
        token = token.strip()
        if not token:
            continue
        values.append(int(token))
    if not values:
        raise PipelineError(f"CPU list is empty: {cpu_list}")
    return tuple(values)


def _read_cpu_topology_key(cpu_id: int, key_name: str) -> Optional[int]:
    path = Path(f"/sys/devices/system/cpu/cpu{cpu_id}/topology/{key_name}")
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except (FileNotFoundError, ValueError):
        return None


def plan_cpu_slots(threads_per_job: int) -> List[Tuple[int, ...]]:
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


def run_baseml_job(
    baseml_bin: str,
    ctlfile: Path,
    job_dir: Path,
    cpu_ids: Optional[Sequence[int]] = None,
) -> subprocess.CompletedProcess:
    child_affinity = tuple(int(cpu_id) for cpu_id in cpu_ids) if cpu_ids else None

    def _apply_affinity():
        if child_affinity is not None:
            os.sched_setaffinity(0, set(child_affinity))

    return subprocess.run(
        [baseml_bin, ctlfile.name],
        cwd=str(job_dir),
        check=False,
        text=True,
        capture_output=True,
        preexec_fn=_apply_affinity if child_affinity is not None else None,
    )


def load_baseml_records(summary_path: Path):
    return load_baseml_summary(summary_path)

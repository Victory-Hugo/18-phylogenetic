#!/usr/bin/env python3
"""Utilities for preparing, running, and parsing baseml jobs."""

from __future__ import annotations

import csv
import io
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

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
            source_id = tip_name.split("_", 1)[0]
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


def write_paml_treefile(source_tree_path: Path, expected_tip_count: int, destination: Path) -> None:
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
    branch_tags_index = next_nonempty_index(lines, lnl_index + 1)
    branch_params_index = next_nonempty_index(lines, branch_tags_index + 1)
    se_anchor_index = next((idx for idx, line in enumerate(lines) if "SEs for parameters" in line), None)
    if se_anchor_index is None:
        raise PipelineError(f"Could not locate SEs section in baseml result: {result_file}")
    if branch_params_index == se_anchor_index:
        branch_tags_line = ""
        branch_params_line = lines[branch_tags_index].strip()
    else:
        branch_tags_line = lines[branch_tags_index].strip()
        branch_params_line = lines[branch_params_index].strip()
    se_line_index = next_nonempty_index(lines, se_anchor_index + 1)

    tree_anchor_index = next((idx for idx, line in enumerate(lines) if "tree length" in line), None)
    if tree_anchor_index is None:
        raise PipelineError(f"Could not locate tree length section in baseml result: {result_file}")

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

    return ParsedBasemlOutput(
        branch_tags_line=branch_tags_line,
        branch_params_line=branch_params_line,
        se_line=lines[se_line_index].strip(),
        named_tree_line=named_tree_line,
        indexed_tree_line=indexed_tree_line,
    )


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
        n_rows = min(len(branch_params), len(se_values))
        if n_rows == 0:
            raise PipelineError("No parameter rows could be parsed from baseml output.")
        rows: List[Dict[str, object]] = []
        for idx, (branch_param, se_value) in enumerate(zip(branch_params[:n_rows], se_values[:n_rows]), start=1):
            rows.append(
                {
                    "Branch": f"PARAM_{idx}",
                    "Param": branch_param,
                    "SE": se_value,
                    "Parent": "",
                    "Child": "",
                    "Parent_ID": "",
                    "Child_ID": "",
                }
            )
        return rows

    indexed_tip_map = parse_indexed_tip_mapping(parsed_output.indexed_tree_line)
    n_rows = min(len(branch_tags), len(branch_params), len(se_values))
    if n_rows == 0:
        raise PipelineError("No branch parameter rows could be parsed from baseml output.")

    rows: List[Dict[str, object]] = []
    for branch_tag, branch_param, se_value in zip(branch_tags[:n_rows], branch_params[:n_rows], se_values[:n_rows]):
        parent_str, child_str = branch_tag.split("..", 1)
        parent_id = int(parent_str)
        child_id = int(child_str)
        rows.append(
            {
                "Branch": branch_tag,
                "Param": branch_param,
                "SE": se_value,
                "Parent": parent_id,
                "Child": child_id,
                "Parent_ID": indexed_tip_map.get(parent_id, ""),
                "Child_ID": indexed_tip_map.get(child_id, ""),
            }
        )
    return rows


def run_baseml_job(baseml_bin: str, ctlfile: Path, job_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [baseml_bin, ctlfile.name],
        cwd=str(job_dir),
        check=False,
        text=True,
        capture_output=True,
    )


def load_baseml_records(summary_path: Path):
    return load_baseml_summary(summary_path)

#!/usr/bin/env python3
"""Shared utilities for phylogenetic core partition and baseml subtree construction.

Dependencies:
  pip install biopython
  conda run -n BigLin gotree --help
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import logging
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree


class PipelineError(RuntimeError):
    """Raised when the pipeline encounters a recoverable domain error."""


GLOBAL_OUTGROUP_TIP = "RSRS"

CORE_SUMMARY_COLUMNS = [
    "core_subtree_id",
    "core_root_node_id",
    "parent_node_id",
    "core_n_tips",
    "core_tip_names",
    "tip_hash",
    "core_tree_file",
]

BASEML_SUMMARY_COLUMNS = [
    "baseml_subtree_id",
    "core_subtree_id",
    "core_root_node_id",
    "outgroup_tip",
    "core_n_tips",
    "anchor_n_tips",
    "total_n_tips",
    "core_tip_names",
    "anchor_tip_names",
    "total_tip_names",
    "anchor_source_node_ids",
    "tip_hash",
    "baseml_tree_file",
]

MANIFEST_COLUMNS = [
    "baseml_subtree_id",
    "tip_name",
    "tip_role",
    "source_core_subtree_id",
    "source_node_id",
]

OVERLAP_COLUMNS = [
    "tip_name",
    "owner_core_subtree_id",
    "appearing_baseml_subtree_ids",
    "roles",
]

VALIDATION_COLUMNS = ["check_name", "status", "details"]


@dataclass
class CoreSubtreeRecord:
    core_subtree_id: str
    core_root_node_id: str
    parent_node_id: str
    core_n_tips: int
    core_tip_names: List[str]
    tip_hash: str
    core_tree_file: str
    clade: Clade

    def to_row(self) -> Dict[str, object]:
        return {
            "core_subtree_id": self.core_subtree_id,
            "core_root_node_id": self.core_root_node_id,
            "parent_node_id": self.parent_node_id,
            "core_n_tips": self.core_n_tips,
            "core_tip_names": encode_json_list(self.core_tip_names),
            "tip_hash": self.tip_hash,
            "core_tree_file": self.core_tree_file,
        }


@dataclass
class AnchorUnit:
    source_node_id: str
    tip_names: List[str]
    n_tips: int
    ancestor_distance: int
    unit_type: str


@dataclass
class BasemlSubtreeRecord:
    baseml_subtree_id: str
    core_subtree_id: str
    core_root_node_id: str
    outgroup_tip: str
    core_n_tips: int
    anchor_n_tips: int
    total_n_tips: int
    core_tip_names: List[str]
    anchor_tip_names: List[str]
    total_tip_names: List[str]
    anchor_source_node_ids: List[str]
    tip_hash: str
    baseml_tree_file: str

    def to_row(self) -> Dict[str, object]:
        return {
            "baseml_subtree_id": self.baseml_subtree_id,
            "core_subtree_id": self.core_subtree_id,
            "core_root_node_id": self.core_root_node_id,
            "outgroup_tip": self.outgroup_tip,
            "core_n_tips": self.core_n_tips,
            "anchor_n_tips": self.anchor_n_tips,
            "total_n_tips": self.total_n_tips,
            "core_tip_names": encode_json_list(self.core_tip_names),
            "anchor_tip_names": encode_json_list(self.anchor_tip_names),
            "total_tip_names": encode_json_list(self.total_tip_names),
            "anchor_source_node_ids": encode_json_list(self.anchor_source_node_ids),
            "tip_hash": self.tip_hash,
            "baseml_tree_file": self.baseml_tree_file,
        }


def setup_logger(name: str, log_path: Optional[Path], level: str) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, str(level).upper(), logging.INFO))
    logger.propagate = False
    while logger.handlers:
        handler = logger.handlers.pop()
        handler.close()

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def read_newick_tree(tree_path: Path) -> Tree:
    try:
        return Phylo.read(str(tree_path), "newick")
    except Exception as exc:  # pragma: no cover
        raise PipelineError(f"Failed to read Newick tree: {tree_path}") from exc


def write_tree(tree: Tree, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tree.rooted = True
    Phylo.write(tree, str(destination), "newick")


def get_tip_names_from_tree(tree: Tree) -> List[str]:
    names = []
    for tip in tree.get_terminals():
        if not tip.name:
            raise PipelineError("Found unnamed terminal tip in tree.")
        names.append(tip.name)
    return names


def parse_tip_file(path: Path) -> List[str]:
    if not path.exists():
        raise PipelineError(f"Tip file not found: {path}")
    tips: List[str] = []
    seen = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line in seen:
                raise PipelineError(f"Duplicate tip in file {path}: {line}")
            seen.add(line)
            tips.append(line)
    if not tips:
        raise PipelineError(f"Tip file is empty after filtering comments/blanks: {path}")
    return tips


def validate_outgroup_tips(outgroup_tips: Sequence[str], tree_tip_names: Iterable[str]) -> None:
    tree_tip_set = set(tree_tip_names)
    missing = [tip for tip in outgroup_tips if tip not in tree_tip_set]
    if missing:
        preview = ", ".join(missing[:10])
        raise PipelineError(f"Outgroup tips not found in input tree ({len(missing)} missing): {preview}")


def load_backbone_tips(backbone_tree: Optional[Path], tree_tip_names: Iterable[str]) -> List[str]:
    if backbone_tree is None:
        return []
    tree_tip_set = set(tree_tip_names)
    backbone = read_newick_tree(backbone_tree)
    backbone_tips = get_tip_names_from_tree(backbone)
    missing = [tip for tip in backbone_tips if tip not in tree_tip_set]
    if missing:
        preview = ", ".join(missing[:10])
        raise PipelineError(f"Backbone tips not found in input tree ({len(missing)} missing): {preview}")
    return backbone_tips


def run_command(command: Sequence[str], logger: logging.Logger) -> subprocess.CompletedProcess:
    logger.info("Running command: %s", " ".join(shlex.quote(part) for part in command))
    completed = subprocess.run(list(command), check=False, text=True, capture_output=True)
    if completed.stdout:
        logger.info("Command stdout:\n%s", completed.stdout.rstrip())
    if completed.stderr:
        logger.info("Command stderr:\n%s", completed.stderr.rstrip())
    if completed.returncode != 0:
        raise PipelineError(
            f"Command failed with exit code {completed.returncode}: "
            f"{' '.join(shlex.quote(part) for part in command)}"
        )
    return completed


def run_gotree_stats(
    input_tree: Path,
    conda_env: str,
    gotree_bin: str,
    threads: int,
    logger: logging.Logger,
) -> str:
    command = [
        "conda",
        "run",
        "-n",
        conda_env,
        gotree_bin,
        "stats",
        "-i",
        str(input_tree),
        "-t",
        str(threads),
    ]
    completed = run_command(command, logger)
    return completed.stdout.strip()


def parse_rooted_flag(stats_output: str) -> str:
    lines = [line.strip() for line in stats_output.splitlines() if line.strip()]
    if len(lines) < 2:
        raise PipelineError("Unexpected gotree stats output; could not determine rooted status.")
    headers = lines[0].split("\t")
    values = lines[1].split("\t")
    if "rooted" not in headers:
        raise PipelineError("gotree stats output does not contain a rooted column.")
    rooted_index = headers.index("rooted")
    if rooted_index >= len(values):
        raise PipelineError("gotree stats rooted column is missing a value.")
    rooted_value = values[rooted_index].strip().lower()
    if rooted_value not in {"rooted", "unrooted"}:
        raise PipelineError(f"Unexpected rooted value from gotree stats: {rooted_value}")
    return rooted_value


def detect_tree_rooted_status(
    input_tree: Path,
    conda_env: str,
    gotree_bin: str,
    threads: int,
    logger: logging.Logger,
) -> Tuple[str, str]:
    stats_output = run_gotree_stats(input_tree, conda_env, gotree_bin, threads, logger)
    rooted_status = parse_rooted_flag(stats_output)
    logger.info("gotree stats summary: %s", stats_output.replace("\n", " | "))
    return rooted_status, stats_output


def reroot_with_outgroup(
    input_tree: Path,
    outgroup_tip_file: Path,
    rooted_tree: Path,
    conda_env: str,
    gotree_bin: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    rooted_tree.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "conda",
        "run",
        "-n",
        conda_env,
        gotree_bin,
        "reroot",
        "outgroup",
        "--strict",
        "-i",
        str(input_tree),
        "-l",
        str(outgroup_tip_file),
        "-o",
        str(rooted_tree),
        "-t",
        str(threads),
    ]
    run_command(command, logger)


def prepare_rooted_tree(
    input_tree: Path,
    rooted_tree: Path,
    conda_env: str,
    gotree_bin: str,
    threads: int,
    logger: logging.Logger,
    outgroup_tip_file: Optional[Path] = None,
    rooted_status: Optional[str] = None,
) -> str:
    if rooted_status is None:
        rooted_status, _ = detect_tree_rooted_status(input_tree, conda_env, gotree_bin, threads, logger)
    rooted_tree.parent.mkdir(parents=True, exist_ok=True)
    if rooted_status == "rooted":
        shutil.copy2(input_tree, rooted_tree)
        logger.info("Input tree is already rooted; copied directly to %s", rooted_tree)
        return rooted_status
    if outgroup_tip_file is None:
        raise PipelineError("Input tree is unrooted and no outgroup tip file was provided.")
    reroot_with_outgroup(input_tree, outgroup_tip_file, rooted_tree, conda_env, gotree_bin, threads, logger)
    return rooted_status


def assign_node_ids(tree: Tree) -> Tuple[Dict[Clade, str], Dict[str, Clade], Dict[Clade, Optional[Clade]]]:
    node_id_map: Dict[Clade, str] = {}
    reverse_node_id_map: Dict[str, Clade] = {}
    parent_map: Dict[Clade, Optional[Clade]] = {tree.root: None}
    internal_counter = 1
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            node_id = f"TIP::{clade.name}"
        else:
            node_id = f"NODE_{internal_counter:06d}"
            internal_counter += 1
        node_id_map[clade] = node_id
        reverse_node_id_map[node_id] = clade
        for child in clade.clades:
            parent_map[child] = clade
    return node_id_map, reverse_node_id_map, parent_map


def compute_tip_counts(tree: Tree) -> Dict[Clade, int]:
    tip_counts: Dict[Clade, int] = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            tip_counts[clade] = 1
        else:
            tip_counts[clade] = sum(tip_counts[child] for child in clade.clades)
    return tip_counts


def collect_selected_roots_from_nodes(
    start_nodes: Sequence[Clade],
    tip_counts: Dict[Clade, int],
    max_tips: int,
) -> List[Clade]:
    selected: List[Clade] = []

    def recurse(clade: Clade) -> None:
        if tip_counts[clade] <= max_tips:
            selected.append(clade)
            return
        for child in clade.clades:
            recurse(child)

    for node in start_nodes:
        recurse(node)
    return selected


def collect_selected_roots(tree: Tree, tip_counts: Dict[Clade, int], max_tips: int) -> List[Clade]:
    if tip_counts[tree.root] <= max_tips:
        return [tree.root]
    return collect_selected_roots_from_nodes(tree.root.clades, tip_counts, max_tips)


def is_descendant(node: Clade, ancestor: Clade, parent_map: Dict[Clade, Optional[Clade]]) -> bool:
    current = node
    while current is not None:
        if current is ancestor:
            return True
        current = parent_map[current]
    return False


def merge_small_blocks(
    tree: Tree,
    selected_roots: Sequence[Clade],
    tip_counts: Dict[Clade, int],
    parent_map: Dict[Clade, Optional[Clade]],
    max_tips: int,
    logger: logging.Logger,
) -> List[Clade]:
    selected = set(selected_roots)
    if not selected:
        return []

    changed = True
    while changed:
        changed = False
        for clade in tree.find_clades(order="postorder"):
            if clade in selected or tip_counts[clade] > max_tips:
                continue
            descendants = [root for root in selected if root is not clade and is_descendant(root, clade, parent_map)]
            if not descendants:
                continue
            covered = sum(tip_counts[root] for root in descendants)
            if covered != tip_counts[clade]:
                continue
            for descendant in descendants:
                selected.remove(descendant)
            selected.add(clade)
            changed = True
            logger.info(
                "Merged %d descendant block(s) into %s (%d tips).",
                len(descendants),
                clade.name or "<internal>",
                tip_counts[clade],
            )
    return list(selected)


def get_clade_tip_names(clade: Clade) -> List[str]:
    return [tip.name for tip in clade.get_terminals()]


def compute_tip_hash(tip_names: Sequence[str]) -> str:
    payload = "\n".join(sorted(tip_names)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def encode_json_list(items: Sequence[str]) -> str:
    return json.dumps(list(items), ensure_ascii=False)


def decode_json_list(value: str) -> List[str]:
    parsed = json.loads(value)
    if not isinstance(parsed, list):
        raise PipelineError("Expected a JSON list in table.")
    return [str(item) for item in parsed]


def subtree_sort_key(clade: Clade, node_id_map: Dict[Clade, str]) -> Tuple[int, str]:
    node_id = node_id_map[clade]
    if node_id.startswith("NODE_"):
        return (0, node_id)
    return (1, node_id)


def validate_selected_clades(
    selected_roots: Sequence[Clade],
    tip_counts: Dict[Clade, int],
    max_tips: int,
) -> Tuple[set, set]:
    all_seen = set()
    duplicates = set()
    for root in selected_roots:
        if tip_counts[root] > max_tips:
            raise PipelineError(f"Selected subtree exceeds max_tips: {tip_counts[root]} > {max_tips}")
        for tip_name in get_clade_tip_names(root):
            if tip_name in all_seen:
                duplicates.add(tip_name)
            all_seen.add(tip_name)
    return all_seen, duplicates


def write_clade_tree(clade: Clade, destination: Path) -> None:
    subtree = Tree(root=copy.deepcopy(clade), rooted=True)
    write_tree(subtree, destination)


def write_table(rows: Sequence[Dict[str, object]], fieldnames: Sequence[str], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_table(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise PipelineError(f"Required output table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_root_list(root_ids: Sequence[str], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        for root_id in root_ids:
            handle.write(f"{root_id}\n")


def load_root_list(path: Path) -> List[str]:
    if not path.exists():
        raise PipelineError(f"Required root list not found: {path}")
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def write_validation_report(rows: Sequence[Tuple[str, str, str]], destination: Path) -> None:
    write_table(
        [{"check_name": a, "status": b, "details": c} for a, b, c in rows],
        VALIDATION_COLUMNS,
        destination,
    )


def clean_generated_outputs(output_dir: Path) -> None:
    for pattern in [
        "subtrees/subtree_*.nwk",
        "core_subtrees/core_subtree_*.nwk",
        "baseml_subtrees/baseml_subtree_*.nwk",
    ]:
        for candidate in output_dir.glob(pattern):
            candidate.unlink(missing_ok=True)
    for relative_path in [
        "subtree_summary.tsv",
        "subtree_roots.txt",
        "validation_report.tsv",
        "core_subtree_summary.tsv",
        "core_subtree_roots.txt",
        "baseml_subtree_summary.tsv",
        "baseml_tree_manifest.tsv",
        "baseml_validation_report.tsv",
        "overlap_report.tsv",
        "split_tree.log",
        "intermediate/rooted.tree",
        "intermediate/rooted.validation.tree",
    ]:
        (output_dir / relative_path).unlink(missing_ok=True)


def find_tip_clade(tree: Tree, tip_name: str) -> Clade:
    for clade in tree.get_terminals():
        if clade.name == tip_name:
            return clade
    raise PipelineError(f"Tip not found in tree: {tip_name}")


def get_root_children_for_outgroup(tree: Tree, outgroup_tip: str) -> Tuple[Clade, List[Clade]]:
    if len(tree.root.clades) != 2:
        raise PipelineError("Rooted input tree must have exactly 2 root children for baseml.")
    outgroup_children = [child for child in tree.root.clades if outgroup_tip in set(get_clade_tip_names(child))]
    if len(outgroup_children) != 1:
        raise PipelineError(f"Could not identify unique root child containing outgroup tip {outgroup_tip}.")
    outgroup_child = outgroup_children[0]
    if set(get_clade_tip_names(outgroup_child)) != {outgroup_tip}:
        raise PipelineError(
            f"Root child containing outgroup tip {outgroup_tip} must be a singleton clade for this workflow."
        )
    ingroup_children = [child for child in tree.root.clades if child is not outgroup_child]
    if not ingroup_children:
        raise PipelineError("No ingroup child found after identifying the outgroup root child.")
    return outgroup_child, ingroup_children


def build_core_partition(
    tree: Tree,
    node_id_map: Dict[Clade, str],
    parent_map: Dict[Clade, Optional[Clade]],
    tip_counts: Dict[Clade, int],
    effective_core_max_tips: int,
    outgroup_tip: str,
    enable_merge_small_blocks: bool,
    logger: logging.Logger,
) -> List[CoreSubtreeRecord]:
    _, ingroup_children = get_root_children_for_outgroup(tree, outgroup_tip)
    selected_roots = collect_selected_roots_from_nodes(ingroup_children, tip_counts, effective_core_max_tips)
    before_merge = len(selected_roots)
    if enable_merge_small_blocks:
        selected_roots = merge_small_blocks(
            tree=tree,
            selected_roots=selected_roots,
            tip_counts=tip_counts,
            parent_map=parent_map,
            max_tips=effective_core_max_tips,
            logger=logger,
        )
    logger.info("Selected core subtree count before merge: %d", before_merge)
    logger.info("Selected core subtree count after merge: %d", len(selected_roots))

    ingroup_tip_names = set()
    for node in ingroup_children:
        ingroup_tip_names.update(get_clade_tip_names(node))
    seen_tips, duplicates = validate_selected_clades(selected_roots, tip_counts, effective_core_max_tips)
    if duplicates:
        preview = ", ".join(sorted(duplicates)[:10])
        raise PipelineError(f"Duplicate tip assignments detected in core partition: {preview}")
    missing = sorted(ingroup_tip_names - seen_tips)
    if missing:
        preview = ", ".join(missing[:10])
        raise PipelineError(f"Missing tips in core partition ({len(missing)} missing): {preview}")

    records: List[CoreSubtreeRecord] = []
    for index, clade in enumerate(sorted(selected_roots, key=lambda c: subtree_sort_key(c, node_id_map)), start=1):
        tips = get_clade_tip_names(clade)
        parent = parent_map[clade]
        parent_node_id = "ROOT" if parent is None or parent is tree.root else node_id_map[parent]
        record = CoreSubtreeRecord(
            core_subtree_id=f"core_subtree_{index:04d}",
            core_root_node_id=node_id_map[clade],
            parent_node_id=parent_node_id,
            core_n_tips=len(tips),
            core_tip_names=tips,
            tip_hash=compute_tip_hash(tips),
            core_tree_file=Path("core_subtrees", f"core_subtree_{index:04d}.nwk").as_posix(),
            clade=clade,
        )
        records.append(record)
    return records


def build_tip_owner_map(core_records: Sequence[CoreSubtreeRecord]) -> Dict[str, str]:
    owner: Dict[str, str] = {}
    for record in core_records:
        for tip in record.core_tip_names:
            if tip in owner:
                raise PipelineError(f"Tip assigned to more than one core partition: {tip}")
            owner[tip] = record.core_subtree_id
    return owner


def build_ordered_tip_list(master_tip_order: Sequence[str], selected_tip_set: set) -> List[str]:
    return [tip for tip in master_tip_order if tip in selected_tip_set]


def collect_anchor_units_for_sister(
    sister: Clade,
    tip_counts: Dict[Clade, int],
    node_id_map: Dict[Clade, str],
    ancestor_distance: int,
    budget: int,
) -> List[AnchorUnit]:
    if budget <= 0:
        return []
    if tip_counts[sister] <= budget:
        tips = get_clade_tip_names(sister)
        return [
            AnchorUnit(
                source_node_id=node_id_map[sister],
                tip_names=tips,
                n_tips=len(tips),
                ancestor_distance=ancestor_distance,
                unit_type="whole_clade",
            )
        ]
    units = []
    selected = collect_selected_roots_from_nodes([sister], tip_counts, budget)
    for root in sorted(selected, key=lambda c: subtree_sort_key(c, node_id_map)):
        tips = get_clade_tip_names(root)
        units.append(
            AnchorUnit(
                source_node_id=node_id_map[root],
                tip_names=tips,
                n_tips=len(tips),
                ancestor_distance=ancestor_distance,
                unit_type="split_clade",
            )
        )
    return units


def gather_candidate_anchor_units(
    core_clade: Clade,
    parent_map: Dict[Clade, Optional[Clade]],
    tip_counts: Dict[Clade, int],
    node_id_map: Dict[Clade, str],
    outgroup_tip: str,
    budget: int,
) -> List[AnchorUnit]:
    units: List[AnchorUnit] = []
    current = core_clade
    ancestor_distance = 1
    while parent_map[current] is not None:
        parent = parent_map[current]
        for sibling in parent.clades:
            if sibling is current:
                continue
            sibling_tips = set(get_clade_tip_names(sibling))
            if outgroup_tip in sibling_tips:
                continue
            units.extend(
                collect_anchor_units_for_sister(
                    sister=sibling,
                    tip_counts=tip_counts,
                    node_id_map=node_id_map,
                    ancestor_distance=ancestor_distance,
                    budget=budget,
                )
            )
        current = parent
        ancestor_distance += 1
    units.sort(key=lambda unit: (unit.ancestor_distance, -unit.n_tips, unit.source_node_id))
    return units


def build_manifest_rows(
    record: BasemlSubtreeRecord,
    anchor_tip_source_node: Dict[str, str],
    tip_owner_map: Dict[str, str],
    outgroup_tip: str,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for tip in record.core_tip_names:
        rows.append(
            {
                "baseml_subtree_id": record.baseml_subtree_id,
                "tip_name": tip,
                "tip_role": "core",
                "source_core_subtree_id": record.core_subtree_id,
                "source_node_id": record.core_root_node_id,
            }
        )
    for tip in record.anchor_tip_names:
        rows.append(
            {
                "baseml_subtree_id": record.baseml_subtree_id,
                "tip_name": tip,
                "tip_role": "anchor",
                "source_core_subtree_id": tip_owner_map.get(tip, ""),
                "source_node_id": anchor_tip_source_node[tip],
            }
        )
    rows.append(
        {
            "baseml_subtree_id": record.baseml_subtree_id,
            "tip_name": outgroup_tip,
            "tip_role": "outgroup",
            "source_core_subtree_id": "",
            "source_node_id": f"TIP::{outgroup_tip}",
        }
    )
    return rows


def build_baseml_subtrees(
    tree: Tree,
    core_records: Sequence[CoreSubtreeRecord],
    parent_map: Dict[Clade, Optional[Clade]],
    tip_counts: Dict[Clade, int],
    node_id_map: Dict[Clade, str],
    outgroup_tip: str,
    min_baseml_tips: int,
    max_tips: int,
    logger: logging.Logger,
) -> Tuple[List[BasemlSubtreeRecord], List[Dict[str, object]]]:
    master_tip_order = get_tip_names_from_tree(tree)
    input_smaller_than_min = len(master_tip_order) < min_baseml_tips
    tip_owner_map = build_tip_owner_map(core_records)
    records: List[BasemlSubtreeRecord] = []
    manifest_rows: List[Dict[str, object]] = []

    for index, core_record in enumerate(core_records, start=1):
        selected_tip_set = set(core_record.core_tip_names)
        selected_tip_set.add(outgroup_tip)
        anchor_tip_source_node: Dict[str, str] = {}
        anchor_source_node_ids: List[str] = []

        if len(selected_tip_set) < min_baseml_tips:
            budget = max_tips - len(selected_tip_set)
            candidate_units = gather_candidate_anchor_units(
                core_clade=core_record.clade,
                parent_map=parent_map,
                tip_counts=tip_counts,
                node_id_map=node_id_map,
                outgroup_tip=outgroup_tip,
                budget=budget,
            )
            for unit in candidate_units:
                if len(selected_tip_set) >= min_baseml_tips:
                    break
                if unit.n_tips > max_tips - len(selected_tip_set):
                    continue
                new_tips = [tip for tip in unit.tip_names if tip not in selected_tip_set]
                if not new_tips:
                    continue
                if len(selected_tip_set) + len(new_tips) > max_tips:
                    continue
                for tip in new_tips:
                    selected_tip_set.add(tip)
                    anchor_tip_source_node[tip] = unit.source_node_id
                if unit.source_node_id not in anchor_source_node_ids:
                    anchor_source_node_ids.append(unit.source_node_id)

        if len(selected_tip_set) < min_baseml_tips and not input_smaller_than_min:
            raise PipelineError(
                f"Could not construct baseml subtree >= {min_baseml_tips} tips for {core_record.core_subtree_id}."
            )
        if len(selected_tip_set) > max_tips:
            raise PipelineError(f"Constructed baseml subtree exceeds max_tips for {core_record.core_subtree_id}.")

        total_tip_names = build_ordered_tip_list(master_tip_order, selected_tip_set)
        anchor_tip_names = [tip for tip in total_tip_names if tip not in set(core_record.core_tip_names) and tip != outgroup_tip]
        record = BasemlSubtreeRecord(
            baseml_subtree_id=f"baseml_subtree_{index:04d}",
            core_subtree_id=core_record.core_subtree_id,
            core_root_node_id=core_record.core_root_node_id,
            outgroup_tip=outgroup_tip,
            core_n_tips=core_record.core_n_tips,
            anchor_n_tips=len(anchor_tip_names),
            total_n_tips=len(total_tip_names),
            core_tip_names=list(core_record.core_tip_names),
            anchor_tip_names=anchor_tip_names,
            total_tip_names=total_tip_names,
            anchor_source_node_ids=anchor_source_node_ids,
            tip_hash=compute_tip_hash(total_tip_names),
            baseml_tree_file=Path("baseml_subtrees", f"baseml_subtree_{index:04d}.nwk").as_posix(),
        )
        records.append(record)
        manifest_rows.extend(build_manifest_rows(record, anchor_tip_source_node, tip_owner_map, outgroup_tip))
        logger.info(
            "%s core=%d anchor=%d total=%d trigger_anchor=%s source_nodes=%s",
            record.baseml_subtree_id,
            record.core_n_tips,
            record.anchor_n_tips,
            record.total_n_tips,
            "yes" if record.anchor_n_tips else "no",
            ",".join(record.anchor_source_node_ids) if record.anchor_source_node_ids else "-",
        )
    return records, manifest_rows


def normalize_tree_binary(tree: Tree) -> Tree:
    def collapse(clade: Clade) -> Clade:
        if clade.is_terminal():
            return clade
        clade.clades = [collapse(child) for child in clade.clades]
        while len(clade.clades) == 1 and not clade.is_terminal():
            child = clade.clades[0]
            if child.branch_length is None:
                child.branch_length = clade.branch_length
            elif clade.branch_length is not None:
                child.branch_length += clade.branch_length
            if clade.name and not child.name:
                child.name = clade.name
            clade = child
            if not clade.is_terminal():
                clade.clades = [collapse(grandchild) for grandchild in clade.clades]
        return clade

    tree.root = collapse(tree.root)
    tree.root.branch_length = None
    tree.rooted = True
    return tree


def is_binary_rooted_with_outgroup(tree: Tree, outgroup_tip: str) -> Tuple[bool, str]:
    if len(tree.root.clades) != 2:
        return False, "Root does not have exactly 2 children."
    root_child_tip_sets = [set(get_clade_tip_names(child)) for child in tree.root.clades]
    if {outgroup_tip} not in root_child_tip_sets:
        return False, f"Root does not contain singleton outgroup child {outgroup_tip}."
    for clade in tree.get_nonterminals(order="preorder"):
        if len(clade.clades) != 2:
            return False, f"Internal node is not bifurcating: {clade.name or '<internal>'}"
    return True, "Tree is rooted and binary."


def build_induced_tree_file(
    rooted_master_tree: Path,
    selected_tip_names: Sequence[str],
    outgroup_tip: str,
    destination: Path,
    intermediate_dir: Path,
    conda_env: str,
    gotree_bin: str,
    threads: int,
    logger: logging.Logger,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=intermediate_dir, delete=False, suffix=".tips.txt") as handle:
        keep_tip_file = Path(handle.name)
        handle.write("\n".join(selected_tip_names))
        handle.write("\n")
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=intermediate_dir, delete=False, suffix=".outgroup.txt") as handle:
        rsrs_file = Path(handle.name)
        handle.write(f"{outgroup_tip}\n")
    pruned_tree = intermediate_dir / f"{destination.stem}.pruned.nwk"
    rerooted_tree = intermediate_dir / f"{destination.stem}.rerooted.nwk"
    try:
        run_command(
            [
                "conda",
                "run",
                "-n",
                conda_env,
                gotree_bin,
                "prune",
                "-r",
                "-i",
                str(rooted_master_tree),
                "-f",
                str(keep_tip_file),
                "-o",
                str(pruned_tree),
            ],
            logger,
        )
        run_command(
            [
                "conda",
                "run",
                "-n",
                conda_env,
                gotree_bin,
                "reroot",
                "outgroup",
                "--strict",
                "-i",
                str(pruned_tree),
                "-l",
                str(rsrs_file),
                "-o",
                str(rerooted_tree),
                "-t",
                str(threads),
            ],
            logger,
        )
        tree = read_newick_tree(rerooted_tree)
        tree = normalize_tree_binary(tree)
        ok, reason = is_binary_rooted_with_outgroup(tree, outgroup_tip)
        if not ok:
            raise PipelineError(f"Constructed baseml tree is not rooted binary with outgroup: {reason}")
        write_tree(tree, destination)
    finally:
        keep_tip_file.unlink(missing_ok=True)
        rsrs_file.unlink(missing_ok=True)
        pruned_tree.unlink(missing_ok=True)
        rerooted_tree.unlink(missing_ok=True)


def build_overlap_rows(manifest_rows: Sequence[Dict[str, object]], tip_owner_map: Dict[str, str]) -> List[Dict[str, object]]:
    grouped: Dict[str, List[Dict[str, object]]] = {}
    for row in manifest_rows:
        grouped.setdefault(str(row["tip_name"]), []).append(row)
    rows: List[Dict[str, object]] = []
    for tip_name in sorted(grouped):
        entries = grouped[tip_name]
        rows.append(
            {
                "tip_name": tip_name,
                "owner_core_subtree_id": tip_owner_map.get(tip_name, ""),
                "appearing_baseml_subtree_ids": encode_json_list([str(entry["baseml_subtree_id"]) for entry in entries]),
                "roles": encode_json_list([str(entry["tip_role"]) for entry in entries]),
            }
        )
    return rows


def summarize_context(
    rooted_tree: Tree,
    outgroup_tip: str,
    max_tips: int,
    reserve_slots_for_outgroup: int,
) -> Dict[str, object]:
    effective_core_max_tips = max_tips - reserve_slots_for_outgroup
    if effective_core_max_tips < 1:
        raise PipelineError("reserve_slots_for_outgroup leaves no room for core partition tips.")
    total_tips = get_tip_names_from_tree(rooted_tree)
    if outgroup_tip not in total_tips:
        raise PipelineError(f"Required outgroup tip {outgroup_tip} is not present in rooted tree.")
    return {
        "outgroup_tip": outgroup_tip,
        "all_tip_names": total_tips,
        "effective_core_max_tips": effective_core_max_tips,
    }

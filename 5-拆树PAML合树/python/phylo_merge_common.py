#!/usr/bin/env python3
"""Shared utilities for simulated baseml outputs and heuristic tree merging."""

from __future__ import annotations

import hashlib
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.Phylo.BaseTree import Clade, Tree

from phylo_split_common import (
    GLOBAL_OUTGROUP_TIP,
    PipelineError,
    assign_node_ids,
    compute_tip_hash,
    decode_json_list,
    encode_json_list,
    get_tip_names_from_tree,
    is_binary_rooted_with_outgroup,
    load_table,
    read_newick_tree,
    write_table,
    write_tree,
)

SIMULATION_MANIFEST_COLUMNS = [
    "baseml_subtree_id",
    "source_tree_file",
    "simulated_tree_file",
    "tree_seed",
    "n_edges_randomized",
    "min_multiplier",
    "max_multiplier",
    "mean_multiplier",
]

MERGE_REPORT_COLUMNS = [
    "core_subtree_id",
    "core_root_node_id",
    "baseml_subtree_id",
    "analysis_tree_file",
    "core_n_tips",
    "analysis_total_n_tips",
    "matched_internal_edges",
    "matched_tip_edges",
    "replaced_core_root_edge",
    "topology_match",
    "merge_status",
    "details",
]

EDGE_UPDATE_COLUMNS = [
    "edge_role",
    "core_subtree_id",
    "baseml_subtree_id",
    "node_id",
    "descendant_tip_hash",
    "descendant_n_tips",
    "old_branch_length",
    "new_branch_length",
    "scale_multiplier",
    "scale_source_core_ids",
    "changed",
]


@dataclass
class CoreSummaryRecord:
    core_subtree_id: str
    core_root_node_id: str
    parent_node_id: str
    core_n_tips: int
    core_tip_names: List[str]
    tip_hash: str
    core_tree_file: str


@dataclass
class BasemlSummaryRecord:
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


def load_core_summary(path: Path) -> List[CoreSummaryRecord]:
    rows = load_table(path)
    records: List[CoreSummaryRecord] = []
    for row in rows:
        records.append(
            CoreSummaryRecord(
                core_subtree_id=row["core_subtree_id"],
                core_root_node_id=row["core_root_node_id"],
                parent_node_id=row["parent_node_id"],
                core_n_tips=int(row["core_n_tips"]),
                core_tip_names=decode_json_list(row["core_tip_names"]),
                tip_hash=row["tip_hash"],
                core_tree_file=row["core_tree_file"],
            )
        )
    return records


def load_baseml_summary(path: Path) -> List[BasemlSummaryRecord]:
    rows = load_table(path)
    records: List[BasemlSummaryRecord] = []
    for row in rows:
        records.append(
            BasemlSummaryRecord(
                baseml_subtree_id=row["baseml_subtree_id"],
                core_subtree_id=row["core_subtree_id"],
                core_root_node_id=row["core_root_node_id"],
                outgroup_tip=row["outgroup_tip"],
                core_n_tips=int(row["core_n_tips"]),
                anchor_n_tips=int(row["anchor_n_tips"]),
                total_n_tips=int(row["total_n_tips"]),
                core_tip_names=decode_json_list(row["core_tip_names"]),
                anchor_tip_names=decode_json_list(row["anchor_tip_names"]),
                total_tip_names=decode_json_list(row["total_tip_names"]),
                anchor_source_node_ids=decode_json_list(row["anchor_source_node_ids"]),
                tip_hash=row["tip_hash"],
                baseml_tree_file=row["baseml_tree_file"],
            )
        )
    return records


def stable_hash_int(value: str) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:16], 16)


def ensure_positive_branch_length(branch_length: Optional[float], min_branch_length: float) -> float:
    if branch_length is None or branch_length <= 0:
        return float(min_branch_length)
    return float(branch_length)


def format_float(value: float) -> str:
    return f"{float(value):.12g}"


def write_rows(rows: Sequence[Dict[str, object]], columns: Sequence[str], destination: Path) -> None:
    write_table(rows, columns, destination)


def iter_nonroot_clades(tree: Tree) -> Iterable[Clade]:
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        yield clade


def randomize_tree_branch_lengths(
    tree: Tree,
    sigma: float,
    seed: int,
    min_branch_length: float,
) -> Tuple[int, float, float, float]:
    rng = random.Random(seed)
    mu = -(sigma ** 2) / 2.0
    multipliers: List[float] = []
    n_edges = 0
    for clade in iter_nonroot_clades(tree):
        baseline = ensure_positive_branch_length(clade.branch_length, min_branch_length)
        multiplier = math.exp(rng.gauss(mu, sigma))
        clade.branch_length = baseline * multiplier
        multipliers.append(multiplier)
        n_edges += 1
    if not multipliers:
        return 0, 1.0, 1.0, 1.0
    return n_edges, min(multipliers), max(multipliers), sum(multipliers) / len(multipliers)


def validate_tip_hash(expected_tip_names: Sequence[str], expected_hash: str, label: str) -> None:
    actual_hash = compute_tip_hash(expected_tip_names)
    if actual_hash != expected_hash:
        raise PipelineError(f"{label} tip hash mismatch.")


def validate_tree_against_expected_tip_set(tree: Tree, expected_tip_names: Sequence[str], label: str) -> None:
    actual_tip_names = get_tip_names_from_tree(tree)
    if actual_tip_names.count(GLOBAL_OUTGROUP_TIP) != expected_tip_names.count(GLOBAL_OUTGROUP_TIP):
        raise PipelineError(f"{label}: outgroup tip count mismatch.")
    if set(actual_tip_names) != set(expected_tip_names):
        missing = sorted(set(expected_tip_names) - set(actual_tip_names))
        extra = sorted(set(actual_tip_names) - set(expected_tip_names))
        raise PipelineError(f"{label}: tip set mismatch. missing={missing[:5]} extra={extra[:5]}")
    validate_tip_hash(actual_tip_names, compute_tip_hash(expected_tip_names), label)


def get_tip_lookup(tree: Tree) -> Dict[str, Clade]:
    lookup: Dict[str, Clade] = {}
    for tip in tree.get_terminals():
        if tip.name in lookup:
            raise PipelineError(f"Duplicate tip name found in tree: {tip.name}")
        lookup[str(tip.name)] = tip
    return lookup


def find_core_mrca(tree: Tree, core_tip_names: Sequence[str]) -> Clade:
    tip_lookup = get_tip_lookup(tree)
    try:
        targets = [tip_lookup[tip_name] for tip_name in core_tip_names]
    except KeyError as exc:
        raise PipelineError(f"Result tree is missing a core tip: {exc}") from exc
    core_mrca = tree.common_ancestor(targets)
    descendant_tips = get_tip_names_from_clade(core_mrca)
    if set(descendant_tips) != set(core_tip_names):
        raise PipelineError("Core tips are not monophyletic in the analysis tree.")
    return core_mrca


def get_tip_names_from_clade(clade: Clade) -> List[str]:
    names: List[str] = []
    for tip in clade.get_terminals():
        if not tip.name:
            raise PipelineError("Found unnamed terminal tip in clade.")
        names.append(str(tip.name))
    return names


def build_clade_signature_maps(root_clade: Clade) -> Tuple[Dict[str, Clade], Dict[str, int]]:
    signature_to_clade: Dict[str, Clade] = {}
    signature_to_count: Dict[str, int] = {}

    def recurse(clade: Clade) -> Tuple[List[str], str]:
        if clade.is_terminal():
            tip_names = [str(clade.name)]
        else:
            tip_names = []
            for child in clade.clades:
                child_tip_names, _ = recurse(child)
                tip_names.extend(child_tip_names)
            tip_names.sort()
        signature = compute_tip_hash(tip_names)
        signature_to_clade[signature] = clade
        signature_to_count[signature] = len(tip_names)
        return tip_names, signature

    recurse(root_clade)
    return signature_to_clade, signature_to_count


def build_tree_signature_set(tree: Tree) -> set:
    signature_map, _ = build_clade_signature_maps(tree.root)
    return set(signature_map)


def build_analysis_tree_path(
    merge_output_dir: Path,
    baseml_subtree_id: str,
    analysis_tree_source: str,
    external_result_dir: Optional[Path],
) -> Path:
    if analysis_tree_source == "simulated":
        return merge_output_dir / "simulated_baseml_subtrees" / f"{baseml_subtree_id}.nwk"
    if analysis_tree_source == "external":
        if external_result_dir is None:
            raise PipelineError("analysis_tree_source=external requires external_result_dir.")
        return external_result_dir / f"{baseml_subtree_id}.nwk"
    raise PipelineError(f"Unsupported analysis_tree_source: {analysis_tree_source}")


def compute_multiplier(old_value: Optional[float], new_value: Optional[float], min_branch_length: float) -> float:
    old_positive = ensure_positive_branch_length(old_value, min_branch_length)
    new_positive = ensure_positive_branch_length(new_value, min_branch_length)
    return new_positive / old_positive


def weighted_geometric_mean(values: Sequence[Tuple[float, int]]) -> float:
    if not values:
        raise PipelineError("Cannot compute weighted geometric mean of an empty set.")
    numerator = 0.0
    denominator = 0.0
    for value, weight in values:
        if value <= 0:
            raise PipelineError("Weighted geometric mean requires strictly positive values.")
        if weight <= 0:
            continue
        numerator += weight * math.log(value)
        denominator += weight
    if denominator <= 0:
        raise PipelineError("Weighted geometric mean received non-positive total weight.")
    return math.exp(numerator / denominator)


def get_outgroup_child(tree: Tree, outgroup_tip: str = GLOBAL_OUTGROUP_TIP) -> Clade:
    matching_children: List[Clade] = []
    for child in tree.root.clades:
        child_tip_names = set(get_tip_names_from_clade(child))
        if child_tip_names == {outgroup_tip}:
            matching_children.append(child)
    if len(matching_children) != 1:
        raise PipelineError(f"Could not find a unique singleton outgroup child for {outgroup_tip}.")
    return matching_children[0]


def build_descendant_core_ids(
    tree: Tree,
    core_root_node_to_core_id: Dict[str, str],
    node_id_map: Dict[Clade, str],
) -> Dict[Clade, List[str]]:
    descendant_core_ids: Dict[Clade, List[str]] = {}
    for clade in tree.find_clades(order="postorder"):
        node_id = node_id_map[clade]
        if node_id in core_root_node_to_core_id:
            descendant_core_ids[clade] = [core_root_node_to_core_id[node_id]]
            continue
        core_ids = set()
        for child in clade.clades:
            core_ids.update(descendant_core_ids[child])
        descendant_core_ids[clade] = sorted(core_ids)
    return descendant_core_ids


def build_node_id_branch_map(tree: Tree) -> Dict[str, Optional[float]]:
    node_id_map, _, _ = assign_node_ids(tree)
    return {node_id: clade.branch_length for clade, node_id in node_id_map.items() if clade is not tree.root}


def edge_role_for_core(node_id: str, clade: Clade, core_root_node_id: str) -> str:
    if node_id == core_root_node_id:
        return "core_root"
    if clade.is_terminal():
        return "core_tip"
    return "core_internal"


def edge_role_for_scaffold(clade: Clade) -> str:
    if clade.is_terminal():
        return "scaffold_tip"
    return "scaffold_internal"


def build_edge_update_row(
    edge_role: str,
    core_subtree_id: str,
    baseml_subtree_id: str,
    node_id: str,
    descendant_tip_hash: str,
    descendant_n_tips: int,
    old_branch_length: float,
    new_branch_length: float,
    scale_multiplier: float,
    scale_source_core_ids: Sequence[str],
) -> Dict[str, object]:
    return {
        "edge_role": edge_role,
        "core_subtree_id": core_subtree_id,
        "baseml_subtree_id": baseml_subtree_id,
        "node_id": node_id,
        "descendant_tip_hash": descendant_tip_hash,
        "descendant_n_tips": descendant_n_tips,
        "old_branch_length": format_float(old_branch_length),
        "new_branch_length": format_float(new_branch_length),
        "scale_multiplier": format_float(scale_multiplier),
        "scale_source_core_ids": encode_json_list(scale_source_core_ids),
        "changed": "true" if abs(new_branch_length - old_branch_length) > 0 else "false",
    }


def validate_branch_lengths_complete(tree: Tree) -> None:
    for clade in iter_nonroot_clades(tree):
        if clade.branch_length is None or clade.branch_length < 0:
            raise PipelineError("Merged tree contains missing or negative branch lengths.")


def validate_rooted_binary_tree(tree: Tree, outgroup_tip: str = GLOBAL_OUTGROUP_TIP) -> None:
    ok, reason = is_binary_rooted_with_outgroup(tree, outgroup_tip)
    if not ok:
        raise PipelineError(reason)


def write_tree_file(tree: Tree, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_tree(tree, destination)


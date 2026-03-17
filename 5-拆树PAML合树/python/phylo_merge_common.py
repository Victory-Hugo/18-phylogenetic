#!/usr/bin/env python3
"""Shared utilities for backbone-graft merge operations."""

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
    PipelineError,
    assign_node_ids,
    collapse_unary_tree,
    clone_clade,
    clone_tree,
    compute_tip_hash,
    decode_json_list,
    encode_json_list,
    get_tip_names_from_clade,
    get_tip_names_from_tree,
    is_rooted_with_singleton_outgroup,
    is_binary_rooted_with_outgroup,
    load_table,
    normalize_tree_binary,
    read_newick_tree,
    write_table,
    write_tree,
)

SIMULATION_MANIFEST_COLUMNS = [
    "paml_subtree_id",
    "source_tree_file",
    "simulated_tree_file",
    "tree_seed",
    "n_edges_randomized",
    "min_multiplier",
    "max_multiplier",
    "mean_multiplier",
]

BACKBONE_EDGE_ESTIMATE_COLUMNS = [
    "signature_key",
    "descendant_backbone_n_tips",
    "supporting_paml_subtree_ids",
    "n_estimates",
    "aggregated_branch_length",
    "details",
]

GRAFT_REPORT_COLUMNS = [
    "target_subtree_id",
    "target_root_node_id",
    "paml_subtree_id",
    "analysis_tree_file",
    "target_nonbackbone_n_tips",
    "graft_status",
    "attachment_placeholder",
    "topology_match",
    "details",
]

EDGE_UPDATE_COLUMNS = [
    "edge_role",
    "target_subtree_id",
    "paml_subtree_id",
    "node_id",
    "descendant_tip_hash",
    "descendant_n_tips",
    "old_branch_length",
    "new_branch_length",
    "aggregation_source",
    "changed",
]


@dataclass
class BackboneSummaryRecord:
    backbone_tip_name: str
    selection_source: str
    selection_rank: int
    is_user_supplied: bool
    patristic_seed_distance: Optional[float]


@dataclass
class TargetSummaryRecord:
    target_subtree_id: str
    target_root_node_id: str
    parent_node_id: str
    target_nonbackbone_n_tips: int
    target_nonbackbone_tip_names: List[str]
    target_tip_hash: str
    target_tree_file: str


@dataclass
class PamlSummaryRecord:
    paml_subtree_id: str
    target_subtree_id: str
    target_root_node_id: str
    outgroup_tip: str
    backbone_n_tips: int
    target_nonbackbone_n_tips: int
    total_n_tips: int
    backbone_tip_names: List[str]
    target_nonbackbone_tip_names: List[str]
    total_tip_names: List[str]
    tip_hash: str
    paml_tree_file: str

    @property
    def baseml_subtree_id(self) -> str:
        return self.paml_subtree_id

    @property
    def baseml_tree_file(self) -> str:
        return self.paml_tree_file


def load_backbone_summary(path: Path) -> List[BackboneSummaryRecord]:
    rows = load_table(path)
    records: List[BackboneSummaryRecord] = []
    for row in rows:
        records.append(
            BackboneSummaryRecord(
                backbone_tip_name=row["backbone_tip_name"],
                selection_source=row["selection_source"],
                selection_rank=int(row["selection_rank"]),
                is_user_supplied=str(row["is_user_supplied"]).lower() == "true",
                patristic_seed_distance=float(row["patristic_seed_distance"]) if row["patristic_seed_distance"] else None,
            )
        )
    return records


def load_target_summary(path: Path) -> List[TargetSummaryRecord]:
    rows = load_table(path)
    records: List[TargetSummaryRecord] = []
    for row in rows:
        records.append(
            TargetSummaryRecord(
                target_subtree_id=row["target_subtree_id"],
                target_root_node_id=row["target_root_node_id"],
                parent_node_id=row["parent_node_id"],
                target_nonbackbone_n_tips=int(row["target_nonbackbone_n_tips"]),
                target_nonbackbone_tip_names=decode_json_list(row["target_nonbackbone_tip_names"]),
                target_tip_hash=row["target_tip_hash"],
                target_tree_file=row["target_tree_file"],
            )
        )
    return records


def load_paml_summary(path: Path) -> List[PamlSummaryRecord]:
    rows = load_table(path)
    records: List[PamlSummaryRecord] = []
    for row in rows:
        records.append(
            PamlSummaryRecord(
                paml_subtree_id=row["paml_subtree_id"],
                target_subtree_id=row["target_subtree_id"],
                target_root_node_id=row["target_root_node_id"],
                outgroup_tip=row["outgroup_tip"],
                backbone_n_tips=int(row["backbone_n_tips"]),
                target_nonbackbone_n_tips=int(row["target_nonbackbone_n_tips"]),
                total_n_tips=int(row["total_n_tips"]),
                backbone_tip_names=decode_json_list(row["backbone_tip_names"]),
                target_nonbackbone_tip_names=decode_json_list(row["target_nonbackbone_tip_names"]),
                total_tip_names=decode_json_list(row["total_tip_names"]),
                tip_hash=row["tip_hash"],
                paml_tree_file=row["paml_tree_file"],
            )
        )
    return records


def load_baseml_summary(path: Path) -> List[PamlSummaryRecord]:
    return load_paml_summary(path)


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


def validate_tip_hash(expected_tip_names: Sequence[str], expected_hash: str, label: str) -> None:
    actual_hash = compute_tip_hash(expected_tip_names)
    if actual_hash != expected_hash:
        raise PipelineError(f"{label} tip hash mismatch.")


def validate_tree_against_expected_tip_set(
    tree: Tree,
    expected_tip_names: Sequence[str],
    label: str,
    outgroup_tip: str,
) -> None:
    actual_tip_names = get_tip_names_from_tree(tree)
    if actual_tip_names.count(outgroup_tip) != expected_tip_names.count(outgroup_tip):
        raise PipelineError(f"{label}: outgroup tip count mismatch.")
    if set(actual_tip_names) != set(expected_tip_names):
        missing = sorted(set(expected_tip_names) - set(actual_tip_names))
        extra = sorted(set(actual_tip_names) - set(expected_tip_names))
        raise PipelineError(f"{label}: tip set mismatch. missing={missing[:5]} extra={extra[:5]}")
    validate_tip_hash(actual_tip_names, compute_tip_hash(expected_tip_names), label)


def get_tip_lookup(tree: Tree) -> Dict[str, Clade]:
    lookup: Dict[str, Clade] = {}
    for tip in tree.get_terminals():
        tip_name = str(tip.name)
        if tip_name in lookup:
            raise PipelineError(f"Duplicate tip name found in tree: {tip_name}")
        lookup[tip_name] = tip
    return lookup


def build_clade_signature_maps(root_clade: Clade, allowed_tip_set: Optional[set[str]] = None) -> Tuple[Dict[str, Clade], Dict[str, int]]:
    signature_to_clade: Dict[str, Clade] = {}
    signature_to_count: Dict[str, int] = {}

    def recurse(clade: Clade) -> List[str]:
        if clade.is_terminal():
            tip_names = [str(clade.name)] if allowed_tip_set is None or str(clade.name) in allowed_tip_set else []
        else:
            tip_names = []
            for child in clade.clades:
                tip_names.extend(recurse(child))
            tip_names.sort()
        if tip_names:
            signature = compute_tip_hash(tip_names)
            signature_to_clade[signature] = clade
            signature_to_count[signature] = len(tip_names)
        return tip_names

    recurse(root_clade)
    return signature_to_clade, signature_to_count


def build_tree_signature_set(tree: Tree) -> set[str]:
    signature_map, _ = build_clade_signature_maps(tree.root)
    return set(signature_map)


def build_analysis_tree_path(
    merge_output_dir: Path,
    paml_subtree_id: str,
    analysis_tree_source: str,
    external_result_dir: Optional[Path],
) -> Path:
    if analysis_tree_source == "simulated":
        return merge_output_dir / "simulated_baseml_subtrees" / f"{paml_subtree_id}.nwk"
    if analysis_tree_source == "external":
        if external_result_dir is None:
            raise PipelineError("analysis_tree_source=external requires external_result_dir.")
        return external_result_dir / f"{paml_subtree_id}.nwk"
    raise PipelineError(f"Unsupported analysis_tree_source: {analysis_tree_source}")


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


def get_outgroup_child(tree: Tree, outgroup_tip: str) -> Clade:
    matching_children: List[Clade] = []
    for child in tree.root.clades:
        child_tip_names = set(get_tip_names_from_clade(child))
        if child_tip_names == {outgroup_tip}:
            matching_children.append(child)
    if len(matching_children) != 1:
        raise PipelineError(f"Could not find a unique singleton outgroup child for {outgroup_tip}.")
    return matching_children[0]


def validate_branch_lengths_complete(tree: Tree) -> None:
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        if clade.branch_length is None or clade.branch_length < 0:
            raise PipelineError("Merged tree contains missing or negative branch lengths.")


def validate_rooted_binary_tree(tree: Tree, outgroup_tip: str) -> None:
    ok, reason = is_binary_rooted_with_outgroup(tree, outgroup_tip)
    if not ok:
        raise PipelineError(reason)


def validate_rooted_tree_with_outgroup(tree: Tree, outgroup_tip: str) -> None:
    ok, reason = is_rooted_with_singleton_outgroup(tree, outgroup_tip)
    if not ok:
        raise PipelineError(reason)


def write_tree_file(tree: Tree, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_tree(tree, destination)


def iter_nonroot_clades(tree: Tree) -> Iterable[Clade]:
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        yield clade


def randomize_tree_branch_lengths(tree: Tree, sigma: float, seed: int, min_branch_length: float) -> Tuple[int, float, float, float]:
    rng = random.Random(seed)
    mu = -(sigma**2) / 2.0
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


def build_scaffold_tree(
    master_tree: Tree,
    target_records: Sequence[TargetSummaryRecord],
    backbone_tip_names: Sequence[str],
    outgroup_tip_name: str,
) -> Tuple[Tree, Dict[str, str]]:
    del backbone_tip_names
    placeholder_to_target_id: Dict[str, str] = {}
    target_root_map: Dict[str, Tuple[TargetSummaryRecord, str]] = {}
    for record in target_records:
        placeholder = f"TARGET_{record.target_subtree_id.split('_')[-1]}"
        placeholder_to_target_id[placeholder] = record.target_subtree_id
        if record.target_root_node_id in target_root_map:
            raise PipelineError(f"Duplicate target_root_node_id in target summary: {record.target_root_node_id}")
        target_root_map[record.target_root_node_id] = (record, placeholder)

    node_id_map, _, _ = assign_node_ids(master_tree)

    def prune_target_tips(clade: Clade, target_tip_set: set[str]) -> Optional[Clade]:
        if clade.is_terminal():
            tip_name = str(clade.name)
            if tip_name in target_tip_set:
                return None
            return Clade(name=tip_name, branch_length=clade.branch_length)

        kept_children = []
        for child in clade.clades:
            pruned_child = prune_target_tips(child, target_tip_set)
            if pruned_child is not None:
                kept_children.append(pruned_child)
        if not kept_children:
            return None
        return Clade(name=clade.name, branch_length=clade.branch_length, clades=kept_children)

    def recurse(clade: Clade) -> Clade:
        node_id = node_id_map[clade]
        target_entry = target_root_map.get(node_id)
        if target_entry is not None:
            record, placeholder = target_entry
            target_tip_set = set(record.target_nonbackbone_tip_names)
            retained_children = []
            for child in clade.clades:
                pruned_child = prune_target_tips(child, target_tip_set)
                if pruned_child is not None:
                    retained_children.append(pruned_child)
            if not retained_children:
                return Clade(name=placeholder, branch_length=clade.branch_length)
            retained_children.append(Clade(name=placeholder, branch_length=None))
            return Clade(name=clade.name, branch_length=clade.branch_length, clades=retained_children)

        if clade.is_terminal():
            return Clade(name=str(clade.name), branch_length=clade.branch_length)
        return Clade(
            name=clade.name,
            branch_length=clade.branch_length,
            clades=[recurse(child) for child in clade.clades],
        )

    scaffold_tree = Tree(root=recurse(master_tree.root), rooted=True)
    scaffold_tree = collapse_unary_tree(scaffold_tree)
    validate_rooted_tree_with_outgroup(scaffold_tree, outgroup_tip_name)
    return scaffold_tree, placeholder_to_target_id


def build_backbone_signature_key(clade: Clade, backbone_tip_set: set[str], outgroup_tip_name: str) -> Optional[str]:
    descendant_backbone = sorted(tip for tip in get_tip_names_from_clade(clade) if tip in backbone_tip_set)
    if descendant_backbone:
        return compute_tip_hash(descendant_backbone)
    if set(get_tip_names_from_clade(clade)) == {outgroup_tip_name}:
        return "__OUTGROUP__"
    return None


def extract_monophyletic_target_clade(tree: Tree, target_tip_names: Sequence[str]) -> Clade:
    tip_lookup = get_tip_lookup(tree)
    try:
        target_tips = [tip_lookup[tip_name] for tip_name in target_tip_names]
    except KeyError as exc:
        raise PipelineError(f"Analysis tree is missing target tip: {exc}") from exc
    target_mrca = tree.common_ancestor(target_tips)
    extracted = Tree(root=clone_clade(target_mrca), rooted=True)
    target_tip_set = set(target_tip_names)
    for tip in list(extracted.get_terminals()):
        if str(tip.name) not in target_tip_set:
            extracted.prune(target=tip)
    extracted = collapse_unary_tree(extracted)
    extracted_tip_names = sorted(get_tip_names_from_tree(extracted))
    if extracted_tip_names != sorted(target_tip_names):
        raise PipelineError("Could not extract the exact target tip set from the analysis tree.")
    extracted.root.branch_length = target_mrca.branch_length
    return extracted.root


def replace_placeholder_with_clade(tree: Tree, placeholder_name: str, replacement: Clade) -> bool:
    if tree.root.is_terminal():
        return False
    stack = [tree.root]
    while stack:
        parent = stack.pop()
        for idx, child in enumerate(parent.clades):
            if child.is_terminal() and str(child.name) == placeholder_name:
                parent.clades[idx] = replacement
                return True
            if not child.is_terminal():
                stack.append(child)
    return False


def compute_branch_by_signature(tree: Tree, backbone_tip_set: set[str], outgroup_tip_name: str, min_branch_length: float) -> Dict[str, float]:
    signature_map: Dict[str, float] = {}
    for clade in iter_nonroot_clades(tree):
        signature_key = build_backbone_signature_key(clade, backbone_tip_set, outgroup_tip_name)
        if signature_key is None:
            continue
        signature_map[signature_key] = ensure_positive_branch_length(clade.branch_length, min_branch_length)
    return signature_map


def build_edge_update_row(
    edge_role: str,
    target_subtree_id: str,
    paml_subtree_id: str,
    node_id: str,
    descendant_tip_hash: str,
    descendant_n_tips: int,
    old_branch_length: float,
    new_branch_length: float,
    aggregation_source: str,
) -> Dict[str, object]:
    return {
        "edge_role": edge_role,
        "target_subtree_id": target_subtree_id,
        "paml_subtree_id": paml_subtree_id,
        "node_id": node_id,
        "descendant_tip_hash": descendant_tip_hash,
        "descendant_n_tips": descendant_n_tips,
        "old_branch_length": format_float(old_branch_length),
        "new_branch_length": format_float(new_branch_length),
        "aggregation_source": aggregation_source,
        "changed": "true" if abs(new_branch_length - old_branch_length) > 0 else "false",
    }


def copy_tree(tree: Tree) -> Tree:
    return clone_tree(tree)

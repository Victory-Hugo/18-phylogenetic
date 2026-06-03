#!/usr/bin/env python3
"""Shared utilities for ultrastandard backbone-anchored ultrametric grafting."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.Phylo.BaseTree import Clade, Tree

from phylo_merge_common import (
    build_backbone_signature_key,
    compute_branch_by_signature,
    ensure_positive_branch_length,
    load_paml_summary,
    write_rows,
    write_tree_file,
)
from phylo_split_common import (
    PipelineError,
    collapse_unary_tree,
    clone_tree,
    compute_tip_hash,
    get_tip_names_from_clade,
    get_tip_names_from_tree,
    read_newick_tree,
)

ROOT_TO_TIP_REPORT_COLUMNS = [
    "label",
    "n_tips",
    "min_depth",
    "max_depth",
    "delta",
    "target_height",
    "is_ultrametric",
]

BACKBONE_EDGE_ASSIGNMENT_COLUMNS = [
    "node_id",
    "signature_key",
    "descendant_backbone_n_tips",
    "old_branch_length",
    "new_branch_length",
    "status",
    "details",
]

TARGET_SCALING_REPORT_COLUMNS = [
    "target_subtree_id",
    "beast_subtree_id",
    "attachment_placeholder",
    "anchor_node_id",
    "anchor_height",
    "tree_height",
    "available_depth",
    "relative_target_height",
    "scaling_factor",
    "status",
    "details",
]

GRAFT_REPORT_COLUMNS = [
    "target_subtree_id",
    "beast_subtree_id",
    "attachment_placeholder",
    "anchor_node_id",
    "graft_status",
    "details",
]

BACKBONE_ULTRAMETRIC_REPORT_COLUMNS = [
    "label",
    "n_tips",
    "min_depth",
    "max_depth",
    "delta",
    "target_height",
    "is_ultrametric",
]

BACKBONE_ANALYSIS_MANIFEST_COLUMNS = [
    "backbone_job_id",
    "tip_name",
    "tip_role",
    "source_id",
]


def infer_unique_outgroup(paml_summary_tsv: Path) -> str:
    records = load_paml_summary(paml_summary_tsv)
    unique_outgroups = sorted({record.outgroup_tip for record in records})
    if len(unique_outgroups) != 1:
        raise PipelineError(f"Could not infer a unique outgroup tip from {paml_summary_tsv}: {unique_outgroups}")
    return unique_outgroups[0]


def compute_root_to_tip_distances(tree: Tree) -> List[float]:
    if not tree.get_terminals():
        raise PipelineError("Tree contains no terminal tips.")
    return [float(tree.distance(tree.root, tip)) for tip in tree.get_terminals()]


def build_root_to_tip_row(label: str, tree: Tree, tolerance: float) -> Dict[str, str]:
    depths = compute_root_to_tip_distances(tree)
    min_depth = min(depths)
    max_depth = max(depths)
    delta = max_depth - min_depth
    return {
        "label": label,
        "n_tips": str(len(depths)),
        "min_depth": f"{min_depth:.12g}",
        "max_depth": f"{max_depth:.12g}",
        "delta": f"{delta:.12g}",
        "target_height": f"{max_depth:.12g}",
        "is_ultrametric": "true" if delta <= tolerance else "false",
    }


def validate_ultrametric(tree: Tree, tolerance: float, label: str) -> None:
    row = build_root_to_tip_row(label, tree, tolerance)
    if row["is_ultrametric"] != "true":
        raise PipelineError(
            f"{label} is not ultrametric within tolerance={tolerance}: "
            f"min={row['min_depth']} max={row['max_depth']} delta={row['delta']}"
        )


def extend_terminal_branches_to_height(
    tree: Tree,
    min_branch_length: float,
    tolerance: float,
) -> Tuple[Tree, float, Dict[str, str]]:
    tree = clone_tree(tree)
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            clade.branch_length = None
            continue
        clade.branch_length = ensure_positive_branch_length(clade.branch_length, min_branch_length)
    depths = compute_root_to_tip_distances(tree)
    target_height = max(depths)
    pre_row = build_root_to_tip_row("pre_normalization", tree, tolerance)
    for tip in tree.get_terminals():
        current_depth = float(tree.distance(tree.root, tip))
        deficit = target_height - current_depth
        if deficit < -tolerance:
            raise PipelineError(
                f"Tip {tip.name} exceeds target ultrametric height by {abs(deficit):.12g}; "
                "cannot normalize by terminal extension only."
            )
        tip.branch_length = ensure_positive_branch_length(tip.branch_length, min_branch_length) + max(0.0, deficit)
    tree.root.branch_length = None
    post_row = build_root_to_tip_row("post_normalization", tree, tolerance)
    return tree, target_height, {"pre": pre_row, "post": post_row}


def project_tree_by_constrained_height_fit(
    tree: Tree,
    min_branch_length: float,
    tolerance: float,
) -> Tuple[Tree, float, Dict[str, str]]:
    projected_tree = clone_tree(tree)
    for clade in projected_tree.find_clades(order="preorder"):
        if clade is projected_tree.root:
            clade.branch_length = None
            continue
        clade.branch_length = ensure_positive_branch_length(clade.branch_length, min_branch_length)

    root_depths: Dict[Clade, float] = {projected_tree.root: 0.0}
    parent_map: Dict[Clade, Optional[Clade]] = {projected_tree.root: None}
    stack: List[Clade] = [projected_tree.root]
    while stack:
        parent = stack.pop()
        parent_depth = root_depths[parent]
        for child in parent.clades:
            branch_length = ensure_positive_branch_length(child.branch_length, min_branch_length)
            root_depths[child] = parent_depth + branch_length
            parent_map[child] = parent
            stack.append(child)

    pre_row = build_root_to_tip_row("pre_normalization", projected_tree, tolerance)
    node_heights: Dict[Clade, float] = {}

    def postorder_height(clade: Clade) -> Tuple[int, float]:
        if clade.is_terminal():
            node_heights[clade] = 0.0
            return 1, root_depths[clade]

        total_tips = 0
        total_tip_depth = 0.0
        child_min_height = 0.0
        for child in clade.clades:
            child_tip_count, child_depth_sum = postorder_height(child)
            total_tips += child_tip_count
            total_tip_depth += child_depth_sum
            child_min_height = max(child_min_height, node_heights[child] + min_branch_length)

        mean_tip_depth = total_tip_depth / float(total_tips)
        raw_height = max(0.0, mean_tip_depth - root_depths[clade])
        node_heights[clade] = max(raw_height, child_min_height)
        return total_tips, total_tip_depth

    postorder_height(projected_tree.root)

    for clade in projected_tree.find_clades(order="preorder"):
        if clade is projected_tree.root:
            clade.branch_length = None
            continue
        parent = parent_map.get(clade)
        if parent is None:
            raise PipelineError("Failed to determine parent while projecting ultrametric heights.")
        projected_length = max(min_branch_length, node_heights[parent] - node_heights[clade])
        clade.branch_length = projected_length

    post_row = build_root_to_tip_row("post_normalization", projected_tree, tolerance)
    return projected_tree, node_heights[projected_tree.root], {"pre": pre_row, "post": post_row}


def project_tree_to_ultrametric(
    tree: Tree,
    projection_mode: str,
    min_branch_length: float,
    tolerance: float,
) -> Tuple[Tree, float, Dict[str, str]]:
    if projection_mode == "extend_terminal_to_max_depth":
        return extend_terminal_branches_to_height(
            tree,
            min_branch_length=min_branch_length,
            tolerance=tolerance,
        )
    if projection_mode == "constrained_height_fit":
        return project_tree_by_constrained_height_fit(
            tree,
            min_branch_length=min_branch_length,
            tolerance=tolerance,
        )
    raise PipelineError(f"Unsupported ultrametric projection mode: {projection_mode}")


def normalize_relative_ultrametric_tree(
    tree: Tree,
    min_branch_length: float,
    tolerance: float,
) -> Tuple[Tree, float, Dict[str, str]]:
    ultrametric_tree, target_height, report_rows = extend_terminal_branches_to_height(
        tree,
        min_branch_length=min_branch_length,
        tolerance=tolerance,
    )
    if target_height <= 0:
        tip_names = get_tip_names_from_tree(ultrametric_tree)
        if len(tip_names) != 1:
            raise PipelineError("Target tree height must be positive for relative ultrametric scaling.")
        synthetic_tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        name=tip_names[0],
                        branch_length=1.0,
                    )
                ]
            ),
            rooted=True,
        )
        synthetic_tree.root.branch_length = None
        return synthetic_tree, 0.0, report_rows
    for clade in ultrametric_tree.find_clades(order="preorder"):
        if clade is ultrametric_tree.root:
            clade.branch_length = None
            continue
        branch_length = ensure_positive_branch_length(clade.branch_length, min_branch_length)
        clade.branch_length = branch_length / target_height
    ultrametric_tree.root.branch_length = None
    return ultrametric_tree, target_height, report_rows


def scale_tree_nonroot(tree: Tree, factor: float, min_branch_length: float) -> Tree:
    if factor <= 0:
        raise PipelineError(f"Scaling factor must be positive, got {factor}")
    scaled_tree = clone_tree(tree)
    for clade in scaled_tree.find_clades(order="preorder"):
        if clade is scaled_tree.root:
            clade.branch_length = None
            continue
        clade.branch_length = max(min_branch_length, ensure_positive_branch_length(clade.branch_length, min_branch_length) * factor)
    scaled_tree.root.branch_length = None
    return scaled_tree


def build_backbone_signature_length_map(
    backbone_tree: Tree,
    backbone_tip_set: set[str],
    outgroup_tip_name: str,
    min_branch_length: float,
) -> Dict[str, float]:
    signature_lengths = compute_branch_by_signature(
        backbone_tree,
        backbone_tip_set=backbone_tip_set,
        outgroup_tip_name=outgroup_tip_name,
        min_branch_length=min_branch_length,
    )
    if not signature_lengths:
        raise PipelineError("Could not compute any signature-based backbone branch lengths.")
    return signature_lengths


def find_placeholder_parent(tree: Tree, placeholder_name: str) -> Tuple[Clade, Clade]:
    stack: List[Clade] = [tree.root]
    while stack:
        parent = stack.pop()
        for child in parent.clades:
            if child.is_terminal() and str(child.name) == placeholder_name:
                return parent, child
            if not child.is_terminal():
                stack.append(child)
    raise PipelineError(f"Could not locate placeholder {placeholder_name} in scaffold tree.")


def count_backbone_descendants(tree: Tree, backbone_tip_set: set[str]) -> int:
    return sum(1 for tip_name in get_tip_names_from_tree(tree) if tip_name in backbone_tip_set)


def build_signature_tip_hash(clade: Clade, backbone_tip_set: set[str]) -> str:
    descendant_backbone = sorted(tip for tip in get_tip_names_from_clade(clade) if tip in backbone_tip_set)
    if not descendant_backbone:
        return compute_tip_hash(["__NONE__"])
    return compute_tip_hash(descendant_backbone)


def assign_backbone_lengths_to_scaffold(
    scaffold_tree: Tree,
    backbone_tree: Tree,
    backbone_tip_set: set[str],
    outgroup_tip_name: str,
    min_branch_length: float,
    node_id_map: Dict[Clade, str],
) -> List[Dict[str, str]]:
    signature_lengths = build_backbone_signature_length_map(
        backbone_tree=backbone_tree,
        backbone_tip_set=backbone_tip_set,
        outgroup_tip_name=outgroup_tip_name,
        min_branch_length=min_branch_length,
    )
    report_rows: List[Dict[str, str]] = []
    assigned_count = 0
    for clade, node_id in node_id_map.items():
        if clade is scaffold_tree.root:
            continue
        signature_key = build_backbone_signature_key(clade, backbone_tip_set, outgroup_tip_name)
        if signature_key is None:
            continue
        old_branch_length = ensure_positive_branch_length(clade.branch_length, min_branch_length)
        new_branch_length = signature_lengths.get(signature_key)
        if new_branch_length is None:
            report_rows.append(
                {
                    "node_id": node_id,
                    "signature_key": signature_key,
                    "descendant_backbone_n_tips": str(
                        len([tip for tip in get_tip_names_from_clade(clade) if tip in backbone_tip_set])
                    ),
                    "old_branch_length": f"{old_branch_length:.12g}",
                    "new_branch_length": "",
                    "status": "missing_signature",
                    "details": "Signature not present in backbone ultrametric tree.",
                }
            )
            continue
        clade.branch_length = new_branch_length
        assigned_count += 1
        report_rows.append(
            {
                "node_id": node_id,
                "signature_key": signature_key,
                "descendant_backbone_n_tips": str(
                    len([tip for tip in get_tip_names_from_clade(clade) if tip in backbone_tip_set])
                ),
                "old_branch_length": f"{old_branch_length:.12g}",
                "new_branch_length": f"{new_branch_length:.12g}",
                "status": "assigned",
                "details": "Assigned from backbone ultrametric signature map.",
            }
        )
    if assigned_count == 0:
        raise PipelineError("No scaffold backbone edges were assigned from the backbone ultrametric tree.")
    failures = [row for row in report_rows if row["status"] != "assigned"]
    if failures:
        preview = "; ".join(f"{row['node_id']}: {row['signature_key']}" for row in failures[:5])
        raise PipelineError(f"{len(failures)} scaffold backbone edge(s) could not be assigned: {preview}")
    return report_rows


def tree_height(tree: Tree) -> float:
    return max(compute_root_to_tip_distances(tree))


def replace_placeholder_with_scaled_tree(tree: Tree, placeholder_name: str, replacement_root: Clade) -> None:
    parent, placeholder = find_placeholder_parent(tree, placeholder_name)
    for idx, child in enumerate(parent.clades):
        if child is placeholder:
            parent.clades[idx] = replacement_root
            return
    raise PipelineError(f"Failed to replace placeholder {placeholder_name}.")


def resolve_backbone_ultrametric_tree_path(paml_output_dir: Path, configured_path: Optional[str]) -> Path:
    if configured_path in (None, "", "null"):
        return paml_output_dir / "backbone_analysis" / "backbone_ultrametric_tree.nwk"
    path = Path(configured_path)
    return path if path.is_absolute() else (paml_output_dir / path).resolve()


def write_root_to_tip_report(rows: Sequence[Dict[str, str]], destination: Path) -> None:
    write_rows(rows, ROOT_TO_TIP_REPORT_COLUMNS, destination)


def write_simple_tree(destination: Path, tree: Tree) -> None:
    write_tree_file(tree, destination)


def read_tree(path: Path) -> Tree:
    return read_newick_tree(path)


def collapse_tree(tree: Tree) -> Tree:
    return collapse_unary_tree(tree)


def iter_internal_nonroot_clades(tree: Tree) -> Iterable[Clade]:
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        if not clade.is_terminal():
            yield clade

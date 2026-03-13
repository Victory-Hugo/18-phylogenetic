#!/usr/bin/env python3
"""Deterministic backbone tip sampling for phylogenetic trees."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from Bio.Phylo.BaseTree import Clade, Tree


def _pipeline_error(message: str):
    from phylo_split_common import PipelineError

    return PipelineError(message)


@dataclass
class BackboneSelectionStep:
    tip_name: str
    selection_source: str
    selection_rank: int
    patristic_seed_distance: float


def _build_tip_lookup(tree: Tree) -> Dict[str, Clade]:
    lookup: Dict[str, Clade] = {}
    for tip in tree.get_terminals():
        tip_name = str(tip.name)
        if tip_name in lookup:
            raise _pipeline_error(f"Duplicate tip name found in tree: {tip_name}")
        lookup[tip_name] = tip
    return lookup


def _tip_names(clade: Clade) -> List[str]:
    return sorted(str(tip.name) for tip in clade.get_terminals())


def _clade_sort_key(clade: Clade) -> Tuple[int, str]:
    tip_names = _tip_names(clade)
    return (-len(tip_names), tip_names[0])


def _build_frontier_clades(root_clade: Clade, backbone_size: int) -> List[Clade]:
    frontier: List[Clade] = [root_clade]
    while len(frontier) < backbone_size:
        split_index = None
        split_clade = None
        for idx, clade in enumerate(sorted(frontier, key=_clade_sort_key)):
            if clade.is_terminal():
                continue
            split_clade = clade
            split_index = frontier.index(clade)
            break
        if split_clade is None or split_index is None:
            raise _pipeline_error("Could not split frontier clades enough to satisfy backbone_size.")
        frontier.pop(split_index)
        frontier.extend(split_clade.clades)
    return sorted(frontier, key=_clade_sort_key)


def greedy_farthest_patristic_sampling(
    tree: Tree,
    ingroup_tip_names: Sequence[str],
    backbone_size: int,
) -> List[BackboneSelectionStep]:
    if backbone_size < 1:
        raise _pipeline_error("runtime.backbone_size must be at least 1.")
    unique_tip_names = sorted(set(str(name) for name in ingroup_tip_names))
    if len(unique_tip_names) < backbone_size:
        raise _pipeline_error(
            f"Requested backbone_size={backbone_size}, but only {len(unique_tip_names)} ingroup tips are available."
        )

    tip_lookup = _build_tip_lookup(tree)
    mrca = tree.common_ancestor([tip_lookup[tip_name] for tip_name in unique_tip_names])
    frontier = _build_frontier_clades(mrca, backbone_size)
    steps: List[BackboneSelectionStep] = []
    for rank, frontier_clade in enumerate(frontier, start=1):
        candidates = _tip_names(frontier_clade)
        best_tip = sorted(
            candidates,
            key=lambda tip_name: (
                -float(tree.distance(mrca, tip_lookup[tip_name])),
                tip_name,
            ),
        )[0]
        score = float(tree.distance(mrca, tip_lookup[best_tip]))
        steps.append(
            BackboneSelectionStep(
                tip_name=best_tip,
                selection_source="auto_frontier_deepest_tip",
                selection_rank=rank,
                patristic_seed_distance=score,
            )
        )
    return steps

#!/usr/bin/env python3
"""Shared utilities for backbone-target phylogenetic splitting."""

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

from backbone_sampling import BackboneSelectionStep, greedy_farthest_patristic_sampling


class PipelineError(RuntimeError):
    """Raised when the pipeline encounters a recoverable domain error."""


DEFAULT_OUTGROUP_TIP_NAME = "RSRS_MRCA"

BACKBONE_SUMMARY_COLUMNS = [
    "backbone_tip_name",
    "selection_source",
    "selection_rank",
    "is_user_supplied",
    "patristic_seed_distance",
]

TARGET_SUMMARY_COLUMNS = [
    "target_subtree_id",
    "target_root_node_id",
    "parent_node_id",
    "target_nonbackbone_n_tips",
    "target_nonbackbone_tip_names",
    "target_tip_hash",
    "target_tree_file",
]

TARGET_MANIFEST_COLUMNS = [
    "target_subtree_id",
    "tip_name",
    "tip_role",
    "source_node_id",
]

PAML_SUMMARY_COLUMNS = [
    "paml_subtree_id",
    "target_subtree_id",
    "target_root_node_id",
    "outgroup_tip",
    "backbone_n_tips",
    "target_nonbackbone_n_tips",
    "total_n_tips",
    "backbone_tip_names",
    "target_nonbackbone_tip_names",
    "total_tip_names",
    "tip_hash",
    "paml_tree_file",
]

PAML_MANIFEST_COLUMNS = [
    "paml_subtree_id",
    "tip_name",
    "tip_role",
    "target_subtree_id",
    "source_node_id",
]

VALIDATION_COLUMNS = ["check_name", "status", "details"]


@dataclass
class BackboneTipRecord:
    backbone_tip_name: str
    selection_source: str
    selection_rank: int
    is_user_supplied: bool
    patristic_seed_distance: Optional[float]

    def to_row(self) -> Dict[str, object]:
        return {
            "backbone_tip_name": self.backbone_tip_name,
            "selection_source": self.selection_source,
            "selection_rank": self.selection_rank,
            "is_user_supplied": "true" if self.is_user_supplied else "false",
            "patristic_seed_distance": ""
            if self.patristic_seed_distance is None
            else format_float(self.patristic_seed_distance),
        }


@dataclass
class TargetSubtreeRecord:
    target_subtree_id: str
    target_root_node_id: str
    parent_node_id: str
    target_nonbackbone_n_tips: int
    target_nonbackbone_tip_names: List[str]
    target_tip_hash: str
    target_tree_file: str
    clade: Clade

    def to_row(self) -> Dict[str, object]:
        return {
            "target_subtree_id": self.target_subtree_id,
            "target_root_node_id": self.target_root_node_id,
            "parent_node_id": self.parent_node_id,
            "target_nonbackbone_n_tips": self.target_nonbackbone_n_tips,
            "target_nonbackbone_tip_names": encode_json_list(self.target_nonbackbone_tip_names),
            "target_tip_hash": self.target_tip_hash,
            "target_tree_file": self.target_tree_file,
        }


@dataclass
class PamlSubtreeRecord:
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

    def to_row(self) -> Dict[str, object]:
        return {
            "paml_subtree_id": self.paml_subtree_id,
            "target_subtree_id": self.target_subtree_id,
            "target_root_node_id": self.target_root_node_id,
            "outgroup_tip": self.outgroup_tip,
            "backbone_n_tips": self.backbone_n_tips,
            "target_nonbackbone_n_tips": self.target_nonbackbone_n_tips,
            "total_n_tips": self.total_n_tips,
            "backbone_tip_names": encode_json_list(self.backbone_tip_names),
            "target_nonbackbone_tip_names": encode_json_list(self.target_nonbackbone_tip_names),
            "total_tip_names": encode_json_list(self.total_tip_names),
            "tip_hash": self.tip_hash,
            "paml_tree_file": self.paml_tree_file,
        }


def format_float(value: float) -> str:
    return f"{float(value):.12g}"


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
    if not tree_path.exists():
        raise PipelineError(f"Newick tree file not found: {tree_path}")
    if tree_path.stat().st_size == 0:
        raise PipelineError(f"Newick tree file is empty: {tree_path}")
    try:
        return Phylo.read(str(tree_path), "newick")
    except Exception as exc:  # pragma: no cover
        raise PipelineError(f"Failed to read Newick tree: {tree_path}") from exc


def estimate_tree_depth(tree: Tree) -> int:
    max_depth = 0
    stack: List[Tuple[Clade, int]] = [(tree.root, 1)]
    while stack:
        clade, depth = stack.pop()
        if depth > max_depth:
            max_depth = depth
        for child in clade.clades:
            stack.append((child, depth + 1))
    return max_depth


def clone_clade(clade: Clade) -> Clade:
    cloned_root = Clade()
    stack: List[Tuple[Clade, Clade]] = [(clade, cloned_root)]
    while stack:
        source, target = stack.pop()
        for key, value in vars(source).items():
            if key == "clades":
                continue
            setattr(target, key, copy.copy(value))
        target.clades = []
        for child in source.clades:
            cloned_child = Clade()
            target.clades.append(cloned_child)
            stack.append((child, cloned_child))
    return cloned_root


def clone_tree(tree: Tree) -> Tree:
    cloned = Tree(root=clone_clade(tree.root), rooted=tree.rooted)
    for key, value in vars(tree).items():
        if key == "root":
            continue
        setattr(cloned, key, copy.copy(value))
    return cloned


def write_tree(tree: Tree, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tree.rooted = True
    depth = estimate_tree_depth(tree)
    previous_limit = sys.getrecursionlimit()
    required_limit = max(previous_limit, depth * 4 + 500)
    if required_limit > previous_limit:
        sys.setrecursionlimit(required_limit)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=destination.parent,
        delete=False,
        suffix=destination.suffix or ".nwk",
    ) as handle:
        tmp_path = Path(handle.name)
    try:
        Phylo.write(tree, str(tmp_path), "newick", format_branch_length="%1.12g")
        if tmp_path.stat().st_size == 0:
            raise PipelineError(f"Generated Newick tree is empty: {tmp_path}")
        tmp_path.replace(destination)
    finally:
        tmp_path.unlink(missing_ok=True)
        if required_limit > previous_limit:
            sys.setrecursionlimit(previous_limit)


def get_tip_names_from_tree(tree: Tree) -> List[str]:
    names = []
    for tip in tree.get_terminals():
        if not tip.name:
            raise PipelineError("Found unnamed terminal tip in tree.")
        names.append(str(tip.name))
    return names


def get_tip_names_from_clade(clade: Clade) -> List[str]:
    return [str(tip.name) for tip in clade.get_terminals()]


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


def resolve_outgroup_tip_name(outgroup_tip_file: Optional[Path], outgroup_tip_name: Optional[str]) -> str:
    if outgroup_tip_name not in (None, "", "null"):
        return str(outgroup_tip_name)
    if outgroup_tip_file is None:
        raise PipelineError("Either --outgroup-tip-name or --outgroup-tip-file must be provided.")
    outgroup_tips = parse_tip_file(outgroup_tip_file)
    if len(outgroup_tips) != 1:
        raise PipelineError(
            "--outgroup-tip-file must contain exactly one tip when --outgroup-tip-name is omitted."
        )
    return outgroup_tips[0]


def validate_outgroup_tips(outgroup_tips: Sequence[str], tree_tip_names: Iterable[str]) -> None:
    tree_tip_set = set(tree_tip_names)
    missing = [tip for tip in outgroup_tips if tip not in tree_tip_set]
    if missing:
        raise PipelineError(f"Outgroup tips not found in input tree: {', '.join(missing[:10])}")


def load_backbone_tip_ids_from_file(path: Path, tree_tip_names: Iterable[str]) -> List[str]:
    backbone_tips = parse_tip_file(path)
    tree_tip_set = set(tree_tip_names)
    missing = [tip for tip in backbone_tips if tip not in tree_tip_set]
    if missing:
        raise PipelineError(f"Backbone tips not found in input tree: {', '.join(missing[:10])}")
    return backbone_tips


def load_backbone_tips_from_tree(backbone_tree: Path, tree_tip_names: Iterable[str]) -> List[str]:
    tree_tip_set = set(tree_tip_names)
    backbone = read_newick_tree(backbone_tree)
    backbone_tips = get_tip_names_from_tree(backbone)
    missing = [tip for tip in backbone_tips if tip not in tree_tip_set]
    if missing:
        raise PipelineError(f"Backbone tips not found in input tree: {', '.join(missing[:10])}")
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
    outgroup_tip_name: Optional[str] = None,
    pre_binarize_rooted_tree: bool = True,
) -> str:
    if rooted_status is None:
        rooted_status, _ = detect_tree_rooted_status(input_tree, conda_env, gotree_bin, threads, logger)
    rooted_tree.parent.mkdir(parents=True, exist_ok=True)
    if rooted_status == "rooted":
        shutil.copy2(input_tree, rooted_tree)
        logger.info("Input tree is already rooted; copied directly to %s", rooted_tree)
    else:
        if outgroup_tip_file is None:
            raise PipelineError("Input tree is unrooted and no outgroup tip file was provided.")
        reroot_with_outgroup(input_tree, outgroup_tip_file, rooted_tree, conda_env, gotree_bin, threads, logger)
    if pre_binarize_rooted_tree:
        tree = read_newick_tree(rooted_tree)
        if outgroup_tip_name is not None:
            ok, reason = is_rooted_with_singleton_outgroup(tree, outgroup_tip_name)
            if not ok:
                raise PipelineError(f"Rooted tree cannot be binarized with singleton outgroup {outgroup_tip_name}: {reason}")
            was_binary, _ = is_binary_rooted_with_outgroup(tree, outgroup_tip_name)
        else:
            was_binary = False
        tree = normalize_tree_binary(tree)
        if outgroup_tip_name is not None:
            ok, reason = is_binary_rooted_with_outgroup(tree, outgroup_tip_name)
            if not ok:
                raise PipelineError(f"Binarized rooted tree failed rooted-binary validation: {reason}")
        write_tree(tree, rooted_tree)
        if was_binary:
            logger.info("Intermediate rooted tree already satisfied rooted-binary constraints: %s", rooted_tree)
        else:
            logger.info("Normalized intermediate rooted tree to rooted binary form: %s", rooted_tree)
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


def encode_json_list(items: Sequence[str]) -> str:
    return json.dumps(list(items), ensure_ascii=False)


def decode_json_list(value: str) -> List[str]:
    parsed = json.loads(value)
    if not isinstance(parsed, list):
        raise PipelineError("Expected a JSON list in table.")
    return [str(item) for item in parsed]


def compute_tip_hash(tip_names: Sequence[str]) -> str:
    payload = "\n".join(sorted(tip_names)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


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


def write_validation_report(rows: Sequence[Tuple[str, str, str]], destination: Path) -> None:
    write_table(
        [{"check_name": name, "status": status, "details": details} for name, status, details in rows],
        VALIDATION_COLUMNS,
        destination,
    )


def clean_generated_outputs(output_dir: Path) -> None:
    for pattern in [
        "target_subtrees/target_subtree_*.nwk",
        "paml_subtrees/paml_subtree_*.nwk",
    ]:
        for candidate in output_dir.glob(pattern):
            candidate.unlink(missing_ok=True)
    for relative_path in [
        "backbone_tips.txt",
        "backbone_tree.nwk",
        "backbone_summary.tsv",
        "target_subtree_summary.tsv",
        "target_tree_manifest.tsv",
        "paml_subtree_summary.tsv",
        "paml_tree_manifest.tsv",
        "split_validation_report.tsv",
        "split_tree.log",
        "intermediate/rooted.tree",
        "intermediate/rooted.validation.tree",
    ]:
        (output_dir / relative_path).unlink(missing_ok=True)


def get_root_children_for_outgroup(tree: Tree, outgroup_tip: str) -> Tuple[Clade, Clade]:
    if len(tree.root.clades) != 2:
        raise PipelineError("Rooted input tree must have exactly 2 root children.")
    outgroup_children = [child for child in tree.root.clades if set(get_tip_names_from_clade(child)) == {outgroup_tip}]
    if len(outgroup_children) != 1:
        raise PipelineError(f"Could not identify unique singleton root child for outgroup {outgroup_tip}.")
    outgroup_child = outgroup_children[0]
    ingroup_child = next(child for child in tree.root.clades if child is not outgroup_child)
    return outgroup_child, ingroup_child


def collapse_unary_tree(tree: Tree) -> Tree:
    def collapse(clade: Clade) -> Clade:
        if clade.is_terminal():
            return clade
        clade.clades = [collapse(child) for child in clade.clades]
        while len(clade.clades) == 1 and not clade.is_terminal():
            child = clade.clades[0]
            child.branch_length = (
                (child.branch_length or 0.0) + (clade.branch_length or 0.0)
                if child.branch_length is not None or clade.branch_length is not None
                else None
            )
            clade = child
        return clade

    tree.root = collapse(tree.root)
    tree.root.branch_length = None
    tree.rooted = True
    return tree


def normalize_tree_binary(tree: Tree) -> Tree:
    tree = collapse_unary_tree(tree)

    def bifurcate(clade: Clade) -> Clade:
        if clade.is_terminal():
            return clade
        clade.clades = [bifurcate(child) for child in clade.clades]
        while len(clade.clades) > 2:
            right_child = clade.clades.pop()
            left_child = clade.clades.pop()
            synthetic = Clade(branch_length=0.0, clades=[left_child, right_child])
            clade.clades.append(synthetic)
        return clade

    tree.root = bifurcate(tree.root)
    tree.root.branch_length = None
    tree.rooted = True
    return tree


def is_rooted_with_singleton_outgroup(tree: Tree, outgroup_tip: str) -> Tuple[bool, str]:
    if len(tree.root.clades) < 2:
        return False, "Root does not have at least 2 children."
    root_child_tip_sets = [set(get_tip_names_from_clade(child)) for child in tree.root.clades]
    if root_child_tip_sets.count({outgroup_tip}) != 1:
        return False, f"Root does not contain exactly one singleton outgroup child {outgroup_tip}."
    for clade in tree.get_nonterminals(order="preorder"):
        if clade is tree.root:
            continue
        if len(clade.clades) < 2:
            return False, f"Internal node has fewer than 2 children: {clade.name or '<internal>'}"
    return True, "Tree is rooted with a singleton outgroup."


def is_binary_rooted_with_outgroup(tree: Tree, outgroup_tip: str) -> Tuple[bool, str]:
    if len(tree.root.clades) != 2:
        return False, "Root does not have exactly 2 children."
    root_child_tip_sets = [set(get_tip_names_from_clade(child)) for child in tree.root.clades]
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
        outgroup_file = Path(handle.name)
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
                str(outgroup_file),
                "-o",
                str(rerooted_tree),
                "-t",
                str(threads),
            ],
            logger,
        )
        if not pruned_tree.exists() or pruned_tree.stat().st_size == 0:
            raise PipelineError(f"gotree prune produced an empty tree: {pruned_tree}")
        if not rerooted_tree.exists() or rerooted_tree.stat().st_size == 0:
            raise PipelineError(f"gotree reroot produced an empty tree: {rerooted_tree}")
        tree = read_newick_tree(rerooted_tree)
        tree = normalize_tree_binary(tree)
        ok, reason = is_binary_rooted_with_outgroup(tree, outgroup_tip)
        if not ok:
            raise PipelineError(f"Constructed subtree is not rooted binary with outgroup: {reason}")
        write_tree(tree, destination)
        if not destination.exists() or destination.stat().st_size == 0:
            raise PipelineError(f"Failed to write induced subtree: {destination}")
    finally:
        keep_tip_file.unlink(missing_ok=True)
        outgroup_file.unlink(missing_ok=True)
        pruned_tree.unlink(missing_ok=True)
        rerooted_tree.unlink(missing_ok=True)


def resolve_backbone_selection(
    rooted_tree: Tree,
    outgroup_tip_name: str,
    backbone_tree: Optional[Path],
    backbone_tip_id_file: Optional[Path],
    backbone_size: int,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[str], List[BackboneTipRecord], str]:
    tree_tip_names = get_tip_names_from_tree(rooted_tree)
    if outgroup_tip_name not in tree_tip_names:
        raise PipelineError(f"Required outgroup tip {outgroup_tip_name} is missing from rooted tree.")
    if backbone_tip_id_file is not None:
        if backbone_tree is not None and logger is not None:
            logger.info("Using backbone_tip_id_file as backbone source; backbone_tree ignored.")
        backbone_tips = load_backbone_tip_ids_from_file(backbone_tip_id_file, tree_tip_names)
        if outgroup_tip_name in backbone_tips:
            raise PipelineError("backbone_tip_id_file must not include the singleton outgroup tip.")
        records = [
            BackboneTipRecord(
                backbone_tip_name=tip_name,
                selection_source="tip_id_file",
                selection_rank=index,
                is_user_supplied=True,
                patristic_seed_distance=None,
            )
            for index, tip_name in enumerate(backbone_tips, start=1)
        ]
        return backbone_tips, records, "tip_id_file"
    if backbone_tree is not None:
        backbone_tips = load_backbone_tips_from_tree(backbone_tree, tree_tip_names)
        if outgroup_tip_name in backbone_tips:
            raise PipelineError("backbone_tree must not include the singleton outgroup tip.")
        records = [
            BackboneTipRecord(
                backbone_tip_name=tip_name,
                selection_source="backbone_tree",
                selection_rank=index,
                is_user_supplied=True,
                patristic_seed_distance=None,
            )
            for index, tip_name in enumerate(backbone_tips, start=1)
        ]
        return backbone_tips, records, "backbone_tree"

    _, ingroup_child = get_root_children_for_outgroup(rooted_tree, outgroup_tip_name)
    ingroup_tip_names = sorted(set(get_tip_names_from_clade(ingroup_child)) - {outgroup_tip_name})
    steps = greedy_farthest_patristic_sampling(rooted_tree, ingroup_tip_names, backbone_size)
    records = [
        BackboneTipRecord(
            backbone_tip_name=step.tip_name,
            selection_source=step.selection_source,
            selection_rank=step.selection_rank,
            is_user_supplied=False,
            patristic_seed_distance=step.patristic_seed_distance,
        )
        for step in steps
    ]
    return [step.tip_name for step in steps], records, "auto"


def build_ordered_tip_list(tip_order_index: Dict[str, int], selected_tip_set: set[str]) -> List[str]:
    return sorted(selected_tip_set, key=tip_order_index.__getitem__)


def write_tip_list(tip_names: Sequence[str], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(tip_names) + "\n", encoding="utf-8")


def compute_target_partition_profiles(
    tree: Tree,
    backbone_tip_set: set[str],
    outgroup_tip: str,
) -> Tuple[Dict[Clade, int], Dict[Clade, int], Dict[Clade, List[str]]]:
    nonbackbone_counts: Dict[Clade, int] = {}
    backbone_counts: Dict[Clade, int] = {}
    ordered_nonbackbone_tips: Dict[Clade, List[str]] = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            tip_name = str(clade.name)
            if tip_name == outgroup_tip:
                nonbackbone_counts[clade] = 0
                backbone_counts[clade] = 0
                ordered_nonbackbone_tips[clade] = []
            elif tip_name in backbone_tip_set:
                nonbackbone_counts[clade] = 0
                backbone_counts[clade] = 1
                ordered_nonbackbone_tips[clade] = []
            else:
                nonbackbone_counts[clade] = 1
                backbone_counts[clade] = 0
                ordered_nonbackbone_tips[clade] = [tip_name]
            continue
        nonbackbone_counts[clade] = sum(nonbackbone_counts[child] for child in clade.clades)
        backbone_counts[clade] = sum(backbone_counts[child] for child in clade.clades)
        tip_names: List[str] = []
        for child in clade.clades:
            tip_names.extend(ordered_nonbackbone_tips[child])
        ordered_nonbackbone_tips[clade] = tip_names
    return nonbackbone_counts, backbone_counts, ordered_nonbackbone_tips


def build_target_partition(
    tree: Tree,
    node_id_map: Dict[Clade, str],
    parent_map: Dict[Clade, Optional[Clade]],
    backbone_tip_set: set[str],
    outgroup_tip: str,
    target_capacity: int,
    logger: logging.Logger,
) -> List[TargetSubtreeRecord]:
    if target_capacity < 1:
        raise PipelineError("runtime.max_tips - runtime.backbone_size - 1 leaves no room for target tips.")

    _, ingroup_child = get_root_children_for_outgroup(tree, outgroup_tip)
    nonbackbone_counts, backbone_counts, ordered_nonbackbone_tips = compute_target_partition_profiles(
        tree,
        backbone_tip_set,
        outgroup_tip,
    )
    selected: List[Clade] = []

    def recurse(clade: Clade) -> None:
        nonbackbone_n = nonbackbone_counts[clade]
        if nonbackbone_n == 0:
            return
        if nonbackbone_n <= target_capacity:
            selected.append(clade)
            return
        if clade.is_terminal():
            if str(clade.name) == outgroup_tip or str(clade.name) in backbone_tip_set:
                return
            selected.append(clade)
            return
        for child in clade.clades:
            recurse(child)

    recurse(ingroup_child)

    records: List[TargetSubtreeRecord] = []
    seen_nonbackbone_tips: set[str] = set()
    expected_nonbackbone_tips = set(ordered_nonbackbone_tips[ingroup_child])

    for index, clade in enumerate(selected, start=1):
        tip_names = list(ordered_nonbackbone_tips[clade])
        if not tip_names:
            continue
        overlap = seen_nonbackbone_tips.intersection(tip_names)
        if overlap:
            raise PipelineError(f"Target partition overlaps on tips: {', '.join(sorted(overlap)[:10])}")
        seen_nonbackbone_tips.update(tip_names)
        parent = parent_map[clade]
        parent_node_id = "ROOT" if parent is None or parent is tree.root else node_id_map[parent]
        record = TargetSubtreeRecord(
            target_subtree_id=f"target_subtree_{index:04d}",
            target_root_node_id=node_id_map[clade],
            parent_node_id=parent_node_id,
            target_nonbackbone_n_tips=len(tip_names),
            target_nonbackbone_tip_names=tip_names,
            target_tip_hash=compute_tip_hash(tip_names),
            target_tree_file=Path("target_subtrees", f"target_subtree_{index:04d}.nwk").as_posix(),
            clade=clade,
        )
        records.append(record)
        logger.info(
            "%s root=%s nonbackbone_tips=%d backbone_descendants=%d",
            record.target_subtree_id,
            record.target_root_node_id,
            record.target_nonbackbone_n_tips,
            backbone_counts[clade],
        )

    missing = expected_nonbackbone_tips - seen_nonbackbone_tips
    extra = seen_nonbackbone_tips - expected_nonbackbone_tips
    if missing or extra:
        raise PipelineError(
            f"Target partition coverage mismatch. missing={sorted(missing)[:10]} extra={sorted(extra)[:10]}"
        )
    return records


def build_target_manifest_rows(records: Sequence[TargetSubtreeRecord]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for record in records:
        for tip_name in record.target_nonbackbone_tip_names:
            rows.append(
                {
                    "target_subtree_id": record.target_subtree_id,
                    "tip_name": tip_name,
                    "tip_role": "target",
                    "source_node_id": record.target_root_node_id,
                }
            )
    return rows


def build_paml_subtrees(
    tree: Tree,
    backbone_tip_names: Sequence[str],
    target_records: Sequence[TargetSubtreeRecord],
    outgroup_tip: str,
    max_tips: int,
    logger: logging.Logger,
) -> Tuple[List[PamlSubtreeRecord], List[Dict[str, object]]]:
    master_tip_order = get_tip_names_from_tree(tree)
    tip_order_index = {tip_name: index for index, tip_name in enumerate(master_tip_order)}
    backbone_tip_set = set(backbone_tip_names)
    records: List[PamlSubtreeRecord] = []
    manifest_rows: List[Dict[str, object]] = []
    for index, target_record in enumerate(target_records, start=1):
        selected_tip_set = set(backbone_tip_set)
        selected_tip_set.update(target_record.target_nonbackbone_tip_names)
        selected_tip_set.add(outgroup_tip)
        ordered_tip_names = build_ordered_tip_list(tip_order_index, selected_tip_set)
        if len(ordered_tip_names) > max_tips:
            raise PipelineError(
                f"{target_record.target_subtree_id} produces {len(ordered_tip_names)} tips, exceeding max_tips={max_tips}"
            )
        record = PamlSubtreeRecord(
            paml_subtree_id=f"paml_subtree_{index:04d}",
            target_subtree_id=target_record.target_subtree_id,
            target_root_node_id=target_record.target_root_node_id,
            outgroup_tip=outgroup_tip,
            backbone_n_tips=len(backbone_tip_names),
            target_nonbackbone_n_tips=target_record.target_nonbackbone_n_tips,
            total_n_tips=len(ordered_tip_names),
            backbone_tip_names=list(backbone_tip_names),
            target_nonbackbone_tip_names=list(target_record.target_nonbackbone_tip_names),
            total_tip_names=ordered_tip_names,
            tip_hash=compute_tip_hash(ordered_tip_names),
            paml_tree_file=Path("paml_subtrees", f"paml_subtree_{index:04d}.nwk").as_posix(),
        )
        records.append(record)
        for tip_name in record.backbone_tip_names:
            manifest_rows.append(
                {
                    "paml_subtree_id": record.paml_subtree_id,
                    "tip_name": tip_name,
                    "tip_role": "backbone",
                    "target_subtree_id": record.target_subtree_id,
                    "source_node_id": f"TIP::{tip_name}",
                }
            )
        for tip_name in record.target_nonbackbone_tip_names:
            manifest_rows.append(
                {
                    "paml_subtree_id": record.paml_subtree_id,
                    "tip_name": tip_name,
                    "tip_role": "target",
                    "target_subtree_id": record.target_subtree_id,
                    "source_node_id": record.target_root_node_id,
                }
            )
        manifest_rows.append(
            {
                "paml_subtree_id": record.paml_subtree_id,
                "tip_name": outgroup_tip,
                "tip_role": "outgroup",
                "target_subtree_id": record.target_subtree_id,
                "source_node_id": f"TIP::{outgroup_tip}",
            }
        )
        logger.info(
            "%s backbone=%d target=%d total=%d",
            record.paml_subtree_id,
            record.backbone_n_tips,
            record.target_nonbackbone_n_tips,
            record.total_n_tips,
        )
    return records, manifest_rows


def write_clade_tree(clade: Clade, destination: Path) -> None:
    subtree = Tree(root=clone_clade(clade), rooted=True)
    write_tree(subtree, destination)

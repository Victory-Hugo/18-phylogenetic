#!/usr/bin/env python3
"""Select nested backbone sample tiers from mtDNA variation profiles."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Sequence

import numpy as np


STRICT_FULL_LENGTH_RANGE = "1-16569"
AUTO_DEEP_LINEAGE_PREFIXES = [f"L{i}" for i in range(7)] + ["M", "N", "R"]
AUTO_BASAL_MIN_COUNT = 100
DOT = ord(".")
ZERO = ord("0")
TAB = ord("\t")

QC_COLUMNS = ["metric", "value", "note"]
CLEAN_POOL_COLUMNS = [
    "SampleID",
    "Haplogroup",
    "root_haplogroup",
    "level2_haplogroup",
    "deep_lineage_label",
    "Quality",
    "Range",
    "is_in_vcf",
    "is_strict_full_length",
    "dedupe_status",
    "vcf_column_index",
    "call_n",
    "missing_n",
    "missing_rate",
]
SEED_SUMMARY_COLUMNS = [
    "deep_lineage_label",
    "candidate_n",
    "selected_sample_id",
    "selected_haplogroup",
    "selected_quality",
    "selected_missing_rate",
    "seed_rank_global",
]
DISTANCE_SUMMARY_COLUMNS = [
    "tier",
    "selected_n",
    "mean_pairwise_distance_estimate",
    "min_pairwise_distance_estimate",
    "deep_lineage_covered_n",
    "unique_haplogroup_label_n",
    "notes",
]
SELECTION_COLUMNS = [
    "global_selection_rank",
    "tier_first_included",
    "SampleID",
    "Haplogroup",
    "root_haplogroup",
    "level2_haplogroup",
    "deep_lineage_label",
    "Quality",
    "Range",
    "selection_stage",
    "selection_reason",
    "distance_to_selected_min",
    "distance_to_selected_mean",
    "comparable_site_n",
    "missing_rate",
]


@dataclass
class Candidate:
    sample_id: str
    haplogroup: str
    root_haplogroup: str
    level2_haplogroup: str
    deep_lineage_label: str
    quality: float
    range_value: str
    dedupe_status: str
    vcf_column_index: int
    alt_bits: int = 0
    missing_bits: int = 0
    call_n: int = 0
    missing_n: int = 0

    @property
    def missing_rate(self) -> float:
        denominator = self.call_n + self.missing_n
        if denominator == 0:
            return 1.0
        return self.missing_n / denominator


def normalize_whitespace(text: str) -> str:
    return " ".join((text or "").strip().split())


def normalize_range(text: str) -> str:
    return normalize_whitespace(text)


def parse_root_haplogroup(haplogroup: str) -> str:
    text = (haplogroup or "").strip()
    match = re.match(r"^([A-Z]+)", text)
    return match.group(1) if match else text


def parse_level2_haplogroup(haplogroup: str) -> str:
    text = (haplogroup or "").strip()
    root = parse_root_haplogroup(text)
    suffix = text[len(root) :]
    match = re.match(r"^(\d+)", suffix)
    return root + match.group(1) if match else root


def parse_deep_lineage_label(haplogroup: str) -> str:
    return parse_level2_haplogroup(haplogroup)


def conflict_sort_key(row: dict) -> tuple:
    haplogroup = normalize_whitespace(row["Haplogroup"])
    return (-float(row["Quality"]), len(haplogroup), haplogroup, row["SampleID"])


def write_tsv(rows: Iterable[dict], columns: Sequence[str], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_float(value: float | None, digits: int = 6) -> str:
    if value is None or (isinstance(value, float) and (math.isnan(value) or math.isinf(value))):
        return ""
    return f"{value:.{digits}f}"


def read_meta_rows(meta_path: Path) -> List[dict]:
    with meta_path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_basic_info_sources(basic_info_path: Path) -> Dict[str, str]:
    with basic_info_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ID" not in reader.fieldnames or "Source" not in reader.fieldnames:
            raise ValueError(f"Expected ID and Source columns in basic info TSV: {basic_info_path}")
        sources = {}
        for row in reader:
            sample_id = normalize_whitespace(row["ID"])
            if sample_id:
                sources[sample_id] = normalize_whitespace(row["Source"])
        return sources


def read_vcf_sample_columns(vcf_path: Path) -> tuple[List[str], Dict[str, int]]:
    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                sample_names = columns[9:]
                sample_to_column = {sample_id: index for index, sample_id in enumerate(columns)}
                return sample_names, sample_to_column
    raise ValueError(f"Could not find #CHROM header in VCF: {vcf_path}")


def group_rows_by_sample(meta_rows: List[dict]) -> Dict[str, List[dict]]:
    by_id: Dict[str, List[dict]] = defaultdict(list)
    for row in meta_rows:
        by_id[row["SampleID"]].append(row)
    return by_id


def build_clean_candidate_pool(
    meta_rows: List[dict],
    sample_to_vcf_column: Dict[str, int],
    basic_info_sources: Dict[str, str],
    min_quality: float,
) -> tuple[List[Candidate], dict]:
    raw_unique_ids = {normalize_whitespace(row["SampleID"]) for row in meta_rows}
    vcf_samples = set(sample_to_vcf_column)
    overlap_ids = raw_unique_ids & vcf_samples
    missing_ids = raw_unique_ids - vcf_samples

    by_id: Dict[str, List[dict]] = defaultdict(list)
    strict_full_rows = 0
    quality_filtered_rows = 0
    quality_not_one_rows = 0
    missing_basic_info_rows = 0
    llt_source_rows = 0
    for row in meta_rows:
        sample_id = normalize_whitespace(row["SampleID"])
        if normalize_range(row["Range"]) == STRICT_FULL_LENGTH_RANGE:
            strict_full_rows += 1
        if sample_id in vcf_samples and normalize_range(row["Range"]) == STRICT_FULL_LENGTH_RANGE:
            quality = float(row["Quality"])
            if not math.isclose(quality, 1.0, rel_tol=0.0, abs_tol=1e-12):
                quality_not_one_rows += 1
                continue
            source = basic_info_sources.get(sample_id)
            if source is None:
                missing_basic_info_rows += 1
                continue
            if source.upper().startswith("LLT"):
                llt_source_rows += 1
                continue
            if quality >= min_quality:
                normalized_row = dict(row)
                normalized_row["SampleID"] = sample_id
                by_id[sample_id].append(normalized_row)
            else:
                quality_filtered_rows += 1

    candidates: List[Candidate] = []
    identical_duplicate_ids = 0
    conflicting_duplicate_ids = 0

    for sample_id in sorted(by_id):
        rows = by_id[sample_id]
        unique_signatures = {
            (
                normalize_whitespace(row["Haplogroup"]),
                f"{float(row['Quality']):.10f}",
                normalize_range(row["Range"]),
            )
            for row in rows
        }
        if len(rows) == 1:
            chosen = rows[0]
            dedupe_status = "unique"
        elif len(unique_signatures) == 1:
            chosen = rows[0]
            dedupe_status = "deduped_identical"
            identical_duplicate_ids += 1
        else:
            chosen = sorted(rows, key=conflict_sort_key)[0]
            dedupe_status = "deduped_conflict_resolved"
            conflicting_duplicate_ids += 1

        haplogroup = normalize_whitespace(chosen["Haplogroup"])
        candidates.append(
            Candidate(
                sample_id=sample_id,
                haplogroup=haplogroup,
                root_haplogroup=parse_root_haplogroup(haplogroup),
                level2_haplogroup=parse_level2_haplogroup(haplogroup),
                deep_lineage_label=parse_deep_lineage_label(haplogroup),
                quality=float(chosen["Quality"]),
                range_value=normalize_range(chosen["Range"]),
                dedupe_status=dedupe_status,
                vcf_column_index=sample_to_vcf_column[sample_id],
            )
        )

    stats = {
        "raw_row_count": len(meta_rows),
        "raw_unique_sample_ids": len(raw_unique_ids),
        "vcf_overlap_sample_ids": len(overlap_ids),
        "missing_in_vcf_sample_ids": len(missing_ids),
        "duplicate_sample_ids_total": sum(1 for rows in group_rows_by_sample(meta_rows).values() if len(rows) > 1),
        "strict_full_length_rows": strict_full_rows,
        "quality_not_one_rows": quality_not_one_rows,
        "quality_filtered_rows": quality_filtered_rows,
        "missing_basic_info_rows": missing_basic_info_rows,
        "llt_source_rows": llt_source_rows,
        "clean_candidate_pool_size": len(candidates),
        "identical_duplicate_ids": identical_duplicate_ids,
        "conflicting_duplicate_ids": conflicting_duplicate_ids,
    }
    return candidates, stats


def _parse_gt_tail_fast(tail_bytes: bytes, sample_count: int) -> np.ndarray | None:
    tail = tail_bytes.rstrip(b"\n\r")
    flat = np.frombuffer(tail + b"\t", dtype=np.uint8)
    if flat.size != sample_count * 4:
        return None
    matrix = flat.reshape(sample_count, 4)
    if not np.all(matrix[:, 3] == TAB):
        return None
    return matrix[:, :3]


def _apply_variant_bit(candidates: Sequence[Candidate], candidate_indices: np.ndarray, mask: np.ndarray, bit: int, attr: str) -> None:
    if not np.any(mask):
        return
    for candidate_index in candidate_indices[mask].tolist():
        current = getattr(candidates[candidate_index], attr)
        setattr(candidates[candidate_index], attr, current | bit)


def encode_candidate_variants_from_vcf(vcf_path: Path, sample_names: Sequence[str], candidates: Sequence[Candidate]) -> int:
    if not candidates:
        raise ValueError("No candidates remain after cleaning.")

    ordered = sorted(enumerate(candidates), key=lambda item: item[1].vcf_column_index)
    candidate_indices = np.array([item[0] for item in ordered], dtype=np.int64)
    sample_relative_indices = np.array([item[1].vcf_column_index - 9 for item in ordered], dtype=np.int64)
    sample_count = len(sample_names)
    total_variants = 0
    header_seen = False

    with gzip.open(vcf_path, "rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b"##"):
                continue
            if raw_line.startswith(b"#CHROM"):
                header_seen = True
                continue
            if not header_seen:
                raise ValueError(f"VCF header missing before variant rows in {vcf_path}")
            parts = raw_line.rstrip(b"\n").split(b"\t", 9)
            if len(parts) < 10:
                raise ValueError(f"Malformed VCF row in {vcf_path}: {raw_line[:120]!r}")

            total_variants += 1
            bit = 1 << (total_variants - 1)
            format_value = parts[8]
            fast_matrix = _parse_gt_tail_fast(parts[9], sample_count) if format_value == b"GT" else None

            if fast_matrix is not None:
                subset = fast_matrix[sample_relative_indices]
                missing_mask = (subset[:, 0] == DOT) | (subset[:, 2] == DOT)
                alt_mask = (~missing_mask) & ((subset[:, 0] != ZERO) | (subset[:, 2] != ZERO))
                _apply_variant_bit(candidates, candidate_indices, alt_mask, bit, "alt_bits")
                _apply_variant_bit(candidates, candidate_indices, missing_mask, bit, "missing_bits")
                continue

            sample_fields = parts[9].split(b"\t")
            if len(sample_fields) != sample_count:
                raise ValueError(
                    f"Expected {sample_count} sample fields, found {len(sample_fields)} at variant {total_variants} in {vcf_path}"
                )
            for candidate_index, sample_relative_index in zip(candidate_indices.tolist(), sample_relative_indices.tolist()):
                gt = sample_fields[sample_relative_index].split(b":", 1)[0]
                if b"." in gt:
                    candidates[candidate_index].missing_bits |= bit
                    continue
                normalized = gt.replace(b"|", b"/")
                alleles = normalized.split(b"/")
                if any(allele != b"0" for allele in alleles):
                    candidates[candidate_index].alt_bits |= bit

    for candidate in candidates:
        candidate.missing_n = candidate.missing_bits.bit_count()
        candidate.call_n = total_variants - candidate.missing_n
    return total_variants


def infer_deep_lineage_labels(candidates: Sequence[Candidate], label_spec: str) -> List[str]:
    counts = Counter(candidate.deep_lineage_label for candidate in candidates)
    if label_spec and label_spec != "auto":
        labels = []
        for token in label_spec.split(","):
            label = token.strip()
            if label and label in counts and label not in labels:
                labels.append(label)
        return labels

    labels: List[str] = []
    for label in AUTO_DEEP_LINEAGE_PREFIXES:
        if counts[label] > 0:
            labels.append(label)

    basal_candidates = [
        label
        for label, count in counts.items()
        if label == parse_root_haplogroup(label) and count >= AUTO_BASAL_MIN_COUNT and label not in labels
    ]
    basal_candidates.sort(key=lambda label: (-counts[label], label))
    labels.extend(basal_candidates)
    return labels


def compute_distance(candidate_a: Candidate, candidate_b: Candidate, total_variants: int) -> tuple[float, int]:
    excluded = candidate_a.missing_bits | candidate_b.missing_bits
    comparable = total_variants - excluded.bit_count()
    if comparable <= 0:
        return -1.0, 0
    diff = ((candidate_a.alt_bits ^ candidate_b.alt_bits) & ~excluded).bit_count()
    return diff / comparable, comparable


def choose_initial_seed(candidates: Sequence[Candidate]) -> int:
    return min(
        range(len(candidates)),
        key=lambda idx: (-candidates[idx].quality, candidates[idx].missing_n, candidates[idx].sample_id),
    )


def average_distance_within_group(candidate_index: int, group_indices: Sequence[int], candidates: Sequence[Candidate], total_variants: int) -> float:
    distances = []
    for other_index in group_indices:
        if other_index == candidate_index:
            continue
        distance, _ = compute_distance(candidates[candidate_index], candidates[other_index], total_variants)
        if distance >= 0:
            distances.append(distance)
    if not distances:
        return -1.0
    return mean(distances)


def choose_deep_lineage_seed(group_indices: Sequence[int], candidates: Sequence[Candidate], total_variants: int) -> int:
    sorted_indices = sorted(
        group_indices,
        key=lambda idx: (-candidates[idx].quality, candidates[idx].missing_n, candidates[idx].sample_id),
    )
    top = [sorted_indices[0]]
    first = candidates[sorted_indices[0]]
    for candidate_index in sorted_indices[1:]:
        candidate = candidates[candidate_index]
        if candidate.quality == first.quality and candidate.missing_n == first.missing_n:
            top.append(candidate_index)
        else:
            break
    if len(top) == 1:
        return top[0]
    return max(
        top,
        key=lambda idx: (
            average_distance_within_group(idx, group_indices, candidates, total_variants),
            -candidates[idx].missing_n,
            candidates[idx].sample_id,
        ),
    )


def tier_for_rank(rank: int, tiers: Sequence[int]) -> int:
    for tier in tiers:
        if rank <= tier:
            return tier
    raise ValueError(f"Rank {rank} exceeds largest tier {tiers[-1]}")


def compute_mean_distance_to_selected(candidate_index: int, selected_indices: Sequence[int], candidates: Sequence[Candidate], total_variants: int) -> tuple[float | None, int]:
    distances = []
    comparable_sites = 0
    for selected_index in selected_indices:
        distance, comparable = compute_distance(candidates[candidate_index], candidates[selected_index], total_variants)
        if distance >= 0:
            distances.append(distance)
            comparable_sites += comparable
    if not distances:
        return None, 0
    return mean(distances), comparable_sites // len(distances)


def update_min_distance_trackers(
    new_selected_index: int,
    candidates: Sequence[Candidate],
    total_variants: int,
    selected_mask: Sequence[bool],
    current_min_distance: List[float],
    current_min_comparable_sites: List[int],
) -> None:
    selected_candidate = candidates[new_selected_index]
    for candidate_index, candidate in enumerate(candidates):
        if selected_mask[candidate_index]:
            continue
        distance, comparable = compute_distance(candidate, selected_candidate, total_variants)
        if distance < 0:
            continue
        if current_min_distance[candidate_index] < 0 or distance < current_min_distance[candidate_index]:
            current_min_distance[candidate_index] = distance
            current_min_comparable_sites[candidate_index] = comparable


def build_selection_row(
    candidate: Candidate,
    rank: int,
    tiers: Sequence[int],
    selection_stage: str,
    selection_reason: str,
    distance_to_selected_min: float | None,
    distance_to_selected_mean: float | None,
    comparable_site_n: int,
) -> dict:
    return {
        "global_selection_rank": rank,
        "tier_first_included": tier_for_rank(rank, tiers),
        "SampleID": candidate.sample_id,
        "Haplogroup": candidate.haplogroup,
        "root_haplogroup": candidate.root_haplogroup,
        "level2_haplogroup": candidate.level2_haplogroup,
        "deep_lineage_label": candidate.deep_lineage_label,
        "Quality": f"{candidate.quality:.4f}",
        "Range": candidate.range_value,
        "selection_stage": selection_stage,
        "selection_reason": selection_reason,
        "distance_to_selected_min": format_float(distance_to_selected_min),
        "distance_to_selected_mean": format_float(distance_to_selected_mean),
        "comparable_site_n": comparable_site_n,
        "missing_rate": format_float(candidate.missing_rate),
    }


def build_selection(
    candidates: Sequence[Candidate],
    tiers: Sequence[int],
    total_variants: int,
    deep_lineage_labels: Sequence[str],
    seed_mode: str,
) -> tuple[List[dict], Dict[str, dict], List[int]]:
    if seed_mode != "deep_lineage_cover":
        raise ValueError(f"Unsupported --seed-mode: {seed_mode}")
    if not candidates:
        raise ValueError("No candidates remain after cleaning.")

    max_tier = max(tiers)
    if len(candidates) < max_tier:
        raise ValueError(f"Only {len(candidates)} clean candidates remain, fewer than requested tier {max_tier}.")

    selected_indices: List[int] = []
    selected_mask = [False] * len(candidates)
    selection_rows: List[dict] = []
    deep_label_to_seed: Dict[str, dict] = {}
    lineage_selection_counts = Counter()
    current_min_distance = [-1.0] * len(candidates)
    current_min_comparable_sites = [0] * len(candidates)
    label_to_candidate_indices: Dict[str, List[int]] = defaultdict(list)
    for candidate_index, candidate in enumerate(candidates):
        label_to_candidate_indices[candidate.deep_lineage_label].append(candidate_index)

    def add_selected(candidate_index: int, stage: str, reason: str, label_for_seed: str | None = None) -> None:
        rank = len(selected_indices) + 1
        candidate = candidates[candidate_index]
        mean_distance = None
        comparable_sites = 0
        min_distance = None
        if selected_indices:
            min_value = current_min_distance[candidate_index]
            min_distance = min_value if min_value >= 0 else None
            comparable_sites = current_min_comparable_sites[candidate_index]
            mean_distance, _ = compute_mean_distance_to_selected(candidate_index, selected_indices, candidates, total_variants)
        selected_mask[candidate_index] = True
        selected_indices.append(candidate_index)
        lineage_selection_counts[candidate.deep_lineage_label] += 1
        selection_rows.append(
            build_selection_row(
                candidate=candidate,
                rank=rank,
                tiers=tiers,
                selection_stage=stage,
                selection_reason=reason,
                distance_to_selected_min=min_distance,
                distance_to_selected_mean=mean_distance,
                comparable_site_n=comparable_sites,
            )
        )
        if label_for_seed is not None:
            deep_label_to_seed[label_for_seed] = {
                "deep_lineage_label": label_for_seed,
                "candidate_n": len(label_to_candidate_indices[label_for_seed]),
                "selected_sample_id": candidate.sample_id,
                "selected_haplogroup": candidate.haplogroup,
                "selected_quality": f"{candidate.quality:.4f}",
                "selected_missing_rate": format_float(candidate.missing_rate),
                "seed_rank_global": rank,
            }
        update_min_distance_trackers(
            new_selected_index=candidate_index,
            candidates=candidates,
            total_variants=total_variants,
            selected_mask=selected_mask,
            current_min_distance=current_min_distance,
            current_min_comparable_sites=current_min_comparable_sites,
        )

    initial_index = choose_initial_seed(candidates)
    add_selected(initial_index, "seed_cover", "initial_seed_best_quality")
    if candidates[initial_index].deep_lineage_label in deep_lineage_labels:
        label = candidates[initial_index].deep_lineage_label
        deep_label_to_seed[label] = {
            "deep_lineage_label": label,
            "candidate_n": len(label_to_candidate_indices[label]),
            "selected_sample_id": candidates[initial_index].sample_id,
            "selected_haplogroup": candidates[initial_index].haplogroup,
            "selected_quality": f"{candidates[initial_index].quality:.4f}",
            "selected_missing_rate": format_float(candidates[initial_index].missing_rate),
            "seed_rank_global": 1,
        }

    for label in deep_lineage_labels:
        if label in deep_label_to_seed:
            continue
        group_indices = [idx for idx in label_to_candidate_indices.get(label, []) if not selected_mask[idx]]
        if not group_indices:
            continue
        selected_index = choose_deep_lineage_seed(group_indices, candidates, total_variants)
        add_selected(selected_index, "seed_cover", f"deep_lineage_seed:{label}", label_for_seed=label)

    while len(selected_indices) < max_tier:
        best_index = None
        best_score = None
        best_sample_id = None
        for candidate_index, candidate in enumerate(candidates):
            if selected_mask[candidate_index]:
                continue
            distance_score = current_min_distance[candidate_index]
            score = (
                distance_score,
                candidate.quality,
                -candidate.missing_n,
                -lineage_selection_counts[candidate.deep_lineage_label],
            )
            if (
                best_index is None
                or score > best_score
                or (score == best_score and candidate.sample_id < best_sample_id)
            ):
                best_index = candidate_index
                best_score = score
                best_sample_id = candidate.sample_id
        if best_index is None:
            raise ValueError("No candidate available during diversity expansion.")
        add_selected(best_index, "diversity_expansion", "max_min_distance")

    return selection_rows, deep_label_to_seed, selected_indices


def build_clean_rows(candidates: Sequence[Candidate]) -> List[dict]:
    rows = []
    for candidate in sorted(candidates, key=lambda item: (-item.quality, item.sample_id)):
        rows.append(
            {
                "SampleID": candidate.sample_id,
                "Haplogroup": candidate.haplogroup,
                "root_haplogroup": candidate.root_haplogroup,
                "level2_haplogroup": candidate.level2_haplogroup,
                "deep_lineage_label": candidate.deep_lineage_label,
                "Quality": f"{candidate.quality:.4f}",
                "Range": candidate.range_value,
                "is_in_vcf": 1,
                "is_strict_full_length": 1,
                "dedupe_status": candidate.dedupe_status,
                "vcf_column_index": candidate.vcf_column_index,
                "call_n": candidate.call_n,
                "missing_n": candidate.missing_n,
                "missing_rate": format_float(candidate.missing_rate),
            }
        )
    return rows


def build_qc_summary(
    meta_rows: Sequence[dict],
    sample_names: Sequence[str],
    candidates: Sequence[Candidate],
    stats: dict,
    total_variants: int,
    deep_lineage_labels: Sequence[str],
    tiers: Sequence[int],
    distance_mode: str,
    seed_mode: str,
) -> List[dict]:
    rows = [
        {"metric": "raw_row_count", "value": len(meta_rows), "note": "Raw metadata rows before any filtering."},
        {"metric": "raw_unique_sample_ids", "value": stats["raw_unique_sample_ids"], "note": "Unique SampleID values in metadata."},
        {"metric": "vcf_sample_count", "value": len(sample_names), "note": "Total samples present in the VCF header."},
        {"metric": "vcf_overlap_sample_ids", "value": stats["vcf_overlap_sample_ids"], "note": "Metadata SampleIDs present in the VCF header."},
        {"metric": "missing_in_vcf_sample_ids", "value": stats["missing_in_vcf_sample_ids"], "note": "Metadata SampleIDs absent from the VCF header."},
        {"metric": "duplicate_sample_ids", "value": stats["duplicate_sample_ids_total"], "note": "SampleIDs with more than one metadata row."},
        {"metric": "strict_full_length_rows", "value": stats["strict_full_length_rows"], "note": "Rows whose normalized Range equals 1-16569."},
        {"metric": "quality_not_one_rows", "value": stats["quality_not_one_rows"], "note": "Rows removed because Quality was not exactly 1."},
        {"metric": "quality_filtered_rows", "value": stats["quality_filtered_rows"], "note": "Rows removed because Quality was below --min-quality."},
        {"metric": "missing_basic_info_rows", "value": stats["missing_basic_info_rows"], "note": "Rows removed because the SampleID was absent from basic info."},
        {"metric": "llt_source_rows", "value": stats["llt_source_rows"], "note": "Rows removed because Source in basic info started with LLT."},
        {"metric": "identical_duplicate_ids", "value": stats["identical_duplicate_ids"], "note": "Duplicate SampleIDs collapsed because all retained fields matched."},
        {"metric": "conflicting_duplicate_ids", "value": stats["conflicting_duplicate_ids"], "note": "Duplicate SampleIDs resolved by quality and haplogroup tie-break rules."},
        {"metric": "clean_candidate_pool_size", "value": len(candidates), "note": "Unique SampleIDs retained after filtering and deduplication."},
        {"metric": "variant_site_count", "value": total_variants, "note": "Variant rows scanned from the VCF body for distance encoding."},
        {"metric": "distance_mode", "value": distance_mode, "note": "Distance function used during diversity expansion."},
        {"metric": "seed_mode", "value": seed_mode, "note": "Seed selection mode used before diversity expansion."},
        {"metric": "deep_lineage_label_count", "value": len(deep_lineage_labels), "note": "Deep-lineage labels considered during seed coverage."},
        {"metric": "requested_tiers", "value": ",".join(str(tier) for tier in tiers), "note": "Tier cutoffs sliced from the global ordered backbone master list."},
    ]
    return rows


def build_seed_summary(deep_lineage_labels: Sequence[str], deep_label_to_seed: Dict[str, dict]) -> List[dict]:
    rows = []
    for label in deep_lineage_labels:
        rows.append(
            deep_label_to_seed.get(
                label,
                {
                    "deep_lineage_label": label,
                    "candidate_n": 0,
                    "selected_sample_id": "",
                    "selected_haplogroup": "",
                    "selected_quality": "",
                    "selected_missing_rate": "",
                    "seed_rank_global": "",
                },
            )
        )
    return rows


def build_distance_summary(
    tiers: Sequence[int],
    selected_indices: Sequence[int],
    candidates: Sequence[Candidate],
    total_variants: int,
    deep_lineage_labels: Sequence[str],
    distance_mode: str,
    seed_mode: str,
) -> List[dict]:
    rows = []
    selected_deep_label_set = set(deep_lineage_labels)
    for tier in tiers:
        prefix = selected_indices[:tier]
        pairwise_distances = []
        for left_index in range(len(prefix)):
            for right_index in range(left_index + 1, len(prefix)):
                distance, _ = compute_distance(candidates[prefix[left_index]], candidates[prefix[right_index]], total_variants)
                if distance >= 0:
                    pairwise_distances.append(distance)
        deep_lineage_covered_n = len({candidates[idx].deep_lineage_label for idx in prefix if candidates[idx].deep_lineage_label in selected_deep_label_set})
        unique_label_n = len({candidates[idx].deep_lineage_label for idx in prefix})
        rows.append(
            {
                "tier": tier,
                "selected_n": len(prefix),
                "mean_pairwise_distance_estimate": format_float(mean(pairwise_distances) if pairwise_distances else None),
                "min_pairwise_distance_estimate": format_float(min(pairwise_distances) if pairwise_distances else None),
                "deep_lineage_covered_n": deep_lineage_covered_n,
                "unique_haplogroup_label_n": unique_label_n,
                "notes": f"distance_mode={distance_mode};seed_mode={seed_mode}",
            }
        )
    return rows


def build_readme(
    tiers: Sequence[int],
    deep_lineage_labels: Sequence[str],
    clean_candidates: Sequence[Candidate],
    selection_rows: Sequence[dict],
    distance_summary_rows: Sequence[dict],
) -> str:
    seed_rows = [row for row in selection_rows if row["selection_stage"] == "seed_cover"]
    lines = [
        "# 筛选说明",
        "",
        "## 数据清洗规则",
        "",
        "- 仅保留出现在当前 VCF 头部中的 SampleID。",
        "- 仅保留严格全长样本，即 `Range` 规范化后严格等于 `1-16569`。",
        "- 仅保留 `Quality == 1` 的样本。",
        "- 仅保留在基础信息表中存在且 `Source` 不以 `LLT` 开头的样本。",
        "- 对重复 SampleID 去重；完全相同的重复只保留 1 条，冲突重复按 Quality、更浅层 Haplogroup、SampleID 顺序决议。",
        "- 可额外通过 `--min-quality` 再收紧阈值，但在当前规则下必须先满足 `Quality == 1`。",
        "",
        "## 变异距离编码",
        "",
        "- 主流程直接扫描 VCF 正文，而不是只检查 VCF header。",
        "- 默认距离模式为 `alt_hamming`：每个位点仅记录是否携带任一 ALT。",
        "- `0/0` 记为参考型，任意非 0 等位基因组合如 `1/1`、`0/1`、`1/2` 统一记为 ALT presence。",
        "- `./.` 或含 `.` 的 GT 视为缺失，该位点不进入两样本距离分母。",
        "- 当前实现依赖 `numpy` 做 GT 向量解析，但不构建全量两两距离矩阵。",
        "",
        "## 种子覆盖与多样性扩展",
        "",
        "- 第一阶段先加入 1 个全局初始种子样本，优先 Quality 高且缺失少者。",
        "- 随后对关键深分支做最低覆盖约束，默认自动识别 `L0-L6`、`M`、`N`、`R` 及显著 basal 标签。",
        "- 例如 `L0/L1/L2/L3` 会在种子阶段分别考虑，不再像旧版那样被统一压缩进单一 root=`L`。",
        "- 第二阶段使用增量式 farthest-first / max-min diversity 扩展，每轮选择与当前已选集合最远的候选样本。",
        "- `100/150/200/250/300` 五档名单均由同一全局排序列表前缀切片得到，因此严格嵌套。",
        "",
        "## 与旧版差异",
        "",
        "- 旧版以 `root_haplogroup` / `level2_haplogroup` 覆盖为主；新版以真实变异差异为主。",
        "- 旧版只用 VCF header 验证样本存在；新版直接用 VCF 正文编码距离。",
        "- 旧版输出的 `03_root_haplogroup_summary.tsv` 与 `04_level2_haplogroup_summary.tsv` 已由深分支与距离统计表替代。",
        "",
        "## 当前运行概览",
        "",
        f"- 清洗后候选池大小：`{len(clean_candidates)}`",
        f"- 自动/手工识别的关键深分支标签数：`{len(deep_lineage_labels)}`",
        f"- 实际种子样本数：`{len(seed_rows)}`",
        "",
        "## 各档距离统计",
        "",
        "| Tier | 样本数 | 平均成对距离 | 最小成对距离 | 深分支覆盖数 |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in distance_summary_rows:
        lines.append(
            f"| {row['tier']} | {row['selected_n']} | {row['mean_pairwise_distance_estimate'] or 'NA'} | {row['min_pairwise_distance_estimate'] or 'NA'} | {row['deep_lineage_covered_n']} |"
        )
    lines.extend(
        [
            "",
            "## 结果文件",
            "",
            "- `01_input_qc_summary.tsv`：输入、过滤与运行参数统计。",
            "- `02_clean_candidate_pool.tsv`：最终可参与筛选的样本池，包含 VCF 列索引与缺失统计。",
            "- `03_deep_lineage_seed_summary.tsv`：关键深分支种子覆盖结果。",
            "- `04_distance_selection_summary.tsv`：各 tier 的距离与覆盖统计。",
            "- `05_backbone_selection_master.tsv`：到最大 tier 为止的全局有序主名单。",
            "- `backbone_100.tsv` 至 `backbone_300.tsv`：按 tier 切片得到的严格嵌套结果。",
        ]
    )
    return "\n".join(lines) + "\n"


def remove_legacy_outputs(output_dir: Path) -> None:
    for legacy_name in ("03_root_haplogroup_summary.tsv", "04_level2_haplogroup_summary.tsv"):
        legacy_path = output_dir / legacy_name
        if legacy_path.exists():
            legacy_path.unlink()


def parse_tiers(raw_tiers: str) -> List[int]:
    tiers = sorted({int(token.strip()) for token in raw_tiers.split(",") if token.strip()})
    if not tiers:
        raise ValueError("At least one tier must be provided.")
    return tiers


def run(
    meta_path: Path,
    vcf_path: Path,
    basic_info_path: Path,
    output_dir: Path,
    tiers: Sequence[int],
    distance_mode: str,
    seed_mode: str,
    deep_lineage_labels_spec: str,
    min_quality: float,
) -> None:
    if distance_mode != "alt_hamming":
        raise ValueError(f"Unsupported --distance-mode: {distance_mode}")

    meta_rows = read_meta_rows(meta_path)
    basic_info_sources = read_basic_info_sources(basic_info_path)
    sample_names, sample_to_vcf_column = read_vcf_sample_columns(vcf_path)
    candidates, stats = build_clean_candidate_pool(meta_rows, sample_to_vcf_column, basic_info_sources, min_quality)
    total_variants = encode_candidate_variants_from_vcf(vcf_path, sample_names, candidates)
    deep_lineage_labels = infer_deep_lineage_labels(candidates, deep_lineage_labels_spec)
    selection_rows, deep_label_to_seed, selected_indices = build_selection(
        candidates=candidates,
        tiers=tiers,
        total_variants=total_variants,
        deep_lineage_labels=deep_lineage_labels,
        seed_mode=seed_mode,
    )

    qc_summary_rows = build_qc_summary(
        meta_rows=meta_rows,
        sample_names=sample_names,
        candidates=candidates,
        stats=stats,
        total_variants=total_variants,
        deep_lineage_labels=deep_lineage_labels,
        tiers=tiers,
        distance_mode=distance_mode,
        seed_mode=seed_mode,
    )
    clean_rows = build_clean_rows(candidates)
    seed_summary_rows = build_seed_summary(deep_lineage_labels, deep_label_to_seed)
    distance_summary_rows = build_distance_summary(
        tiers=tiers,
        selected_indices=selected_indices,
        candidates=candidates,
        total_variants=total_variants,
        deep_lineage_labels=deep_lineage_labels,
        distance_mode=distance_mode,
        seed_mode=seed_mode,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    remove_legacy_outputs(output_dir)
    write_tsv(qc_summary_rows, QC_COLUMNS, output_dir / "01_input_qc_summary.tsv")
    write_tsv(clean_rows, CLEAN_POOL_COLUMNS, output_dir / "02_clean_candidate_pool.tsv")
    write_tsv(seed_summary_rows, SEED_SUMMARY_COLUMNS, output_dir / "03_deep_lineage_seed_summary.tsv")
    write_tsv(distance_summary_rows, DISTANCE_SUMMARY_COLUMNS, output_dir / "04_distance_selection_summary.tsv")
    write_tsv(selection_rows, SELECTION_COLUMNS, output_dir / "05_backbone_selection_master.tsv")

    for tier in tiers:
        tier_rows = [row for row in selection_rows if int(row["global_selection_rank"]) <= tier]
        write_tsv(tier_rows, SELECTION_COLUMNS, output_dir / f"backbone_{tier}.tsv")

    readme_text = build_readme(
        tiers=tiers,
        deep_lineage_labels=deep_lineage_labels,
        clean_candidates=candidates,
        selection_rows=selection_rows,
        distance_summary_rows=distance_summary_rows,
    )
    (output_dir / "README_筛选说明.md").write_text(readme_text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Select nested backbone sample tiers from mtDNA variation profiles.")
    parser.add_argument("--meta", required=True, help="Input metadata TSV with SampleID/Haplogroup/Quality/Range columns.")
    parser.add_argument("--vcf", required=True, help="Input VCF.gz used to compute sample-level variation distances.")
    parser.add_argument("--basic-info", required=True, help="Basic info TSV with ID and Source columns used for LLT filtering.")
    parser.add_argument("--output-dir", required=True, help="Output directory for selection artifacts.")
    parser.add_argument("--tiers", default="100,150,200,250,300", help="Comma-separated nested tier sizes.")
    parser.add_argument("--distance-mode", default="alt_hamming", help="Distance definition. Currently only alt_hamming is supported.")
    parser.add_argument("--seed-mode", default="deep_lineage_cover", help="Seed selection mode. Currently only deep_lineage_cover is supported.")
    parser.add_argument("--deep-lineage-labels", default="auto", help="Comma-separated deep-lineage labels or 'auto'.")
    parser.add_argument("--min-quality", type=float, default=0.0, help="Drop candidates whose Quality is below this threshold.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    run(
        meta_path=Path(args.meta).resolve(),
        vcf_path=Path(args.vcf).resolve(),
        basic_info_path=Path(args.basic_info).resolve(),
        output_dir=Path(args.output_dir).resolve(),
        tiers=parse_tiers(args.tiers),
        distance_mode=args.distance_mode,
        seed_mode=args.seed_mode,
        deep_lineage_labels_spec=args.deep_lineage_labels,
        min_quality=args.min_quality,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

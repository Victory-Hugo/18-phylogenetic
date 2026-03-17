#!/usr/bin/env python3
"""Select nested backbone sample tiers from mtDNA haplogroup metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


STRICT_FULL_LENGTH_RANGE = "1-16569"
QC_COLUMNS = ["metric", "value", "note"]
CLEAN_POOL_COLUMNS = [
    "SampleID",
    "Haplogroup",
    "root_haplogroup",
    "level2_haplogroup",
    "Quality",
    "Range",
    "is_in_vcf",
    "is_strict_full_length",
    "dedupe_status",
]
ROOT_SUMMARY_COLUMNS = [
    "root_haplogroup",
    "candidate_n",
    "best_quality",
    "selected_in_100",
    "selected_in_150",
    "selected_in_200",
    "selected_in_250",
    "selected_in_300",
]
LEVEL2_SUMMARY_COLUMNS = [
    "root_haplogroup",
    "level2_haplogroup",
    "candidate_n",
    "best_quality",
    "representative_sample_id",
    "selected_rank_global",
    "selected_in_100",
    "selected_in_150",
    "selected_in_200",
    "selected_in_250",
    "selected_in_300",
]
SELECTION_COLUMNS = [
    "global_selection_rank",
    "tier_first_included",
    "SampleID",
    "Haplogroup",
    "root_haplogroup",
    "level2_haplogroup",
    "Quality",
    "Range",
    "selection_stage",
    "selection_reason",
]


@dataclass(frozen=True)
class Candidate:
    sample_id: str
    haplogroup: str
    root_haplogroup: str
    level2_haplogroup: str
    quality: float
    range_value: str
    dedupe_status: str


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


def representative_sort_key(candidate: Candidate, group_label: str) -> tuple:
    exact_penalty = 0 if candidate.haplogroup == group_label else 1
    suffix_length = len(candidate.haplogroup[len(group_label) :]) if candidate.haplogroup.startswith(group_label) else len(candidate.haplogroup)
    return (exact_penalty, suffix_length, -candidate.quality, candidate.sample_id)


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


def read_vcf_samples(vcf_path: Path) -> set[str]:
    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return set(line.rstrip("\n").split("\t")[9:])
    raise ValueError(f"Could not find #CHROM header in VCF: {vcf_path}")


def read_meta_rows(meta_path: Path) -> List[dict]:
    with meta_path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_clean_candidate_pool(meta_rows: List[dict], vcf_samples: set[str]) -> tuple[List[Candidate], List[dict], dict]:
    raw_unique_ids = {row["SampleID"] for row in meta_rows}
    overlap_ids = raw_unique_ids & vcf_samples
    missing_ids = raw_unique_ids - vcf_samples

    strict_full_rows = []
    for row in meta_rows:
        if normalize_range(row["Range"]) == STRICT_FULL_LENGTH_RANGE:
            strict_full_rows.append(row)

    by_id: Dict[str, List[dict]] = defaultdict(list)
    for row in meta_rows:
        sample_id = row["SampleID"]
        if sample_id in vcf_samples and normalize_range(row["Range"]) == STRICT_FULL_LENGTH_RANGE:
            by_id[sample_id].append(row)

    candidates: List[Candidate] = []
    clean_rows: List[dict] = []
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
        candidate = Candidate(
            sample_id=sample_id,
            haplogroup=haplogroup,
            root_haplogroup=parse_root_haplogroup(haplogroup),
            level2_haplogroup=parse_level2_haplogroup(haplogroup),
            quality=float(chosen["Quality"]),
            range_value=normalize_range(chosen["Range"]),
            dedupe_status=dedupe_status,
        )
        candidates.append(candidate)
        clean_rows.append(
            {
                "SampleID": candidate.sample_id,
                "Haplogroup": candidate.haplogroup,
                "root_haplogroup": candidate.root_haplogroup,
                "level2_haplogroup": candidate.level2_haplogroup,
                "Quality": f"{candidate.quality:.4f}",
                "Range": candidate.range_value,
                "is_in_vcf": 1,
                "is_strict_full_length": 1,
                "dedupe_status": candidate.dedupe_status,
            }
        )

    stats = {
        "raw_row_count": len(meta_rows),
        "raw_unique_sample_ids": len(raw_unique_ids),
        "vcf_overlap_sample_ids": len(overlap_ids),
        "missing_in_vcf_sample_ids": len(missing_ids),
        "duplicate_sample_ids": sum(1 for sample_id, rows in defaultdict(list, ((row["SampleID"], []) for row in [])).items()),
        "strict_full_length_rows": len(strict_full_rows),
        "clean_candidate_pool_size": len(candidates),
        "identical_duplicate_ids": identical_duplicate_ids,
        "conflicting_duplicate_ids": conflicting_duplicate_ids,
        "duplicate_sample_ids_total": sum(1 for rows in _group_rows_by_sample(meta_rows).values() if len(rows) > 1),
    }
    return candidates, clean_rows, stats


def _group_rows_by_sample(meta_rows: List[dict]) -> Dict[str, List[dict]]:
    by_id: Dict[str, List[dict]] = defaultdict(list)
    for row in meta_rows:
        by_id[row["SampleID"]].append(row)
    return by_id


def build_group_maps(candidates: Sequence[Candidate]) -> tuple[dict, dict]:
    root_groups: Dict[str, List[Candidate]] = defaultdict(list)
    level2_groups: Dict[str, List[Candidate]] = defaultdict(list)
    for candidate in candidates:
        root_groups[candidate.root_haplogroup].append(candidate)
        level2_groups[candidate.level2_haplogroup].append(candidate)
    return root_groups, level2_groups


def choose_representatives(candidates: Sequence[Candidate], group_key: str) -> Dict[str, Candidate]:
    grouped: Dict[str, List[Candidate]] = defaultdict(list)
    for candidate in candidates:
        label = getattr(candidate, group_key)
        grouped[label].append(candidate)
    representatives = {}
    for label, rows in grouped.items():
        representatives[label] = sorted(rows, key=lambda row: representative_sort_key(row, label))[0]
    return representatives


def tier_for_rank(rank: int, tiers: Sequence[int]) -> int:
    for tier in tiers:
        if rank <= tier:
            return tier
    raise ValueError(f"Rank {rank} exceeds largest tier {tiers[-1]}")


def build_selection(
    candidates: Sequence[Candidate],
    tiers: Sequence[int],
) -> tuple[List[dict], Dict[str, int], Dict[str, Candidate], Dict[str, Candidate], Dict[str, int]]:
    if not candidates:
        raise ValueError("No candidates remain after cleaning.")

    max_tier = max(tiers)
    root_groups, level2_groups = build_group_maps(candidates)
    root_reps = choose_representatives(candidates, "root_haplogroup")
    level2_reps = choose_representatives(candidates, "level2_haplogroup")

    root_order = sorted(root_groups, key=lambda root: (-len(root_groups[root]), root))
    selection_rows: List[dict] = []
    selected_ids: set[str] = set()
    selected_rank_by_level2: Dict[str, int] = {}

    for root in root_order:
        candidate = root_reps[root]
        if candidate.sample_id in selected_ids:
            continue
        rank = len(selection_rows) + 1
        selected_ids.add(candidate.sample_id)
        selected_rank_by_level2[candidate.level2_haplogroup] = rank
        selection_rows.append(
            {
                "global_selection_rank": rank,
                "tier_first_included": tier_for_rank(rank, tiers),
                "SampleID": candidate.sample_id,
                "Haplogroup": candidate.haplogroup,
                "root_haplogroup": candidate.root_haplogroup,
                "level2_haplogroup": candidate.level2_haplogroup,
                "Quality": f"{candidate.quality:.4f}",
                "Range": candidate.range_value,
                "selection_stage": "root_coverage",
                "selection_reason": f"best root representative for {root}",
            }
        )

    root_to_level2_labels: Dict[str, List[str]] = defaultdict(list)
    for level2_label, rows in level2_groups.items():
        root_to_level2_labels[rows[0].root_haplogroup].append(level2_label)

    root_level2_queues: Dict[str, List[str]] = {}
    for root in root_order:
        labels = sorted(
            root_to_level2_labels[root],
            key=lambda label: (-len(level2_groups[label]), label),
        )
        queue = []
        for label in labels:
            representative = level2_reps[label]
            if representative.sample_id in selected_ids:
                continue
            queue.append(label)
        root_level2_queues[root] = queue

    while len(selection_rows) < max_tier:
        progress = False
        for root in root_order:
            queue = root_level2_queues[root]
            if not queue:
                continue
            label = queue.pop(0)
            representative = level2_reps[label]
            if representative.sample_id in selected_ids:
                continue
            rank = len(selection_rows) + 1
            selected_ids.add(representative.sample_id)
            selected_rank_by_level2[label] = rank
            selection_rows.append(
                {
                    "global_selection_rank": rank,
                    "tier_first_included": tier_for_rank(rank, tiers),
                    "SampleID": representative.sample_id,
                    "Haplogroup": representative.haplogroup,
                    "root_haplogroup": representative.root_haplogroup,
                    "level2_haplogroup": representative.level2_haplogroup,
                    "Quality": f"{representative.quality:.4f}",
                    "Range": representative.range_value,
                    "selection_stage": "level2_coverage",
                    "selection_reason": f"first {label} representative added during level2 breadth expansion",
                }
            )
            progress = True
            if len(selection_rows) >= max_tier:
                break
        if not progress:
            break

    if len(selection_rows) < max_tier:
        raise ValueError(f"Could only select {len(selection_rows)} samples, fewer than required max tier {max_tier}.")

    selected_rank_by_sample = {row["SampleID"]: int(row["global_selection_rank"]) for row in selection_rows}
    return selection_rows, selected_rank_by_sample, root_reps, level2_reps, selected_rank_by_level2


def build_root_summary(
    root_groups: Dict[str, List[Candidate]],
    selected_rank_by_sample: Dict[str, int],
    tiers: Sequence[int],
) -> List[dict]:
    rows = []
    for root in sorted(root_groups, key=lambda label: (-len(root_groups[label]), label)):
        candidates = root_groups[root]
        best_quality = max(candidate.quality for candidate in candidates)
        sample_ids = {candidate.sample_id for candidate in candidates}
        row = {
            "root_haplogroup": root,
            "candidate_n": len(candidates),
            "best_quality": f"{best_quality:.4f}",
        }
        for tier in tiers:
            row[f"selected_in_{tier}"] = 1 if any(selected_rank_by_sample.get(sample_id, tier + 1) <= tier for sample_id in sample_ids) else 0
        rows.append(row)
    return rows


def build_level2_summary(
    level2_groups: Dict[str, List[Candidate]],
    level2_reps: Dict[str, Candidate],
    selected_rank_by_level2: Dict[str, int],
    tiers: Sequence[int],
) -> List[dict]:
    rows = []
    for label in sorted(level2_groups, key=lambda key: (-len(level2_groups[key]), key)):
        candidates = level2_groups[label]
        best_quality = max(candidate.quality for candidate in candidates)
        representative = level2_reps[label]
        selected_rank = selected_rank_by_level2.get(label, "")
        row = {
            "root_haplogroup": representative.root_haplogroup,
            "level2_haplogroup": label,
            "candidate_n": len(candidates),
            "best_quality": f"{best_quality:.4f}",
            "representative_sample_id": representative.sample_id,
            "selected_rank_global": selected_rank,
        }
        for tier in tiers:
            row[f"selected_in_{tier}"] = 1 if selected_rank and selected_rank <= tier else 0
        rows.append(row)
    return rows


def build_qc_summary(meta_rows: List[dict], vcf_samples: set[str], clean_candidates: Sequence[Candidate], stats: dict) -> List[dict]:
    missing_ids = {row["SampleID"] for row in meta_rows if row["SampleID"] not in vcf_samples}
    strict_full_rows = sum(1 for row in meta_rows if normalize_range(row["Range"]) == STRICT_FULL_LENGTH_RANGE)
    return [
        {"metric": "raw_row_count", "value": len(meta_rows), "note": "Raw metadata rows before any filtering."},
        {"metric": "raw_unique_sample_ids", "value": len({row['SampleID'] for row in meta_rows}), "note": "Unique SampleID values in metadata."},
        {"metric": "vcf_overlap_sample_ids", "value": len({row['SampleID'] for row in meta_rows} & vcf_samples), "note": "Metadata SampleIDs present in the VCF header."},
        {"metric": "missing_in_vcf_sample_ids", "value": len(missing_ids), "note": "Metadata SampleIDs absent from the VCF header."},
        {"metric": "duplicate_sample_ids", "value": stats["duplicate_sample_ids_total"], "note": "SampleIDs with more than one metadata row."},
        {"metric": "strict_full_length_rows", "value": strict_full_rows, "note": "Rows whose normalized Range equals 1-16569."},
        {"metric": "identical_duplicate_ids", "value": stats["identical_duplicate_ids"], "note": "Duplicate SampleIDs collapsed because all retained fields matched."},
        {"metric": "conflicting_duplicate_ids", "value": stats["conflicting_duplicate_ids"], "note": "Duplicate SampleIDs resolved by quality/shallowness tie-break rules."},
        {"metric": "clean_candidate_pool_size", "value": len(clean_candidates), "note": "Unique SampleIDs retained after VCF, strict full-length, and deduplication filters."},
    ]


def build_readme(
    tiers: Sequence[int],
    clean_candidates: Sequence[Candidate],
    root_summary: Sequence[dict],
    level2_summary: Sequence[dict],
) -> str:
    root_coverages = {
        tier: sum(int(row[f"selected_in_{tier}"]) for row in root_summary)
        for tier in tiers
    }
    level2_coverages = {
        tier: sum(int(row[f"selected_in_{tier}"]) for row in level2_summary)
        for tier in tiers
    }
    lines = [
        "# 筛选说明",
        "",
        "## 数据清洗规则",
        "",
        "- 仅保留出现在当前 VCF 头部中的 SampleID。",
        "- 仅保留严格全长样本，即 `Range` 规范化后严格等于 `1-16569`。",
        "- 对重复 SampleID 去重；完全相同的重复只保留 1 条，冲突重复按 Quality、更浅层 Haplogroup、SampleID 顺序决议。",
        "",
        "## 单倍群层级解析",
        "",
        "- `root_haplogroup`：从 Haplogroup 左侧提取连续大写字母。",
        "- `level2_haplogroup`：`root_haplogroup` 加 root 后首段数字；若无数字，则等于 root。",
        "",
        "## 嵌套策略",
        "",
        "- 第一阶段先覆盖每个 root haplogroup 的 1 个代表样本。",
        "- 第二阶段按 root 候选规模轮转，为新的 level-2 haplogroup 逐步补充代表样本。",
        "- `100/150/200/250/300` 五档名单均由同一全局排序列表前缀切片得到，因此严格嵌套。",
        "",
        "## 候选池规模",
        "",
        f"- 清洗后可参与筛选的唯一样本数：`{len(clean_candidates)}`",
        f"- 清洗后 root haplogroup 数：`{len(root_summary)}`",
        f"- 清洗后 level-2 haplogroup 数：`{len(level2_summary)}`",
        "",
        "## 各档覆盖统计",
        "",
        "| Tier | 样本数 | root 覆盖 | level-2 覆盖 |",
        "| --- | ---: | ---: | ---: |",
    ]
    for tier in tiers:
        lines.append(f"| {tier} | {tier} | {root_coverages[tier]} | {level2_coverages[tier]} |")
    lines.extend(
        [
            "",
            "## 结果文件",
            "",
            "- `01_input_qc_summary.tsv`：输入与清洗统计。",
            "- `02_clean_candidate_pool.tsv`：最终可参与筛选的样本池。",
            "- `03_root_haplogroup_summary.tsv`：root 层级覆盖统计。",
            "- `04_level2_haplogroup_summary.tsv`：level-2 层级覆盖统计。",
            "- `05_backbone_selection_master.tsv`：到 300 档为止的全局有序主名单。",
            "- `backbone_100.tsv` 至 `backbone_300.tsv`：五档嵌套结果。",
        ]
    )
    return "\n".join(lines) + "\n"


def run(meta_path: Path, vcf_path: Path, output_dir: Path, tiers: Sequence[int]) -> None:
    meta_rows = read_meta_rows(meta_path)
    vcf_samples = read_vcf_samples(vcf_path)
    clean_candidates, clean_rows, stats = build_clean_candidate_pool(meta_rows, vcf_samples)
    if len(clean_candidates) < max(tiers):
        raise ValueError(f"Only {len(clean_candidates)} clean candidates remain, fewer than requested tier {max(tiers)}.")

    root_groups, level2_groups = build_group_maps(clean_candidates)
    selection_rows, selected_rank_by_sample, root_reps, level2_reps, selected_rank_by_level2 = build_selection(clean_candidates, tiers)
    root_summary_rows = build_root_summary(root_groups, selected_rank_by_sample, tiers)
    level2_summary_rows = build_level2_summary(level2_groups, level2_reps, selected_rank_by_level2, tiers)
    qc_summary_rows = build_qc_summary(meta_rows, vcf_samples, clean_candidates, stats)

    write_tsv(qc_summary_rows, QC_COLUMNS, output_dir / "01_input_qc_summary.tsv")
    write_tsv(clean_rows, CLEAN_POOL_COLUMNS, output_dir / "02_clean_candidate_pool.tsv")
    write_tsv(root_summary_rows, ROOT_SUMMARY_COLUMNS, output_dir / "03_root_haplogroup_summary.tsv")
    write_tsv(level2_summary_rows, LEVEL2_SUMMARY_COLUMNS, output_dir / "04_level2_haplogroup_summary.tsv")
    write_tsv(selection_rows, SELECTION_COLUMNS, output_dir / "05_backbone_selection_master.tsv")

    for tier in tiers:
        tier_rows = [row for row in selection_rows if int(row["global_selection_rank"]) <= tier]
        write_tsv(tier_rows, SELECTION_COLUMNS, output_dir / f"backbone_{tier}.tsv")

    readme_text = build_readme(tiers, clean_candidates, root_summary_rows, level2_summary_rows)
    (output_dir / "README_筛选说明.md").write_text(readme_text, encoding="utf-8")


def parse_tiers(raw_tiers: str) -> List[int]:
    tiers = sorted({int(token.strip()) for token in raw_tiers.split(",") if token.strip()})
    if not tiers:
        raise ValueError("At least one tier must be provided.")
    return tiers


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Select nested backbone sample tiers from mtDNA haplogroup metadata.")
    parser.add_argument("--meta", required=True, help="Input metadata TSV with SampleID/Haplogroup/Quality/Range columns.")
    parser.add_argument("--vcf", required=True, help="Input VCF.gz used to verify sample availability.")
    parser.add_argument("--output-dir", required=True, help="Output directory for selection artifacts.")
    parser.add_argument("--tiers", default="100,150,200,250,300", help="Comma-separated nested tier sizes.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    run(
        meta_path=Path(args.meta).resolve(),
        vcf_path=Path(args.vcf).resolve(),
        output_dir=Path(args.output_dir).resolve(),
        tiers=parse_tiers(args.tiers),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

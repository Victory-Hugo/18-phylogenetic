"""
s11_compare_literature_ages.py
统一解析 download/ 中的文献单倍群年龄表，并与本项目 Stage 2 的
ml_vs_rho_comparison.tsv 逐单倍群比较。
"""

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


NA = ""


def clean_haplogroup(value: str) -> str:
    return value.strip().replace("’", "'").replace("`", "'")


def to_float(value: str) -> float | None:
    if value is None:
        return None
    value = value.strip().replace(",", "")
    if value in {"", "-", "NA", "nan"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_range_years(value: str) -> tuple[float | None, float | None, float | None]:
    value = value.strip().replace(",", "").replace("—", "-").replace("–", "-")
    match = re.search(r"(-?\d+(?:\.\d+)?)\s*-\s*(-?\d+(?:\.\d+)?)", value)
    if not match:
        return None, None, None
    lo = float(match.group(1)) / 1000.0
    hi = float(match.group(2)) / 1000.0
    return (lo + hi) / 2.0, lo, hi


def parse_ci_kya(value: str) -> tuple[float | None, float | None, float | None]:
    value = value.strip()
    match = re.search(
        r"(-?\d+(?:\.\d+)?)\s*\(\s*(-?\d+(?:\.\d+)?)\s*[,;]\s*(-?\d+(?:\.\d+)?)\s*\)",
        value,
    )
    if not match:
        return None, None, None
    return tuple(float(match.group(i)) for i in range(1, 4))


def parse_ci_years(value: str) -> tuple[float | None, float | None, float | None]:
    age, lo, hi = parse_ci_kya(value.replace(",", ""))
    if age is None:
        return None, None, None
    return age / 1000.0, lo / 1000.0, hi / 1000.0


def parse_mean_minmax_kya(value: str) -> tuple[float | None, float | None, float | None]:
    value = value.strip()
    match = re.search(
        r"(-?\d+(?:\.\d+)?)\s*\(\s*(-?\d+(?:\.\d+)?)\s*;\s*(-?\d+(?:\.\d+)?)\s*\)",
        value,
    )
    if not match:
        return None, None, None
    return tuple(float(match.group(i)) for i in range(1, 4))


def row_template(source_file: str, reference: str, method: str, haplogroup: str) -> dict:
    return {
        "source_file": source_file,
        "reference": reference,
        "literature_method": method,
        "haplogroup": clean_haplogroup(haplogroup),
        "n_literature": NA,
        "literature_age_kya": NA,
        "literature_ci95_lower_kya": NA,
        "literature_ci95_upper_kya": NA,
        "literature_sd_kya": NA,
        "literature_period": NA,
    }


def add_age(row: dict, age, lo=None, hi=None, sd=None) -> dict | None:
    if age is None:
        return None
    row["literature_age_kya"] = round(float(age), 6)
    if lo is not None:
        row["literature_ci95_lower_kya"] = round(float(lo), 6)
    if hi is not None:
        row["literature_ci95_upper_kya"] = round(float(hi), 6)
    if sd is not None:
        row["literature_sd_kya"] = round(float(sd), 6)
    return row


def parse_study_1(path: Path) -> list[dict]:
    rows = []
    reference = ""
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("# Reference:"):
                reference = line.split(":", 1)[1].strip()
                continue
            if not line or line.startswith("Haplogroup"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            hap = parts[0]
            age, lo, hi = parse_range_years(parts[1])
            row = add_age(row_template(path.name, reference, "ML_range_midpoint", hap), age, lo, hi)
            if row:
                rows.append(row)
            age, lo, hi = parse_ci_years(f"{parts[2]} {parts[3]}")
            row = add_age(row_template(path.name, reference, "MCMC", hap), age, lo, hi)
            if row:
                rows.append(row)
    return rows


def parse_study_2(path: Path) -> list[dict]:
    rows = []
    reference = ""
    methods = [
        ("coding_region_rho", 3),
        ("complete_genome_rho", 6),
        ("synonymous_positions_rho", 9),
        ("bayesian", 12),
    ]
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("# Reference:"):
                reference = line.split(":", 1)[1].strip()
                continue
            parts = line.split("\t")
            if len(parts) < 14:
                continue
            hap = clean_haplogroup(parts[0])
            if not hap or hap in {"Haplogroup", "Coalescent ages based on all publised whole mtDNA sequences"}:
                continue
            n = to_float(parts[1])
            for method, idx in methods:
                age, lo, hi = parse_mean_minmax_kya(parts[idx])
                row = add_age(row_template(path.name, reference, method, hap), age, lo, hi)
                if row:
                    row["n_literature"] = int(n) if n is not None else NA
                    row["literature_period"] = parts[14].strip() if len(parts) > 14 else NA
                    rows.append(row)
            avg = to_float(parts[13])
            row = add_age(row_template(path.name, reference, "average", hap), avg)
            if row:
                row["n_literature"] = int(n) if n is not None else NA
                row["literature_period"] = parts[14].strip() if len(parts) > 14 else NA
                rows.append(row)
    return rows


def parse_study_3(path: Path) -> list[dict]:
    rows = []
    reference = ""
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("# Reference:"):
                reference = line.split(":", 1)[1].strip()
                continue
            if not line or line.startswith("HG") or line.startswith("\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            hap = parts[0]
            for method, value in [("rho_based", parts[1]), ("ML", parts[2])]:
                age, lo, hi = parse_ci_kya(value)
                row = add_age(row_template(path.name, reference, method, hap), age, lo, hi)
                if row:
                    rows.append(row)
    return rows


def parse_study_4(path: Path) -> list[dict]:
    rows = []
    reference = ""
    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if not parts:
                continue
            if parts[0].startswith("# Reference:"):
                reference = parts[0].split(":", 1)[1].strip()
                continue
            if parts[0] in {"Haplogroup", ""} or len(parts) < 3:
                continue
            age_years = to_float(parts[1])
            sd_years = to_float(parts[2])
            if age_years is None:
                continue
            age = age_years / 1000.0
            sd = sd_years / 1000.0 if sd_years is not None else None
            lo = max(0.0, age - 1.96 * sd) if sd is not None else None
            hi = age + 1.96 * sd if sd is not None else None
            row = add_age(row_template(path.name, reference, "purifying_selection_corrected", parts[0]), age, lo, hi, sd)
            if row:
                rows.append(row)
    return rows


def load_literature(download_dir: Path) -> list[dict]:
    parsers = {
        "以往研究-1.tsv": parse_study_1,
        "以往研究-2.tsv": parse_study_2,
        "以往研究-3.tsv": parse_study_3,
        "以往研究-4.tsv": parse_study_4,
    }
    rows = []
    for name, parser in parsers.items():
        path = download_dir / name
        if path.exists():
            rows.extend(parser(path))
    return rows


def load_project(path: Path) -> dict[str, dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {clean_haplogroup(row["haplogroup"]): row for row in reader}


def fnum(value) -> float | None:
    if value in {None, ""}:
        return None
    return float(value)


def interval_overlap(lo_a: float | None, hi_a: float | None, lo_b: float | None, hi_b: float | None) -> str:
    """Return yes/no/NA for whether two closed intervals overlap."""
    if None in {lo_a, hi_a, lo_b, hi_b}:
        return NA
    return "yes" if max(lo_a, lo_b) <= min(hi_a, hi_b) else "no"


def compare_rows(literature_rows: list[dict], project_rows: dict[str, dict], large_diff: float, large_rel: float) -> list[dict]:
    out = []
    for lit in literature_rows:
        hap = lit["haplogroup"]
        project = project_rows.get(hap)
        matched = project is not None
        lit_age = fnum(lit["literature_age_kya"])
        lit_lo = fnum(lit["literature_ci95_lower_kya"])
        lit_hi = fnum(lit["literature_ci95_upper_kya"])
        rho_age = fnum(project.get("rho_age_kya")) if project else None
        rho_lo = fnum(project.get("rho_ci95_lower_kya")) if project else None
        rho_hi = fnum(project.get("rho_ci95_upper_kya")) if project else None
        ml_age = fnum(project.get("ml_age_kya")) if project else None
        ml_lo = fnum(project.get("ml_ci95_lower_kya")) if project else None
        ml_hi = fnum(project.get("ml_ci95_upper_kya")) if project else None
        diff_rho = lit_age - rho_age if matched and lit_age is not None and rho_age is not None else None
        diff_ml = lit_age - ml_age if matched and lit_age is not None and ml_age is not None else None
        rel_rho = abs(diff_rho) / rho_age * 100.0 if diff_rho is not None and rho_age not in {None, 0.0} else None
        rel_ml = abs(diff_ml) / ml_age * 100.0 if diff_ml is not None and ml_age not in {None, 0.0} else None
        nearest = "rho" if diff_rho is not None and (diff_ml is None or abs(diff_rho) <= abs(diff_ml)) else "ml"
        min_abs = min([abs(x) for x in [diff_rho, diff_ml] if x is not None], default=None)
        large_flag = False
        if diff_rho is not None and abs(diff_rho) >= large_diff:
            large_flag = True
        if diff_ml is not None and abs(diff_ml) >= large_diff:
            large_flag = True
        if rel_rho is not None and rel_rho >= large_rel:
            large_flag = True
        if rel_ml is not None and rel_ml >= large_rel:
            large_flag = True
        overlap_rho = interval_overlap(lit_lo, lit_hi, rho_lo, rho_hi) if matched else NA
        overlap_ml = interval_overlap(lit_lo, lit_hi, ml_lo, ml_hi) if matched else NA
        if not matched:
            reasonable_flag = NA
        elif overlap_rho == "yes" or overlap_ml == "yes":
            reasonable_flag = "yes"
        elif overlap_rho == "no" or overlap_ml == "no":
            reasonable_flag = "no"
        else:
            reasonable_flag = NA
        row = dict(lit)
        row.update({
            "matched_in_project": "yes" if matched else "no",
            "project_n_samples": project.get("n_samples", NA) if project else NA,
            "project_confidence": project.get("confidence", NA) if project else NA,
            "project_calibration_method": project.get("calibration_method", NA) if project else NA,
            "project_rho_age_kya": project.get("rho_age_kya", NA) if project else NA,
            "project_rho_ci95_lower_kya": project.get("rho_ci95_lower_kya", NA) if project else NA,
            "project_rho_ci95_upper_kya": project.get("rho_ci95_upper_kya", NA) if project else NA,
            "project_ml_age_kya": project.get("ml_age_kya", NA) if project else NA,
            "project_ml_ci95_lower_kya": project.get("ml_ci95_lower_kya", NA) if project else NA,
            "project_ml_ci95_upper_kya": project.get("ml_ci95_upper_kya", NA) if project else NA,
            "diff_lit_minus_project_rho_kya": round(diff_rho, 6) if diff_rho is not None else NA,
            "diff_lit_minus_project_ml_kya": round(diff_ml, 6) if diff_ml is not None else NA,
            "abs_diff_vs_rho_kya": round(abs(diff_rho), 6) if diff_rho is not None else NA,
            "abs_diff_vs_ml_kya": round(abs(diff_ml), 6) if diff_ml is not None else NA,
            "rel_diff_vs_rho_pct": round(rel_rho, 6) if rel_rho is not None else NA,
            "rel_diff_vs_ml_pct": round(rel_ml, 6) if rel_ml is not None else NA,
            "nearest_project_method": nearest if matched else NA,
            "min_abs_diff_to_project_kya": round(min_abs, 6) if min_abs is not None else NA,
            "overlap_with_project_rho_ci95": overlap_rho,
            "overlap_with_project_ml_ci95": overlap_ml,
            "interval_reasonable_flag": reasonable_flag,
            "large_difference_flag": "yes" if matched and large_flag else ("no" if matched else NA),
        })
        out.append(row)
    return out


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def median(values: list[float]) -> float:
    if not values:
        return math.nan
    values = sorted(values)
    n = len(values)
    mid = n // 2
    if n % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2.0


def summarize(comparison_rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in comparison_rows:
        grouped[(row["source_file"], row["reference"], row["literature_method"])].append(row)
    summary = []
    for (source_file, reference, method), rows in sorted(grouped.items()):
        matched = [r for r in rows if r["matched_in_project"] == "yes"]
        abs_rho = [float(r["abs_diff_vs_rho_kya"]) for r in matched if r["abs_diff_vs_rho_kya"] != NA]
        abs_ml = [float(r["abs_diff_vs_ml_kya"]) for r in matched if r["abs_diff_vs_ml_kya"] != NA]
        large = sum(1 for r in matched if r["large_difference_flag"] == "yes")
        interval_evaluable = [r for r in matched if r["interval_reasonable_flag"] != NA]
        reasonable = sum(1 for r in interval_evaluable if r["interval_reasonable_flag"] == "yes")
        non_overlap = sum(1 for r in interval_evaluable if r["interval_reasonable_flag"] == "no")
        nearest_counts = Counter(r["nearest_project_method"] for r in matched if r["nearest_project_method"])
        summary.append({
            "source_file": source_file,
            "reference": reference,
            "literature_method": method,
            "n_literature_records": len(rows),
            "n_matched_project": len(matched),
            "match_rate_pct": round(len(matched) / len(rows) * 100.0, 2) if rows else 0,
            "mean_abs_diff_vs_rho_kya": round(mean(abs_rho), 4) if abs_rho else NA,
            "median_abs_diff_vs_rho_kya": round(median(abs_rho), 4) if abs_rho else NA,
            "mean_abs_diff_vs_ml_kya": round(mean(abs_ml), 4) if abs_ml else NA,
            "median_abs_diff_vs_ml_kya": round(median(abs_ml), 4) if abs_ml else NA,
            "n_interval_evaluable": len(interval_evaluable),
            "n_interval_reasonable": reasonable,
            "pct_interval_reasonable": round(reasonable / len(interval_evaluable) * 100.0, 2) if interval_evaluable else NA,
            "n_interval_non_overlap": non_overlap,
            "n_large_difference": large,
            "pct_large_difference": round(large / len(matched) * 100.0, 2) if matched else NA,
            "nearest_rho_count": nearest_counts.get("rho", 0),
            "nearest_ml_count": nearest_counts.get("ml", 0),
        })
    return summary


def top_differences(comparison_rows: list[dict], limit: int = 30) -> list[dict]:
    matched = [r for r in comparison_rows if r["matched_in_project"] == "yes" and r["min_abs_diff_to_project_kya"] != NA]
    return sorted(matched, key=lambda r: float(r["min_abs_diff_to_project_kya"]), reverse=True)[:limit]


def write_report(path: Path, literature_rows: list[dict], comparison_rows: list[dict], summary_rows: list[dict], top_rows: list[dict], project_rows: dict[str, dict]) -> None:
    matched_haps = {r["haplogroup"] for r in comparison_rows if r["matched_in_project"] == "yes"}
    literature_haps = {r["haplogroup"] for r in literature_rows}
    interval_evaluable = [r for r in comparison_rows if r["interval_reasonable_flag"] != NA]
    reasonable_rows = [r for r in interval_evaluable if r["interval_reasonable_flag"] == "yes"]
    non_overlap_rows = [r for r in interval_evaluable if r["interval_reasonable_flag"] == "no"]
    exact_same_rho = [
        r for r in comparison_rows
        if r["matched_in_project"] == "yes" and r["abs_diff_vs_rho_kya"] != NA and float(r["abs_diff_vs_rho_kya"]) < 0.01
    ]
    exact_same_ml = [
        r for r in comparison_rows
        if r["matched_in_project"] == "yes" and r["abs_diff_vs_ml_kya"] != NA and float(r["abs_diff_vs_ml_kya"]) < 0.01
    ]

    lines = [
        "# 文献单倍群年龄与本项目年龄估计比较报告",
        "",
        "## 1. 数据来源与处理",
        "",
        f"- 文献输入目录：`download/`，共解析 {len(literature_rows)} 条文献年龄记录，覆盖 {len(literature_haps)} 个唯一单倍群。",
        f"- 本项目输入：`output/stage_2/results/ml_vs_rho_comparison.tsv`，共 {len(project_rows)} 个单倍群。",
        "- 所有年龄统一换算为 kya；`以往研究-1.tsv` 与 `以往研究-4.tsv` 原始单位为 years，已除以 1000。",
        "- 对只有区间的 ML 年龄，使用区间中点作为文献年龄，同时保留上下界；对仅给出 SD 的表，使用 age ± 1.96 SD 近似 95% CI。",
        "- 合理性主判定采用区间覆盖规则：文献年龄上下限与本项目 ρ 或 ML 的 95% CI 只要有任一交集，即认定为合理。",
        "",
        "## 2. 总体一致性",
        "",
        f"- 与本项目存在同名单倍群匹配的文献记录：{sum(1 for r in comparison_rows if r['matched_in_project'] == 'yes')} 条。",
        f"- 与本项目存在同名单倍群匹配的唯一单倍群：{len(matched_haps)} 个。",
        f"- 文献有记录但本项目未匹配的记录：{sum(1 for r in comparison_rows if r['matched_in_project'] == 'no')} 条。",
        f"- 可进行区间覆盖判定的记录：{len(interval_evaluable)} 条。",
        f"- 区间与本项目 ρ 或 ML 至少一个结果有覆盖、因此认定合理的记录：{len(reasonable_rows)} 条。",
        f"- 区间与本项目 ρ 和 ML 均不覆盖的记录：{len(non_overlap_rows)} 条。",
        f"- 与本项目 ρ 年龄几乎相同（绝对差异 < 0.01 kya）的记录：{len(exact_same_rho)} 条。",
        f"- 与本项目 ML 年龄几乎相同（绝对差异 < 0.01 kya）的记录：{len(exact_same_ml)} 条。",
        "- 绝对差异和相对差异仍保留在结果表中，但只作为参考指标，不作为合理性主判定。",
        "",
        "## 3. 按文献和方法汇总",
        "",
        "| 文献文件 | 方法 | 文献记录 | 匹配记录 | 匹配率(%) | 可判定区间 | 合理记录 | 合理率(%) | 区间不重叠 | ρ中位绝对差(kya) | ML中位绝对差(kya) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['source_file']} | {row['literature_method']} | {row['n_literature_records']} | "
            f"{row['n_matched_project']} | {row['match_rate_pct']} | {row['n_interval_evaluable']} | "
            f"{row['n_interval_reasonable']} | {row['pct_interval_reasonable']} | {row['n_interval_non_overlap']} | "
            f"{row['median_abs_diff_vs_rho_kya']} | {row['median_abs_diff_vs_ml_kya']} |"
        )

    lines.extend([
        "",
        "## 4. 最大差异记录",
        "",
        "| 单倍群 | 文献文件 | 方法 | 文献年龄(kya) | 本项目ρ(kya) | 本项目ML(kya) | 最小绝对差(kya) | 更接近 |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in top_rows[:20]:
        lines.append(
            f"| {row['haplogroup']} | {row['source_file']} | {row['literature_method']} | "
            f"{row['literature_age_kya']} | {row['project_rho_age_kya']} | {row['project_ml_age_kya']} | "
            f"{row['min_abs_diff_to_project_kya']} | {row['nearest_project_method']} |"
        )

    lines.extend([
        "",
        "## 5. 输出文件",
        "",
        "- `output/stage_3_compare/results/literature_age_normalized.tsv`：文献年龄统一格式表。",
        "- `output/stage_3_compare/results/literature_vs_project_comparison.tsv`：逐条文献记录与本项目 ρ/ML 年龄比较。",
        "- `output/stage_3_compare/results/literature_vs_project_summary.tsv`：按文献与方法汇总的匹配率和差异统计。",
        "- `output/stage_3_compare/results/literature_unmatched.tsv`：文献中存在但本项目结果中未匹配到的记录。",
        "- `output/stage_3_compare/results/top_literature_project_differences.tsv`：最大差异记录。",
    ])

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-dir", required=True)
    parser.add_argument("--project-comparison", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--large-diff-kya", type=float, default=5.0)
    parser.add_argument("--large-rel-diff-pct", type=float, default=30.0)
    args = parser.parse_args()

    download_dir = Path(args.download_dir)
    results_dir = Path(args.results_dir)
    project_rows = load_project(Path(args.project_comparison))
    literature_rows = load_literature(download_dir)
    comparison_rows = compare_rows(literature_rows, project_rows, args.large_diff_kya, args.large_rel_diff_pct)
    summary_rows = summarize(comparison_rows)
    top_rows = top_differences(comparison_rows)
    unmatched_rows = [row for row in comparison_rows if row["matched_in_project"] == "no"]

    lit_fields = list(literature_rows[0].keys()) if literature_rows else []
    cmp_fields = list(comparison_rows[0].keys()) if comparison_rows else []
    sum_fields = list(summary_rows[0].keys()) if summary_rows else []

    write_tsv(results_dir / "literature_age_normalized.tsv", literature_rows, lit_fields)
    write_tsv(results_dir / "literature_vs_project_comparison.tsv", comparison_rows, cmp_fields)
    write_tsv(results_dir / "literature_vs_project_summary.tsv", summary_rows, sum_fields)
    write_tsv(results_dir / "literature_unmatched.tsv", unmatched_rows, cmp_fields)
    write_tsv(results_dir / "top_literature_project_differences.tsv", top_rows, cmp_fields)
    write_report(Path(args.report), literature_rows, comparison_rows, summary_rows, top_rows, project_rows)


if __name__ == "__main__":
    main()

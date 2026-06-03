#!/usr/bin/env python3
"""
1-2-vcf_fps_select.py — 第二层：VCF 成对差异 + FPS 降采样（双模块）
====================================================================
输入：第一层硬过滤后的样本列表（含 major_haplogroup 列）。
流程：
  1. 用 cyvcf2 一次读取全部过滤样本的基因型
  2. 按主单倍群分块并行计算成对差异矩阵
  3. 严格折叠基因型向量完全一致的单倍型
  4. 自适应 FPS 降采样
  5. 用块级日志支持断点续跑
"""

import argparse
import fcntl
import hashlib
import json
import logging
import os
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from cyvcf2 import VCF
from joblib import Parallel, delayed
from numba import njit, prange, set_num_threads


LOG = logging.getLogger(__name__)
CACHE_COLUMNS = ["ID", "Haplogroup_17.2", "QC_Haplogrep", "Reason"]


# ─────────────────────────────────────────────────────────────────────────────
# VCF 基因型读取
# ─────────────────────────────────────────────────────────────────────────────

def load_gt_matrix(vcf_path: str, sample_ids: list[str]) -> tuple[np.ndarray, list[str]]:
    """
    一次读取指定样本的 mtDNA 基因型。

    返回：
      G           : ndarray, shape (n_sites, n_samples), dtype int8
                    0=参考态, k=ALT 等位序号, -1=缺失
      vcf_samples : 与 G 列一一对应的样本 ID

    mtDNA 共识序列应使用同型二倍体形式记录。杂合位点无法可靠转换为
    单倍型，因此直接报错。
    """
    vcf = VCF(vcf_path)
    available = set(vcf.samples)
    vcf_samples = [sample_id for sample_id in sample_ids if sample_id in available]
    if not vcf_samples:
        return np.empty((0, 0), dtype=np.int8), []

    vcf.set_samples(vcf_samples)
    rows = []
    for variant in vcf:
        alleles = variant.genotype.array()[:, :2]
        missing = np.any(alleles < 0, axis=1)
        heterozygous = (~missing) & (alleles[:, 0] != alleles[:, 1])
        if np.any(heterozygous):
            sample = vcf_samples[int(np.flatnonzero(heterozygous)[0])]
            raise ValueError(
                f"检测到 mtDNA 杂合基因型: {variant.CHROM}:{variant.POS}, sample={sample}"
            )
        first_allele = alleles[:, 0]
        if np.any(first_allele > np.iinfo(np.int8).max):
            raise ValueError(f"ALT 等位序号超出 int8 范围: {variant.CHROM}:{variant.POS}")
        rows.append(first_allele.astype(np.int8, copy=True))

    if not rows:
        return np.empty((0, len(vcf_samples)), dtype=np.int8), vcf_samples
    return np.stack(rows), vcf_samples


def get_segregating_mask(G: np.ndarray) -> np.ndarray:
    """返回块内分离位点掩码；缺失值 -1 不参与最小值和最大值计算。"""
    if G.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    valid = G >= 0
    row_max = np.where(valid, G, -1).max(axis=1)
    row_min = np.where(valid, G, np.iinfo(np.int8).max).min(axis=1)
    return valid.any(axis=1) & (row_max > row_min)


# ─────────────────────────────────────────────────────────────────────────────
# 成对差异与严格折叠
# ─────────────────────────────────────────────────────────────────────────────

@njit(cache=False, nogil=True, parallel=True)
def _pairwise_diff_numba(Gt: np.ndarray) -> np.ndarray:
    """Numba 内核：任一侧缺失时不计该位点差异。"""
    n_samples, n_sites = Gt.shape
    D = np.zeros((n_samples, n_samples), dtype=np.int16)
    for i in prange(n_samples):
        for j in range(i + 1, n_samples):
            diff = 0
            for site in range(n_sites):
                left = Gt[i, site]
                right = Gt[j, site]
                if left >= 0 and right >= 0 and left != right:
                    diff += 1
            D[i, j] = diff
            D[j, i] = diff
    return D


def pairwise_diff(Gt: np.ndarray) -> np.ndarray:
    """
    计算成对差异矩阵。

    参数 Gt shape: (n_samples, n_var_sites), dtype int8, -1=缺失。
    返回 shape: (n_samples, n_samples), dtype int16。
    """
    return _pairwise_diff_numba(np.ascontiguousarray(Gt, dtype=np.int8))


def collapse_identical(Gt: np.ndarray) -> tuple[np.ndarray, dict[int, int]]:
    """
    严格折叠基因型向量完全一致的样本。

    缺失模式也必须一致，避免把“忽略缺失后距离为 0”误当成等价关系。
    返回代表索引，以及代表索引到簇大小的映射。
    """
    if Gt.shape[0] == 0:
        return np.empty(0, dtype=int), {}
    _, first_indices, counts = np.unique(
        Gt, axis=0, return_index=True, return_counts=True
    )
    order = np.argsort(first_indices)
    reps = first_indices[order]
    cluster_sizes = {
        int(rep): int(count) for rep, count in zip(reps, counts[order])
    }
    return reps, cluster_sizes


# ─────────────────────────────────────────────────────────────────────────────
# 自适应 FPS
# ─────────────────────────────────────────────────────────────────────────────

def fps_adaptive(D_reps: np.ndarray, min_dist: int, max_tips: int) -> tuple[list, list]:
    """执行自适应最远点采样，返回代表索引和选中时的最近距离。"""
    n_reps = D_reps.shape[0]
    if n_reps == 0:
        return [], []
    if n_reps == 1:
        return [0], [-1]

    selected = [0]
    distances = [-1]
    nearest = D_reps[0].astype(np.int32).copy()
    nearest[0] = 0

    while True:
        best_dist = int(nearest.max())
        if best_dist <= min_dist:
            break
        if max_tips > 0 and len(selected) >= max_tips:
            break
        next_index = int(np.argmax(nearest))
        selected.append(next_index)
        distances.append(best_dist)
        nearest = np.minimum(nearest, D_reps[next_index].astype(np.int32))
        nearest[next_index] = 0
    return selected, distances


# ─────────────────────────────────────────────────────────────────────────────
# 日志驱动断点续跑
# ─────────────────────────────────────────────────────────────────────────────

def _write_json_atomic(path: Path, data: dict) -> None:
    """原子写 JSON，避免中断后留下半截文件。"""
    temp_path = path.with_suffix(f"{path.suffix}.tmp.{os.getpid()}")
    temp_path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temp_path, path)


def _append_locked(log_dir: Path, filename: str, text: str) -> None:
    """在 flock 保护下追加日志。"""
    lock_path = log_dir / ".resume.lock"
    with lock_path.open("a", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        with (log_dir / filename).open("a", encoding="utf-8") as log_handle:
            log_handle.write(text)
        fcntl.flock(lock_handle, fcntl.LOCK_UN)


def prepare_resume(log_dir: Path, temp_dir: Path, fingerprint: dict,
                   overwrite: bool) -> set[str]:
    """
    初始化日志并读取成功块集合。

    完成状态只来自 success.log 第一列。运行指纹变化时必须显式覆盖，
    防止不同参数的块缓存静默混用。
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    success_log = log_dir / "success.log"
    fail_log = log_dir / "fail.log"
    fingerprint_path = log_dir / "run_fingerprint.json"

    if overwrite:
        success_log.write_text("", encoding="utf-8")
        fail_log.write_text("", encoding="utf-8")
        for cache in temp_dir.glob("block_*_selected.tsv"):
            cache.unlink()
        for cache in temp_dir.glob("block_*_summary.json"):
            cache.unlink()
        _write_json_atomic(fingerprint_path, fingerprint)
        return set()

    if fingerprint_path.exists():
        previous = json.loads(fingerprint_path.read_text(encoding="utf-8"))
        if previous != fingerprint:
            raise RuntimeError(
                "运行指纹与现有断点状态不一致；请确认参数后使用 overwrite=true 重算。"
            )
    else:
        _write_json_atomic(fingerprint_path, fingerprint)

    success_log.touch()
    fail_log.touch()
    return {
        line.split("\t", 1)[0]
        for line in success_log.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def _file_sha256(path: str) -> str:
    """计算小型输入文件哈希，用于识别第一层筛选结果变化。"""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_fingerprint(filtered_list: str, vcf_path: str,
                      min_dist: int, max_tips: int) -> dict:
    """生成影响块结果的运行指纹。"""
    vcf = Path(vcf_path).resolve()
    stat = vcf.stat()
    return {
        "filtered_list_sha256": _file_sha256(filtered_list),
        "vcf_path": str(vcf),
        "vcf_size": stat.st_size,
        "vcf_mtime_ns": stat.st_mtime_ns,
        "min_dist": min_dist,
        "max_tips": max_tips,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 单块处理
# ─────────────────────────────────────────────────────────────────────────────

def _write_tsv_atomic(data: pd.DataFrame, path: Path) -> None:
    """原子写块 TSV。"""
    temp_path = path.with_suffix(f"{path.suffix}.tmp.{os.getpid()}")
    data.to_csv(temp_path, sep="\t", index=False)
    os.replace(temp_path, path)


def process_block(block_name: str, block_df: pd.DataFrame,
                  G: np.ndarray, sample_indices: np.ndarray,
                  min_dist: int, max_tips: int,
                  temp_dir: str, numba_threads: int) -> tuple[str, pd.DataFrame, dict]:
    """处理一个主单倍群块，并原子写入块缓存。"""
    cache = Path(temp_dir) / f"block_{block_name}_selected.tsv"
    summary_cache = Path(temp_dir) / f"block_{block_name}_summary.json"
    n_input = len(block_df)
    set_num_threads(numba_threads)
    LOG.info("[块 %s] 输入样本: %s", block_name, f"{n_input:,}")

    block_G = G[:, sample_indices]
    seg_mask = get_segregating_mask(block_G)
    n_seg = int(seg_mask.sum())
    LOG.info("[块 %s] 块内分离位点: %s", block_name, f"{n_seg:,}")

    Gt = np.ascontiguousarray(block_G[seg_mask].T)
    reps, cluster_sizes = collapse_identical(np.ascontiguousarray(block_G.T))
    n_unique = len(reps)

    if n_seg == 0:
        selected_global = np.array([reps[0]])
        selected_local = [0]
        distances = [-1]
    else:
        LOG.info("[块 %s] 计算成对差异矩阵 (%s x %s)",
                 block_name, f"{n_input:,}", f"{n_input:,}")
        D = pairwise_diff(Gt)
        Dr = D[np.ix_(reps, reps)]
        selected_local, distances = fps_adaptive(Dr, min_dist, max_tips)
        selected_global = reps[selected_local]

    rows = block_df.iloc[selected_global].copy()
    reasons = []
    for rank, (local_index, distance) in enumerate(
            zip(selected_local, distances), start=1):
        rep = int(reps[local_index])
        cluster_size = cluster_sizes[rep]
        if rank == 1:
            reasons.append(f"fps_seed;block={block_name};cluster={cluster_size}")
        else:
            reasons.append(
                f"fps_diverse;block={block_name};rank={rank};"
                f"dist={distance};cluster={cluster_size}"
            )
    rows["Reason"] = reasons
    rows = rows[CACHE_COLUMNS]
    _write_tsv_atomic(rows, cache)

    summary = {
        "major_haplogroup": block_name,
        "n_layer1_input": n_input,
        "n_unique_haplotypes": n_unique,
        "n_fps_selected": len(rows),
        "n_segregating_sites": n_seg,
    }
    _write_json_atomic(summary_cache, summary)
    LOG.info("[块 %s] 唯一单倍型: %s; FPS 选中: %s",
             block_name, f"{n_unique:,}", f"{len(rows):,}")
    return block_name, rows, summary


def _process_block_safe(*args, **kwargs) -> tuple:
    """捕获单块异常，让其他块继续运行。"""
    block_name = args[0]
    try:
        return True, process_block(*args, **kwargs), ""
    except Exception:
        return False, block_name, traceback.format_exc()


# ─────────────────────────────────────────────────────────────────────────────
# 主流程与 CLI
# ─────────────────────────────────────────────────────────────────────────────

def run(filtered_list: str, vcf_path: str, output_path: str,
        summary_path: str, temp_dir: str, log_dir: str,
        min_dist: int, max_tips: int, jobs: int, overwrite: bool) -> int:
    """执行第二层流程，返回 0 表示成功。"""
    LOG.info("[1-2] 读取硬过滤结果: %s", filtered_list)
    df = pd.read_csv(filtered_list, sep="\t", low_memory=False)
    if "major_haplogroup" not in df.columns:
        df["major_haplogroup"] = (
            df["Haplogroup_17.2"].astype(str).str[0].str.upper().fillna("Unknown")
        )

    temp_path = Path(temp_dir)
    log_path = Path(log_dir)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
    fingerprint = build_fingerprint(filtered_list, vcf_path, min_dist, max_tips)
    completed = prepare_resume(log_path, temp_path, fingerprint, overwrite)

    sample_ids = df["ID"].astype(str).tolist()
    LOG.info("[1-2] 用 cyvcf2 一次读取 %s 个过滤样本", f"{len(sample_ids):,}")
    G, vcf_samples = load_gt_matrix(vcf_path, sample_ids)
    if len(vcf_samples) != len(sample_ids):
        missing = sorted(set(sample_ids) - set(vcf_samples))
        raise RuntimeError(f"{len(missing)} 个过滤样本不在 VCF 中，例如: {missing[:3]}")
    index_by_sample = {sample_id: index for index, sample_id in enumerate(vcf_samples)}

    blocks = sorted(df["major_haplogroup"].unique())
    results = []
    summaries = []
    pending = []
    numba_threads = max(1, (os.cpu_count() or jobs) // jobs)
    for block_name in blocks:
        cache = temp_path / f"block_{block_name}_selected.tsv"
        summary_cache = temp_path / f"block_{block_name}_summary.json"
        block_df = df[df["major_haplogroup"] == block_name].copy()
        if block_name in completed:
            if not cache.exists():
                raise RuntimeError(f"成功日志存在但块缓存缺失: {cache}")
            if not summary_cache.exists():
                raise RuntimeError(f"成功日志存在但块摘要缺失: {summary_cache}")
            LOG.info("[块 %s] 使用成功日志确认的缓存", block_name)
            cached = pd.read_csv(cache, sep="\t")
            results.append(cached)
            summaries.append(json.loads(summary_cache.read_text(encoding="utf-8")))
            continue
        indices = np.array(
            [index_by_sample[str(sample_id)] for sample_id in block_df["ID"]],
            dtype=np.int64,
        )
        pending.append(
            (block_name, block_df, G, indices, min_dist, max_tips, temp_dir, numba_threads)
        )

    LOG.info("[1-2] 待计算块: %s; 并行 jobs=%s; 每块 numba_threads=%s",
             len(pending), jobs, numba_threads)
    outcomes = Parallel(n_jobs=jobs, backend="threading")(
        delayed(_process_block_safe)(*args) for args in pending
    )
    failed = []
    for success, payload, details in outcomes:
        if success:
            block_name, rows, summary = payload
            results.append(rows)
            summaries.append(summary)
            _append_locked(log_path, "success.log",
                           f"{block_name}\t{datetime.now().astimezone().isoformat()}\n")
        else:
            block_name = payload
            failed.append(block_name)
            _append_locked(log_path, "fail.log", f"{block_name}\n")
            LOG.error("[块 %s] 失败:\n%s", block_name, details)

    if failed:
        LOG.error("[1-2] %s 个块失败: %s", len(failed), ", ".join(failed))
        return 1
    if not results:
        LOG.error("[1-2] 所有块均为空，无输出")
        return 1

    combined = pd.concat(results, ignore_index=True)
    combined[CACHE_COLUMNS].to_csv(output_path, sep="\t", index=False)
    summary_df = pd.DataFrame(summaries).sort_values("major_haplogroup")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    LOG.info("[1-2] 最终选中样本: %s -> %s", f"{len(combined):,}", output_path)
    LOG.info("[1-2] 汇总报告: %s", summary_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """构建 CLI 参数解析器。"""
    parser = argparse.ArgumentParser(
        description="Step 1-2: VCF 成对差异 + FPS 降采样",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--filtered-list", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--temp", required=True)
    parser.add_argument("--log-dir", required=True)
    parser.add_argument("--min-dist", type=int, default=1)
    parser.add_argument("--max-tips", type=int, default=0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--overwrite", default="false")
    return parser


def main(argv=None) -> int:
    """CLI 入口。"""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    if args.jobs < 1:
        raise ValueError("--jobs 必须大于等于 1")
    overwrite = args.overwrite.lower() in ("true", "1", "yes")
    return run(
        filtered_list=args.filtered_list,
        vcf_path=args.vcf,
        output_path=args.output,
        summary_path=args.summary,
        temp_dir=args.temp,
        log_dir=args.log_dir,
        min_dist=args.min_dist,
        max_tips=args.max_tips,
        jobs=args.jobs,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""按作业清单并行运行 BEAST + TreeAnnotator。

对应原 PAML 版本的 ``run_paml_jobs.py``：复用基于 success.log / fail.log 的断点续跑
追踪器与 CPU 亲和并行框架。每个作业先运行 ``beast`` 产出后验 ``.trees``，再运行
``treeannotator`` 汇总为 MCC 树。

断点续跑遵循规范：完成判定**只**依据 success.log，绝不扫描输出目录或输出文件判定。
"""

from __future__ import annotations

import argparse
import concurrent.futures
import queue
import sys
import time
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from beast_common import (
    BEAST_RUN_STATUS_COLUMNS,
    ResumeTracker,
    load_rows,
    plan_cpu_slots,
    run_beast_job,
    run_treeannotator,
    write_status_rows,
)


def _mcc_looks_complete(mcc_file: Path) -> bool:
    """轻量检查 MCC 文件是否存在且包含 tree 行（仅用于 legacy skip，不用于续跑判定）。"""
    if not mcc_file.exists() or mcc_file.stat().st_size == 0:
        return False
    try:
        with mcc_file.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.lstrip().lower().startswith("tree "):
                    return True
    except OSError:
        return False
    return False


def _run_single_job(
    row,
    beast_bin,
    treeannotator_bin,
    burnin_percent,
    mcc_heights,
    base_seed,
    legacy_skip_existing,
    resume_tracker,
    logger,
    job_index,
    total_jobs,
    cpu_slot_queue=None,
):
    job_dir = Path(row["job_dir"])
    xmlfile = Path(row["xmlfile"])
    trees_file = Path(row["trees_file"])
    mcc_file = Path(row["mcc_file"])
    subtree_id = row["beast_subtree_id"]
    start_time = time.monotonic()

    cpu_ids = None
    cpu_label = "system-scheduled"
    if cpu_slot_queue is not None:
        cpu_ids = cpu_slot_queue.get()
        cpu_label = ",".join(str(cpu_id) for cpu_id in cpu_ids)

    try:
        logger.info(
            "Starting %s (%d/%d) job_dir=%s cpus=%s",
            subtree_id, job_index, total_jobs, job_dir.as_posix(), cpu_label,
        )

        if legacy_skip_existing and _mcc_looks_complete(mcc_file):
            logger.info("Finished %s status=skipped_existing", subtree_id)
            return {
                "beast_subtree_id": subtree_id,
                "job_dir": job_dir.as_posix(),
                "trees_file": trees_file.as_posix(),
                "mcc_file": mcc_file.as_posix(),
                "return_code": "0",
                "status": "skipped_existing",
                "detail": "Existing MCC retained",
            }

        # 各作业用 base_seed + 子树名哈希偏移，保证可复现且互不相同
        job_seed = int(base_seed) + (abs(hash(subtree_id)) % 1_000_000)

        # 1) 运行 BEAST
        beast_completed = run_beast_job(beast_bin, xmlfile, job_dir, seed=job_seed, cpu_ids=cpu_ids)
        (job_dir / "beast.stdout.log").write_text(beast_completed.stdout or "", encoding="utf-8")
        (job_dir / "beast.stderr.log").write_text(beast_completed.stderr or "", encoding="utf-8")
        if beast_completed.returncode != 0 or not trees_file.exists():
            detail = (beast_completed.stderr or beast_completed.stdout or "beast failed").strip()
            resume_tracker.record_fail(subtree_id)
            logger.error("Finished %s status=failed (beast) detail=%s", subtree_id, detail[:200])
            return {
                "beast_subtree_id": subtree_id,
                "job_dir": job_dir.as_posix(),
                "trees_file": trees_file.as_posix(),
                "mcc_file": mcc_file.as_posix(),
                "return_code": str(beast_completed.returncode),
                "status": "failed",
                "detail": f"beast: {detail[:300]}",
            }

        # 2) 运行 TreeAnnotator
        ta_completed = run_treeannotator(
            treeannotator_bin, trees_file, mcc_file, job_dir,
            burnin_percent=burnin_percent, heights=mcc_heights,
        )
        (job_dir / "treeannotator.stdout.log").write_text(ta_completed.stdout or "", encoding="utf-8")
        (job_dir / "treeannotator.stderr.log").write_text(ta_completed.stderr or "", encoding="utf-8")
        if ta_completed.returncode != 0 or not _mcc_looks_complete(mcc_file):
            detail = (ta_completed.stderr or ta_completed.stdout or "treeannotator failed").strip()
            resume_tracker.record_fail(subtree_id)
            logger.error("Finished %s status=failed (treeannotator) detail=%s", subtree_id, detail[:200])
            return {
                "beast_subtree_id": subtree_id,
                "job_dir": job_dir.as_posix(),
                "trees_file": trees_file.as_posix(),
                "mcc_file": mcc_file.as_posix(),
                "return_code": str(ta_completed.returncode),
                "status": "failed",
                "detail": f"treeannotator: {detail[:300]}",
            }

        elapsed = time.monotonic() - start_time
        resume_tracker.record_success(subtree_id)
        logger.info("Finished %s (%d/%d) status=success elapsed=%.1fs cpus=%s",
                    subtree_id, job_index, total_jobs, elapsed, cpu_label)
        return {
            "beast_subtree_id": subtree_id,
            "job_dir": job_dir.as_posix(),
            "trees_file": trees_file.as_posix(),
            "mcc_file": mcc_file.as_posix(),
            "return_code": "0",
            "status": "success",
            "detail": "beast + treeannotator finished successfully",
        }
    except Exception as exc:  # pragma: no cover - 防御性兜底
        resume_tracker.record_fail(subtree_id)
        logger.error("Finished %s status=failed detail=%s", subtree_id, exc)
        return {
            "beast_subtree_id": subtree_id,
            "job_dir": job_dir.as_posix(),
            "trees_file": trees_file.as_posix(),
            "mcc_file": mcc_file.as_posix(),
            "return_code": "1",
            "status": "failed",
            "detail": str(exc),
        }
    finally:
        if cpu_slot_queue is not None and cpu_ids is not None:
            cpu_slot_queue.put(cpu_ids)


def run(
    manifest_path,
    beast_bin,
    treeannotator_bin,
    parallel_jobs,
    threads_per_job=1,
    burnin_percent=10.0,
    mcc_heights="mean",
    base_seed=42,
    bind_cpu_affinity=False,
    resume_log_dir=None,
    legacy_skip_existing=False,
    log_level="INFO",
):
    manifest_path = Path(manifest_path).resolve()
    rows = load_rows(manifest_path)
    if not rows:
        raise PipelineError(f"BEAST job manifest is empty: {manifest_path}")

    output_dir = manifest_path.parent
    resume_tracker = ResumeTracker(output_dir / "log" if resume_log_dir is None else Path(resume_log_dir))
    resume_tracker.initialize()
    logger = setup_logger("run_beast_jobs", output_dir / "run_beast_jobs.log", log_level)
    logger.info(
        "Running %d BEAST jobs parallel=%s threads/job=%s bind_cpu=%s resume_dir=%s",
        len(rows), parallel_jobs, threads_per_job, bind_cpu_affinity, resume_tracker.log_dir,
    )

    statuses = []
    total_jobs = len(rows)
    pending_rows = []
    for row in rows:
        subtree_id = row["beast_subtree_id"]
        if resume_tracker.is_completed(subtree_id):
            statuses.append({
                "beast_subtree_id": subtree_id,
                "job_dir": row["job_dir"],
                "trees_file": row["trees_file"],
                "mcc_file": row["mcc_file"],
                "return_code": "0",
                "status": "skipped_logged",
                "detail": "Already recorded in success.log",
            })
            logger.info("Skipping %s (recorded in success.log)", subtree_id)
            continue
        pending_rows.append(row)

    logger.info("Resume summary: total=%d skipped_logged=%d pending=%d",
                total_jobs, len(statuses), len(pending_rows))

    cpu_slot_queue = None
    max_workers = max(1, int(parallel_jobs))
    if bind_cpu_affinity and pending_rows:
        cpu_slots = plan_cpu_slots(int(threads_per_job))
        max_workers = min(max_workers, len(cpu_slots), len(pending_rows))
        cpu_slot_queue = queue.Queue()
        for cpu_slot in cpu_slots[:max_workers]:
            cpu_slot_queue.put(cpu_slot)
        logger.info("CPU affinity enabled: slots=%d effective_parallel=%d", len(cpu_slots), max_workers)
    else:
        max_workers = min(max_workers, len(pending_rows)) if pending_rows else 1

    if pending_rows:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(
                    _run_single_job,
                    row, beast_bin, treeannotator_bin, burnin_percent, mcc_heights, base_seed,
                    legacy_skip_existing, resume_tracker, logger, index, total_jobs, cpu_slot_queue,
                ): row["beast_subtree_id"]
                for index, row in enumerate(pending_rows, start=1 + len(statuses))
            }
            completed_count = len(statuses)
            for future in concurrent.futures.as_completed(future_map):
                subtree_id = future_map[future]
                try:
                    status_row = future.result()
                except Exception as exc:  # pragma: no cover
                    resume_tracker.record_fail(subtree_id)
                    status_row = {
                        "beast_subtree_id": subtree_id, "job_dir": "", "trees_file": "", "mcc_file": "",
                        "return_code": "1", "status": "failed", "detail": str(exc),
                    }
                statuses.append(status_row)
                completed_count += 1
                logger.info("Progress %d/%d: %s -> %s",
                            completed_count, total_jobs, status_row["beast_subtree_id"], status_row["status"])

    statuses.sort(key=lambda r: r["beast_subtree_id"])
    write_status_rows(statuses, BEAST_RUN_STATUS_COLUMNS, output_dir / "beast_run_status.tsv")

    failures = [r for r in statuses if r["status"] == "failed"]
    if failures:
        details = "; ".join(f"{r['beast_subtree_id']}: {r['detail']}" for r in failures[:5])
        raise PipelineError(f"{len(failures)} BEAST job(s) failed: {details}")
    logger.info("All BEAST jobs completed successfully.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Run BEAST + TreeAnnotator jobs from a manifest.")
    parser.add_argument("--manifest", required=True, help="Path to beast_job_manifest.tsv.")
    parser.add_argument("--beast-bin", required=True, help="Path/command for the beast executable.")
    parser.add_argument("--treeannotator-bin", required=True, help="Path/command for treeannotator.")
    parser.add_argument("--parallel-jobs", required=True, type=int, help="Maximum parallel jobs.")
    parser.add_argument("--threads-per-job", type=int, default=1, help="CPU cores reserved per job.")
    parser.add_argument("--burnin-percent", type=float, default=10.0, help="TreeAnnotator burn-in percentage.")
    parser.add_argument("--mcc-heights", default="mean", help="TreeAnnotator node heights: mean/median/ca/keep.")
    parser.add_argument("--base-seed", type=int, default=42, help="Base random seed.")
    parser.add_argument("--bind-cpu-affinity", action="store_true", help="Bind each job to an explicit CPU slot.")
    parser.add_argument("--resume-log-dir", default=None, help="Resume state directory. Defaults to <manifest_dir>/log.")
    parser.add_argument("--legacy-skip-existing", action="store_true", dest="legacy_skip_existing",
                        help="Legacy fallback: skip jobs when an existing MCC already looks complete.")
    parser.add_argument("--skip-existing", action="store_true", dest="legacy_skip_existing", help=argparse.SUPPRESS)
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            manifest_path=args.manifest,
            beast_bin=args.beast_bin,
            treeannotator_bin=args.treeannotator_bin,
            parallel_jobs=args.parallel_jobs,
            threads_per_job=args.threads_per_job,
            burnin_percent=args.burnin_percent,
            mcc_heights=args.mcc_heights,
            base_seed=args.base_seed,
            bind_cpu_affinity=args.bind_cpu_affinity,
            resume_log_dir=args.resume_log_dir,
            legacy_skip_existing=args.legacy_skip_existing,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run baseml jobs from a prepared manifest."""

from __future__ import annotations

import argparse
import concurrent.futures
import queue
import re
import sys
import time
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from paml_common import PAML_RUN_STATUS_COLUMNS, load_rows, plan_cpu_slots, run_baseml_job, write_rows


def _has_complete_baseml_output(outfile: Path) -> bool:
    if not outfile.exists() or outfile.stat().st_size == 0:
        return False
    text = outfile.read_text(encoding="utf-8", errors="replace")
    has_lnl = re.search(r"^\s*lnL\(", text, flags=re.MULTILINE) is not None
    return has_lnl and "SEs for parameters" in text and "tree length" in text


def _has_known_numeric_failure(*paths: Path) -> bool:
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if "fh =" in text and "negative" in text:
            return True
    return False


def _rewrite_ctl_file(source_ctl: Path, destination_ctl: Path, overrides):
    lines = source_ctl.read_text(encoding="utf-8").splitlines()
    rendered = []
    seen = set()
    for line in lines:
        updated = False
        for key, value in overrides.items():
            if re.match(rf"^\s*\*?\s*{re.escape(key)}\s*=", line):
                prefix = re.match(r"^(\s*)", line).group(1)
                rendered.append(f"{prefix}{key} = {value}")
                seen.add(key)
                updated = True
                break
        if not updated:
            rendered.append(line)
    for key, value in overrides.items():
        if key not in seen:
            rendered.append(f"{key} = {value}")
    destination_ctl.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def _run_single_job(
    row,
    baseml_bin: str,
    skip_existing: bool,
    logger,
    job_index: int,
    total_jobs: int,
    cpu_slot_queue=None,
):
    job_dir = Path(row["job_dir"])
    ctlfile = Path(row["ctlfile"])
    outfile = Path(row["outfile"])
    stdout_log = job_dir / "baseml.stdout.log"
    stderr_log = job_dir / "baseml.stderr.log"
    baseml_subtree_id = row["baseml_subtree_id"]
    start_time = time.monotonic()
    cpu_ids = None
    cpu_label = "system-scheduled"
    if cpu_slot_queue is not None:
        cpu_ids = cpu_slot_queue.get()
        cpu_label = ",".join(str(cpu_id) for cpu_id in cpu_ids)

    try:
        logger.info(
            "Starting %s (%d/%d) job_dir=%s cpus=%s",
            baseml_subtree_id,
            job_index,
            total_jobs,
            job_dir.as_posix(),
            cpu_label,
        )

        if skip_existing and _has_complete_baseml_output(outfile):
            elapsed_seconds = time.monotonic() - start_time
            logger.info(
                "Finished %s (%d/%d) status=skipped_existing elapsed=%.1fs",
                baseml_subtree_id,
                job_index,
                total_jobs,
                elapsed_seconds,
            )
            return {
                "baseml_subtree_id": baseml_subtree_id,
                "job_dir": job_dir.as_posix(),
                "outfile": outfile.as_posix(),
                "return_code": "0",
                "status": "skipped_existing",
                "detail": "Existing mlb retained",
            }

        attempts = [
            ("primary", ctlfile, {}),
            ("retry_fix_blength2", job_dir / "baseml.retry_fix_blength2.ctl", {"fix_blength": "2"}),
        ]
        if _has_known_numeric_failure(outfile, stdout_log, stderr_log):
            attempts = [attempts[1]]
        log_sections = []
        completed = None
        detail = ""
        status = "failed"

        for attempt_name, attempt_ctl, overrides in attempts:
            logger.info("%s attempt=%s ctl=%s cpus=%s", baseml_subtree_id, attempt_name, attempt_ctl.name, cpu_label)
            if overrides:
                _rewrite_ctl_file(ctlfile, attempt_ctl, overrides)
            if outfile.exists():
                outfile.unlink()
            completed = run_baseml_job(baseml_bin, attempt_ctl, job_dir, cpu_ids=cpu_ids)
            log_sections.append(
                f"===== {attempt_name} stdout =====\n{completed.stdout}\n"
            )
            log_sections.append(
                f"===== {attempt_name} stderr =====\n{completed.stderr}\n"
            )
            if completed.returncode == 0 and _has_complete_baseml_output(outfile):
                status = "success"
                if attempt_name == "primary":
                    detail = "baseml finished successfully"
                else:
                    detail = f"Recovered with fallback profile: {attempt_name}"
                break
            detail = completed.stderr.strip() or completed.stdout.strip() or f"baseml failed during {attempt_name}"

        stdout_log.write_text("".join(log_sections[0::2]), encoding="utf-8")
        stderr_log.write_text("".join(log_sections[1::2]), encoding="utf-8")

        if completed is None:
            raise PipelineError(f"Internal error: no baseml attempt executed for {baseml_subtree_id}")

        elapsed_seconds = time.monotonic() - start_time
        logger.info(
            "Finished %s (%d/%d) status=%s return_code=%s elapsed=%.1fs cpus=%s",
            baseml_subtree_id,
            job_index,
            total_jobs,
            status,
            completed.returncode,
            elapsed_seconds,
            cpu_label,
        )

        return {
            "baseml_subtree_id": baseml_subtree_id,
            "job_dir": job_dir.as_posix(),
            "outfile": outfile.as_posix(),
            "return_code": str(completed.returncode),
            "status": status,
            "detail": detail,
        }
    finally:
        if cpu_slot_queue is not None and cpu_ids is not None:
            cpu_slot_queue.put(cpu_ids)


def run(
    manifest_path,
    baseml_bin,
    parallel_jobs,
    threads_per_job=1,
    bind_cpu_affinity=False,
    skip_existing=True,
    log_level="INFO",
):
    manifest_path = Path(manifest_path).resolve()
    rows = load_rows(manifest_path)
    if not rows:
        raise PipelineError(f"PAML job manifest is empty: {manifest_path}")

    output_dir = manifest_path.parent
    logger = setup_logger("run_paml_jobs", output_dir / "run_paml_jobs.log", log_level)
    logger.info(
        "Running %d baseml jobs with requested_parallel_jobs=%s threads_per_job=%s bind_cpu_affinity=%s",
        len(rows),
        parallel_jobs,
        threads_per_job,
        bind_cpu_affinity,
    )

    statuses = []
    total_jobs = len(rows)
    cpu_slot_queue = None
    max_workers = max(1, int(parallel_jobs))
    if bind_cpu_affinity:
        cpu_slots = plan_cpu_slots(int(threads_per_job))
        max_workers = min(max_workers, len(cpu_slots), total_jobs)
        cpu_slot_queue = queue.Queue()
        for cpu_slot in cpu_slots[:max_workers]:
            cpu_slot_queue.put(cpu_slot)
        logger.info(
            "CPU affinity enabled: available_slots=%d effective_parallel_jobs=%d slots=%s",
            len(cpu_slots),
            max_workers,
            ";".join(",".join(str(cpu_id) for cpu_id in slot) for slot in cpu_slots[:max_workers]),
        )
    else:
        max_workers = min(max_workers, total_jobs)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(
                _run_single_job,
                row,
                baseml_bin,
                skip_existing,
                logger,
                index,
                total_jobs,
                cpu_slot_queue,
            ): row["baseml_subtree_id"]
            for index, row in enumerate(rows, start=1)
        }
        completed_count = 0
        for future in concurrent.futures.as_completed(future_map):
            status_row = future.result()
            statuses.append(status_row)
            completed_count += 1
            logger.info(
                "Progress %d/%d: %s -> %s",
                completed_count,
                total_jobs,
                status_row["baseml_subtree_id"],
                status_row["status"],
            )

    statuses.sort(key=lambda row: row["baseml_subtree_id"])
    status_path = output_dir / "paml_run_status.tsv"
    write_rows(statuses, PAML_RUN_STATUS_COLUMNS, status_path)

    failures = [row for row in statuses if row["status"] == "failed"]
    if failures:
        details = "; ".join(f"{row['baseml_subtree_id']}: {row['detail']}" for row in failures[:5])
        raise PipelineError(f"{len(failures)} baseml job(s) failed: {details}")
    logger.info("All baseml jobs completed successfully.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Run baseml jobs from a manifest.")
    parser.add_argument("--manifest", required=True, help="Path to paml_job_manifest.tsv.")
    parser.add_argument("--baseml-bin", required=True, help="Path to the baseml executable.")
    parser.add_argument("--parallel-jobs", required=True, type=int, help="Maximum parallel baseml jobs.")
    parser.add_argument("--threads-per-job", type=int, default=1, help="OpenMP threads reserved for each baseml job.")
    parser.add_argument("--bind-cpu-affinity", action="store_true", help="Bind each baseml job to an explicit CPU slot.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip jobs with an existing non-empty mlb file.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            manifest_path=args.manifest,
            baseml_bin=args.baseml_bin,
            parallel_jobs=args.parallel_jobs,
            threads_per_job=args.threads_per_job,
            bind_cpu_affinity=args.bind_cpu_affinity,
            skip_existing=args.skip_existing,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

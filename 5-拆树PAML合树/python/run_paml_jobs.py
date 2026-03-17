#!/usr/bin/env python3
"""Run baseml jobs from a prepared manifest."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import fcntl
import queue
import re
import sys
import threading
import time
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from paml_common import (
    PAML_RUN_STATUS_COLUMNS,
    baseml_output_looks_complete,
    load_rows,
    plan_cpu_slots,
    run_baseml_job,
    write_rows,
)


class ResumeTracker:
    def __init__(self, resume_log_dir: Path):
        self.log_dir = Path(resume_log_dir).resolve()
        self.lock_file = self.log_dir / ".resume.lock"
        self.success_log = self.log_dir / "success.log"
        self.fail_log = self.log_dir / "fail.log"
        self._state_lock = threading.Lock()
        self.completed_jobs = set()

    def initialize(self) -> None:
        self.log_dir.mkdir(parents=True, exist_ok=True)
        for path in (self.lock_file, self.success_log, self.fail_log):
            path.touch(exist_ok=True)

        completed_jobs = set()
        with self.success_log.open("r", encoding="utf-8") as handle:
            for line in handle:
                job_id = line.split("\t", 1)[0].strip()
                if job_id:
                    completed_jobs.add(job_id)

        with self._state_lock:
            self.completed_jobs = completed_jobs

    def is_completed(self, job_id: str) -> bool:
        with self._state_lock:
            return job_id in self.completed_jobs

    def record_success(self, job_id: str) -> None:
        timestamp = dt.datetime.now().astimezone().isoformat(timespec="seconds")
        with self._state_lock:
            self._append_line(self.success_log, f"{job_id}\t{timestamp}\n")
            self.completed_jobs.add(job_id)

    def record_fail(self, job_id: str) -> None:
        with self._state_lock:
            self._append_line(self.fail_log, f"{job_id}\n")

    def _append_line(self, destination: Path, line: str) -> None:
        with self.lock_file.open("a", encoding="utf-8") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            try:
                with destination.open("a", encoding="utf-8") as destination_handle:
                    destination_handle.write(line)
                    destination_handle.flush()
            finally:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _has_complete_baseml_output(outfile: Path) -> bool:
    return baseml_output_looks_complete(outfile)


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
    legacy_skip_existing: bool,
    resume_tracker: ResumeTracker,
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

        if legacy_skip_existing and _has_complete_baseml_output(outfile):
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

        if status == "success":
            resume_tracker.record_success(baseml_subtree_id)
        else:
            resume_tracker.record_fail(baseml_subtree_id)

        return {
            "baseml_subtree_id": baseml_subtree_id,
            "job_dir": job_dir.as_posix(),
            "outfile": outfile.as_posix(),
            "return_code": str(completed.returncode),
            "status": status,
            "detail": detail,
        }
    except Exception as exc:
        detail = str(exc)
        resume_tracker.record_fail(baseml_subtree_id)
        logger.error("Finished %s (%d/%d) status=failed detail=%s", baseml_subtree_id, job_index, total_jobs, detail)
        return {
            "baseml_subtree_id": baseml_subtree_id,
            "job_dir": job_dir.as_posix(),
            "outfile": outfile.as_posix(),
            "return_code": "1",
            "status": "failed",
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
    resume_log_dir=None,
    legacy_skip_existing=False,
    log_level="INFO",
):
    manifest_path = Path(manifest_path).resolve()
    rows = load_rows(manifest_path)
    if not rows:
        raise PipelineError(f"PAML job manifest is empty: {manifest_path}")

    output_dir = manifest_path.parent
    resume_tracker = ResumeTracker(output_dir / "log" if resume_log_dir is None else Path(resume_log_dir))
    resume_tracker.initialize()
    logger = setup_logger("run_paml_jobs", output_dir / "run_paml_jobs.log", log_level)
    logger.info(
        "Running %d baseml jobs with requested_parallel_jobs=%s threads_per_job=%s bind_cpu_affinity=%s resume_log_dir=%s legacy_skip_existing=%s",
        len(rows),
        parallel_jobs,
        threads_per_job,
        bind_cpu_affinity,
        resume_tracker.log_dir,
        legacy_skip_existing,
    )

    statuses = []
    total_jobs = len(rows)
    pending_rows = []
    for row in rows:
        baseml_subtree_id = row["baseml_subtree_id"]
        if resume_tracker.is_completed(baseml_subtree_id):
            statuses.append(
                {
                    "baseml_subtree_id": baseml_subtree_id,
                    "job_dir": row["job_dir"],
                    "outfile": row["outfile"],
                    "return_code": "0",
                    "status": "skipped_logged",
                    "detail": "Already recorded in success.log",
                }
            )
            logger.info("Skipping %s because it is already recorded in %s", baseml_subtree_id, resume_tracker.success_log)
            continue
        pending_rows.append(row)

    logger.info(
        "Resume summary: total_jobs=%d skipped_logged=%d pending=%d",
        total_jobs,
        len(statuses),
        len(pending_rows),
    )

    cpu_slot_queue = None
    max_workers = max(1, int(parallel_jobs))
    if bind_cpu_affinity and pending_rows:
        cpu_slots = plan_cpu_slots(int(threads_per_job))
        max_workers = min(max_workers, len(cpu_slots), len(pending_rows))
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
        max_workers = min(max_workers, len(pending_rows)) if pending_rows else 1

    if pending_rows:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(
                    _run_single_job,
                    row,
                    baseml_bin,
                    legacy_skip_existing,
                    resume_tracker,
                    logger,
                    index,
                    total_jobs,
                    cpu_slot_queue,
                ): row["baseml_subtree_id"]
                for index, row in enumerate(pending_rows, start=1 + len(statuses))
            }
            completed_count = len(statuses)
            for future in concurrent.futures.as_completed(future_map):
                baseml_subtree_id = future_map[future]
                try:
                    status_row = future.result()
                except Exception as exc:  # pragma: no cover - defensive fallback
                    resume_tracker.record_fail(baseml_subtree_id)
                    status_row = {
                        "baseml_subtree_id": baseml_subtree_id,
                        "job_dir": "",
                        "outfile": "",
                        "return_code": "1",
                        "status": "failed",
                        "detail": str(exc),
                    }
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
    parser.add_argument("--resume-log-dir", default=None, help="Resume state directory. Defaults to <manifest_dir>/log.")
    parser.add_argument(
        "--legacy-skip-existing",
        action="store_true",
        dest="legacy_skip_existing",
        help="Legacy fallback: skip jobs when an existing mlb already looks complete.",
    )
    parser.add_argument("--skip-existing", action="store_true", dest="legacy_skip_existing", help=argparse.SUPPRESS)
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
            resume_log_dir=args.resume_log_dir,
            legacy_skip_existing=args.legacy_skip_existing,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

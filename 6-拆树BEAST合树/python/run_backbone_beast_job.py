#!/usr/bin/env python3
"""运行单个 backbone BEAST 作业（beast + treeannotator）。

对应原 PAML 版本的 ``run_backbone_paml_job.py``。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from beast_common import run_beast_job, run_treeannotator


def _mcc_looks_complete(mcc_file: Path) -> bool:
    if not mcc_file.exists() or mcc_file.stat().st_size == 0:
        return False
    with mcc_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.lstrip().lower().startswith("tree "):
                return True
    return False


def run(
    job_dir,
    beast_bin,
    treeannotator_bin,
    burnin_percent=10.0,
    mcc_heights="mean",
    seed=42,
    skip_existing=False,
    bind_cpu_affinity=False,
    threads=16,
    log_level="INFO",
):
    job_dir = Path(job_dir).resolve()
    logger = setup_logger("run_backbone_beast_job", job_dir / "run_backbone_beast_job.log", log_level)

    xmlfile = job_dir / "beast.xml"
    trees_file = job_dir / "backbone_only.trees"
    mcc_file = job_dir / "backbone_only.mcc.tree"
    if not xmlfile.exists():
        raise PipelineError(f"Backbone BEAST XML not found: {xmlfile}")

    if skip_existing and _mcc_looks_complete(mcc_file):
        logger.info("Backbone MCC already exists; skip running.")
        return 0

    cpu_ids = None
    if bind_cpu_affinity:
        import os
        # backbone 在各子树之后单独运行，可独占多核：绑定前 threads 个核，
        # run_beast_job 会据此自动加 -beagle_instances 实现多线程似然计算。
        cpu_ids = sorted(os.sched_getaffinity(0))[:max(1, int(threads))]

    logger.info("Running backbone BEAST job in %s", job_dir)
    beast_completed = run_beast_job(beast_bin, xmlfile, job_dir, seed=int(seed), cpu_ids=cpu_ids)
    (job_dir / "beast.stdout.log").write_text(beast_completed.stdout or "", encoding="utf-8")
    (job_dir / "beast.stderr.log").write_text(beast_completed.stderr or "", encoding="utf-8")
    if beast_completed.returncode != 0 or not trees_file.exists():
        detail = (beast_completed.stderr or beast_completed.stdout or "beast failed").strip()
        raise PipelineError(f"Backbone beast failed: {detail[:300]}")

    ta_completed = run_treeannotator(
        treeannotator_bin, trees_file, mcc_file, job_dir,
        burnin_percent=burnin_percent, heights=mcc_heights,
    )
    (job_dir / "treeannotator.stdout.log").write_text(ta_completed.stdout or "", encoding="utf-8")
    (job_dir / "treeannotator.stderr.log").write_text(ta_completed.stderr or "", encoding="utf-8")
    if ta_completed.returncode != 0 or not _mcc_looks_complete(mcc_file):
        detail = (ta_completed.stderr or ta_completed.stdout or "treeannotator failed").strip()
        raise PipelineError(f"Backbone treeannotator failed: {detail[:300]}")

    logger.info("Backbone BEAST job finished successfully.")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description="Run a single backbone BEAST job.")
    parser.add_argument("--job-dir", required=True, help="Backbone BEAST job directory.")
    parser.add_argument("--beast-bin", required=True, help="Path/command for beast.")
    parser.add_argument("--treeannotator-bin", required=True, help="Path/command for treeannotator.")
    parser.add_argument("--burnin-percent", type=float, default=10.0, help="TreeAnnotator burn-in percentage.")
    parser.add_argument("--mcc-heights", default="mean", help="TreeAnnotator node heights.")
    parser.add_argument("--seed", type=int, default=42, help="BEAST random seed.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip when MCC already looks complete.")
    parser.add_argument("--bind-cpu-affinity", action="store_true", help="Bind the job to the first N available cores.")
    parser.add_argument("--threads", type=int, default=16, help="Cores for the backbone job (=> -beagle_instances).")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            job_dir=args.job_dir,
            beast_bin=args.beast_bin,
            treeannotator_bin=args.treeannotator_bin,
            burnin_percent=args.burnin_percent,
            mcc_heights=args.mcc_heights,
            seed=args.seed,
            skip_existing=args.skip_existing,
            bind_cpu_affinity=args.bind_cpu_affinity,
            threads=args.threads,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

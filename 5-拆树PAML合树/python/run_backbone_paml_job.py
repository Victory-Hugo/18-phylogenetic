#!/usr/bin/env python3
"""Run a single backbone-only baseml job."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from phylo_split_common import PipelineError, setup_logger
from paml_common import run_baseml_job


def _has_complete_baseml_output(outfile: Path) -> bool:
    if not outfile.exists() or outfile.stat().st_size == 0:
        return False
    text = outfile.read_text(encoding="utf-8", errors="replace")
    has_lnl = re.search(r"^\s*lnL\(", text, flags=re.MULTILINE) is not None
    return has_lnl and "SEs for parameters" in text and "tree length" in text


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


def run(job_dir, baseml_bin, skip_existing=True, log_level="INFO"):
    job_dir = Path(job_dir).resolve()
    ctlfile = job_dir / "baseml.ctl"
    outfile = job_dir / "mlb"
    logger = setup_logger("run_backbone_paml_job", job_dir / "backbone_paml.log", log_level)

    if not ctlfile.exists():
        raise PipelineError(f"Backbone baseml ctl file not found: {ctlfile}")

    if skip_existing and _has_complete_baseml_output(outfile):
        logger.info("Skipping backbone-only baseml job because a complete mlb already exists.")
        return 0

    attempts = [
        ("primary", ctlfile, {}),
        ("retry_fix_blength2", job_dir / "baseml.retry_fix_blength2.ctl", {"fix_blength": "2"}),
    ]
    completed = None
    for attempt_name, attempt_ctl, overrides in attempts:
        logger.info("Running backbone-only baseml attempt=%s", attempt_name)
        if overrides:
            _rewrite_ctl_file(ctlfile, attempt_ctl, overrides)
        if outfile.exists():
            outfile.unlink()
        completed = run_baseml_job(baseml_bin, attempt_ctl, job_dir)
        (job_dir / f"{attempt_name}.stdout.log").write_text(completed.stdout, encoding="utf-8")
        (job_dir / f"{attempt_name}.stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode == 0 and _has_complete_baseml_output(outfile):
            logger.info("Backbone-only baseml finished successfully with attempt=%s", attempt_name)
            return 0

    detail = "" if completed is None else completed.stderr.strip() or completed.stdout.strip()
    raise PipelineError(f"Backbone-only baseml failed: {detail}")


def build_parser():
    parser = argparse.ArgumentParser(description="Run a single backbone-only baseml job.")
    parser.add_argument("--job-dir", required=True, help="Backbone baseml job directory.")
    parser.add_argument("--baseml-bin", required=True, help="Path to baseml.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip when mlb already looks complete.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(
            job_dir=args.job_dir,
            baseml_bin=args.baseml_bin,
            skip_existing=args.skip_existing,
            log_level=args.log_level,
        )
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

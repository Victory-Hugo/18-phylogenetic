#!/usr/bin/env python3
"""汇总 IQ2MC 流水线运行信息。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, List

LOG = logging.getLogger(__name__)


def _split_csv(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def run(
    summary_file: str,
    run_name: str,
    config_file: str,
    steps_executed: Iterable[str],
    commands_file: str,
    key_outputs: Iterable[str],
    exit_status: int,
) -> int:
    """生成 Markdown 摘要文件。"""
    out = Path(summary_file)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# IQ2MC Pipeline Summary",
        "",
        f"- Run name: `{run_name}`",
        f"- Config file: `{config_file}`",
        f"- Exit status: `{exit_status}`",
        "",
        "## Steps Executed",
    ]

    for step in steps_executed:
        lines.append(f"- {step}")

    lines.extend(["", "## Command Log", f"- `{commands_file}`", "", "## Key Outputs"])

    for item in key_outputs:
        lines.append(f"- `{item}`")

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    LOG.info("摘要文件已生成: %s", out)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Collect IQ2MC run summary")
    parser.add_argument("--summary-file", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--config-file", required=True)
    parser.add_argument("--steps-executed", required=True, help="Comma-separated list")
    parser.add_argument("--commands-file", required=True)
    parser.add_argument("--key-outputs", required=True, help="Comma-separated list")
    parser.add_argument("--exit-status", required=True, type=int)
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    return run(
        summary_file=args.summary_file,
        run_name=args.run_name,
        config_file=args.config_file,
        steps_executed=_split_csv(args.steps_executed),
        commands_file=args.commands_file,
        key_outputs=_split_csv(args.key_outputs),
        exit_status=args.exit_status,
    )


if __name__ == "__main__":
    raise SystemExit(main())

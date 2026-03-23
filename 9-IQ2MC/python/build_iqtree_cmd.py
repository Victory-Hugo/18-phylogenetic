#!/usr/bin/env python3
"""构建 IQ-TREE 运行命令。"""

from __future__ import annotations

import argparse
import json
import logging
import shlex
from pathlib import Path
from typing import Any, Dict, List

LOG = logging.getLogger(__name__)


def _as_list(value: str, sep: str = ",") -> List[str]:
    return [x.strip() for x in value.split(sep) if x.strip()]


def run(
    iqtree3_bin: str,
    msa_file: str,
    output_prefix: str,
    threads: int,
    substitution_model: str,
    rooted_tree: str = "",
    partition_file: str = "",
    partition_mode: str = "Q",
    dating_mode: str = "",
    mcmc_iter: str = "",
    mcmc_bds: str = "",
    mcmc_clock: str = "",
    redo: bool = False,
) -> Dict[str, Any]:
    """生成 IQ-TREE 命令数组与 shell 字符串。"""
    cmd: List[str] = [
        str(Path(iqtree3_bin)),
        "-s",
        str(Path(msa_file)),
        "-m",
        substitution_model,
        "-nt",
        str(threads),
        "--prefix",
        str(Path(output_prefix)),
    ]

    if partition_file:
        mode = partition_mode.strip().upper()
        if mode not in {"P", "Q"}:
            raise ValueError("partition_mode 仅支持 P 或 Q")
        cmd.extend([f"-{mode}", str(Path(partition_file))])

    if rooted_tree:
        cmd.extend(["-te", str(Path(rooted_tree))])

    if dating_mode:
        cmd.extend(["--dating", dating_mode])

    if mcmc_iter:
        parts = _as_list(mcmc_iter)
        if len(parts) != 3:
            raise ValueError("mcmc_iter 需为 burnin,samplefreq,nsample")
        cmd.extend(["--mcmc-iter", ",".join(parts)])

    if mcmc_bds:
        parts = _as_list(mcmc_bds)
        if len(parts) != 3:
            raise ValueError("mcmc_bds 需为 birth,death,sampling_fraction")
        cmd.extend(["--mcmc-bds", ",".join(parts)])

    if mcmc_clock:
        if mcmc_clock.upper() not in {"EQUAL", "IND", "CORR"}:
            raise ValueError("mcmc_clock 仅支持 EQUAL/IND/CORR")
        cmd.extend(["--mcmc-clock", mcmc_clock.upper()])
    if redo:
        cmd.append("-redo")

    result = {
        "cmd_list": cmd,
        "cmd_shell": " ".join(shlex.quote(x) for x in cmd),
    }
    LOG.info("IQ-TREE 命令构建完成")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build IQ-TREE command")
    parser.add_argument("--iqtree3-bin", required=True)
    parser.add_argument("--msa-file", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--substitution-model", required=True)
    parser.add_argument("--rooted-tree", default="")
    parser.add_argument("--partition-file", default="")
    parser.add_argument("--partition-mode", default="Q")
    parser.add_argument("--dating-mode", default="")
    parser.add_argument("--mcmc-iter", default="")
    parser.add_argument("--mcmc-bds", default="")
    parser.add_argument("--mcmc-clock", default="")
    parser.add_argument("--redo", action="store_true")
    parser.add_argument("--shell-only", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    result = run(
        iqtree3_bin=args.iqtree3_bin,
        msa_file=args.msa_file,
        output_prefix=args.output_prefix,
        threads=args.threads,
        substitution_model=args.substitution_model,
        rooted_tree=args.rooted_tree,
        partition_file=args.partition_file,
        partition_mode=args.partition_mode,
        dating_mode=args.dating_mode,
        mcmc_iter=args.mcmc_iter,
        mcmc_bds=args.mcmc_bds,
        mcmc_clock=args.mcmc_clock,
        redo=args.redo,
    )

    if args.shell_only:
        print(result["cmd_shell"])
    else:
        print(json.dumps(result, ensure_ascii=False, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

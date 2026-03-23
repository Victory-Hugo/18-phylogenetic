#!/usr/bin/env python3
"""根据配置修订 MCMCTree 控制文件。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict

LOG = logging.getLogger(__name__)


KEY_VALUE_FIELDS = {
    "seqfile",
    "treefile",
    "mcmcfile",
    "outfile",
    "ckpfile",
    "hessianfile",
    "clock",
    "BDparas",
    "burnin",
    "sampfreq",
    "nsample",
    "checkpoint",
}

CLOCK_MAP: Dict[str, str] = {
    "EQUAL": "1",
    "IND": "2",
    "CORR": "3",
}


def _split_comment(line: str) -> tuple[str, str]:
    if "*" in line:
        left, right = line.split("*", 1)
        return left.rstrip(), "*" + right
    return line.rstrip(), ""


def run(
    input_ctl: str,
    output_ctl: str,
    seqfile: str,
    treefile: str,
    mcmcfile: str,
    outfile: str,
    ckpfile: str,
    hessianfile: str,
    clock_model: str,
    mcmc_iter: str,
    mcmc_bds: str,
    checkpoint: int,
    resume: bool,
) -> int:
    """更新 mcmctree.ctl 关键参数。"""
    in_path = Path(input_ctl)
    out_path = Path(output_ctl)
    if not in_path.is_file():
        raise FileNotFoundError(f"输入 ctl 文件不存在: {in_path}")

    iter_parts = [x.strip() for x in mcmc_iter.split(",") if x.strip()]
    if len(iter_parts) != 3:
        raise ValueError("mcmc_iter 必须为 burnin,sampfreq,nsample")

    bds_parts = [x.strip() for x in mcmc_bds.split(",") if x.strip()]
    if len(bds_parts) != 3:
        raise ValueError("mcmc_bds 必须为 birth,death,sampling_fraction")

    clock_code = CLOCK_MAP.get(clock_model.upper())
    if clock_code is None:
        raise ValueError("clock_model 仅支持 EQUAL/IND/CORR")

    values = {
        "seqfile": seqfile,
        "treefile": treefile,
        "mcmcfile": mcmcfile,
        "outfile": outfile,
        "ckpfile": ckpfile,
        "hessianfile": hessianfile,
        "clock": clock_code,
        "BDparas": " ".join(bds_parts),
        "burnin": iter_parts[0],
        "sampfreq": iter_parts[1],
        "nsample": iter_parts[2],
        "checkpoint": "2" if resume else str(checkpoint),
    }

    lines_out = []
    seen = set()

    for raw in in_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("*") or "=" not in stripped:
            lines_out.append(raw)
            continue

        key = stripped.split("=", 1)[0].strip()
        if key not in KEY_VALUE_FIELDS:
            lines_out.append(raw)
            continue

        left, comment = _split_comment(raw)
        if key in values:
            new_line = f"{key} = {values[key]}"
            if comment:
                new_line = f"{new_line}  {comment}"
            lines_out.append(new_line)
            seen.add(key)
        else:
            lines_out.append(raw)

    for key in [
        "seqfile",
        "treefile",
        "mcmcfile",
        "outfile",
        "ckpfile",
        "hessianfile",
        "checkpoint",
        "clock",
        "BDparas",
        "burnin",
        "sampfreq",
        "nsample",
    ]:
        if key not in seen:
            lines_out.append(f"{key} = {values[key]}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines_out) + "\n", encoding="utf-8")

    LOG.info("ctl 文件已更新: %s", out_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare mcmctree ctl")
    parser.add_argument("--input-ctl", required=True)
    parser.add_argument("--output-ctl", required=True)
    parser.add_argument("--seqfile", required=True)
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--mcmcfile", required=True)
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--ckpfile", required=True)
    parser.add_argument("--hessianfile", required=True)
    parser.add_argument("--clock-model", required=True)
    parser.add_argument("--mcmc-iter", required=True)
    parser.add_argument("--mcmc-bds", required=True)
    parser.add_argument("--checkpoint", required=True, type=int)
    parser.add_argument("--resume", action="store_true")
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
        input_ctl=args.input_ctl,
        output_ctl=args.output_ctl,
        seqfile=args.seqfile,
        treefile=args.treefile,
        mcmcfile=args.mcmcfile,
        outfile=args.outfile,
        ckpfile=args.ckpfile,
        hessianfile=args.hessianfile,
        clock_model=args.clock_model,
        mcmc_iter=args.mcmc_iter,
        mcmc_bds=args.mcmc_bds,
        checkpoint=args.checkpoint,
        resume=args.resume,
    )


if __name__ == "__main__":
    raise SystemExit(main())

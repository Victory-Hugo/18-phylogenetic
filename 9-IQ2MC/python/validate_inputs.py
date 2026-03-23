#!/usr/bin/env python3
"""验证 IQ2MC 流水线输入文件并准备输出目录。"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Dict

LOG = logging.getLogger(__name__)


def _resolve_path(project_root: Path, value: str) -> Path:
    p = Path(value)
    if p.is_absolute():
        return p
    return (project_root / p).resolve()


def run(
    project_root: str,
    msa_file: str,
    rooted_calibrated_tree: str,
    partition_file: str,
    output_dir: str,
) -> Dict[str, Any]:
    """执行输入校验，返回标准化路径。"""
    root = Path(project_root).resolve()
    if not root.exists():
        raise FileNotFoundError(f"项目根目录不存在: {root}")

    msa_path = _resolve_path(root, msa_file)
    tree_path = _resolve_path(root, rooted_calibrated_tree)
    output_path = _resolve_path(root, output_dir)

    if not msa_path.is_file():
        raise FileNotFoundError(f"MSA 文件不存在: {msa_path}")
    if not tree_path.is_file():
        raise FileNotFoundError(f"有根校准树文件不存在: {tree_path}")

    partition_path = None
    partition_file = partition_file.strip()
    if partition_file:
        partition_path = _resolve_path(root, partition_file)
        if not partition_path.is_file():
            raise FileNotFoundError(f"分区文件不存在: {partition_path}")

    output_path.mkdir(parents=True, exist_ok=True)

    result = {
        "project_root": str(root),
        "msa_file": str(msa_path),
        "rooted_calibrated_tree": str(tree_path),
        "partition_file": str(partition_path) if partition_path else "",
        "output_dir": str(output_path),
    }
    LOG.info("输入校验通过")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate IQ2MC pipeline inputs")
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--msa-file", required=True)
    parser.add_argument("--rooted-calibrated-tree", required=True)
    parser.add_argument("--partition-file", default="")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--json-out", default="")
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
        project_root=args.project_root,
        msa_file=args.msa_file,
        rooted_calibrated_tree=args.rooted_calibrated_tree,
        partition_file=args.partition_file,
        output_dir=args.output_dir,
    )

    payload = json.dumps(result, ensure_ascii=False, indent=2)
    if args.json_out:
        Path(args.json_out).write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env bash
set -euo pipefail

# 中文说明：检查基础命令与关键二进制是否可用。
require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd" >&2
        exit 1
    fi
}

require_exec_file() {
    local path="$1"
    local label="$2"
    if [[ ! -f "$path" ]]; then
        echo "[ERROR] Missing $label file: $path" >&2
        exit 1
    fi
    if [[ ! -x "$path" ]]; then
        echo "[ERROR] $label is not executable: $path" >&2
        exit 1
    fi
}

CONDA_BIN="${1:-}"
IQTREE3_BIN="${2:-}"
MCMCTREE_BIN="${3:-}"

require_cmd bash
require_cmd awk
require_cmd sed

if [[ -z "$CONDA_BIN" || -z "$IQTREE3_BIN" || -z "$MCMCTREE_BIN" ]]; then
    echo "[ERROR] Usage: check_env.sh <conda_bin> <iqtree3_bin> <mcmctree_bin>" >&2
    exit 1
fi

require_exec_file "$CONDA_BIN" "conda"
require_exec_file "$IQTREE3_BIN" "iqtree3"
require_exec_file "$MCMCTREE_BIN" "mcmctree"

echo "[OK] Environment check passed."

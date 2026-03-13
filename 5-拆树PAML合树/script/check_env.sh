#!/usr/bin/env bash
set -euo pipefail

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd" >&2
        exit 1
    fi
}

require_cmd bash
require_cmd python3
require_cmd conda

if ! conda run -n BigLin bash -lc 'command -v gotree >/dev/null 2>&1'; then
    echo "[ERROR] gotree is not available in conda environment: BigLin" >&2
    exit 1
fi

if ! conda run -n BigLin gotree stats --help >/dev/null 2>&1; then
    echo "[ERROR] Failed to execute 'conda run -n BigLin gotree stats --help'" >&2
    exit 1
fi

echo "[OK] Environment check passed."

#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN=""
CONDA_ENV=""
GOTREE_BIN=""
BASEML_BIN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --gotree)
            GOTREE_BIN="$2"
            shift 2
            ;;
        --baseml)
            BASEML_BIN="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd" >&2
        exit 1
    fi
}

require_executable() {
    local path_or_cmd="$1"
    if [[ -z "$path_or_cmd" ]]; then
        echo "[ERROR] Executable path is empty." >&2
        exit 1
    fi
    if [[ "$path_or_cmd" == */* ]]; then
        if [[ ! -x "$path_or_cmd" ]]; then
            echo "[ERROR] Executable is not runnable: $path_or_cmd" >&2
            exit 1
        fi
        return
    fi
    require_cmd "$path_or_cmd"
}

require_cmd bash

if [[ -z "$PYTHON_BIN" ]]; then
    PYTHON_BIN="python3"
fi
require_executable "$PYTHON_BIN"

if [[ -n "$CONDA_ENV" || -n "$GOTREE_BIN" ]]; then
    require_cmd conda
    if [[ -z "$CONDA_ENV" ]]; then
        echo "[ERROR] --conda-env is required when --gotree is provided." >&2
        exit 1
    fi
    if [[ -z "$GOTREE_BIN" ]]; then
        echo "[ERROR] --gotree is required when --conda-env is provided." >&2
        exit 1
    fi

    if ! conda run -n "$CONDA_ENV" bash -lc "command -v '$GOTREE_BIN' >/dev/null 2>&1"; then
        echo "[ERROR] gotree is not available in conda environment: $CONDA_ENV" >&2
        exit 1
    fi

    if ! conda run -n "$CONDA_ENV" "$GOTREE_BIN" stats --help >/dev/null 2>&1; then
        echo "[ERROR] Failed to execute 'conda run -n $CONDA_ENV $GOTREE_BIN stats --help'" >&2
        exit 1
    fi
fi

if [[ -n "$BASEML_BIN" && "$BASEML_BIN" != "null" ]]; then
    require_executable "$BASEML_BIN"
    set +e
    "$BASEML_BIN" </dev/null >/dev/null 2>&1
    status=$?
    set -e
    if [[ "$status" -ne 0 && "$status" -ne 1 && "$status" -ne 255 ]]; then
        echo "[ERROR] baseml probe exited unexpectedly with status $status" >&2
        exit 1
    fi
fi

echo "[OK] Environment check passed."

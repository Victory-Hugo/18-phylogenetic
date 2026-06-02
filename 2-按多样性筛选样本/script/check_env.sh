#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="${1:-$PROJECT_ROOT/conf/Config.yaml}"

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd" >&2
        exit 1
    fi
}

require_executable() {
    local exe_path="$1"
    if [[ ! -x "$exe_path" ]]; then
        echo "[ERROR] Required executable is not runnable: $exe_path" >&2
        exit 1
    fi
}

require_cmd bash
require_cmd awk

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

if [[ -z "${CFG_TOOLS_PYTHON_BIN:-}" ]]; then
    echo "[ERROR] Missing required config value: CFG_TOOLS_PYTHON_BIN" >&2
    exit 1
fi

if [[ "$CFG_TOOLS_PYTHON_BIN" = /* ]]; then
    require_executable "$CFG_TOOLS_PYTHON_BIN"
else
    require_cmd "$CFG_TOOLS_PYTHON_BIN"
fi

echo "[OK] Environment check passed."

#!/usr/bin/env bash
set -euo pipefail

# 中文说明：自动检查并创建 conda 环境 IQ2MC。
CONDA_BIN=""
ENV_NAME="IQ2MC"
CREATE_IF_MISSING="true"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --conda-bin)
            CONDA_BIN="$2"
            shift 2
            ;;
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --create-if-missing)
            CREATE_IF_MISSING="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$CONDA_BIN" ]]; then
    echo "[ERROR] Missing --conda-bin" >&2
    exit 1
fi

if [[ ! -x "$CONDA_BIN" ]]; then
    echo "[ERROR] conda binary is not executable: $CONDA_BIN" >&2
    exit 1
fi

if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[OK] Conda environment exists: $ENV_NAME"
    exit 0
fi

if [[ "${CREATE_IF_MISSING,,}" != "true" ]]; then
    echo "[ERROR] Conda environment '$ENV_NAME' not found and auto-create is disabled." >&2
    exit 1
fi

echo "[INFO] Creating conda environment: $ENV_NAME"
"$CONDA_BIN" create -y -n "$ENV_NAME" python=3.11 pip

echo "[OK] Conda environment created: $ENV_NAME"

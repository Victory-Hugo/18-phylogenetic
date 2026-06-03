#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/console_ui.sh"
ui_init

PYTHON_BIN=""
CONDA_ENV=""
GOTREE_BIN=""
BEAST_BIN=""
TREEANNOTATOR_BIN=""

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
        --beast)
            BEAST_BIN="$2"
            shift 2
            ;;
        --treeannotator)
            TREEANNOTATOR_BIN="$2"
            shift 2
            ;;
        *)
            ui_error "参数错误 | Unknown argument: $1"
            exit 1
            ;;
    esac
done

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        ui_error "环境缺失 | Required command not found: $cmd"
        exit 1
    fi
}

require_executable() {
    local path_or_cmd="$1"
    if [[ -z "$path_or_cmd" ]]; then
        ui_error "环境缺失 | Executable path is empty."
        exit 1
    fi
    if [[ "$path_or_cmd" == */* ]]; then
        if [[ ! -x "$path_or_cmd" ]]; then
            ui_error "环境缺失 | Executable is not runnable: $path_or_cmd"
            exit 1
        fi
        ui_ok "Executable ready | $path_or_cmd"
        return
    fi
    require_cmd "$path_or_cmd"
    ui_ok "Command ready | $path_or_cmd"
}

ui_section "Checking environment" "Python / conda / gotree / beast / treeannotator"

require_cmd bash
ui_ok "Command ready | bash"

if [[ -z "$PYTHON_BIN" ]]; then
    PYTHON_BIN="python3"
fi
require_executable "$PYTHON_BIN"

if [[ -n "$CONDA_ENV" || -n "$GOTREE_BIN" ]]; then
    require_cmd conda
    ui_ok "Command ready | conda"
    if [[ -z "$CONDA_ENV" ]]; then
        ui_error "参数错误 | --conda-env is required when --gotree is provided."
        exit 1
    fi
    if [[ -z "$GOTREE_BIN" ]]; then
        ui_error "参数错误 | --gotree is required when --conda-env is provided."
        exit 1
    fi

    # 用非登录 shell（bash -c）探测：部分服务器的登录 shell（bash -lc）会在 profile 中
    # 重置 PATH，丢掉 conda 环境的 bin 目录，导致已安装的 gotree 误判为缺失。
    if ! conda run -n "$CONDA_ENV" bash -c "command -v '$GOTREE_BIN' >/dev/null 2>&1"; then
        ui_error "环境缺失 | gotree is not available in conda environment: $CONDA_ENV"
        exit 1
    fi
    ui_ok "Conda env ready | $CONDA_ENV provides $GOTREE_BIN"

    if ! conda run -n "$CONDA_ENV" "$GOTREE_BIN" stats --help >/dev/null 2>&1; then
        ui_error "环境探测失败 | Failed to execute 'conda run -n $CONDA_ENV $GOTREE_BIN stats --help'"
        exit 1
    fi
    ui_ok "Probe passed | gotree stats --help"
fi

if [[ -n "$BEAST_BIN" && "$BEAST_BIN" != "null" ]]; then
    require_executable "$BEAST_BIN"
    set +e
    "$BEAST_BIN" -version </dev/null >/dev/null 2>&1
    status=$?
    set -e
    if [[ "$status" -ne 0 && "$status" -ne 1 && "$status" -ne 255 ]]; then
        ui_error "环境探测失败 | beast probe exited unexpectedly with status $status"
        exit 1
    fi
    ui_ok "Probe passed | beast -version exited with accepted status $status"
fi

if [[ -n "$TREEANNOTATOR_BIN" && "$TREEANNOTATOR_BIN" != "null" ]]; then
    require_executable "$TREEANNOTATOR_BIN"
    set +e
    "$TREEANNOTATOR_BIN" -version </dev/null >/dev/null 2>&1
    status=$?
    set -e
    if [[ "$status" -ne 0 && "$status" -ne 1 && "$status" -ne 255 ]]; then
        ui_error "环境探测失败 | treeannotator probe exited unexpectedly with status $status"
        exit 1
    fi
    ui_ok "Probe passed | treeannotator -version exited with accepted status $status"
fi

ui_ok "Completed | Environment check passed"

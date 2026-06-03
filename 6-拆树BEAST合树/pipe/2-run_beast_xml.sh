#!/usr/bin/env bash

# 2-run_beast_xml.sh
# 运行已准备好的 BEAST XML，并解析 TreeAnnotator MCC 输出。

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo

PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/2-beast.yaml"
CONFIG_PATH="$DEFAULT_CONFIG_PATH"
BOOTSTRAP_PYTHON="python3"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        *)
            ui_error "参数错误 | Unknown argument: $1"
            exit 1
            ;;
    esac
done

resolve_path() {
    local value="$1"
    if [[ -z "$value" || "$value" == "null" ]]; then
        printf '%s\n' "$value"
        return
    fi
    if [[ "$value" == /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$PATH_ROOT/$value"
    fi
}

resolve_command_or_path() {
    local value="$1"
    if [[ "$value" == */* || "$value" == .* ]]; then
        resolve_path "$value"
    else
        printf '%s\n' "$value"
    fi
}

if [[ ! -f "$CONFIG_PATH" ]]; then
    ui_error "配置缺失 | Config file not found: $CONFIG_PATH"
    exit 1
fi

CONFIG_LOADER="$PROJECT_ROOT/python/config_loader.py"
config_get() { "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1"; }
config_get_default() { "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1" --default "$2"; }

CONFIG_PROJECT_ROOT=$(config_get_default projectpath "$PROJECT_ROOT")
PATH_ROOT=$(resolve_path "$CONFIG_PROJECT_ROOT")

PYTHON_BIN=$(resolve_command_or_path "$(config_get tools.python)")
BEAST_BIN=$(resolve_command_or_path "$(config_get tools.beast)")
TREEANNOTATOR_BIN=$(resolve_command_or_path "$(config_get tools.treeannotator)")
BEAST_EXECUTION_MODE=$(config_get_default beast.execution_mode real)

if [[ "$BEAST_EXECUTION_MODE" == "fake" ]]; then
    ui_error "模式不匹配 | execution_mode=fake does not run XML. Set beast.execution_mode=real."
    exit 1
fi

bash "$PROJECT_ROOT/script/check_env.sh" \
    --python "$PYTHON_BIN" \
    --beast "$BEAST_BIN" \
    --treeannotator "$TREEANNOTATOR_BIN"

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)
BEAST_OUTPUT_DIR=$(resolve_path "$(config_get beast.output_dir)")
BURNIN_PERCENT=$(config_get beast.burnin_percent)
MCC_HEIGHTS=$(config_get beast.mcc_heights)
BEAST_SEED=$(config_get_default beast.seed 42)
BEAST_PARALLEL_JOBS=$(config_get beast.parallel_jobs)
BEAST_THREADS_PER_JOB=$(config_get beast.threads_per_job)
BEAST_BIND_CPU_AFFINITY=$(config_get_default beast.bind_cpu_affinity true)
BEAST_SKIP_EXISTING=$(config_get beast.skip_existing)
BACKBONE_ONLY_ENABLED=$(config_get beast.backbone_only_enabled)
BACKBONE_JOB_DIR=$(resolve_path "$(config_get beast.backbone_job_dir)")
BACKBONE_ANALYSIS_DIR=$(resolve_path "$(config_get beast.backbone_analysis_dir)")
BACKBONE_SKIP_EXISTING=$(config_get beast.backbone_skip_existing)
ULTRAMETRIC_TOLERANCE=$(config_get beast.ultrametric_tolerance)
MIN_BRANCH_LENGTH=$(config_get_default beast.min_branch_length 1e-8)

BEAST_SUBTREE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/beast_subtree_summary.tsv"
MANIFEST_PATH="$BEAST_OUTPUT_DIR/beast_job_manifest.tsv"
BEAST_RESUME_LOG_DIR="$BEAST_OUTPUT_DIR/log"
BACKBONE_TREE="$SPLIT_OUTPUT_DIR/backbone_tree.nwk"

if ! [[ "$BEAST_THREADS_PER_JOB" =~ ^[1-9][0-9]*$ ]]; then
    ui_error "参数错误 | beast.threads_per_job must be a positive integer: $BEAST_THREADS_PER_JOB"
    exit 1
fi

for required in "$BEAST_SUBTREE_SUMMARY_TSV" "$MANIFEST_PATH" "$BACKBONE_TREE"; do
    if [[ ! -f "$required" ]]; then
        ui_error "输入缺失 | Required BEAST run input not found: $required"
        exit 1
    fi
done

if [[ "$BACKBONE_ONLY_ENABLED" == "True" || "$BACKBONE_ONLY_ENABLED" == "true" ]]; then
    if [[ ! -f "$BACKBONE_JOB_DIR/beast.xml" ]]; then
        ui_error "输入缺失 | Backbone XML not found; run prepare stage first: $BACKBONE_JOB_DIR/beast.xml"
        exit 1
    fi
fi

ui_section "BEAST XML run" "Step 2B/5 | 运行 BEAST XML 并解析 MCC"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Output dir" "$BEAST_OUTPUT_DIR"
ui_kv "Manifest" "$MANIFEST_PATH"
ui_kv "Parallel jobs" "$BEAST_PARALLEL_JOBS"
ui_kv "Threads/job" "$BEAST_THREADS_PER_JOB"

RUN_ARGS=(
    --manifest "$MANIFEST_PATH"
    --beast-bin "$BEAST_BIN"
    --treeannotator-bin "$TREEANNOTATOR_BIN"
    --parallel-jobs "$BEAST_PARALLEL_JOBS"
    --threads-per-job "$BEAST_THREADS_PER_JOB"
    --burnin-percent "$BURNIN_PERCENT"
    --mcc-heights "$MCC_HEIGHTS"
    --base-seed "$BEAST_SEED"
    --resume-log-dir "$BEAST_RESUME_LOG_DIR"
    --log-level "$LOG_LEVEL"
)
if [[ "$BEAST_SKIP_EXISTING" == "True" || "$BEAST_SKIP_EXISTING" == "true" ]]; then
    RUN_ARGS+=(--legacy-skip-existing)
fi
if [[ "$BEAST_BIND_CPU_AFFINITY" == "True" || "$BEAST_BIND_CPU_AFFINITY" == "true" ]]; then
    RUN_ARGS+=(--bind-cpu-affinity)
fi

ui_stage_start "Stage 1/4" "run_beast_jobs | 执行子树 BEAST + TreeAnnotator"
"$PYTHON_BIN" "$PROJECT_ROOT/python/run_beast_jobs.py" "${RUN_ARGS[@]}"
ui_stage_end "Stage 1/4" "run_beast_jobs | 执行子树 BEAST + TreeAnnotator"

ui_stage_start "Stage 2/4" "parse_beast_outputs | 解析子树 MCC 树"
"$PYTHON_BIN" "$PROJECT_ROOT/python/parse_beast_outputs.py" \
    --beast-summary-tsv "$BEAST_SUBTREE_SUMMARY_TSV" \
    --manifest "$MANIFEST_PATH" \
    --output-dir "$BEAST_OUTPUT_DIR" \
    --min-branch-length "$MIN_BRANCH_LENGTH" \
    --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE" \
    --log-level "$LOG_LEVEL"
ui_stage_end "Stage 2/4" "parse_beast_outputs | 解析子树 MCC 树"

if [[ "$BACKBONE_ONLY_ENABLED" == "True" || "$BACKBONE_ONLY_ENABLED" == "true" ]]; then
    BACKBONE_RUN_ARGS=(
        --job-dir "$BACKBONE_JOB_DIR"
        --beast-bin "$BEAST_BIN"
        --treeannotator-bin "$TREEANNOTATOR_BIN"
        --burnin-percent "$BURNIN_PERCENT"
        --mcc-heights "$MCC_HEIGHTS"
        --seed "$BEAST_SEED"
        --threads "$BEAST_THREADS_PER_JOB"
        --log-level "$LOG_LEVEL"
    )
    if [[ "$BACKBONE_SKIP_EXISTING" == "True" || "$BACKBONE_SKIP_EXISTING" == "true" ]]; then
        BACKBONE_RUN_ARGS+=(--skip-existing)
    fi
    if [[ "$BEAST_BIND_CPU_AFFINITY" == "True" || "$BEAST_BIND_CPU_AFFINITY" == "true" ]]; then
        BACKBONE_RUN_ARGS+=(--bind-cpu-affinity)
    fi

    ui_stage_start "Stage 3/4" "run_backbone_beast_job | 执行 backbone BEAST"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/run_backbone_beast_job.py" "${BACKBONE_RUN_ARGS[@]}"
    ui_stage_end "Stage 3/4" "run_backbone_beast_job | 执行 backbone BEAST"

    ui_stage_start "Stage 4/4" "parse_backbone_beast_output | 解析 backbone 输出"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/parse_backbone_beast_output.py" \
        --backbone-tree "$BACKBONE_TREE" \
        --beast-summary-tsv "$BEAST_SUBTREE_SUMMARY_TSV" \
        --job-dir "$BACKBONE_JOB_DIR" \
        --analysis-dir "$BACKBONE_ANALYSIS_DIR" \
        --min-branch-length "$MIN_BRANCH_LENGTH" \
        --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE" \
        --log-level "$LOG_LEVEL"
    ui_stage_end "Stage 4/4" "parse_backbone_beast_output | 解析 backbone 输出"
fi

ui_ok "Completed | BEAST XML run/parse finished"

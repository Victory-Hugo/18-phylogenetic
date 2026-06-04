#!/usr/bin/env bash

# 脚本名称: 6-run_multinode_calib_pipeline.sh
# 描述: 多节点时间校准 + 基于子树离散度的 Monte-Carlo 年代区间（不重算 PAML）
#
# 功能:
#   - 从配置文件读取多节点校准参数（N/M 单系目标年龄、root 年龄、MC 参数）
#   - 以拆树/合并阶段已有的子树尺度离散度为基础，估计每个节点年龄的置信区间
#   - 输出点估计年化树、节点年龄区间表与校准摘要
#
# 用法:
#   ./6-run_multinode_calib_pipeline.sh [--config CONFIG_PATH]
#
# 依赖:
#   - python3 / config_loader.py / calibrate_multinode_interval.py / check_env.sh
#   - 配置文件: conf/6-multinode_calib.yaml
#   - 输入: output/ultrastandard/merged_ultrametric_tree.nwk,
#           output/merge/subtree_scale_report.tsv,
#           output/split/target_subtree_summary.tsv,
#           data/1-ID-Hap.tsv, data/2-phylotree_index_v17.2.json
#
# 输出文件 (output/multinode_calib/):
#   - merged_multinode_tree_years.nwk     多节点校准后的点估计年化树
#   - node_age_ci.tsv                     每个内部节点的年龄与 95% 区间
#   - multinode_calibration_summary.tsv   校准与 MC 摘要
#   - multinode_calibration.log           日志

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/6-multinode_calib.yaml"
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

config_get() {
    "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1"
}

config_get_optional() {
    "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1" --default "${2:-}"
}

CONFIG_PROJECT_ROOT=$(config_get_optional projectpath "$PROJECT_ROOT")
PATH_ROOT=$(resolve_path "$CONFIG_PROJECT_ROOT")

PYTHON_BIN=$(resolve_command_or_path "$(config_get tools.python)")
bash "$PROJECT_ROOT/script/check_env.sh" --python "$PYTHON_BIN"

ULTRA_OUTPUT_DIR=$(resolve_path "$(config_get paths.ultrastandard_output_dir)")
MERGE_OUTPUT_DIR=$(resolve_path "$(config_get paths.merge_output_dir)")
SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

MC_INPUT_TREE=$(config_get_optional multinode_calib.input_tree)
ID_HAP_TSV=$(resolve_path "$(config_get multinode_calib.id_hap_tsv)")
PHYLOTREE_JSON=$(resolve_path "$(config_get multinode_calib.phylotree_json)")
MC_OUTPUT_DIR=$(resolve_path "$(config_get multinode_calib.output_dir)")
CLEAN_OUTPUT_DIR=$(config_get multinode_calib.clean_output_dir)
CALIBRATIONS=$(config_get multinode_calib.calibrations)
ROOT_AGE=$(config_get multinode_calib.root_age)
REPLICATES=$(config_get_optional multinode_calib.replicates 1000)
SEED=$(config_get_optional multinode_calib.seed 42)
USE_OFFSET=$(config_get_optional multinode_calib.use_systematic_offset false)
CI_LOWER=$(config_get_optional multinode_calib.ci_lower 2.5)
CI_UPPER=$(config_get_optional multinode_calib.ci_upper 97.5)
CHUNK_SIZE=$(config_get_optional multinode_calib.chunk_size 250)
OUTPUT_TREE_NAME=$(config_get_optional multinode_calib.output_tree_name merged_multinode_tree_years.nwk)

if [[ "$MC_INPUT_TREE" == "null" || -z "$MC_INPUT_TREE" ]]; then
    INPUT_TREE="$ULTRA_OUTPUT_DIR/merged_ultrametric_tree.nwk"
else
    INPUT_TREE=$(resolve_path "$MC_INPUT_TREE")
fi

if [[ ! -f "$INPUT_TREE" ]]; then
    ui_error "输入缺失 | Required input tree not found: $INPUT_TREE"
    exit 1
fi

ui_section "Multinode calibration pipeline" "多节点校准 + 子树离散度区间"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Input tree" "$INPUT_TREE"
ui_kv "Output dir" "$MC_OUTPUT_DIR"
ui_kv "Calibrations" "$CALIBRATIONS; root=$ROOT_AGE"
ui_kv "MC replicates" "$REPLICATES"

mkdir -p "$MC_OUTPUT_DIR"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$MC_OUTPUT_DIR"/multinode_calibration.log
    rm -f "$MC_OUTPUT_DIR"/multinode_calibration_summary.tsv
    rm -f "$MC_OUTPUT_DIR"/node_age_ci.tsv
    rm -f "$MC_OUTPUT_DIR"/"$OUTPUT_TREE_NAME"
fi

OFFSET_FLAG=()
if [[ "$USE_OFFSET" == "True" || "$USE_OFFSET" == "true" ]]; then
    OFFSET_FLAG=(--use-systematic-offset)
fi

ui_stage_start "Stage 1/1" "calibrate_multinode_interval | 多节点校准与区间"
"$PYTHON_BIN" "$PROJECT_ROOT/python/calibrate_multinode_interval.py" \
    --input-tree "$INPUT_TREE" \
    --id-hap-tsv "$ID_HAP_TSV" \
    --phylotree-json "$PHYLOTREE_JSON" \
    --merge-output-dir "$MERGE_OUTPUT_DIR" \
    --split-output-dir "$SPLIT_OUTPUT_DIR" \
    --output-dir "$MC_OUTPUT_DIR" \
    --calibrations "$CALIBRATIONS" \
    --root-age "$ROOT_AGE" \
    --replicates "$REPLICATES" \
    --seed "$SEED" \
    --ci-lower "$CI_LOWER" \
    --ci-upper "$CI_UPPER" \
    --chunk-size "$CHUNK_SIZE" \
    --output-tree-name "$OUTPUT_TREE_NAME" \
    "${OFFSET_FLAG[@]}" \
    --log-level "$LOG_LEVEL"
ui_stage_end "Stage 1/1" "calibrate_multinode_interval | 多节点校准与区间"

ui_ok "Completed | Multinode calibration pipeline finished"

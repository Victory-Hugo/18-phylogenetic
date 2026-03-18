#!/usr/bin/env bash

# 脚本名称: 5-run_time_calib_pipeline.sh
# 描述: 系统发育树时间校准管道脚本
# 
# 功能:
#   - 从配置文件中读取时间校准相关参数
#   - 验证输入树文件的有效性
#   - 执行树的分子钟校准
#   - 生成时间校准后的超度量树
#
# 用法:
#   ./5-run_time_calib_pipeline.sh [--config CONFIG_PATH]
#
# 参数:
#   --config CONFIG_PATH      指定配置文件路径（可选，默认为 conf/5-time_calib.yaml）
#
# 依赖:
#   - python3: Python运行环境
#   - config_loader.py: 配置文件解析工具
#   - time_calibrate_tree.py: 时间校准核心模块
#   - check_env.sh: 环境检查脚本
#   - 配置文件: conf/5-time_calib.yaml
#
# 环境变量:
#   - BOOTSTRAP_PYTHON: Python3解释器路径（默认: python3）
#
# 配置参数说明:
#   - projectpath: 项目根路径
#   - tools.python: Python可执行文件路径
#   - paths.ultrastandard_output_dir: 超度量树输出目录
#   - runtime.log_level: 日志级别
#   - time_calib.input_tree: 输入树文件路径（可选）
#   - time_calib.output_dir: 输出目录
#   - time_calib.clean_output_dir: 是否清空输出目录
#   - time_calib.method: 校准方法（默认: molecular_clock）
#   - time_calib.substitution_rate_per_site_per_year: 每位点每年替代率
#   - time_calib.sequence_length: 序列长度
#   - time_calib.branch_length_unit: 分支长度单位（默认: substitutions_per_site）
#   - time_calib.node_calibration_tip_name: 校准节点名称（默认: RSRS）
#   - time_calib.node_calibration_divergence_years: 校准分化时间年数（默认: 180000）
#   - time_calib.output_tree_name: 输出树文件名（默认: merged_ultrametric_tree_years.nwk）
#
# 输出文件:
#   - time_calibration.log: 校准日志文件
#   - time_calibration_summary.tsv: 校准摘要表
#   - time_calibration_edge_report.tsv: 分支报告表
#   - merged_ultrametric_tree_years.nwk: 时间校准后的树文件
#
# 错误处理:
#   - 脚本在任何命令失败时立即退出
#   - 配置文件不存在时报错并退出
#   - 输入树文件不存在时报错并退出

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/5-time_calib.yaml"
CONFIG_PATH="$DEFAULT_CONFIG_PATH"
BOOTSTRAP_PYTHON="python3"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
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
    echo "[ERROR] Config file not found: $CONFIG_PATH" >&2
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
LOG_LEVEL=$(config_get runtime.log_level)

TIME_CALIB_INPUT_TREE=$(config_get_optional time_calib.input_tree)
TIME_CALIB_OUTPUT_DIR=$(resolve_path "$(config_get time_calib.output_dir)")
CLEAN_OUTPUT_DIR=$(config_get time_calib.clean_output_dir)
TIME_CALIB_METHOD=$(config_get_optional time_calib.method molecular_clock)
SUBSTITUTION_RATE=$(config_get time_calib.substitution_rate_per_site_per_year)
SEQUENCE_LENGTH=$(config_get time_calib.sequence_length)
BRANCH_LENGTH_UNIT=$(config_get_optional time_calib.branch_length_unit substitutions_per_site)
NODE_CALIBRATION_TIP_NAME=$(config_get_optional time_calib.node_calibration_tip_name RSRS)
NODE_CALIBRATION_DIVERGENCE_YEARS=$(config_get_optional time_calib.node_calibration_divergence_years 180000)
OUTPUT_TREE_NAME=$(config_get_optional time_calib.output_tree_name merged_ultrametric_tree_years.nwk)

if [[ "$TIME_CALIB_INPUT_TREE" == "null" || -z "$TIME_CALIB_INPUT_TREE" ]]; then
    INPUT_TREE="$ULTRA_OUTPUT_DIR/merged_ultrametric_tree.nwk"
else
    INPUT_TREE=$(resolve_path "$TIME_CALIB_INPUT_TREE")
fi

if [[ ! -f "$INPUT_TREE" ]]; then
    echo "[ERROR] Required time calibration input not found: $INPUT_TREE" >&2
    exit 1
fi

mkdir -p "$TIME_CALIB_OUTPUT_DIR"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration.log
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration_summary.tsv
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration_edge_report.tsv
    rm -f "$TIME_CALIB_OUTPUT_DIR"/"$OUTPUT_TREE_NAME"
fi

echo "[INFO] Stage 1/1: time_calibrate_tree"
"$PYTHON_BIN" "$PROJECT_ROOT/python/time_calibrate_tree.py" \
    --input-tree "$INPUT_TREE" \
    --output-dir "$TIME_CALIB_OUTPUT_DIR" \
    --output-tree-name "$OUTPUT_TREE_NAME" \
    --method "$TIME_CALIB_METHOD" \
    --substitution-rate-per-site-per-year "$SUBSTITUTION_RATE" \
    --sequence-length "$SEQUENCE_LENGTH" \
    --branch-length-unit "$BRANCH_LENGTH_UNIT" \
    --node-calibration-tip-name "$NODE_CALIBRATION_TIP_NAME" \
    --node-calibration-divergence-years "$NODE_CALIBRATION_DIVERGENCE_YEARS" \
    --log-level "$LOG_LEVEL"

echo "[OK] Time calibration pipeline finished."

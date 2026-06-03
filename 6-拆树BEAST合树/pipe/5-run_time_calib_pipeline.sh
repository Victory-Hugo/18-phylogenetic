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
#   - time_calib.beast_output_dir: stage2 BEAST 输出目录（用于定位节点 HPD 表，默认: output/beast）
#   - time_calib.nexus_output_name: BEAST 风格注释 NEXUS 文件名（默认: merged_time_tree.nexus）
#   - time_calib.node_table_name: 全节点时间信息表文件名（默认: node_time_table.tsv）
#
# 输出文件:
#   - time_calibration.log: 校准日志文件
#   - time_calibration_summary.tsv: 校准摘要表
#   - time_calibration_edge_report.tsv: 分支报告表
#   - merged_ultrametric_tree_years.nwk: 时间校准后的树文件（纯 Newick，年）
#   - merged_time_tree.nexus: BEAST/FigTree 风格年代树（每节点带 height 与 95% HPD 注释）
#   - node_time_table.tsv: 全节点时间信息表（年龄/根深/HPD/来源/枝长/代表 tip）
#
# 错误处理:
#   - 脚本在任何命令失败时立即退出
#   - 配置文件不存在时报错并退出
#   - 输入树文件不存在时报错并退出

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo
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
BEAST_OUTPUT_DIR=$(resolve_path "$(config_get_optional time_calib.beast_output_dir output/beast)")
NEXUS_OUTPUT_NAME=$(config_get_optional time_calib.nexus_output_name merged_time_tree.nexus)
NODE_TABLE_NAME=$(config_get_optional time_calib.node_table_name node_time_table.tsv)
# 节点高度 95% HPD 表（substitutions/site）：backbone 在前、各 target 子树在后，
# 使局部采样更密的 target 子树覆盖 backbone。stage5 读取后按 clade 线性传播到年。
BACKBONE_NODE_HPD="$BEAST_OUTPUT_DIR/backbone_analysis/backbone_node_hpd.tsv"
SUBTREE_NODE_HPD="$BEAST_OUTPUT_DIR/beast_node_hpd.tsv"

if [[ "$TIME_CALIB_INPUT_TREE" == "null" || -z "$TIME_CALIB_INPUT_TREE" ]]; then
    INPUT_TREE="$ULTRA_OUTPUT_DIR/merged_ultrametric_tree.nwk"
else
    INPUT_TREE=$(resolve_path "$TIME_CALIB_INPUT_TREE")
fi

if [[ ! -f "$INPUT_TREE" ]]; then
    ui_error "输入缺失 | Required time calibration input not found: $INPUT_TREE"
    exit 1
fi

ui_section "Time calibration pipeline" "Step 5/5 | 分子钟校准与年代尺度输出"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Input tree" "$INPUT_TREE"
ui_kv "Output dir" "$TIME_CALIB_OUTPUT_DIR"
ui_kv "Method" "$TIME_CALIB_METHOD"
ui_kv "Output tree" "$OUTPUT_TREE_NAME"

mkdir -p "$TIME_CALIB_OUTPUT_DIR"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration.log
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration_summary.tsv
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration_edge_report.tsv
    rm -f "$TIME_CALIB_OUTPUT_DIR"/"$OUTPUT_TREE_NAME"
    rm -f "$TIME_CALIB_OUTPUT_DIR"/"$NEXUS_OUTPUT_NAME"
    rm -f "$TIME_CALIB_OUTPUT_DIR"/"$NODE_TABLE_NAME"
fi

# 仅当存在的 HPD 表才传入（fake 模式无 MCC 时这些文件不存在，stage5 自动退化为无 HPD NEXUS）。
NODE_HPD_ARGS=()
if [[ -f "$BACKBONE_NODE_HPD" ]]; then
    NODE_HPD_ARGS+=(--beast-node-hpd "$BACKBONE_NODE_HPD")
fi
if [[ -f "$SUBTREE_NODE_HPD" ]]; then
    NODE_HPD_ARGS+=(--beast-node-hpd "$SUBTREE_NODE_HPD")
fi

ui_stage_start "Stage 1/1" "time_calibrate_tree | 时间校准"
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
    --nexus-output-name "$NEXUS_OUTPUT_NAME" \
    --node-table-name "$NODE_TABLE_NAME" \
    "${NODE_HPD_ARGS[@]}" \
    --log-level "$LOG_LEVEL"
ui_stage_end "Stage 1/1" "time_calibrate_tree | 时间校准"

ui_ok "Completed | Time calibration pipeline finished"

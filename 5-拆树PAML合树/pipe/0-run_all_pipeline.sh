#!/usr/bin/env bash

# 0-run_all_pipeline.sh - 全流程调度脚本（文档）
#
# 用途：
#   按顺序运行本项目中预定义的 5 个子流程：
#     1) split
#     2) paml
#     3) merge
#     4) ultrastandard
#     5) time calibration
#   每一步均通过独立的脚本执行，脚本运行失败则整个流程中止。
#
# 工作方式：
#   - 脚本自动确定自身所在目录并把项目根目录设为 SCRIPT_DIR/..
#   - 使用一组 YAML 配置文件指定各步骤的参数（默认路径在 conf/ 下）。
#   - 逐步校验所有配置文件是否存在，若任一文件缺失则退出（非零状态码）。
#   - 以 run_stage 函数包装每一步：打印开始/结束信息并执行对应子脚本。
#   - set -euo pipefail：遇到错误、未定义变量或管道失败即退出。
#
# 默认配置文件（可通过命令行覆盖）：
#   conf/1-split.yaml           - 第 1 步 split 配置
#   conf/2-paml.yaml            - 第 2 步 paml 配置
#   conf/3-merge.yaml           - 第 3 步 merge 配置
#   conf/4-ultrastandard.yaml   - 第 4 步 ultrastandard 配置
#   conf/5-time_calib.yaml      - 第 5 步时间校准配置
#
# 命令行选项：
#   --split-config PATH
#   --paml-config PATH
#   --merge-config PATH
#   --ultrastandard-config PATH
#   --time-calib-config PATH
#   --help, -h
#   未知参数会导致脚本打印错误并以非零状态退出。
#
# 前提与依赖：
#   - 需要 bash（POSIX bash 支持数组与函数）。
#   - 项目目录下存在 pipe/1-run_split_pipeline.sh ... pipe/5-run_time_calib_pipeline.sh。
#   - 配置文件为合法的 YAML，且各子脚本能正确解析对应配置。
#
# 输出与日志：
#   - 在标准输出打印每一步的开始与完成标记（[INFO] / [OK]）。
#   - 错误信息打印到标准错误并伴随非零退出码。
#
# 使用示例：
#   bash pipe/0-run_all_pipeline.sh
#   bash pipe/0-run_all_pipeline.sh --paml-config /path/to/custom_paml.yaml
#
# 退出码：
#   0   - 所有步骤成功完成
#   1   - 找不到配置文件、参数错误或任一步骤失败
#
# 维护说明：
#   - 若增加或删除子流程，请同步修改此脚本中步骤顺序与配置项。
#   - 可在子脚本中增加更细粒度的日志或并行化策略，但要保持本脚本对失败的即时响应。

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo

SPLIT_CONFIG="$PROJECT_ROOT/conf/1-split.yaml"
PAML_CONFIG="$PROJECT_ROOT/conf/2-paml.yaml"
MERGE_CONFIG="$PROJECT_ROOT/conf/3-merge.yaml"
ULTRASTANDARD_CONFIG="$PROJECT_ROOT/conf/4-ultrastandard.yaml"
TIME_CALIB_CONFIG="$PROJECT_ROOT/conf/5-time_calib.yaml"


usage() {
    cat <<EOF
Usage:
  bash pipe/0-run_all_pipeline.sh [options]

Options:
  --split-config PATH          Config for step 1 split
  --paml-config PATH           Config for step 2 paml
  --merge-config PATH          Config for step 3 merge
  --ultrastandard-config PATH  Config for step 4 ultrastandard
  --time-calib-config PATH     Config for step 5 time calibration
  --help                       Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --split-config)
            SPLIT_CONFIG="$2"
            shift 2
            ;;
        --paml-config)
            PAML_CONFIG="$2"
            shift 2
            ;;
        --merge-config)
            MERGE_CONFIG="$2"
            shift 2
            ;;
        --ultrastandard-config)
            ULTRASTANDARD_CONFIG="$2"
            shift 2
            ;;
        --time-calib-config)
            TIME_CALIB_CONFIG="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            ui_error "参数错误 | Unknown argument: $1"
            usage >&2
            exit 1
            ;;
    esac
done

for config_path in \
    "$SPLIT_CONFIG" \
    "$PAML_CONFIG" \
    "$MERGE_CONFIG" \
    "$ULTRASTANDARD_CONFIG" \
    "$TIME_CALIB_CONFIG"; do
    if [[ ! -f "$config_path" ]]; then
        ui_error "配置缺失 | Config file not found: $config_path"
        exit 1
    fi
done

run_stage() {
    local stage_step="$1"
    local stage_title="$2"
    shift
    shift
    ui_stage_start "$stage_step" "$stage_title"
    "$@"
    ui_stage_end "$stage_step" "$stage_title"
}

ui_section "Phylogenetic Pipeline Console UI" "Full workflow | split -> paml -> merge -> ultrastandard -> time calibration"
ui_kv "Split config" "$SPLIT_CONFIG"
ui_kv "PAML config" "$PAML_CONFIG"
ui_kv "Merge config" "$MERGE_CONFIG"
ui_kv "Ultra config" "$ULTRASTANDARD_CONFIG"
ui_kv "Time config" "$TIME_CALIB_CONFIG"

run_stage "Step 1/5" "split | 拆树与校验" \
    bash "$PROJECT_ROOT/pipe/1-run_split_pipeline.sh" \
    --config "$SPLIT_CONFIG"

run_stage "Step 2/5" "paml | 准备与运行 baseml" \
    bash "$PROJECT_ROOT/pipe/2-run_paml_pipeline.sh" \
    --config "$PAML_CONFIG"

run_stage "Step 3/5" "merge | 合并子树结果" \
    bash "$PROJECT_ROOT/pipe/3-run_merge_pipeline.sh" \
    --config "$MERGE_CONFIG"

run_stage "Step 4/5" "ultrastandard | 构建超度量树" \
    bash "$PROJECT_ROOT/pipe/4-run_ultrastandard_pipeline.sh" \
    --config "$ULTRASTANDARD_CONFIG"

run_stage "Step 5/5" "time calibration | 时间校准" \
    bash "$PROJECT_ROOT/pipe/5-run_time_calib_pipeline.sh" \
    --config "$TIME_CALIB_CONFIG"

ui_section "Workflow completed" "All configured stages finished successfully."
ui_summary_line "Result" "Full pipeline finished"

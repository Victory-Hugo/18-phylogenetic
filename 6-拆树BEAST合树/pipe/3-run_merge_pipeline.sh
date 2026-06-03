#!/usr/bin/env bash

# 3-run_merge_pipeline.sh
#
# 功能
#   根据配置文件执行合并（merge）流水线：校验环境、可选的模拟基线结果、合并 PAML 子树结果并验证合并后的树。
#
# 用法
#   ./3-run_merge_pipeline.sh [--config /path/to/config.yaml]
#
# 参数
#   --config    指定配置文件路径（可为相对或绝对路径）。默认：conf/3-merge.yaml（相对于项目根目录）。
#
# 配置项（来自 YAML，使用 python/config_loader.py 读取）
#   projectpath                       （可选）项目根路径，用于解析相对路径
#   tools.python                      Python 可执行文件（可为命令或路径）
#   paths.split_output_dir            拆分步骤输出目录（必需，包含多种 summary 和子树）
#   paths.beast_output_dir             PAML 输出基础目录
#   merge.output_dir                  合并结果输出目录（将被创建）
#   merge.analysis_tree_source        one of: external|simulated|...（决定是否模拟）
#   merge.external_result_dir         （可选）外部结果目录
#   merge.clean_output_dir            是否清理旧输出（"True"/"true" 表示清理）
#   merge.mode                        合并模式（默认 "backbone_graft"）
#   merge.backbone_tree               主干树文件（可为自定义路径）
#   merge.redline_tolerance           redline 检查容忍度（默认 1e-6）
#   merge.parallel_jobs               并行作业数（默认 8）
#   merge.randomization_model         模拟时使用的随机化模型
#   merge.randomization_sigma         模拟时的 sigma 参数
#   merge.randomization_seed          模拟时的随机种子
#   merge.min_branch_length           最短分支长度阈值
#   merge.outgroup_tip_name           （可选）外群 tip 名称，用于模拟或验证
#   runtime.log_level                 日志等级（传递给下游脚本）
#
# 主要行为
#   1. 解析并规范化路径（相对路径以 projectpath 或脚本根路径为基准）。
#   2. 使用 python/config_loader.py 读取配置；默认使用 python3 启动该 loader。
#   3. 调用 script/check_env.sh 验证运行时环境（传递 --python）。
#   4. 检查必需输入文件是否存在：
#        - $SPLIT_OUTPUT_DIR/intermediate/rooted.tree
#        - $SPLIT_OUTPUT_DIR/backbone_summary.tsv
#        - $SPLIT_OUTPUT_DIR/target_subtree_summary.tsv
#        - $SPLIT_OUTPUT_DIR/beast_subtree_summary.tsv
#        - $SPLIT_OUTPUT_DIR/beast_tree_manifest.tsv
#        - $BACKBONE_TREE
#      若缺失则以错误退出。
#   5. 根据配置创建输出目录：merge.output_dir 及其子目录 intermediate、simulated_baseml_subtrees。
#   6. 若 merge.clean_output_dir 为真，则删除若干旧的合并产物文件以便重跑。
#   7. 若 merge.analysis_tree_source == "simulated"，调用 python/simulate_baseml_results.py 生成模拟的 Baseml 结果（可传 outgroup）。
#   8. 调用 python/merge_baseml_subtrees_redline.py 执行合并（传入分裂输出目录、backbone_tree、redline 容忍度、并行数等）。
#   9. 调用 python/validate_merged_tree.py 验证合并结果（传入 merge.mode、outgroup 等）。
#  10. 成功完成后输出 "[OK] Merge pipeline finished."
#
# 退出码
#   0   成功完成
#   1   常见错误：找不到配置文件、必需的拆分输出文件缺失、未知脚本参数或下游脚本返回错误（脚本使用 set -euo pipefail）
#
# 依赖（项目内脚本/模块）
#   - python/config_loader.py
#   - script/check_env.sh
#   - python/simulate_baseml_results.py
#   - python/merge_baseml_subtrees_redline.py
#   - python/validate_merged_tree.py
#
# 注意
#   - 相对路径会基于配置中的 projectpath（若存在）或脚本所在的项目根目录进行解析。
#   - config_loader 使用独立的 Python（脚本内默认为 python3）加载配置；实际 pipeline 脚本使用配置中指定的 tools.python。
#   - 本脚本严格模式运行（set -euo pipefail），遇到未捕获错误即退出。

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/3-merge.yaml"
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

bash "$PROJECT_ROOT/script/check_env.sh" \
    --python "$PYTHON_BIN"

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
PAML_OUTPUT_DIR=$(resolve_path "$(config_get paths.beast_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

MERGE_OUTPUT_DIR=$(resolve_path "$(config_get merge.output_dir)")
ANALYSIS_TREE_SOURCE=$(config_get merge.analysis_tree_source)
EXTERNAL_RESULT_DIR=$(config_get_optional merge.external_result_dir)
CLEAN_OUTPUT_DIR=$(config_get merge.clean_output_dir)
MERGE_MODE=$(config_get_optional merge.mode "backbone_graft")
BACKBONE_TREE=$(resolve_path "$(config_get_optional merge.backbone_tree "$PAML_OUTPUT_DIR/backbone_analysis/backbone_ultrametric_tree.nwk")")
REDLINE_TOLERANCE=$(config_get_optional merge.redline_tolerance 1e-6)
PARALLEL_JOBS=$(config_get_optional merge.parallel_jobs 8)
RANDOMIZATION_MODEL=$(config_get merge.randomization_model)
RANDOMIZATION_SIGMA=$(config_get merge.randomization_sigma)
RANDOMIZATION_SEED=$(config_get merge.randomization_seed)
MIN_BRANCH_LENGTH=$(config_get merge.min_branch_length)
OUTGROUP_TIP_NAME=$(config_get_optional merge.outgroup_tip_name)

if [[ "$EXTERNAL_RESULT_DIR" != "null" && -n "$EXTERNAL_RESULT_DIR" ]]; then
    EXTERNAL_RESULT_DIR_ABS=$(resolve_path "$EXTERNAL_RESULT_DIR")
elif [[ "$ANALYSIS_TREE_SOURCE" == "external" ]]; then
    EXTERNAL_RESULT_DIR_ABS="$PAML_OUTPUT_DIR/analysis_trees"
else
    EXTERNAL_RESULT_DIR_ABS=""
fi

for required in \
    "$SPLIT_OUTPUT_DIR/intermediate/rooted.tree" \
    "$SPLIT_OUTPUT_DIR/backbone_summary.tsv" \
    "$SPLIT_OUTPUT_DIR/target_subtree_summary.tsv" \
    "$SPLIT_OUTPUT_DIR/beast_subtree_summary.tsv" \
    "$SPLIT_OUTPUT_DIR/beast_tree_manifest.tsv" \
    "$BACKBONE_TREE"; do
    if [[ ! -f "$required" ]]; then
        ui_error "输入缺失 | Required split output not found: $required"
        exit 1
    fi
done

ui_section "Merge pipeline" "Step 3/5 | 合并 PAML 子树并验证结果"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Source" "$ANALYSIS_TREE_SOURCE"
ui_kv "Split dir" "$SPLIT_OUTPUT_DIR"
ui_kv "PAML dir" "$PAML_OUTPUT_DIR"
ui_kv "Output dir" "$MERGE_OUTPUT_DIR"

mkdir -p "$MERGE_OUTPUT_DIR" "$MERGE_OUTPUT_DIR/intermediate" "$MERGE_OUTPUT_DIR/simulated_baseml_subtrees"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$MERGE_OUTPUT_DIR"/assembly_scaffold.nwk
    rm -f "$MERGE_OUTPUT_DIR"/backbone_edge_estimates.tsv
    rm -f "$MERGE_OUTPUT_DIR"/graft_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/edge_update_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/subtree_scale_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/merge_validation_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/merged_ml_tree.nwk
    rm -f "$MERGE_OUTPUT_DIR"/merged_tree.nwk
    rm -f "$MERGE_OUTPUT_DIR"/redline_failed_tip_depths.tsv
    rm -f "$MERGE_OUTPUT_DIR"/simulated_baseml_subtrees/*.nwk
fi

if [[ "$ANALYSIS_TREE_SOURCE" == "simulated" ]]; then
    ui_stage_start "Stage 1/3" "simulate_baseml_results | 生成模拟分析树"
    SIMULATE_ARGS=(
        --paml-summary-tsv "$SPLIT_OUTPUT_DIR/beast_subtree_summary.tsv"
        --paml-tree-dir "$SPLIT_OUTPUT_DIR/beast_subtrees"
        --output-dir "$MERGE_OUTPUT_DIR"
        --randomization-model "$RANDOMIZATION_MODEL"
        --randomization-sigma "$RANDOMIZATION_SIGMA"
        --randomization-seed "$RANDOMIZATION_SEED"
        --min-branch-length "$MIN_BRANCH_LENGTH"
        --log-level "$LOG_LEVEL"
    )
    if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
        SIMULATE_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
    fi
    "$PYTHON_BIN" "$PROJECT_ROOT/python/simulate_baseml_results.py" "${SIMULATE_ARGS[@]}"
    ui_stage_end "Stage 1/3" "simulate_baseml_results | 生成模拟分析树"
else
    ui_info "Stage 1/3 | skip simulated results because analysis_tree_source=$ANALYSIS_TREE_SOURCE"
fi

MERGE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --merge-output-dir "$MERGE_OUTPUT_DIR"
    --analysis-tree-source "$ANALYSIS_TREE_SOURCE"
    --backbone-tree "$BACKBONE_TREE"
    --min-branch-length "$MIN_BRANCH_LENGTH"
    --redline-tolerance "$REDLINE_TOLERANCE"
    --parallel-jobs "$PARALLEL_JOBS"
)

if [[ -n "$EXTERNAL_RESULT_DIR_ABS" ]]; then
    MERGE_ARGS+=(--external-result-dir "$EXTERNAL_RESULT_DIR_ABS")
fi

ui_stage_start "Stage 2/3" "merge_baseml_subtrees_redline | 合并子树"
"$PYTHON_BIN" "$PROJECT_ROOT/python/merge_baseml_subtrees_redline.py" "${MERGE_ARGS[@]}"
ui_stage_end "Stage 2/3" "merge_baseml_subtrees_redline | 合并子树"

VALIDATE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --merge-output-dir "$MERGE_OUTPUT_DIR"
    --merge-mode "$MERGE_MODE"
    --log-level "$LOG_LEVEL"
)

if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
    VALIDATE_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
fi

ui_stage_start "Stage 3/3" "validate_merged_tree | 校验合并树"
"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_merged_tree.py" "${VALIDATE_ARGS[@]}"
ui_stage_end "Stage 3/3" "validate_merged_tree | 校验合并树"

ui_ok "Completed | Merge pipeline finished"

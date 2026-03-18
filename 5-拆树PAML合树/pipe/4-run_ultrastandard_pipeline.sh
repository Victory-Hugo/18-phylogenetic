#!/usr/bin/env bash
################################################################################
# 脚本名称: 4-run_ultrastandard_pipeline.sh
# 功能描述: 
#   执行超度量树构建的完整流程，包括三个主要阶段：
#   1. 合并超度量树 (merge_ultrastandard_tree)
#   2. 验证超度量树 (validate_ultrastandard_tree)
#   3. 比较最终拓扑结构 (check_tree_topology_match)
#
# 用法:
#   ./4-run_ultrastandard_pipeline.sh [--config CONFIG_FILE]
#
# 参数:
#   --config CONFIG_FILE  指定配置文件路径（可选，默认为 conf/4-ultrastandard.yaml）
#
# 依赖文件:
#   - conf/4-ultrastandard.yaml         超度量树配置文件
#   - python/config_loader.py           配置加载器
#   - python/merge_ultrastandard_tree.py 树合并模块
#   - python/validate_ultrastandard_tree.py 树验证模块
#   - python/check_tree_topology_match.py 拓扑比较模块
#   - script/check_env.sh               环境检查脚本
#
# 必需输入文件:
#   - split_output_dir/intermediate/rooted.tree
#   - split_output_dir/paml_subtree_summary.tsv
#   - merge_output_dir/merged_ml_tree.nwk
#
# 输出文件:
#   - ultrastandard_output_dir/assembly_scaffold.nwk
#   - ultrastandard_output_dir/backbone_assigned_scaffold.nwk
#   - ultrastandard_output_dir/merged_ultrametric_tree.nwk
#   - ultrastandard_output_dir/*.tsv （各类报告文件）
#
# 配置项:
#   tools.python                        Python 解释器路径
#   paths.split_output_dir              拆分输出目录
#   paths.merge_output_dir              合并输出目录
#   runtime.log_level                   日志级别
#   ultrastandard.output_dir            超度量树输出目录
#   ultrastandard.clean_output_dir      是否清理输出目录
#   ultrastandard.min_branch_length     最小分支长度
#   ultrastandard.ultrametric_tolerance 超度量容差值
#   ultrastandard.projection_mode       投影模式
#   ultrastandard.primary_tree_input    主要树输入（可选）
#
# 返回值:
#   0  成功完成所有阶段
#   1  配置文件缺失或必需输入文件缺失
#   2  拓扑结构不匹配但继续执行
#   其他  执行过程中出现错误
#
# 错误处理:
#   - 使用 set -euo pipefail 确保错误立即退出
#   - 验证所有必需输入文件存在
#   - 捕获拓扑比较的返回状态进行判断
################################################################################
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/4-ultrastandard.yaml"
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

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
MERGE_OUTPUT_DIR=$(resolve_path "$(config_get paths.merge_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

ULTRA_OUTPUT_DIR=$(resolve_path "$(config_get ultrastandard.output_dir)")
CLEAN_OUTPUT_DIR=$(config_get ultrastandard.clean_output_dir)
MIN_BRANCH_LENGTH=$(config_get ultrastandard.min_branch_length)
ULTRAMETRIC_TOLERANCE=$(config_get ultrastandard.ultrametric_tolerance)
PROJECTION_MODE=$(config_get ultrastandard.projection_mode)
PRIMARY_TREE_INPUT=$(config_get_optional ultrastandard.primary_tree_input)

for required in \
    "$SPLIT_OUTPUT_DIR/intermediate/rooted.tree" \
    "$SPLIT_OUTPUT_DIR/paml_subtree_summary.tsv" \
    "$MERGE_OUTPUT_DIR/merged_ml_tree.nwk"; do
    if [[ ! -f "$required" ]]; then
        ui_error "输入缺失 | Required ultrastandard input not found: $required"
        exit 1
    fi
done

ui_section "Ultrastandard pipeline" "Step 4/5 | 合并与验证超度量树"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Split dir" "$SPLIT_OUTPUT_DIR"
ui_kv "Merge dir" "$MERGE_OUTPUT_DIR"
ui_kv "Output dir" "$ULTRA_OUTPUT_DIR"
ui_kv "Projection" "$PROJECTION_MODE"

mkdir -p "$ULTRA_OUTPUT_DIR"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$ULTRA_OUTPUT_DIR"/assembly_scaffold.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/backbone_assigned_scaffold.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/backbone_edge_assignment_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/graft_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/projection_input_tree.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/target_scaling_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/ultrametric_projection_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/root_to_tip_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/merged_ultrametric_tree.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/merged_tree.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/ultrastandard_validation_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/ultrastandard.log
    rm -rf "$ULTRA_OUTPUT_DIR/relative_target_trees"
fi

MERGE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --merge-output-dir "$MERGE_OUTPUT_DIR"
    --ultrastandard-output-dir "$ULTRA_OUTPUT_DIR"
    --min-branch-length "$MIN_BRANCH_LENGTH"
    --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE"
    --projection-mode "$PROJECTION_MODE"
    --log-level "$LOG_LEVEL"
)

if [[ "$PRIMARY_TREE_INPUT" != "null" && -n "$PRIMARY_TREE_INPUT" ]]; then
    MERGE_ARGS+=(--primary-tree-input "$(resolve_path "$PRIMARY_TREE_INPUT")")
fi

ui_stage_start "Stage 1/3" "merge_ultrastandard_tree | 构建超度量树"
"$PYTHON_BIN" "$PROJECT_ROOT/python/merge_ultrastandard_tree.py" "${MERGE_ARGS[@]}"
ui_stage_end "Stage 1/3" "merge_ultrastandard_tree | 构建超度量树"

VALIDATE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --ultrastandard-output-dir "$ULTRA_OUTPUT_DIR"
    --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE"
    --log-level "$LOG_LEVEL"
)

ui_stage_start "Stage 2/3" "validate_ultrastandard_tree | 校验超度量树"
"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_ultrastandard_tree.py" "${VALIDATE_ARGS[@]}"
ui_stage_end "Stage 2/3" "validate_ultrastandard_tree | 校验超度量树"

ui_stage_start "Stage 3/3" "compare_final_topology | 比较最终拓扑"
if "$PYTHON_BIN" "$PROJECT_ROOT/python/check_tree_topology_match.py" \
    --reference-tree "$SPLIT_OUTPUT_DIR/intermediate/rooted.tree" \
    --query-tree "$ULTRA_OUTPUT_DIR/merged_ultrametric_tree.nwk" \
    --log-level "$LOG_LEVEL"; then
    ui_ok "Topology matched | Ultrastandard merged tree topology matches the original rooted tree"
else
    compare_status=$?
    if [[ "$compare_status" -eq 2 ]]; then
        ui_warn "拓扑不一致 | Ultrastandard merged tree topology does not match the original rooted tree"
    else
        exit "$compare_status"
    fi
fi
ui_stage_end "Stage 3/3" "compare_final_topology | 比较最终拓扑"

ui_ok "Completed | Ultrastandard pipeline finished"

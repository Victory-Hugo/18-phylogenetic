#!/usr/bin/env bash

# 脚本说明：
#   1-run_split_pipeline.sh
#   用途：基于输入系统发育树拆分为若干子树，生成用于 PAML 等下游分析的子树与校验报告。
#
# 用法：
#   ./1-run_split_pipeline.sh [--config /path/to/config.yaml]
#   默认配置文件：conf/1-split.yaml（相对于项目根目录）
#
# 主要工作流程（3 个阶段）：
#   1) precompute_split_context：预计算拆分所需的上下文（precompute_split_context.py）
#   2) split_phylo_tree：按照配置拆分主树并生成子树与清单（split_phylo_tree.py）
#   3) validate_phylo_split：对拆分结果进行校验并输出报告（validate_phylo_split.py）
#
# 关键行为：
#   - 解析并规范化配置与路径（支持相对/绝对路径）。
#   - 调用 script/check_env.sh 验证 python、conda 环境与 gotree 工具可用性。
#   - 创建输出目录：output_dir、output_dir/target_subtrees、output_dir/beast_subtrees、output_dir/intermediate。
#   - 当配置中 clean_output_dir 为 true 时，清理历史输出文件（若存在）。
#   - 将配置中的各项通过命令行参数传递给三个 Python 程序，并使用指定的 python 可执行文件运行。
#
# 主要配置键（config.yaml 中常用项）：
#   projectpath
#   tools.python                 -> Python 可执行路径或命令
#   tools.conda_env              -> Conda 环境名称（用于 check_env.sh）
#   tools.gotree                 -> gotree 可执行路径或命令
#   paths.input_tree             -> 输入树文件（Newick）
#   paths.outgroup_tip_file      -> 外群（outgroup）tip id 文件
#   paths.output_dir             -> 输出根目录
#   paths.backbone_tree          -> 可选：骨干树文件
#   paths.backbone_tip_id_file   -> 可选：骨干 tip id 文件
#   runtime.max_tips
#   runtime.backbone_size
#   runtime.threads
#   runtime.outgroup_tip_name
#   runtime.clean_output_dir
#   runtime.log_level
#   runtime.backbone_sampling_strategy
#   runtime.target_partition_mode
#   runtime.local_anchor_count
#   runtime.anchor_selection_strategy
#   runtime.benchmark_tree_tips
#   runtime.pre_binarize_rooted_tree
#   runtime.validation_mode
#
# 主要输出（位于 output_dir）：
#   - target_subtrees/target_subtree_*.nwk
#   - beast_subtrees/beast_subtree_*.nwk
#   - intermediate/split_context.json
#   - backbone_tree.nwk、backbone_tips.txt（若采用骨干抽样）
#   - 多个 summary/manifest/validation tsv 文件（subtree_design_summary.tsv、beast_subtree_summary.tsv 等）
#   - 运行日志：split_tree.log、split_precompute.log、split_validation_report.tsv
#
# 错误与退出：
#   - 当缺少配置文件、必要工具不可用或任一步骤失败时脚本会以非 0 状态退出（set -euo pipefail）。
#
# 备注：
#   - 默认使用 python3；可在配置中通过 tools.python 指定不同的 python 可执行文件或路径。
#   - 所有相对路径会基于配置中的 projectpath 解析为绝对路径。
#   - 若需调试，调整 config 中 runtime.log_level 或观察生成的日志/TSV 文件。

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/1-split.yaml"
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
CONDA_ENV=$(config_get tools.conda_env)
GOTREE_BIN=$(resolve_command_or_path "$(config_get tools.gotree)")

bash "$PROJECT_ROOT/script/check_env.sh" \
    --python "$PYTHON_BIN" \
    --conda-env "$CONDA_ENV" \
    --gotree "$GOTREE_BIN"

INPUT_TREE=$(resolve_path "$(config_get paths.input_tree)")
OUTGROUP_TIP_FILE=$(resolve_path "$(config_get paths.outgroup_tip_file)")
OUTPUT_DIR=$(resolve_path "$(config_get paths.output_dir)")
BACKBONE_TREE=$(config_get_optional paths.backbone_tree)
BACKBONE_TIP_ID_FILE=$(config_get_optional paths.backbone_tip_id_file)
MAX_TIPS=$(config_get runtime.max_tips)
BACKBONE_SIZE=$(config_get runtime.backbone_size)
THREADS=$(config_get runtime.threads)
OUTGROUP_TIP_NAME=$(config_get_optional runtime.outgroup_tip_name)
CLEAN_SPLIT_OUTPUT_DIR=$(config_get runtime.clean_output_dir)
LOG_LEVEL=$(config_get runtime.log_level)
BACKBONE_SAMPLING_STRATEGY=$(config_get runtime.backbone_sampling_strategy)
TARGET_PARTITION_MODE=$(config_get runtime.target_partition_mode)
LOCAL_ANCHOR_COUNT=$(config_get_optional runtime.local_anchor_count 16)
ANCHOR_SELECTION_STRATEGY=$(config_get_optional runtime.anchor_selection_strategy nearest_boundary_patristic)
BENCHMARK_TREE_TIPS=$(config_get_optional runtime.benchmark_tree_tips 300)
PRE_BINARIZE_ROOTED_TREE=$(config_get_optional runtime.pre_binarize_rooted_tree true)
VALIDATION_MODE=$(config_get_optional runtime.validation_mode fast)

ui_section "Split pipeline" "Step 1/5 | 拆树、上下文预计算与校验"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Input tree" "$INPUT_TREE"
ui_kv "Output dir" "$OUTPUT_DIR"
ui_kv "Threads" "$THREADS"
ui_kv "Max tips" "$MAX_TIPS"

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/target_subtrees" "$OUTPUT_DIR/beast_subtrees" "$OUTPUT_DIR/intermediate"

if [[ "$CLEAN_SPLIT_OUTPUT_DIR" == "True" || "$CLEAN_SPLIT_OUTPUT_DIR" == "true" ]]; then
    rm -f "$OUTPUT_DIR"/target_subtrees/target_subtree_*.nwk
    rm -f "$OUTPUT_DIR"/beast_subtrees/beast_subtree_*.nwk
    rm -f "$OUTPUT_DIR"/backbone_tips.txt
    rm -f "$OUTPUT_DIR"/backbone_tree.nwk
    rm -f "$OUTPUT_DIR"/backbone_summary.tsv
    rm -f "$OUTPUT_DIR"/target_subtree_summary.tsv
    rm -f "$OUTPUT_DIR"/target_tree_manifest.tsv
    rm -f "$OUTPUT_DIR"/subtree_design_summary.tsv
    rm -f "$OUTPUT_DIR"/anchor_manifest.tsv
    rm -f "$OUTPUT_DIR"/beast_subtree_summary.tsv
    rm -f "$OUTPUT_DIR"/beast_tree_manifest.tsv
    rm -f "$OUTPUT_DIR"/split_tree.log
    rm -f "$OUTPUT_DIR"/split_precompute.log
    rm -f "$OUTPUT_DIR"/split_validation_report.tsv
    rm -f "$OUTPUT_DIR"/intermediate/rooted.tree
    rm -f "$OUTPUT_DIR"/intermediate/rooted.validation.tree
    rm -f "$OUTPUT_DIR"/intermediate/split_context.json
fi

PRECOMPUTE_ARGS=(
    --input-tree "$INPUT_TREE"
    --outgroup-tip-file "$OUTGROUP_TIP_FILE"
    --output-dir "$OUTPUT_DIR"
    --threads "$THREADS"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --pre-binarize-rooted-tree "$PRE_BINARIZE_ROOTED_TREE"
    --log-level "$LOG_LEVEL"
)

if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
    PRECOMPUTE_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
fi

ui_stage_start "Stage 1/3" "precompute_split_context | 预计算拆树上下文"
"$PYTHON_BIN" "$PROJECT_ROOT/python/precompute_split_context.py" "${PRECOMPUTE_ARGS[@]}"
ui_stage_end "Stage 1/3" "precompute_split_context | 预计算拆树上下文"

SPLIT_ARGS=(
    --input-tree "$INPUT_TREE"
    --outgroup-tip-file "$OUTGROUP_TIP_FILE"
    --output-dir "$OUTPUT_DIR"
    --max-tips "$MAX_TIPS"
    --backbone-size "$BACKBONE_SIZE"
    --threads "$THREADS"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --backbone-sampling-strategy "$BACKBONE_SAMPLING_STRATEGY"
    --target-partition-mode "$TARGET_PARTITION_MODE"
    --local-anchor-count "$LOCAL_ANCHOR_COUNT"
    --anchor-selection-strategy "$ANCHOR_SELECTION_STRATEGY"
    --benchmark-tree-tips "$BENCHMARK_TREE_TIPS"
    --pre-binarize-rooted-tree "$PRE_BINARIZE_ROOTED_TREE"
    --precomputed-context-json "$OUTPUT_DIR/intermediate/split_context.json"
    --log-level "$LOG_LEVEL"
)

if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
    SPLIT_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
fi

if [[ "$BACKBONE_TREE" != "null" && -n "$BACKBONE_TREE" ]]; then
    SPLIT_ARGS+=(--backbone-tree "$(resolve_path "$BACKBONE_TREE")")
fi

if [[ "$BACKBONE_TIP_ID_FILE" != "null" && -n "$BACKBONE_TIP_ID_FILE" ]]; then
    SPLIT_ARGS+=(--backbone-tip-id-file "$(resolve_path "$BACKBONE_TIP_ID_FILE")")
fi

ui_stage_start "Stage 2/3" "split_phylo_tree | 生成 target/beast 子树"
"$PYTHON_BIN" "$PROJECT_ROOT/python/split_phylo_tree.py" "${SPLIT_ARGS[@]}"
ui_stage_end "Stage 2/3" "split_phylo_tree | 生成 target/beast 子树"

VALIDATE_ARGS=(
    --input-tree "$INPUT_TREE"
    --outgroup-tip-file "$OUTGROUP_TIP_FILE"
    --output-dir "$OUTPUT_DIR"
    --max-tips "$MAX_TIPS"
    --backbone-size "$BACKBONE_SIZE"
    --threads "$THREADS"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --backbone-sampling-strategy "$BACKBONE_SAMPLING_STRATEGY"
    --target-partition-mode "$TARGET_PARTITION_MODE"
    --local-anchor-count "$LOCAL_ANCHOR_COUNT"
    --anchor-selection-strategy "$ANCHOR_SELECTION_STRATEGY"
    --benchmark-tree-tips "$BENCHMARK_TREE_TIPS"
    --validation-mode "$VALIDATION_MODE"
    --precomputed-context-json "$OUTPUT_DIR/intermediate/split_context.json"
    --pre-binarize-rooted-tree "$PRE_BINARIZE_ROOTED_TREE"
    --log-level "$LOG_LEVEL"
)

if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
    VALIDATE_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
fi

if [[ "$BACKBONE_TREE" != "null" && -n "$BACKBONE_TREE" ]]; then
    VALIDATE_ARGS+=(--backbone-tree "$(resolve_path "$BACKBONE_TREE")")
fi

if [[ "$BACKBONE_TIP_ID_FILE" != "null" && -n "$BACKBONE_TIP_ID_FILE" ]]; then
    VALIDATE_ARGS+=(--backbone-tip-id-file "$(resolve_path "$BACKBONE_TIP_ID_FILE")")
fi

ui_stage_start "Stage 3/3" "validate_phylo_split | 校验拆树结果"
"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_phylo_split.py" "${VALIDATE_ARGS[@]}"
ui_stage_end "Stage 3/3" "validate_phylo_split | 校验拆树结果"

ui_ok "Completed | Split pipeline finished"

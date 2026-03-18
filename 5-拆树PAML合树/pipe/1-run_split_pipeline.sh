#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
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

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/target_subtrees" "$OUTPUT_DIR/paml_subtrees" "$OUTPUT_DIR/intermediate"

if [[ "$CLEAN_SPLIT_OUTPUT_DIR" == "True" || "$CLEAN_SPLIT_OUTPUT_DIR" == "true" ]]; then
    rm -f "$OUTPUT_DIR"/target_subtrees/target_subtree_*.nwk
    rm -f "$OUTPUT_DIR"/paml_subtrees/paml_subtree_*.nwk
    rm -f "$OUTPUT_DIR"/backbone_tips.txt
    rm -f "$OUTPUT_DIR"/backbone_tree.nwk
    rm -f "$OUTPUT_DIR"/backbone_summary.tsv
    rm -f "$OUTPUT_DIR"/target_subtree_summary.tsv
    rm -f "$OUTPUT_DIR"/target_tree_manifest.tsv
    rm -f "$OUTPUT_DIR"/subtree_design_summary.tsv
    rm -f "$OUTPUT_DIR"/anchor_manifest.tsv
    rm -f "$OUTPUT_DIR"/paml_subtree_summary.tsv
    rm -f "$OUTPUT_DIR"/paml_tree_manifest.tsv
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

echo "[INFO] Stage 0/2: precompute_split_context"
"$PYTHON_BIN" "$PROJECT_ROOT/python/precompute_split_context.py" "${PRECOMPUTE_ARGS[@]}"

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

"$PYTHON_BIN" "$PROJECT_ROOT/python/split_phylo_tree.py" "${SPLIT_ARGS[@]}"

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

"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_phylo_split.py" "${VALIDATE_ARGS[@]}"

echo "[OK] Split pipeline finished."

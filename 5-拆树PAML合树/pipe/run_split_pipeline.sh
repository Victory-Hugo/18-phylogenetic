#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="$PROJECT_ROOT/conf/Config.json"
BOOTSTRAP_PYTHON="python3"

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "[ERROR] Config file not found: $CONFIG_PATH" >&2
    exit 1
fi

bash "$PROJECT_ROOT/script/check_env.sh"

CONFIG_LOADER="$PROJECT_ROOT/python/config_loader.py"
PYTHON_BIN=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key tools.python)
INPUT_TREE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.input_tree)
OUTGROUP_TIP_FILE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.outgroup_tip_file)
OUTPUT_DIR=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.output_dir)
BACKBONE_TREE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.backbone_tree)
MAX_TIPS=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.max_tips)
MIN_BASEML_TIPS=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.min_baseml_tips)
THREADS=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.threads)
ENABLE_MERGE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.enable_merge_small_blocks)
CLEAN_SUBTREE_DIR=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.clean_subtree_dir)
LOG_LEVEL=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.log_level)
CONSTRUCT_BASEML=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.construct_baseml_subtrees)
ALWAYS_INCLUDE_OUTGROUP=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.always_include_outgroup)
ANCHOR_STRATEGY=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.anchor_strategy)
RESERVE_SLOTS_FOR_OUTGROUP=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.reserve_slots_for_outgroup)
CONDA_ENV=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key tools.conda_env)
GOTREE_BIN=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key tools.gotree)
BACKBONE_MODE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key backbone.mode)

mkdir -p "$PROJECT_ROOT/$OUTPUT_DIR" "$PROJECT_ROOT/$OUTPUT_DIR/core_subtrees" "$PROJECT_ROOT/$OUTPUT_DIR/baseml_subtrees" "$PROJECT_ROOT/$OUTPUT_DIR/intermediate"

if [[ "$CLEAN_SUBTREE_DIR" == "True" || "$CLEAN_SUBTREE_DIR" == "true" ]]; then
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/subtrees/subtree_*.nwk
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/core_subtrees/core_subtree_*.nwk
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/baseml_subtrees/baseml_subtree_*.nwk
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/subtree_summary.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/subtree_roots.txt
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/core_subtree_summary.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/core_subtree_roots.txt
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/baseml_subtree_summary.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/baseml_tree_manifest.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/overlap_report.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/split_tree.log
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/validation_report.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/baseml_validation_report.tsv
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/intermediate/rooted.tree
    rm -f "$PROJECT_ROOT/$OUTPUT_DIR"/intermediate/rooted.validation.tree
fi

SPLIT_ARGS=(
    --input-tree "$PROJECT_ROOT/$INPUT_TREE"
    --outgroup-tip-file "$PROJECT_ROOT/$OUTGROUP_TIP_FILE"
    --output-dir "$PROJECT_ROOT/$OUTPUT_DIR"
    --max-tips "$MAX_TIPS"
    --min-baseml-tips "$MIN_BASEML_TIPS"
    --threads "$THREADS"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --backbone-mode "$BACKBONE_MODE"
    --log-level "$LOG_LEVEL"
    --anchor-strategy "$ANCHOR_STRATEGY"
    --reserve-slots-for-outgroup "$RESERVE_SLOTS_FOR_OUTGROUP"
)

if [[ "$ENABLE_MERGE" == "True" || "$ENABLE_MERGE" == "true" ]]; then
    SPLIT_ARGS+=(--enable-merge-small-blocks)
fi

if [[ "$CLEAN_SUBTREE_DIR" == "True" || "$CLEAN_SUBTREE_DIR" == "true" ]]; then
    SPLIT_ARGS+=(--clean-subtree-dir)
fi

if [[ "$CONSTRUCT_BASEML" == "True" || "$CONSTRUCT_BASEML" == "true" ]]; then
    SPLIT_ARGS+=(--construct-baseml-subtrees)
fi

if [[ "$ALWAYS_INCLUDE_OUTGROUP" == "True" || "$ALWAYS_INCLUDE_OUTGROUP" == "true" ]]; then
    SPLIT_ARGS+=(--always-include-outgroup)
fi

if [[ "$BACKBONE_TREE" != "null" && -n "$BACKBONE_TREE" ]]; then
    SPLIT_ARGS+=(--backbone-tree "$PROJECT_ROOT/$BACKBONE_TREE")
fi

"$PYTHON_BIN" "$PROJECT_ROOT/python/split_phylo_tree.py" "${SPLIT_ARGS[@]}"

VALIDATE_ARGS=(
    --input-tree "$PROJECT_ROOT/$INPUT_TREE"
    --outgroup-tip-file "$PROJECT_ROOT/$OUTGROUP_TIP_FILE"
    --output-dir "$PROJECT_ROOT/$OUTPUT_DIR"
    --max-tips "$MAX_TIPS"
    --min-baseml-tips "$MIN_BASEML_TIPS"
    --threads "$THREADS"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --backbone-mode "$BACKBONE_MODE"
    --log-level "$LOG_LEVEL"
    --anchor-strategy "$ANCHOR_STRATEGY"
    --reserve-slots-for-outgroup "$RESERVE_SLOTS_FOR_OUTGROUP"
)

if [[ "$CONSTRUCT_BASEML" == "True" || "$CONSTRUCT_BASEML" == "true" ]]; then
    VALIDATE_ARGS+=(--construct-baseml-subtrees)
fi

if [[ "$ALWAYS_INCLUDE_OUTGROUP" == "True" || "$ALWAYS_INCLUDE_OUTGROUP" == "true" ]]; then
    VALIDATE_ARGS+=(--always-include-outgroup)
fi

if [[ "$ENABLE_MERGE" == "True" || "$ENABLE_MERGE" == "true" ]]; then
    VALIDATE_ARGS+=(--enable-merge-small-blocks)
fi

if [[ "$BACKBONE_TREE" != "null" && -n "$BACKBONE_TREE" ]]; then
    VALIDATE_ARGS+=(--backbone-tree "$PROJECT_ROOT/$BACKBONE_TREE")
fi

"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_phylo_split.py" "${VALIDATE_ARGS[@]}"

echo "[OK] Pipeline finished."

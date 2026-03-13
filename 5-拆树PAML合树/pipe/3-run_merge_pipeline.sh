#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
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
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

resolve_path() {
    local value="$1"
    if [[ "$value" == /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$PROJECT_ROOT/$value"
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

PYTHON_BIN=$(config_get tools.python)

bash "$PROJECT_ROOT/script/check_env.sh" \
    --python "$PYTHON_BIN"

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
PAML_OUTPUT_DIR=$(resolve_path "$(config_get paths.paml_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

MERGE_OUTPUT_DIR=$(resolve_path "$(config_get merge.output_dir)")
ANALYSIS_TREE_SOURCE=$(config_get merge.analysis_tree_source)
EXTERNAL_RESULT_DIR=$(config_get_optional merge.external_result_dir)
CLEAN_OUTPUT_DIR=$(config_get merge.clean_output_dir)
MERGE_MODE=$(config_get merge.mode)
BACKBONE_EDGE_AGGREGATION=$(config_get merge.backbone_edge_aggregation)
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
    "$SPLIT_OUTPUT_DIR/paml_subtree_summary.tsv" \
    "$SPLIT_OUTPUT_DIR/paml_tree_manifest.tsv"; do
    if [[ ! -f "$required" ]]; then
        echo "[ERROR] Required split output not found: $required" >&2
        exit 1
    fi
done

mkdir -p "$MERGE_OUTPUT_DIR" "$MERGE_OUTPUT_DIR/intermediate" "$MERGE_OUTPUT_DIR/simulated_baseml_subtrees"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$MERGE_OUTPUT_DIR"/simulation_manifest.tsv
    rm -f "$MERGE_OUTPUT_DIR"/assembly_scaffold.nwk
    rm -f "$MERGE_OUTPUT_DIR"/backbone_edge_estimates.tsv
    rm -f "$MERGE_OUTPUT_DIR"/graft_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/edge_update_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/merge_validation_report.tsv
    rm -f "$MERGE_OUTPUT_DIR"/merged_tree.nwk
    rm -f "$MERGE_OUTPUT_DIR"/merge.log
    rm -f "$MERGE_OUTPUT_DIR"/simulated_baseml_subtrees/*.nwk
fi

if [[ "$ANALYSIS_TREE_SOURCE" == "simulated" ]]; then
    SIMULATE_ARGS=(
        --paml-summary-tsv "$SPLIT_OUTPUT_DIR/paml_subtree_summary.tsv"
        --paml-tree-dir "$SPLIT_OUTPUT_DIR/paml_subtrees"
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
fi

MERGE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --merge-output-dir "$MERGE_OUTPUT_DIR"
    --analysis-tree-source "$ANALYSIS_TREE_SOURCE"
    --merge-mode "$MERGE_MODE"
    --backbone-edge-aggregation "$BACKBONE_EDGE_AGGREGATION"
    --min-branch-length "$MIN_BRANCH_LENGTH"
    --log-level "$LOG_LEVEL"
)

if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
    MERGE_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
fi

if [[ -n "$EXTERNAL_RESULT_DIR_ABS" ]]; then
    MERGE_ARGS+=(--external-result-dir "$EXTERNAL_RESULT_DIR_ABS")
fi

"$PYTHON_BIN" "$PROJECT_ROOT/python/merge_baseml_subtrees.py" "${MERGE_ARGS[@]}"

VALIDATE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --merge-output-dir "$MERGE_OUTPUT_DIR"
    --merge-mode "$MERGE_MODE"
    --backbone-edge-aggregation "$BACKBONE_EDGE_AGGREGATION"
    --log-level "$LOG_LEVEL"
)

if [[ -n "$OUTGROUP_TIP_NAME" ]]; then
    VALIDATE_ARGS+=(--outgroup-tip-name "$OUTGROUP_TIP_NAME")
fi

"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_merged_tree.py" "${VALIDATE_ARGS[@]}"

echo "[OK] Merge pipeline finished."

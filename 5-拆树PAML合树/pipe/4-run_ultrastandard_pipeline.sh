#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
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

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
PAML_OUTPUT_DIR=$(resolve_path "$(config_get paths.paml_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

ULTRA_OUTPUT_DIR=$(resolve_path "$(config_get ultrastandard.output_dir)")
CLEAN_OUTPUT_DIR=$(config_get ultrastandard.clean_output_dir)
MIN_BRANCH_LENGTH=$(config_get ultrastandard.min_branch_length)
ULTRAMETRIC_TOLERANCE=$(config_get ultrastandard.ultrametric_tolerance)
RELATIVE_ULTRA_METHOD=$(config_get ultrastandard.relative_ultrametric_method)
REQUIRE_BACKBONE_DESCENDANT_ANCHOR=$(config_get ultrastandard.require_backbone_descendant_anchor)
ALLOW_MULTIFURCATION=$(config_get ultrastandard.allow_multifurcation)
BACKBONE_ULTRA_TREE=$(config_get_optional ultrastandard.backbone_ultrametric_tree)

for required in \
    "$SPLIT_OUTPUT_DIR/intermediate/rooted.tree" \
    "$SPLIT_OUTPUT_DIR/backbone_summary.tsv" \
    "$SPLIT_OUTPUT_DIR/target_subtree_summary.tsv" \
    "$SPLIT_OUTPUT_DIR/paml_subtree_summary.tsv" \
    "$PAML_OUTPUT_DIR/backbone_analysis/backbone_ultrametric_tree.nwk"; do
    if [[ ! -f "$required" ]]; then
        echo "[ERROR] Required ultrastandard input not found: $required" >&2
        exit 1
    fi
done

if [[ ! -d "$PAML_OUTPUT_DIR/analysis_trees" ]]; then
    echo "[ERROR] Required analysis tree directory not found: $PAML_OUTPUT_DIR/analysis_trees" >&2
    exit 1
fi

mkdir -p "$ULTRA_OUTPUT_DIR" "$ULTRA_OUTPUT_DIR/relative_target_trees"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$ULTRA_OUTPUT_DIR"/assembly_scaffold.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/backbone_assigned_scaffold.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/target_scaling_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/graft_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/backbone_edge_assignment_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/root_to_tip_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/merged_tree.nwk
    rm -f "$ULTRA_OUTPUT_DIR"/ultrastandard_validation_report.tsv
    rm -f "$ULTRA_OUTPUT_DIR"/ultrastandard.log
    rm -rf "$ULTRA_OUTPUT_DIR/relative_target_trees"
    mkdir -p "$ULTRA_OUTPUT_DIR/relative_target_trees"
fi

MERGE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --paml-output-dir "$PAML_OUTPUT_DIR"
    --ultrastandard-output-dir "$ULTRA_OUTPUT_DIR"
    --min-branch-length "$MIN_BRANCH_LENGTH"
    --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE"
    --relative-ultrametric-method "$RELATIVE_ULTRA_METHOD"
    --require-backbone-descendant-anchor "$REQUIRE_BACKBONE_DESCENDANT_ANCHOR"
    --allow-multifurcation "$ALLOW_MULTIFURCATION"
    --log-level "$LOG_LEVEL"
)

if [[ "$BACKBONE_ULTRA_TREE" != "null" && -n "$BACKBONE_ULTRA_TREE" ]]; then
    MERGE_ARGS+=(--backbone-ultrametric-tree "$(resolve_path "$BACKBONE_ULTRA_TREE")")
fi

echo "[INFO] Stage 1/3: merge_ultrastandard_tree"
"$PYTHON_BIN" "$PROJECT_ROOT/python/merge_ultrastandard_tree.py" "${MERGE_ARGS[@]}"

VALIDATE_ARGS=(
    --split-output-dir "$SPLIT_OUTPUT_DIR"
    --paml-output-dir "$PAML_OUTPUT_DIR"
    --ultrastandard-output-dir "$ULTRA_OUTPUT_DIR"
    --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE"
    --log-level "$LOG_LEVEL"
)

if [[ "$BACKBONE_ULTRA_TREE" != "null" && -n "$BACKBONE_ULTRA_TREE" ]]; then
    VALIDATE_ARGS+=(--backbone-ultrametric-tree "$(resolve_path "$BACKBONE_ULTRA_TREE")")
fi

echo "[INFO] Stage 2/3: validate_ultrastandard_tree"
"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_ultrastandard_tree.py" "${VALIDATE_ARGS[@]}"

echo "[INFO] Stage 3/3: compare_final_topology"
if "$PYTHON_BIN" "$PROJECT_ROOT/python/check_tree_topology_match.py" \
    --reference-tree "$SPLIT_OUTPUT_DIR/intermediate/rooted.tree" \
    --query-tree "$ULTRA_OUTPUT_DIR/merged_tree.nwk" \
    --log-level "$LOG_LEVEL"; then
    echo "[OK] Ultrastandard merged tree topology matches the original rooted tree."
else
    compare_status=$?
    if [[ "$compare_status" -eq 2 ]]; then
        echo "[WARN] Ultrastandard merged tree topology does not match the original rooted tree."
    else
        exit "$compare_status"
    fi
fi

echo "[OK] Ultrastandard pipeline finished."

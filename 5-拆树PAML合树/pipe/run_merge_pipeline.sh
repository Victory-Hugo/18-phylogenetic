#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="$PROJECT_ROOT/conf/Config.json"
BOOTSTRAP_PYTHON="python3"

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

bash "$PROJECT_ROOT/script/check_env.sh"

CONFIG_LOADER="$PROJECT_ROOT/python/config_loader.py"
PYTHON_BIN=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key tools.python)
INPUT_TREE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.input_tree)
OUTGROUP_TIP_FILE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.outgroup_tip_file)
SPLIT_OUTPUT_DIR=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key paths.output_dir)
CONDA_ENV=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key tools.conda_env)
GOTREE_BIN=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key tools.gotree)
THREADS=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.threads)
LOG_LEVEL=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key runtime.log_level)

MERGE_ENABLED=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.enabled)
MERGE_OUTPUT_DIR=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.output_dir)
ANALYSIS_TREE_SOURCE=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.analysis_tree_source)
EXTERNAL_RESULT_DIR=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.external_result_dir)
CLEAN_OUTPUT_DIR=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.clean_output_dir)
RANDOMIZATION_MODEL=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.randomization_model)
RANDOMIZATION_SIGMA=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.randomization_sigma)
RANDOMIZATION_SEED=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.randomization_seed)
MIN_BRANCH_LENGTH=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.min_branch_length)
PRESERVE_MASTER_TOPOLOGY=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.preserve_master_topology)
STRICT_CORE_TOPOLOGY_MATCH=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.strict_core_topology_match)
SCAFFOLD_SCALE_AGGREGATION=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.scaffold_scale_aggregation)
SCAFFOLD_SCALE_CLAMP=$("$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key merge.scaffold_scale_clamp)

if [[ "$MERGE_ENABLED" != "True" && "$MERGE_ENABLED" != "true" ]]; then
    echo "[INFO] merge.enabled=false; skipping merge pipeline."
    exit 0
fi

INPUT_TREE_ABS=$(resolve_path "$INPUT_TREE")
OUTGROUP_TIP_FILE_ABS=$(resolve_path "$OUTGROUP_TIP_FILE")
SPLIT_OUTPUT_DIR_ABS=$(resolve_path "$SPLIT_OUTPUT_DIR")
MERGE_OUTPUT_DIR_ABS=$(resolve_path "$MERGE_OUTPUT_DIR")

if [[ "$EXTERNAL_RESULT_DIR" != "null" && -n "$EXTERNAL_RESULT_DIR" ]]; then
    EXTERNAL_RESULT_DIR_ABS=$(resolve_path "$EXTERNAL_RESULT_DIR")
else
    EXTERNAL_RESULT_DIR_ABS=""
fi

for required in \
    "$SPLIT_OUTPUT_DIR_ABS/intermediate/rooted.tree" \
    "$SPLIT_OUTPUT_DIR_ABS/core_subtree_summary.tsv" \
    "$SPLIT_OUTPUT_DIR_ABS/baseml_subtree_summary.tsv" \
    "$SPLIT_OUTPUT_DIR_ABS/baseml_tree_manifest.tsv"; do
    if [[ ! -f "$required" ]]; then
        echo "[ERROR] Required split output not found: $required" >&2
        exit 1
    fi
done

mkdir -p "$MERGE_OUTPUT_DIR_ABS" "$MERGE_OUTPUT_DIR_ABS/intermediate" "$MERGE_OUTPUT_DIR_ABS/simulated_baseml_subtrees"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$MERGE_OUTPUT_DIR_ABS"/simulation_manifest.tsv
    rm -f "$MERGE_OUTPUT_DIR_ABS"/merge_report.tsv
    rm -f "$MERGE_OUTPUT_DIR_ABS"/edge_update_report.tsv
    rm -f "$MERGE_OUTPUT_DIR_ABS"/merge_validation_report.tsv
    rm -f "$MERGE_OUTPUT_DIR_ABS"/merged_tree.nwk
    rm -f "$MERGE_OUTPUT_DIR_ABS"/merge.log
    rm -f "$MERGE_OUTPUT_DIR_ABS"/simulated_baseml_subtrees/*.nwk
fi

if [[ "$ANALYSIS_TREE_SOURCE" == "simulated" ]]; then
    "$PYTHON_BIN" "$PROJECT_ROOT/python/simulate_baseml_results.py" \
        --baseml-summary-tsv "$SPLIT_OUTPUT_DIR_ABS/baseml_subtree_summary.tsv" \
        --baseml-tree-dir "$SPLIT_OUTPUT_DIR_ABS/baseml_subtrees" \
        --output-dir "$MERGE_OUTPUT_DIR_ABS" \
        --randomization-model "$RANDOMIZATION_MODEL" \
        --randomization-sigma "$RANDOMIZATION_SIGMA" \
        --randomization-seed "$RANDOMIZATION_SEED" \
        --min-branch-length "$MIN_BRANCH_LENGTH" \
        --log-level "$LOG_LEVEL"
fi

MERGE_ARGS=(
    --input-tree "$INPUT_TREE_ABS"
    --outgroup-tip-file "$OUTGROUP_TIP_FILE_ABS"
    --split-output-dir "$SPLIT_OUTPUT_DIR_ABS"
    --merge-output-dir "$MERGE_OUTPUT_DIR_ABS"
    --analysis-tree-source "$ANALYSIS_TREE_SOURCE"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --threads "$THREADS"
    --scaffold-scale-aggregation "$SCAFFOLD_SCALE_AGGREGATION"
    --min-branch-length "$MIN_BRANCH_LENGTH"
    --log-level "$LOG_LEVEL"
)

if [[ "$PRESERVE_MASTER_TOPOLOGY" == "True" || "$PRESERVE_MASTER_TOPOLOGY" == "true" ]]; then
    MERGE_ARGS+=(--preserve-master-topology)
fi

if [[ "$STRICT_CORE_TOPOLOGY_MATCH" == "True" || "$STRICT_CORE_TOPOLOGY_MATCH" == "true" ]]; then
    MERGE_ARGS+=(--strict-core-topology-match)
fi

if [[ -n "$EXTERNAL_RESULT_DIR_ABS" ]]; then
    MERGE_ARGS+=(--external-result-dir "$EXTERNAL_RESULT_DIR_ABS")
fi

if [[ "$SCAFFOLD_SCALE_CLAMP" != "null" && -n "$SCAFFOLD_SCALE_CLAMP" ]]; then
    MERGE_ARGS+=(--scaffold-scale-clamp "$SCAFFOLD_SCALE_CLAMP")
fi

"$PYTHON_BIN" "$PROJECT_ROOT/python/merge_baseml_subtrees.py" "${MERGE_ARGS[@]}"

VALIDATE_ARGS=(
    --input-tree "$INPUT_TREE_ABS"
    --outgroup-tip-file "$OUTGROUP_TIP_FILE_ABS"
    --split-output-dir "$SPLIT_OUTPUT_DIR_ABS"
    --merge-output-dir "$MERGE_OUTPUT_DIR_ABS"
    --conda-env "$CONDA_ENV"
    --gotree-bin "$GOTREE_BIN"
    --threads "$THREADS"
    --scaffold-scale-aggregation "$SCAFFOLD_SCALE_AGGREGATION"
    --log-level "$LOG_LEVEL"
)

if [[ "$STRICT_CORE_TOPOLOGY_MATCH" == "True" || "$STRICT_CORE_TOPOLOGY_MATCH" == "true" ]]; then
    VALIDATE_ARGS+=(--strict-core-topology-match)
fi

"$PYTHON_BIN" "$PROJECT_ROOT/python/validate_merged_tree.py" "${VALIDATE_ARGS[@]}"

echo "[OK] Merge pipeline finished."

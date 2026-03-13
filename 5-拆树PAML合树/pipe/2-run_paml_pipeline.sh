#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/2-paml.yaml"
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

PYTHON_BIN=$(config_get tools.python)
BASEML_BIN=$(config_get tools.baseml)

bash "$PROJECT_ROOT/script/check_env.sh" \
    --python "$PYTHON_BIN" \
    --baseml "$BASEML_BIN"

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
INPUT_FASTA=$(resolve_path "$(config_get paths.input_fasta)")
TMP_DIR=$(resolve_path "$(config_get paths.tmp_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

CTL_TEMPLATE=$(resolve_path "$(config_get paml.ctl_template)")
PAML_OUTPUT_DIR=$(resolve_path "$(config_get paml.output_dir)")
PAML_CLEAN_OUTPUT_DIR=$(config_get paml.clean_output_dir)
PAML_PARALLEL_JOBS=$(config_get paml.parallel_jobs)
SEQ_ID_STRATEGY=$(config_get paml.seq_id_strategy)
PAML_SKIP_EXISTING=$(config_get paml.skip_existing)
BACKBONE_ONLY_ENABLED=$(config_get paml.backbone_only_enabled)
BACKBONE_JOB_DIR=$(resolve_path "$(config_get paml.backbone_job_dir)")
BACKBONE_ANALYSIS_DIR=$(resolve_path "$(config_get paml.backbone_analysis_dir)")
BACKBONE_SKIP_EXISTING=$(config_get paml.backbone_skip_existing)
ULTRAMETRIC_NORMALIZATION=$(config_get paml.ultrametric_normalization)
ULTRAMETRIC_TOLERANCE=$(config_get paml.ultrametric_tolerance)

PAML_SUBTREE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/paml_subtree_summary.tsv"
PAML_TREE_DIR="$SPLIT_OUTPUT_DIR/paml_subtrees"
MANIFEST_PATH="$PAML_OUTPUT_DIR/paml_job_manifest.tsv"
BACKBONE_TREE="$SPLIT_OUTPUT_DIR/backbone_tree.nwk"
BACKBONE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/backbone_summary.tsv"

for required in \
    "$PAML_SUBTREE_SUMMARY_TSV" \
    "$INPUT_FASTA" \
    "$CTL_TEMPLATE" \
    "$BACKBONE_TREE" \
    "$BACKBONE_SUMMARY_TSV"; do
    if [[ ! -f "$required" ]]; then
        echo "[ERROR] Required paml input not found: $required" >&2
        exit 1
    fi
done

if [[ ! -d "$PAML_TREE_DIR" ]]; then
    echo "[ERROR] PAML subtree directory not found: $PAML_TREE_DIR" >&2
    exit 1
fi

mkdir -p "$PAML_OUTPUT_DIR" "$TMP_DIR"

if [[ "$PAML_CLEAN_OUTPUT_DIR" == "True" || "$PAML_CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -rf "$PAML_OUTPUT_DIR/jobs"
    rm -rf "$PAML_OUTPUT_DIR/analysis_trees"
    rm -rf "$PAML_OUTPUT_DIR/tables"
    rm -rf "$BACKBONE_JOB_DIR"
    rm -rf "$BACKBONE_ANALYSIS_DIR"
    rm -f "$PAML_OUTPUT_DIR/paml_job_manifest.tsv"
    rm -f "$PAML_OUTPUT_DIR/paml_run_status.tsv"
    rm -f "$PAML_OUTPUT_DIR/paml_parse_status.tsv"
    rm -f "$PAML_OUTPUT_DIR/prepare_paml_inputs.log"
    rm -f "$PAML_OUTPUT_DIR/run_paml_jobs.log"
    rm -f "$PAML_OUTPUT_DIR/parse_paml_outputs.log"
    rm -rf "$TMP_DIR/paml_parse"
fi

echo "[INFO] Stage 1/6: prepare_paml_inputs"
"$PYTHON_BIN" "$PROJECT_ROOT/python/prepare_paml_inputs.py" \
    --paml-summary-tsv "$PAML_SUBTREE_SUMMARY_TSV" \
    --paml-tree-dir "$PAML_TREE_DIR" \
    --input-fasta "$INPUT_FASTA" \
    --ctl-template "$CTL_TEMPLATE" \
    --output-dir "$PAML_OUTPUT_DIR" \
    --seq-id-strategy "$SEQ_ID_STRATEGY" \
    --log-level "$LOG_LEVEL"

RUN_ARGS=(
    --manifest "$MANIFEST_PATH"
    --baseml-bin "$BASEML_BIN"
    --parallel-jobs "$PAML_PARALLEL_JOBS"
    --log-level "$LOG_LEVEL"
)

if [[ "$PAML_SKIP_EXISTING" == "True" || "$PAML_SKIP_EXISTING" == "true" ]]; then
    RUN_ARGS+=(--skip-existing)
fi

echo "[INFO] Stage 2/6: run_paml_jobs"
"$PYTHON_BIN" "$PROJECT_ROOT/python/run_paml_jobs.py" "${RUN_ARGS[@]}"

echo "[INFO] Stage 3/6: parse_paml_outputs"
"$PYTHON_BIN" "$PROJECT_ROOT/python/parse_paml_outputs.py" \
    --paml-summary-tsv "$PAML_SUBTREE_SUMMARY_TSV" \
    --manifest "$MANIFEST_PATH" \
    --output-dir "$PAML_OUTPUT_DIR" \
    --tmp-dir "$TMP_DIR" \
    --log-level "$LOG_LEVEL"

if [[ "$BACKBONE_ONLY_ENABLED" == "True" || "$BACKBONE_ONLY_ENABLED" == "true" ]]; then
    if [[ "$ULTRAMETRIC_NORMALIZATION" != "extend_terminal_to_max_depth" ]]; then
        echo "[ERROR] Unsupported paml.ultrametric_normalization: $ULTRAMETRIC_NORMALIZATION" >&2
        exit 1
    fi

    echo "[INFO] Stage 4/6: prepare_backbone_paml_input"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/prepare_backbone_paml_input.py" \
        --backbone-tree "$BACKBONE_TREE" \
        --backbone-summary-tsv "$BACKBONE_SUMMARY_TSV" \
        --paml-summary-tsv "$PAML_SUBTREE_SUMMARY_TSV" \
        --input-fasta "$INPUT_FASTA" \
        --ctl-template "$CTL_TEMPLATE" \
        --job-dir "$BACKBONE_JOB_DIR" \
        --analysis-dir "$BACKBONE_ANALYSIS_DIR" \
        --seq-id-strategy "$SEQ_ID_STRATEGY" \
        --log-level "$LOG_LEVEL"

    BACKBONE_RUN_ARGS=(
        --job-dir "$BACKBONE_JOB_DIR"
        --baseml-bin "$BASEML_BIN"
        --log-level "$LOG_LEVEL"
    )
    if [[ "$BACKBONE_SKIP_EXISTING" == "True" || "$BACKBONE_SKIP_EXISTING" == "true" ]]; then
        BACKBONE_RUN_ARGS+=(--skip-existing)
    fi

    echo "[INFO] Stage 5/6: run_backbone_paml_job"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/run_backbone_paml_job.py" "${BACKBONE_RUN_ARGS[@]}"

    echo "[INFO] Stage 6/6: parse_backbone_paml_output"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/parse_backbone_paml_output.py" \
        --backbone-tree "$BACKBONE_TREE" \
        --paml-summary-tsv "$PAML_SUBTREE_SUMMARY_TSV" \
        --job-dir "$BACKBONE_JOB_DIR" \
        --analysis-dir "$BACKBONE_ANALYSIS_DIR" \
        --min-branch-length 1e-8 \
        --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE" \
        --log-level "$LOG_LEVEL"
fi

echo "[OK] PAML pipeline finished."
echo "[INFO] Cleaning up temporary files..."
rm -f ./rst ./rst1 ./rub

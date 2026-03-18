#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PATH_ROOT="$PROJECT_ROOT"
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

config_get_default() {
    "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1" --default "$2"
}

CONFIG_PROJECT_ROOT=$(config_get_default projectpath "$PROJECT_ROOT")
PATH_ROOT=$(resolve_path "$CONFIG_PROJECT_ROOT")

PYTHON_BIN=$(resolve_command_or_path "$(config_get tools.python)")
BASEML_BIN=$(resolve_command_or_path "$(config_get tools.baseml)")
PAML_EXECUTION_MODE=$(config_get_default paml.execution_mode real)

if [[ "$PAML_EXECUTION_MODE" == "fake" ]]; then
    bash "$PROJECT_ROOT/script/check_env.sh" --python "$PYTHON_BIN"
else
    bash "$PROJECT_ROOT/script/check_env.sh" \
        --python "$PYTHON_BIN" \
        --baseml "$BASEML_BIN"
fi

SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
INPUT_FASTA=$(resolve_path "$(config_get paths.input_fasta)")
TMP_DIR=$(resolve_path "$(config_get paths.tmp_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

CTL_TEMPLATE=$(resolve_path "$(config_get paml.ctl_template)")
PAML_OUTPUT_DIR=$(resolve_path "$(config_get paml.output_dir)")
PAML_CLEAN_OUTPUT_DIR=$(config_get paml.clean_output_dir)
PAML_PARALLEL_JOBS=$(config_get paml.parallel_jobs)
PAML_THREADS_PER_JOB=$(config_get paml.threads_per_job)
PAML_ENABLE_PAR3=$(config_get paml.enable_par3)
PAML_BIND_CPU_AFFINITY=$(config_get_default paml.bind_cpu_affinity true)
SEQ_ID_STRATEGY=$(config_get paml.seq_id_strategy)
PAML_SKIP_EXISTING=$(config_get paml.skip_existing)
BACKBONE_ONLY_ENABLED=$(config_get paml.backbone_only_enabled)
BACKBONE_JOB_DIR=$(resolve_path "$(config_get paml.backbone_job_dir)")
BACKBONE_ANALYSIS_DIR=$(resolve_path "$(config_get paml.backbone_analysis_dir)")
BACKBONE_SKIP_EXISTING=$(config_get paml.backbone_skip_existing)
ULTRAMETRIC_NORMALIZATION=$(config_get paml.ultrametric_normalization)
ULTRAMETRIC_TOLERANCE=$(config_get paml.ultrametric_tolerance)
FAKE_BRANCH_LENGTH_MODEL=$(config_get_default paml.fake_branch_length_model lognormal_then_extend_terminal)
FAKE_BRANCH_LENGTH_SIGMA=$(config_get_default paml.fake_branch_length_sigma 0.25)
FAKE_RANDOM_SEED=$(config_get_default paml.fake_random_seed 42)
FAKE_MIN_BRANCH_LENGTH=$(config_get_default paml.fake_min_branch_length 1e-8)

PAML_SUBTREE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/paml_subtree_summary.tsv"
PAML_TREE_DIR="$SPLIT_OUTPUT_DIR/paml_subtrees"
MANIFEST_PATH="$PAML_OUTPUT_DIR/paml_job_manifest.tsv"
PAML_RESUME_LOG_DIR="$PAML_OUTPUT_DIR/log"
PAML_RESUME_SUCCESS_LOG="$PAML_RESUME_LOG_DIR/success.log"
BACKBONE_TREE="$SPLIT_OUTPUT_DIR/backbone_tree.nwk"
BACKBONE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/backbone_summary.tsv"

if ! [[ "$PAML_THREADS_PER_JOB" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] paml.threads_per_job must be a positive integer: $PAML_THREADS_PER_JOB" >&2
    exit 1
fi

PAML_RUNTIME_ENV=(
    "OMP_NUM_THREADS=$PAML_THREADS_PER_JOB"
)

if [[ "$PAML_ENABLE_PAR3" == "True" || "$PAML_ENABLE_PAR3" == "true" ]]; then
    PAML_RUNTIME_ENV+=("BASEML_PAR3=1")
fi

for required in \
    "$PAML_SUBTREE_SUMMARY_TSV" \
    "$BACKBONE_TREE" \
    "$BACKBONE_SUMMARY_TSV"; do
    if [[ ! -f "$required" ]]; then
        echo "[ERROR] Required paml input not found: $required" >&2
        exit 1
    fi
done

if [[ "$PAML_EXECUTION_MODE" != "fake" ]]; then
    for required in "$INPUT_FASTA" "$CTL_TEMPLATE"; do
        if [[ ! -f "$required" ]]; then
            echo "[ERROR] Required paml input not found: $required" >&2
            exit 1
        fi
    done
fi

if [[ ! -d "$PAML_TREE_DIR" ]]; then
    echo "[ERROR] PAML subtree directory not found: $PAML_TREE_DIR" >&2
    exit 1
fi

mkdir -p "$PAML_OUTPUT_DIR" "$TMP_DIR"

if [[ "$PAML_CLEAN_OUTPUT_DIR" == "True" || "$PAML_CLEAN_OUTPUT_DIR" == "true" ]]; then
    if [[ -s "$PAML_RESUME_SUCCESS_LOG" ]]; then
        echo "[INFO] Resume logs detected at $PAML_RESUME_SUCCESS_LOG; skip destructive cleanup even though paml.clean_output_dir=true."
    else
        rm -rf "$PAML_OUTPUT_DIR/jobs"
        rm -rf "$PAML_OUTPUT_DIR/analysis_trees"
        rm -rf "$PAML_OUTPUT_DIR/tables"
        rm -rf "$PAML_OUTPUT_DIR/log"
        rm -rf "$BACKBONE_JOB_DIR"
        rm -rf "$BACKBONE_ANALYSIS_DIR"
        rm -f "$PAML_OUTPUT_DIR/paml_job_manifest.tsv"
        rm -f "$PAML_OUTPUT_DIR/paml_run_status.tsv"
        rm -f "$PAML_OUTPUT_DIR/paml_parse_status.tsv"
        rm -f "$PAML_OUTPUT_DIR/fake_paml_manifest.tsv"
        rm -f "$PAML_OUTPUT_DIR/prepare_paml_inputs.log"
        rm -f "$PAML_OUTPUT_DIR/run_paml_jobs.log"
        rm -f "$PAML_OUTPUT_DIR/parse_paml_outputs.log"
        rm -f "$PAML_OUTPUT_DIR/fake_paml_outputs.log"
        rm -rf "$TMP_DIR/paml_parse"
    fi
fi

if [[ "$PAML_EXECUTION_MODE" == "fake" ]]; then
    echo "[INFO] Stage 1/1: fake_paml_outputs"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/fake_paml_outputs.py" \
        --paml-summary-tsv "$PAML_SUBTREE_SUMMARY_TSV" \
        --paml-tree-dir "$PAML_TREE_DIR" \
        --backbone-tree "$BACKBONE_TREE" \
        --output-dir "$PAML_OUTPUT_DIR" \
        --backbone-analysis-dir "$BACKBONE_ANALYSIS_DIR" \
        --branch-length-model "$FAKE_BRANCH_LENGTH_MODEL" \
        --branch-length-sigma "$FAKE_BRANCH_LENGTH_SIGMA" \
        --random-seed "$FAKE_RANDOM_SEED" \
        --min-branch-length "$FAKE_MIN_BRANCH_LENGTH" \
        --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE" \
        --backbone-only-enabled "$BACKBONE_ONLY_ENABLED" \
        --log-level "$LOG_LEVEL"
    echo "[OK] PAML pipeline finished in fake mode."
    exit 0
fi

RUN_ARGS=(
    --manifest "$MANIFEST_PATH"
    --baseml-bin "$BASEML_BIN"
    --parallel-jobs "$PAML_PARALLEL_JOBS"
    --threads-per-job "$PAML_THREADS_PER_JOB"
    --resume-log-dir "$PAML_RESUME_LOG_DIR"
    --log-level "$LOG_LEVEL"
)

if [[ "$PAML_SKIP_EXISTING" == "True" || "$PAML_SKIP_EXISTING" == "true" ]]; then
    RUN_ARGS+=(--legacy-skip-existing)
fi
if [[ "$PAML_BIND_CPU_AFFINITY" == "True" || "$PAML_BIND_CPU_AFFINITY" == "true" ]]; then
    RUN_ARGS+=(--bind-cpu-affinity)
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

echo "[INFO] Stage 2/6: run_paml_jobs"
env "${PAML_RUNTIME_ENV[@]}" \
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

    BACKBONE_RUN_ARGS=(
        --job-dir "$BACKBONE_JOB_DIR"
        --baseml-bin "$BASEML_BIN"
        --threads-per-job "$PAML_THREADS_PER_JOB"
        --log-level "$LOG_LEVEL"
    )
    if [[ "$BACKBONE_SKIP_EXISTING" == "True" || "$BACKBONE_SKIP_EXISTING" == "true" ]]; then
        BACKBONE_RUN_ARGS+=(--skip-existing)
    fi
    if [[ "$PAML_BIND_CPU_AFFINITY" == "True" || "$PAML_BIND_CPU_AFFINITY" == "true" ]]; then
        BACKBONE_RUN_ARGS+=(--bind-cpu-affinity)
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

    echo "[INFO] Stage 5/6: run_backbone_paml_job"
    env "${PAML_RUNTIME_ENV[@]}" \
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

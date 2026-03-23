#!/usr/bin/env bash
set -euo pipefail

# 中文说明：IQ2MC 主控脚本，负责配置读取、环境检查、流程调度与结果汇总。
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="${1:-$PROJECT_ROOT/conf/Config.yaml}"

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "[ERROR] Config file not found: $CONFIG_PATH" >&2
    exit 1
fi

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

require_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        echo "[ERROR] Missing required config value: $var_name" >&2
        exit 1
    fi
}

to_bool() {
    local raw="${1:-}"
    raw="${raw,,}"
    case "$raw" in
        true|1|yes|y) echo "true" ;;
        false|0|no|n|"") echo "false" ;;
        *)
            echo "[ERROR] Invalid bool value: $1" >&2
            exit 1
            ;;
    esac
}

resolve_path() {
    local value="$1"
    if [[ "$value" = /* ]]; then
        echo "$value"
    else
        echo "$PROJECT_ROOT/$value"
    fi
}

require_var CFG_PATHS_MSA_FILE
require_var CFG_PATHS_ROOTED_CALIBRATED_TREE
require_var CFG_PATHS_OUTPUT_DIR
require_var CFG_TOOLS_CONDA_BIN
require_var CFG_TOOLS_IQTREE3_BIN
require_var CFG_TOOLS_MCMCTREE_BIN
require_var CFG_RUNTIME_THREADS
require_var CFG_RUNTIME_CREATE_ENV_IF_MISSING
require_var CFG_RUNTIME_RUN_STEP1
require_var CFG_RUNTIME_RUN_STEP2
require_var CFG_RUNTIME_RUN_STEP3
require_var CFG_MODELING_SUBSTITUTION_MODEL
require_var CFG_MODELING_CLOCK_MODEL
require_var CFG_MODELING_MCMC_ITER
require_var CFG_MODELING_MCMC_BDS
require_var CFG_MCMCTREE_CHECKPOINT
require_var CFG_MCMCTREE_RESUME

CONDA_BIN=$(resolve_path "$CFG_TOOLS_CONDA_BIN")
IQTREE3_BIN=$(resolve_path "$CFG_TOOLS_IQTREE3_BIN")
MCMCTREE_BIN=$(resolve_path "$CFG_TOOLS_MCMCTREE_BIN")
MSA_FILE="$CFG_PATHS_MSA_FILE"
ROOTED_TREE="$CFG_PATHS_ROOTED_CALIBRATED_TREE"
PARTITION_FILE="${CFG_PATHS_PARTITION_FILE:-}"
OUTPUT_DIR_REL="$CFG_PATHS_OUTPUT_DIR"
THREADS="$CFG_RUNTIME_THREADS"
CREATE_ENV_IF_MISSING=$(to_bool "$CFG_RUNTIME_CREATE_ENV_IF_MISSING")
RUN_STEP1=$(to_bool "$CFG_RUNTIME_RUN_STEP1")
RUN_STEP2=$(to_bool "$CFG_RUNTIME_RUN_STEP2")
RUN_STEP3=$(to_bool "$CFG_RUNTIME_RUN_STEP3")
SUBSTITUTION_MODEL="$CFG_MODELING_SUBSTITUTION_MODEL"
PARTITION_MODE="${CFG_MODELING_PARTITION_MODE:-Q}"
CLOCK_MODEL="$CFG_MODELING_CLOCK_MODEL"
MCMC_ITER="$CFG_MODELING_MCMC_ITER"
MCMC_BDS="$CFG_MODELING_MCMC_BDS"
CHECKPOINT="$CFG_MCMCTREE_CHECKPOINT"
RESUME=$(to_bool "$CFG_MCMCTREE_RESUME")
ENV_NAME="IQ2MC"

OUTPUT_DIR=$(resolve_path "$OUTPUT_DIR_REL")
RUN_LOG_DIR="$OUTPUT_DIR/logs"
STEP1_DIR="$OUTPUT_DIR/step1"
STEP2_DIR="$OUTPUT_DIR/step2"
STEP3_DIR="$OUTPUT_DIR/step3"
mkdir -p "$RUN_LOG_DIR" "$STEP1_DIR" "$STEP2_DIR" "$STEP3_DIR"

COMMANDS_LOG="$RUN_LOG_DIR/commands.log"
: > "$COMMANDS_LOG"

CONDA_PYTHON=("$CONDA_BIN" run --no-capture-output -n "$ENV_NAME" python)

"$PROJECT_ROOT/script/check_env.sh" "$CONDA_BIN" "$IQTREE3_BIN" "$MCMCTREE_BIN"

"$PROJECT_ROOT/script/setup_conda_env.sh" \
    --conda-bin "$CONDA_BIN" \
    --env-name "$ENV_NAME" \
    --create-if-missing "$CREATE_ENV_IF_MISSING"

# 中文说明：先做输入校验并标准化路径，避免后续步骤运行时才报错。
"${CONDA_PYTHON[@]}" "$PROJECT_ROOT/python/validate_inputs.py" \
    --project-root "$PROJECT_ROOT" \
    --msa-file "$MSA_FILE" \
    --rooted-calibrated-tree "$ROOTED_TREE" \
    --partition-file "$PARTITION_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --json-out "$RUN_LOG_DIR/validated_inputs.json"

STEPS_EXECUTED=()
KEY_OUTPUTS=()

if [[ "$RUN_STEP1" == "true" ]]; then
    STEP1_PREFIX="$STEP1_DIR/step1"
    STEP1_ARGS=(
        --iqtree3-bin "$IQTREE3_BIN"
        --msa-file "$(resolve_path "$MSA_FILE")"
        --output-prefix "$STEP1_PREFIX"
        --threads "$THREADS"
        --substitution-model "$SUBSTITUTION_MODEL"
        --partition-mode "$PARTITION_MODE"
        --shell-only
    )
    if [[ -n "$PARTITION_FILE" ]]; then
        STEP1_ARGS+=(--partition-file "$(resolve_path "$PARTITION_FILE")")
    fi

    STEP1_CMD=$("${CONDA_PYTHON[@]}" "$PROJECT_ROOT/python/build_iqtree_cmd.py" \
        "${STEP1_ARGS[@]}")

    echo "$STEP1_CMD" | tee -a "$COMMANDS_LOG"
    eval "$STEP1_CMD" | tee "$RUN_LOG_DIR/step1.log"
    STEPS_EXECUTED+=("step1")
    KEY_OUTPUTS+=("$STEP1_PREFIX.treefile")
fi

if [[ "$RUN_STEP2" == "true" ]]; then
    STEP2_PREFIX="$STEP2_DIR/iq2mc"
    ROOTED_TREE_ABS=$(resolve_path "$ROOTED_TREE")

    STEP2_ARGS=(
        --iqtree3-bin "$IQTREE3_BIN"
        --msa-file "$(resolve_path "$MSA_FILE")"
        --output-prefix "$STEP2_PREFIX"
        --threads "$THREADS"
        --substitution-model "$SUBSTITUTION_MODEL"
        --rooted-tree "$ROOTED_TREE_ABS"
        --partition-mode "$PARTITION_MODE"
        --dating-mode mcmctree
        --mcmc-iter "$MCMC_ITER"
        --mcmc-bds "$MCMC_BDS"
        --mcmc-clock "$CLOCK_MODEL"
        --redo
        --shell-only
    )
    if [[ -n "$PARTITION_FILE" ]]; then
        STEP2_ARGS+=(--partition-file "$(resolve_path "$PARTITION_FILE")")
    fi

    STEP2_CMD=$("${CONDA_PYTHON[@]}" "$PROJECT_ROOT/python/build_iqtree_cmd.py" \
        "${STEP2_ARGS[@]}")

    echo "$STEP2_CMD" | tee -a "$COMMANDS_LOG"
    eval "$STEP2_CMD" | tee "$RUN_LOG_DIR/step2.log"

    HESSIAN_FILE="$STEP2_PREFIX.mcmctree.hessian"
    CTL_FILE="$STEP2_PREFIX.mcmctree.ctl"
    ROOTED_OUT="$STEP2_PREFIX.rooted.nwk"
    DUMMY_ALN=""
    if [[ -f "$STEP2_PREFIX.dummy.aln" ]]; then
        DUMMY_ALN="$STEP2_PREFIX.dummy.aln"
    elif [[ -f "$STEP2_PREFIX.dummy.phy" ]]; then
        DUMMY_ALN="$STEP2_PREFIX.dummy.phy"
    fi

    [[ -f "$HESSIAN_FILE" ]] || { echo "[ERROR] Missing hessian file: $HESSIAN_FILE" >&2; exit 1; }
    [[ -f "$CTL_FILE" ]] || { echo "[ERROR] Missing ctl file: $CTL_FILE" >&2; exit 1; }
    [[ -f "$ROOTED_OUT" ]] || { echo "[ERROR] Missing rooted tree output: $ROOTED_OUT" >&2; exit 1; }
    [[ -n "$DUMMY_ALN" ]] || { echo "[ERROR] Missing dummy alignment: $STEP2_PREFIX.dummy.aln or $STEP2_PREFIX.dummy.phy" >&2; exit 1; }

    STEPS_EXECUTED+=("step2")
    KEY_OUTPUTS+=("$HESSIAN_FILE" "$CTL_FILE" "$ROOTED_OUT" "$DUMMY_ALN")
fi

if [[ "$RUN_STEP3" == "true" ]]; then
    [[ -n "${CTL_FILE:-}" ]] || { echo "[ERROR] Step3 requires Step2 output ctl file." >&2; exit 1; }

    FINAL_CTL="$STEP3_DIR/mcmctree.run.ctl"
    MCMC_LOG="$STEP3_DIR/mcmc.txt"
    MCMC_OUT="$STEP3_DIR/mcmctree.out"
    MCMC_CKP="$STEP3_DIR/mcmctree.ckp"

    PREP_ARGS=(
        --input-ctl "$CTL_FILE"
        --output-ctl "$FINAL_CTL"
        --seqfile "$DUMMY_ALN"
        --treefile "$ROOTED_OUT"
        --mcmcfile "$MCMC_LOG"
        --outfile "$MCMC_OUT"
        --ckpfile "$MCMC_CKP"
        --hessianfile "$HESSIAN_FILE"
        --clock-model "$CLOCK_MODEL"
        --mcmc-iter "$MCMC_ITER"
        --mcmc-bds "$MCMC_BDS"
        --checkpoint "$CHECKPOINT"
    )

    if [[ "$RESUME" == "true" ]]; then
        PREP_ARGS+=(--resume)
    fi

    "${CONDA_PYTHON[@]}" "$PROJECT_ROOT/python/prepare_mcmctree_ctl.py" "${PREP_ARGS[@]}"

    STEP3_CMD="$MCMCTREE_BIN $FINAL_CTL"
    echo "$STEP3_CMD" | tee -a "$COMMANDS_LOG"
    eval "$STEP3_CMD" | tee "$RUN_LOG_DIR/step3.log"

    [[ -f "$MCMC_OUT" ]] || { echo "[ERROR] Missing MCMC output: $MCMC_OUT" >&2; exit 1; }

    STEPS_EXECUTED+=("step3")
    KEY_OUTPUTS+=("$FINAL_CTL" "$MCMC_LOG" "$MCMC_OUT" "$MCMC_CKP")
fi

STEPS_CSV=$(IFS=, ; echo "${STEPS_EXECUTED[*]}")
OUTPUTS_CSV=$(IFS=, ; echo "${KEY_OUTPUTS[*]}")

"${CONDA_PYTHON[@]}" "$PROJECT_ROOT/python/collect_summary.py" \
    --summary-file "$RUN_LOG_DIR/summary.md" \
    --run-name "$(basename "$OUTPUT_DIR")" \
    --config-file "$CONFIG_PATH" \
    --steps-executed "$STEPS_CSV" \
    --commands-file "$COMMANDS_LOG" \
    --key-outputs "$OUTPUTS_CSV" \
    --exit-status 0

echo "[OK] IQ2MC pipeline finished. Summary: $RUN_LOG_DIR/summary.md"

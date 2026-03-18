#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/5-time_calib.yaml"
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

ULTRA_OUTPUT_DIR=$(resolve_path "$(config_get paths.ultrastandard_output_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

TIME_CALIB_INPUT_TREE=$(config_get_optional time_calib.input_tree)
TIME_CALIB_OUTPUT_DIR=$(resolve_path "$(config_get time_calib.output_dir)")
CLEAN_OUTPUT_DIR=$(config_get time_calib.clean_output_dir)
TIME_CALIB_METHOD=$(config_get_optional time_calib.method molecular_clock)
SUBSTITUTION_RATE=$(config_get time_calib.substitution_rate_per_site_per_year)
SEQUENCE_LENGTH=$(config_get time_calib.sequence_length)
BRANCH_LENGTH_UNIT=$(config_get_optional time_calib.branch_length_unit substitutions_per_site)
NODE_CALIBRATION_TIP_NAME=$(config_get_optional time_calib.node_calibration_tip_name RSRS)
NODE_CALIBRATION_DIVERGENCE_YEARS=$(config_get_optional time_calib.node_calibration_divergence_years 180000)
OUTPUT_TREE_NAME=$(config_get_optional time_calib.output_tree_name merged_ultrametric_tree_years.nwk)

if [[ "$TIME_CALIB_INPUT_TREE" == "null" || -z "$TIME_CALIB_INPUT_TREE" ]]; then
    INPUT_TREE="$ULTRA_OUTPUT_DIR/merged_ultrametric_tree.nwk"
else
    INPUT_TREE=$(resolve_path "$TIME_CALIB_INPUT_TREE")
fi

if [[ ! -f "$INPUT_TREE" ]]; then
    echo "[ERROR] Required time calibration input not found: $INPUT_TREE" >&2
    exit 1
fi

mkdir -p "$TIME_CALIB_OUTPUT_DIR"

if [[ "$CLEAN_OUTPUT_DIR" == "True" || "$CLEAN_OUTPUT_DIR" == "true" ]]; then
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration.log
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration_summary.tsv
    rm -f "$TIME_CALIB_OUTPUT_DIR"/time_calibration_edge_report.tsv
    rm -f "$TIME_CALIB_OUTPUT_DIR"/"$OUTPUT_TREE_NAME"
fi

echo "[INFO] Stage 1/1: time_calibrate_tree"
"$PYTHON_BIN" "$PROJECT_ROOT/python/time_calibrate_tree.py" \
    --input-tree "$INPUT_TREE" \
    --output-dir "$TIME_CALIB_OUTPUT_DIR" \
    --output-tree-name "$OUTPUT_TREE_NAME" \
    --method "$TIME_CALIB_METHOD" \
    --substitution-rate-per-site-per-year "$SUBSTITUTION_RATE" \
    --sequence-length "$SEQUENCE_LENGTH" \
    --branch-length-unit "$BRANCH_LENGTH_UNIT" \
    --node-calibration-tip-name "$NODE_CALIBRATION_TIP_NAME" \
    --node-calibration-divergence-years "$NODE_CALIBRATION_DIVERGENCE_YEARS" \
    --log-level "$LOG_LEVEL"

echo "[OK] Time calibration pipeline finished."

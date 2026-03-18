#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

SPLIT_CONFIG="$PROJECT_ROOT/conf/1-split.yaml"
PAML_CONFIG="$PROJECT_ROOT/conf/2-paml.yaml"
MERGE_CONFIG="$PROJECT_ROOT/conf/3-merge.yaml"
ULTRASTANDARD_CONFIG="$PROJECT_ROOT/conf/4-ultrastandard.yaml"
TIME_CALIB_CONFIG="$PROJECT_ROOT/conf/5-time_calib.yaml"

usage() {
    cat <<EOF
Usage:
  bash pipe/0-run_all_pipeline.sh [options]

Options:
  --split-config PATH          Config for step 1 split
  --paml-config PATH           Config for step 2 paml
  --merge-config PATH          Config for step 3 merge
  --ultrastandard-config PATH  Config for step 4 ultrastandard
  --time-calib-config PATH     Config for step 5 time calibration
  --help                       Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --split-config)
            SPLIT_CONFIG="$2"
            shift 2
            ;;
        --paml-config)
            PAML_CONFIG="$2"
            shift 2
            ;;
        --merge-config)
            MERGE_CONFIG="$2"
            shift 2
            ;;
        --ultrastandard-config)
            ULTRASTANDARD_CONFIG="$2"
            shift 2
            ;;
        --time-calib-config)
            TIME_CALIB_CONFIG="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

for config_path in \
    "$SPLIT_CONFIG" \
    "$PAML_CONFIG" \
    "$MERGE_CONFIG" \
    "$ULTRASTANDARD_CONFIG" \
    "$TIME_CALIB_CONFIG"; do
    if [[ ! -f "$config_path" ]]; then
        echo "[ERROR] Config file not found: $config_path" >&2
        exit 1
    fi
done

run_stage() {
    local stage_label="$1"
    shift
    echo "[INFO] ===== $stage_label ====="
    "$@"
    echo "[OK] ===== $stage_label ====="
}

run_stage "Step 1/5: split" \
    bash "$PROJECT_ROOT/pipe/1-run_split_pipeline.sh" \
    --config "$SPLIT_CONFIG"

run_stage "Step 2/5: paml" \
    bash "$PROJECT_ROOT/pipe/2-run_paml_pipeline.sh" \
    --config "$PAML_CONFIG"

run_stage "Step 3/5: merge" \
    bash "$PROJECT_ROOT/pipe/3-run_merge_pipeline.sh" \
    --config "$MERGE_CONFIG"

run_stage "Step 4/5: ultrastandard" \
    bash "$PROJECT_ROOT/pipe/4-run_ultrastandard_pipeline.sh" \
    --config "$ULTRASTANDARD_CONFIG"

run_stage "Step 5/5: time calibration" \
    bash "$PROJECT_ROOT/pipe/5-run_time_calib_pipeline.sh" \
    --config "$TIME_CALIB_CONFIG"

echo "[OK] Full pipeline finished."

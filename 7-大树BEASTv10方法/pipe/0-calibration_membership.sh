#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

CONFIG="conf/0-calibration_membership.yaml"
source script/load_config.sh "$CONFIG"

mkdir -p "$OUTPUT_TABLE" "$OUTPUT_REPORT" "$TEMP_DIR" "$LOG_DIR"

"$PYTHON_BIN" python/0-build_calibration_membership.py \
    --id-hap "$ID_HAP" \
    --phylotree-json "$PHYLOTREE_JSON" \
    --calibration-nodes "$CALIBRATION_NODES" \
    --output-table "$OUTPUT_TABLE" \
    --output-report "$OUTPUT_REPORT" \
    --temp-dir "$TEMP_DIR" \
    --overwrite "$OVERWRITE" \
    2>&1 | tee "$LOG_DIR/calibration_membership.log"

#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

CONFIG="conf/1-lsd2.yaml"
source script/load_config.sh "$CONFIG"

mkdir -p "$OUTPUT_TABLE" "$OUTPUT_REPORT" "$TEMP_DIR" "$LOG_DIR"

"$PYTHON_BIN" python/1-1-run_lsd2.py \
    --fasta "$FASTA" \
    --ml-tree "$ML_TREE" \
    --calibration-membership "$CALIBRATION_MEMBERSHIP" \
    --calibration-nodes "$CALIBRATION_NODES" \
    --output-table "$OUTPUT_TABLE" \
    --output-report "$OUTPUT_REPORT" \
    --temp-dir "$TEMP_DIR" \
    --iqtree-bin "$IQTREE_BIN" \
    --iqtree-model "$IQTREE_MODEL" \
    --iqtree-threads "$IQTREE_THREADS" \
    --iqtree-outgroup "$IQTREE_OUTGROUP" \
    --iqtree-date-options "$IQTREE_DATE_OPTIONS" \
    --overwrite "$OVERWRITE" \
    2>&1 | tee "$LOG_DIR/lsd2.log"

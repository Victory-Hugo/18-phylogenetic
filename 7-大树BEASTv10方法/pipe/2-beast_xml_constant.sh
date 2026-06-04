#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

CONFIG="conf/2-beast_xml_constant.yaml"
source script/load_config.sh "$CONFIG"

mkdir -p "$OUTPUT_TABLE" "$OUTPUT_REPORT" "$TEMP_DIR" "$LOG_DIR"

"$PYTHON_BIN" python/2-1-build_beast_xml.py \
    --ml-tree "$ML_TREE" \
    --starting-tree "$STARTING_TREE" \
    --calibration-membership "$CALIBRATION_MEMBERSHIP" \
    --input-validation "$INPUT_VALIDATION" \
    --template-xml "$TEMPLATE_XML" \
    --output-table "$OUTPUT_TABLE" \
    --output-report "$OUTPUT_REPORT" \
    --temp-dir "$TEMP_DIR" \
    --output-prefix "$BEAST_OUTPUT_PREFIX" \
    --scale "$BEAST_SCALE" \
    --chain-length "$BEAST_CHAIN_LENGTH" \
    --log-every "$BEAST_LOG_EVERY" \
    --tree-log-every "$BEAST_TREE_LOG_EVERY" \
    --clock-rate "$BEAST_CLOCK_RATE" \
    --skygrid-points "$BEAST_SKYGRID_POINTS" \
    --skygrid-init-logpop "$BEAST_SKYGRID_INIT_LOGPOP" \
    --calibration-nodes "$CALIBRATION_NODES" \
    --overwrite "$OVERWRITE" \
    2>&1 | tee "$LOG_DIR/build_beast_xml.log"

# ---------------------------------------------------------------------------
# BEAST 解析校验：对极短 chain 的校验 XML 实跑几步，确认能被 BEAST v10.5.0
# 解析、模型可构建、posterior 非 -Inf（捕获 Unrecognised element / parser 报错）。
# ---------------------------------------------------------------------------
VALIDATE_XML="$TEMP_DIR/beast_validate/${BEAST_OUTPUT_PREFIX}.validate.xml"
echo "[validate] 用 BEAST 实跑校验 $VALIDATE_XML"
if "$BEAST_BIN" -overwrite -seed 1 -working "$VALIDATE_XML" \
    > "$LOG_DIR/beast_validate.log" 2>&1; then
    echo "[validate] BEAST 校验通过"
else
    echo "[validate] BEAST 校验失败，详见 $LOG_DIR/beast_validate.log" >&2
    tail -40 "$LOG_DIR/beast_validate.log" >&2
    exit 1
fi

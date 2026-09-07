#!/usr/bin/env bash
# 步骤 1：Y 染色体单倍群分层频率统计与分级堆叠柱状图
set -euo pipefail

STEP="1-haplogroup_level_frequency"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

source script/load_config.sh
eval "$(parse_yaml "conf/${STEP}.yaml")"

mkdir -p "$OUTPUT_RESULT" "$OUTPUT_FIGURE"

"$PYTHON_BIN" script/1-1-haplogroup_level_frequency.py \
    --input-path        "$INPUT" \
    --output-result     "$OUTPUT_RESULT" \
    --id-column         "$ID_COLUMN" \
    --haplogroup-column "$HAPLOGROUP_COLUMN" \
    --population-column "$POPULATION_COLUMN" \
    --label-column      "$LABEL_COLUMN" \
    --max-level         "$MAX_LEVEL" \
    --min-sample-size   "$MIN_SAMPLE_SIZE" \
    --population-order  "$POPULATION_ORDER" \
    --cluster-level     "$CLUSTER_LEVEL"

"$RSCRIPT_BIN" script/1-2-plot_haplogroup_level_frequency.R \
    --input-tsv            "${OUTPUT_RESULT}/⭐1-1-Haplogroup-Frequency-By-Level.tsv" \
    --color-tsv            "$COLOR" \
    --label-color-tsv      "$LABEL_COLOR" \
    --output-figure        "$OUTPUT_FIGURE" \
    --width-per-population "$WIDTH_PER_POPULATION" \
    --panel-height         "$PANEL_HEIGHT" \
    --legend-max-rows      "$LEGEND_MAX_ROWS" \
    --legend-max-categories "$LEGEND_MAX_CATEGORIES" \
    --dpi                  "$DPI"

echo "步骤 ${STEP} 完成"

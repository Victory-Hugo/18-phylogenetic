#!/usr/bin/env bash
# 步骤 2：交叉验证选择最优系统发育分辨率并绘制对应柱状图
set -euo pipefail

STEP="2-optimal_depth"
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

source script/load_config.sh
eval "$(parse_yaml "conf/${STEP}.yaml")"

mkdir -p "$OUTPUT_RESULT" "$OUTPUT_FIGURE"

"$PYTHON_BIN" script/2-1-optimal_depth.py --self-test

"$PYTHON_BIN" script/2-1-optimal_depth.py \
    --input-path          "$INPUT" \
    --output-result       "$OUTPUT_RESULT" \
    --output-figure       "$OUTPUT_FIGURE" \
    --id-column           "$ID_COLUMN" \
    --haplogroup-column   "$HAPLOGROUP_COLUMN" \
    --population-column   "$POPULATION_COLUMN" \
    --label-column        "$LABEL_COLUMN" \
    --max-level           "$MAX_LEVEL" \
    --min-sample-size     "$MIN_SAMPLE_SIZE" \
    --smoothing-alpha     "$SMOOTHING_ALPHA" \
    --cv-folds            "$CV_FOLDS" \
    --cv-repeats          "$CV_REPEATS" \
    --random-seed         "$RANDOM_SEED" \
    --lambda-min          "$LAMBDA_MIN" \
    --lambda-max          "$LAMBDA_MAX" \
    --lambda-steps        "$LAMBDA_STEPS" \
    --min-resolution-rate "$MIN_RESOLUTION_RATE" \
    --cluster-level       "$CLUSTER_LEVEL"

"$RSCRIPT_BIN" script/2-2-plot_optimal_depth.R \
    --frequency-tsv        "${OUTPUT_FIGURE}/⭐2-4-Optimal-Depth-Frequency.tsv" \
    --summary-tsv          "${OUTPUT_RESULT}/⭐2-1-Optimal-Depth-Selection.tsv" \
    --color-tsv            "$COLOR" \
    --label-color-tsv      "$LABEL_COLOR" \
    --output-figure        "$OUTPUT_FIGURE" \
    --width-per-population "$WIDTH_PER_POPULATION" \
    --panel-height         "$PANEL_HEIGHT" \
    --legend-max-rows      "$LEGEND_MAX_ROWS" \
    --dpi                  "$DPI"

echo "步骤 ${STEP} 完成"

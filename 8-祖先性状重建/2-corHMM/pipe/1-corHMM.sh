#!/usr/bin/env bash
# 第1步：比较离散性状模型、执行随机映射并输出逐事件位置。

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONF="$BASE_DIR/conf/1-corHMM.yaml"
source "$BASE_DIR/script/load_config.sh"

RSCRIPT_BIN="$(yaml_get "$CONF" tools rscript_bin)"
INPUT_TREE="$BASE_DIR/$(yaml_get "$CONF" paths input_tree)"
INPUT_TRAITS="$BASE_DIR/$(yaml_get "$CONF" paths input_traits)"
OUTPUT_MODEL="$BASE_DIR/$(yaml_get "$CONF" paths output_model)"
OUTPUT_SIMMAP="$BASE_DIR/$(yaml_get "$CONF" paths output_simmap)"
OUTPUT_TRANS="$BASE_DIR/$(yaml_get "$CONF" paths output_trans)"
TEMP_DIR="$BASE_DIR/$(yaml_get "$CONF" paths temp)"
LOG_PATH="$BASE_DIR/$(yaml_get "$CONF" paths log)"

mkdir -p "$OUTPUT_MODEL" "$OUTPUT_SIMMAP" "$OUTPUT_TRANS" "$TEMP_DIR" "$(dirname "$LOG_PATH")"

"$RSCRIPT_BIN" --vanilla "$BASE_DIR/R/1-corHMM.R" \
  --input-tree "$INPUT_TREE" \
  --input-traits "$INPUT_TRAITS" \
  --output-model "$OUTPUT_MODEL" \
  --output-simmap "$OUTPUT_SIMMAP" \
  --output-trans "$OUTPUT_TRANS" \
  --temp-dir "$TEMP_DIR" \
  --id-column "$(yaml_get "$CONF" columns id)" \
  --trait-column "$(yaml_get "$CONF" columns trait)" \
  --models "$(yaml_get "$CONF" models candidates)" \
  --criterion "$(yaml_get "$CONF" models selection_criterion)" \
  --root-prior "$(yaml_get "$CONF" models root_prior)" \
  --rate-categories "$(yaml_get "$CONF" models rate_categories)" \
  --node-states "$(yaml_get "$CONF" models node_states)" \
  --n-starts "$(yaml_get "$CONF" models n_starts)" \
  --model-cores "$(yaml_get "$CONF" models n_cores)" \
  --n-simulations "$(yaml_get "$CONF" simmap n_simulations)" \
  --simmap-cores "$(yaml_get "$CONF" simmap n_cores)" \
  --max-attempt "$(yaml_get "$CONF" simmap max_attempt)" \
  --resolve-polytomies "$(yaml_get "$CONF" tree resolve_polytomies)" \
  --branch-floor-ratio "$(yaml_get "$CONF" tree branch_floor_ratio)" \
  --scale-branch-lengths "$(yaml_get "$CONF" tree scale_branch_lengths)" \
  --seed "$(yaml_get "$CONF" runtime seed)" \
  --overwrite "$(yaml_get "$CONF" runtime overwrite)" \
  2>&1 | tee "$LOG_PATH"

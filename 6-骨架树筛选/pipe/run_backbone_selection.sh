#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/Config.yaml"
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
META_TSV=$(resolve_path "$(config_get paths.meta_tsv)")
BASIC_INFO_TSV=$(resolve_path "$(config_get paths.basic_info_tsv)")
VCF_GZ=$(resolve_path "$(config_get paths.vcf_gz)")
OUTPUT_DIR=$(resolve_path "$(config_get paths.output_dir)")
TIERS=$(config_get runtime.tiers)
DISTANCE_MODE=$(config_get_default runtime.distance_mode "alt_hamming")
SEED_MODE=$(config_get_default runtime.seed_mode "deep_lineage_cover")
DEEP_LINEAGE_LABELS=$(config_get_default runtime.deep_lineage_labels "auto")
MIN_QUALITY=$(config_get_default runtime.min_quality "0.0")

for required in "$META_TSV" "$BASIC_INFO_TSV" "$VCF_GZ"; do
    if [[ ! -f "$required" ]]; then
        echo "[ERROR] Required input not found: $required" >&2
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

"$PYTHON_BIN" "$PROJECT_ROOT/python/select_backbone_by_variation.py" \
    --meta "$META_TSV" \
    --vcf "$VCF_GZ" \
    --basic-info "$BASIC_INFO_TSV" \
    --output-dir "$OUTPUT_DIR" \
    --tiers "$TIERS" \
    --distance-mode "$DISTANCE_MODE" \
    --seed-mode "$SEED_MODE" \
    --deep-lineage-labels "$DEEP_LINEAGE_LABELS" \
    --min-quality "$MIN_QUALITY"

echo "[OK] Backbone sample selection finished."

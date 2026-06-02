#!/usr/bin/env bash
# 脚本说明：
#   1-sample_select.sh
#   用途：对全部 mtDNA 样本进行两层筛选，输出带有筛选理由的样本列表。
#
# 用法：
#   bash pipe/1-sample_select.sh [--config /path/to/config.yaml]
#   默认配置文件：conf/1-sample_select.yaml
#
# 两步流程：
#   Stage 1 — 1-1-tsv_filter.py
#     读取 meta/总表.tsv，应用三项硬过滤（QC_Other、QC_Haplogrep、Sequence_method），
#     输出通过硬过滤的样本列表到 temp/。
#   Stage 2 — 1-2-vcf_fps_select.py
#     按主单倍群分块，从 VCF 提取基因型（bcftools），计算成对差异，
#     折叠相同单倍型后执行自适应 FPS 降采样，输出带理由的最终样本列表。
#
# 主要输出：
#   output/1-sample_select/1-table/selected_samples.tsv  ← 最终筛选结果
#   output/1-sample_select/3-report/selection_summary.tsv ← 各块统计摘要
#   temp/1-sample_select/tsv_filtered.tsv                 ← 第一层中间结果

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
DEFAULT_CONFIG="$PROJECT_ROOT/conf/1-sample_select.yaml"
CONFIG_PATH="$DEFAULT_CONFIG"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG_PATH="$2"; shift 2 ;;
        *) echo "[ERROR] 未知参数: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "[ERROR] 配置文件不存在: $CONFIG_PATH" >&2
    exit 1
fi

source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"
bash "$PROJECT_ROOT/script/check_env.sh" "$CONFIG_PATH"

# ── 路径解析 ───────────────────────────────────────────────────────────────────
resolve_path() {
    local p="$1"
    if [[ "$p" = /* ]]; then printf '%s\n' "$p"
    else printf '%s\n' "$PROJECT_ROOT/$p"
    fi
}

require_var() {
    local v="$1"
    if [[ -z "${!v:-}" ]]; then
        echo "[ERROR] 缺少配置项: $v" >&2; exit 1
    fi
}

require_var CFG_PATHS_TSV
require_var CFG_PATHS_VCF
require_var CFG_PATHS_OUTPUT_TABLE
require_var CFG_PATHS_OUTPUT_REPORT
require_var CFG_PATHS_TEMP
require_var CFG_TOOLS_PYTHON_BIN
require_var CFG_TOOLS_BCFTOOLS_BIN
require_var CFG_RUNTIME_QC_MIN
require_var CFG_RUNTIME_MIN_DIST
require_var CFG_RUNTIME_MAX_TIPS

PYTHON_BIN=$(resolve_path   "$CFG_TOOLS_PYTHON_BIN")
BCFTOOLS_BIN=$(resolve_path "$CFG_TOOLS_BCFTOOLS_BIN")
TSV=$(resolve_path          "$CFG_PATHS_TSV")
VCF=$(resolve_path          "$CFG_PATHS_VCF")
OUTPUT_TABLE=$(resolve_path  "$CFG_PATHS_OUTPUT_TABLE")
OUTPUT_REPORT=$(resolve_path "$CFG_PATHS_OUTPUT_REPORT")
TEMP_DIR=$(resolve_path      "$CFG_PATHS_TEMP")

QC_MIN="$CFG_RUNTIME_QC_MIN"
MIN_DIST="$CFG_RUNTIME_MIN_DIST"
MAX_TIPS="$CFG_RUNTIME_MAX_TIPS"
OVERWRITE="${CFG_RUNTIME_OVERWRITE:-false}"

mkdir -p "$OUTPUT_TABLE" "$OUTPUT_REPORT" "$TEMP_DIR"

TSV_FILTERED="$TEMP_DIR/tsv_filtered.tsv"
SELECTED_OUT="$OUTPUT_TABLE/selected_samples.tsv"
SUMMARY_OUT="$OUTPUT_REPORT/selection_summary.tsv"

# ── Stage 1：TSV 硬过滤 ────────────────────────────────────────────────────────
echo "[Stage 1/2] TSV 硬过滤..."
"$PYTHON_BIN" "$PROJECT_ROOT/python/1-1-tsv_filter.py" \
    --tsv      "$TSV" \
    --output   "$TSV_FILTERED" \
    --qc-min   "$QC_MIN"

# ── Stage 2：VCF 成对差异 + FPS 降采样 ────────────────────────────────────────
echo "[Stage 2/2] VCF 成对差异 + FPS 降采样..."
"$PYTHON_BIN" "$PROJECT_ROOT/python/1-2-vcf_fps_select.py" \
    --filtered-list "$TSV_FILTERED" \
    --vcf           "$VCF" \
    --output        "$SELECTED_OUT" \
    --summary       "$SUMMARY_OUT" \
    --temp          "$TEMP_DIR" \
    --bcftools      "$BCFTOOLS_BIN" \
    --min-dist      "$MIN_DIST" \
    --max-tips      "$MAX_TIPS" \
    --overwrite     "$OVERWRITE"

echo "[OK] 1-sample_select 完成。结果: $SELECTED_OUT"

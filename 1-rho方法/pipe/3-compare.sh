#!/usr/bin/env bash
# =============================================================================
# 3-compare.sh — Stage 3：文献单倍群年龄与本项目 ρ/ML 年龄比较
# 用法：bash pipe/3-compare.sh [conf/stage_3_compare.yaml]
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG="${1:-$BASE_DIR/conf/stage_3_compare.yaml}"

source "$BASE_DIR/script/load_config.sh" "$CONFIG"

BASE="$CFG_PROJECT_BASE_DIR"
PYTHON_DIR="$BASE/$CFG_PATHS_PYTHON_DIR"
DOWNLOAD_DIR="$BASE/$CFG_PATHS_DOWNLOAD_DIR"
PROJECT_COMPARISON="$BASE/$CFG_PATHS_PROJECT_COMPARISON"
RESULTS="$BASE/$CFG_PATHS_RESULTS_DIR"
LOG_DIR="$BASE/$CFG_PATHS_LOG_DIR"
REPORT_DIR="$BASE/$CFG_PATHS_REPORT_DIR"

CONDA_BIN="$CFG_TOOLS_CONDA_BIN"
ENV_PREFIX="$CFG_TOOLS_PYTHON_ENV_PREFIX"
LARGE_DIFF_KYA="$CFG_RUNTIME_LARGE_DIFF_KYA"
LARGE_REL_DIFF_PCT="$CFG_RUNTIME_LARGE_REL_DIFF_PCT"

PYTHON="$CONDA_BIN run -p $ENV_PREFIX python3"
REPORT="$REPORT_DIR/文献与本项目单倍群年龄比较报告.md"

mkdir -p "$RESULTS" "$LOG_DIR" "$REPORT_DIR"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

[[ -d "$DOWNLOAD_DIR" ]] || die "文献输入目录不存在：$DOWNLOAD_DIR"
[[ -f "$PROJECT_COMPARISON" ]] || die "本项目比较表不存在：$PROJECT_COMPARISON"
[[ -f "$PYTHON_DIR/s11_compare_literature_ages.py" ]] || die "Python 模块不存在：$PYTHON_DIR/s11_compare_literature_ages.py"

log "=========================================="
log "Stage 3 文献年龄比较 pipeline 启动"
log "BASE_DIR:           $BASE"
log "DOWNLOAD_DIR:       $DOWNLOAD_DIR"
log "PROJECT_COMPARISON: $PROJECT_COMPARISON"
log "RESULTS:            $RESULTS"
log "REPORT:             $REPORT"
log "=========================================="

log "[Step 11] 解析文献年龄，比较本项目 ρ/ML 年龄，并生成中文报告..."
$PYTHON "$PYTHON_DIR/s11_compare_literature_ages.py" \
    --download-dir "$DOWNLOAD_DIR" \
    --project-comparison "$PROJECT_COMPARISON" \
    --results-dir "$RESULTS" \
    --report "$REPORT" \
    --large-diff-kya "$LARGE_DIFF_KYA" \
    --large-rel-diff-pct "$LARGE_REL_DIFF_PCT" \
    2>&1 | tee "$LOG_DIR/s11_compare_literature_ages.log"
log "[Step 11] 完成"

log "=========================================="
log "Stage 3 比较完成！"
log "主要结果文件："
log "  文献统一表:     $RESULTS/literature_age_normalized.tsv"
log "  逐条比较表:     $RESULTS/literature_vs_project_comparison.tsv"
log "  汇总统计表:     $RESULTS/literature_vs_project_summary.tsv"
log "  未匹配记录:     $RESULTS/literature_unmatched.tsv"
log "  最大差异记录:   $RESULTS/top_literature_project_differences.tsv"
log "  中文报告:       $REPORT"
log "=========================================="

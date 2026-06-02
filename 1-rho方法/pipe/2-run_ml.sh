#!/usr/bin/env bash
# =============================================================================
# run_ml.sh — Stage 2：IQ-TREE 3 ML 单倍群年龄估计总控脚本
# 用法：bash pipe/run_ml.sh [conf/stage_2_ml.yaml]
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG="${1:-$BASE_DIR/conf/stage_2_ml.yaml}"

source "$BASE_DIR/script/load_config.sh" "$CONFIG"

BASE="$CFG_PROJECT_BASE_DIR"
PYTHON_DIR="$BASE/$CFG_PATHS_PYTHON_DIR"
INTERMEDIATE="$BASE/$CFG_PATHS_INTERMEDIATE_DIR"
RESULTS="$BASE/$CFG_PATHS_RESULTS_DIR"
LOG_DIR="$BASE/$CFG_PATHS_LOG_DIR"
ML_DIR="$INTERMEDIATE/ml"

ANC_FASTA="$BASE/$CFG_PATHS_ANC_FASTA"
HAP_GRAPH="$BASE/$CFG_PATHS_HAP_GRAPH"
RHO_DATING="$BASE/$CFG_PATHS_RHO_DATING"
HALIGN4="$BASE/$CFG_PATHS_HALIGN4_BIN"

CONDA_BIN="$CFG_TOOLS_CONDA_BIN"
ENV_PREFIX="$CFG_TOOLS_PYTHON_ENV_PREFIX"
IQTREE3="$CFG_TOOLS_IQTREE3_BIN"
THREADS="$CFG_RUNTIME_THREADS"
ROOT_AGE="$CFG_RUNTIME_ROOT_AGE"
ROOT_NODE="$CFG_RUNTIME_ROOT_NODE"

PYTHON="$CONDA_BIN run -p $ENV_PREFIX python3"

mkdir -p "$ML_DIR" "$RESULTS" "$LOG_DIR"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

log "=========================================="
log "Stage 2 ML 单倍群年龄估计 pipeline 启动"
log "BASE_DIR:    $BASE"
log "ML_DIR:      $ML_DIR"
log "RESULTS:     $RESULTS"
log "THREADS:     $THREADS"
log "ROOT_AGE:    $ROOT_AGE 年"
log "=========================================="

# =============================================================================
# Step 8：净化祖先序列名称，提取叶节点 FASTA
# =============================================================================
LEAVES_FASTA="$ML_DIR/ancestors_sanitized_leaves.fasta"
NAME_MAP="$ML_DIR/name_sanitize_map.tsv"
if [[ ! -f "$NAME_MAP" ]]; then
    log "[Step 8] 净化名称，提取叶节点序列..."
    $PYTHON "$PYTHON_DIR/s08_prepare_ml_inputs.py" \
        --fasta        "$ANC_FASTA" \
        --hap-graph    "$HAP_GRAPH" \
        --output-fasta "$LEAVES_FASTA" \
        --output-map   "$NAME_MAP" \
        2>&1 | tee "$LOG_DIR/s08_prepare_ml_inputs.log"
    log "[Step 8] 完成"
else
    log "[Step 8] 跳过（已存在 $NAME_MAP）"
fi

# =============================================================================
# Step 8b：多序列比对（halign4，以根节点序列最近的叶为参考，约 1 分钟）
# =============================================================================
ALIGNED_FASTA="$ML_DIR/ancestors_aligned_leaves.fasta"
if [[ ! -f "$ALIGNED_FASTA" ]]; then
    log "[Step 8b] 多序列比对（halign4，叶节点序列）..."
    "$HALIGN4" \
        "$LEAVES_FASTA" \
        "$ALIGNED_FASTA" \
        -t "$THREADS" \
        2>&1 | tee "$LOG_DIR/s08b_halign4.log"
    log "[Step 8b] 完成"
else
    log "[Step 8b] 跳过（已存在 $ALIGNED_FASTA）"
fi

# =============================================================================
# Step 9：生成固定拓扑 Newick
# =============================================================================
NEWICK="$ML_DIR/phylotree_fixed.nwk"
if [[ ! -f "$NEWICK" ]]; then
    log "[Step 9] 生成 PhyloTree 固定拓扑 Newick..."
    $PYTHON "$PYTHON_DIR/s09_phylotree_to_newick.py" \
        --graph  "$HAP_GRAPH" \
        --map    "$NAME_MAP" \
        --output "$NEWICK" \
        --root   "$ROOT_NODE" \
        2>&1 | tee "$LOG_DIR/s09_phylotree_to_newick.log"
    log "[Step 9] 完成"
else
    log "[Step 9] 跳过（已存在 $NEWICK）"
fi

# =============================================================================
# Step 10：IQ-TREE 3（固定拓扑 ML 枝长估计，约 1-3 小时）
# =============================================================================
TREEFILE="$ML_DIR/anc_tree.treefile"
PREFIX="$ML_DIR/anc_tree"
if [[ ! -f "$TREEFILE" ]]; then
    log "[Step 10] 运行 IQ-TREE 3（固定拓扑，GTR+G，约 1-3 小时）..."
    "$IQTREE3" \
        -s  "$ALIGNED_FASTA" \
        -te "$NEWICK" \
        -m  GTR+G \
        -T  "$THREADS" \
        --prefix "$PREFIX" \
        --redo \
        2>&1 | tee "$LOG_DIR/s10_iqtree3.log"
    log "[Step 10] 完成"
else
    log "[Step 10] 跳过（已存在 $TREEFILE）"
fi

# =============================================================================
# Step 11：分子钟校准 + 与 ρ 结果比较
# =============================================================================
ML_DATING="$RESULTS/ml_dating.tsv"
ML_COMPARISON="$RESULTS/ml_vs_rho_comparison.tsv"
log "[Step 11] 分子钟校准与 ρ 结果比较..."
$PYTHON "$PYTHON_DIR/s10_ml_dating.py" \
    --treefile          "$TREEFILE" \
    --map               "$NAME_MAP" \
    --rho               "$RHO_DATING" \
    --hap-graph         "$HAP_GRAPH" \
    --root-age          "$ROOT_AGE" \
    --output-dating     "$ML_DATING" \
    --output-comparison "$ML_COMPARISON" \
    2>&1 | tee "$LOG_DIR/s11_ml_dating.log"
log "[Step 11] 完成"

# =============================================================================
# Step 11b：祖先-后代 CI 一致性检查（ML 结果）
# =============================================================================
log "[Step 11b] ML dating CI 一致性检查（ml_ci bounds → 降级 confidence）..."
$PYTHON "$PYTHON_DIR/s05b_ci_consistency.py" \
    --mode      ml \
    --input     "$ML_DATING" \
    --hap-graph "$HAP_GRAPH" \
    --output    "$ML_DATING" \
    2>&1 | tee "$LOG_DIR/s11b_ci_consistency.log"
log "[Step 11b] 完成"

# =============================================================================
# 完成汇报
# =============================================================================
log "=========================================="
log "Stage 2 ML pipeline 完成！"
log "主要结果文件："
log "  ML dating:      $ML_DATING"
log "  rho 比较:       $ML_COMPARISON"
log "中间产物："
log "  叶节点比对 FASTA: $ALIGNED_FASTA"
log "  固定拓扑 Newick: $NEWICK"
log "  IQ-TREE treefile: $TREEFILE"
log "=========================================="

#!/usr/bin/env bash
# =============================================================================
# run_all.sh — Stage 1：mtDNA单倍群年龄估计（ρ统计量 + Soares 2009）总控脚本
# 用法：bash pipe/run_all.sh [conf/stage_1_rho.yaml]
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG="${1:-$BASE_DIR/conf/stage_1_rho.yaml}"

# 加载配置
source "$BASE_DIR/script/load_config.sh" "$CONFIG"

# 从配置变量构建所有路径（均相对于base_dir展开）
BASE="$CFG_PROJECT_BASE_DIR"
PYTHON_DIR="$BASE/$CFG_PATHS_PYTHON_DIR"
INTERMEDIATE="$BASE/$CFG_PATHS_INTERMEDIATE_DIR"
RESULTS="$BASE/$CFG_PATHS_RESULTS_DIR"
LOG_DIR="$BASE/$CFG_PATHS_LOG_DIR"

VCF="$BASE/$CFG_PATHS_VCF"
PHYLOTREE_JSON="$BASE/$CFG_PATHS_PHYLOTREE_JSON"
BACK_MUTATION="$BASE/$CFG_PATHS_BACK_MUTATION_TSV"
HAPLOGREP="$BASE/$CFG_PATHS_HAPLOGREP_TSV"

CONDA_BIN="$CFG_TOOLS_CONDA_BIN"
ENV_PREFIX="$CFG_TOOLS_PYTHON_ENV_PREFIX"
BCFTOOLS="$CFG_TOOLS_BCFTOOLS_BIN"
VALIDATE_NODES="$CFG_RUNTIME_VALIDATE_NODES"
HAPLOGREP_QUALITY_MIN="$CFG_RUNTIME_HAPLOGREP_QUALITY_MIN"

PYTHON="$CONDA_BIN run -p $ENV_PREFIX python3"

# 创建输出目录
mkdir -p "$INTERMEDIATE" "$RESULTS" "$LOG_DIR"

# 时间戳日志函数
log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

log "=========================================="
log "Stage 1 rho 单倍群年龄估计 pipeline 启动"
log "BASE_DIR:     $BASE"
log "INTERMEDIATE: $INTERMEDIATE"
log "RESULTS:      $RESULTS"
log "Haplogrep Quality >= $HAPLOGREP_QUALITY_MIN"
log "=========================================="

# =============================================================================
# Step 0：提取VCF全样本列表
# =============================================================================
ALL_SAMPLES="$INTERMEDIATE/all_samples.txt"
if [[ ! -f "$ALL_SAMPLES" ]]; then
    log "[Step 0] 提取VCF全样本列表..."
    "$BCFTOOLS" query -l "$VCF" > "$ALL_SAMPLES"
    log "  样本数: $(wc -l < "$ALL_SAMPLES")"
else
    log "[Step 0] 跳过（已存在 $ALL_SAMPLES）"
fi

# =============================================================================
# Step 1：构建位点掩码
# =============================================================================
SITE_MASKS="$INTERMEDIATE/site_masks.json"
if [[ ! -f "$SITE_MASKS" ]]; then
    log "[Step 1] 构建各区域位点掩码..."
    $PYTHON "$PYTHON_DIR/s01_build_site_masks.py" \
        --vcf             "$VCF" \
        --back-mutation   "$BACK_MUTATION" \
        --output          "$SITE_MASKS" \
        --bcftools        "$BCFTOOLS" \
        2>&1 | tee "$LOG_DIR/s01_build_site_masks.log"
else
    log "[Step 1] 跳过（已存在 $SITE_MASKS）"
fi

# =============================================================================
# Step 2：构建单倍群图
# =============================================================================
HAP_GRAPH="$INTERMEDIATE/haplogroup_graph.json"
if [[ ! -f "$HAP_GRAPH" ]]; then
    log "[Step 2] 构建单倍群祖先-后代图..."
    $PYTHON "$PYTHON_DIR/s02_build_haplogroup_graph.py" \
        --phylotree "$PHYLOTREE_JSON" \
        --output    "$HAP_GRAPH" \
        2>&1 | tee "$LOG_DIR/s02_build_haplogroup_graph.log"
else
    log "[Step 2] 跳过（已存在 $HAP_GRAPH）"
fi

# =============================================================================
# Step 3：建立样本-单倍群-后代映射
# =============================================================================
MODERN_SAMPLES="$INTERMEDIATE/modern_samples.txt"
ANCESTOR_SAMPLES="$INTERMEDIATE/ancestor_samples.txt"
SAMPLE_HAP="$INTERMEDIATE/sample_to_haplogroup.tsv"
HAP_DESC="$INTERMEDIATE/haplogroup_to_descendants.json"
HAPLOGREP_FILTERED="$INTERMEDIATE/haplogrep_quality_filtered.tsv"
if [[ ! -f "$HAP_DESC" ]]; then
    log "[Step 3] 按Quality阈值筛选Haplogrep3结果..."
    awk -F'\t' -v OFS='\t' -v qmin="$HAPLOGREP_QUALITY_MIN" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col = tolower($i)
                gsub(/^[[:space:]"]+|[[:space:]"]+$/, "", col)
                if (col == "quality") quality_col = i
            }
            if (!quality_col) {
                print "[ERROR] Haplogrep3结果缺少Quality列，无法按Quality阈值筛选" > "/dev/stderr"
                exit 1
            }
            print
            next
        }
        NR > 1 {
            total++
            q = $quality_col
            gsub(/"/, "", q)
            if ((q + 0) >= qmin) { print; kept++ }
        }
        END {
            if (quality_col) {
                printf "  Quality >= %s: 保留 %d/%d 条记录\n", qmin, kept, total > "/dev/stderr"
            }
        }
    ' "$HAPLOGREP" > "$HAPLOGREP_FILTERED"
    log "[Step 3] 建立样本-单倍群-后代映射..."
    $PYTHON "$PYTHON_DIR/s03_build_sample_mapping.py" \
        --haplogrep         "$HAPLOGREP_FILTERED" \
        --haplogroup-graph  "$HAP_GRAPH" \
        --vcf-sample-list   "$ALL_SAMPLES" \
        --output-modern     "$MODERN_SAMPLES" \
        --output-ancestor   "$ANCESTOR_SAMPLES" \
        --output-sample-hap "$SAMPLE_HAP" \
        --output-hap-desc   "$HAP_DESC" \
        2>&1 | tee "$LOG_DIR/s03_build_sample_mapping.log"
else
    log "[Step 3] 跳过（已存在 $HAP_DESC）"
fi

# =============================================================================
# Step 4：提取基因型矩阵（最耗时步骤，约5-15分钟）
# =============================================================================
GENO_MATRIX="$INTERMEDIATE/geno_matrix.npz"
GENO_POS="$INTERMEDIATE/geno_positions.npy"
GENO_SAMPLES="$INTERMEDIATE/geno_samples.txt"
if [[ ! -f "$GENO_MATRIX" ]]; then
    log "[Step 4] 从VCF提取基因型矩阵（需要数分钟）..."
    $PYTHON "$PYTHON_DIR/s04_extract_genotype_matrix.py" \
        --vcf              "$VCF" \
        --site-masks       "$SITE_MASKS" \
        --output-matrix    "$GENO_MATRIX" \
        --output-positions "$GENO_POS" \
        --output-samples   "$GENO_SAMPLES" \
        2>&1 | tee "$LOG_DIR/s04_extract_genotype_matrix.log"
else
    log "[Step 4] 跳过（已存在 $GENO_MATRIX）"
fi

# =============================================================================
# Step 5：计算ρ dating（主计算步骤）
# =============================================================================
RHO_DATING="$RESULTS/rho_dating.tsv"
if [[ ! -f "$RHO_DATING" ]]; then
    log "[Step 5] 计算全树ρ dating..."
    $PYTHON "$PYTHON_DIR/s05_calc_rho_dating.py" \
        --geno-matrix    "$GENO_MATRIX" \
        --geno-positions "$GENO_POS" \
        --geno-samples   "$GENO_SAMPLES" \
        --site-masks     "$SITE_MASKS" \
        --hap-desc       "$HAP_DESC" \
        --output         "$RHO_DATING" \
        2>&1 | tee "$LOG_DIR/s05_calc_rho_dating.log"
else
    log "[Step 5] 跳过（已存在 $RHO_DATING）"
fi

# =============================================================================
# Step 5b：祖先-后代 CI 一致性检查（complete 区域）
# =============================================================================
log "[Step 5b] CI 一致性检查（complete 区域 → 降级 confidence）..."
$PYTHON "$PYTHON_DIR/s05b_ci_consistency.py" \
    --mode      rho \
    --input     "$RHO_DATING" \
    --hap-graph "$HAP_GRAPH" \
    --output    "$RHO_DATING" \
    2>&1 | tee "$LOG_DIR/s05b_ci_consistency.log"

# =============================================================================
# Step 6：验证关键节点
# =============================================================================
VALIDATION="$RESULTS/validation_key_nodes.tsv"
log "[Step 6] 验证关键单倍群节点年龄..."
$PYTHON "$PYTHON_DIR/s06_validate_key_nodes.py" \
    --rho-dating     "$RHO_DATING" \
    --validate-nodes "$VALIDATE_NODES" \
    --output         "$VALIDATION" \
    2>&1 | tee "$LOG_DIR/s06_validate_key_nodes.log"

# =============================================================================
# Step 7：合并输出为宽格式
# =============================================================================
RHO_WIDE="$RESULTS/rho_dating_wide.tsv"
RHO_SUMMARY="$RESULTS/rho_dating_summary.tsv"
log "[Step 7] 生成宽格式结果和汇总统计..."
$PYTHON "$PYTHON_DIR/s07_merge_output.py" \
    --rho-dating     "$RHO_DATING" \
    --output-wide    "$RHO_WIDE" \
    --output-summary "$RHO_SUMMARY" \
    2>&1 | tee "$LOG_DIR/s07_merge_output.log"

# =============================================================================
# 完成汇报
# =============================================================================
log "=========================================="
log "Stage 1 pipeline 完成！"
log "主要结果文件："
log "  ρ dating（长格式）: $RHO_DATING"
log "  ρ dating（宽格式）: $RHO_WIDE"
log "  汇总统计:           $RHO_SUMMARY"
log "  关键节点验证:       $VALIDATION"
log "=========================================="

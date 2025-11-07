#!/bin/bash

# 退出条件：任何命令失败时立即退出
set -euo pipefail

# ================== 配置参数 ==================
BASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/2-Rho/"
# 存放所有 .vcf.gz 文件的目录
VCF_DIR="${BASE_DIR}/data"
# FASTA 输出根目录（脚本会自动按样本名创建子目录和前缀）
OUTPUT_DIR="${BASE_DIR}/data"
# Python 可执行文件路径
PYTHON_PATH="/home/luolintao/miniconda3/bin/python3"
# 脚本
SCRIPT_PATH="${BASE_DIR}/script/3-vcf→fasta_通过plink.py"
# 0转N修复脚本路径
FIX_ZEROS_SCRIPT="${BASE_DIR}/script/3-5-vcf2fasta_通过plink_0toN.py"
#? ID 与单倍群映射表可选可不选
# MAPPING_FILE="${BASE_DIR}/data/质量控制_ID_Hap.tsv"
# 并行线程数
THREADS=32
# 并行处理的最大任务数
PARALLEL_JOBS=8
# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# ================== 函数定义 ==================

function run_plink() {
    echo "[INFO] 使用 PLINK 处理 ${VCF_PATH} ..."
    /mnt/e/Scientifc_software/plink_linux_x86_64_20241022/plink \
        --vcf "${VCF_PATH}" \
        --threads "${THREADS}" \
        --recode \
        --double-id \
        --out "${OUT_PREFIX}"
    echo "[INFO] PLINK 处理完成：${OUT_PREFIX}.ped"
}

function ped_to_fasta() {
    local ped_file="${OUT_PREFIX}.ped"
    local fasta_file="${OUT_PREFIX}.fasta"

    echo "[INFO] 将 PED 转为 FASTA：${fasta_file} ..."
    awk '{
        printf(">%s\n", $1);
        seq = "";
        for (i = 7; i <= NF; i += 2) {
            seq = seq $i;
        }
        printf("%s\n", seq);
    }' "${ped_file}" > "${fasta_file}"
    echo "[INFO] PED 转 FASTA 完成：${fasta_file}"
}

function rename_fasta_ids() {
    local input_fasta="${OUT_PREFIX}.fasta"
    local output_fasta="${OUT_PREFIX}_rename.fasta"

    echo "[INFO] 调用 Python 脚本重命名：${output_fasta} ..."
    "${PYTHON_PATH}" "${SCRIPT_PATH}" \
        --mapping "${MAPPING_FILE}" \
        --input_fasta "${input_fasta}" \
        --output_fasta "${output_fasta}"
    echo "[INFO] 重命名完成：${output_fasta}"
}

function cleanup() {
    echo "[INFO] 清理中间文件：${OUT_PREFIX}.*"
    rm -f "${OUT_PREFIX}.ped" \
          "${OUT_PREFIX}.map" \
          "${OUT_PREFIX}.log" \
          "${OUT_PREFIX}.nosex"
}

function fix_zeros_in_fasta() {
    local fasta_dir="$1"
    echo "[INFO] 开始修复 FASTA 文件中的'0'字符（使用 parallel）..."
    "${PYTHON_PATH}" "${FIX_ZEROS_SCRIPT}" "${fasta_dir}"
    echo "[INFO] FASTA 修复完成！"
}

# ================== 主流程循环 ==================

echo "[START] 开始处理所有 VCF 文件..."

# 收集所有输出目录以便后续批量处理
FASTA_DIRS=()

for vcf in "${VCF_DIR}"/*.vcf.gz; do
    # 样本名前缀（去掉 .vcf.gz）
    SAMPLE=$(basename "${vcf}" .vcf.gz)
    VCF_PATH="${vcf}"
    OUT_PREFIX="${OUTPUT_DIR}/${SAMPLE}"

    echo "=============================================="
    echo "[START] 处理样本：${SAMPLE}"
    
    run_plink
    ped_to_fasta
    # 如果需要重命名 ID，请取消下一行注释
    # rename_fasta_ids
    cleanup

    echo "[DONE] 样本 ${SAMPLE} 完成"
    echo "=============================================="
    
    # 将输出目录加入列表
    FASTA_DIRS+=("${OUTPUT_DIR}")
done

echo "[ALL DONE] 所有 VCF 文件处理完毕。"

# ================== 使用 parallel 批量修复 FASTA 文件 ==================

if [ ${#FASTA_DIRS[@]} -gt 0 ]; then
    echo ""
    echo "=============================================="
    echo "[START] 使用 parallel 批量修复 FASTA 文件中的'0'字符..."
    
    # 去重
    UNIQUE_DIRS=($(printf '%s\n' "${FASTA_DIRS[@]}" | sort -u))
    
    # 使用 parallel 并行处理每个目录
    printf '%s\n' "${UNIQUE_DIRS[@]}" | parallel -j "${PARALLEL_JOBS}" \
        "${PYTHON_PATH}" "${FIX_ZEROS_SCRIPT}" {}
    
    echo "[DONE] 所有 FASTA 文件修复完成！"
    echo "=============================================="
fi

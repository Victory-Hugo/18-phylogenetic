# #!/bin/bash

# # 退出条件：任何命令失败时立即退出
# set -euo pipefail

# # ================== 配置参数 ==================
# VCF_PATH="/mnt/f/0_现代DNA处理流程/output/merged_biallelic_complete_annot_去除hots_fixed.vcf.gz"
# OUT_PREFIX="/mnt/f/0_现代DNA处理流程/output/Merge_FASTA/merged_biallelic_complete_annot_去除hots"
# PYTHON_PATH="/home/luolintao/miniconda3/bin/python3"
# SCRIPT_PATH="/mnt/f/0_现代DNA处理流程/script/6_vcf2fasta_通过plink.py"
# MAPPING_FILE="/home/luolintao/07_20K_CPGDP/1_单倍群分型/data/质量控制_ID_Hap.tsv"
# THREADS=32

# # ================== 函数定义 ==================

# function run_plink() {
#     echo "[INFO] 正在使用 PLINK 处理 VCF 文件..."
#     /mnt/e/Scientifc_software/plink_linux_x86_64_20241022/plink --vcf "$VCF_PATH" \
#         --threads "$THREADS" \
#         --recode \
#         --double-id \
#         --out "$OUT_PREFIX"
#     echo "[INFO] PLINK 处理完成！"
# }

# function ped_to_fasta() {
#     local ped_file="${OUT_PREFIX}.ped"
#     local fasta_file="${OUT_PREFIX}.fasta"

#     echo "[INFO] 正在将 PED 文件转为 FASTA 格式..."
#     awk '{
#         printf(">%s\n", $1);
#         seq = "";
#         for (i = 7; i <= NF; i += 2) {
#             seq = seq $i;
#         }
#         printf("%s\n", seq);
#     }' "$ped_file" > "$fasta_file"
#     echo "[INFO] PED 转换为 FASTA 成功！"
# }

# function rename_fasta_ids() {
#     local input_fasta="${OUT_PREFIX}.fasta"
#     local output_fasta="${OUT_PREFIX}_rename.fasta"

#     echo "[INFO] 正在调用 Python 脚本重命名 FASTA 序列..."
#     "$PYTHON_PATH" "$SCRIPT_PATH" \
#         --mapping "$MAPPING_FILE" \
#         --input_fasta "$input_fasta" \
#         --output_fasta "$output_fasta"
#     echo "[INFO] FASTA 序列重命名完成！"
# }

# # ================== 主流程执行 ==================

# run_plink
# ped_to_fasta
# # rename_fasta_ids
# rm "$OUT_PREFIX".ped
# rm "$OUT_PREFIX".map
# rm "$OUT_PREFIX".log
# rm "$OUT_PREFIX".nosex
# echo "[DONE] 所有步骤执行完毕。"
#!/bin/bash

# 退出条件：任何命令失败时立即退出
set -euo pipefail

# ================== 配置参数 ==================
# 存放所有 .vcf.gz 文件的目录
VCF_DIR="/mnt/f/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/ρ方法/2-从vcf/data"

# FASTA 输出根目录（脚本会自动按样本名创建子目录和前缀）
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/ρ方法/2-从vcf/data"

# Python 可执行文件路径
PYTHON_PATH="/home/luolintao/miniconda3/bin/python3"

# 重命名脚本
SCRIPT_PATH="/mnt/f/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/ρ方法/2-从vcf/script/6_vcf2fasta_通过plink.py"

# ID 与单倍群映射表
MAPPING_FILE="/home/luolintao/07_20K_CPGDP/1_单倍群分型/data/质量控制_ID_Hap.tsv"

# 并行线程数
THREADS=32

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

# ================== 主流程循环 ==================

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
done

echo "[ALL DONE] 所有 VCF 文件处理完毕。"

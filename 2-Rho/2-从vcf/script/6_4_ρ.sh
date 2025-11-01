#!/bin/bash
set -euo pipefail

# 并行线程数
THREADS=32

# 基础目录
BASE="/mnt/f/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/ρ方法/2-从vcf/data"

# 输入文件
INPUT="$BASE/merged_biallelic_complete_fixed_anno.vcf.gz"

# 仍然可选地生成一个仅含 SNP 的文件
NOINDEL="$BASE/merged_biallelic_complete_fixed_anno_SNP.vcf.gz"

# 下游各类派生文件 —— 统统直接基于 $INPUT
OUT_COMPLETE="$BASE/RHO_Complete_sequence.vcf.gz"
OUT_CONTROL="$BASE/RHO_Control_Region.vcf.gz"
OUT_HVSI="$BASE/RHO_HVS-I.vcf.gz"
OUT_HVSII="$BASE/RHO_HVS-II.vcf.gz"
OUT_HVSI_TRANS="$BASE/RHO_HVS-I_transitions.vcf.gz"

########################################
echo "[0] （可选）提取所有 SNP（去除 InDel）"
bcftools view \
  --threads $THREADS \
  -r chrM \
  -v snps \
  -Oz -o "$NOINDEL" \
  -i '(REF="A" | REF="C" | REF="G" | REF="T") & (ALT="A" | ALT="C" | ALT="G" | ALT="T")' \
  "$INPUT"
bcftools index --threads $THREADS "$NOINDEL"
########################################
echo "[1] 导出完整线粒体序列区段"
bcftools view \
  --threads $THREADS \
  -r chrM:1-16181,chrM:16184-16193,chrM:16195-16569 \
  -Oz -o "$OUT_COMPLETE" \
  "$NOINDEL"

echo "[2] 导出控制区（16024-576）"
bcftools view \
  --threads $THREADS \
  -r chrM:1-576,chrM:16024-16569 \
  -Oz -o "$OUT_CONTROL" \
  "$NOINDEL"

echo "[3] 导出 HVS-I（16051-16400）"
bcftools view \
  --threads $THREADS \
  -r chrM:16051-16181,chrM:16184-16193,chrM:16195-16400, \
  -Oz -o "$OUT_HVSI" \
  "$NOINDEL"

echo "[4] 导出 HVS-II（68-263）"
bcftools view \
  --threads $THREADS \
  -r chrM:68-263 \
  -Oz -o "$OUT_HVSII" \
  "$NOINDEL"

echo "[5] 导出 HVS-I transitions（16090-16365，且仅 transitions）"
# -------- [5] HVS-I transitions 修正版（方案 B） --------
bcftools view --threads $THREADS \
  -r chrM:16090-16365 -v snps \
  -Ou  "$NOINDEL" | \
bcftools view --threads $THREADS \
  -i '(REF="A" & ALT="G") | (REF="G" & ALT="A") | (REF="C" & ALT="T") | (REF="T" & ALT="C")' \
  -Oz -o "$OUT_HVSI_TRANS"



# 索引所有生成的 VCF
for f in "$OUT_COMPLETE" "$OUT_CONTROL" "$OUT_HVSI" "$OUT_HVSII" "$OUT_HVSI_TRANS"  ; do
  bcftools index --threads $THREADS -f "$f"
done

echo "全部处理完成。"

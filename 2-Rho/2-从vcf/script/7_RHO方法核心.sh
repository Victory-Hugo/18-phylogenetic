#!/usr/bin/env bash
set -euo pipefail

# —————————————— 配置区 ——————————————
BASE_DIR='/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/2-Rho/2-从vcf'
INPUT="$BASE_DIR/input/rho.csv"
TEMPDIR="$BASE_DIR/temp"
OUTDIR="$BASE_DIR/output/rho_number"
LOGFILE="$BASE_DIR/process.log"
PYTHON_BIN="/home/luolintao/miniconda3/envs/BigLin/bin/python3"
#todo 选择使用vcf计算还是使用fasta计算
#PY_SCRIPT="$BASE_DIR/script/7_RHO方法.py"
PY_SCRIPT="$BASE_DIR/script/7_RHO方法_fasta.py"
THREADS=24     # 传给 Python 脚本的 -t 参数
# —————————————— 结束配置 ——————————————

# 创建目录和清空日志
mkdir -p "$TEMPDIR" "$OUTDIR" "$(dirname "$LOGFILE")"
: > "$LOGFILE"

echo "[$(date '+%F %T')] === 开始拆分 ${INPUT} ===" | tee -a "$LOGFILE"
# 拆分 rho.csv，根据第二列（hap）分文件，并只写一次表头
{
  read header
  tail -n +2 | awk -F',' -v hdr="$header" -v outdir="$TEMPDIR" '
  {
    hap = $2
    gsub(/^[ \t]+|[ \t]+$/, "", hap)
    fname = outdir "/" hap ".csv"
    if (!(hap in seen)) {
      print hdr > fname
      seen[hap] = 1
    }
    print $0 >> fname
  }'
} < "$INPUT"

echo "[$(date '+%F %T')] === 拆分完成，临时目录：${TEMPDIR} ===" | tee -a "$LOGFILE"

# 遍历每个子文件，顺序调用 Python，跳过已存在的输出
for csv in "$TEMPDIR"/*.csv; do
  name="$(basename "$csv" .csv)"
  out="$OUTDIR/${name}.csv"

  if [[ -f "$out" ]]; then
    echo "[$(date '+%F %T')] 跳过 ${name}（已存在）" | tee -a "$LOGFILE"
    continue
  fi

  echo "[$(date '+%F %T')] 开始处理 ${name}" | tee -a "$LOGFILE"
  # 调用 Python，-t 参数留给脚本内部并行使用
  "$PYTHON_BIN" "$PY_SCRIPT" \
    "$csv" \
    "$BASE_DIR/data" \
    "$out" \
    -t "$THREADS" -v
  echo "[$(date '+%F %T')] 完成 ${name}" | tee -a "$LOGFILE"
done

echo "[$(date '+%F %T')] === 全部任务完成 ===" | tee -a "$LOGFILE"


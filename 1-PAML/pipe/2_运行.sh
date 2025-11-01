#!/bin/bash
set -euo pipefail

# --- OpenMP threads ---
export OMP_NUM_THREADS=18
export OMP_DYNAMIC=false
export OMP_PROC_BIND=true
export OMP_PLACES=cores
# 调试可开：export OMP_DISPLAY_ENV=true

# 路径
BASEML_PATH="/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/Linux/bin/baseml_LLST"
CTL_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/conf/baseml_mtDNA.ctl"
OUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/output"

# 准备输出目录
mkdir -p "$OUT_DIR"

# 拷贝 ctl 到输出目录并在该目录运行（保证相对路径一致）
cp "$CTL_FILE" "$OUT_DIR/"
cd "$OUT_DIR"

# 运行 baseml，使用当前目录下的 ctl
"${BASEML_PATH}" "baseml_mtDNA.ctl"

# 可选：清理临时/中间文件（-f 避免 set -e 下因不存在而报错）
# rm -f "lnf" "rst" "rst1" "rst2" "rub" "2base.t"
# 如果不想保留 ctl：
# rm -f "baseml.ctl"

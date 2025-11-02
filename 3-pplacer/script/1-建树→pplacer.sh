#!/usr/bin/env bash
set -euo pipefail

# Author: BigLin
# Requirements:
#   - FastTree (for --method fasttree)
#   - raxmlHPC-PTHREADS (for --method)
#   - taxit
#   - pplacer
#   - perl
#   - grep
#   - tput (optional, for colored output)

DATA_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/example/data"
QUERY_ALN="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/example/input/placer.aln.snp-sites.fasta"
REF_ALN="$DATA_DIR/example.aln.snp-sites.fasta"

FASTTREE_REF_PKG="$DATA_DIR/example.refpkg"
FASTTREE_OUTPUT_JPLACE="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/placer.aln.snp-sites.jplace"

RAXML_REF_PKG="$DATA_DIR/example.raxml.refpkg"
RAXML_OUTPUT_JPLACE="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/placer.aln.snp-sites.raxml.jplace"

# 固定执行策略：fasttree 或 raxml
METHOD="raxmlHPC-PTHREADS"
# 若使用 RAxML，可在此设置线程数（为空时自动检测）
THREADS_OVERRIDE=""

# pplacer 参数
PPLACER_THREADS=8
PPLACER_KEEP=5

METHOD="${METHOD,,}"
case "$METHOD" in
  fasttree)
    METHOD="fasttree"
    ;;
  raxml|raxmlhpc|raxmlhpc-pthreads|raxmlhpc_pthreads)
    METHOD="raxml"
    ;;
  *)
    echo "Unsupported METHOD value: $METHOD" >&2
    exit 64
    ;;
esac

if command -v tput >/dev/null 2>&1; then
  _bold="$(tput bold || true)"
  _sgr0="$(tput sgr0 || true)"
else
  _bold=""; _sgr0=""
fi
RED="${_bold}\033[31m"
GREEN="${_bold}\033[32m"
YELLOW="${_bold}\033[33m"
BLUE="${_bold}\033[34m"
RESET="\033[0m${_sgr0}"

cecho() { printf "%b%s%b\n" "${1}" "${2}" "${RESET}"; }

TOTAL_STEPS=6
CURRENT_STEP=0

progress_step() {
  local message="$1"
  (( ++CURRENT_STEP ))
  local done="$CURRENT_STEP"
  local total="$TOTAL_STEPS"
  if (( done > total )); then
    done="$total"
  fi
  local width=24
  local filled=$(( done * width / total ))
  local empty=$(( width - filled ))
  local bar=""
  if (( filled > 0 )); then
    bar="$(printf "%${filled}s" "" | tr ' ' '#')"
  fi
  if (( empty > 0 )); then
    bar+=$(printf "%${empty}s" "" | tr ' ' '-')
  fi
  cecho "${BLUE}" "[${bar}] ${done}/${total} ${message}"
}

progress_step "Checking dependencies"
need_cmds=(taxit pplacer perl grep)
for cmd in "${need_cmds[@]}"; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    cecho "${RED}" "[x] 缺少依赖：$cmd"
    exit 127
  fi
done

RAXML_CMD=""
RAXML_SUPPORTS_THREADS=0
if [[ "$METHOD" == "fasttree" ]]; then
  if ! command -v FastTree >/dev/null 2>&1; then
    cecho "${RED}" "[x] 缺少依赖：FastTree"
    exit 127
  fi
else
  if command -v raxmlHPC-PTHREADS >/dev/null 2>&1; then
    RAXML_CMD="raxmlHPC-PTHREADS"
    RAXML_SUPPORTS_THREADS=1
  elif command -v raxmlHPC >/dev/null 2>&1; then
    RAXML_CMD="raxmlHPC"
  else
    cecho "${RED}" "[x] 缺少依赖：raxmlHPC 或 raxmlHPC-PTHREADS"
    exit 127
  fi
  cecho "${GREEN}" "[✓] 使用 RAxML 命令：${RAXML_CMD}"
fi
cecho "${GREEN}" "[✓] 依赖检查通过。"

progress_step "Validating & sanitizing alignments"
for aln in "$REF_ALN" "$QUERY_ALN"; do
  if [[ ! -f $aln ]]; then
    cecho "${RED}" "[x] 未找到对齐文件：$aln"
    exit 2
  fi
done

sanitize_alignment() {
  local file="$1"
  if grep -q '\.' "$file"; then
    cecho "${YELLOW}" "[*] 将 $file 中的 '.' 替换为 '-'"
    perl -i -pe 's/\./-/g unless /^>/' "$file"
  fi
}

cecho "${BLUE}" "[*] Working directory: $DATA_DIR"
cd "$DATA_DIR"

sanitize_alignment "$REF_ALN"
sanitize_alignment "$QUERY_ALN"

BASE="$(basename "$REF_ALN")"
BASE_NOEXT="${BASE%.fasta}"

progress_step "Estimating reference tree (${METHOD})"
if [[ "$METHOD" == "fasttree" ]]; then
  TREE_FILE="${BASE_NOEXT}.tree"
  STATS_FILE="${BASE_NOEXT}.log"
  REF_PKG="$FASTTREE_REF_PKG"
  OUTPUT_JPLACE="$FASTTREE_OUTPUT_JPLACE"

  rm -f "$TREE_FILE" "$STATS_FILE"
  cecho "${YELLOW}" "[*] 使用 FastTree 估计参考树..."
  FastTree -nt -gtr -gamma -fastest -log "$STATS_FILE" < "$REF_ALN" > "$TREE_FILE"
else
  RUN_NAME="ref_tree"
  TREE_FILE="RAxML_bestTree.${RUN_NAME}"
  STATS_FILE="RAxML_info.${RUN_NAME}"
  REF_PKG="$RAXML_REF_PKG"
  OUTPUT_JPLACE="$RAXML_OUTPUT_JPLACE"

  rm -f RAxML_*."$RUN_NAME" 2>/dev/null || true
  cecho "${YELLOW}" "[*] 使用 ${RAXML_CMD} 估计参考树..."
  THREADS_VALUE="$THREADS_OVERRIDE"
  if [[ -z "${THREADS_VALUE}" && $RAXML_SUPPORTS_THREADS -eq 1 ]]; then
    if command -v nproc >/dev/null 2>&1; then
      THREADS_VALUE="$(nproc)"
    else
      THREADS_VALUE=2
    fi
  fi
  if [[ $RAXML_SUPPORTS_THREADS -eq 1 ]]; then
    cecho "${BLUE}" "[*] RAxML 线程数：${THREADS_VALUE}"
    "$RAXML_CMD" -T "${THREADS_VALUE}" -m GTRGAMMA -p 12345 -s "$REF_ALN" -n "$RUN_NAME"
  else
    "$RAXML_CMD" -m GTRGAMMA -p 12345 -s "$REF_ALN" -n "$RUN_NAME"
  fi
fi

if [[ ! -f $TREE_FILE || ! -f $STATS_FILE ]]; then
  cecho "${RED}" "[x] 未生成预期的树或统计文件：$TREE_FILE / $STATS_FILE"
  exit 4
fi
cecho "${GREEN}" "[✓] 树构建完成：TREE=$TREE_FILE  INFO=$STATS_FILE"

progress_step "Building reference package"
STATS_TYPE="FastTree"
PKG_LABEL="example_ref"
if [[ "$METHOD" == "raxml" ]]; then
  STATS_TYPE="RAxML"
  PKG_LABEL="example_raxml"
fi

taxit create \
  --clobber \
  -l "$PKG_LABEL" \
  -P "$REF_PKG" \
  --aln-fasta "$REF_ALN" \
  --tree-file "$TREE_FILE" \
  --tree-stats "$STATS_FILE" \
  --stats-type "$STATS_TYPE"
cecho "${GREEN}" "[✓] refpkg 构建完成：$REF_PKG"

progress_step "Running pplacer"
pplacer \
  -j "${PPLACER_THREADS}" \
  --keep-at-most "${PPLACER_KEEP}" \
  -p \
  -c "$REF_PKG" \
  "$QUERY_ALN"

progress_step "Finalizing outputs"
DEFAULT_JPLACE="$(basename "$QUERY_ALN" .fasta).jplace"
if [[ -f "$DEFAULT_JPLACE" ]]; then
  mkdir -p "$(dirname "$OUTPUT_JPLACE")"
  mv -f "$DEFAULT_JPLACE" "$OUTPUT_JPLACE"
  cecho "${GREEN}" "[✓] 输出文件已移动到：$OUTPUT_JPLACE"
else
  cecho "${RED}" "[x] 未找到默认输出文件 (${DEFAULT_JPLACE})，请检查 pplacer 运行是否成功。"
  exit 5
fi

cecho "${GREEN}" "[✓] 全部完成，方法：$METHOD"

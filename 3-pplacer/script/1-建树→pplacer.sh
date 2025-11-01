#!/usr/bin/env bash
set -euo pipefail

# =================== 固定路径 ===================
DATA_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/example/data"
QUERY_ALN="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/example/input/placer.aln.snp-sites.fasta"
REF_ALN="$DATA_DIR/example.aln.snp-sites.fasta"
REF_PKG="$DATA_DIR/example.refpkg"

# ✅ 输出路径（保持）
OUTPUT_JPLACE="/mnt/f/OneDrive/文档（科研）/脚本/Download/18-phylogenetic/3-pplacer/output/placer.aln.snp-sites.jplace"

# =================== 彩色打印 ===================
if command -v tput >/dev/null 2>&1; then
  _bold="$(tput bold || true)"
  _sgr0="$(tput sgr0 || true)"
else
  _bold=""; _sgr0=""
fi
RED="${_bold}\033[31m"; GREEN="${_bold}\033[32m"; YELLOW="${_bold}\033[33m"; BLUE="${_bold}\033[34m"; RESET="\033[0m${_sgr0}"
cecho() { printf "%b%s%b\n" "${1}" "${2}" "${RESET}"; }

# =================== 依赖检查 ===================
need_cmds=(FastTree taxit pplacer perl grep)
for cmd in "${need_cmds[@]}"; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    cecho "${RED}" "[x] 缺少依赖：$cmd"
    exit 127
  fi
done
cecho "${GREEN}" "[✓] 依赖检查通过。"

# =================== 文件验证 ===================
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

# =================== 基于输入派生默认文件名 ===================
BASE="$(basename "$REF_ALN")"
BASE_NOEXT="${BASE%.fasta}"            # example.aln.snp-sites
TREE_FILE="${BASE_NOEXT}.tree"         # example.aln.snp-sites.tree
STATS_FILE="${BASE_NOEXT}.log"         # example.aln.snp-sites.log（FastTree 日志）

# =================== 建树（不再强行命名为 refpkg.bestTree/refpkg.info） ===================
cecho "${YELLOW}" "[*] 使用 FastTree 估计参考树..."
rm -f "$TREE_FILE" "$STATS_FILE"
FastTree -nt -gtr -gamma -fastest -log "$STATS_FILE" < "$REF_ALN" > "$TREE_FILE"
cecho "${GREEN}" "[✓] FastTree 完成：TREE=$TREE_FILE  LOG=$STATS_FILE"

# =================== 构建 refpkg ===================
cecho "${YELLOW}" "[*] 构建 reference package..."
taxit create \
  --clobber \
  -l example_ref \
  -P "$REF_PKG" \
  --aln-fasta "$REF_ALN" \
  --tree-file "$TREE_FILE" \
  --tree-stats "$STATS_FILE" \
  --stats-type FastTree
cecho "${GREEN}" "[✓] refpkg 构建完成：$REF_PKG"

# =================== pplacer 放置 ===================
cecho "${YELLOW}" "[*] 运行 pplacer..."
pplacer -p --keep-at-most 5 -c "$REF_PKG" "$QUERY_ALN"

# ✅ 自动将结果移动到指定输出位置
DEFAULT_JPLACE="$(basename "$QUERY_ALN" .fasta).jplace"
if [[ -f "$DEFAULT_JPLACE" ]]; then
  mkdir -p "$(dirname "$OUTPUT_JPLACE")"
  mv -f "$DEFAULT_JPLACE" "$OUTPUT_JPLACE"
  cecho "${GREEN}" "[✓] 输出文件已移动到：$OUTPUT_JPLACE"
else
  cecho "${RED}" "[x] 未找到默认输出文件 (${DEFAULT_JPLACE})，请检查 pplacer 运行是否成功。"
  exit 3
fi

cecho "${BLUE}" "[*] Done."

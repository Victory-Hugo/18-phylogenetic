#!/usr/bin/env bash
# 时间树累计节点山脊图：准备输入 → 节点年龄 → 山脊密度 → 出图
# 用法：bash pipe/1-pipe.sh [conf/1-conf.yaml]
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

CONF="${1:-conf/1-conf.yaml}"
[[ -f "$CONF" ]] || { echo "找不到配置文件：$CONF" >&2; exit 1; }
mkdir -p log

# conf 的标量字段展开成大写下划线的 shell 变量；groupings 等嵌套结构不展开，由 Python 直接读 conf
eval "$(python3 - "$CONF" <<'PY'
import sys, yaml
conf = yaml.safe_load(open(sys.argv[1], encoding="utf-8")) or {}
for section in ("paths", "analysis", "plot", "tools"):
    for key, value in (conf.get(section) or {}).items():
        if isinstance(value, (dict, list)):
            if isinstance(value, list):
                print(f"export {key.upper()}={','.join(map(str, value))!r}")
            continue
        print(f"export {key.upper()}={str(value)!r}")
PY
)"

# conf 登记的解释器路径只是本机默认值，不可执行时回落到 PATH 上的同名命令
resolve_bin() {
  local var="$1" fallback="$2" value="${!1:-}"
  if [[ -z "$value" || ! -x "$value" ]]; then
    value="$(command -v "$fallback")" || { echo "找不到可用的 $fallback" >&2; exit 1; }
    echo "提示：$var 不可用，改用 $value" >&2
  fi
  printf -v "$var" '%s' "$value"
}
resolve_bin PYTHON_BIN python3
resolve_bin R_BIN Rscript

RESULT_1="output/result/1-prepare_input"
RESULT_2="output/result/2-node_ages"
RESULT_3="output/result/3-ridge_density"

#* =====1 准备输入=====
"$PYTHON_BIN" python/1-1-prepare_input.py \
  --conf "$CONF" \
  --output-result  "$RESULT_1" \
  2>&1 | tee log/1-prepare_input.log

#* =====2 节点年龄=====
"$PYTHON_BIN" python/2-1-node_ages.py \
  --conf "$CONF" \
  --sample-groups "$RESULT_1/sample-groups.tsv" \
  --output-result "$RESULT_2" \
  --temp          "temp/2-node_ages" \
  2>&1 | tee log/2-node_ages.log

#* =====3 山脊密度=====
"$PYTHON_BIN" python/3-1-ridge_density.py \
  --conf "$CONF" \
  --node-ages          "$RESULT_2/⭐节点年龄长表.tsv.gz" \
  --rarefied-node-ages "temp/2-node_ages/rarefied-node-ages.tsv.gz" \
  --sample-counts      "$RESULT_1/sample-counts.tsv" \
  --group-design       "$RESULT_1/group-design.tsv" \
  --output-result      "$RESULT_3" \
  2>&1 | tee log/3-ridge_density.log

#* =====4 出图=====
"$R_BIN" R/4-1-plot_ridges.R \
  --ridge-density "$RESULT_3/⭐山脊密度长表.tsv" \
  --group-design  "$RESULT_1/group-design.tsv" \
  --output-figure "output/figure/4-plot_ridges" \
  --fill-mode              "$FILL_MODE" \
  --colour-ramp            "$COLOUR_RAMP" \
  --ridge-scale            "$RIDGE_SCALE" \
  --ridge-alpha            "$RIDGE_ALPHA" \
  --outline-width          "$OUTLINE_WIDTH" \
  --min-height-fraction    "$MIN_HEIGHT_FRACTION" \
  --baseline-colour        "$BASELINE_COLOUR" \
  --outline-colour         "$OUTLINE_COLOUR" \
  --reference-colour       "$REFERENCE_COLOUR" \
  --interval-colour        "$INTERVAL_COLOUR" \
  --interval-alpha         "$INTERVAL_ALPHA" \
  --reference-lines-recent "$REFERENCE_LINES_RECENT" \
  --reference-lines-full   "$REFERENCE_LINES_FULL" \
  --font-size              "$FONT_SIZE" \
  --panel-mm               "$PANEL_MM" \
  --label-mm               "$LABEL_MM" \
  --min-row-mm             "$MIN_ROW_MM" \
  2>&1 | tee log/4-plot_ridges.log

#* =====字体验收=====
if command -v pdffonts >/dev/null; then
  echo "PDF 字体检查："
  for pdf in output/figure/4-plot_ridges/*.pdf; do
    printf '  %-58s %s\n' "$(basename "$pdf")" \
      "$(pdffonts "$pdf" | awk 'NR>2 && $1 != "" {print $1}' | sort -u | paste -sd, -)"
  done
else
  echo "未安装 pdffonts，跳过 PDF 字体检查" >&2
fi

#!/usr/bin/env bash
# run.sh — 唯一入口：特征提取 → 聚类 → 四类图
#
# 用法：
#   bash pipe/run.sh                      # 用 conf/config.yaml
#   NE_CONFIG=conf/other.yaml bash pipe/run.sh
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

CONFIG="${NE_CONFIG:-conf/config.yaml}"
[[ -f "$CONFIG" ]] || { echo "找不到配置文件 $CONFIG" >&2; exit 1; }


# 只从配置里取脚本自己要用的几个标量，其余全部交给 Python 解析
yaml_scalar() {
    awk -v target="$1" '
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        /^[^[:space:]#][^:]*:/ { section = $0; sub(/:.*/, "", section); next }
        /^[[:space:]]{2}[A-Za-z0-9_]+:/ {
            line = $0; sub(/^[[:space:]]+/, "", line)
            key = line; sub(/:.*/, "", key)
            val = line; sub(/^[^:]+:[[:space:]]*/, "", val)
            sub(/[[:space:]]+#.*$/, "", val)
            gsub(/^["'\''"]+|["'\''"]+$/, "", val)
            if (section "." key == target) { print val; exit }
        }' "$CONFIG"
}

# 解释器：环境变量 > 配置 > 当前 conda 环境 > PATH。裸 python3 在很多机器上指向系统
# 自带的那个（没有 numpy），激活了 conda 环境也未必排在 PATH 最前面
PYTHON_BIN="${NE_PYTHON:-$(yaml_scalar tools.python_bin)}"
if [[ -z "$PYTHON_BIN" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python3" ]]; then
        PYTHON_BIN="$CONDA_PREFIX/bin/python3"
    else
        PYTHON_BIN="python3"
    fi
fi
RSCRIPT_BIN="${NE_RSCRIPT:-$(yaml_scalar tools.rscript_bin)}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

# 依赖不全时当场说清楚该怎么办，而不是让 Python 在第 24 行抛 ModuleNotFoundError
MISSING="$("$PYTHON_BIN" - <<'PYCHECK' 2>/dev/null || echo "python-unusable"
import importlib.util
need = {"numpy": "numpy", "scipy": "scipy", "sklearn": "scikit-learn",
        "pandas": "pandas", "yaml": "pyyaml"}
print(" ".join(pkg for mod, pkg in need.items()
               if importlib.util.find_spec(mod) is None))
PYCHECK
)"
if [[ -n "$MISSING" ]]; then
    echo "[run] 解释器 $PYTHON_BIN 不可用或缺少依赖：$MISSING" >&2
    echo "      装依赖：conda install -c conda-forge $MISSING" >&2
    echo "      或指定解释器：NE_PYTHON=/path/to/python3 bash pipe/run.sh" >&2
    echo "      也可以把路径写进 $CONFIG 的 tools.python_bin" >&2
    exit 1
fi
command -v "$RSCRIPT_BIN" >/dev/null || { echo "[run] 找不到 $RSCRIPT_BIN" >&2; exit 1; }
echo "[run] Python: $PYTHON_BIN"

OUT_DIR="$(yaml_scalar io.output_dir)";        OUT_DIR="${OUT_DIR:-output}"
DEBUG_MODE="$(yaml_scalar runtime.debug_mode)"; DEBUG_MODE="${DEBUG_MODE:-false}"
SCATTER="$(yaml_scalar figure.scatter_method)"; SCATTER="${SCATTER:-pca}"
THRESHOLD="$(yaml_scalar figure.chord_threshold)"; THRESHOLD="${THRESHOLD:-0.5}"
HEAT_PAL="$(yaml_scalar figure.heatmap_palette)"; HEAT_PAL="${HEAT_PAL:-diverging}"
SCORE_PAL="$(yaml_scalar figure.score_palette)"; SCORE_PAL="${SCORE_PAL:-diverging}"
ELLIPSE="$(yaml_scalar figure.scatter_ellipse)"; ELLIPSE="${ELLIPSE:-true}"
N_COL="$(yaml_scalar figure.panels_per_row)"; N_COL="${N_COL:-3}"

mkdir -p "$OUT_DIR"
SUCCESS_LOG="$OUT_DIR/success.log"
FAIL_LOG="$OUT_DIR/fail.log"
touch "$SUCCESS_LOG" "$FAIL_LOG"
mark_ok()   { echo "$1" >> "$SUCCESS_LOG"; }
mark_fail() { printf '%s\t%s\n' "$(date '+%F %T')" "$1" >> "$FAIL_LOG"; }

run_step() {
    local name="$1"; shift
    if "$@"; then mark_ok "$name"; else mark_fail "$name"; exit 1; fi
}

run_step "1-extract_features" "$PYTHON_BIN" python/1-extract_features.py \
    --config "$CONFIG" --project-root "$PROJECT_ROOT"
run_step "2-cluster"          "$PYTHON_BIN" python/2-cluster.py \
    --config "$CONFIG" --project-root "$PROJECT_ROOT"

FEATURE_TABLE="$OUT_DIR/⭐1-Feature-Summary.tsv"
CLUSTER_TABLE="$OUT_DIR/⭐2-Cluster-And-Similarity.tsv"

run_step "1-plot_period_heatmap" "$RSCRIPT_BIN" R/1-plot_period_heatmap.R \
    "$FEATURE_TABLE" "$OUT_DIR" "$HEAT_PAL" "$N_COL"
run_step "2-plot_all_time_cluster_heatmap" "$RSCRIPT_BIN" R/2-plot_all_time_cluster_heatmap.R \
    "$CLUSTER_TABLE" "$OUT_DIR" "$SCORE_PAL"
run_step "3-plot_cluster_scatter" "$RSCRIPT_BIN" R/3-plot_cluster_scatter.R \
    "$CLUSTER_TABLE" "$OUT_DIR" "$SCATTER" "$ELLIPSE" "$N_COL"
run_step "4-plot_chord_diagram" "$RSCRIPT_BIN" R/4-plot_chord_diagram.R \
    "$CLUSTER_TABLE" "$OUT_DIR" "$THRESHOLD" "$N_COL"

# 曲线缓存只是两个 Python 步骤之间的接力，默认不留
if [[ "${DEBUG_MODE,,}" != "true" ]]; then
    rm -rf "$OUT_DIR/.cache"
fi
rm -f Rplots.pdf

echo "[run] 完成，结果见 $OUT_DIR"

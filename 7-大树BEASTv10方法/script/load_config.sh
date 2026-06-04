#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:?Usage: source script/load_config.sh conf/<step>.yaml}"

yaml_value() {
    local key="$1"
    awk -v target="$key" '
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        /^[^[:space:]][^:]*:/ {
            gsub(":", "", $1)
            section=$1
            next
        }
        /^[[:space:]]+[A-Za-z0-9_]+:/ {
            line=$0
            sub(/^[[:space:]]+/, "", line)
            split(line, parts, ":")
            name=parts[1]
            value=line
            sub(/^[^:]+:[[:space:]]*/, "", value)
            gsub(/^["'\''"]|["'\''"]$/, "", value)
            path=section "." name
            if (path == target) {
                print value
                exit
            }
        }
    ' "$CONFIG_FILE"
}

# 解析 calibration.nodes 列表（任意数量、支持嵌套），输出紧凑串：
#   name:age:stdev:is_root  以 ';' 分隔多个节点；is_root 用 1/0。
# 下游 python 模块通过单个 --calibration-nodes 参数接收并解析（模块不读 YAML）。
yaml_calibration_nodes() {
    awk '
        function flush() {
            if (name != "") {
                out = out (out == "" ? "" : ";") name ":" age ":" stdev ":" isroot
            }
            name=""; age=""; stdev=""; isroot="0"
        }
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        /^[^[:space:]]/ {
            flush()
            in_cal = ($1 == "calibration:")
            in_nodes = 0
            next
        }
        in_cal && /^[[:space:]]{2}nodes:[[:space:]]*$/ { in_nodes=1; next }
        in_nodes && /^[[:space:]]{4}-[[:space:]]*name:/ {
            flush()
            line=$0
            sub(/^[[:space:]]*-[[:space:]]*name:[[:space:]]*/, "", line)
            gsub(/"/, "", line)
            name=line
            next
        }
        in_nodes && /^[[:space:]]{6}[A-Za-z_]+:/ {
            line=$0
            sub(/^[[:space:]]+/, "", line)
            key=line; sub(/:.*/, "", key)
            val=line; sub(/^[^:]+:[[:space:]]*/, "", val)
            gsub(/"/, "", val)
            if (key == "age") age=val
            else if (key == "stdev") stdev=val
            else if (key == "is_root") { if (val == "true" || val == "1" || val == "True") isroot="1" }
            next
        }
        END { flush(); print out }
    ' "$CONFIG_FILE"
}

BASE_DIR="$(yaml_value project.base_dir)"
FASTA="$(yaml_value paths.fasta)"
ML_TREE="$(yaml_value paths.ml_tree)"
ID_HAP="$(yaml_value paths.id_hap)"
PHYLOTREE_JSON="$(yaml_value paths.phylotree_json)"
CALIBRATION_MEMBERSHIP="$(yaml_value paths.calibration_membership)"
STARTING_TREE="$(yaml_value paths.starting_tree)"
INPUT_VALIDATION="$(yaml_value paths.input_validation)"
TEMPLATE_XML="$(yaml_value paths.template_xml)"
OUTPUT_TABLE="$(yaml_value paths.output_table)"
OUTPUT_REPORT="$(yaml_value paths.output_report)"
TEMP_DIR="$(yaml_value paths.temp)"
LOG_DIR="$(yaml_value paths.log)"
PYTHON_BIN="$(yaml_value tools.python_bin)"
IQTREE_BIN="$(yaml_value tools.iqtree_bin)"
BEAST_BIN="$(yaml_value tools.beast_bin)"
OVERWRITE="$(yaml_value runtime.overwrite)"
IQTREE_MODEL="$(yaml_value dating.model)"
IQTREE_THREADS="$(yaml_value dating.threads)"
IQTREE_OUTGROUP="$(yaml_value dating.outgroup)"
IQTREE_DATE_OPTIONS="$(yaml_value dating.date_options)"
BEAST_OUTPUT_PREFIX="$(yaml_value beast.output_prefix)"
BEAST_SCALE="$(yaml_value beast.scale)"
BEAST_CHAIN_LENGTH="$(yaml_value beast.chain_length)"
BEAST_LOG_EVERY="$(yaml_value beast.log_every)"
BEAST_TREE_LOG_EVERY="$(yaml_value beast.tree_log_every)"
BEAST_CLOCK_RATE="$(yaml_value beast.clock_rate)"
BEAST_SKYGRID_POINTS="$(yaml_value beast.skygrid_points)"
BEAST_SKYGRID_INIT_LOGPOP="$(yaml_value beast.skygrid_init_logpop)"
CALIBRATION_NODES="$(yaml_calibration_nodes)"

export BASE_DIR FASTA ML_TREE ID_HAP PHYLOTREE_JSON CALIBRATION_MEMBERSHIP STARTING_TREE INPUT_VALIDATION TEMPLATE_XML OUTPUT_TABLE OUTPUT_REPORT TEMP_DIR LOG_DIR
export PYTHON_BIN IQTREE_BIN BEAST_BIN OVERWRITE
export IQTREE_MODEL IQTREE_THREADS IQTREE_OUTGROUP IQTREE_DATE_OPTIONS
export BEAST_OUTPUT_PREFIX BEAST_SCALE BEAST_CHAIN_LENGTH BEAST_LOG_EVERY BEAST_TREE_LOG_EVERY BEAST_CLOCK_RATE BEAST_SKYGRID_POINTS BEAST_SKYGRID_INIT_LOGPOP
export CALIBRATION_NODES

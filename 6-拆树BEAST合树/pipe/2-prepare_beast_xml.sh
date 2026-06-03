#!/usr/bin/env bash

# 2-prepare_beast_xml.sh
# 只准备 BEAST XML/FASTA/manifest，不运行 BEAST MCMC。

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo

PATH_ROOT="$PROJECT_ROOT"
DEFAULT_CONFIG_PATH="$PROJECT_ROOT/conf/2-beast.yaml"
CONFIG_PATH="$DEFAULT_CONFIG_PATH"
BOOTSTRAP_PYTHON="python3"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        *)
            ui_error "参数错误 | Unknown argument: $1"
            exit 1
            ;;
    esac
done

resolve_path() {
    local value="$1"
    if [[ -z "$value" || "$value" == "null" ]]; then
        printf '%s\n' "$value"
        return
    fi
    if [[ "$value" == /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$PATH_ROOT/$value"
    fi
}

resolve_command_or_path() {
    local value="$1"
    if [[ "$value" == */* || "$value" == .* ]]; then
        resolve_path "$value"
    else
        printf '%s\n' "$value"
    fi
}

if [[ ! -f "$CONFIG_PATH" ]]; then
    ui_error "配置缺失 | Config file not found: $CONFIG_PATH"
    exit 1
fi

CONFIG_LOADER="$PROJECT_ROOT/python/config_loader.py"
config_get() { "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1"; }
config_get_default() { "$BOOTSTRAP_PYTHON" "$CONFIG_LOADER" --config "$CONFIG_PATH" --key "$1" --default "$2"; }

CONFIG_PROJECT_ROOT=$(config_get_default projectpath "$PROJECT_ROOT")
PATH_ROOT=$(resolve_path "$CONFIG_PROJECT_ROOT")

PYTHON_BIN=$(resolve_command_or_path "$(config_get tools.python)")
bash "$PROJECT_ROOT/script/check_env.sh" --python "$PYTHON_BIN"

BEAST_EXECUTION_MODE=$(config_get_default beast.execution_mode real)
SPLIT_OUTPUT_DIR=$(resolve_path "$(config_get paths.split_output_dir)")
INPUT_FASTA=$(resolve_path "$(config_get paths.input_fasta)")
TMP_DIR=$(resolve_path "$(config_get paths.tmp_dir)")
LOG_LEVEL=$(config_get runtime.log_level)

XML_TEMPLATE=$(resolve_path "$(config_get beast.xml_template)")
BEAST_OUTPUT_DIR=$(resolve_path "$(config_get beast.output_dir)")
BEAST_CLEAN_OUTPUT_DIR=$(config_get beast.clean_output_dir)
SEQ_ID_STRATEGY=$(config_get beast.seq_id_strategy)
CHAIN_LENGTH=$(config_get beast.chain_length)
LOG_EVERY=$(config_get beast.log_every)
CLOCK_RATE=$(config_get_default beast.clock_rate 1.0)
BACKBONE_ONLY_ENABLED=$(config_get beast.backbone_only_enabled)
BACKBONE_JOB_DIR=$(resolve_path "$(config_get beast.backbone_job_dir)")
BACKBONE_ANALYSIS_DIR=$(resolve_path "$(config_get beast.backbone_analysis_dir)")
CALIBRATION_ENABLED=$(config_get_default beast.calibration.enabled false)
CALIBRATION_TEMPLATE=$(resolve_path "$(config_get_default beast.calibration.template null)")
CALIBRATION_ID_HAP_TSV=$(resolve_path "$(config_get_default beast.calibration.id_hap_tsv null)")
CALIBRATION_PHYLOTREE_JSON=$(resolve_path "$(config_get_default beast.calibration.phylotree_json null)")
CALIBRATION_WARM_START=$(config_get_default beast.calibration.warm_start true)
CALIBRATION_POPSIZE_START=$(config_get_default beast.calibration.popsize_start 100000.0)
CALIBRATION_ROOT_AGE=$(config_get_default beast.calibration.nodes.root.age 180000)
CALIBRATION_ROOT_STDEV=$(config_get_default beast.calibration.nodes.root.stdev 20000)
CALIBRATION_M_AGE=$(config_get_default beast.calibration.nodes.M.age 50000)
CALIBRATION_M_STDEV=$(config_get_default beast.calibration.nodes.M.stdev 5000)
CALIBRATION_N_AGE=$(config_get_default beast.calibration.nodes.N.age 57000)
CALIBRATION_N_STDEV=$(config_get_default beast.calibration.nodes.N.stdev 5000)
ULTRAMETRIC_TOLERANCE=$(config_get beast.ultrametric_tolerance)
MIN_BRANCH_LENGTH=$(config_get_default beast.min_branch_length 1e-8)

BEAST_SUBTREE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/beast_subtree_summary.tsv"
BEAST_TREE_DIR="$SPLIT_OUTPUT_DIR/beast_subtrees"
BEAST_RESUME_SUCCESS_LOG="$BEAST_OUTPUT_DIR/log/success.log"
BACKBONE_TREE="$SPLIT_OUTPUT_DIR/backbone_tree.nwk"
BACKBONE_SUMMARY_TSV="$SPLIT_OUTPUT_DIR/backbone_summary.tsv"

for required in \
    "$BEAST_SUBTREE_SUMMARY_TSV" \
    "$BACKBONE_TREE" \
    "$BACKBONE_SUMMARY_TSV" \
    "$INPUT_FASTA" \
    "$XML_TEMPLATE"; do
    if [[ ! -f "$required" ]]; then
        ui_error "输入缺失 | Required BEAST prepare input not found: $required"
        exit 1
    fi
done

if [[ ! -d "$BEAST_TREE_DIR" ]]; then
    ui_error "输入缺失 | Subtree directory not found: $BEAST_TREE_DIR"
    exit 1
fi

if [[ "$CALIBRATION_ENABLED" == "True" || "$CALIBRATION_ENABLED" == "true" ]]; then
    for required in "$CALIBRATION_TEMPLATE" "$CALIBRATION_ID_HAP_TSV" "$CALIBRATION_PHYLOTREE_JSON"; do
        if [[ ! -f "$required" ]]; then
            ui_error "输入缺失 | Required calibration input not found: $required"
            exit 1
        fi
    done
fi

ui_section "BEAST XML preparation" "Step 2A/5 | 生成 BEAST XML/FASTA/manifest"
ui_kv "Config" "$CONFIG_PATH"
ui_kv "Mode" "$BEAST_EXECUTION_MODE"
ui_kv "Split dir" "$SPLIT_OUTPUT_DIR"
ui_kv "Output dir" "$BEAST_OUTPUT_DIR"
ui_kv "Chain length" "$CHAIN_LENGTH"
ui_kv "Calibration XML" "$CALIBRATION_ENABLED"

if [[ "$BEAST_EXECUTION_MODE" == "fake" ]]; then
    ui_error "模式不匹配 | execution_mode=fake does not prepare XML. Set beast.execution_mode=real."
    exit 1
fi

mkdir -p "$BEAST_OUTPUT_DIR" "$TMP_DIR"

if [[ "$BEAST_CLEAN_OUTPUT_DIR" == "True" || "$BEAST_CLEAN_OUTPUT_DIR" == "true" ]]; then
    if [[ -s "$BEAST_RESUME_SUCCESS_LOG" ]]; then
        ui_warn "检测到断点续跑日志 | Resume logs detected at $BEAST_RESUME_SUCCESS_LOG; skip destructive cleanup."
    else
        rm -rf "$BEAST_OUTPUT_DIR/jobs"
        rm -rf "$BEAST_OUTPUT_DIR/analysis_trees"
        rm -rf "$BEAST_OUTPUT_DIR/log"
        rm -rf "$BACKBONE_JOB_DIR"
        rm -rf "$BACKBONE_ANALYSIS_DIR"
        rm -f "$BEAST_OUTPUT_DIR/beast_job_manifest.tsv"
        rm -f "$BEAST_OUTPUT_DIR/beast_run_status.tsv"
        rm -f "$BEAST_OUTPUT_DIR/beast_parse_status.tsv"
        rm -f "$BEAST_OUTPUT_DIR/fake_beast_manifest.tsv"
        rm -f "$BEAST_OUTPUT_DIR/prepare_beast_inputs.log"
        rm -f "$BEAST_OUTPUT_DIR/run_beast_jobs.log"
        rm -f "$BEAST_OUTPUT_DIR/parse_beast_outputs.log"
        rm -f "$BEAST_OUTPUT_DIR/fake_beast_outputs.log"
    fi
fi

CALIBRATION_PREPARE_ARGS=()
if [[ "$CALIBRATION_ENABLED" == "True" || "$CALIBRATION_ENABLED" == "true" ]]; then
    CALIBRATION_PREPARE_ARGS=(
        --calibration-enabled true
        --calibrated-xml-template "$CALIBRATION_TEMPLATE"
        --calibration-id-hap-tsv "$CALIBRATION_ID_HAP_TSV"
        --calibration-phylotree-json "$CALIBRATION_PHYLOTREE_JSON"
        --calibration-warm-start "$CALIBRATION_WARM_START"
        --calibration-popsize-start "$CALIBRATION_POPSIZE_START"
        --calibration-root-age "$CALIBRATION_ROOT_AGE"
        --calibration-root-stdev "$CALIBRATION_ROOT_STDEV"
        --calibration-m-age "$CALIBRATION_M_AGE"
        --calibration-m-stdev "$CALIBRATION_M_STDEV"
        --calibration-n-age "$CALIBRATION_N_AGE"
        --calibration-n-stdev "$CALIBRATION_N_STDEV"
    )
fi

ui_stage_start "Stage 1/2" "prepare_beast_inputs | 生成子树 XML"
"$PYTHON_BIN" "$PROJECT_ROOT/python/prepare_beast_inputs.py" \
    --beast-summary-tsv "$BEAST_SUBTREE_SUMMARY_TSV" \
    --beast-tree-dir "$BEAST_TREE_DIR" \
    --input-fasta "$INPUT_FASTA" \
    --xml-template "$XML_TEMPLATE" \
    --output-dir "$BEAST_OUTPUT_DIR" \
    --seq-id-strategy "$SEQ_ID_STRATEGY" \
    --chain-length "$CHAIN_LENGTH" \
    --log-every "$LOG_EVERY" \
    --clock-rate "$CLOCK_RATE" \
    --min-branch-length "$MIN_BRANCH_LENGTH" \
    --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE" \
    --log-level "$LOG_LEVEL" \
    "${CALIBRATION_PREPARE_ARGS[@]}"
ui_stage_end "Stage 1/2" "prepare_beast_inputs | 生成子树 XML"

if [[ "$BACKBONE_ONLY_ENABLED" == "True" || "$BACKBONE_ONLY_ENABLED" == "true" ]]; then
    ui_stage_start "Stage 2/2" "prepare_backbone_beast_input | 生成 backbone XML"
    "$PYTHON_BIN" "$PROJECT_ROOT/python/prepare_backbone_beast_input.py" \
        --backbone-tree "$BACKBONE_TREE" \
        --backbone-summary-tsv "$BACKBONE_SUMMARY_TSV" \
        --beast-summary-tsv "$BEAST_SUBTREE_SUMMARY_TSV" \
        --input-fasta "$INPUT_FASTA" \
        --xml-template "$XML_TEMPLATE" \
        --job-dir "$BACKBONE_JOB_DIR" \
        --analysis-dir "$BACKBONE_ANALYSIS_DIR" \
        --seq-id-strategy "$SEQ_ID_STRATEGY" \
        --chain-length "$CHAIN_LENGTH" \
        --log-every "$LOG_EVERY" \
        --clock-rate "$CLOCK_RATE" \
        --min-branch-length "$MIN_BRANCH_LENGTH" \
        --ultrametric-tolerance "$ULTRAMETRIC_TOLERANCE" \
        --log-level "$LOG_LEVEL" \
        "${CALIBRATION_PREPARE_ARGS[@]}"
    ui_stage_end "Stage 2/2" "prepare_backbone_beast_input | 生成 backbone XML"
fi

ui_ok "Completed | BEAST XML preparation finished"

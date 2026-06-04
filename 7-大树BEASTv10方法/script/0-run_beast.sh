#!/usr/bin/env bash
set -euo pipefail

DEFAULT_XML='/mnt/f/Onedrive/文档（科研）/脚本/Download/18-phylogenetic/7-大树BEASTv10方法/output/2-beast_skygrid_xml/1-table/mtDNA.beast.v10.thorney_skygrid.xml'
DEFAULT_BEAST_BIN='/home/luolintao/opt/BEASTv10.5.0/bin/beast'

XML_PATH="$DEFAULT_XML"
BEAST_BIN="$DEFAULT_BEAST_BIN"
RUNS=3
THREADS=8
SAVE_EVERY=100000
STATE_STEM='mtDNA.state'
STDOUT_LOG='mtDNA.beast.stdout.log'
LEGACY_SEED=777
TIMESTAMP=''
JOBS_ROOT=''
ORIGINAL_ARGC="$#"

usage() {
    cat <<'EOF'
Usage:
  bash script/0-run_beast.sh [options]

Options:
  --xml PATH          BEAST XML 路径。默认运行 Skygrid XML。
  --runs N           同一个 XML 同时运行的随机种子任务数。默认: 1。
  --threads N        每个 BEAST 任务使用的 CPU 线程数。默认: 16。
  --beast-bin PATH   BEAST 可执行文件路径。默认使用本项目 BEAST v10.5.0。
  --save-every N     BEAST state 保存间隔。默认: 100000。
  --jobs-root DIR    批量任务根目录。默认: XML 所在目录/jobs。
  --timestamp ID     批量任务目录名；复用同一 ID 可对已有 job 独立续跑。
  --help             显示帮助。

Examples:
  bash script/0-run_beast.sh
  bash script/0-run_beast.sh --xml output/2-beast_skygrid_xml/1-table/mtDNA.beast.v10.thorney_skygrid.xml --runs 3 --threads 8
EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

log_msg() {
    local log_file="$1"
    local message="$2"
    echo "[$(date '+%F %T')] $message" | tee -a "$log_file"
}

absolute_path() {
    local path="$1"
    local dir
    local base

    if [[ "$path" == /* ]]; then
        dir="$(dirname "$path")"
        base="$(basename "$path")"
    else
        dir="$(dirname "$PWD/$path")"
        base="$(basename "$path")"
    fi

    (cd "$dir" && printf '%s/%s\n' "$PWD" "$base")
}

require_positive_integer() {
    local name="$1"
    local value="$2"

    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        die "$name 必须是正整数: $value"
    fi
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --xml)
            [[ "$#" -ge 2 ]] || die "--xml 缺少参数"
            XML_PATH="$2"
            shift 2
            ;;
        --runs)
            [[ "$#" -ge 2 ]] || die "--runs 缺少参数"
            RUNS="$2"
            shift 2
            ;;
        --threads)
            [[ "$#" -ge 2 ]] || die "--threads 缺少参数"
            THREADS="$2"
            shift 2
            ;;
        --beast-bin)
            [[ "$#" -ge 2 ]] || die "--beast-bin 缺少参数"
            BEAST_BIN="$2"
            shift 2
            ;;
        --save-every)
            [[ "$#" -ge 2 ]] || die "--save-every 缺少参数"
            SAVE_EVERY="$2"
            shift 2
            ;;
        --jobs-root)
            [[ "$#" -ge 2 ]] || die "--jobs-root 缺少参数"
            JOBS_ROOT="$2"
            shift 2
            ;;
        --timestamp)
            [[ "$#" -ge 2 ]] || die "--timestamp 缺少参数"
            TIMESTAMP="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            die "未知参数: $1"
            ;;
    esac
done

require_positive_integer "--runs" "$RUNS"
require_positive_integer "--threads" "$THREADS"
require_positive_integer "--save-every" "$SAVE_EVERY"

XML_PATH="$(absolute_path "$XML_PATH")"
[[ -f "$XML_PATH" ]] || die "找不到 XML: $XML_PATH"
[[ -x "$BEAST_BIN" || -n "$(command -v "$BEAST_BIN" 2>/dev/null || true)" ]] || die "找不到可执行 BEAST: $BEAST_BIN"

XML_DIR="$(dirname "$XML_PATH")"
XML_NAME="$(basename "$XML_PATH")"
XML_PREFIX="${XML_NAME%.xml}"

if [[ -z "$JOBS_ROOT" ]]; then
    JOBS_ROOT="$XML_DIR/jobs"
else
    JOBS_ROOT="$(absolute_path "$JOBS_ROOT")"
fi

if [[ -z "$TIMESTAMP" ]]; then
    TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
fi

COMMON_ARGS=(
    -working
    -beagle
    -beagle_double
    -threads "$THREADS"
    -save_every "$SAVE_EVERY"
    -save_stem "$STATE_STEM"
)

latest_state_file() {
    local latest

    latest="$(
        for state_file in "${STATE_STEM}"_*; do
            [[ -e "$state_file" ]] || continue
            [[ -s "$state_file" ]] || continue
            state_index="${state_file#${STATE_STEM}_}"
            [[ "$state_index" =~ ^[0-9]+$ ]] || continue
            printf '%s\t%s\n' "$state_index" "$state_file"
        done | sort -n | tail -n 1 | cut -f2-
    )"

    printf '%s\n' "$latest"
}

backup_existing_outputs() {
    local state_file="$1"
    local state_tag
    local backup_dir
    local output_files
    local file

    state_tag="${state_file#${STATE_STEM}_}"
    backup_dir="resume_backup/state_${state_tag}_$(date '+%Y%m%d_%H%M%S')"
    mkdir -p "$backup_dir"

    # BEAST 续跑需要 -overwrite；先备份旧输出，避免历史日志被静默覆盖。
    output_files=(
        "$STDOUT_LOG"
        "${XML_PREFIX}.log"
        "${XML_PREFIX}.trees"
        "${XML_PREFIX}.citations.txt"
        "${XML_PREFIX}.ops"
    )

    for file in "${output_files[@]}"; do
        if [[ -e "$file" ]]; then
            cp -p "$file" "$backup_dir/"
        fi
    done

    echo "$backup_dir"
}

run_beast_job() {
    local job_dir="$1"
    local seed="$2"
    local source_xml="$3"
    local xml_name="$4"
    local stdout_log="$5"
    local latest_state
    local backup_dir

    mkdir -p "$job_dir"
    cp -p "$source_xml" "$job_dir/$xml_name"

    cd "$job_dir"

    log_msg "$stdout_log" "Working directory: $job_dir"
    log_msg "$stdout_log" "XML: $xml_name"
    log_msg "$stdout_log" "Seed: $seed"
    log_msg "$stdout_log" "Threads: $THREADS"

    latest_state="$(latest_state_file)"
    if [[ -n "$latest_state" ]]; then
        log_msg "$stdout_log" "检测到断点文件，续跑: $latest_state"
        backup_dir="$(backup_existing_outputs "$latest_state")"
        log_msg "$stdout_log" "已备份续跑前输出: $backup_dir"

        "$BEAST_BIN" \
            "${COMMON_ARGS[@]}" \
            -overwrite \
            -seed "$seed" \
            -load_state "$latest_state" \
            -force_resume \
            "$xml_name" \
            >> "$stdout_log" 2>&1
    else
        log_msg "$stdout_log" "未检测到断点文件，从头开始运行"

        "$BEAST_BIN" \
            "${COMMON_ARGS[@]}" \
            -seed "$seed" \
            "$xml_name" \
            >> "$stdout_log" 2>&1
    fi

    log_msg "$stdout_log" "BEAST 运行结束"
}

seed_for_index() {
    local index="$1"
    local digits

    digits="$(printf '%s' "$TIMESTAMP" | tr -cd '0-9')"
    if [[ -z "$digits" ]]; then
        digits="$(date '+%s')"
    fi

    printf '%s%02d\n' "$digits" "$index"
}

run_single_in_place() {
    run_beast_job "$XML_DIR" "$LEGACY_SEED" "$XML_PATH" "$XML_NAME" "$STDOUT_LOG"
}

run_batch_jobs() {
    local batch_dir="$1"
    local manifest="$batch_dir/job_manifest.tsv"
    local -a pids=()
    local -a job_dirs=()
    local -a seeds=()
    local -a start_times=()
    local index
    local seed
    local job_dir
    local status_file
    local pid
    local failures=0
    local status
    local exit_code
    local start_time
    local end_time

    mkdir -p "$batch_dir"
    printf 'job_index\tseed\tjob_dir\txml\tthreads\tstatus\tstart_time\tend_time\texit_code\n' > "$manifest"

    for ((index = 1; index <= RUNS; index++)); do
        seed="$(seed_for_index "$index")"
        job_dir="$batch_dir/seed_${seed}"
        status_file="$job_dir/.status"
        mkdir -p "$job_dir"
        start_time="$(date '+%F %T')"
        printf '%s\t%s\t%s\t%s\t%s\tRUNNING\t%s\t\t\n' \
            "$index" "$seed" "$job_dir" "$XML_PATH" "$THREADS" "$start_time" \
            >> "$manifest"

        (
            set +e
            run_beast_job "$job_dir" "$seed" "$XML_PATH" "$XML_NAME" "$STDOUT_LOG"
            exit_code="$?"
            end_time="$(date '+%F %T')"
            if [[ "$exit_code" -eq 0 ]]; then
                status='DONE'
            else
                status='FAILED'
            fi
            printf '%s\t%s\t%s\n' "$status" "$end_time" "$exit_code" > "$status_file"
            exit "$exit_code"
        ) &

        pids+=("$!")
        job_dirs+=("$job_dir")
        seeds+=("$seed")
        start_times+=("$start_time")
    done

    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failures=$((failures + 1))
        fi
    done

    printf 'job_index\tseed\tjob_dir\txml\tthreads\tstatus\tstart_time\tend_time\texit_code\n' > "$manifest"
    for ((index = 1; index <= RUNS; index++)); do
        job_dir="${job_dirs[$((index - 1))]}"
        seed="${seeds[$((index - 1))]}"
        status_file="$job_dir/.status"
        start_time="${start_times[$((index - 1))]}"
        if [[ -f "$status_file" ]]; then
            IFS=$'\t' read -r status end_time exit_code < "$status_file"
        else
            status='FAILED'
            end_time="$(date '+%F %T')"
            exit_code='missing_status'
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$index" "$seed" "$job_dir" "$XML_PATH" "$THREADS" "$status" "${start_time:-}" "$end_time" "$exit_code" \
            >> "$manifest"
    done

    if [[ "$failures" -gt 0 ]]; then
        die "$failures 个 BEAST job 失败，详见: $manifest"
    fi

    echo "[INFO] 批量 BEAST 运行结束: $batch_dir"
    echo "[INFO] 任务清单: $manifest"
}

if [[ "$RUNS" -eq 1 && "$XML_PATH" == "$DEFAULT_XML" && "$ORIGINAL_ARGC" -eq 0 ]]; then
    run_single_in_place
else
    run_batch_jobs "$JOBS_ROOT/$TIMESTAMP"
fi

#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="$PROJECT_ROOT/script/1-merge_result.sh"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

assert_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "[FAIL] 缺少文件: $path" >&2
        exit 1
    fi
}

assert_contains() {
    local path="$1"
    local pattern="$2"
    if ! grep -Fq -- "$pattern" "$path"; then
        echo "[FAIL] $path 未包含: $pattern" >&2
        exit 1
    fi
}

write_log() {
    local path="$1"
    local first_state="$2"
    local second_state="$3"

    cat > "$path" <<EOF
# fake BEAST log
state	posterior
$first_state	-$first_state
$second_state	-$second_state
EOF
}

write_tree() {
    local path="$1"
    local first_state="$2"
    local second_state="$3"

    cat > "$path" <<EOF
#NEXUS
Begin trees;
    Translate
        1 A,
        2 B
    ;
    tree STATE_$first_state = [&R] (1:0.1,2:0.1);
    tree STATE_$second_state = [&R] (1:0.2,2:0.2);
End;
EOF
}

BATCH_DIR="$TMPDIR/jobs/20260604_010203"
mkdir -p "$BATCH_DIR/seed_2026060401020301" "$BATCH_DIR/seed_2026060401020302/resume_backup/state_100_20260604"

write_log "$BATCH_DIR/seed_2026060401020301/mtDNA.beast.v10.thorney_skygrid.log" 0 100
write_tree "$BATCH_DIR/seed_2026060401020301/mtDNA.beast.v10.thorney_skygrid.trees" 0 100
write_log "$BATCH_DIR/seed_2026060401020302/resume_backup/state_100_20260604/mtDNA.beast.v10.thorney_skygrid.log" 0 100
write_tree "$BATCH_DIR/seed_2026060401020302/resume_backup/state_100_20260604/mtDNA.beast.v10.thorney_skygrid.trees" 0 100
write_log "$BATCH_DIR/seed_2026060401020302/mtDNA.beast.v10.thorney_skygrid.log" 100 200
write_tree "$BATCH_DIR/seed_2026060401020302/mtDNA.beast.v10.thorney_skygrid.trees" 100 200

"$SCRIPT" \
    --batch-dir "$BATCH_DIR" \
    --outdir "$BATCH_DIR/merged" \
    > "$TMPDIR/merge.out" 2> "$TMPDIR/merge.err"

MERGE_LOG="$BATCH_DIR/merged/merge.log"
MERGE_TREE="$BATCH_DIR/merged/merge.tree"

assert_file "$MERGE_LOG"
assert_file "$MERGE_TREE"
assert_contains "$TMPDIR/merge.out" "合并完成"

log_states="$(awk 'BEGIN{n=0} /^[0-9]+[[:space:]]/{print $1; n++} END{if (n == 0) exit 1}' "$MERGE_LOG" | paste -sd ',')"
if [[ "$log_states" != "0,100,200" ]]; then
    echo "[FAIL] log state 合并错误: $log_states" >&2
    exit 1
fi

tree_states="$(awk '/^[[:space:]]*tree STATE_[0-9]+/ {state=$2; sub(/^STATE_/, "", state); print state}' "$MERGE_TREE" | paste -sd ',')"
if [[ "$tree_states" != "0,100,200" ]]; then
    echo "[FAIL] tree state 合并错误: $tree_states" >&2
    exit 1
fi

end_count="$(grep -c '^End;$' "$MERGE_TREE")"
if [[ "$end_count" -ne 1 ]]; then
    echo "[FAIL] merge.tree End 数量错误: $end_count" >&2
    exit 1
fi

WORKDIR="$TMPDIR/workdir"
LATEST_BATCH="$WORKDIR/jobs/20260604_020304"
mkdir -p "$LATEST_BATCH/seed_2026060402030401"
write_log "$LATEST_BATCH/seed_2026060402030401/mtDNA.beast.v10.thorney_skygrid.log" 0 100
write_tree "$LATEST_BATCH/seed_2026060402030401/mtDNA.beast.v10.thorney_skygrid.trees" 0 100

"$SCRIPT" \
    --workdir "$WORKDIR" \
    > "$TMPDIR/auto_latest.out" 2> "$TMPDIR/auto_latest.err"

assert_file "$LATEST_BATCH/merge.log"
assert_file "$LATEST_BATCH/merge.tree"
assert_contains "$TMPDIR/auto_latest.out" "自动选择最新批量目录"

echo "[PASS] 1-merge_result.sh 批量结果合并测试通过"

#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="$PROJECT_ROOT/script/0-run_beast.sh"
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

XML_DIR="$TMPDIR/xml"
mkdir -p "$XML_DIR"
cat > "$XML_DIR/test.xml" <<'XML'
<beast>
</beast>
XML

FAKE_BEAST="$TMPDIR/fake_beast.sh"
cat > "$FAKE_BEAST" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> fake_beast.args
for arg in "$@"; do
    case "$arg" in
        *.xml)
            xml="$arg"
            ;;
    esac
done
prefix="${xml%.xml}"
printf '# fake log\nstate\tposterior\n0\t0\n' > "${prefix}.log"
printf '#NEXUS\nBegin trees;\n tree STATE_0 = [&R] (A:0.1);\nEnd;\n' > "${prefix}.trees"
printf 'citation\n' > "${prefix}.citations.txt"
SH
chmod +x "$FAKE_BEAST"

"$SCRIPT" \
    --xml "$XML_DIR/test.xml" \
    --runs 3 \
    --threads 8 \
    --beast-bin "$FAKE_BEAST" \
    --timestamp 20260604_010203 \
    > "$TMPDIR/run.out" 2> "$TMPDIR/run.err"

JOBS_DIR="$XML_DIR/jobs/20260604_010203"
assert_file "$JOBS_DIR/job_manifest.tsv"
assert_contains "$JOBS_DIR/job_manifest.tsv" $'job_index\tseed\tjob_dir\txml\tthreads\tstatus'

job_count="$(find "$JOBS_DIR" -maxdepth 1 -type d -name 'seed_*' | wc -l)"
if [[ "$job_count" -ne 3 ]]; then
    echo "[FAIL] seed job 数量错误: $job_count" >&2
    exit 1
fi

for job_dir in "$JOBS_DIR"/seed_*; do
    assert_file "$job_dir/test.xml"
    assert_file "$job_dir/mtDNA.beast.stdout.log"
    assert_file "$job_dir/fake_beast.args"
    assert_contains "$job_dir/fake_beast.args" "-threads 8"
    assert_contains "$job_dir/fake_beast.args" "-seed"
done

first_job="$(find "$JOBS_DIR" -maxdepth 1 -type d -name 'seed_*' | sort | head -n 1)"
printf 'state\n' > "$first_job/mtDNA.state_100000"
printf 'old stdout\n' > "$first_job/mtDNA.beast.stdout.log"

"$SCRIPT" \
    --xml "$XML_DIR/test.xml" \
    --runs 1 \
    --threads 8 \
    --beast-bin "$FAKE_BEAST" \
    --timestamp 20260604_010203 \
    > "$TMPDIR/resume.out" 2> "$TMPDIR/resume.err"

assert_contains "$first_job/fake_beast.args" "-load_state mtDNA.state_100000"
if ! find "$first_job/resume_backup" -mindepth 1 -maxdepth 1 -type d | grep -q .; then
    echo "[FAIL] 续跑时未创建 resume_backup" >&2
    exit 1
fi

echo "[PASS] 0-run_beast.sh 多 seed 批量运行测试通过"

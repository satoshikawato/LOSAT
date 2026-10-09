#!/usr/bin/env bash
# SFb (E2h) gate: every stage of the full gate, run in sequence, each heavy stage under
# flock "$BUILD_ROOT/oracle.lock". From docs/evidence/losat_web_e2d/gates/s11_gates.sh, changed
# for stage E2h (see NOTES.md). Gate A (sf_gate_a.sh) and V-PERF (sf_perf.sh) are separate runs.
#
#   sf_gates.sh                start a new run (refuses when a run of this GATE exists: use
#                              sf_gates_resume.sh, or FRESH=1 to start another run)
#   sf_gates_resume.sh         continue after a stop: finished stages (marker files) are skipped
#   STAGES="a b" sf_gates.sh   only the named stages (markers still apply; FORCE=1 reruns them)
#
# Build directories   $BUILD_ROOT/$P/{clippy,test,adapter,wasm32-test,native,native-serial,wasi,
#                     wasi-artifacts,reactors}   P=${GATE:-gate-sf}   (one frozen gate build)
# Run output          $BUILD_ROOT/sfb-e2h/gate-<UTC>/   (everything; markers in .done/)
# Tracked run record  docs/evidence/losat_web_e2h/run-<UTC>/   (written by stage "collect": status
#                     table, summaries, hashes, logs of failures; no bulky output)
# Gates record failures in status.tsv and go on; a build failure or a dirty tree stops the run.
set -u
WORK_ROOT=${WORK_ROOT:-/home/kawato/losat-work}
BUILD_ROOT=${BUILD_ROOT:-/home/kawato/.cache/losat-work}
W=${WT:-$WORK_ROOT/.worktrees/web-gui}
NCBI=${NCBI_BIN:-/home/kawato/micromamba/bin}
NCBI_CPP=${NCBI_SRC:-/home/kawato/tools/ncbi-blast}/c++
OLD=/home/kawato/.cache/losat-web-gui-target          # baselines kept there: sf/, s08p/ (api-oracle, verify_refs.py)
API=$OLD/s08p/api-oracle/tblastn_stage_e_local_oracle
BASE=${BASE:-f3048ffde}                                # SF's base commit (sf/base-commit.txt)
P=${GATE:-gate-sf}
SB=$BUILD_ROOT/sfb-e2h
LOCK=$BUILD_ROOT/oracle.lock
SWEEPS=$SB/sweeps                                      # frozen NCBI side of the two SF sweeps
FIXNCBI=$BUILD_ROOT/sf-e2h/fixtures/ncbi               # frozen NCBI side of fasta_input_fixtures.py
E2G_RUN=docs/evidence/losat_web_e2e/run-20261004T163746Z   # E2e's gate: option-sweep classes
E2D_RUN=docs/evidence/losat_web_e2d/run-20261005T133711Z   # S11's gate: v1-requests
STAGE_G_ROOT=/mnt/c/Users/genom/GitHub/LOSAT           # named by the TBLASTN Stage G commands (see NOTES.md)
export RUSTUP_TOOLCHAIN=1.92.0
export PATH=/home/kawato/.local/bin:$PATH
export PYTHONDONTWRITEBYTECODE=1
SELF=$(readlink -f "$0")
POINTER=$SB/$P.run
STAGE_LIST="prep lint clippy tests build pychecks quick-fixtures regression-fixtures sf-fixtures sf-sweeps input-sweeps fast-all v1-wasi vabi-quick option-blastp option-tblastn option-tblastx capture vabi-blastn vabi-blastp vabi-tblastn vabi-tblastx collect"
UNLOCKED="prep lint collect"

step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
# rec NAME LOG cmd...   stdout and stderr to LOG, "exit N" appended, a row in status.tsv; never stops
rec() { local name=$1 log=$2; shift 2; mkdir -p "$(dirname "$log")"; "$@" > "$log" 2>&1; local rc=$?
  echo "exit $rc" >> "$log"; printf '%s\t%s\t%s\t%s\n' "$STAGE" "$name" "$rc" "${log#$OUT/}" >> "$OUT/status.tsv"; return 0; }
# rec2 NAME OUTFILE ERRFILE cmd...   stdout to OUTFILE, stderr to ERRFILE (the sweeps)
rec2() { local name=$1 o=$2 e=$3; shift 3; mkdir -p "$(dirname "$o")" "$(dirname "$e")"; "$@" > "$o" 2> "$e"; local rc=$?
  echo "exit $rc" >> "$e"; printf '%s\t%s\t%s\t%s\n' "$STAGE" "$name" "$rc" "${e#$OUT/}" >> "$OUT/status.tsv"; return 0; }
in_dir() { local d=$1; shift; (cd "$W/$d" && "$@"); }
need_load() { # a pool starts below load 12 (skill losat-oracle-runs)
  for _ in $(seq 60); do l=$(cut -d. -f1 /proc/loadavg); [ "$l" -lt 12 ] && return 0; sleep 30; done
  echo "load average still >= 12 after 30 minutes; going on" >&2; }

load_state() {
  [ -n "${SF_TS:-}" ] || { [ -f "$POINTER" ] && . "$POINTER"; }
  [ -n "${SF_TS:-}" ] || fail "no run for GATE=$P ($POINTER missing)"
  TS=$SF_TS; OUT=$SB/gate-$TS; RUN=$W/docs/evidence/losat_web_e2h/run-$TS; DONE=$OUT/.done
  N=$BUILD_ROOT/$P/native/release/LOSAT
  NS=$BUILD_ROOT/$P/native-serial/release/LOSAT
  WA=$BUILD_ROOT/$P/wasi-artifacts
  R=$BUILD_ROOT/$P/reactors
  mkdir -p "$OUT" "$DONE"
}

# ---------------------------------------------------------------- stages
# Run before every stage (also after a WSL restart): the lexical fixture of Gate A, the capture and the
# fast regressions lives in /tmp/losat-pr5-runtime-cert-*, which a restart wipes.
guard() {
  cd "$W" || fail worktree
  (cd LOSAT/tests && python3 -c 'import ci_fast_regressions as c; c.stage_lexical_fixtures()') || fail lexical-fixtures
  # Only inside the private namespace (SF_NS=1), where $STAGE_G_ROOT is the bound worktree; outside it the path is on
  # /mnt/c (9p), which the gate never reads.
  if [ -n "${SF_NS:-}" ]; then [ -d "$STAGE_G_ROOT/LOSAT" ] || fail "$STAGE_G_ROOT is not bound to the worktree; see NOTES.md"; fi
  [ -d "$FIXNCBI" ] || fail "frozen NCBI fixtures missing: $FIXNCBI"
  if [ -f "$OUT/head.txt" ]; then [ "$(git rev-parse HEAD)" = "$(cat "$OUT/head.txt")" ] || fail "HEAD moved since prep ($(cat "$OUT/head.txt"))"; fi
}
st_prep() {
  cd "$W" || fail worktree
  git rev-parse HEAD > "$OUT/head.txt"; git status --short > "$OUT/worktree-status.txt"
  if git status --porcelain --untracked-files=no | grep -q . && [ "${ALLOW_DIRTY:-0}" != 1 ]; then
    fail "tracked files are modified: the gate builds one committed tree (ALLOW_DIRTY=1 overrides; the record names it)"; fi
  uptime > "$OUT/uptime-start.txt"; free -g >> "$OUT/uptime-start.txt"
}

st_lint() {
  cd "$W" || fail worktree
  rec verify-refs "$OUT/verify-refs.log" python3 "$OLD/s08p/verify_refs.py" $(git diff --name-only --diff-filter=d $BASE..HEAD -- '*.rs')
  rec verify-added "$OUT/verify-refs-session-added.txt" python3 docs/evidence/losat_web_e2h/gates/verify_added.py "$OUT/verify-refs.log" $BASE
  rec protein-tables "$OUT/protein-tables-check.log" python3 docs/evidence/losat_web_e2e/gen_protein_tables.py --ncbi-src "$NCBI_CPP" --check
  rec fmt "$OUT/fmt.log" in_dir LOSAT cargo fmt --check
  rec fmt-adapter "$OUT/fmt-adapter.log" in_dir web/adapter cargo fmt --check
  rec pure-rust-boundary "$OUT/pure-rust-boundary.log" python3 LOSAT/tests/check_pure_rust_runtime_boundary.py --root .
  rec unittest-boundary "$OUT/unittest-boundary.log" python3 -m unittest discover -s LOSAT/tests -p 'test_check_pure_rust_runtime_boundary.py'
  rec unittest-fast "$OUT/unittest-fast.log" python3 -m unittest discover -s LOSAT/tests -p 'test_ci_fast_regressions.py'
}

clippy_cfg() { local dir=$1 d=$2 bad=0 r; shift 2   # one "exit N" per configuration; non-zero when any failed
  for cfg in "$@"; do echo "== cargo clippy --locked $cfg -- -D warnings"
    # shellcheck disable=SC2086
    (cd "$W/$dir" && cargo clippy --locked $cfg --target-dir "$d" -- -D warnings) 2>&1; r=$?; echo "exit $r"; [ $r = 0 ] || bad=1; done
  return $bad; }
st_clippy() {
  cd "$W" || fail worktree
  rec clippy "$OUT/clippy.log" clippy_cfg LOSAT "$BUILD_ROOT/$P/clippy" "--all-targets --all-features" "--all-targets --no-default-features" \
    "--lib --target wasm32-wasip1 --no-default-features" "--lib --target wasm32-wasip1-threads --features wasm-threads"
  rec adapter-clippy "$OUT/adapter-clippy.log" clippy_cfg web/adapter "$BUILD_ROOT/$P/adapter" "--all-targets" \
    "--lib --target wasm32-wasip1" "--lib --target wasm32-wasip1-threads --features threads"
}

st_tests() {
  cd "$W" || fail worktree
  # cargo test --all-features runs every unit test, including the reader's (blastinput::fasta_reader)
  rec cargo-test "$OUT/cargo-test.log" in_dir LOSAT env LOSAT_BLASTX_WORKER_LOG="$BUILD_ROOT/blastx-worker.log" CARGO_PROFILE_TEST_OPT_LEVEL=1 \
    cargo test --locked --all-features --no-fail-fast --target-dir "$BUILD_ROOT/$P/test"
  rec reader-tests "$OUT/reader-tests.log" in_dir LOSAT env CARGO_PROFILE_TEST_OPT_LEVEL=1 \
    cargo test --locked --all-features --target-dir "$BUILD_ROOT/$P/test" blastinput::fasta_reader
  # adapter: all tests, which include tests/scan_properties.rs and tests/scan_ncbi_properties.rs
  rec adapter-test "$OUT/adapter-test.log" in_dir web/adapter cargo test --locked --target-dir "$BUILD_ROOT/$P/adapter" -- --nocapture
  for t in scan_properties scan_ncbi_properties; do
    if [ -f web/adapter/tests/$t.rs ]; then
      rec "adapter-$t" "$OUT/adapter-$t.log" in_dir web/adapter cargo test --locked --target-dir "$BUILD_ROOT/$P/adapter" --test $t
    else
      printf '%s\t%s\t%s\t%s\n' "$STAGE" "adapter-$t" 99 "web/adapter/tests/$t.rs is missing" >> "$OUT/status.tsv"
    fi
  done
  rec wasm32-web-api-tests "$OUT/wasm32-web-api-tests.log" in_dir LOSAT env CARGO_TARGET_WASM32_WASIP1_RUNNER="node $W/web/tools/wasi-test-runner.mjs" \
    cargo test --locked --lib --target wasm32-wasip1 --no-default-features --target-dir "$BUILD_ROOT/$P/wasm32-test" -- web_api::
}

st_build() {
  cd "$W" || fail worktree
  (cd LOSAT && cargo build --release --locked --target-dir "$BUILD_ROOT/$P/native") > "$OUT/build-native.log" 2>&1 || fail native
  (cd LOSAT && cargo build --release --locked --no-default-features --target-dir "$BUILD_ROOT/$P/native-serial") > "$OUT/build-native-serial.log" 2>&1 || fail native-serial
  python3 LOSAT/tests/build_wasi_artifacts.py --target-dir "$BUILD_ROOT/$P/wasi" --output-dir "$WA" --include-serial > "$OUT/build-wasi.log" 2>&1 || fail wasi
  python3 web/adapter/tools/build_reactors.py --target-dir "$BUILD_ROOT/$P/adapter" --output-dir "$R" > "$OUT/build-reactors.log" 2>&1 || fail reactors
  mkdir -p "$OUT/reactors" && cp "$R"/*.json "$OUT/reactors/"
  python3 web/adapter/tools/check_build_identity.py --out "$OUT/reactors/build-identity.json" > /dev/null || fail identity
  sha256sum "$N" "$NS" "$WA"/*.wasm "$R"/*.wasm > "$OUT/artifacts.sha256"
}

st_pychecks() { # engine/python checks that need the binary but not the oracle
  cd "$W" || fail worktree
  for n in 1 2 4; do rec2 check-losat-n$n "$OUT/check-losat-n$n.tsv" "$OUT/check-losat-n$n.err" python3 docs/evidence/losat_web_e2a/check_losat.py --losat "$N" --threads $n; done
  # TMPDIR=/tmp: the -db_gencode stand-in report names its BLAST database path (# Database: $TMPDIR/losat_outfmt0_db/...),
  # frozen under /tmp (SFc, 2026-10-10).
  rec oracle-check "$OUT/oracle-check-gate.log" env TMPDIR=/tmp python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir "$NCBI" --out "$OUT/work/oracle-check"
}
st_quick-fixtures() { st_pychecks; }

st_regression-fixtures() {
  cd "$W" || fail worktree
  need_load
  rec range-fixtures "$OUT/range-fixtures.log" in_dir LOSAT python3 tests/range_regression_fixtures.py check --losat "$N" --jobs 4 --out "$OUT/range-fixtures.tsv"
  rec tblastx-fixtures "$OUT/tblastx-fixtures.log" in_dir LOSAT python3 tests/tblastx_regression_fixtures.py check --losat "$N" --jobs 4 --out "$OUT/tblastx-fixtures.tsv"
  rec blastn-fixtures "$OUT/blastn-fixtures.log" in_dir LOSAT python3 tests/blastn_regression_fixtures.py check --losat "$N" --jobs 4 --out "$OUT/blastn-fixtures.tsv"
}

st_sf-fixtures() { # SF's NCBI-frozen FASTA-input fixtures, every program and tier
  cd "$W" || fail worktree
  need_load
  rec sf-fasta-input-fixtures "$OUT/sf-fixtures.log" python3 LOSAT/tests/fasta_input_fixtures.py check --losat "$N" --jobs 3 \
    --expected "$FIXNCBI" --out "$OUT/sf-fixtures.tsv"
}

st_sf-sweeps() { # the two SF sweeps; inputs and frozen NCBI side are made once (deterministic) under $SWEEPS
  cd "$W" || fail worktree
  need_load
  local E=docs/evidence/losat_web_e2h
  mkdir -p "$SWEEPS" "$OUT/sweeps"
  [ -f "$SWEEPS/fasta_sweep/cases.tsv" ] || python3 $E/fasta_sweep.py generate --dir "$SWEEPS/fasta_sweep" > "$OUT/sweeps/fasta-sweep-generate.log" 2>&1 || fail fasta-sweep-generate
  [ -f "$SWEEPS/ncbi/fasta_sweep/manifest.tsv" ] || (cd "$SWEEPS" && python3 "$W/$E/fasta_sweep.py" freeze-ncbi --dir "$SWEEPS/fasta_sweep" \
      --ncbi-bin "$NCBI" --jobs 3 --out "$SWEEPS/ncbi/fasta_sweep") > "$OUT/sweeps/fasta-sweep-freeze.log" 2>&1 || fail fasta-sweep-freeze
  [ -f "$SWEEPS/ncbi/check_inputs/manifest.tsv" ] || (cd "$SWEEPS" && python3 "$W/$E/check_inputs.py" freeze-ncbi --ncbi-bin "$NCBI" \
      --work "$SWEEPS/work/check_inputs" --out "$SWEEPS/ncbi/check_inputs" --jobs 3) > "$OUT/sweeps/check-inputs-freeze.log" 2>&1 || fail check-inputs-freeze
  # both run from the sweep directory (NCBI's side was frozen with relative paths)
  rec fasta-sweep "$OUT/sweeps/fasta-sweep.log" bash -c 'cd "$1" && python3 "$2/docs/evidence/losat_web_e2h/fasta_sweep.py" check --dir "$1/fasta_sweep" --losat "$3" --jobs 3 --ncbi "$1/ncbi/fasta_sweep" --out "$4"' _ "$SWEEPS" "$W" "$N" "$OUT/sweeps/fasta-sweep.tsv"
  rec check-inputs "$OUT/sweeps/check-inputs.log" bash -c 'cd "$1" && python3 "$2/docs/evidence/losat_web_e2h/check_inputs.py" check --losat "$3" --work "$1/work/check_inputs" --ncbi "$1/ncbi/check_inputs" --out "$4" --jobs 3' _ "$SWEEPS" "$W" "$N" "$OUT/sweeps/check-inputs.tsv"
  rec sweep-summary "$OUT/sweeps/summary.md" python3 $E/fasta_sweep.py summary "$OUT/sweeps/fasta-sweep.tsv" "$OUT/sweeps/check-inputs.tsv"
}

st_input-sweeps() { # E2b-E2d sweeps against NCBI that the SF reader can change
  cd "$W" || fail worktree
  need_load
  rec ctoolkit-compare "$OUT/ctoolkit-compare.tsv" python3 docs/evidence/losat_web_e2b/ctoolkit_compare.py --bin-dir "$NCBI" --losat "$N" --jobs 3
  rm -rf "$OUT/work/punct"; rec punct-defline "$OUT/punct-defline.tsv" python3 docs/evidence/losat_web_e2b/punct_defline.py --bin-dir "$NCBI" --losat "$N" --work "$OUT/work/punct"
  rm -rf "$OUT/work/titles"; rec title-sweep "$OUT/title-sweep.tsv" python3 docs/evidence/losat_web_e2e/title_sweep.py --bin-dir "$NCBI" --losat "$N" --work "$OUT/work/titles" --jobs 3
  rm -rf "$OUT/work/ptitles"; rec protein-title-sweep "$OUT/protein-title-sweep.tsv" python3 docs/evidence/losat_web_e2e/protein_title_sweep.py --bin-dir "$NCBI" --losat "$N" --work "$OUT/work/ptitles" --jobs 3
  rec2 range-sweep "$OUT/range-sweep.tsv" "$OUT/range-sweep.err" python3 docs/evidence/losat_web_e2d/range_sweep.py --bin-dir "$NCBI" --losat "$N" --jobs 3 --out "$OUT/range-sweep.tsv"
  rm -rf "$OUT/work/punct" "$OUT/work/titles" "$OUT/work/ptitles"
}

st_fast-all() { # S02 baseline hashes, Gate A and TLOSAN Stage G frozen hashes, fixtures of every program
  cd "$W" || fail worktree
  need_load
  rm -rf "$OUT/fast"
  rec fast-regressions-all "$OUT/fast-regressions-all.log" python3 LOSAT/tests/ci_fast_regressions.py --losat "$N" --out "$OUT/fast" --all-cases --jobs 4
  cp "$OUT/fast/summary.json" "$OUT/fast-regressions-all-summary.json" 2>/dev/null
}

st_v1-wasi() { # the v1 WASI matrix and v1-requests (ABI v1 is frozen: bytes equal S11's)
  cd "$W" || fail worktree
  rm -rf "$OUT/work/wasm-threading"
  rec wasm-threading "$OUT/wasm-threading.log" python3 LOSAT/tests/check_wasm_threading.py --native "$N" --native-serial "$NS" \
    --serial "$WA/losat-serial-command.wasm" --threaded "$WA/losat-threaded-command.wasm" \
    --reactor "$WA/losat-threaded-reactor.wasm" --serial-reactor "$WA/losat-serial-reactor.wasm" \
    --oracle-dir "$NCBI" --output-dir "$OUT/work/wasm-threading"
  cp "$OUT/work/wasm-threading/metadata.json" "$OUT/wasm-threading-metadata.json" 2>/dev/null
  cp "$OUT/work/wasm-threading/runs.json" "$OUT/wasm-threading-runs.json" 2>/dev/null
  node docs/evidence/losat_web_e1a/v1_requests.js "$WA/losat-threaded-reactor.wasm" "$OUT/work/wasm-threading/fixtures/aa3.fasta" > "$OUT/v1-requests-after.jsonl" 2>/dev/null
  rec v1-requests-compare "$OUT/v1-requests-compare.txt" cmp "$E2D_RUN/v1-requests-after.jsonl" "$OUT/v1-requests-after.jsonl"
}

st_vabi-quick() {
  cd "$W" || fail worktree
  mkdir -p "$OUT/v-abi-quick"
  python3 web/adapter/tools/v_abi_cases.py --suite quick --out "$OUT/v-abi-quick/cases.json" > /dev/null || fail cases-quick
  rec v-abi-quick "$OUT/v-abi-quick/v-abi.log" node web/adapter/tests/v_abi.js --native "$N" --serial "$R/losat-web-serial.wasm" \
    --threads "$R/losat-web-threads.wasm" --cases "$OUT/v-abi-quick/cases.json" --out "$OUT/v-abi-quick"
}

sweep_compare() { # sweep_compare PROGRAM : classes of this run against E2e's (same classification)
  python3 - "$W/$E2G_RUN/sweeps/after-$1.tsv" "$OUT/sweeps/after-$1.tsv" <<'PY'
import sys
def rows(p):
    d, n = {}, 0
    for line in open(p, encoding="utf-8", errors="replace"):
        if line.startswith("#") or line.startswith("set\t"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4: d[tuple(f[:3])] = f[3]
    return d
a, b = rows(sys.argv[1]), rows(sys.argv[2])
chg = [(k, a[k], b[k]) for k in a if k in b and a[k] != b[k]]
print(f"E2e rows {len(a)}, this run {len(b)}; removed {len([k for k in a if k not in b])}, added {len([k for k in b if k not in a])}, changed {len(chg)}")
for k in a:
    if k not in b: print("removed", *k, a[k], sep="\t")
for k in b:
    if k not in a: print("added", *k, b[k], sep="\t")
for k, x, y in chg: print("changed", *k, x, "->", y, sep="\t")
sys.exit(1 if chg or len(a) != len(b) or any(v.startswith("DIFF") or v == "timeout" for v in b.values()) else 0)
PY
}
option_sweep() { # option_sweep PROGRAM JOBS
  local p=$1 j=$2
  cd "$W" || fail worktree
  need_load
  rm -rf "$OUT/work/sweep-$p"; mkdir -p "$OUT/work/sweep-$p" "$OUT/sweeps"
  rec2 "option-sweep-$p" "$OUT/sweeps/after-$p.tsv" "$OUT/sweeps/after-$p.err" python3 docs/evidence/losat_web_e2e/option_sweep.py --program $p \
    --bin-dir "$NCBI" --losat "$N" --ncbi-src "$NCBI_CPP" --jobs $j --timeout 1200 --work "$OUT/work/sweep-$p" --api "$API"
  rec "option-sweep-$p-vs-e2e" "$OUT/sweeps/compare-$p.txt" sweep_compare $p
  rm -rf "$OUT/work/sweep-$p"
}
st_option-blastp() { option_sweep blastp 3; }
st_option-tblastn() { option_sweep tblastn 2; }
st_option-tblastx() {
  option_sweep tblastx 2
  rec gencode-api-check "$OUT/gencode-api-check.tsv" python3 docs/evidence/losat_web_e2e/gencode_api_check.py --api "$API" --bin-dir "$NCBI" --losat "$N"
}

st_capture() { # all 236 cases (jobs 4); compared with the S02 baseline and with S11's build (sf/capture-before)
  cd "$W" || fail worktree
  need_load
  rm -rf "$OUT/capture"
  # capture_outputs.py exits 1 on any frozen-hash mismatch; the one known mismatch blastn/Sakai.MG1655.megablast is allowed
  # (ci_fast_regressions.py's allowlist, S02); the two comparisons below are the gate (SFc, 2026-10-10).
  rec capture "$OUT/capture.log" bash -c 'python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat "$1" --out "$2" --jobs 4 > "$2.run.log" 2>&1; rc=$?
    cat "$2.run.log"; [ $rc = 0 ] && exit 0
    other=$(grep "mismatch:" "$2.run.log" | grep -v "^expected-hash mismatch: blastn Sakai.MG1655.megablast$" | wc -l)
    [ "$other" = 0 ] && grep -q "^236 cases, 1 frozen-hash mismatches" "$2.run.log" && { echo "allowed known mismatch: blastn/Sakai.MG1655.megablast"; exit 0; }
    exit $rc' _ "$N" "$OUT/capture"
  rec capture-compare-s02 "$OUT/capture-compare-s02.txt" python3 docs/evidence/losat_web_e1a/capture_outputs.py compare docs/evidence/losat_web_e1a/baseline/hashes.tsv "$OUT/capture/hashes.tsv"
  rec capture-compare-before "$OUT/capture-compare-before.txt" python3 docs/evidence/losat_web_e1a/capture_outputs.py compare "$OLD/sf/capture-before/hashes.tsv" "$OUT/capture/hashes.tsv"
  wc -l < "$OUT/capture/hashes.tsv" > "$OUT/capture-rows.txt"   # 237 = header + 236 cases
}

vabi_program() { # V-ABI full for one program, alone, --jobs 4; cases.json of the full suite is split by program
  local prog=$1
  cd "$W" || fail worktree
  need_load
  local d=$OUT/v-abi-full/$prog; mkdir -p "$d"
  [ -s "$OUT/v-abi-full/cases.json" ] || python3 web/adapter/tools/v_abi_cases.py --suite full --out "$OUT/v-abi-full/cases.json" > /dev/null || fail cases-full
  python3 - "$OUT/v-abi-full/cases.json" "$d/cases.json" "$prog" > "$d/count.txt" <<'PY' || fail split-cases
import json, sys
keep = [s for s in json.load(open(sys.argv[1])) if s["program"] == sys.argv[3]]
json.dump(keep, open(sys.argv[2], "w"), indent=1)
print(len(keep))
PY
  rec "v-abi-full-$prog" "$d/run.log" python3 web/adapter/tools/run_v_abi_parallel.py --native "$N" --serial "$R/losat-web-serial.wasm" \
    --threads "$R/losat-web-threads.wasm" --cases "$d/cases.json" --out "$d" --jobs 4
}
st_vabi-blastn() { vabi_program blastn; }
st_vabi-blastp() { vabi_program blastp; }
st_vabi-tblastn() { vabi_program tblastn; }
st_vabi-tblastx() { vabi_program tblastx; }

# Copy the small results into the tracked run directory; large files keep a hash and a tail.
keep_file() { # keep_file SRC DEST_NAME
  [ -f "$1" ] || return 0
  mkdir -p "$(dirname "$RUN/$2")"
  if [ "$(stat -c %s "$1")" -le 300000 ]; then cp "$1" "$RUN/$2"
  else
    sha256sum "$1" > "$RUN/$2.sha256"
    case "$2" in
      *.tsv) gzip -9 -n -c "$1" > "$RUN/$2.gz"; [ "$(stat -c %s "$RUN/$2.gz")" -le 300000 ] || { rm -f "$RUN/$2.gz"; tail -n 60 "$1" > "$RUN/$2.tail.txt"; } ;;
      *) tail -n 60 "$1" > "$RUN/$2.tail.txt" ;;
    esac
  fi
}
st_collect() {
  cd "$W" || fail worktree
  mkdir -p "$RUN"
  cp "$OUT/head.txt" "$OUT/worktree-status.txt" "$OUT/artifacts.sha256" "$OUT/uptime-start.txt" "$RUN/" 2>/dev/null
  cp "$OUT/stage-times.tsv" "$RUN/" 2>/dev/null
  mkdir -p "$RUN/reactors" && cp "$OUT"/reactors/*.json "$RUN/reactors/" 2>/dev/null
  # last row per (stage, name); the V-ABI parts must add up to the full case list
  python3 - "$OUT" "$RUN" <<'PY'
import json, sys
from pathlib import Path
out, run = Path(sys.argv[1]), Path(sys.argv[2])
last = {}
for line in (out / "status.tsv").read_text().splitlines():
    f = line.split("\t")
    if len(f) == 4: last[(f[0], f[1])] = f
rows = list(last.values())
total = len(json.load(open(out / "v-abi-full/cases.json"))) if (out / "v-abi-full/cases.json").exists() else None
parts = sum(int((out / f"v-abi-full/{p}/count.txt").read_text() or 0) for p in ("blastn", "blastp", "tblastn", "tblastx") if (out / f"v-abi-full/{p}/count.txt").exists())
if total is not None and parts != total: rows.append(["vabi", "v-abi-full-cases-split", "98", f"{parts} of {total} searches"])
with open(run / "status.tsv", "w") as h:
    h.write("stage\tname\trc\tlog\n")
    for f in rows: h.write("\t".join(f) + "\n")
bad = [f for f in rows if f[2] != "0"]
(run / "status-failed.txt").write_text("".join("\t".join(f) + "\n" for f in bad) or "none\n")
print(f"{len(rows)} checks, {len(bad)} with a non-zero exit")
for f in bad: print("  ", *f)
PY
  # logs: every small log; the failed ones in full up to 1 MB, else their tails
  local f
  for f in "$OUT"/*.log "$OUT"/*.txt "$OUT"/*.tsv "$OUT"/*.err "$OUT"/*.json "$OUT"/*.jsonl "$OUT"/sweeps/* ; do
    [ -f "$f" ] || continue
    case "$f" in */status.tsv|*/head.txt|*/worktree-status.txt|*/uptime-start.txt|*/stage-times.tsv) continue;; esac
    keep_file "$f" "${f#$OUT/}"
  done
  mkdir -p "$RUN/sweeps" "$RUN/capture" "$RUN/v-abi-full" "$RUN/v-abi-quick"
  for f in "$OUT"/sweeps/*; do [ -f "$f" ] && keep_file "$f" "sweeps/$(basename "$f")"; done
  cp "$OUT/capture/hashes.tsv" "$OUT/capture/binary.json" "$RUN/capture/" 2>/dev/null
  for p in blastn blastp tblastn tblastx; do
    mkdir -p "$RUN/v-abi-full/$p"
    for f in run.log count.txt v-abi-results.json summary.json; do keep_file "$OUT/v-abi-full/$p/$f" "v-abi-full/$p/$f"; done
  done
  for f in "$OUT"/v-abi-quick/*.log "$OUT"/v-abi-quick/*.json; do keep_file "$f" "v-abi-quick/$(basename "$f")"; done
  # the full log of every failed check (up to 1 MB; else its tail)
  mkdir -p "$RUN/failed"
  while IFS=$'\t' read -r _ _ _ lg; do
    [ -f "$OUT/$lg" ] || continue
    if [ "$(stat -c %s "$OUT/$lg")" -le 1000000 ]; then cp "$OUT/$lg" "$RUN/failed/$(echo "$lg" | tr / _)"
    else tail -n 200 "$OUT/$lg" > "$RUN/failed/$(echo "$lg" | tr / _).tail.txt"; fi
  done < "$RUN/status-failed.txt"
  step "run record: $RUN"
}

# ---------------------------------------------------------------- driver
if [ "${1:-}" = --stage ]; then
  STAGE=$2; load_state
  grep -v "^$STAGE	" "$OUT/status.tsv" > "$OUT/status.tsv.new" 2>/dev/null; mv "$OUT/status.tsv.new" "$OUT/status.tsv" 2>/dev/null || : > "$OUT/status.tsv"
  # STAGE_G_BIND=1 (default): the stages that run the TLOSAN Stage G cases (their commands name $STAGE_G_ROOT/...) see
  # the worktree at that path through a private mount namespace (unshare -rm, no sudo): a tmpfs on /mnt hides the 9p
  # mounts, the path is created in it and the worktree is bound there, so nothing is read from /mnt/c (SFc, 2026-10-09).
  if [ "${STAGE_G_BIND:-1}" = 1 ] && [ -z "${SF_NS:-}" ]; then
    case " capture fast-all vabi-quick vabi-blastn vabi-blastp vabi-tblastn vabi-tblastx " in *" $STAGE "*)
      exec unshare -rm bash -c 'mount -t tmpfs none /mnt && mkdir -p "$2" && mount --bind "$1" "$2" && SF_NS=1 exec "$3" --stage "$4"' _ "$W" "$STAGE_G_ROOT" "$SELF" "$STAGE";; esac
  fi
  start=$(date -u +%FT%TZ)
  [ "$STAGE" = prep ] || [ "$STAGE" = collect ] || guard
  "st_$STAGE"; rc=$?
  printf '%s\t%s\t%s\t%s\n' "$STAGE" "$start" "$(date -u +%FT%TZ)" "$rc" >> "$OUT/stage-times.tsv"
  [ $rc = 0 ] && touch "$DONE/$STAGE"
  exit $rc
fi

if [ "${RESUME:-0}" != 1 ]; then
  if [ -f "$POINTER" ] && [ "${FRESH:-0}" != 1 ]; then
    echo "a run of GATE=$P exists ($(cat "$POINTER")): use sf_gates_resume.sh, or FRESH=1 for another run" >&2; exit 2
  fi
  mkdir -p "$SB"; echo "SF_TS=$(date -u +%Y%m%dT%H%M%SZ)" > "$POINTER"
fi
. "$POINTER"; export SF_TS GATE
load_state
touch "$OUT/status.tsv"
step "gate $P run $TS: out $OUT, record $RUN"
for s in ${STAGES:-$STAGE_LIST}; do
  if [ -e "$DONE/$s" ] && [ "${FORCE:-0}" != 1 ]; then step "skip $s (done)"; continue; fi
  step "stage $s"
  if case " $UNLOCKED " in *" $s "*) true;; *) false;; esac; then "$SELF" --stage "$s"; else flock "$LOCK" "$SELF" --stage "$s"; fi
  rc=$?
  [ $rc = 0 ] || { step "stage $s stopped with rc $rc; fix it and run sf_gates_resume.sh"; exit 1; }
done
step "done (Gate A: sf_gate_a.sh alone; V-PERF: sf_perf.sh)"

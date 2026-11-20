#!/bin/bash
# NP WW0 gate (local; ws has no Node): #136 / #137 on dev. Steps of web.yml (adapter, app, browser-engine),
# the host tests of wasm-threading.yml, the v1 WASI matrix, V-ABI quick (normal and with the heap above 2 GiB),
# and the section comparison of the dev (base) and WW0 reactors.
# Usage: ww0_gate.sh [steps...]   default: all steps in order. Logs: WW/logs/gate_<head>/, one line per step in summary.txt.
set -u
BR=/home/kawato/.cache/losat-work
W=/home/kawato/losat-work/.worktrees/np-ww0
T=/home/kawato/losat-baselines/ncbi-parity-20261019/WW
BASE=$(cat "$T/BASE_DEV_ww0.txt")
TD=$BR/wt-np-ww0
R=$BR/np-ww0
ORACLE=/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin
export RUSTUP_TOOLCHAIN=1.92.0 TMPDIR=$T/tmp
export PATH=/home/kawato/.local/bin:$PATH
unset CARGO_BUILD_JOBS
for v in $(env | grep -o '^LOSAT_[A-Z0-9_]*'); do unset $v; done
unset NODE_ENV VITE_USER_NODE_ENV
cd "$W" || exit 1
HEAD=$(git rev-parse --short HEAD)
G=$T/logs/gate_$HEAD; mkdir -p "$G" "$R" "$TMPDIR"
S=$G/summary.txt
LOCK="flock -w 3600 $BR/oracle.lock"
rec() { echo "$(date +%H:%M:%S) $1 exit=$2 | $3" | tee -a "$S"; }
last() { tail -n 1 "$1" 2>/dev/null | cut -c1-200; }
tally() { grep -E '^\s+[0-9]+ (passed|failed|skipped|flaky|did not run)' "$1" | tr -s ' ' | tr '\n' ';'; }
steps=${*:-"base_reactors reactors sections identity adapter native host_tests v_abi v_abi_high wasi_matrix app_check engine_runtime e2e_engine e2e_fake"}
echo "$(date +%H:%M:%S) head $HEAD base $BASE steps: $steps (uncommitted paths: $(git status --short | wc -l))" | tee -a "$S"
for step in $steps; do case $step in
  base_reactors)
    # dev's reactors from the same target directory: switch the worktree to the base, build, switch back.
    trap 'git -C "$W" switch -q feature/np-ww0' EXIT
    git switch -q --detach "$BASE" && \
      python3 web/adapter/tools/build_reactors.py --target-dir "$TD/adapter" --output-dir "$R/reactors-base" > "$G/base-reactors.log" 2>&1
    rc=$?; git switch -q feature/np-ww0; trap - EXIT
    rec base-reactors $rc "$(python3 -c "import json;d=json.load(open('$R/reactors-base/losat-web-threads.json'));print('threads memory',d['memory'],d['target_rustflags'])" 2>&1 | cut -c1-200)";;
  reactors)
    python3 web/adapter/tools/build_reactors.py --target-dir "$TD/adapter" --output-dir "$R/reactors" > "$G/reactors.log" 2>&1
    rec reactors $? "$(python3 -c "import json;d=json.load(open('$R/reactors/losat-web-threads.json'));print('threads memory',d['memory'],d['target_rustflags'])" 2>&1 | cut -c1-200)";;
  sections)
    for k in serial threads; do echo "== $k"; python3 docs/evidence/wasm_memory_max_20261012/scripts/wasm_sections.py \
      "$R/reactors-base/losat-web-$k.wasm" "$R/reactors/losat-web-$k.wasm"; done > "$G/sections.log" 2>&1
    rec sections $? "$(grep -ci differ "$G/sections.log") lines with 'differ'";;
  identity)
    python3 web/adapter/tools/check_build_identity.py > "$G/identity.log" 2>&1; rec build-identity $? "$(last "$G/identity.log")";;
  adapter)
    (cd web/adapter && cargo fmt --check) > "$G/adapter-fmt.log" 2>&1; rec adapter-fmt $? ""
    (cd web/adapter && cargo test --locked --target-dir "$TD/adapter") > "$G/adapter-test.log" 2>&1
    rec adapter-cargo-test $? "$(grep '^test result' "$G/adapter-test.log" | awk '{p+=$4; f+=$6} END {print p" passed "f" failed"}')";;
  native)
    (cd LOSAT && cargo build --release --locked --target-dir "$TD") > "$G/native.log" 2>&1; rc=$?
    mkdir -p "$R/native" && cp "$TD/release/losat" "$R/native/losat"
    rec native $rc "sha256 $(sha256sum "$R/native/losat" | cut -c1-16)";;
  host_tests)
    node --test LOSAT/tests/test_wasi_runners.js LOSAT/tests/test_wasi_shared_memory.js > "$G/host-node.log" 2>&1
    rec host-node-tests $? "$(grep -E '^(#|ℹ) (tests|pass|fail) ' "$G/host-node.log" | tr '\n' ' ')"
    python3 -m unittest discover -s LOSAT/tests -p test_wasm_performance.py > "$G/host-py.log" 2>&1
    rec host-py-wasm-performance $? "$(last "$G/host-py.log")";;
  v_abi)
    python3 web/adapter/tools/v_abi_cases.py --suite quick --out "$R/v-abi-cases.json" > "$G/v-abi.log" 2>&1
    $LOCK node web/adapter/tests/v_abi.js --native "$R/native/losat" \
      --serial "$R/reactors/losat-web-serial.wasm" --threads "$R/reactors/losat-web-threads.wasm" \
      --cases "$R/v-abi-cases.json" >> "$G/v-abi.log" 2>&1
    rec v-abi-quick $? "$(last "$G/v-abi.log")";;
  v_abi_high)
    # The threaded rows start from a shared memory of 40,000 pages (2.44 GiB): every allocation lies above 2 GiB.
    LOSAT_WASM_MEMORY_INITIAL_PAGES=40000 $LOCK node web/adapter/tests/v_abi.js --native "$R/native/losat" \
      --serial "$R/reactors/losat-web-serial.wasm" --threads "$R/reactors/losat-web-threads.wasm" \
      --cases "$R/v-abi-cases.json" > "$G/v-abi-high.log" 2>&1
    rec v-abi-quick-heap-2.44GiB $? "$(last "$G/v-abi-high.log")";;
  wasi_matrix)
    python3 LOSAT/tests/build_wasi_artifacts.py --target-dir "$TD/wasi" --output-dir "$R/wasi-artifacts" > "$G/wasi-artifacts.log" 2>&1
    rec wasi-artifacts $? "$(last "$G/wasi-artifacts.log")"
    $LOCK python3 LOSAT/tests/check_wasm_threading.py --native "$R/native/losat" \
      --threaded "$R/wasi-artifacts/losat-threaded-command.wasm" --reactor "$R/wasi-artifacts/losat-threaded-reactor.wasm" \
      --oracle-dir "$ORACLE" --output-dir "$R/wasi-gates" > "$G/wasi-matrix.log" 2>&1
    rec wasi-matrix $? "$(last "$G/wasi-matrix.log")";;
  app_check)
    (cd web/app && npm ci) > "$G/app-npm-ci.log" 2>&1; rec app-npm-ci $? "$(last "$G/app-npm-ci.log")"
    (cd web/app && $LOCK npm run check) > "$G/app-check.log" 2>&1
    rec app-check $? "$(grep -E 'Test Files|Tests  ' "$G/app-check.log" | tr -s ' ' | tr '\n' ';')";;
  engine_runtime)
    (cd web/app && LOSAT_WEB_REACTORS="$R/reactors" npx vitest run tests/unit/engine-runtime.test.ts tests/unit/engine-memory.test.ts) > "$G/app-engine-runtime.log" 2>&1
    rec app-engine-runtime-vitest $? "$(grep -E 'Tests|Test Files' "$G/app-engine-runtime.log" | tr -s ' ' | tr '\n' ' ' | cut -c1-200)";;
  e2e_engine)
    mkdir -p "$T/results/evidence-$HEAD"
    (cd web/app && LOSAT_WEB_REACTORS="$R/reactors" LOSAT_WEB_NATIVE="$R/native/losat" LOSAT_WEB_EVIDENCE="$T/results/evidence-$HEAD" \
      LOSAT_WEB_WEBKIT_EXECUTABLE="$R/webkit-deps/pw_run.sh" LOSAT_WEB_E2E_PORT=4397 $LOCK npm run e2e) > "$G/e2e-engine.log" 2>&1
    rec e2e-engine $? "$(tally "$G/e2e-engine.log")";;
  e2e_fake)
    (cd web/app && LOSAT_WEB_WEBKIT_EXECUTABLE="$R/webkit-deps/pw_run.sh" LOSAT_WEB_E2E_PORT=4397 $LOCK npm run e2e) > "$G/e2e-fake.log" 2>&1
    rec e2e-fake $? "$(tally "$G/e2e-fake.log")";;
esac; done
echo "$(date +%H:%M:%S) done" | tee -a "$S"

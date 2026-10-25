#!/usr/bin/env bash
# Engine-side verification of the worktree's tree with the "new" build (out-new).
# Steps: native build, host unit tests, V-ABI quick (normal and with the heap above
# 2.4 GiB), the v1 WASI matrix. Logs: $T/logs/verify/. Usage: verify_engine.sh [steps...]
set -uo pipefail
WT=/home/kawato/losat-work/.worktrees/wasm-mem
BR=/home/kawato/.cache/losat-work/wasm-mem
T=/home/kawato/losat-baselines/wasm-mem-20261012
NODE=/home/kawato/.local/bin/node
ORACLE=/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin
LOCK=/home/kawato/.cache/losat-work/oracle.lock
export PATH=$HOME/.cargo/bin:/home/kawato/.local/bin:$PATH RUSTUP_TOOLCHAIN=1.92.0 TMPDIR=$T/tmp NODE_NO_WARNINGS=1
for v in $(env | grep -oE '^(LOSAT_[A-Z0-9_]*|RUSTFLAGS|CARGO_ENCODED_RUSTFLAGS|CARGO_BUILD_[A-Z_]*|CARGO_PROFILE_[A-Z_]*|CARGO_TARGET_[A-Z_]*)' || true); do unset "$v"; done
L=$T/logs/verify; mkdir -p "$L"
steps=${*:-"native unit vabi vabi_high matrix"}
N=$BR/target-native/release/LOSAT
cd "$WT"
status() { printf '%s %s rc=%s (%s)\n' "$(date +%T)" "$1" "$2" "$3" | tee -a "$L/summary.txt"; }
for step in $steps; do case $step in
  native)
    (cd LOSAT && cargo build --release --locked --target-dir "$BR/target-native") > "$L/native.log" 2>&1
    status native $? "$(sha256sum "$N" | cut -c1-12)";;
  unit)
    $NODE --test LOSAT/tests/test_wasi_runners.js LOSAT/tests/test_wasi_shared_memory.js > "$L/unit_node.log" 2>&1; rc=$?
    status unit_node $rc "$(grep -E '^ℹ (pass|fail)' "$L/unit_node.log" | tr '\n' ' ')"
    python3 -m unittest discover -s LOSAT/tests -p test_wasm_performance.py > "$L/unit_py1.log" 2>&1; rc=$?
    status unit_py_wasm_performance $rc "$(tail -1 "$L/unit_py1.log")"
    python3 -m unittest discover -s LOSAT/tests -p test_simple_comparison.py > "$L/unit_py2.log" 2>&1; rc=$?
    status unit_py_simple_comparison $rc "$(tail -1 "$L/unit_py2.log")"
    python3 web/adapter/tools/check_build_identity.py --out "$L/build_identity.json" > /dev/null 2>&1
    status build_identity $? "$(python3 -c "import json;d=json.load(open('$L/build_identity.json'));print('identical', d['identical'], d['failures'])")";;
  vabi|vabi_high)
    python3 web/adapter/tools/v_abi_cases.py --suite quick --out "$L/v-abi-cases.json" > "$L/$step.cases.log" 2>&1
    if [ $step = vabi_high ]; then export LOSAT_WASM_MEMORY_INITIAL_PAGES=40000; fi
    flock "$LOCK" $NODE web/adapter/tests/v_abi.js --native "$N" --serial "$BR/out-new/adapter/losat-web-serial.wasm" \
      --threads "$BR/out-new/adapter/losat-web-threads.wasm" --cases "$L/v-abi-cases.json" > "$L/$step.log" 2>&1; rc=$?
    unset LOSAT_WASM_MEMORY_INITIAL_PAGES
    status $step $rc "$(grep -E '^ℹ (tests|pass|fail)|^# (tests|pass|fail)' "$L/$step.log" | tr '\n' ' ')";;
  matrix)
    rm -rf "$BR/matrix"
    flock "$LOCK" python3 LOSAT/tests/check_wasm_threading.py --native "$N" \
      --threaded "$BR/out-new/engine/losat-threaded-command.wasm" --reactor "$BR/out-new/engine/losat-threaded-reactor.wasm" \
      --oracle-dir "$ORACLE" --node $NODE --output-dir "$BR/matrix" > "$L/matrix.log" 2>&1; rc=$?
    status matrix $rc "$(tail -2 "$L/matrix.log" | tr '\n' ' ' | cut -c1-200)";;
esac; done
status done 0 "$steps"

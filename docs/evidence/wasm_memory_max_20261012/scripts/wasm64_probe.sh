#!/usr/bin/env bash
# Can the LOSAT library be compiled for wasm64-unknown-unknown as it is?
# Uses the installed 1.97.0 toolchain with rust-src and RUSTC_BOOTSTRAP=1 for -Zbuild-std
# (no nightly is installed; the pinned 1.92.0 toolchain is not touched). No source change.
# Usage: wasm64_probe.sh [extra RUSTFLAGS]   Logs: $T/logs/wasm64/
set -uo pipefail
T=/home/kawato/losat-baselines/wasm-mem-20261012
WT=/home/kawato/losat-work/.worktrees/wasm-mem
BR=/home/kawato/.cache/losat-work/wasm-mem
export PATH=$HOME/.cargo/bin:$PATH TMPDIR=$T/tmp
L=$T/logs/wasm64; mkdir -p "$L"
TC=1.97.0
rustup component add rust-src --toolchain $TC > "$L/rust-src.log" 2>&1; echo "rust-src rc=$?" | tee "$L/summary.txt"
cd "$WT/LOSAT"
for variant in "serial --no-default-features" ; do
  set -- $variant; name=$1; shift
  RUSTC_BOOTSTRAP=1 RUSTFLAGS="${EXTRA_RUSTFLAGS:-}" cargo +$TC build --release --lib "$@" --target wasm64-unknown-unknown \
    -Z build-std=std,panic_abort --target-dir "$BR/target-wasm64" > "$L/build_$name.log" 2>&1
  rc=$?
  echo "build $name rc=$rc errors=$(grep -c '^error' "$L/build_$name.log") warnings=$(grep -c '^warning' "$L/build_$name.log")" | tee -a "$L/summary.txt"
  grep -E '^error(\[E[0-9]+\])?:' "$L/build_$name.log" | sort | uniq -c | sort -rn | head -15 | tee -a "$L/summary.txt"
done
ls -la "$BR/target-wasm64/wasm64-unknown-unknown/release/"*.wasm 2>/dev/null | tee -a "$L/summary.txt"

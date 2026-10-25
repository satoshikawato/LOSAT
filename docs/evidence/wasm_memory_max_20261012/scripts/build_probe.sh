#!/usr/bin/env bash
# Build memprobe for wasm32-wasip1-threads with a given --max-memory (bytes) or "default".
# Usage: build_probe.sh <label> <max-memory bytes | default>
set -euo pipefail
label=$1; max=$2
export PATH=$HOME/.cargo/bin:$PATH RUSTUP_TOOLCHAIN=1.92.0
BR=/home/kawato/.cache/losat-work/wasm-mem
cd /home/kawato/losat-baselines/wasm-mem-20261012/scripts/memprobe
flags="-C target-feature=+simd128"
[ "$max" = default ] || flags="$flags -C link-arg=--max-memory=$max"
RUSTFLAGS="$flags" cargo build --release --target wasm32-wasip1-threads --target-dir "$BR/target-probe-$label" 2>&1 | tail -3
cp "$BR/target-probe-$label/wasm32-wasip1-threads/release/memprobe.wasm" "$BR/memprobe-$label.wasm"
/home/kawato/.local/bin/node /home/kawato/losat-work/.worktrees/wasm-mem/LOSAT/tests/wasi_artifact.js "$BR/memprobe-$label.wasm" threaded-command | python3 -c "import json,sys; d=json.load(sys.stdin); print('$label', d['kind'], d['memory'])"

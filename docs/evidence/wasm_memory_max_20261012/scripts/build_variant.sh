#!/usr/bin/env bash
# Build the threaded command, the threaded reactor (engine) and the two adapter reactors
# of the worktree's current tree into $BUILD_ROOT/wasm-mem/out-<label>.
# Usage: build_variant.sh <label>
set -euo pipefail
label=$1
WT=/home/kawato/losat-work/.worktrees/wasm-mem
BR=/home/kawato/.cache/losat-work/wasm-mem
export PATH=$HOME/.cargo/bin:/home/kawato/.local/bin:$PATH
export RUSTUP_TOOLCHAIN=1.92.0 TMPDIR=/home/kawato/losat-baselines/wasm-mem-20261012/tmp
for v in $(env | grep -oE '^(RUSTFLAGS|CARGO_ENCODED_RUSTFLAGS|CARGO_BUILD_[A-Z_]*|CARGO_PROFILE_[A-Z_]*|CARGO_TARGET_[A-Z_]*)' || true); do unset "$v"; done
out=$BR/out-$label
rm -rf "$out"; mkdir -p "$out/engine" "$out/adapter"
cd "$WT"
git rev-parse HEAD > "$out/head.txt"; git status --short > "$out/status.txt"
date +"%F %T start engine"
python3 LOSAT/tests/build_wasi_artifacts.py --output-dir "$out/engine" --target-dir "$BR/target-engine" --node /home/kawato/.local/bin/node
date +"%F %T start adapter"
python3 web/adapter/tools/build_reactors.py --target-dir "$BR/target-adapter" --output-dir "$out/adapter" --node /home/kawato/.local/bin/node
date +"%F %T done"

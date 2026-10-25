#!/usr/bin/env bash
# Run memprobe for each declared maximum and block size with the worktree's Node host.
# Usage: run_probe_matrix.sh <host dir> <out prefix> [labels...]
set -uo pipefail
host=$1; prefix=$2; shift 2
labels=${*:-"def p65535 p65536"}
BR=/home/kawato/.cache/losat-work/wasm-mem
T=/home/kawato/losat-baselines/wasm-mem-20261012
for l in $labels; do for mib in 256 64 8 4 1; do
  log=$T/logs/${prefix}_${l}_b${mib}.log
  NODE_NO_WARNINGS=1 timeout 120 /usr/bin/time -f 'wall %e s maxrss %M KB' /home/kawato/.local/bin/node $host/run_losat_wasi_threads.js $BR/memprobe-$l.wasm $mib 4 $T/tmp/scratch-$l-$mib.bin > $log 2>&1
  rc=$?
  printf '%s b%-3s rc=%s | %s | %s | %s | %s | %s | %s\n' "$l" "$mib" "$rc" \
    "$(grep -o 'blocks=[0-9]* total_mib=[0-9]*' $log)" "$(grep -o 'small_chunks=.*' $log)" \
    "$(grep -o 'bulk_ok=[a-z]*' $log) $(grep -o 'tag_mismatches=[0-9]*' $log | tr '\n' ' ')" \
    "$(grep -o 'stacks_above_2g=[0-9]*' $log)" "$(grep -o 'file_ok=[a-z]*' $log) $(grep -o 'end pages=[0-9]*' $log)" "$(grep -o 'wall .*' $log)"
done; done

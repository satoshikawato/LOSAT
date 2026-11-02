#!/usr/bin/env bash
# Line WO modes: does the widened window trap in every mode, and does the fix hold there? The widened
# dev module, 16 threads, the 9 daily cases, the base host against the fix, arms in turn.
# usage: wo_modes.sh <cpus> <memnode> <rounds>
set -u
source $HOME/losat-bench/v4prog37/scripts/vp_cases.sh
CPUS=$1; MEMNODE=$2; ROUNDS=$3
B=$HOME/losat-bench; W=$B/wo; OUT=$W/results/modes; NODE=$B/tools/node-v26.8.2-linux-x64/bin/node
mkdir -p $OUT/err
TSV=$OUT/runs.tsv
[ -f "$TSV" ] || printf 'time\tcase\tround\tarm\twall\tmaxrss_kb\trc\tsha16\n' > $TSV
for ((r = 0; r < ROUNDS; r++)); do
  for c in $CASES; do
    for a in $([ $((r % 2)) = 0 ] && echo base fix || echo fix base); do
      H=$W/hostclean_wg; [ $a = fix ] && H=$W/hostclean_cell
      env -i HOME="$HOME" PATH="$PATH" XDG_CACHE_HOME=$W/wcache_m$a NODE_NO_WARNINGS=1 taskset -c $CPUS numactl --membind=$MEMNODE \
        /usr/bin/time -o $OUT/t -f '%e %M' $NODE $H/run_losat_wasi_threads.js $W/bin/d41amp100k.wasm $(largs $c) -num_threads 16 -out $OUT/o 2> $OUT/e
      RC=$?
      read -r WALL RSS < <(grep -E '^[0-9.]+ [0-9]+$' $OUT/t | tail -n 1)
      printf '%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\n' "$(date +%T)" $c $r $a "${WALL:-?}" "${RSS:-?}" $RC "$(sha256sum $OUT/o | cut -c1-16)" >> $TSV
      [ "$RC" -ne 0 ] && cp $OUT/e $OUT/err/$c.$a.r$r.err
    done
  done
done
rm -f $OUT/o $OUT/e $OUT/t

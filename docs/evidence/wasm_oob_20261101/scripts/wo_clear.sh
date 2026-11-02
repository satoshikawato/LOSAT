#!/usr/bin/env bash
# Line WO: the clearing runs of release DECISIONS 54 on x64 (line DECISIONS 15): the widened module, dc-megablast and
# blastp at 16 threads, N=40 per arm, untimed. after = merged dev's host; before = #693's host, in the same run, because
# the widened rate depends on the CPUs (dc-megablast 7/40 on one set, 0/10 on another).
set -u
source $HOME/losat-bench/v4prog37/scripts/vp_cases.sh
B=$HOME/losat-bench; W=$B/wo; OUT=$W/results/clear; NODE=$B/tools/node-v26.8.2-linux-x64/bin/node
CPUS=16-23,80-87; MEMNODE=0; ROUNDS=40
M=$B/state/NODE0_SHARED_wo-clear.txt
mkdir -p $OUT/err
for ((i = 0; i < 180; i++)); do
  python3 $W/scripts/cpuguard.py "$CPUS" --timed > /dev/null && break
  sleep 60
done
python3 $W/scripts/cpuguard.py "$CPUS" --timed || { echo GUARD_REFUSED; echo CLEAR_DONE; exit 3; }
echo "shared by line WO (local losat-baselines-4c, /home/kawato/losat-baselines/wasm-oob-20261101) since $(date '+%F %T'): node 0 CPUs 16-23, siblings 80-87 ONLY, UNTIMED trap counts (~40 min; wo/scripts/wo_clear.sh); driver PID $$; removed at its end" > "$M"
sleep 3
python3 $W/scripts/cpuguard.py "$CPUS" || { rm -f "$M"; echo GUARD_REFUSED; echo CLEAR_DONE; exit 3; }
TSV=$OUT/runs.tsv
[ -f "$TSV" ] || printf 'time\tcase\tround\tarm\twall\tmaxrss_kb\trc\tsha16\n' > $TSV
for ((r = 0; r < ROUNDS; r++)); do
  for c in bp dc; do
    for a in $([ $((r % 2)) = 0 ] && echo before after || echo after before); do
      H=$W/hostclean_wg; [ $a = after ] && H=$W/hostmerged
      env -i HOME="$HOME" PATH="$PATH" XDG_CACHE_HOME=$W/wcache_c$a NODE_NO_WARNINGS=1 taskset -c $CPUS numactl --membind=$MEMNODE \
        /usr/bin/time -o $OUT/t -f '%e %M' $NODE $H/run_losat_wasi_threads.js $W/bin/d41amp100k.wasm $(largs $c) -num_threads 16 -out $OUT/o 2> $OUT/e
      RC=$?
      read -r WALL RSS < <(grep -E '^[0-9.]+ [0-9]+$' $OUT/t | tail -n 1)
      printf '%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\n' "$(date +%T)" $c $r $a "${WALL:-?}" "${RSS:-?}" $RC "$(sha256sum $OUT/o | cut -c1-16)" >> $TSV
      [ "$RC" -ne 0 ] && cp $OUT/e $OUT/err/$c.$a.r$r.err
    done
  done
done
rm -f "$M" $OUT/o $OUT/e $OUT/t
echo CLEAR_DONE

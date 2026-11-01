#!/usr/bin/env bash
# Line WO ab2: time and max RSS of the fixed guard against its base (line WG's host), in one run: the 9 daily cases at
# THREADS threads, arms in turn with the order rotated per round, a null arm (dev again, another
# copy of the host), 3 s warm-up before each run. Each arm has its own cache directory, because the
# guard cache index of one directory goes stale when arms with different guards take turns. Round 0
# fills the caches and is not used.
# usage: wo_ab2.sh <tag> <cpus for taskset> <cpus text> <siblings text> <memnode> <threads> <rounds> <cases>
set -u
source $HOME/losat-bench/v4prog37/scripts/vp_cases.sh
TAG=$1; CPUS=$2; CT=$3; ST=$4; MEMNODE=$5; THREADS=$6; ROUNDS=$7; shift 7; CASELIST="$*"
B=$HOME/losat-bench; W=$B/wo; OUT=$W/results/$TAG; NODE=$B/tools/node-v26.8.2-linux-x64/bin/node
mkdir -p $OUT/err
M=$B/state/NODE${MEMNODE}_IN_USE_wo-$TAG.txt
echo "in use by line WO (local losat-baselines-4c, /home/kawato/losat-baselines/wasm-oob-20261101) since $(date '+%F %T'): node $MEMNODE CPUs $CT, siblings $ST (kept idle), TIMED A/B of the Wasm host, $THREADS threads (wo/scripts/wo_ab2.sh); driver PID $$; removed at its end" > "$M"
sleep 3
python3 $W/scripts/cpuguard.py "$CPUS" --timed || { rm -f "$M"; echo GUARD_REFUSED; echo AB_DONE; exit 3; }
TSV=$OUT/runs.tsv
[ -f "$TSV" ] || printf 'time\tcase\tround\tarm\twall\tuser\tsys\tmaxrss_kb\trc\tsha16\tload1\n' > $TSV
ARMS=(dev null fix)
host() { case $1 in dev) echo $W/hostclean_wg ;; null) echo $W/hostclean_null2 ;; fix) echo $W/hostclean_fix ;; esac; }
for ((r = 0; r < ROUNDS; r++)); do
  for c in $CASELIST; do
    for ((k = 0; k < 3; k++)); do
      a=${ARMS[$(((k + r) % 3))]}
      timeout 3 taskset -c $CPUS bash -c 'while :; do :; done' & 
      wait $! 2>/dev/null
      L=$(cut -d' ' -f1 /proc/loadavg)
      env -i HOME="$HOME" PATH="$PATH" XDG_CACHE_HOME=$W/wcache_$a NODE_NO_WARNINGS=1 taskset -c $CPUS numactl --membind=$MEMNODE \
        /usr/bin/time -o $OUT/t -f '%e %U %S %M' $NODE $(host $a)/run_losat_wasi_threads.js $W/bin/d41.wasm $(largs $c) -num_threads $THREADS -out $OUT/o 2> $OUT/e
      RC=$?
      read -r WALL USER SYS RSS < <(grep -E '^[0-9.]+ [0-9.]+ [0-9.]+ [0-9]+$' $OUT/t | tail -n 1)
      printf '%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n' "$(date +%T)" $c $r $a "${WALL:-?}" "${USER:-?}" "${SYS:-?}" "${RSS:-?}" $RC "$(sha256sum $OUT/o | cut -c1-16)" "$L" >> $TSV
      [ "$RC" -ne 0 ] && cp $OUT/e $OUT/err/$c.$a.r$r.err
    done
  done
done
rm -f "$M" $OUT/o $OUT/e $OUT/t
echo AB_DONE

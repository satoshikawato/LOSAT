#!/usr/bin/env bash
# Line WO: runs arms in turn, each on the same command. An arm = name, host dir, module, WO_PTR.
# usage: wo_arms.sh <out dir> <rounds> <cpulist> <node for membind> <threads> <armfile>
# armfile lines: <name> <host dir> <module> <WO_PTR>
# Stops early when the file <out dir>/STOP exists. A failed run keeps its stderr.
set -u
OUT=$1; ROUNDS=$2; CPUS=$3; MEMNODE=$4; THREADS=$5; ARMS=$6
B=$HOME/losat-bench
NODE=$B/tools/node-v26.8.2-linux-x64/bin/node
Q=$B/scaling/in/dmel_cdna_first31168.fa
S=$B/scaling/in/dmel_dna_3M.fa
mkdir -p "$OUT/err" "$OUT/tmp"
export XDG_CACHE_HOME=$B/wo/wcache
mapfile -t LINES < <(grep -v '^#' "$ARMS")
TSV=$OUT/runs.tsv
[ -f "$TSV" ] || printf 'time\tround\tarm\twall\tuser\tsys\tmaxrss_kb\trc\tsha16\n' > "$TSV"
for ((r = 0; r < ROUNDS; r++)); do
  for ((i = 0; i < ${#LINES[@]}; i++)); do
    [ -e "$OUT/STOP" ] && exit 0
    read -r NAME HOST MODULE PTR <<< "${LINES[$(((i + r) % ${#LINES[@]}))]}"
    E=$OUT/tmp/$NAME.err
    O=$OUT/tmp/$NAME.out
    T=$OUT/tmp/$NAME.time
    env -i HOME="$HOME" PATH="$PATH" XDG_CACHE_HOME="$XDG_CACHE_HOME" NODE_NO_WARNINGS=1 WO_PTR="$PTR" \
      taskset -c "$CPUS" numactl --membind="$MEMNODE" \
      /usr/bin/time -o "$T" -f '%e %U %S %M' \
      "$NODE" "$HOST/run_losat_wasi_threads.js" "$MODULE" blastn -task dc-megablast -query "$Q" -subject "$S" \
      -outfmt 6 -num_threads "$THREADS" -out "$O" 2> "$E"
    RC=$?
    read -r WALL USER SYS RSS < <(grep -E '^[0-9.]+ [0-9.]+ [0-9.]+ [0-9]+$' "$T" | tail -n 1)
    SHA=$(sha256sum "$O" 2>/dev/null | cut -c1-16)
    printf '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n' "$(date +%T)" "$r" "$NAME" "${WALL:-?}" "${USER:-?}" "${SYS:-?}" "${RSS:-?}" "$RC" "$SHA" >> "$TSV"
    if [ "$RC" -ne 0 ] || [ -s "$E" ]; then cp "$E" "$OUT/err/$NAME.r$r.err"; fi
  done
done
rm -rf "$OUT/tmp"
echo ARMS_DONE

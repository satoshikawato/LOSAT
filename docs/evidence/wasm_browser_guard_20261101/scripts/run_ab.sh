#!/bin/bash
# usage: run_ab.sh ROUNDS THREADS [CASES] [EXTRA node-script options, e.g. "--debug 1"]
# -> timing/results/ab_n<THREADS>_r<ROUNDS>_<stamp>.jsonl (+ .log = guest stderr and progress). Summarize: python3 summarize_ab.py FILE.jsonl
set -u
ROUNDS=${1:?rounds}; THREADS=${2:?threads}; CASES=${3:-}; shift 3 2>/dev/null || shift $#
B=/home/kawato/losat-baselines/wg-wasm-guard-20261101
for v in $(env | grep -o '^LOSAT_[A-Za-z0-9_]*'); do unset $v; done
unset GLIBC_TUNABLES MALLOC_ARENA_MAX
export TMPDIR=$B/tmp NODE_NO_WARNINGS=1
mkdir -p $B/tmp $B/timing/results
OUT=$B/timing/results/ab_n${THREADS}_r${ROUNDS}_$(date +%m%d_%H%M%S).jsonl
ARGS=(--rounds $ROUNDS --threads $THREADS --out $OUT --outdir $B/tmp)
[ -n "$CASES" ] && ARGS+=(--cases $CASES)
# A trap or a thread-start timeout ends the node process (SIGTERM, as in the shipped runner). Record the search it was
# running (OUT.current) as a failed sample and go on with the next round in a new process.
FIRST=0
while :; do
  node $B/timing/wasm_guard_ab.js "${ARGS[@]}" --first-round $FIRST "$@" >> ${OUT%.jsonl}.stdout 2>> ${OUT%.jsonl}.log
  rc=$?
  echo "exit $rc first-round $FIRST -> $OUT"
  [ -f $OUT.current ] || break
  FIRST=$(python3 - "$OUT" "$rc" "$FIRST" <<'PY'
import json, sys
out, rc, first = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
cur = json.load(open(out + ".current"))
cur.update(wall_s=None, exit=rc, error=f"node process ended (exit {rc}) during the search; see the .log", out_sha256=None)
open(out, "a").write(json.dumps(cur) + "\n")
print(max(cur["round"] + 1, first))
PY
)
  rm -f $OUT.current
  [ "$FIRST" -lt "$ROUNDS" ] || break
done
python3 $B/timing/summarize_ab.py $OUT

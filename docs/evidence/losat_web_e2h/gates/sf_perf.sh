#!/usr/bin/env bash
# SFb (E2h) V-PERF: before (S11's last gate artifacts, ~/.cache/losat-web-gui-target/sf/bin/) vs after
# (this run's gate build), alternating on every repetition (E2b perf_cases.py -> measure_perf.py). From
# s11_perf.sh. Needs a quiet machine: start it after sf_gates.sh and sf_gate_a.sh have finished.
#
#   sf_perf.sh                       standard cases + read-heavy cases, --repeat 3 (SUFFIX=1; REPEAT, CASES, READ_CASES override)
#   a case over +5% is measured again with --repeat 5 (Owner 2026-10-01; not 10) into perf-<n>-repeat5.json
#
# Read-heavy cases (perf_cases.py, inputs from gen_perf_inputs.py into $BUILD_ROOT/sfb-e2h/perf-inputs/):
#   blastn-q100k, blastn-genome-1line, blastn-genome-80col, blastp-many.
# PERF_READ_MODES (default native,serial-wasi,threaded-wasi) restricts the modes of the read-heavy run.
# Standard-input cases (native only; STDIN_CASES overrides): blastn-q100k-stdin-file, blastn-q100k-stdin-pipe.
set -u
WORK_ROOT=${WORK_ROOT:-/home/kawato/losat-work}
BUILD_ROOT=${BUILD_ROOT:-/home/kawato/.cache/losat-work}
W=${WT:-$WORK_ROOT/.worktrees/web-gui}
OLD=/home/kawato/.cache/losat-web-gui-target
P=${GATE:-gate-sf}
LOCK=$BUILD_ROOT/oracle.lock
if [ -z "${SF_LOCKED:-}" ]; then SF_LOCKED=1 exec flock "$LOCK" "$0" "$@"; fi
. "$BUILD_ROOT/sfb-e2h/$P.run"
OUT=$BUILD_ROOT/sfb-e2h/gate-$SF_TS
RUN=$W/docs/evidence/losat_web_e2h/run-$SF_TS
REPEAT=${REPEAT:-3}
SUFFIX=${SUFFIX:-1}
CASES=${CASES:-blastp,blastp-fmt0,tblastn,tblastn-fmt0,tblastx,tblastx-multi,tblastx-many,blastn,blastn-large,blastn-large-fmt0,blastn-many}
READ_CASES=${READ_CASES:-blastn-q100k,blastn-genome-1line,blastn-genome-80col,blastp-many}
STDIN_CASES=${STDIN_CASES:-blastn-q100k-stdin-file,blastn-q100k-stdin-pipe}
PERF=docs/evidence/losat_web_e2h/gates/perf_cases.py
export PYTHONDONTWRITEBYTECODE=1
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd "$W" || fail worktree
B=$OLD/sf/bin; BEFORE=$B/LOSAT-native,$B/losat-serial-command.wasm,$B/losat-threaded-command.wasm
A=$BUILD_ROOT/$P; AFTER=$A/native/release/LOSAT,$A/wasi-artifacts/losat-serial-command.wasm,$A/wasi-artifacts/losat-threaded-command.wasm
for f in ${BEFORE//,/ } ${AFTER//,/ }; do [ -f "$f" ] || fail "missing $f"; done
mkdir -p "$OUT/perf" "$RUN"

# deterministic read-heavy inputs (script, not committed data)
PI=$BUILD_ROOT/sfb-e2h/perf-inputs
[ -f "$PI/inputs.sha256" ] || python3 docs/evidence/losat_web_e2h/gates/gen_perf_inputs.py --out "$PI" > "$OUT/perf/gen-inputs.log" 2>&1 || fail gen-inputs
export PERF_INPUTS=$PI
cp "$PI/inputs.sha256" "$RUN/perf-inputs.sha256"

# the V-PERF lock (plan DW-7): the app track pauses while it exists. vperf_lock.sh's APP_WORKTREE names the
# old /mnt/c path; point it at the app track's worktree on the Linux clone.
source "$OLD/s07p-resume/vperf_lock.sh"
APP_WORKTREE=$WORK_ROOT/.worktrees/web-gui-app
trap vperf_lock_release EXIT
vperf_lock_take
# a quiet machine: 1-minute load average under 1 (wait up to 20 minutes inside the script)
for _ in $(seq 40); do l=$(cut -d' ' -f1 /proc/loadavg); awk -v l="$l" 'BEGIN{exit !(l < 1)}' && break; sleep 30; done
echo "load average at start: $(cat /proc/loadavg)" | tee "$OUT/perf/load-start-$SUFFIX.txt"

step "perf standard cases (repeat $REPEAT)"
python3 $PERF run --before "$BEFORE" --after "$AFTER" --repeat "$REPEAT" --cases "$CASES" --out "$OUT/perf/perf-$SUFFIX.json" > "$OUT/perf/perf-$SUFFIX.log" 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check "$OUT/perf/perf-$SUFFIX.json" > "$OUT/perf/perf-check-$SUFFIX.txt" 2>&1; rc1=$?
step "perf read-heavy cases (repeat $REPEAT)"
PERF_MODES=${PERF_READ_MODES:-native,serial-wasi,threaded-wasi} python3 $PERF run --before "$BEFORE" --after "$AFTER" --repeat "$REPEAT" --cases "$READ_CASES" \
  --out "$OUT/perf/perf-read-$SUFFIX.json" > "$OUT/perf/perf-read-$SUFFIX.log" 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check "$OUT/perf/perf-read-$SUFFIX.json" > "$OUT/perf/perf-read-check-$SUFFIX.txt" 2>&1; rc2=$?
step "perf standard-input cases, native (repeat $REPEAT)"
PERF_MODES=native python3 $PERF run --before "$BEFORE" --after "$AFTER" --repeat "$REPEAT" --cases "$STDIN_CASES" \
  --out "$OUT/perf/perf-stdin-$SUFFIX.json" > "$OUT/perf/perf-stdin-$SUFFIX.log" 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check "$OUT/perf/perf-stdin-$SUFFIX.json" > "$OUT/perf/perf-stdin-check-$SUFFIX.txt" 2>&1; rc3=$?

# settle: cases over +5% (or with different output) are measured again with --repeat 5, only those cases
for kind in "" read- stdin-; do
  f=$OUT/perf/perf-${kind}$SUFFIX.json
  [ -f "$f" ] || continue
  over=$(python3 - "$f" <<'PY'
import json, sys
cases = json.load(open(sys.argv[1]))["cases"]
bad = {c["case"] for c in cases if c["after"]["median_s"] / c["before"]["median_s"] > 1.05 or c["after"]["output_sha256"] != c["before"]["output_sha256"]}
print(",".join(sorted(bad)))
PY
)
  if [ -n "$over" ]; then
    step "settling ${kind}cases over +5%: $over (repeat 5)"
    modes=""; [ "$kind" = read- ] && modes=${PERF_READ_MODES:-native,serial-wasi,threaded-wasi}; [ "$kind" = stdin- ] && modes=native
    PERF_MODES=$modes python3 $PERF run --before "$BEFORE" --after "$AFTER" --repeat 5 --cases "$over" --out "$OUT/perf/perf-${kind}$SUFFIX-repeat5.json" \
      > "$OUT/perf/perf-${kind}$SUFFIX-repeat5.log" 2>&1
    python3 docs/evidence/losat_web_e1a/measure_perf.py check "$OUT/perf/perf-${kind}$SUFFIX-repeat5.json" > "$OUT/perf/perf-${kind}check-$SUFFIX-repeat5.txt" 2>&1
  fi
done
vperf_lock_release
cp "$OUT"/perf/perf-*.json "$OUT"/perf/perf-*.txt "$OUT"/perf/perf-*.log "$OUT"/perf/load-start-*.txt "$RUN/" 2>/dev/null
printf 'perf\tperf-standard\t%s\tperf/perf-check-%s.txt\nperf\tperf-read-heavy\t%s\tperf/perf-read-check-%s.txt\nperf\tperf-stdin\t%s\tperf/perf-stdin-check-%s.txt\n' "$rc1" "$SUFFIX" "$rc2" "$SUFFIX" "$rc3" "$SUFFIX" >> "$OUT/status.tsv"
step "done (standard rc $rc1, read-heavy rc $rc2, standard input rc $rc3; the settled cases are in perf-*-repeat5.*)"

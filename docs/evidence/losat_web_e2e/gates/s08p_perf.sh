#!/usr/bin/env bash
# S08+ V-PERF: before (SD's final engine, the session's start) vs after, alternating.
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/s08p
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $S/rundir)
REPEAT=${REPEAT:-3}
SUFFIX=${SUFFIX:-1}
step() { echo "$(date -u +%H:%M:%S) $*"; }
cd $W
N=${N_AFTER:-$A/s08p-gate-native/release/LOSAT}
WA=${WA_AFTER:-$A/s08p-gate-wasi-artifacts}
B=$A/sd-final-native/release/LOSAT; BW=$A/sd-final-wasi-artifacts
step "perf (V-PERF lock: app track paused)"
source $A/s07p-resume/vperf_lock.sh
trap vperf_lock_release EXIT
vperf_lock_take
python3 docs/evidence/losat_web_e2b/perf_cases.py run --before $B,$BW/losat-serial-command.wasm,$BW/losat-threaded-command.wasm --repeat $REPEAT \
  --after $N,$WA/losat-serial-command.wasm,$WA/losat-threaded-command.wasm \
  --cases ${CASES:-blastp,blastp-fmt0,tblastn,tblastn-fmt0,tblastx,tblastx-multi,tblastx-many,blastn-large} --out $RUN/perf-$SUFFIX.json > $RUN/perf-$SUFFIX.log 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check $RUN/perf-$SUFFIX.json > $RUN/perf-check-$SUFFIX.txt 2>&1
vperf_lock_release
step "done"

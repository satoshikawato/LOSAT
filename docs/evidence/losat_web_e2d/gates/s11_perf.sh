#!/usr/bin/env bash
# S11 (E2d) V-PERF: before (E2e's last gate artifacts, the session's start, native 6f070575...)
# vs after (this gate's), alternating. From docs/evidence/losat_web_e2e/gates/s08p_perf.sh.
set -u
A=/home/kawato/.cache/losat-web-gui-target
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $A/s11/rundir)
REPEAT=${REPEAT:-3}
SUFFIX=${SUFFIX:-1}
P=${GATE:-s11-gate}
step() { echo "$(date -u +%H:%M:%S) $*"; }
cd $W
N=$A/$P-native/release/LOSAT
WA=$A/$P-wasi-artifacts
B=$A/s08pb2-gate-native/release/LOSAT; BW=$A/s08pb2-gate-wasi-artifacts
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

#!/usr/bin/env bash
# S08 V-PERF: the cases over the threshold in perf-1, measured again with --repeat 5.
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/s08
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $S/rundir)
step() { echo "$(date -u +%H:%M:%S) $*"; }
cd $W
N=$A/s08-gate-native/release/LOSAT
WA=$A/s08-gate-wasi-artifacts
B=$S/bin/LOSAT-base; BW=$A/e2g-gate-wasi-artifacts
step "perf rerun (V-PERF lock: app track paused)"
source $A/s07p-resume/vperf_lock.sh
trap vperf_lock_release EXIT
vperf_lock_take
python3 docs/evidence/losat_web_e2b/perf_cases.py run --before $B,$BW/losat-serial-command.wasm,$BW/losat-threaded-command.wasm --repeat 5 \
  --after $N,$WA/losat-serial-command.wasm,$WA/losat-threaded-command.wasm \
  --cases ${CASES:-tblastx,tblastx-multi,blastp-fmt0} --out $RUN/perf-${TAG:-2}.json > $RUN/perf-${TAG:-2}.log 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check $RUN/perf-${TAG:-2}.json > $RUN/perf-check-${TAG:-2}.txt 2>&1
vperf_lock_release
step "done"

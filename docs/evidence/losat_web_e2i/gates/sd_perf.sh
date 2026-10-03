#!/usr/bin/env bash
# SD V-PERF: megablast and blastn before (S08b's engine, copied to sd/bin at the session
# start) vs after, alternating; then dc-megablast and blastn-short against NCBI (new paths).
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/sd
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $S/rundir)
NCBI=/home/kawato/micromamba/bin
REPEAT=${REPEAT:-3}
SUFFIX=${SUFFIX:-1}
step() { echo "$(date -u +%H:%M:%S) $*"; }
cd $W
N=$A/sd-gate-native/release/LOSAT
WA=$A/sd-gate-wasi-artifacts
B=$S/bin/native/LOSAT; BW=$S/bin/wasi
step "perf (V-PERF lock: app track paused)"
source $A/s07p-resume/vperf_lock.sh
trap vperf_lock_release EXIT
vperf_lock_take
python3 docs/evidence/losat_web_e2c/perf_cases.py run --before $B,$BW/losat-serial-command.wasm,$BW/losat-threaded-command.wasm --repeat $REPEAT \
  --after $N,$WA/losat-serial-command.wasm,$WA/losat-threaded-command.wasm \
  --cases ${CASES:-blastn,blastn-large,blastn-large-fmt0,blastn-many} --out $RUN/perf-$SUFFIX.json > $RUN/perf-$SUFFIX.log 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check $RUN/perf-$SUFFIX.json > $RUN/perf-check-$SUFFIX.txt 2>&1
if [ "${NCBI_RATIO:-1}" = 1 ]; then
  step "dc-megablast and blastn-short against NCBI (-db), native"
  rm -rf $A/sd-perf-ncbi
  python3 docs/evidence/losat_web_e2i/perf_ncbi.py --losat $N --bin-dir $NCBI --work $A/sd-perf-ncbi --out $RUN/perf-ncbi.json > $RUN/perf-ncbi.log 2>&1
fi
vperf_lock_release
step "done"

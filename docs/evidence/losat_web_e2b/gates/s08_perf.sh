#!/usr/bin/env bash
# S08 V-PERF: before (HEAD of the session start, 2bcb86b1f's engine) vs after, alternating.
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
step "perf (V-PERF lock: app track paused)"
source $A/s07p-resume/vperf_lock.sh
trap vperf_lock_release EXIT
vperf_lock_take
python3 docs/evidence/losat_web_e2b/perf_cases.py run --before $B,$BW/losat-serial-command.wasm,$BW/losat-threaded-command.wasm --repeat 3 \
  --after $N,$WA/losat-serial-command.wasm,$WA/losat-threaded-command.wasm \
  --cases tblastx,tblastx-multi,tblastx-many,tblastn,tblastn-fmt0,blastp,blastp-fmt0,blastn-large-fmt0,blastn-large --out $RUN/perf-1.json > $RUN/perf-1.log 2>&1
python3 docs/evidence/losat_web_e1a/measure_perf.py check $RUN/perf-1.json > $RUN/perf-check-1.txt 2>&1
step "format cost (after only, native, alternating 0/6/7, 3 each)"
python3 - "$N" "$W/LOSAT/tests/fasta" > $RUN/perf-formats.txt 2>&1 <<'PY'
import statistics, subprocess, sys, time
losat, fasta = sys.argv[1:]
cases = {"LC738884 x LC741431": [f"{fasta}/LC738884.fasta", f"{fasta}/LC741431.fasta"],
         "LC738874 x LC738875": [f"{fasta}/LC738874.fasta", f"{fasta}/LC738875.fasta"]}
for name, (q, s) in cases.items():
    times = {f: [] for f in ("6", "0", "7")}
    for f in times:
        subprocess.run([losat, "tblastx", "-query", q, "-subject", s, "-outfmt", f], capture_output=True)
    for rep in range(3):
        for f in (("6", "0", "7") if rep % 2 == 0 else ("7", "0", "6")):
            t = time.perf_counter()
            subprocess.run([losat, "tblastx", "-query", q, "-subject", s, "-outfmt", f], capture_output=True, check=True)
            times[f].append(time.perf_counter() - t)
    print(name, "; ".join(f"outfmt {f}: median {statistics.median(v):.3f}s {[round(x, 3) for x in v]}" for f, v in times.items()))
PY
vperf_lock_release
step "done"

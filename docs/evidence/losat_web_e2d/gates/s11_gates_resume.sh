#!/usr/bin/env bash
# S11 (E2d): the rest of s11_gates.sh after a WSL restart stopped it at the TBLASTX option
# sweep (2026-10-05), with the same run directory and build artifacts (checked against
# artifacts.sha256). The steps left out are recorded in the gate record (README.md, "ゲート"):
# the TBLASTX option sweep (E2e's; S11 does not change the options), the capture (the same
# 236 cases and hashes as the CI fast regressions of this run), and Gate A (its frozen
# hashes are in the fast regressions). V-ABI full had finished BLASTN, BLASTP and TBLASTN
# (v-abi-full/parts); this runs V-ABI for the TBLASTX searches of the E2d fixtures, V-ABI
# quick and the v1 WASI matrix.
set -u
A=/home/kawato/.cache/losat-web-gui-target
B=$A/s11
P=${GATE:-s11-gate}
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $B/rundir)
NCBI=/home/kawato/micromamba/bin
export RUSTUP_TOOLCHAIN=1.92.0
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd $W
(cd LOSAT/tests && python3 -c 'import ci_fast_regressions as c; c.stage_lexical_fixtures()') || fail lexical-fixtures

N=$A/$P-native/release/LOSAT
WA=$A/$P-wasi-artifacts
R=$A/$P-reactors
git rev-parse HEAD > $RUN/head-resume.txt
sha256sum $N $A/$P-native-serial/release/LOSAT $WA/*.wasm $R/*.wasm > $RUN/artifacts-resume.sha256
cmp $RUN/artifacts.sha256 $RUN/artifacts-resume.sha256 > /dev/null || fail artifacts-changed

step "v-abi: tblastx searches of the e2d fixtures"
mkdir -p $RUN/v-abi-tblastx-e2d
python3 -c '
import json, sys
cases = json.load(open(sys.argv[1]))
keep = [s for s in cases if s["program"] == "tblastx" and any(c.startswith("e2d.") for c in s["cases"])]
json.dump(keep, open(sys.argv[2], "w"), indent=1)
print(len(keep))
' $RUN/v-abi-full/cases.json $RUN/v-abi-tblastx-e2d/cases.json > $RUN/v-abi-tblastx-e2d/count.txt || fail cases-e2d
python3 web/adapter/tools/run_v_abi_parallel.py --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-tblastx-e2d/cases.json --out $RUN/v-abi-tblastx-e2d --jobs 6 > $RUN/v-abi-tblastx-e2d/run.log 2>&1
echo "exit $?" >> $RUN/v-abi-tblastx-e2d/run.log

step "v-abi quick"
mkdir -p $RUN/v-abi-quick
python3 web/adapter/tools/v_abi_cases.py --suite quick --out $RUN/v-abi-quick/cases.json > /dev/null
node web/adapter/tests/v_abi.js --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-quick/cases.json --out $RUN/v-abi-quick > $RUN/v-abi-quick/v-abi.log 2>&1
echo "exit $?" >> $RUN/v-abi-quick/v-abi.log

step "v1 wasi matrix"
rm -rf $A/wasm-threading-$P
python3 LOSAT/tests/check_wasm_threading.py --native $N --native-serial $A/$P-native-serial/release/LOSAT \
  --serial $WA/losat-serial-command.wasm --threaded $WA/losat-threaded-command.wasm \
  --reactor $WA/losat-threaded-reactor.wasm --serial-reactor $WA/losat-serial-reactor.wasm \
  --oracle-dir $NCBI --output-dir $A/wasm-threading-$P > $RUN/wasm-threading.log 2>&1
echo "exit $?" >> $RUN/wasm-threading.log
cp $A/wasm-threading-$P/metadata.json $RUN/wasm-threading-metadata.json 2>/dev/null
cp $A/wasm-threading-$P/runs.json $RUN/wasm-threading-runs.json 2>/dev/null
node docs/evidence/losat_web_e1a/v1_requests.js $WA/losat-threaded-reactor.wasm $A/wasm-threading-$P/fixtures/aa3.fasta > $RUN/v1-requests-after.jsonl 2>/dev/null
cmp $W/docs/evidence/losat_web_e2e/run-20261004T163746Z/v1-requests-after.jsonl $RUN/v1-requests-after.jsonl > $RUN/v1-requests-compare.txt 2>&1; echo "cmp exit $?" >> $RUN/v1-requests-compare.txt

step "done (V-PERF separately)"

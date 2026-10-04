#!/usr/bin/env bash
# S08+b (E2e close) gate runs, from s08p_gates.sh with its own build and work directories
# (s08pb-gate-*) and the TBLASTX option sweep at 3 jobs (-threshold +inf needs about 9 GB):
# lint, tests, builds,
# then the gates; logs go into the run directory (docs/evidence/losat_web_e2e/run-*).
# Made from docs/evidence/losat_web_e2i/gates/sd_gates.sh and the S08 gates.
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/s08p
B=$A/s08pb
# Build and work directories $A/$P-*: GATE=s08pb2-gate gives a second gate its own while
# another gate's artifacts are in use (Gate A after the R2A2-1 fix).
P=${GATE:-s08pb-gate}
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $B/rundir)
NCBI=/home/kawato/micromamba/bin
BASE=78c06fe61
export RUSTUP_TOOLCHAIN=1.92.0
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd $W
(cd LOSAT/tests && python3 -c 'import ci_fast_regressions as c; c.stage_lexical_fixtures()') || fail lexical-fixtures

mkdir -p $RUN
step "lint and tests"
git rev-parse HEAD > $RUN/head.txt
git status --short > $RUN/worktree-status.txt
python3 $S/verify_refs.py $(git diff --name-only --diff-filter=d $BASE..HEAD -- '*.rs') > $RUN/verify-refs.log 2>&1; echo "exit $?" >> $RUN/verify-refs.log
python3 $W/docs/evidence/losat_web_e2e/gates/verify_added.py $RUN/verify-refs.log > $RUN/verify-refs-session-added.txt 2>&1; echo "exit $?" >> $RUN/verify-refs-session-added.txt
python3 docs/evidence/losat_web_e2e/gen_protein_tables.py --ncbi-src $S/ncbi/c++ --check > $RUN/protein-tables-check.log 2>&1; echo "exit $?" >> $RUN/protein-tables-check.log
(cd LOSAT && cargo fmt --check) > $RUN/fmt.log 2>&1; echo "exit $?" >> $RUN/fmt.log
(cd web/adapter && cargo fmt --check) >> $RUN/fmt.log 2>&1; echo "adapter exit $?" >> $RUN/fmt.log
: > $RUN/clippy.log
for cfg in "--all-targets --all-features" "--all-targets --no-default-features" "--lib --target wasm32-wasip1 --no-default-features" "--lib --target wasm32-wasip1-threads --features wasm-threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/clippy.log
  (cd LOSAT && cargo clippy --locked $cfg --target-dir $A/$P-clippy -- -D warnings) >> $RUN/clippy.log 2>&1; echo "exit $?" >> $RUN/clippy.log
done
: > $RUN/adapter-clippy.log
for cfg in "--all-targets" "--lib --target wasm32-wasip1" "--lib --target wasm32-wasip1-threads --features threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/adapter-clippy.log
  (cd web/adapter && cargo clippy --locked $cfg --target-dir $A/$P-adapter -- -D warnings) >> $RUN/adapter-clippy.log 2>&1; echo "exit $?" >> $RUN/adapter-clippy.log
done
(cd LOSAT && LOSAT_BLASTX_WORKER_LOG=$A/blastx-worker.log CARGO_PROFILE_TEST_OPT_LEVEL=1 cargo test --locked --all-features --no-fail-fast --target-dir $A/$P-test) > $RUN/cargo-test.log 2>&1; echo "exit $?" >> $RUN/cargo-test.log
(cd web/adapter && cargo test --locked --target-dir $A/$P-adapter -- --nocapture) > $RUN/adapter-test.log 2>&1; echo "exit $?" >> $RUN/adapter-test.log
(cd LOSAT && CARGO_TARGET_WASM32_WASIP1_RUNNER="node $W/web/tools/wasi-test-runner.mjs" cargo test --locked --lib --target wasm32-wasip1 --no-default-features --target-dir $A/$P-wasm32-test -- web_api::tests) > $RUN/wasm32-web-api-tests.log 2>&1; echo "exit $?" >> $RUN/wasm32-web-api-tests.log
{ python3 LOSAT/tests/check_pure_rust_runtime_boundary.py --root .; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_check_pure_rust_runtime_boundary.py'; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_ci_fast_regressions.py'; echo "exit $?"; } > $RUN/ci-python-checks.log 2>&1

step "build native"
(cd LOSAT && cargo build --release --locked --target-dir $A/$P-native) > $RUN/build-native.log 2>&1 || fail native
(cd LOSAT && cargo build --release --locked --no-default-features --target-dir $A/$P-native-serial) > $RUN/build-native-serial.log 2>&1 || fail native-serial
step "build wasi"
python3 LOSAT/tests/build_wasi_artifacts.py --target-dir $A/$P-wasi --output-dir $A/$P-wasi-artifacts --include-serial > $RUN/build-wasi.log 2>&1 || fail wasi
step "build reactors"
python3 web/adapter/tools/build_reactors.py --target-dir $A/$P-adapter --output-dir $A/$P-reactors > $RUN/build-reactors.log 2>&1 || fail reactors
mkdir -p $RUN/reactors && cp $A/$P-reactors/*.json $RUN/reactors/
python3 web/adapter/tools/check_build_identity.py --out $RUN/reactors/build-identity.json > /dev/null || fail identity

N=$A/$P-native/release/LOSAT
WA=$A/$P-wasi-artifacts
R=$A/$P-reactors
sha256sum $N $A/$P-native-serial/release/LOSAT $WA/*.wasm $R/*.wasm > $RUN/artifacts.sha256

step "v-abi full (background)"
mkdir -p $RUN/v-abi-full $RUN/v-abi-quick
python3 web/adapter/tools/v_abi_cases.py --suite full --out $RUN/v-abi-full/cases.json > /dev/null || fail cases-full
python3 web/adapter/tools/run_v_abi_parallel.py --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-full/cases.json --out $RUN/v-abi-full --jobs 6 > $RUN/v-abi-full/run.log 2>&1 &
VABI=$!

step "fixtures"
for n in 1 2 4; do python3 docs/evidence/losat_web_e2a/check_losat.py --losat $N --threads $n > $RUN/check-losat-n$n.tsv 2>&1; echo "exit $?" >> $RUN/check-losat-n$n.tsv; done
python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir $NCBI --out $A/$P-oracle-check > $RUN/oracle-check-gate.log 2>&1; echo "exit $?" >> $RUN/oracle-check-gate.log
(cd LOSAT && python3 tests/tblastx_regression_fixtures.py check --losat $N --jobs 6 --out $RUN/tblastx-fixtures.tsv) > $RUN/tblastx-fixtures.log 2>&1; echo "exit $?" >> $RUN/tblastx-fixtures.log
(cd LOSAT && python3 tests/blastn_regression_fixtures.py check --losat $N --jobs 6 --out $RUN/blastn-fixtures.tsv) > $RUN/blastn-fixtures.log 2>&1; echo "exit $?" >> $RUN/blastn-fixtures.log
python3 docs/evidence/losat_web_e2b/ctoolkit_compare.py --bin-dir $NCBI --losat $N --jobs 6 > $RUN/ctoolkit-compare.tsv 2>&1; echo "exit $?" >> $RUN/ctoolkit-compare.tsv
rm -rf $A/$P-punct
python3 docs/evidence/losat_web_e2b/punct_defline.py --bin-dir $NCBI --losat $N --work $A/$P-punct > $RUN/punct-defline.tsv 2>&1; echo "exit $?" >> $RUN/punct-defline.tsv
rm -rf $A/$P-titles
python3 docs/evidence/losat_web_e2e/title_sweep.py --bin-dir $NCBI --losat $N --work $A/$P-titles --jobs 6 > $RUN/title-sweep.tsv 2>&1; echo "exit $?" >> $RUN/title-sweep.tsv
rm -rf $A/$P-inputs
python3 docs/evidence/losat_web_e2g/check_inputs.py --bin-dir $NCBI --losat $N --work $A/$P-inputs > $RUN/blastn-check-inputs.tsv 2>&1; echo "exit $?" >> $RUN/blastn-check-inputs.tsv
rm -rf $A/$P-fast
python3 LOSAT/tests/ci_fast_regressions.py --losat $N --out $A/$P-fast --all-cases --jobs 6 > $RUN/fast-regressions-all.log 2>&1; echo "exit $?" >> $RUN/fast-regressions-all.log
cp $A/$P-fast/summary.json $RUN/fast-regressions-all-summary.json 2>/dev/null

step "option sweeps (after)"
mkdir -p $RUN/sweeps
for p in blastp tblastn tblastx; do
  rm -rf $A/$P-sweep-$p; mkdir -p $A/$P-sweep-$p
  jobs=6; [ $p = tblastx ] && jobs=3
  python3 docs/evidence/losat_web_e2e/option_sweep.py --program $p --bin-dir $NCBI --losat $N --ncbi-src $S/ncbi/c++ --jobs $jobs --timeout 1200 \
    --work $A/$P-sweep-$p --api $S/api-oracle/tblastn_stage_e_local_oracle > $RUN/sweeps/after-$p.tsv 2> $RUN/sweeps/after-$p.err; echo "exit $?" >> $RUN/sweeps/after-$p.err
done
python3 docs/evidence/losat_web_e2e/gencode_api_check.py --api $S/api-oracle/tblastn_stage_e_local_oracle --bin-dir $NCBI --losat $N > $RUN/gencode-api-check.tsv 2>&1; echo "exit $?" >> $RUN/gencode-api-check.tsv
rm -rf $A/$P-ptitles
python3 docs/evidence/losat_web_e2e/protein_title_sweep.py --bin-dir $NCBI --losat $N --work $A/$P-ptitles --jobs 6 > $RUN/protein-title-sweep.tsv 2>&1; echo "exit $?" >> $RUN/protein-title-sweep.tsv

step "capture"
rm -rf $A/capture-$P
python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat $N --out $A/capture-$P --jobs 6 > $RUN/capture.log 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare docs/evidence/losat_web_e1a/baseline/hashes.tsv $A/capture-$P/hashes.tsv > $RUN/capture-compare.txt 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare $A/capture-sd-final/hashes.tsv $A/capture-$P/hashes.tsv > $RUN/capture-compare-before.txt 2>&1
mkdir -p $RUN/capture && cp $A/capture-$P/hashes.tsv $A/capture-$P/binary.json $RUN/capture/ 2>/dev/null

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
cmp $W/docs/evidence/losat_web_e1d/run-20260928T221453Z/v1-requests-after.jsonl $RUN/v1-requests-after.jsonl > $RUN/v1-requests-compare.txt 2>&1; echo "cmp exit $?" >> $RUN/v1-requests-compare.txt

step "v-abi quick"
python3 web/adapter/tools/v_abi_cases.py --suite quick --out $RUN/v-abi-quick/cases.json > /dev/null
node web/adapter/tests/v_abi.js --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-quick/cases.json --out $RUN/v-abi-quick > $RUN/v-abi-quick/v-abi.log 2>&1
echo "exit $?" >> $RUN/v-abi-quick/v-abi.log

step "wait for v-abi full"
wait $VABI
echo "exit $?" >> $RUN/v-abi-full/run.log

step "done (V-PERF separately)"

#!/usr/bin/env bash
# S08 (E2b, TBLASTX outfmt 0/7) gate runs: lint, tests, builds, then the gates; logs go
# into the S08 run directory (docs/evidence/losat_web_e2b/run-*).
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/s08
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $S/rundir)
NCBI=/home/kawato/micromamba/bin
export RUSTUP_TOOLCHAIN=1.92.0
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd $W
# The Gate A lexical root that frozen hashes name (/tmp is emptied when WSL restarts):
# staged as the CI fast job stages it (ci_fast_regressions.py stage_lexical_fixtures).
(cd LOSAT/tests && python3 -c 'import ci_fast_regressions as c; c.stage_lexical_fixtures()') || fail lexical-fixtures

mkdir -p $RUN
step "lint and tests"
git rev-parse HEAD > $RUN/head.txt
git status --short > $RUN/worktree-status.txt
python3 $A/s07p-resume/verify_refs.py $(git diff --name-only 2bcb86b1f..HEAD -- '*.rs') > $RUN/verify-refs.log 2>&1; echo "exit $?" >> $RUN/verify-refs.log
python3 $S/verify_added.py $RUN/verify-refs.log > $RUN/verify-refs-session-added.txt 2>&1; echo "exit $?" >> $RUN/verify-refs-session-added.txt
(cd LOSAT && cargo fmt --check) > $RUN/fmt.log 2>&1; echo "exit $?" >> $RUN/fmt.log
(cd web/adapter && cargo fmt --check) >> $RUN/fmt.log 2>&1; echo "adapter exit $?" >> $RUN/fmt.log
: > $RUN/clippy.log
for cfg in "--all-targets --all-features" "--all-targets --no-default-features" "--lib --target wasm32-wasip1 --no-default-features" "--lib --target wasm32-wasip1-threads --features wasm-threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/clippy.log
  (cd LOSAT && cargo clippy --locked $cfg --target-dir $A/s08-gate-clippy -- -D warnings) >> $RUN/clippy.log 2>&1; echo "exit $?" >> $RUN/clippy.log
done
: > $RUN/adapter-clippy.log
for cfg in "--all-targets" "--lib --target wasm32-wasip1" "--lib --target wasm32-wasip1-threads --features threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/adapter-clippy.log
  (cd web/adapter && cargo clippy --locked $cfg --target-dir $A/s08-gate-adapter -- -D warnings) >> $RUN/adapter-clippy.log 2>&1; echo "exit $?" >> $RUN/adapter-clippy.log
done
(cd LOSAT && LOSAT_BLASTX_WORKER_LOG=$A/blastx-worker.log cargo test --locked --all-features --no-fail-fast --target-dir $A/s08-gate-test) > $RUN/cargo-test.log 2>&1; echo "exit $?" >> $RUN/cargo-test.log
(cd web/adapter && cargo test --locked --target-dir $A/s08-gate-adapter -- --nocapture) > $RUN/adapter-test.log 2>&1; echo "exit $?" >> $RUN/adapter-test.log
(cd LOSAT && CARGO_TARGET_WASM32_WASIP1_RUNNER="node $W/web/tools/wasi-test-runner.mjs" cargo test --locked --lib --target wasm32-wasip1 --no-default-features --target-dir $A/s08-gate-wasm32-test -- web_api::tests) > $RUN/wasm32-web-api-tests.log 2>&1; echo "exit $?" >> $RUN/wasm32-web-api-tests.log
{ python3 LOSAT/tests/check_pure_rust_runtime_boundary.py --root .; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_check_pure_rust_runtime_boundary.py'; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_ci_fast_regressions.py'; echo "exit $?"; } > $RUN/ci-python-checks.log 2>&1

step "build native"
(cd LOSAT && cargo build --release --locked --target-dir $A/s08-gate-native) > $RUN/build-native.log 2>&1 || fail native
(cd LOSAT && cargo build --release --locked --no-default-features --target-dir $A/s08-gate-native-serial) > $RUN/build-native-serial.log 2>&1 || fail native-serial
step "build wasi"
python3 LOSAT/tests/build_wasi_artifacts.py --target-dir $A/s08-gate-wasi --output-dir $A/s08-gate-wasi-artifacts --include-serial > $RUN/build-wasi.log 2>&1 || fail wasi
step "build reactors"
python3 web/adapter/tools/build_reactors.py --target-dir $A/s08-gate-adapter --output-dir $A/s08-gate-reactors > $RUN/build-reactors.log 2>&1 || fail reactors
mkdir -p $RUN/reactors && cp $A/s08-gate-reactors/*.json $RUN/reactors/
python3 web/adapter/tools/check_build_identity.py --out $RUN/reactors/build-identity.json > /dev/null || fail identity

N=$A/s08-gate-native/release/LOSAT
WA=$A/s08-gate-wasi-artifacts
R=$A/s08-gate-reactors
sha256sum $N $A/s08-gate-native-serial/release/LOSAT $WA/*.wasm $R/*.wasm > $RUN/artifacts.sha256

step "v-abi full (background)"
mkdir -p $RUN/v-abi-full $RUN/v-abi-quick
python3 web/adapter/tools/v_abi_cases.py --suite full --out $RUN/v-abi-full/cases.json > /dev/null || fail cases-full
python3 web/adapter/tools/run_v_abi_parallel.py --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-full/cases.json --out $RUN/v-abi-full --jobs 8 > $RUN/v-abi-full/run.log 2>&1 &
VABI=$!

step "fixtures"
for n in 1 2 4; do python3 docs/evidence/losat_web_e2a/check_losat.py --losat $N --threads $n > $RUN/check-losat-n$n.tsv 2>&1; echo "exit $?" >> $RUN/check-losat-n$n.tsv; done
python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir $NCBI --out $A/s08-gate-oracle-check > $RUN/oracle-check-gate.log 2>&1; echo "exit $?" >> $RUN/oracle-check-gate.log
(cd docs/evidence/losat_web_e2a && python3 precheck_hits.py --bin-dir $NCBI --losat $N > $RUN/precheck.tsv 2>&1; echo "exit $?" >> $RUN/precheck.tsv)
(cd LOSAT && python3 tests/tblastx_regression_fixtures.py check --losat $N --jobs 8 --out $RUN/tblastx-fixtures.tsv) > $RUN/tblastx-fixtures.log 2>&1; echo "exit $?" >> $RUN/tblastx-fixtures.log
(cd LOSAT && python3 tests/tblastx_regression_fixtures.py check --losat $S/bin/LOSAT-base --jobs 8 --out $RUN/tblastx-fixtures-before.tsv) > $RUN/tblastx-fixtures-before.log 2>&1; echo "exit $?" >> $RUN/tblastx-fixtures-before.log
(cd LOSAT && python3 tests/tblastx_regression_fixtures.py check --losat $S/bin/LOSAT-0d533ba76 --jobs 8 --out $RUN/tblastx-fixtures-0d533ba76.tsv) > $RUN/tblastx-fixtures-0d533ba76.log 2>&1; echo "exit $?" >> $RUN/tblastx-fixtures-0d533ba76.log
python3 docs/evidence/losat_web_e2b/env_discrimination.py --losat $N --out $RUN/tblastx-env-discrimination.tsv > $RUN/tblastx-env-discrimination.log 2>&1; echo "exit $?" >> $RUN/tblastx-env-discrimination.log
(cd LOSAT && python3 tests/blastn_regression_fixtures.py check --losat $N --jobs 8 --out $RUN/blastn-fixtures.tsv) > $RUN/blastn-fixtures.log 2>&1; echo "exit $?" >> $RUN/blastn-fixtures.log
python3 docs/evidence/losat_web_e2b/ctoolkit_compare.py --bin-dir $NCBI --losat $N --jobs 8 > $RUN/ctoolkit-compare.tsv 2>&1; echo "exit $?" >> $RUN/ctoolkit-compare.tsv
rm -rf $A/s08-gate-punct
python3 docs/evidence/losat_web_e2b/punct_defline.py --bin-dir $NCBI --losat $N --work $A/s08-gate-punct > $RUN/punct-defline.tsv 2>&1; echo "exit $?" >> $RUN/punct-defline.tsv
rm -rf $A/s08-gate-html
python3 docs/evidence/losat_web_e2b/html_titles.py --bin-dir $NCBI --losat $N --work $A/s08-gate-html --jobs 8 > $RUN/html-titles.tsv 2>&1; echo "exit $?" >> $RUN/html-titles.tsv
rm -rf $A/s08-gate-inputs
python3 docs/evidence/losat_web_e2g/check_inputs.py --bin-dir $NCBI --losat $N --work $A/s08-gate-inputs > $RUN/blastn-check-inputs.tsv 2>&1; echo "exit $?" >> $RUN/blastn-check-inputs.tsv
for f in 0 6 7; do python3 docs/evidence/losat_web_e2c/scoring_sweep.py --bin-dir $NCBI --losat $N --jobs 8 --outfmt $f > $RUN/blastn-scoring-sweep-fmt$f.tsv 2>&1; echo "exit $?" >> $RUN/blastn-scoring-sweep-fmt$f.tsv; done
rm -rf $A/s08-gate-titles
python3 docs/evidence/losat_web_e2g/title_sweep.py --bin-dir $NCBI --losat $N --work $A/s08-gate-titles --jobs 8 > $RUN/blastn-title-sweep.tsv 2>&1; echo "exit $?" >> $RUN/blastn-title-sweep.tsv
{
  echo "# Closed pipe (approved exception 5 for outfmt 0; outfmt 6/7 write failures are non-zero)."
  for f in 0 6 7; do
    (cd LOSAT && $N tblastx -query tests/fasta/LC738874.fasta -subject tests/fasta/LC738875.fasta -outfmt $f 2>$A/s08-pipe-err.txt | head -c 100 > /dev/null; echo "outfmt $f: exit ${PIPESTATUS[0]}; stderr: $(tr '\n' '|' < $A/s08-pipe-err.txt)")
  done
  echo "# Standard output closed at the start (>&-): NCBI exit 6 in outfmt 0, abort (134) in 6/7."
  for prog in tblastx blastn; do for f in 0 6 7; do
    (cd LOSAT && sh -c "exec $N $prog -query tests/fasta/LC738874.fasta -subject tests/fasta/LC738875.fasta -outfmt $f >&-" 2>$A/s08-pipe-err.txt; echo "$prog outfmt $f: exit $?; stderr: $(tr '\n' '|' < $A/s08-pipe-err.txt)")
    (cd LOSAT && sh -c "exec $NCBI/$prog -query tests/fasta/LC738874.fasta -subject tests/fasta/LC738875.fasta -outfmt $f >&-" 2>$A/s08-pipe-err.txt; echo "NCBI $prog outfmt $f: exit $?; stderr: $(head -c 60 $A/s08-pipe-err.txt | tr '\n' '|')")
  done; done
} > $RUN/closed-pipe.txt 2>&1
rm -rf $A/s08-gate-fast
python3 LOSAT/tests/ci_fast_regressions.py --losat $N --out $A/s08-gate-fast --all-cases --jobs 6 > $RUN/fast-regressions-all.log 2>&1; echo "exit $?" >> $RUN/fast-regressions-all.log
cp $A/s08-gate-fast/summary.json $RUN/fast-regressions-all-summary.json 2>/dev/null

step "capture"
rm -rf $A/capture-s08
python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat $N --out $A/capture-s08 --jobs 6 > $RUN/capture.log 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare docs/evidence/losat_web_e1a/baseline/hashes.tsv $A/capture-s08/hashes.tsv > $RUN/capture-compare.txt 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare $S/capture-base/hashes.tsv $A/capture-s08/hashes.tsv > $RUN/capture-compare-before.txt 2>&1
mkdir -p $RUN/capture && cp $A/capture-s08/hashes.tsv $A/capture-s08/binary.json $RUN/capture/ 2>/dev/null

step "v1 wasi matrix"
rm -rf $A/wasm-threading-s08
python3 LOSAT/tests/check_wasm_threading.py --native $N --native-serial $A/s08-gate-native-serial/release/LOSAT \
  --serial $WA/losat-serial-command.wasm --threaded $WA/losat-threaded-command.wasm \
  --reactor $WA/losat-threaded-reactor.wasm --serial-reactor $WA/losat-serial-reactor.wasm \
  --oracle-dir $NCBI --output-dir $A/wasm-threading-s08 > $RUN/wasm-threading.log 2>&1
echo "exit $?" >> $RUN/wasm-threading.log
cp $A/wasm-threading-s08/metadata.json $RUN/wasm-threading-metadata.json 2>/dev/null
cp $A/wasm-threading-s08/runs.json $RUN/wasm-threading-runs.json 2>/dev/null
{
  echo "# v1 reactor records (check_wasi_reactor.js, run by check_wasm_threading.py) of the S05 and S08 reactors."
  echo '$ diff -r -x runs.json s05/reactor-records s08/reactor-records'
  diff -r -x runs.json $A/wasm-threading-s05/reactor-records $A/wasm-threading-s08/reactor-records; echo "exit $?"
  echo "records compared: $(ls $A/wasm-threading-s08/reactor-records | grep -c '\.out$')"
  echo '$ diff -r -x runs.json s05/serial-reactor-records s08/serial-reactor-records'
  diff -r -x runs.json $A/wasm-threading-s05/serial-reactor-records $A/wasm-threading-s08/serial-reactor-records; echo "exit $?"
  echo "records compared: $(ls $A/wasm-threading-s08/serial-reactor-records | grep -c '\.out$')"
} > $RUN/v1-reactor-records-compare.txt 2>&1
node docs/evidence/losat_web_e1a/v1_requests.js $WA/losat-threaded-reactor.wasm $A/wasm-threading-s08/fixtures/aa3.fasta > $RUN/v1-requests-after.jsonl 2>/dev/null
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

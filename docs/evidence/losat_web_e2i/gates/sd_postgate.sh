#!/usr/bin/env bash
# SD post-gate check of the final commit. The gate (sd_gates.sh) ran on b99f42493; the later
# commits change comments, an unused constant and its unit test (engine, no behaviour), and
# add fixtures and evidence. This rebuilds the artifacts from the final commit and repeats
# the checks those changes can reach: lint and tests, the fixtures (outfmt 0 manifest, the
# BLASTN and TBLASTX regression fixtures with the new dc.div_* cases), CI's fast checks, the
# E2g and SD input checks, the SD sweeps, the capture and V-ABI quick. V-PERF (sd_perf.sh
# with N_AFTER and WA_AFTER) then measures these artifacts. Logs go into the run directory
# named by ~/.cache/losat-web-gui-target/sd/rundir.
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/sd
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $S/rundir)
NCBI=/home/kawato/micromamba/bin
export RUSTUP_TOOLCHAIN=1.92.0
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd $W
(cd LOSAT/tests && python3 -c 'import ci_fast_regressions as c; c.stage_lexical_fixtures()') || fail lexical-fixtures

mkdir -p $RUN
step "lint and tests"
git rev-parse HEAD > $RUN/head.txt
git status --short > $RUN/worktree-status.txt
python3 $A/s07p-resume/verify_refs.py $(git diff --name-only a92fa902f..HEAD -- '*.rs') > $RUN/verify-refs.log 2>&1; echo "exit $?" >> $RUN/verify-refs.log
python3 $W/docs/evidence/losat_web_e2i/gates/verify_added.py $RUN/verify-refs.log > $RUN/verify-refs-session-added.txt 2>&1; echo "exit $?" >> $RUN/verify-refs-session-added.txt
(cd LOSAT && cargo fmt --check) > $RUN/fmt.log 2>&1; echo "exit $?" >> $RUN/fmt.log
(cd web/adapter && cargo fmt --check) >> $RUN/fmt.log 2>&1; echo "adapter exit $?" >> $RUN/fmt.log
: > $RUN/clippy.log
for cfg in "--all-targets --all-features" "--all-targets --no-default-features" "--lib --target wasm32-wasip1 --no-default-features" "--lib --target wasm32-wasip1-threads --features wasm-threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/clippy.log
  (cd LOSAT && cargo clippy --locked $cfg --target-dir $A/sd-final-clippy -- -D warnings) >> $RUN/clippy.log 2>&1; echo "exit $?" >> $RUN/clippy.log
done
: > $RUN/adapter-clippy.log
for cfg in "--all-targets" "--lib --target wasm32-wasip1" "--lib --target wasm32-wasip1-threads --features threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/adapter-clippy.log
  (cd web/adapter && cargo clippy --locked $cfg --target-dir $A/sd-final-adapter -- -D warnings) >> $RUN/adapter-clippy.log 2>&1; echo "exit $?" >> $RUN/adapter-clippy.log
done
(cd LOSAT && LOSAT_BLASTX_WORKER_LOG=$A/blastx-worker.log cargo test --locked --all-features --no-fail-fast --target-dir $A/sd-final-test) > $RUN/cargo-test.log 2>&1; echo "exit $?" >> $RUN/cargo-test.log
(cd web/adapter && cargo test --locked --target-dir $A/sd-final-adapter -- --nocapture) > $RUN/adapter-test.log 2>&1; echo "exit $?" >> $RUN/adapter-test.log
(cd LOSAT && CARGO_TARGET_WASM32_WASIP1_RUNNER="node $W/web/tools/wasi-test-runner.mjs" cargo test --locked --lib --target wasm32-wasip1 --no-default-features --target-dir $A/sd-final-wasm32-test -- web_api::tests) > $RUN/wasm32-web-api-tests.log 2>&1; echo "exit $?" >> $RUN/wasm32-web-api-tests.log
{ python3 LOSAT/tests/check_pure_rust_runtime_boundary.py --root .; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_check_pure_rust_runtime_boundary.py'; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_ci_fast_regressions.py'; echo "exit $?"; } > $RUN/ci-python-checks.log 2>&1

step "build native, wasi and reactors"
(cd LOSAT && cargo build --release --locked --target-dir $A/sd-final-native) > $RUN/build-native.log 2>&1 || fail native
python3 LOSAT/tests/build_wasi_artifacts.py --target-dir $A/sd-final-wasi --output-dir $A/sd-final-wasi-artifacts --include-serial > $RUN/build-wasi.log 2>&1 || fail wasi
python3 web/adapter/tools/build_reactors.py --target-dir $A/sd-final-adapter --output-dir $A/sd-final-reactors > $RUN/build-reactors.log 2>&1 || fail reactors
mkdir -p $RUN/reactors && cp $A/sd-final-reactors/*.json $RUN/reactors/
python3 web/adapter/tools/check_build_identity.py --out $RUN/reactors/build-identity.json > /dev/null || fail identity
N=$A/sd-final-native/release/LOSAT
WA=$A/sd-final-wasi-artifacts
R=$A/sd-final-reactors
sha256sum $N $WA/*.wasm $R/*.wasm > $RUN/artifacts.sha256

step "fixtures"
for n in 1 2 4; do python3 docs/evidence/losat_web_e2a/check_losat.py --losat $N --threads $n > $RUN/check-losat-n$n.tsv 2>&1; echo "exit $?" >> $RUN/check-losat-n$n.tsv; done
python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir $NCBI --out $A/sd-final-oracle-check > $RUN/oracle-check-gate.log 2>&1; echo "exit $?" >> $RUN/oracle-check-gate.log
(cd LOSAT && python3 tests/blastn_regression_fixtures.py check --losat $N --jobs 8 --out $RUN/blastn-fixtures.tsv) > $RUN/blastn-fixtures.log 2>&1; echo "exit $?" >> $RUN/blastn-fixtures.log
(cd LOSAT && python3 tests/blastn_regression_fixtures.py check --losat $S/bin/native/LOSAT --jobs 8 --out $RUN/blastn-fixtures-before.tsv) > $RUN/blastn-fixtures-before.log 2>&1; echo "exit $?" >> $RUN/blastn-fixtures-before.log
(cd LOSAT && python3 tests/tblastx_regression_fixtures.py check --losat $N --jobs 8 --out $RUN/tblastx-fixtures.tsv) > $RUN/tblastx-fixtures.log 2>&1; echo "exit $?" >> $RUN/tblastx-fixtures.log
python3 docs/evidence/losat_web_e2i/check_authority.py > $RUN/check-authority.log 2>&1; echo "exit $?" >> $RUN/check-authority.log

step "input checks (E2g and SD) and SD sweeps"
rm -rf $A/sd-final-inputs
python3 docs/evidence/losat_web_e2g/check_inputs.py --bin-dir $NCBI --losat $N --work $A/sd-final-inputs > $RUN/check-inputs.tsv 2>&1; echo "exit $?" >> $RUN/check-inputs.tsv
rm -rf $A/sd-final-sweeps
bash docs/evidence/losat_web_e2i/sweeps.sh $N $NCBI $RUN/sd-sweeps $A/sd-final-sweeps > $RUN/sd-sweeps.log 2>&1; echo "exit $?" >> $RUN/sd-sweeps.log

step "ci fast checks"
rm -rf $A/sd-final-fast
python3 LOSAT/tests/ci_fast_regressions.py --losat $N --out $A/sd-final-fast --all-cases --jobs 6 > $RUN/fast-regressions-all.log 2>&1; echo "exit $?" >> $RUN/fast-regressions-all.log
cp $A/sd-final-fast/summary.json $RUN/fast-regressions-all-summary.json 2>/dev/null

step "capture"
rm -rf $A/capture-sd-final
python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat $N --out $A/capture-sd-final --jobs 6 > $RUN/capture.log 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare docs/evidence/losat_web_e1a/baseline/hashes.tsv $A/capture-sd-final/hashes.tsv > $RUN/capture-compare.txt 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare $S/capture-base/hashes.tsv $A/capture-sd-final/hashes.tsv > $RUN/capture-compare-before.txt 2>&1
mkdir -p $RUN/capture && cp $A/capture-sd-final/hashes.tsv $A/capture-sd-final/binary.json $RUN/capture/ 2>/dev/null

step "v-abi quick"
mkdir -p $RUN/v-abi-quick
python3 web/adapter/tools/v_abi_cases.py --suite quick --out $RUN/v-abi-quick/cases.json > /dev/null
node web/adapter/tests/v_abi.js --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-quick/cases.json --out $RUN/v-abi-quick > $RUN/v-abi-quick/v-abi.log 2>&1
echo "exit $?" >> $RUN/v-abi-quick/v-abi.log

step "done (V-PERF separately)"

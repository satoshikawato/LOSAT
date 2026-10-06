#!/usr/bin/env bash
# SD (E2i, BLASTN dc-megablast and blastn-short) gate runs: lint, tests, builds, then the
# gates; logs go into the SD run directory (docs/evidence/losat_web_e2i/run-*). Made from
# docs/evidence/losat_web_e2b/gates/s08_gates.sh and E2g's BLASTN sweeps.
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
# The Gate A lexical root that frozen hashes name (/tmp is emptied when WSL restarts):
# staged as the CI fast job stages it (ci_fast_regressions.py stage_lexical_fixtures).
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
  (cd LOSAT && cargo clippy --locked $cfg --target-dir $A/sd-gate-clippy -- -D warnings) >> $RUN/clippy.log 2>&1; echo "exit $?" >> $RUN/clippy.log
done
: > $RUN/adapter-clippy.log
for cfg in "--all-targets" "--lib --target wasm32-wasip1" "--lib --target wasm32-wasip1-threads --features threads"; do
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $RUN/adapter-clippy.log
  (cd web/adapter && cargo clippy --locked $cfg --target-dir $A/sd-gate-adapter -- -D warnings) >> $RUN/adapter-clippy.log 2>&1; echo "exit $?" >> $RUN/adapter-clippy.log
done
(cd LOSAT && LOSAT_BLASTX_WORKER_LOG=$A/blastx-worker.log cargo test --locked --all-features --no-fail-fast --target-dir $A/sd-gate-test) > $RUN/cargo-test.log 2>&1; echo "exit $?" >> $RUN/cargo-test.log
(cd web/adapter && cargo test --locked --target-dir $A/sd-gate-adapter -- --nocapture) > $RUN/adapter-test.log 2>&1; echo "exit $?" >> $RUN/adapter-test.log
(cd LOSAT && CARGO_TARGET_WASM32_WASIP1_RUNNER="node $W/web/tools/wasi-test-runner.mjs" cargo test --locked --lib --target wasm32-wasip1 --no-default-features --target-dir $A/sd-gate-wasm32-test -- web_api::tests) > $RUN/wasm32-web-api-tests.log 2>&1; echo "exit $?" >> $RUN/wasm32-web-api-tests.log
{ python3 LOSAT/tests/check_pure_rust_runtime_boundary.py --root .; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_check_pure_rust_runtime_boundary.py'; echo "exit $?";
  python3 -m unittest discover -s LOSAT/tests -p 'test_ci_fast_regressions.py'; echo "exit $?"; } > $RUN/ci-python-checks.log 2>&1

step "build native"
(cd LOSAT && cargo build --release --locked --target-dir $A/sd-gate-native) > $RUN/build-native.log 2>&1 || fail native
(cd LOSAT && cargo build --release --locked --no-default-features --target-dir $A/sd-gate-native-serial) > $RUN/build-native-serial.log 2>&1 || fail native-serial
step "build wasi"
python3 LOSAT/tests/build_wasi_artifacts.py --target-dir $A/sd-gate-wasi --output-dir $A/sd-gate-wasi-artifacts --include-serial > $RUN/build-wasi.log 2>&1 || fail wasi
step "build reactors"
python3 web/adapter/tools/build_reactors.py --target-dir $A/sd-gate-adapter --output-dir $A/sd-gate-reactors > $RUN/build-reactors.log 2>&1 || fail reactors
mkdir -p $RUN/reactors && cp $A/sd-gate-reactors/*.json $RUN/reactors/
python3 web/adapter/tools/check_build_identity.py --out $RUN/reactors/build-identity.json > /dev/null || fail identity

N=$A/sd-gate-native/release/LOSAT
WA=$A/sd-gate-wasi-artifacts
R=$A/sd-gate-reactors
sha256sum $N $A/sd-gate-native-serial/release/LOSAT $WA/*.wasm $R/*.wasm > $RUN/artifacts.sha256

step "v-abi full (background)"
mkdir -p $RUN/v-abi-full $RUN/v-abi-quick
python3 web/adapter/tools/v_abi_cases.py --suite full --out $RUN/v-abi-full/cases.json > /dev/null || fail cases-full
python3 web/adapter/tools/run_v_abi_parallel.py --native $N --serial $R/losat-web-serial.wasm --threads $R/losat-web-threads.wasm \
  --cases $RUN/v-abi-full/cases.json --out $RUN/v-abi-full --jobs 8 > $RUN/v-abi-full/run.log 2>&1 &
VABI=$!

step "fixtures"
for n in 1 2 4; do python3 docs/evidence/losat_web_e2a/check_losat.py --losat $N --threads $n > $RUN/check-losat-n$n.tsv 2>&1; echo "exit $?" >> $RUN/check-losat-n$n.tsv; done
python3 docs/evidence/losat_web_e2a/run_oracle.py --bin-dir $NCBI --out $A/sd-gate-oracle-check > $RUN/oracle-check-gate.log 2>&1; echo "exit $?" >> $RUN/oracle-check-gate.log
(cd LOSAT && python3 tests/blastn_regression_fixtures.py check --losat $S/bin/native/LOSAT --jobs 8 --out $RUN/blastn-fixtures-before.tsv) > $RUN/blastn-fixtures-before.log 2>&1; echo "exit $?" >> $RUN/blastn-fixtures-before.log
(cd LOSAT && python3 tests/tblastx_regression_fixtures.py check --losat $N --jobs 8 --out $RUN/tblastx-fixtures.tsv) > $RUN/tblastx-fixtures.log 2>&1; echo "exit $?" >> $RUN/tblastx-fixtures.log
python3 docs/evidence/losat_web_e2b/ctoolkit_compare.py --bin-dir $NCBI --losat $N --jobs 8 > $RUN/ctoolkit-compare.tsv 2>&1; echo "exit $?" >> $RUN/ctoolkit-compare.tsv
python3 docs/evidence/losat_web_e2i/check_authority.py > $RUN/check-authority.log 2>&1; echo "exit $?" >> $RUN/check-authority.log
step "blastn sweeps (E2c, E2f, E2g; megablast and blastn)"
(cd docs/evidence/losat_web_e2a && python3 precheck_hits.py --bin-dir $NCBI --losat $N > $RUN/precheck.tsv 2>&1)
(cd docs/evidence/losat_web_e2a && python3 scoring_sweep.py --bin-dir $NCBI --losat $N > $RUN/scoring-sweep-s06.tsv 2>&1)
for f in 0 6 7; do python3 docs/evidence/losat_web_e2c/scoring_sweep.py --bin-dir $NCBI --losat $N --jobs 8 --outfmt $f > $RUN/scoring-sweep-fmt$f.tsv 2>&1; echo "exit $?" >> $RUN/scoring-sweep-fmt$f.tsv; done
python3 docs/evidence/losat_web_e2c/word_size_sweep.py --bin-dir $NCBI --losat $N --jobs 8 > $RUN/word-size-sweep.tsv 2>&1; echo "exit $?" >> $RUN/word-size-sweep.tsv
for opts in "" "-task blastn" "-task blastn -subject_besthit" "-subject_besthit" "-word_size 8" "-reward 3 -penalty -4 -gapopen 10 -gapextend 3"; do
  name=$(echo "default $opts" | tr ' ' '_' | tr -d '-'); rm -rf $A/sd-gate-slices
  python3 docs/evidence/losat_web_e2c/slice_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-slices --options="$opts" --jobs 8 > $RUN/slice-sweep-$name.tsv 2>&1; echo "exit $?" >> $RUN/slice-sweep-$name.tsv
done
for opts in "" "-task blastn" "-subject_besthit" "-task blastn -word_size 7"; do
  name=$(echo "ambiguity $opts" | tr ' ' '_' | tr -d '-'); rm -rf $A/sd-gate-slices
  python3 docs/evidence/losat_web_e2c/slice_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-slices --pool ambiguity --options="$opts" --cases 100 --jobs 8 > $RUN/slice-sweep-$name.tsv 2>&1; echo "exit $?" >> $RUN/slice-sweep-$name.tsv
done
for opts in "-lcase_masking" "-task blastn -lcase_masking" "-task blastn -word_size 7 -lcase_masking"; do
  name=$(echo "lcase $opts" | tr ' ' '_' | tr -d '-'); rm -rf $A/sd-gate-slices
  python3 docs/evidence/losat_web_e2c/slice_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-slices --pool lcase --options="$opts" --jobs 8 > $RUN/slice-sweep-$name.tsv 2>&1; echo "exit $?" >> $RUN/slice-sweep-$name.tsv
done
rm -rf $A/sd-gate-slices
python3 docs/evidence/losat_web_e2c/slice_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-slices --pool ambiguity --options="-task blastn -word_size 4 -evalue 1e6" --cases 40 --jobs 8 > $RUN/slice-sweep-ambiguity_task_blastn_word_size_4_evalue_1e6.tsv 2>&1; echo "exit $?" >> $RUN/slice-sweep-ambiguity_task_blastn_word_size_4_evalue_1e6.tsv
for opts in "" "-word_size 24" "-task blastn"; do
  name=$(echo "iupac $opts" | tr ' ' '_' | tr -d '-'); rm -rf $A/sd-gate-slices
  python3 docs/evidence/losat_web_e2c/slice_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-slices --pool iupac --options="$opts" --cases 300 --jobs 8 > $RUN/slice-sweep-$name.tsv 2>&1; echo "exit $?" >> $RUN/slice-sweep-$name.tsv
done
rm -rf $A/sd-gate-titles
python3 docs/evidence/losat_web_e2g/title_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-titles --jobs 8 > $RUN/title-sweep.tsv 2>&1; echo "exit $?" >> $RUN/title-sweep.tsv
for opts in "" "-task blastn" "-subject_besthit" "-task blastn -word_size 7" "-task blastn -subject_besthit -max_target_seqs 3"; do
  name=$(echo "batches $opts" | tr ' ' '_' | tr -d '-'); rm -rf $A/sd-gate-bt
  python3 docs/evidence/losat_web_e2f/batch_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-bt --options="$opts" --cases 150 --jobs 8 > $RUN/batch-sweep-$name.tsv 2>&1; echo "exit $?" >> $RUN/batch-sweep-$name.tsv
done
rm -rf $A/sd-gate-split
python3 docs/evidence/losat_web_e2f/split_check.py --bin-dir $NCBI --losat $N --work $A/sd-gate-split --jobs 6 > $RUN/split-check.tsv 2>&1; echo "exit $?" >> $RUN/split-check.tsv
for env in "BATCH_SIZE=1000" "BATCH_SIZE=100000" "CHUNK_SIZE=40000" "CHUNK_SIZE=2000 OVERLAP_CHUNK_SIZE=50"; do
  name=$(echo "batches env $env" | tr ' =' '__'); rm -rf $A/sd-gate-bt
  env $env python3 docs/evidence/losat_web_e2f/batch_sweep.py --bin-dir $NCBI --losat $N --work $A/sd-gate-bt --options="-task blastn" --cases 60 --jobs 8 > $RUN/batch-sweep-$name.tsv 2>&1; echo "exit $?" >> $RUN/batch-sweep-$name.tsv
done
(cd LOSAT && python3 tests/blastn_regression_fixtures.py check --losat $N --jobs 8 --out $RUN/blastn-fixtures.tsv) > $RUN/blastn-fixtures.log 2>&1; echo "exit $?" >> $RUN/blastn-fixtures.log
rm -rf $A/sd-gate-fast
python3 LOSAT/tests/ci_fast_regressions.py --losat $N --out $A/sd-gate-fast --all-cases --jobs 6 > $RUN/fast-regressions-all.log 2>&1; echo "exit $?" >> $RUN/fast-regressions-all.log
cp $A/sd-gate-fast/summary.json $RUN/fast-regressions-all-summary.json 2>/dev/null
rm -rf $A/sd-gate-inputs
python3 docs/evidence/losat_web_e2g/check_inputs.py --bin-dir $NCBI --losat $N --work $A/sd-gate-inputs > $RUN/check-inputs.tsv 2>&1; echo "exit $?" >> $RUN/check-inputs.tsv
step "sd sweeps (dc-megablast and blastn-short)"
rm -rf $A/sd-gate-sweeps
bash docs/evidence/losat_web_e2i/sweeps.sh $N $NCBI $RUN/sd-sweeps $A/sd-gate-sweeps > $RUN/sd-sweeps.log 2>&1; echo "exit $?" >> $RUN/sd-sweeps.log
rm -rf $A/sd-gate-fast
python3 LOSAT/tests/ci_fast_regressions.py --losat $N --out $A/sd-gate-fast --all-cases --jobs 6 > $RUN/fast-regressions-all.log 2>&1; echo "exit $?" >> $RUN/fast-regressions-all.log
cp $A/sd-gate-fast/summary.json $RUN/fast-regressions-all-summary.json 2>/dev/null

step "capture"
rm -rf $A/capture-sd
python3 docs/evidence/losat_web_e1a/capture_outputs.py run --losat $N --out $A/capture-sd --jobs 6 > $RUN/capture.log 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare docs/evidence/losat_web_e1a/baseline/hashes.tsv $A/capture-sd/hashes.tsv > $RUN/capture-compare.txt 2>&1
python3 docs/evidence/losat_web_e1a/capture_outputs.py compare $S/capture-base/hashes.tsv $A/capture-sd/hashes.tsv > $RUN/capture-compare-before.txt 2>&1
mkdir -p $RUN/capture && cp $A/capture-sd/hashes.tsv $A/capture-sd/binary.json $RUN/capture/ 2>/dev/null

step "v1 wasi matrix"
rm -rf $A/wasm-threading-sd
python3 LOSAT/tests/check_wasm_threading.py --native $N --native-serial $A/sd-gate-native-serial/release/LOSAT \
  --serial $WA/losat-serial-command.wasm --threaded $WA/losat-threaded-command.wasm \
  --reactor $WA/losat-threaded-reactor.wasm --serial-reactor $WA/losat-serial-reactor.wasm \
  --oracle-dir $NCBI --output-dir $A/wasm-threading-sd > $RUN/wasm-threading.log 2>&1
echo "exit $?" >> $RUN/wasm-threading.log
cp $A/wasm-threading-sd/metadata.json $RUN/wasm-threading-metadata.json 2>/dev/null
cp $A/wasm-threading-sd/runs.json $RUN/wasm-threading-runs.json 2>/dev/null
{
  echo "# v1 reactor records (check_wasi_reactor.js, run by check_wasm_threading.py) of the S05 and SD reactors."
  echo '$ diff -r -x runs.json s05/reactor-records sd/reactor-records'
  diff -r -x runs.json $A/wasm-threading-s05/reactor-records $A/wasm-threading-sd/reactor-records; echo "exit $?"
  echo "records compared: $(ls $A/wasm-threading-sd/reactor-records | grep -c '\.out$')"
  echo '$ diff -r -x runs.json s05/serial-reactor-records sd/serial-reactor-records'
  diff -r -x runs.json $A/wasm-threading-s05/serial-reactor-records $A/wasm-threading-sd/serial-reactor-records; echo "exit $?"
  echo "records compared: $(ls $A/wasm-threading-sd/serial-reactor-records | grep -c '\.out$')"
} > $RUN/v1-reactor-records-compare.txt 2>&1
node docs/evidence/losat_web_e1a/v1_requests.js $WA/losat-threaded-reactor.wasm $A/wasm-threading-sd/fixtures/aa3.fasta > $RUN/v1-requests-after.jsonl 2>/dev/null
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

#!/usr/bin/env bash
# S08+a lint and tests on the HEAD export (made from docs/evidence/losat_web_e2e/gates/s08p_gates.sh).
set -u
A=/home/kawato/.cache/losat-web-gui-target
H=$A/s08pa/head-src
OUT=$A/s08pa/checks
BASE=78c06fe61
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui-s08pa
export RUSTUP_TOOLCHAIN=1.92.0
waitlock() { while [ -e $A/vperf.lock ]; do echo "$(date -u +%H:%M:%S) vperf.lock present, waiting" >> $OUT/progress.log; sleep 60; done; }
step() { waitlock; echo "$(date -u +%H:%M:%S) $*" >> $OUT/progress.log; }
step verify-refs
(cd $W && python3 $A/s08p/verify_refs.py $(git diff --name-only --diff-filter=d $BASE..HEAD -- '*.rs')) > $OUT/verify-refs.log 2>&1; echo "exit $?" >> $OUT/verify-refs.log
(cd $W && python3 $H/docs/evidence/losat_web_e2e/gates/verify_added.py $OUT/verify-refs.log) > $OUT/verify-refs-session-added.txt 2>&1; echo "exit $?" >> $OUT/verify-refs-session-added.txt
cd $H
step fmt
(cd LOSAT && cargo fmt --check) > $OUT/fmt.log 2>&1; echo "exit $?" >> $OUT/fmt.log
(cd web/adapter && cargo fmt --check) >> $OUT/fmt.log 2>&1; echo "adapter exit $?" >> $OUT/fmt.log
: > $OUT/clippy.log
for cfg in "--all-targets --all-features" "--all-targets --no-default-features" "--lib --target wasm32-wasip1 --no-default-features" "--lib --target wasm32-wasip1-threads --features wasm-threads"; do
  step "clippy $cfg"
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $OUT/clippy.log
  (cd LOSAT && cargo clippy --locked $cfg --target-dir $A/s08pa-clippy -- -D warnings) >> $OUT/clippy.log 2>&1; echo "exit $?" >> $OUT/clippy.log
done
: > $OUT/adapter-clippy.log
for cfg in "--all-targets" "--lib --target wasm32-wasip1" "--lib --target wasm32-wasip1-threads --features threads"; do
  step "adapter clippy $cfg"
  echo "== cargo clippy --locked $cfg -- -D warnings" >> $OUT/adapter-clippy.log
  (cd web/adapter && cargo clippy --locked $cfg --target-dir $A/s08pa-adapter -- -D warnings) >> $OUT/adapter-clippy.log 2>&1; echo "exit $?" >> $OUT/adapter-clippy.log
done
step cargo-test
(cd LOSAT && LOSAT_BLASTX_WORKER_LOG=$A/s08pa/blastx-worker.log CARGO_PROFILE_TEST_OPT_LEVEL=1 cargo test --locked --all-features --no-fail-fast --target-dir $A/s08pa-gate-test) > $OUT/cargo-test.log 2>&1; echo "exit $?" >> $OUT/cargo-test.log
step adapter-test
(cd web/adapter && cargo test --locked --target-dir $A/s08pa-adapter -- --nocapture) > $OUT/adapter-test.log 2>&1; echo "exit $?" >> $OUT/adapter-test.log
step done

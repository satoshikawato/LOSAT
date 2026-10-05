#!/usr/bin/env bash
set -u
A=/home/kawato/.cache/losat-web-gui-target
P=$A/s08pb
SRC=$P/src-merge
export RUSTUP_TOOLCHAIN=1.92.0
echo "$(date -u +%T) build"
(cd $SRC/LOSAT && cargo build --release --locked --target-dir $A/s08pb-native) > $P/merge-build.log 2>&1 || { echo build failed; exit 1; }
N=$A/s08pb-native/release/LOSAT
sha256sum $N > $P/merge-native.sha256
echo "$(date -u +%T) fixtures"
cd $SRC
for n in 1 2 4; do python3 docs/evidence/losat_web_e2a/check_losat.py --losat $N --threads $n > $P/merge-check-n$n.tsv 2>&1; echo "exit $?" >> $P/merge-check-n$n.tsv; done
echo "$(date -u +%T) cargo test"
(cd $SRC/LOSAT && LOSAT_BLASTX_WORKER_LOG=$A/blastx-worker.log CARGO_PROFILE_TEST_OPT_LEVEL=1 cargo test --locked --all-features --no-fail-fast --target-dir $A/s08pb-test) > $P/merge-cargo-test.log 2>&1; echo "exit $?" >> $P/merge-cargo-test.log
(cd $SRC/web/adapter && cargo test --locked --target-dir $A/s08pb-adapter) > $P/merge-adapter-test.log 2>&1; echo "exit $?" >> $P/merge-adapter-test.log
echo "$(date -u +%T) done"

#!/usr/bin/env bash
# App-side verification of the app worktree: npm run check, the FakeEngine E2E, and the
# engine-build E2E with the new (65,535 pages) and the baseline (16,384 pages) reactors.
# Each step holds the machine-wide lock. Logs: $T/logs/app/. Usage: verify_app.sh [steps...]
set -uo pipefail
T=/home/kawato/losat-baselines/wasm-mem-20261012
APP=/home/kawato/losat-work/.worktrees/wasm-mem-app/web/app
BR=/home/kawato/.cache/losat-work/wasm-mem
LOCK=/home/kawato/.cache/losat-work/oracle.lock
export PATH=/home/kawato/.local/bin:$PATH TMPDIR=$T/tmp LOSAT_WEB_E2E_PORT=4391
export LOSAT_WEB_WEBKIT_EXECUTABLE=/home/kawato/.cache/losat-web-gui-target/app-s09-webkit-deps/pw_run.sh
unset NODE_ENV VITE_USER_NODE_ENV LOSAT_WEB_REACTORS LOSAT_WEB_NATIVE LOSAT_WEB_EVIDENCE
L=$T/logs/app; mkdir -p "$L"
steps=${*:-"check e2e_fake e2e_new e2e_base"}
cd "$APP"
status() { printf '%s %s rc=%s (%s)\n' "$(date +%T)" "$1" "$2" "$3" | tee -a "$L/summary.txt"; }
tally() { grep -E '^\s+[0-9]+ (passed|failed|skipped|flaky|did not run)' "$1" | tr -s ' ' | tr '\n' ';'; }
for step in $steps; do case $step in
  check)
    flock "$LOCK" npm run check > "$L/check.log" 2>&1; rc=$?
    status check $rc "$(grep -E 'Test Files|Tests  ' "$L/check.log" | tr -s ' ' | tr '\n' ';')";;
  e2e_fake)
    flock "$LOCK" npm run e2e > "$L/e2e_fake.log" 2>&1; rc=$?
    status e2e_fake $rc "$(tally "$L/e2e_fake.log")";;
  e2e_new|e2e_base)
    v=${step#e2e_}; mkdir -p "$L/evidence-$v"
    LOSAT_WEB_REACTORS=$BR/out-$v/adapter LOSAT_WEB_NATIVE=$BR/target-native/release/LOSAT LOSAT_WEB_EVIDENCE=$L/evidence-$v \
      flock "$LOCK" npm run e2e > "$L/$step.log" 2>&1; rc=$?
    status $step $rc "$(tally "$L/$step.log")";;
esac; done
status done 0 "$steps"

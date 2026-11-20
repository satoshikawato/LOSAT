#!/bin/bash
# NP WW0: the WebKit project of the engine-build and FakeEngine E2E, with the launcher of webkit-deps (the gate's
# first launcher did not reach the bundle's MiniBrowser). Logs next to the gate's: logs/gate_<head>/.
set -u
BR=/home/kawato/.cache/losat-work
W=/home/kawato/losat-work/.worktrees/np-ww0
T=/home/kawato/losat-baselines/ncbi-parity-20261019/WW
R=$BR/np-ww0
export TMPDIR=$T/tmp PATH=/home/kawato/.local/bin:$PATH LOSAT_WEB_WEBKIT_EXECUTABLE=$R/webkit-deps/pw_run.sh LOSAT_WEB_E2E_PORT=4397
for v in $(env | grep -o '^LOSAT_[A-Z0-9_]*' | grep -v 'LOSAT_WEB_WEBKIT_EXECUTABLE\|LOSAT_WEB_E2E_PORT'); do unset $v; done
cd "$W/web/app" || exit 1
HEAD=$(git rev-parse --short HEAD)
G=$T/logs/gate_$HEAD; S=$G/summary.txt; mkdir -p "$G" "$T/results/evidence-$HEAD"
LOCK="flock -w 3600 $BR/oracle.lock"
tally() { grep -E '^\s+[0-9]+ (passed|failed|skipped|flaky|did not run)' "$1" | tr -s ' ' | tr '\n' ';'; }
LOSAT_WEB_REACTORS="$R/reactors" LOSAT_WEB_NATIVE="$R/native/losat" LOSAT_WEB_EVIDENCE="$T/results/evidence-$HEAD" \
  $LOCK npx playwright test --project=webkit > "$G/e2e-engine-webkit.log" 2>&1
echo "$(date +%H:%M:%S) e2e-engine-webkit exit=$? | $(tally "$G/e2e-engine-webkit.log")" | tee -a "$S"
$LOCK npx playwright test --project=webkit > "$G/e2e-fake-webkit.log" 2>&1
echo "$(date +%H:%M:%S) e2e-fake-webkit exit=$? | $(tally "$G/e2e-fake-webkit.log")" | tee -a "$S"

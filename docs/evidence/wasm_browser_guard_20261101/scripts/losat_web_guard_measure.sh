#!/usr/bin/env bash
# Run guard-measure.spec.ts (Chromium) with the reactors of $REACTORS; evidence to $1.
set -uo pipefail
out=$(realpath -m "${1:?evidence dir}"); shift
mkdir -p "$out"
cd /home/kawato/losat-work/.worktrees/wg-guard/web/app
export LOSAT_WEB_REACTORS=${REACTORS:-/home/kawato/.cache/losat-work/webapp/reactors}
export LOSAT_WEB_MEASURE=guard LOSAT_WEB_EVIDENCE="$out" LOSAT_WEB_E2E_PORT=${PORT:-4391}
export TMPDIR=/home/kawato/losat-baselines/wg-wasm-guard-20261101/tmp
mkdir -p "$TMPDIR"
echo "start $(date -Is) load $(cut -d' ' -f1-3 /proc/loadavg)"
npx playwright test tests/e2e/guard-measure.spec.ts --project chromium --workers 1 "$@"
echo "exit $? end $(date -Is) load $(cut -d' ' -f1-3 /proc/loadavg)"

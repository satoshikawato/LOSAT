#!/usr/bin/env bash
# Gate run of LOSAT Web W2 (Session S10). It records the environment, then runs
#   npm ci, npm run check, the unit tests with every case listed, npm run e2e, and the
#   storage and contract E2E specs three more times (stability),
# and writes the logs to docs/evidence/losat_web_w2/run-<UTC>/, which is not changed later.
# Before each step it waits while the engine track measures V-PERF (plan DW-7).
# Usage: docs/evidence/losat_web_w2/run_gate.sh
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
app="$repo/web/app"
lock="${LOSAT_VPERF_LOCK:-/home/kawato/.cache/losat-web-gui-target/vperf.lock}"

wait_for_vperf() {
  until [ ! -e "$lock" ]; do
    echo "waiting: V-PERF is being measured ($lock)" >&2
    sleep 60
  done
}

if curl -s -o /dev/null "http://localhost:4173/"; then
  echo "port 4173 is in use; the E2E tests would reuse another server" >&2
  exit 1
fi

run="$here/run-$(date -u +%Y%m%dT%H%M%SZ)"
mkdir "$run"
{
  echo "commit: $(git -C "$repo" rev-parse HEAD)"
  echo "branch: $(git -C "$repo" branch --show-current)"
  echo "uncommitted files under web/app: $(git -C "$repo" status --porcelain -- web/app | wc -l)"
  echo "node: $(node --version)"
  echo "npm: $(npm --version)"
  echo "playwright: $(cd "$app" && npx playwright --version 2>/dev/null || echo unknown)"
  echo "playwright browsers: $(ls "${PLAYWRIGHT_BROWSERS_PATH:-$HOME/.cache/ms-playwright}" 2>/dev/null | tr '\n' ' ')"
  echo "os: $(uname -srm)"
} > "$run/environment.txt"

cd "$app"
wait_for_vperf
npm ci > "$run/npm-ci.log" 2>&1
wait_for_vperf
npm run check > "$run/npm-check.log" 2>&1
wait_for_vperf
npx vitest run --reporter=verbose > "$run/unit-cases.log" 2>&1
wait_for_vperf
npm run e2e -- --reporter=list > "$run/npm-e2e.log" 2>&1
wait_for_vperf
npx playwright test tests/e2e/storage.spec.ts tests/e2e/contracts.spec.ts --repeat-each 3 --reporter=list \
  > "$run/e2e-repeat.log" 2>&1

echo "gate run passed: $run"

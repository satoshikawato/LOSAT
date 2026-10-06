#!/usr/bin/env bash
# Gate run of LOSAT Web W3 (Session S12, the search screen). It records the environment and
# the engine build, then runs
#   V-ABI quick (Node) on the reactors of this commit, npm ci, npm run check, the unit
#   tests with every case listed (with the reactors: the describe and input-check cases on
#   the serial reactor in Node), npm run e2e without the engine (FakeEngine build) and with
#   it (the search screen's research work and boundaries, S09's V-BR and contracts;
#   Chromium, Firefox and WebKit), the search screen's E2E twice more with the engine
#   (stability), the measurements of many records (tests/e2e/search-measure.spec.ts) in
#   Chromium, Firefox and WebKit, and the screen records for the visual review
#   (tests/e2e/screens.spec.ts) in the three browsers.
# Logs and records go to docs/evidence/losat_web_w3/run-<UTC>/, which is not changed later;
# the screen records (PNG) go to $LOSAT_WEB_SCREENS_OUT (outside the repository; the
# README lists their SHA-256).
# Before each step it waits while the engine track measures V-PERF (plan DW-7).
#
# Needs:
#   LOSAT_WEB_REACTORS  the output directory of web/adapter/tools/build_reactors.py for this commit
#   LOSAT_WEB_NATIVE    the native LOSAT of this commit (cargo +1.92.0 build --release)
#   LOSAT_WEB_SCREENS_OUT  where the screen records go
# Optional:
#   LOSAT_WEB_E2E_PORT           the port of the application under test (4173)
#   LOSAT_WEB_WEBKIT_EXECUTABLE  a launcher for Playwright's WebKit on a host whose system
#                                libraries Playwright cannot install (W1 README)
#   LOSAT_WEB_MEASURE_RECORDS    record counts to measure (1000,10000,100000)
#   LOSAT_WEB_GATE_STEPS         `all` (the gate), or `after-review`: npm ci, check, the unit
#                                tests, both E2E builds and the screen records only (the
#                                run after a fix of a review finding; README)
# Usage: docs/evidence/losat_web_w3/run_gate.sh
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
app="$repo/web/app"
lock="${LOSAT_VPERF_LOCK:-/home/kawato/.cache/losat-web-gui-target/vperf.lock}"
: "${LOSAT_WEB_REACTORS:?set LOSAT_WEB_REACTORS to the reactors of this commit}"
: "${LOSAT_WEB_NATIVE:?set LOSAT_WEB_NATIVE to the native LOSAT of this commit}"
: "${LOSAT_WEB_SCREENS_OUT:?set LOSAT_WEB_SCREENS_OUT to a directory for the screen records}"
reactors="$(cd "$LOSAT_WEB_REACTORS" && pwd)"
native="$(cd "$(dirname "$LOSAT_WEB_NATIVE")" && pwd)/$(basename "$LOSAT_WEB_NATIVE")"
screens="$LOSAT_WEB_SCREENS_OUT"
records="${LOSAT_WEB_MEASURE_RECORDS:-1000,10000,100000}"
steps="${LOSAT_WEB_GATE_STEPS:-all}"
case "$steps" in all | after-review) ;; *) echo "LOSAT_WEB_GATE_STEPS is all or after-review" >&2; exit 2 ;; esac
# Each step sets the engine variables that it needs; the others build with the FakeEngine.
unset LOSAT_WEB_REACTORS LOSAT_WEB_NATIVE LOSAT_WEB_EVIDENCE LOSAT_WEB_MEASURE LOSAT_WEB_SCREENS LOSAT_WEB_SCREENS_OUT

wait_for_vperf() {
  until [ ! -e "$lock" ]; do
    echo "waiting: V-PERF is being measured ($lock)" >&2
    sleep 60
  done
}

port="${LOSAT_WEB_E2E_PORT:-4173}"
if curl -s -o /dev/null "http://localhost:$port/"; then
  echo "port $port is in use (set LOSAT_WEB_E2E_PORT to another port)" >&2
  exit 1
fi

run="$here/run-$(date -u +%Y%m%dT%H%M%SZ)"
[ "$steps" = all ] || run="$run-$steps"
mkdir "$run"
{
  echo "steps: $steps"
  echo "commit: $(git -C "$repo" rev-parse HEAD)"
  echo "branch: $(git -C "$repo" branch --show-current)"
  echo "uncommitted files under web/app: $(git -C "$repo" status --porcelain -- web/app | wc -l)"
  echo "node: $(node --version)"
  echo "npm: $(npm --version)"
  echo "playwright: $(cd "$app" && npx playwright --version 2>/dev/null || echo unknown)"
  echo "playwright browsers: $(ls "${PLAYWRIGHT_BROWSERS_PATH:-$HOME/.cache/ms-playwright}" 2>/dev/null | tr '\n' ' ')"
  echo "webkit launcher: ${LOSAT_WEB_WEBKIT_EXECUTABLE:-playwright default}"
  echo "os: $(uname -srm)"
  echo "processors: $(nproc)"
  echo "native: $native sha256 $(sha256sum "$native" | cut -d' ' -f1)"
  echo "reactors: $reactors"
  (cd "$reactors" && sha256sum losat-web-serial.wasm losat-web-threads.wasm)
  echo "screens: $screens"
} > "$run/environment.txt"
cp "$reactors/artifacts.json" "$run/reactors-artifacts.json"

if [ "$steps" = all ]; then
  wait_for_vperf
  python3 "$repo/web/adapter/tools/v_abi_cases.py" --suite quick --out "$run/v-abi-cases.json"
  node "$repo/web/adapter/tests/v_abi.js" --native "$native" \
    --serial "$reactors/losat-web-serial.wasm" --threads "$reactors/losat-web-threads.wasm" \
    --cases "$run/v-abi-cases.json" --out "$run" > "$run/v-abi.log" 2>&1
fi

cd "$app"
wait_for_vperf
npm ci > "$run/npm-ci.log" 2>&1
wait_for_vperf
npm run check > "$run/npm-check.log" 2>&1
wait_for_vperf
LOSAT_WEB_REACTORS="$reactors" npx vitest run --reporter=verbose > "$run/unit-cases.log" 2>&1
wait_for_vperf
npm run e2e -- --reporter=list > "$run/npm-e2e-fake-engine.log" 2>&1
wait_for_vperf
mkdir "$run/records"
LOSAT_WEB_REACTORS="$reactors" LOSAT_WEB_NATIVE="$native" LOSAT_WEB_EVIDENCE="$run/records" \
  npm run e2e -- --reporter=list > "$run/npm-e2e-engine.log" 2>&1
if [ "$steps" = all ]; then
  wait_for_vperf
  LOSAT_WEB_REACTORS="$reactors" \
    npx playwright test tests/e2e/search.spec.ts --repeat-each 2 --reporter=list > "$run/e2e-search-repeat.log" 2>&1

  mkdir -p "$run/measure"
  for project in chromium firefox webkit; do
    wait_for_vperf
    LOSAT_WEB_MEASURE=records LOSAT_WEB_MEASURE_RECORDS="$records" LOSAT_WEB_REACTORS="$reactors" \
      LOSAT_WEB_EVIDENCE="$run/measure" \
      npx playwright test tests/e2e/search-measure.spec.ts --project="$project" --workers=1 --reporter=list \
      > "$run/measure/$project.log" 2>&1
  done
fi

wait_for_vperf
mkdir -p "$screens"
LOSAT_WEB_SCREENS="$screens" LOSAT_WEB_REACTORS="$reactors" \
  npx playwright test tests/e2e/screens.spec.ts --reporter=list > "$run/screens.log" 2>&1
(cd "$screens" && find . -name '*.png' | sort | xargs sha256sum) > "$run/screens.sha256"

echo "gate run passed: $run"

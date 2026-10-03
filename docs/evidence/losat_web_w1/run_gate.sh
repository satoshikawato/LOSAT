#!/usr/bin/env bash
# Gate run of LOSAT Web W1 (Session S09). It records the environment and the engine build,
# then runs
#   V-ABI quick (Node) on the reactors, npm ci, npm run check, the unit tests with every
#   case listed (with the reactors: the RecordScanner contract on the serial reactor in
#   Node), npm run e2e without the engine (FakeEngine build) and with it (V-BR, R1, cancel,
#   renewal, the serial fallback, the S10 contracts; Chromium, Firefox and WebKit), the
#   engine E2E twice more (stability), and the shared-memory growth probe,
# and, with LOSAT_WEB_MEASURE (for example `all`), the measurements of
# tests/e2e/measure.spec.ts in Chromium, and with LOSAT_WEB_MEASURE_OTHERS (for example
# `dw8,cancel`) those measurements in Firefox and WebKit.
# Logs and records go to docs/evidence/losat_web_w1/run-<UTC>/, which is not changed later.
# Before each step it waits while the engine track measures V-PERF (plan DW-7).
#
# Needs:
#   LOSAT_WEB_REACTORS  the output directory of web/adapter/tools/build_reactors.py for this commit
#   LOSAT_WEB_NATIVE    the native LOSAT of this commit (cargo +1.92.0 build --release)
# Optional:
#   LOSAT_WEB_E2E_PORT           the port of the application under test (4173)
#   LOSAT_WEB_WEBKIT_EXECUTABLE  a launcher for Playwright's WebKit on a host whose system
#                                libraries Playwright cannot install (README.md)
# Usage: docs/evidence/losat_web_w1/run_gate.sh
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
app="$repo/web/app"
lock="${LOSAT_VPERF_LOCK:-/home/kawato/.cache/losat-web-gui-target/vperf.lock}"
: "${LOSAT_WEB_REACTORS:?set LOSAT_WEB_REACTORS to the reactors of this commit}"
: "${LOSAT_WEB_NATIVE:?set LOSAT_WEB_NATIVE to the native LOSAT of this commit}"
reactors="$(cd "$LOSAT_WEB_REACTORS" && pwd)"
native="$(cd "$(dirname "$LOSAT_WEB_NATIVE")" && pwd)/$(basename "$LOSAT_WEB_NATIVE")"
# Each step sets the engine variables that it needs; the others build with the FakeEngine.
measure="${LOSAT_WEB_MEASURE:-}"
measure_others="${LOSAT_WEB_MEASURE_OTHERS:-}"
unset LOSAT_WEB_REACTORS LOSAT_WEB_NATIVE LOSAT_WEB_EVIDENCE LOSAT_WEB_MEASURE LOSAT_WEB_MEASURE_OTHERS

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
mkdir "$run"
{
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
} > "$run/environment.txt"
cp "$reactors/artifacts.json" "$run/reactors-artifacts.json"

wait_for_vperf
python3 "$repo/web/adapter/tools/v_abi_cases.py" --suite quick --out "$run/v-abi-cases.json"
node "$repo/web/adapter/tests/v_abi.js" --native "$native" \
  --serial "$reactors/losat-web-serial.wasm" --threads "$reactors/losat-web-threads.wasm" \
  --cases "$run/v-abi-cases.json" --out "$run" > "$run/v-abi.log" 2>&1

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
wait_for_vperf
LOSAT_WEB_REACTORS="$reactors" LOSAT_WEB_NATIVE="$native" \
  npx playwright test tests/e2e/engine.spec.ts --repeat-each 2 --reporter=list > "$run/e2e-engine-repeat.log" 2>&1
wait_for_vperf
node "$here/shared_memory_growth_probe.mjs" > "$run/shared-memory-growth-probe.log" 2>&1

measure_in() {
  local project="$1" set="$2"
  mkdir -p "$run/measure"
  wait_for_vperf
  LOSAT_WEB_MEASURE="$set" LOSAT_WEB_REACTORS="$reactors" LOSAT_WEB_EVIDENCE="$run/measure" \
    npx playwright test tests/e2e/measure.spec.ts --project="$project" --workers=1 --reporter=list \
    > "$run/measure/$project.log" 2>&1
}
if [ -n "$measure" ]; then measure_in chromium "$measure"; fi
if [ -n "$measure_others" ]; then
  measure_in firefox "$measure_others"
  measure_in webkit "$measure_others"
fi

echo "gate run passed: $run"

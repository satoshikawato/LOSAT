#!/usr/bin/env bash
# Gate run of LOSAT Web W7 (Session S16, delivery: the site layout with each release under
# /r/<release-id>/ and the boot page that starts the pinned release, the release manifest; the
# service worker, "Offline ready" and the explicit update, checked by V-OFF; usage statistics after
# consent through the site's /api/visit; V-PRIV with TD-11, no research data on the network and no
# machine path in a served file; onboarding: the examples, "Try an example", "Getting started" and
# the browser self-test; the Cloudflare Pages preview check; docs/losat_web_gui_sessions/session_s16_w7_delivery.md).
# Copied from docs/evidence/losat_web_w6/run_gate.sh with the same steps and these changes: Node
# >= 22.18 (build/site-server.ts and build/assemble-site.ts run under Node's TypeScript support);
# before the first step, the release IDs of this tree are made by build/release.ts's rule, and a
# change that is not committed in any input of the build (build/inputs.ts: web/app and the files
# of the repository that the build reads) ends the gate, since its releases would be `-dirty` and
# no record would name them; the expected IDs go into environment.txt; after the engine-build E2E
# and after each measurement, the served site's index and every release manifest are kept in the
# run, and a site of another release ends the gate; the delivery measurement follows the results
# measurement; the screen records include the onboarding states (moved into screens.spec.ts); and,
# for `all`, a dry run of the preview check against the engine-build site.
# It records the environment and the engine build, then runs
#   V-ABI quick (Node) on the reactors of this commit, npm ci, npm run check, the unit
#   tests with every case listed (with the reactors: the HSP records against the rows and
#   sections of every outfmt 0 fixture of BLASTN, BLASTP, TBLASTN and TBLASTX,
#   tests/unit/hsp-correspondence.test.ts, unchanged from W4; the residues extracted for the HSP
#   records of real searches, tests/unit/extraction-engine.test.ts; the Data worker's reads for
#   extraction and the candidate tray, tests/unit/data-extraction.test.ts and
#   tests/unit/candidates.test.ts; W6's files for the Writer contract, CSV, JSON and the report,
#   the dot plot SVG, settings files, reproduction and session files, as in W6's script; W7's
#   files: the release ID, the manifest and the site index, tests/unit/release-id.test.ts,
#   tests/unit/release-manifest.test.ts and tests/unit/site-index.test.ts, the build of the site
#   and its worker, tests/unit/site-build.test.ts and tests/unit/service-worker-build.test.ts, the
#   assembly with older releases, tests/unit/assemble-site.test.ts, the local site server and its
#   headers, tests/unit/site-server.test.ts and tests/unit/headers.test.ts, the boot page,
#   tests/unit/boot.test.ts, the service worker's preparation, answers and cleanup,
#   tests/unit/sw-prepare.test.ts and tests/unit/sw-serve.test.ts, the release manager, the pin and
#   the worker's client, tests/unit/release.test.ts, tests/unit/release-pin.test.ts and
#   tests/unit/offline-client.test.ts, usage statistics and the /api/visit Function,
#   tests/unit/analytics.test.ts and tests/unit/visit-function.test.ts, onboarding and the
#   self-test, tests/unit/onboarding.test.ts, tests/unit/self-test.test.ts and, with the reactors,
#   tests/unit/self-test-engine.test.ts (the self-test's expected outputs are the engine's), the
#   same bytes from two builds of the commit, the second from another checkout path,
#   tests/unit/release-bytes.test.ts, and V-PRIV's scanners (the canaries, the request rules and
#   TD-11's machine paths), tests/unit/privacy-scan.test.ts; and the generated verification table,
#   kept in the run record), the NCBI comparison commands of the
#   reproduction panel run with the native CLI and NCBI BLAST+ 2.17.0 and compared byte for byte
#   (docs/evidence/losat_web_w6/check_commands.py, unchanged, under the machine-wide oracle lock),
#   npm run e2e without the engine (FakeEngine build) and with it (every spec, W7's among them:
#   release.spec.ts, the site layout and the boot page; offline.spec.ts, V-OFF; analytics.spec.ts;
#   onboarding.spec.ts, every example searched and the self-test in each browser; privacy.spec.ts,
#   V-PRIV with TD-11; Chromium, Firefox and WebKit), the site that the engine-build E2E served
#   (dist/releases.json and dist/r/*/release-manifest.json, kept under site/ in the run), for `all`
#   a dry run of the preview check (docs/evidence/losat_web_w7/check_preview/, README there) against
#   that site served by build/site-server.ts, the results screen's E2E twice more with the engine
#   (stability), the measurements of the results screen (tests/e2e/results-measure.spec.ts: many
#   queries, and the dot plot and Graphic Summary of a pair with many HSPs, each followed by the
#   candidate tray and W6's files; since W7 each page waits for "Offline ready" first) and of
#   delivery (tests/e2e/delivery-measure.spec.ts: the first, second and offline visits, the
#   preparation until "Offline ready", the self-test, and "Try an example" to the first results of
#   each program) in Chromium, Firefox and WebKit, with the measured site's index and manifests
#   (under measure/site/ in the run), and the screen records for the visual review
#   (tests/e2e/screens.spec.ts: W6's states 01 to 40, which wait for "Offline ready", and W7's 41
#   to 56) in the three browsers.
# Logs and records go to docs/evidence/losat_web_w7/run-<UTC>/, which is not changed later;
# the screen records (PNG) go to $LOSAT_WEB_SCREENS_OUT (outside the repository; the
# README lists their SHA-256).
# Before each step it waits while the engine track measures V-PERF (plan DW-7).
#
# Needs:
#   LOSAT_WEB_REACTORS  the output directory of web/adapter/tools/build_reactors.py for this commit
#   LOSAT_WEB_NATIVE    the native LOSAT of this commit (cargo +1.92.0 build --release)
#   LOSAT_WEB_SCREENS_OUT  where the screen records go
#   LOSAT_WEB_NCBI_BIN  the bin directory of NCBI BLAST+ 2.17.0 (for `all`: the comparison commands)
#   Node >= 22.18
# Optional:
#   LOSAT_WEB_E2E_PORT           the port of the application under test (4173); V-OFF, the screen
#                                records of W7 and the delivery measurement also use this port
#                                + 100, + 200, ... (one per Playwright worker)
#   LOSAT_WEB_PREVIEW_PORT       the port of the preview check's dry run (the E2E port + 50)
#   LOSAT_WEB_WEBKIT_EXECUTABLE  a launcher for Playwright's WebKit on a host whose system
#                                libraries Playwright cannot install (W1 README)
#   LOSAT_WEB_MEASURE_QUERIES    query counts of the results measurement (10000,100000)
#   LOSAT_WEB_MEASURE_COPIES     copies of the repeated unit of its dot-plot case (1500,3000,4000)
#   LOSAT_ORACLE_LOCK            the machine-wide lock of the oracle runs
#                                (/home/kawato/.cache/losat-work/oracle.lock)
#   LOSAT_WEB_COMMANDS_WORK      where the comparison commands run, outside the repository
#                                ($HOME/.cache/losat-work/s16-commands)
#   LOSAT_WEB_GATE_STEPS         `all` (the gate), `after-review`: npm ci, check, the unit
#                                tests, both E2E builds, the site's records and the screen
#                                records only (the run after a fix of a review finding; README),
#                                or `measure`: npm ci and the measurements only (after a fix of
#                                a measurement spec)
# Usage: bash docs/evidence/losat_web_w7/run_gate.sh
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
app="$repo/web/app"
w6="$repo/docs/evidence/losat_web_w6"
lock="${LOSAT_VPERF_LOCK:-/home/kawato/.cache/losat-web-gui-target/vperf.lock}"
: "${LOSAT_WEB_REACTORS:?set LOSAT_WEB_REACTORS to the reactors of this commit}"
: "${LOSAT_WEB_NATIVE:?set LOSAT_WEB_NATIVE to the native LOSAT of this commit}"
: "${LOSAT_WEB_SCREENS_OUT:?set LOSAT_WEB_SCREENS_OUT to a directory for the screen records}"
# build/site-server.ts (the E2E server) and build/assemble-site.ts run under Node's TypeScript support (package.json engines).
if ! node -e 'const [major, minor] = process.versions.node.split(".").map(Number); process.exit(major > 22 || (major === 22 && minor >= 18) ? 0 : 1)'; then
  echo "Node >= 22.18 is needed, not $(node --version)" >&2
  exit 2
fi
reactors="$(cd "$LOSAT_WEB_REACTORS" && pwd)"
native="$(cd "$(dirname "$LOSAT_WEB_NATIVE")" && pwd)/$(basename "$LOSAT_WEB_NATIVE")"
screens="$LOSAT_WEB_SCREENS_OUT"
queries="${LOSAT_WEB_MEASURE_QUERIES:-10000,100000}"
steps="${LOSAT_WEB_GATE_STEPS:-all}"
case "$steps" in all | after-review | measure) ;; *) echo "LOSAT_WEB_GATE_STEPS is all, after-review or measure" >&2; exit 2 ;; esac
oracle_lock="${LOSAT_ORACLE_LOCK:-/home/kawato/.cache/losat-work/oracle.lock}"
commands_work="${LOSAT_WEB_COMMANDS_WORK:-$HOME/.cache/losat-work/s16-commands}"
ncbi_bin=""
if [ "$steps" = all ]; then
  : "${LOSAT_WEB_NCBI_BIN:?set LOSAT_WEB_NCBI_BIN to the bin directory of NCBI BLAST+ 2.17.0}"
  ncbi_bin="$(cd "$LOSAT_WEB_NCBI_BIN" && pwd)"
fi
# Each step sets the engine variables that it needs; the others build with the FakeEngine.
unset LOSAT_WEB_REACTORS LOSAT_WEB_NATIVE LOSAT_WEB_EVIDENCE LOSAT_WEB_MEASURE LOSAT_WEB_SCREENS LOSAT_WEB_SCREENS_OUT
# The site is built with the release ID of this tree (build/release.ts), never one set by hand,
# and without the measurement endpoint (analytics.spec.ts serves its own site index).
unset LOSAT_WEB_RELEASE_ID LOSAT_WEB_ANALYTICS_ENDPOINT
# `vite build` makes a production build only while NODE_ENV is unset (or production): another value
# gives a development build, with each Vue component's absolute path in its bytes. Vite also takes
# VITE_USER_NODE_ENV (exported, or NODE_ENV in a .env file) as NODE_ENV; vite.config.ts stops such
# a build, and the gate does not inherit one.
unset NODE_ENV VITE_USER_NODE_ENV

# The release IDs of this tree, before any step: build/release.ts's rule over the inputs that
# build/inputs.ts lists. A change that is not committed in any of them would make every build of
# the gate a `-dirty` release that no record names, so it ends the gate here.
release_ids="$(cd "$app" && node --input-type=module -e '
import { BUILD_INPUTS } from "./build/inputs.ts";
import { findReactors } from "./build/reactors.ts";
import { releaseId } from "./build/release.ts";
console.log(`fake ${releaseId(undefined, { env: {} })}`);
console.log(`engine ${releaseId(findReactors({ LOSAT_WEB_REACTORS: process.argv[1] }), { env: {} })}`);
console.log(`inputs ${BUILD_INPUTS.join(" ")}`);
' "$reactors")"
fake_release="$(sed -n 's/^fake //p' <<< "$release_ids")"
engine_release="$(sed -n 's/^engine //p' <<< "$release_ids")"
read -r -a build_inputs <<< "$(sed -n 's/^inputs //p' <<< "$release_ids")"
uncommitted="$(git -C "$repo" status --porcelain -- "${build_inputs[@]}")"
case "$fake_release $engine_release" in
  *-dirty*)
    echo "the release IDs of this tree are $fake_release and $engine_release: commit the changes in the build's inputs (build/inputs.ts) and run the gate again" >&2
    [ -z "$uncommitted" ] || echo "$uncommitted" >&2
    exit 1
    ;;
esac

wait_for_vperf() {
  until [ ! -e "$lock" ]; do
    echo "waiting: V-PERF is being measured ($lock)" >&2
    sleep 60
  done
}

port="${LOSAT_WEB_E2E_PORT:-4173}"
preview_port="${LOSAT_WEB_PREVIEW_PORT:-$((port + 50))}"
# The front servers of V-OFF, the W7 screen records and the delivery measurement use the next hundreds.
for used in "$port" "$((port + 100))" "$((port + 200))" "$preview_port"; do
  if curl -s -o /dev/null "http://localhost:$used/"; then
    echo "port $used is in use (set LOSAT_WEB_E2E_PORT or LOSAT_WEB_PREVIEW_PORT to another port)" >&2
    exit 1
  fi
done

run="$here/run-$(date -u +%Y%m%dT%H%M%SZ)"
[ "$steps" = all ] || run="$run-$steps"
mkdir "$run"
{
  echo "steps: $steps"
  echo "commit: $(git -C "$repo" rev-parse HEAD)"
  echo "branch: $(git -C "$repo" branch --show-current)"
  echo "uncommitted files in the build's inputs (build/inputs.ts: ${build_inputs[*]}): $(grep -c . <<< "$uncommitted" || true)"
  echo "release (FakeEngine build): $fake_release"
  echo "release (engine build): $engine_release"
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
  if [ -n "$ncbi_bin" ]; then
    echo "ncbi blast+: $ncbi_bin ($("$ncbi_bin/blastn" -version 2>/dev/null | head -1 || echo unknown))"
  fi
} > "$run/environment.txt"
cp "$reactors/artifacts.json" "$run/reactors-artifacts.json"

# Keeps the index of the site in dist/ (built by Playwright's web server for the step that just ran)
# and the manifest of each of its releases, as the site serves them, under $1; a site whose release
# is not this tree's engine build (checked before the first step) ends the gate.
keep_site() {
  local manifest release latest
  mkdir -p "$1/r"
  cp dist/releases.json "$1/releases.json"
  for manifest in dist/r/*/release-manifest.json; do
    release="$(basename "$(dirname "$manifest")")"
    mkdir -p "$1/r/$release"
    cp "$manifest" "$1/r/$release/release-manifest.json"
  done
  latest="$(node -e 'console.log(JSON.parse(require("node:fs").readFileSync(process.argv[1], "utf8")).latest)' "$1/releases.json")"
  if [ "$latest" != "$engine_release" ]; then
    echo "the site in dist/ is release $latest, not $engine_release (this tree's engine build): did a build input change during the gate?" >&2
    exit 1
  fi
}

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
if [ "$steps" != measure ]; then
  wait_for_vperf
  npm run check > "$run/npm-check.log" 2>&1
  wait_for_vperf
  LOSAT_WEB_REACTORS="$reactors" LOSAT_WEB_VERIFICATION_OUT="$run/verification-table.json" \
    npx vitest run --reporter=verbose > "$run/unit-cases.log" 2>&1
fi
if [ "$steps" = all ]; then
  # The NCBI comparison commands of the reproduction panel (S15 item 5, design §12.3): the commands
  # that tests/unit/reproduce.test.ts writes for its fixed argvs, run with the native CLI and with
  # NCBI BLAST+ and compared byte for byte; an unexpected difference ends the gate. W6's script, unchanged.
  wait_for_vperf
  work="$commands_work/$(basename "$run")"
  flock "$oracle_lock" python3 "$w6/check_commands.py" --losat "$native" --ncbi-bin "$ncbi_bin" \
    --work "$work" --out "$run/check-commands.tsv" --jobs 2 > "$run/check-commands.log" 2>&1
  cp "$work/commands.json" "$run/check-commands.json"
fi
if [ "$steps" != measure ]; then
  wait_for_vperf
  npm run e2e -- --reporter=list > "$run/npm-e2e-fake-engine.log" 2>&1
  wait_for_vperf
  mkdir "$run/records"
  LOSAT_WEB_REACTORS="$reactors" LOSAT_WEB_NATIVE="$native" LOSAT_WEB_EVIDENCE="$run/records" \
    npm run e2e -- --reporter=list > "$run/npm-e2e-engine.log" 2>&1

  # The site that the engine-build E2E served: its index and the manifest of each release.
  keep_site "$run/site"
fi
if [ "$steps" = all ]; then
  # A dry run of the preview check (S16 item 7) against the same site served as Cloudflare Pages
  # serves it (build/site-server.ts with dist/_headers); the check's own records go to the run.
  wait_for_vperf
  node build/site-server.ts dist --port "$preview_port" > "$run/preview-check-server.log" 2>&1 &
  server=$!
  trap 'kill "$server" 2>/dev/null || true' EXIT
  for _ in $(seq 1 30); do
    curl -s -o /dev/null "http://localhost:$preview_port/" && break
    sleep 1
  done
  LOSAT_WEB_PREVIEW_URL="http://localhost:$preview_port" LOSAT_WEB_PREVIEW_HEADERS=dist/_headers \
    LOSAT_WEB_PREVIEW_OUT="$run/preview-check" \
    npx playwright test -c ../../docs/evidence/losat_web_w7/check_preview/playwright.config.ts > "$run/preview-check.log" 2>&1
  kill "$server"
  wait "$server" 2>/dev/null || true
  trap - EXIT
fi
if [ "$steps" = all ]; then
  wait_for_vperf
  LOSAT_WEB_REACTORS="$reactors" \
    npx playwright test tests/e2e/results.spec.ts --repeat-each 2 --reporter=list > "$run/e2e-results-repeat.log" 2>&1
fi
if [ "$steps" != after-review ]; then
  mkdir -p "$run/measure"
  for project in chromium firefox webkit; do
    wait_for_vperf
    LOSAT_WEB_MEASURE=results LOSAT_WEB_MEASURE_QUERIES="$queries" LOSAT_WEB_REACTORS="$reactors" \
      LOSAT_WEB_EVIDENCE="$run/measure" \
      npx playwright test tests/e2e/results-measure.spec.ts --project="$project" --workers=1 --reporter=list \
      > "$run/measure/$project.log" 2>&1
  done
  # The measured site (both measurements serve the engine build in dist/): its index and manifests.
  keep_site "$run/measure/site"
  # Delivery (S16): visits, the preparation until "Offline ready", the self-test and the examples.
  for project in chromium firefox webkit; do
    wait_for_vperf
    LOSAT_WEB_MEASURE=delivery LOSAT_WEB_REACTORS="$reactors" LOSAT_WEB_EVIDENCE="$run/measure" \
      npx playwright test tests/e2e/delivery-measure.spec.ts --project="$project" --workers=1 --reporter=list \
      > "$run/measure/delivery-$project.log" 2>&1
  done
  keep_site "$run/measure/site"
fi

if [ "$steps" != measure ]; then
  # W6's states and W7's, the onboarding states 41 to 45 among them (screens.spec.ts).
  wait_for_vperf
  mkdir -p "$screens"
  LOSAT_WEB_SCREENS="$screens" LOSAT_WEB_REACTORS="$reactors" \
    npx playwright test tests/e2e/screens.spec.ts --reporter=list > "$run/screens.log" 2>&1
  (cd "$screens" && find . -name '*.png' | sort | xargs sha256sum) > "$run/screens.sha256"
fi

echo "gate run passed: $run"

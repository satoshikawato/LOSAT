#!/usr/bin/env bash
# Gate run of LOSAT Web W6 (Session S15, outputs and reproducibility: the Writer contract, the
# compatibility outputs and the extraction written in blocks, CSV, JSON and the static report,
# settings files, Edit Search and the reproduction panel, session files and the re-attachment of
# the original FASTA, the dot plot as SVG; docs/losat_web_gui_sessions/session_s15_w6_export_session.md).
# Copied from docs/evidence/losat_web_w5/run_gate.sh with the same steps and one more, the check of
# the NCBI comparison commands. It records the environment and the engine build, then runs
#   V-ABI quick (Node) on the reactors of this commit, npm ci, npm run check, the unit
#   tests with every case listed (with the reactors: the HSP records against the rows and
#   sections of every outfmt 0 fixture of BLASTN, BLASTP, TBLASTN and TBLASTX,
#   tests/unit/hsp-correspondence.test.ts, unchanged from W4; the residues extracted for the HSP
#   records of real searches, tests/unit/extraction-engine.test.ts; the Data worker's reads for
#   extraction and the candidate tray, tests/unit/data-extraction.test.ts and
#   tests/unit/candidates.test.ts; W6's files: the Writer contract, tests/unit/export-writer.test.ts,
#   the extraction's bounded steps, tests/unit/extraction-steps.test.ts, CSV, JSON and the report,
#   tests/unit/hsp-export.test.ts and tests/unit/result-export.test.ts, the dot plot SVG,
#   tests/unit/plot-svg.test.ts, settings files, the commands and the run's input FASTA,
#   tests/unit/settings-file.test.ts, tests/unit/reproduce.test.ts and tests/unit/run-files.test.ts,
#   and session files, tests/unit/session-file.test.ts, tests/unit/session-data.test.ts and
#   tests/unit/session-save-load.test.ts; and the generated verification table, kept in the run
#   record), the NCBI comparison commands of the reproduction panel run with the native CLI and
#   NCBI BLAST+ 2.17.0 and compared byte for byte (docs/evidence/losat_web_w6/check_commands.py,
#   under the machine-wide oracle lock), npm run e2e without the engine (FakeEngine build) and with
#   it (every spec, the search screen's search.spec.ts, the results screen's results.spec.ts, the
#   candidate tray's candidates.spec.ts, and W6's exports.spec.ts, repro.spec.ts and session.spec.ts
#   - real searches saved in a session file and opened without searching, their outputs compared
#   with the native CLI's - among them; Chromium, Firefox and WebKit), the results screen's E2E
#   twice more with the engine (stability), the measurements of the results screen
#   (tests/e2e/results-measure.spec.ts: many queries, and the dot plot and Graphic Summary of a
#   pair with many HSPs, each followed by the candidate tray filled from them, then by W6's files:
#   CSV, JSON, the report and the outfmt 6 export of the whole run, the dot plot as SVG, and the
#   session file saved and opened in a new page) in Chromium, Firefox and WebKit, and the screen
#   records for the visual review (tests/e2e/screens.spec.ts; W6's states 29 to 40 among them) in
#   the three browsers.
# Logs and records go to docs/evidence/losat_web_w6/run-<UTC>/, which is not changed later;
# the screen records (PNG) go to $LOSAT_WEB_SCREENS_OUT (outside the repository; the
# README lists their SHA-256).
# Before each step it waits while the engine track measures V-PERF (plan DW-7).
#
# Needs:
#   LOSAT_WEB_REACTORS  the output directory of web/adapter/tools/build_reactors.py for this commit
#   LOSAT_WEB_NATIVE    the native LOSAT of this commit (cargo +1.92.0 build --release)
#   LOSAT_WEB_SCREENS_OUT  where the screen records go
#   LOSAT_WEB_NCBI_BIN  the bin directory of NCBI BLAST+ 2.17.0 (for `all`: the comparison commands)
# Optional:
#   LOSAT_WEB_E2E_PORT           the port of the application under test (4173)
#   LOSAT_WEB_WEBKIT_EXECUTABLE  a launcher for Playwright's WebKit on a host whose system
#                                libraries Playwright cannot install (W1 README)
#   LOSAT_WEB_MEASURE_QUERIES    query counts of the results measurement (10000,100000)
#   LOSAT_WEB_MEASURE_COPIES     copies of the repeated unit of its dot-plot case (1500,3000,4000)
#   LOSAT_ORACLE_LOCK            the machine-wide lock of the oracle runs
#                                (/home/kawato/.cache/losat-work/oracle.lock)
#   LOSAT_WEB_COMMANDS_WORK      where the comparison commands run, outside the repository
#                                ($HOME/.cache/losat-work/s15-commands)
#   LOSAT_WEB_GATE_STEPS         `all` (the gate), `after-review`: npm ci, check, the unit
#                                tests, both E2E builds and the screen records only (the
#                                run after a fix of a review finding; README), or `measure`:
#                                npm ci and the measurements only (after a fix of the
#                                measurement spec)
# Usage: bash docs/evidence/losat_web_w6/run_gate.sh
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
queries="${LOSAT_WEB_MEASURE_QUERIES:-10000,100000}"
steps="${LOSAT_WEB_GATE_STEPS:-all}"
case "$steps" in all | after-review | measure) ;; *) echo "LOSAT_WEB_GATE_STEPS is all, after-review or measure" >&2; exit 2 ;; esac
oracle_lock="${LOSAT_ORACLE_LOCK:-/home/kawato/.cache/losat-work/oracle.lock}"
commands_work="${LOSAT_WEB_COMMANDS_WORK:-$HOME/.cache/losat-work/s15-commands}"
ncbi_bin=""
if [ "$steps" = all ]; then
  : "${LOSAT_WEB_NCBI_BIN:?set LOSAT_WEB_NCBI_BIN to the bin directory of NCBI BLAST+ 2.17.0}"
  ncbi_bin="$(cd "$LOSAT_WEB_NCBI_BIN" && pwd)"
fi
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
  if [ -n "$ncbi_bin" ]; then
    echo "ncbi blast+: $ncbi_bin ($("$ncbi_bin/blastn" -version 2>/dev/null | head -1 || echo unknown))"
  fi
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
  # NCBI BLAST+ and compared byte for byte; an unexpected difference ends the gate.
  wait_for_vperf
  work="$commands_work/$(basename "$run")"
  flock "$oracle_lock" python3 "$here/check_commands.py" --losat "$native" --ncbi-bin "$ncbi_bin" \
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
fi

if [ "$steps" != measure ]; then
  wait_for_vperf
  mkdir -p "$screens"
  LOSAT_WEB_SCREENS="$screens" LOSAT_WEB_REACTORS="$reactors" \
    npx playwright test tests/e2e/screens.spec.ts --reporter=list > "$run/screens.log" 2>&1
  (cd "$screens" && find . -name '*.png' | sort | xargs sha256sum) > "$run/screens.sha256"
fi

echo "gate run passed: $run"

#!/usr/bin/env bash
# SFb (E2h): continue sf_gates.sh after a stop (a failed build, a WSL restart, a usage limit) with the
# same run directory and build artifacts. Stages whose marker file ($OUT/.done/<stage>) exists are
# skipped; a stage that was running is run again from its start (its status rows are replaced).
# Before anything else it checks that HEAD and the build artifacts are those of the run
# (artifacts.sha256 of the build stage), as s11_gates_resume.sh did. The lexical fixture in /tmp is
# staged again by every stage (a WSL restart wipes it).
#   sf_gates_resume.sh                  all stages that are not done
#   STAGES="vabi-tblastx collect" sf_gates_resume.sh
#   FORCE=1 STAGES=capture sf_gates_resume.sh     run a finished stage again
set -u
WORK_ROOT=${WORK_ROOT:-/home/kawato/losat-work}
BUILD_ROOT=${BUILD_ROOT:-/home/kawato/.cache/losat-work}
W=${WT:-$WORK_ROOT/.worktrees/web-gui}
P=${GATE:-gate-sf}
POINTER=$BUILD_ROOT/sfb-e2h/$P.run
[ -f "$POINTER" ] || { echo "no run to resume: $POINTER is missing (start with sf_gates.sh)" >&2; exit 2; }
. "$POINTER"
OUT=$BUILD_ROOT/sfb-e2h/gate-$SF_TS
cd "$W" || exit 2
if [ -f "$OUT/head.txt" ]; then
  git rev-parse HEAD > "$OUT/head-resume.txt"
  cmp -s "$OUT/head.txt" "$OUT/head-resume.txt" || { echo "FAILED HEAD is not the run's ($(cat "$OUT/head.txt"))"; exit 1; }
fi
if [ -f "$OUT/artifacts.sha256" ]; then
  N=$BUILD_ROOT/$P/native/release/LOSAT; NS=$BUILD_ROOT/$P/native-serial/release/LOSAT
  WA=$BUILD_ROOT/$P/wasi-artifacts; R=$BUILD_ROOT/$P/reactors
  sha256sum "$N" "$NS" "$WA"/*.wasm "$R"/*.wasm > "$OUT/artifacts-resume.sha256"
  cmp -s "$OUT/artifacts.sha256" "$OUT/artifacts-resume.sha256" || { echo "FAILED the build artifacts changed since the build stage"; exit 1; }
fi
echo "$(date -u +%H:%M:%S) resume run $SF_TS; done: $(ls "$OUT/.done" 2>/dev/null | tr '\n' ' ')"
export RESUME=1 SF_TS
exec "$(dirname "$(readlink -f "$0")")/sf_gates.sh"

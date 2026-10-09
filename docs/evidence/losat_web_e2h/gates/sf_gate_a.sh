#!/usr/bin/env bash
# SFb (E2h) Gate A: TBLASTX v0.1.0 outfmt 6 parity (20 pairs) of the gate's native build. From s11_gate_a.sh.
# Runs ALONE (3-5.6 h): start it after sf_gates.sh has finished, under the machine-wide lock, and start
# nothing else that runs searches meanwhile. The fixture /tmp/losat-pr5-runtime-cert-* (a WSL restart wipes it)
# is staged again here from stage_lexical_fixtures() in LOSAT/tests/ci_fast_regressions.py.
# Logs go to the run's out directory; the summaries are copied to docs/evidence/losat_web_e2h/run-<UTC>/.
set -u
WORK_ROOT=${WORK_ROOT:-/home/kawato/losat-work}
BUILD_ROOT=${BUILD_ROOT:-/home/kawato/.cache/losat-work}
W=${WT:-$WORK_ROOT/.worktrees/web-gui}
P=${GATE:-gate-sf}
LOCK=$BUILD_ROOT/oracle.lock
if [ -z "${SF_LOCKED:-}" ]; then SF_LOCKED=1 exec flock "$LOCK" "$0" "$@"; fi
. "$BUILD_ROOT/sfb-e2h/$P.run"
OUT=$BUILD_ROOT/sfb-e2h/gate-$SF_TS
RUN=$W/docs/evidence/losat_web_e2h/run-$SF_TS
GA=$OUT/gate-a-tblastx-v010
N=$BUILD_ROOT/$P/native/release/LOSAT
export RUSTUP_TOOLCHAIN=1.92.0 PYTHONDONTWRITEBYTECODE=1
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd "$W" || fail worktree
[ -x "$N" ] || fail "no gate build at $N (run sf_gates.sh first)"
if [ -f "$OUT/head.txt" ]; then [ "$(git rev-parse HEAD)" = "$(cat "$OUT/head.txt")" ] || fail "HEAD is not the run's"; fi
if [ -f "$OUT/artifacts.sha256" ]; then
  [ "$(sha256sum "$N" | cut -d' ' -f1)" = "$(grep " $N\$" "$OUT/artifacts.sha256" | cut -d' ' -f1)" ] || fail "native build changed since the build stage"
fi
(cd LOSAT/tests && python3 -c 'import ci_fast_regressions as c; c.stage_lexical_fixtures()') || fail lexical-fixtures
ls -d /tmp/losat-pr5-runtime-cert-* > "$OUT/gate-a-fixture-root.txt" || fail "no /tmp/losat-pr5-runtime-cert-* after staging"
read -r load _ < /proc/loadavg; step "load average $load (Gate A runs alone)"
step "gate a: tblastx v0.1.0 (outfmt 6), 20 pairs"
rm -rf "$GA"
python3 LOSAT/tests/audit_tblastx_v010.py --losat-bin "$N" --output-dir "$GA" > "$OUT/audit-tblastx-v010.log" 2>&1; rc=$?
echo "exit $rc" >> "$OUT/audit-tblastx-v010.log"
mkdir -p "$OUT" "$RUN/audit-tblastx-v010"
cp "$GA"/*.json "$GA"/*.tsv "$RUN/audit-tblastx-v010/" 2>/dev/null
cp "$OUT/audit-tblastx-v010.log" "$RUN/" 2>/dev/null
printf 'gate-a\tgate-a-tblastx-v010\t%s\taudit-tblastx-v010.log\n' "$rc" >> "$OUT/status.tsv"
touch "$OUT/.done/gate-a"
step "done rc $rc"
exit $rc

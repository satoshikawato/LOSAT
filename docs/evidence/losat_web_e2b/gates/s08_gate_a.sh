#!/usr/bin/env bash
# S08 (E2b, TBLASTX outfmt 0/7) gate runs: lint, tests, builds, then the gates; logs go
# into the S08 run directory (docs/evidence/losat_web_e2b/run-*).
set -u
A=/home/kawato/.cache/losat-web-gui-target
S=$A/s08
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $S/rundir)
NCBI=/home/kawato/micromamba/bin
export RUSTUP_TOOLCHAIN=1.92.0
step() { echo "$(date -u +%H:%M:%S) $*"; }
fail() { echo "$(date -u +%H:%M:%S) FAILED $*"; exit 1; }
cd $W
N=$A/s08-gate-native/release/LOSAT
step "gate a: tblastx v0.1.0 (outfmt 6)"
rm -rf /tmp/claude-1000/s08-tblastx-v010
python3 LOSAT/tests/audit_tblastx_v010.py --losat-bin $N --output-dir /tmp/claude-1000/s08-tblastx-v010 > $RUN/audit-tblastx-v010.log 2>&1; echo "exit $?" >> $RUN/audit-tblastx-v010.log
mkdir -p $RUN/audit-tblastx-v010 && cp /tmp/claude-1000/s08-tblastx-v010/*.json /tmp/claude-1000/s08-tblastx-v010/*.tsv $RUN/audit-tblastx-v010/ 2>/dev/null

step "done"

#!/usr/bin/env bash
# S11 (E2d) Gate A: TBLASTX v0.1.0 outfmt 6 parity of the gate's native build; logs go
# into the run directory (docs/evidence/losat_web_e2d/run-*). From s08pb_gate_a.sh (E2e).
set -u
A=/home/kawato/.cache/losat-web-gui-target
B=$A/s11
W=/mnt/c/Users/genom/GitHub/LOSAT-web-gui
RUN=$W/$(cat $B/rundir)
OUT=$A/${GATE:-s11-gate}-tblastx-v010
step() { echo "$(date -u +%H:%M:%S) $*"; }
cd $W
N=$A/${GATE:-s11-gate}-native/release/LOSAT
step "gate a: tblastx v0.1.0 (outfmt 6)"
rm -rf $OUT
python3 LOSAT/tests/audit_tblastx_v010.py --losat-bin $N --output-dir $OUT > $RUN/audit-tblastx-v010.log 2>&1; echo "exit $?" >> $RUN/audit-tblastx-v010.log
mkdir -p $RUN/audit-tblastx-v010 && cp $OUT/*.json $OUT/*.tsv $RUN/audit-tblastx-v010/ 2>/dev/null
step "done"

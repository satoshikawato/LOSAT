#!/usr/bin/env bash
# Gate A once the gate run has rebuilt the native binary (its "build wasi" step).
S=/home/kawato/.cache/losat-web-gui-target/s08
until grep -qE "build wasi|FAILED" $S/gates.log 2>/dev/null; do sleep 20; done
grep -q FAILED $S/gates.log && { echo "gate failed before the native build"; exit 1; }
exec $S/s08_gate_a.sh

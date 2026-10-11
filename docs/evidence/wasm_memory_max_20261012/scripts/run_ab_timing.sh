#!/usr/bin/env bash
# The A/B timing, one case per hold of the machine-wide lock.
set -uo pipefail
T=/home/kawato/losat-baselines/wasm-mem-20261012
cd "$T"; rm -f logs/ab_timing.tsv logs/ab_timing_summary.txt
for c in p03 blastp_ap d04; do
  flock /home/kawato/.cache/losat-work/oracle.lock python3 scripts/ab_time.py timing --rounds 5 --case $c
done

#!/usr/bin/env bash
# WG step 1: gbdraw's threaded host in Chromium (python Playwright), arms interleaved and rotated per round.
# usage: rounds.sh <rounds> <jobs file> <arm>...
set -u
H=/home/kawato/losat-baselines/wg-wasm-guard-20261101/cmdharness
R=$1; J=$2; shift 2; ARMS=("$@")
cd $H
python3 serve.py 9893 > logs/serve.log 2>&1 &
SERVE=$!
trap 'kill $SERVE 2>/dev/null' EXIT
sleep 2
n=${#ARMS[@]}
for ((r = 0; r < R; r++)); do
  for ((k = 0; k < n; k++)); do
    a=${ARMS[$(((k + r) % n))]}
    echo "== round $r arm $a $(date +%T) load $(cut -d' ' -f1 /proc/loadavg)"
    timeout 3600 python3 run.py $a 1 $J 9893 2>&1 | grep -v '^  File\|^    ' | tail -4
  done
done
echo "DONE $(date +%T)"

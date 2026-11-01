#!/usr/bin/env bash
# Line WO: trap counts of the 9 modes on 8 idle cores of node 0, untimed. Waits for the cores to be free of
# other lines' markers and load (at most 3 h), since node 0 had every core named when this started.
set -u
B=$HOME/losat-bench; W=$B/wo; CPUS=4-7,16,18,19,30,68-71,80,82,83,94
M=$B/state/NODE0_SHARED_wo-modes.txt
for ((i = 0; i < 180; i++)); do
  python3 $W/scripts/cpuguard.py "$CPUS" --timed > /dev/null && break
  sleep 60
done
python3 $W/scripts/cpuguard.py "$CPUS" --timed || { echo GUARD_REFUSED; echo MODES_DONE; exit 3; }
echo "shared by line WO (local losat-baselines-4c, /home/kawato/losat-baselines/wasm-oob-20261101) since $(date '+%F %T'): node 0 CPUs 4-7,16,18,19,30, siblings 68-71,80,82,83,94 ONLY, UNTIMED trap counts of the 9 modes at 16 threads (~60 min; wo/scripts/wo_modes_n0.sh); driver PID $$; removed at its end" > "$M"
sleep 3
python3 $W/scripts/cpuguard.py "$CPUS" || { rm -f "$M"; echo GUARD_REFUSED; echo MODES_DONE; exit 3; }
bash $W/scripts/wo_modes.sh "$CPUS" 0 10
rm -f "$M"
echo MODES_DONE

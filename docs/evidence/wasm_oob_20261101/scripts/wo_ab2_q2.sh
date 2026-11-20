#!/usr/bin/env bash
# Line WO: one node 1 slot, ab2 at 8 threads. Replaces wo_ab2_q.sh and keeps its queue line (same tag,
# no add). The modes run on node 0 (wo_modes_n0.sh): they count traps and read no time (Owner 2026-11-02 01:3x).
set -u
B=$HOME/losat-bench; W=$B/wo; TAG=wo-ab2
bash $B/scripts/node1q.sh wait $TAG
bash $W/scripts/wo_ab2.sh ab2 32-39,96-103 32-39 96-103 1 8 9 bp bpff mb bn dc bns tbn bx tbx
bash $B/scripts/node1q.sh done $TAG
echo AB2_SLOT_DONE

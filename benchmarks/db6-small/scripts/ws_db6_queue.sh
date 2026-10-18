#!/bin/bash
# dbsearch-20261012 session 14 (DB-6 small tier, DESIGN_db6.md). Owner 2026-10-17 08:30: decide by per-CPU use, not by
# markers alone. Every minute: pick 16 physical cores of one node whose cores and HT siblings are idle and not named
# by another IN_USE or TIMED marker, not on a node whose TIMED marker asks to keep heavy work off it, and, while an older
# NODE_QUEUED_ marker waits, only on a node that already has a marker (cpu_busy.py pick .. polite: the older queues
# take the next node with no marker), and in the same step write the IN_USE marker naming exactly
# those CPUs; then run ws_db6_bench.sh (TIMED, ~80-100 min) with BIN (built before by ws_db6_build.sh).
#   ~/losat-bench/db3/task/scripts/ws_db6_queue.sh BIN
set -u
N=$1
B=$HOME/losat-bench; ST=$B/state; D3=$B/db3; R=$D3/repo; L=$D3/logs/s14-db6-queue.txt
CB="python3 $D3/task/scripts/cpu_busy.py"
Q=$ST/NODE_QUEUED_dbsearch-s14-db6.txt
echo "dbsearch-20261012 session 14 (local losat-baselines-4d) queued $(date '+%F %T'): TIMED DB-6 benchmark (~80-100 min) on 16 physical cores + siblings that are idle by per-CPU use; yields to older NODE_QUEUED_ markers (takes only a node that already has a marker while they wait); db3/task/scripts/ws_db6_queue.sh; driver PID $$" > $Q
M=""
trap 'rm -f "$Q" ${M:+"$M"}' EXIT
echo "$(date '+%F %T') waiting for 16 idle physical cores + siblings, PID $$" >> $L
while true; do
  pick=$($CB pick 3 16 $ST polite $Q)
  if [ -n "$pick" ]; then
    read -r node CPUS SIBS <<< "$pick"
    M=$ST/NODE${node}_IN_USE_dbsearch-s14-db6.txt
    echo "dbsearch-20261012 session 14 (local losat-baselines-4d) since $(date '+%F %T'): node $node CPUs $CPUS and their siblings $SIBS ONLY (TIMED, membind $node; idle by per-CPU use when taken): DB-6 benchmark NCBI -db vs LOSAT -db n1/n8/n16 (~80-100 min) (db3/task/scripts/ws_db6_queue.sh, ws_db6_bench.sh); driver PID $$; removed at its end. Every other CPU is free." > $M
    sleep 5
    clash=$($CB overlap $ST $CPUS,$SIBS $M)
    if [ -n "$clash" ]; then
      echo "$(date '+%F %T') withdrew: CPUs also named by older marker(s) $clash" >> $L
      rm -f $M
    else
      rm -f $Q
      break
    fi
  fi
  sleep 60
done
echo "$(date '+%F %T') start node $node CPUs $CPUS siblings $SIBS; load $(cut -d' ' -f1-3 /proc/loadavg)" >> $L
CPUS=$CPUS SIBS=$SIBS MEMNODE=$node bash $D3/task/scripts/ws_db6_bench.sh $N >> $L 2>&1
echo "$(date '+%F %T') chain done" >> $L

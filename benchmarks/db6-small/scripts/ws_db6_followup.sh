#!/bin/bash
# dbsearch-20261012 session 14, DB-6 follow-up (decision 151). Question: are the two thread losses of DB-6 (BP-yeast n16
# 1.60, BX-t10 n8 1.73 / n16 2.20 of NCBI -db) the known opt-in fixes' cases? Arms in the same window, 3 alternating
# rounds after a warmup: BP-yeast n16 default / LOSAT_X_BPDBCUT=1 / NCBI; BX-t10 n8 and n16 default / LOSAT_X_BXFLOW=1 /
# NCBI; then one LOSAT_TIMING=1 run of each default case (untimed). Same binary as DB-6 (origin/dev 0b7335ce).
# Picks 16 idle cores + siblings like ws_db6_queue.sh (cpu_busy.py pick .. polite, overlap re-check).
set -u
B=$HOME/losat-bench; ST=$B/state; D3=$B/db3; W=$D3/db6; O=$W/fu; BIN=$D3/bin/LOSAT-db6-0b7335ce
CB="python3 $D3/task/scripts/cpu_busy.py"; NB=$B/tools/env_blast/bin
mkdir -p $O; S=$O/summary.txt; T=$O/runs.tsv
Q=$ST/NODE_QUEUED_dbsearch-s14-db6fu.txt; M=""
echo "dbsearch-20261012 session 14 (local losat-baselines-4d) queued $(date '+%F %T'): TIMED DB-6 follow-up ~25 min on 16 idle cores + siblings (db3/task/scripts/ws_db6_followup.sh); driver PID $$" > $Q
trap 'rm -f "$Q" ${M:+"$M"}' EXIT
while true; do
  pick=$($CB pick 3 16 $ST polite $Q)
  if [ -n "$pick" ]; then
    read -r node CPUS SIBS <<< "$pick"
    M=$ST/NODE${node}_IN_USE_dbsearch-s14-db6fu.txt
    echo "dbsearch-20261012 session 14 (local losat-baselines-4d) since $(date '+%F %T'): node $node CPUs $CPUS and their siblings $SIBS ONLY (TIMED, membind $node): DB-6 follow-up blastp / blastx -db, ~25 min (db3/task/scripts/ws_db6_followup.sh); driver PID $$; removed at its end." > $M
    sleep 5
    if [ -n "$($CB overlap $ST $CPUS,$SIBS $M)" ]; then rm -f $M; M=""; else rm -f $Q; break; fi
  fi
  sleep 60
done
export PATH=/usr/bin:/bin RUST_BACKTRACE=0 TMPDIR=$D3/tmp
for v in $(env | grep -o '^LOSAT_[A-Z0-9_]*'); do unset "$v"; done
IFS=, read -r -a CL <<< "$CPUS"; MINE=$CPUS,$SIBS
cpus_for() { local IFS=,; echo "${CL[*]:0:$1}"; }
echo "$(date '+%F %T') start node $node CPUs $CPUS; load $(cut -d' ' -f1-3 /proc/loadavg)" | tee -a $S
[ -f $T ] || printf 'case\tarm\tthreads\trun\trc\twall_s\tuser_s\tsys_s\tmaxrss_kib\trows\tsha16\tbusy0\tbusy1\n' > $T
one() {  # case arm threads run envs... -- cmd...
  local c=$1 a=$2 t=$3 r=$4; shift 4; local envs=()
  while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
  local o=$O/$c.$a.n$t.$r b0 b1 rc
  b0=$($CB count 0.5 $node $MINE)
  env "${envs[@]}" numactl --physcpubind=$(cpus_for $t) --membind=$node python3 $B/scripts/s7_time.py $o.time "$@" -num_threads $t -outfmt 6 -out $o.out < /dev/null 2> $o.err; rc=$?
  b1=$($CB count 0.5 $node $MINE)
  read wall user sys rss < $o.time
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $c $a $t $r $rc $wall $user $sys $rss $(wc -l < $o.out) $(sha256sum $o.out | cut -c1-16) $b0 $b1 >> $T
}
three() {  # case threads armA envA cmdA... (arms: names in ARMS, envs in ENV_<name>, command CMD_<tool>)
  local c=$1 t=$2 r a
  for a in $ARMS; do one $c $a $t w $(eval echo \$ENV_$a) -- $(eval echo \$CMD_$a); done
  for r in 1 2 3; do
    local order=$ARMS; [ $r = 2 ] && order=$(echo $ARMS | tr ' ' '\n' | tac | tr '\n' ' ')
    for a in $order; do one $c $a $t $r $(eval echo \$ENV_$a) -- $(eval echo \$CMD_$a); done
  done
  echo "$(date '+%F %T') $c n$t: $(awk -F'\t' -v c=$c -v t=$t '$1==c && $3==t && $4!="w" {printf "%s%s=%s ", $2, $4, $6}' $T)" | tee -a $S
}
ARMS="ncbi def cut"
ENV_ncbi=X=1; ENV_def=X=1; ENV_cut=LOSAT_X_BPDBCUT=1
CMD_ncbi="$NB/blastp -query $W/in/yeast.faa -db $W/db/yeast.ncbi5/yeast"
CMD_def="$BIN blastp -query $W/in/yeast.faa -db $W/db/yeast.ncbi5/yeast"; CMD_cut=$CMD_def
three BP-yeast 16
ARMS="ncbi def flow"
ENV_flow=LOSAT_X_BXFLOW=1
CMD_ncbi="$NB/blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp"
CMD_def="$BIN blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp"; CMD_flow=$CMD_def
three BX-t10 8
three BX-t10 16
one BX-t10 timing 16 x LOSAT_TIMING=1 -- $BIN blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp
one BP-yeast timing 16 x LOSAT_TIMING=1 -- $BIN blastp -query $W/in/yeast.faa -db $W/db/yeast.ncbi5/yeast
echo "$(date '+%F %T') shas: $(awk -F'\t' 'NR>1 {print $1, $11}' $T | sort -u | awk '{a[$1]=a[$1] " " $2} END {for (c in a) printf "%s:%s; ", c, a[c]}')" | tee -a $S
echo "$(date '+%F %T') done" | tee -a $S

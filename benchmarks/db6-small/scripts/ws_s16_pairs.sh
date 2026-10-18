#!/bin/bash
# dbsearch-20261012 session 16, item 2: LOSAT_X_BXSUBJFLOW default ON (perf/blastx-db-subject-flow HEAD_SHA) -- paired
# on/off timings in one run (D1 of DEFAULT_ON_POLICY_20261017.md) and V-PERF of the default path against its base on dev.
# Questions and the minimal cases that answer them:
#   1. Does the flow win at every thread count of the DB-6 loss?  BX-t10 x Swiss-Prot -db n16 / n8 / n4 / n2:
#      off (LOSAT_X_BXSUBJFLOW=0) / default (flow) / NCBI.
#   2. Does it lose where BXFLOW (perf-bx #277, default soon) puts the last batches side by side, each running a flow?
#      cdna800 x yeast -db (-evalue 1e-5) n16 / n8 and cdna50 x yeast -db n16: off / default, each with and without
#      LOSAT_X_BXFLOW=1 LOSAT_X_BXSCANBF=1 (#276 / #277's default set, BXDIAG's 5-hit rule is not on dev), and NCBI.
#   3. Does it lose with fewer tasks than threads?  bench5 blastx -subject (MelaMJNV x PemoMJNVA, 115 subjects = 2 tasks) n8.
#   4. V-PERF: head (default = flow) against base (dev BASE_SHA) where the flow does not run or barely does: bench5 pair
#      n1 / n8, cdna50 x yeast -db n1.
# Builds HEAD_SHA and BASE_SHA on 4 idle CPUs (SHARED), then on 16 idle cores + siblings (TIMED; cpu_busy.py pick ..
# polite, overlap re-check) a warmup round and ROUNDS rounds per case, the order reversed in every even round.
#   ~/losat-bench/db3/task/scripts/ws_s16_pairs.sh HEAD_SHA BASE_SHA
set -u
SHA=$1; BASE=$2
B=$HOME/losat-bench; ST=$B/state; D3=$B/db3; R=$D3/repo; W=$D3/db6; O=$D3/s16/pairs; mkdir -p $O
CB="python3 $D3/task/scripts/cpu_busy.py"; NB=$B/tools/env_blast/bin
S=$O/summary.txt; T=$O/runs.tsv; M=""; MK=""
Q=$ST/NODE_QUEUED_dbsearch-s16-pairs.txt
trap 'rm -f "$Q" ${M:+"$M"} ${MK:+"$MK"}' EXIT
export PATH=$HOME/.cargo/bin:/usr/bin:/bin TMPDIR=$D3/tmp RUSTUP_TOOLCHAIN=1.92.0 RUST_BACKTRACE=0
NH=$D3/bin/LOSAT-s16-${SHA:0:8}; NBASE=$D3/bin/LOSAT-s16-${BASE:0:8}
for s in $SHA $BASE; do
  N=$D3/bin/LOSAT-s16-${s:0:8}
  [ -x $N ] && continue
  while true; do
    pick=$($CB pick 3 4 $ST); [ -n "$pick" ] && break; sleep 60
  done
  read -r node CPUS SIBS <<< "$pick"
  MK=$ST/NODE${node}_SHARED_dbsearch-s16-build.txt
  echo "dbsearch-20261012 session 16 (local losat-baselines-92) since $(date '+%F %T'): node $node CPUs $CPUS ONLY (UNTIMED): cargo build --release of $s (~6 min); driver PID $$; removed at its end." > $MK
  git -C $R fetch -q origin perf/blastx-db-subject-flow dev && git -C $R checkout -q --detach "$s" || { echo "FAILED checkout $s" | tee -a $S; exit 1; }
  touch $R/LOSAT/src/main.rs
  (cd $R/LOSAT && CARGO_BUILD_JOBS=4 numactl --physcpubind=$CPUS --membind=$node cargo build --release --locked --bin LOSAT --target-dir $D3/target_db6) > $O/build-${s:0:8}.log 2>&1 || { echo "FAILED build $s" | tee -a $S; exit 1; }
  grep -q "Compiling LOSAT" $O/build-${s:0:8}.log || { echo "FAILED: no Compiling line ($s)" | tee -a $S; exit 1; }
  cp $D3/target_db6/release/LOSAT $N
  echo "$(date '+%F %T') built $N $(sha256sum $N | cut -c1-16) from $(git -C $R rev-parse HEAD) on node $node CPUs $CPUS" | tee -a $S
  rm -f $MK; MK=""
done
echo "dbsearch-20261012 session 16 (local losat-baselines-92) queued $(date '+%F %T'): TIMED ~25 min on 16 idle cores + siblings (db3/task/scripts/ws_s16_pairs.sh); driver PID $$" > $Q
want=16
while true; do
  pick=$($CB pick 3 $want $ST polite $Q)
  if [ -n "$pick" ]; then
    read -r node CPUS SIBS <<< "$pick"
    M=$ST/NODE${node}_IN_USE_dbsearch-s16-pairs.txt
    echo "dbsearch-20261012 session 16 (local losat-baselines-92) since $(date '+%F %T'): node $node CPUs $CPUS and their siblings $SIBS ONLY (TIMED, membind $node): blastx on/off pairs + V-PERF, ~25 min (db3/task/scripts/ws_s16_pairs.sh); driver PID $$; removed at its end." > $M
    sleep 5
    if [ -n "$($CB overlap $ST $CPUS,$SIBS $M)" ]; then rm -f $M; M=""; else rm -f $Q; break; fi
  fi
  sleep 60
done
for v in $(env | grep -o '^LOSAT_[A-Z0-9_]*'); do unset "$v"; done
IFS=, read -r -a CL <<< "$CPUS"; MINE=$CPUS,$SIBS
cpus_for() { local IFS=,; echo "${CL[*]:0:$1}"; }
echo "$(date '+%F %T') start node $node CPUs $CPUS ($want cores); load $(cut -d' ' -f1-3 /proc/loadavg)" | tee -a $S
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
ROUNDS=3
rounds() {  # case threads (arms: names in ARMS, envs in ENV_<name>, command CMD_<name>)
  local c=$1 t=$2 r a order
  for a in $ARMS; do one $c $a $t w $(eval echo \$ENV_$a) -- $(eval echo \$CMD_$a); done
  for r in $(seq 1 $ROUNDS); do
    order=$ARMS; [ $((r % 2)) = 0 ] && order=$(echo $ARMS | tr ' ' '\n' | tac | tr '\n' ' ')
    for a in $order; do one $c $a $t $r $(eval echo \$ENV_$a) -- $(eval echo \$CMD_$a); done
  done
  echo "$(date '+%F %T') $c n$t: $(awk -F'\t' -v c=$c -v t=$t '$1==c && $3==t && $4!="w" {printf "%s%s=%s ", $2, $4, $6}' $T)" | tee -a $S
}
F=$R/LOSAT/tests/fasta; BX=$B/work/bx; Y=$BX/db/scer_pep
# 1. BX-t10
L="$NH blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp"
ARMS="ncbi off def"; ENV_ncbi=X=1; ENV_off=LOSAT_X_BXSUBJFLOW=0; ENV_def=X=1
CMD_ncbi="$NB/blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp"; CMD_off=$L; CMD_def=$L
for t in 16 8 4 2; do rounds BX-t10 $t; done
# 2. many query batches x yeast -db, with and without BXFLOW + BXSCANBF
for c in cdna800 cdna50; do
  L="$NH blastx -query $BX/in/$c.fa -db $Y -evalue 1e-5"
  ARMS="ncbi off def offF defF"; ENV_ncbi=X=1; ENV_off=LOSAT_X_BXSUBJFLOW=0; ENV_def=X=1
  ENV_offF="LOSAT_X_BXSUBJFLOW=0 LOSAT_X_BXFLOW=1 LOSAT_X_BXSCANBF=1"; ENV_defF="LOSAT_X_BXFLOW=1 LOSAT_X_BXSCANBF=1"
  CMD_ncbi="$NB/blastx -query $BX/in/$c.fa -db $Y -evalue 1e-5"; CMD_off=$L; CMD_def=$L; CMD_offF=$L; CMD_defF=$L
  rounds $c 16
  [ $c = cdna800 ] && rounds $c 8
done
# 3. bench5 blastx -subject (2 tasks) n8
L="$NH blastx -query $F/MelaMJNV.fasta -subject $F/PemoMJNVA.faa"
ARMS="ncbi off def"; ENV_ncbi=X=1; ENV_off=LOSAT_X_BXSUBJFLOW=0; ENV_def=X=1
CMD_ncbi="$NB/blastx -query $F/MelaMJNV.fasta -subject $F/PemoMJNVA.faa"; CMD_off=$L; CMD_def=$L
rounds bench5 8
# 4. V-PERF: base (dev) against head (default)
ROUNDS=5
ARMS="base head"; ENV_base=X=1; ENV_head=X=1
CMD_base="$NBASE blastx -query $F/MelaMJNV.fasta -subject $F/PemoMJNVA.faa"; CMD_head="$NH blastx -query $F/MelaMJNV.fasta -subject $F/PemoMJNVA.faa"
rounds vperf-bench5 1
rounds vperf-bench5 8
CMD_base="$NBASE blastx -query $BX/in/cdna50.fa -db $Y -evalue 1e-5"; CMD_head="$NH blastx -query $BX/in/cdna50.fa -db $Y -evalue 1e-5"
rounds vperf-cdna50 1
for t in 16 8; do
  one BX-t10 timing $t x LOSAT_TIMING=1 -- $NH blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp
  echo "$(date '+%T') timing n$t: $(grep -h 'blastx_subject_flow\|blastx_batch_alone:\|blastx_run:' $O/BX-t10.timing.n$t.x.err | cut -c10-300 | paste -sd' ')" >> $S
done
echo "$(date '+%F %T') done; load $(cut -d' ' -f1-3 /proc/loadavg)" | tee -a $S

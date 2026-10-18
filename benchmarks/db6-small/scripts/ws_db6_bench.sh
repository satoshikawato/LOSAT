#!/bin/bash
# dbsearch-20261012 session 14, DB-6 small tier (DESIGN_db6.md). Question: LOSAT -db (default = strict set) against
# NCBI BLAST+ 2.17.0 -db at 1 / 8 / 16 threads: median wall time and min-max of 3 timed runs after 1 warmup (order of
# the two tools alternates by round), max RSS, and SHA-256 of every output against NCBI -num_threads 1. Database
# builds (makeblastdb default and -blastdb_version 4, LOSAT makedb) timed the same way, apart from the searches.
# Owner 2026-10-17 08:30: a TIMED run needs only its own physical cores and their HT siblings idle; every run records
# how many other CPUs of the node were busy (> 50 %) at its start and end, and a timed run where that rose by 3 or more
# is recorded as disturbed (run label "<r>d<try>") and repeated (up to 2 more tries; the last is kept as it is).
#   CPUS=c1,c2,.. SIBS=s1,s2,.. MEMNODE=n ~/losat-bench/db3/task/scripts/ws_db6_bench.sh BIN [STEP...]
# Steps: prep build MB-q100 MB-amp BN-amp BP-q1 BP-q100 BP-yeast BX-t10 summary (default: all, in this order).
set -u
B=$HOME/losat-bench; D3=$B/db3; P=$D3/build_root/dbsearch/probe; W=$D3/db6
BIN=$1; shift
STEPS=${*:-prep build MB-q100 MB-amp BN-amp BP-q1 BP-q100 BP-yeast BX-t10 summary}
CPUS=${CPUS:?}; SIBS=${SIBS:?}; MEMNODE=${MEMNODE:?}
CB="python3 $D3/task/scripts/cpu_busy.py"; MINE=$CPUS,$SIBS
NB=$B/tools/env_blast/bin
YEAST_GZ=$B/data/proteome_yeast/UP000002311_559292.fasta.gz
mkdir -p $W/in $W/db $W/out $W/dbtmp $D3/tmp
S=$W/summary.txt; T=$W/runs.tsv
step() { echo "$(date '+%F %T') $*" | tee -a $S; }
export PATH=/usr/bin:/bin RUST_BACKTRACE=0 TMPDIR=$D3/tmp
for v in $(env | grep -o '^LOSAT_[A-Z0-9_]*'); do unset "$v"; done
IFS=, read -r -a CPU_LIST <<< "$CPUS"
first=${CPU_LIST[0]}
cpus_for() { local IFS=,; echo "${CPU_LIST[*]:0:$1}"; }
TIME="python3 $B/scripts/s7_time.py"
[ -f $T ] || printf 'case\ttool\tthreads\trun\trc\twall_s\tuser_s\tsys_s\tmaxrss_kib\trows\tsha16\tload\tbusy0\tbusy1\n' > $T
row() {  # case tool threads run rc timefile outfile load busy0 busy1
  local wall=NA user=NA sys=NA rss=NA
  [ -f $7 ] || : > $7
  [ -s $6 ] && read wall user sys rss < $6
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $1 $2 $3 $4 $5 $wall $user $sys $rss \
    $(wc -l < $7) $(sha256sum $7 | cut -c1-16) $8 $9 ${10} >> $T
}

search() {  # case tool threads run cmd...
  local c=$1 tool=$2 t=$3 r=$4; shift 4
  local o=$W/out/$c.$tool.n$t.$r load b0 b1 rc try
  for try in 1 2 3; do
    b0=$($CB count 0.5 $MEMNODE $MINE); load=$(cut -d' ' -f1 /proc/loadavg)
    numactl --physcpubind=$(cpus_for $t) --membind=$MEMNODE $TIME $o.time "$@" -num_threads $t -outfmt 6 -out $o.out \
      < /dev/null 2> $o.err
    rc=$?
    b1=$($CB count 0.5 $MEMNODE $MINE)
    if [ "$r" != w ] && [ $try -lt 3 ] && [ $((b1 - b0)) -ge 3 ]; then
      row $c $tool $t ${r}d$try $rc $o.time $o.out $load $b0 $b1
      continue
    fi
    row $c $tool $t $r $rc $o.time $o.out $load $b0 $b1
    break
  done
}

# NCMD / LCMD: the two tools' commands without -num_threads / -outfmt / -out; LDB: LOSAT's -db on its own database.
case_run() {  # case
  local c=$1 t r
  for t in 1 8 16; do
    search $c ncbi $t w "${NCMD[@]}"
    search $c losat $t w "${LCMD[@]}"
    for r in 1 2 3; do
      if [ $r = 2 ]; then
        search $c losat $t $r "${LCMD[@]}"; search $c ncbi $t $r "${NCMD[@]}"
      else
        search $c ncbi $t $r "${NCMD[@]}"; search $c losat $t $r "${LCMD[@]}"
      fi
    done
    step "$c n$t: $(awk -F'\t' -v c=$c -v t=$t '$1==c && $3==t && $4!="w" {printf "%s%s=%s ", $2, $4, $6}' $T)"
  done
  search $c losatdb 8 x "${LCMD[@]:0:${#LCMD[@]}-2}" -db "$LDB"
  step "$c LOSAT on its own database (n8): $(awk -F'\t' -v c=$c '$1==c && $2=="losatdb" {print $11, "rc=" $5}' $T)"
}

mkdb_runs() {  # name dbtype fasta ext(p|n)
  local name=$1 dbt=$2 fa=$3 x=$4 r tool order
  for r in w 1 2 3; do
    if [ $r = 2 ]; then order="losat ncbi4 ncbi5"; else order="ncbi5 ncbi4 losat"; fi
    for tool in $order; do
      local d=$W/dbtmp/$name.$tool o=$W/out/build-$name.$tool.$r load b0 b1 rc
      rm -rf $d; mkdir -p $d
      b0=$($CB count 0.5 $MEMNODE $MINE); load=$(cut -d' ' -f1 /proc/loadavg)
      case $tool in
        ncbi5) numactl --physcpubind=$first --membind=$MEMNODE $TIME $o.time $NB/makeblastdb -in $fa -dbtype $dbt -out $d/$name > $o.log 2>&1;;
        ncbi4) numactl --physcpubind=$first --membind=$MEMNODE $TIME $o.time $NB/makeblastdb -in $fa -dbtype $dbt -out $d/$name -blastdb_version 4 > $o.log 2>&1;;
        losat) numactl --physcpubind=$first --membind=$MEMNODE $TIME $o.time $BIN makedb -in $fa -dbtype $dbt -out $d/$name > $o.log 2>&1;;
      esac
      rc=$?
      b1=$($CB count 0.5 $MEMNODE $MINE)
      row build-$name $tool 1 $r $rc $o.time $d/$name.${x}sq $load $b0 $b1
      if [ $r = 3 ]; then rm -rf $W/db/$name.$tool; mv $d $W/db/$name.$tool; fi
    done
  done
  step "build $name: $(awk -F'\t' -v c=build-$name '$1==c && $4!="w" {printf "%s%s=%s ", $2, $4, $6}' $T)"
  local f
  for f in ${x}sq ${x}hr ${x}in; do
    step "  $name.$f ncbi4 vs losat: $(cmp -l $W/db/$name.ncbi4/$name.$f $W/db/$name.losat/$name.$f 2>&1 | wc -l) differing bytes"
  done
}

for s in $STEPS; do case $s in
prep)
  step "start bin $BIN sha256 $(sha256sum $BIN | cut -c1-16) cpus $CPUS (siblings $SIBS) node $MEMNODE governor $(cat /sys/devices/system/cpu/cpu$first/cpufreq/scaling_governor); load $(cut -d' ' -f1-3 /proc/loadavg)"
  step "NCBI $($NB/blastp -version | head -1) / $($NB/makeblastdb -version | head -1)"
  gzip -dc $YEAST_GZ > $W/in/yeast.faa
  cp $P/16s/16s_all.fna $W/in/16s.fna; cp $P/sp/sp_all.faa $W/in/sp.faa
  cp $P/16s/q100.fna $W/in/16s_q100.fna; cp $P/sp/q1.faa $W/in/sp_q1.faa; cp $P/sp/q100.faa $W/in/sp_q100.faa
  python3 -I - $W/in/16s_q100.fna $W/in/amp100.fna $D3/db5x/in/q100.fna $W/in/t10.fna <<'PY'
import sys
def records(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line, []
        elif line:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)
q100, amp, tx, t10 = sys.argv[1:5]
with open(amp, "w") as f:  # letters 341-800 (1-based) of each full-length 16S query: an amplicon-like 460 nt slice
    for name, seq in records(q100):
        assert len(seq) >= 800, name
        s = seq[340:800]
        f.write(f"{name.split()[0]} slice 341-800\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
with open(t10, "w") as f:
    for i, (name, seq) in enumerate(records(tx)):
        if i == 10:
            break
        f.write(name + "\n" + "".join(seq[j:j + 60] + "\n" for j in range(0, len(seq), 60)))
PY
  head -n $(awk '/^>/{n++} n==11{print NR-1; exit}' $W/in/amp100.fna) $W/in/amp100.fna > $W/in/amp10.fna
  for f in $W/in/*; do
    step "input $(basename $f): $(grep -c '>' $f) records, $(grep -v '>' $f | tr -d '\n' | wc -c) letters, sha256 $(sha256sum $f | cut -c1-16)"
  done
  sha256sum $YEAST_GZ | tee -a $S
  ;;
build)
  mkdb_runs 16s nucl $W/in/16s.fna n
  mkdb_runs sp prot $W/in/sp.faa p
  mkdb_runs yeast prot $W/in/yeast.faa p
  ;;
MB-q100)
  NCMD=($NB/blastn -query $W/in/16s_q100.fna -db $W/db/16s.ncbi5/16s)
  LCMD=($BIN blastn -query $W/in/16s_q100.fna -db $W/db/16s.ncbi5/16s); LDB=$W/db/16s.losat/16s
  case_run MB-q100;;
MB-amp)
  NCMD=($NB/blastn -query $W/in/amp100.fna -db $W/db/16s.ncbi5/16s)
  LCMD=($BIN blastn -query $W/in/amp100.fna -db $W/db/16s.ncbi5/16s); LDB=$W/db/16s.losat/16s
  case_run MB-amp;;
BN-amp)
  search BN-pilot ncbi 16 p $NB/blastn -task blastn -query $W/in/amp10.fna -db $W/db/16s.ncbi5/16s
  search BN-pilot losat 16 p $BIN blastn -task blastn -query $W/in/amp10.fna -db $W/db/16s.ncbi5/16s
  pw=$(awk -F'\t' '$1=="BN-pilot" && $2=="ncbi" {print $6}' $T | tail -1)
  N=10
  for n in 100 50 25; do
    if awk -v w=$pw -v n=$n 'BEGIN { exit !(w / 10 * n * 12 <= 180) }'; then N=$n; break; fi
  done
  step "BN pilot (amp10, n16): $(awk -F'\t' '$1=="BN-pilot" {printf "%s=%s ", $2, $6}' $T)-> N=$N"
  head -n $(awk -v k=$((N + 1)) '/^>/{n++} n==k && !d {print NR-1; d=1; exit} END{if (!d) print NR}' $W/in/amp100.fna) $W/in/amp100.fna > $W/in/amp$N.bn.fna
  NCMD=($NB/blastn -task blastn -query $W/in/amp$N.bn.fna -db $W/db/16s.ncbi5/16s)
  LCMD=($BIN blastn -task blastn -query $W/in/amp$N.bn.fna -db $W/db/16s.ncbi5/16s); LDB=$W/db/16s.losat/16s
  case_run BN-amp$N;;
BP-q1)
  NCMD=($NB/blastp -query $W/in/sp_q1.faa -db $W/db/sp.ncbi5/sp)
  LCMD=($BIN blastp -query $W/in/sp_q1.faa -db $W/db/sp.ncbi5/sp); LDB=$W/db/sp.losat/sp
  case_run BP-q1;;
BP-q100)
  NCMD=($NB/blastp -query $W/in/sp_q100.faa -db $W/db/sp.ncbi5/sp)
  LCMD=($BIN blastp -query $W/in/sp_q100.faa -db $W/db/sp.ncbi5/sp); LDB=$W/db/sp.losat/sp
  case_run BP-q100;;
BP-yeast)
  NCMD=($NB/blastp -query $W/in/yeast.faa -db $W/db/yeast.ncbi5/yeast)
  LCMD=($BIN blastp -query $W/in/yeast.faa -db $W/db/yeast.ncbi5/yeast); LDB=$W/db/yeast.losat/yeast
  case_run BP-yeast;;
BX-t10)
  NCMD=($NB/blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp)
  LCMD=($BIN blastx -query $W/in/t10.fna -db $W/db/sp.ncbi5/sp); LDB=$W/db/sp.losat/sp
  case_run BX-t10;;
summary)
  python3 $D3/task/scripts/db6_summary.py $T | tee -a $S;;
*) step "unknown step $s";;
esac; done
step "done ($STEPS)"

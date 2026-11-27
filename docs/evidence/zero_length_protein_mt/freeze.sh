#!/bin/bash
# NCBI BLAST+ 2.17.0 -num_threads 1 outputs on probe.py's databases (its <work_dir>/db): a_0s (= the test's dba0)
# gave LOSAT/tests/fixtures/zero_length_mt; na_0s and ns_mid are the nucleotide databases of README.md.
# Usage: bash freeze.sh <NCBI bin> <probe work dir> <output dir>
set -u
N=$1; cd "$2/db" || exit 1
O=$3; mkdir -p "$O"
r() { name=$1; shift; prlimit --as=8589934592 --core=1 timeout 120 "$@" -num_threads 1 -outfmt "6 std slen" > "$O/$name.out" 2> "$O/$name.err"; echo $? > "$O/$name.rc"; }
r blastp.dba0        $N/blastp -query q.fa -db a_0s
r blastp_ug.dba0     $N/blastp -query q.fa -db a_0s -ungapped -comp_based_stats 0
r blastx.dba0        $N/blastx -query qn.fa -db a_0s
r blastx_ug.dba0     $N/blastx -query qn.fa -db a_0s -ungapped -comp_based_stats 0
r psiblast.dba0      $N/psiblast -query q.fa -db a_0s -num_iterations 2
for db in na_0s ns_mid; do
  r blastn_mega.$db $N/blastn -query qn.fa -db $db
  r blastn.$db      $N/blastn -task blastn -query qn.fa -db $db
  r tblastn.$db     $N/tblastn -query q.fa -db $db
  r tblastx.$db     $N/tblastx -query qn.fa -db $db
done

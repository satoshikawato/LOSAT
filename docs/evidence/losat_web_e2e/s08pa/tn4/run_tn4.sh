#!/bin/bash
# TN-4 check: time and peak RSS of huge gapped X-drops (pre-change d96412265,
# final 68154c73e, NCBI 2.17.0). Output compared byte for byte with NCBI.
set -u
A=/home/kawato/.cache/losat-web-gui-target/s08pa
F=$A/head-src/LOSAT/tests/fasta/outfmt0
W=$A/tn4final/out
mkdir -p $W
run() { # label prog subject timeout args...
  local label=$1 prog=$2 subj=$3 to=$4; shift 4
  local args="$*"
  for b in ncbi base final; do
    case $b in
      ncbi) cmd=(/home/kawato/micromamba/bin/$prog);;
      base) cmd=($A/bin/LOSAT-base $prog);;
      final) cmd=($A/bin/LOSAT-final $prog);;
    esac
    /usr/bin/time -f "%e %M" -o $W/$label.$b.time timeout $to "${cmd[@]}" -query $F/e2e_protein_query.faa -subject $F/$subj $args -outfmt 6 > $W/$label.$b.out 2> $W/$label.$b.err
    rc=$?
    t=$(tail -1 $W/$label.$b.time)
    printf '%s\t%s\t%s\trc=%s\t%s s / %s KB\n' "$label" "$b" "$args" "$rc" ${t% *} ${t#* }
  done
  for b in base final; do
    if cmp -s $W/$label.ncbi.out $W/$label.$b.out; then s=same; else s=DIFF; fi
    printf '%s\t%s\tstdout vs NCBI: %s (%s lines)\n' "$label" "$b" "$s" "$(wc -l < $W/$label.$b.out)"
  done
}
run tbs_final1e7 tblastn e2e_tblastn_subject.fna 300 -comp_based_stats 0 -xdrop_gap_final 1e7
run tbs_final1e8 tblastn e2e_tblastn_subject.fna 300 -comp_based_stats 0 -xdrop_gap_final 1e8
run tbs_prelim1e8 tblastn e2e_tblastn_subject.fna 300 -xdrop_gap 1e8
run many_final1e7 tblastn e2e_many_subject.fna 300 -comp_based_stats 0 -xdrop_gap_final 1e7

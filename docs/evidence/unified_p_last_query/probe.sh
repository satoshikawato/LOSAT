#!/bin/bash
# Evidence of PD-LOSAT-NCBI-DEFECTS exception 16 and of the rule-2 items of version 1.54:
# BLASTP unified P-values (-comp_based_stats 1u, 2u, 3u) in NCBI BLAST+ 2.17.0.
#
# usage: probe.sh NCBI_BIN REPO WORK [LOSAT]
#   NCBI_BIN  the directory of NCBI BLAST+ 2.17.0's blastp and makeblastdb
#   REPO      a checkout of this repository (the inputs are LOSAT/tests/fasta/*.faa)
#   WORK      a directory for the runs (made if missing)
#   LOSAT     optional: a losat binary, compared with NCBI's reports (cases E and G)
# Writes WORK/results.tsv: case, run, threads, sha256 (12 digits) of the rows of the named
# query, number of those rows. Every run is capped in address space and time.
set -u
B=$1; R=$2; W=$3; L=${4:-}
F=$R/LOSAT/tests/fasta
mkdir -p "$W" && cd "$W" || exit 1
cap() { prlimit --core=1 --as=4000000000 timeout 300 "$@"; }
if ! cap "$B/makeblastdb" -in "$F/AP027131.faa" -dbtype prot -out db > makeblastdb.log 2>&1; then
    echo "makeblastdb failed"; exit 1
fi
python3 - "$F" <<'PY'
import os, sys
F = sys.argv[1]
def records(name):
    out, head, seq = [], None, []
    for line in open(os.path.join(F, name)):
        line = line.strip()
        if line.startswith('>'):
            if head: out.append((head, ''.join(seq)))
            head, seq = line, []
        elif line:
            seq.append(line)
    if head: out.append((head, ''.join(seq)))
    return out
def get(name, acc):
    return next(r for r in records(name) if r[0][1:].split()[0] == acc)
def write(name, items):
    with open(name, 'w') as f:
        for head, seq in items: f.write(head + '\n' + seq + '\n')
qs = get('AP027132.faa', 'BDV03099.1')      # 104 residues
qm = get('SicyWSV.faa', 'BDT63074.1')       # 897 residues
other = get('SicyWSV.faa', 'BDT63079.1')    # another protein, cut to qm's length below
copy = lambda r: ('>COPY copy of ' + r[0][1:].split()[0], r[1])
write('QS.fa', [qs]); write('QM.fa', [qm])
write('QM_copy.fa', [qm, copy(qm)])
write('QM_copy_QS.fa', [qm, copy(qm), qs])
write('QS_QM.fa', [qs, qm])
write('QS_len.fa', [qs, ('>LEN another protein of the length of ' + qm[0][1:].split()[0], other[1][:len(qm[1])])])
write('QS_copy.fa', [qs, copy(qs)])
write('QS_QM_copy.fa', [qs, qm, copy(qm)])
write('QS_allx.fa', [qs, ('>ALLX fifty X', 'X' * 50)])
print('QS', qs[0][1:].split()[0], len(qs[1]), 'QM', qm[0][1:].split()[0], len(qm[1]))
PY
QS=BDV03099.1; QM=BDT63074.1
: > results.tsv
rows_of() { awk -F'\t' -v q="$2" '$1 == q' "$1"; }
record() {  # case run threads file query
    local n s
    n=$(rows_of "$4" "$5" | wc -l)
    s=$(rows_of "$4" "$5" | sha256sum | cut -c1-12)
    printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$s" "$n" >> results.tsv
}
ncbi() {  # out threads query-file options...
    local out=$1 t=$2 q=$3; shift 3
    cap "$B/blastp" -db db -query "$q" -outfmt 6 -num_threads "$t" "$@" > "$out" 2> "$out.err"
    echo $? > "$out.rc"
}
U="-comp_based_stats 2u -evalue 1000"

# A. One query: the report of the batch's last query changes from run to run.
for t in 1 4; do for i in 1 2 3 4 5 6 7 8 9 10; do
    ncbi A.t$t.r$i.out $t QM.fa $U; record A.last_query.QM r$i $t A.t$t.r$i.out $QM
done; done

# B. With a copy after it, the report of the same query is the same in every run.
for i in 1 2 3; do ncbi B.t1.r$i.out 1 QM_copy.fa $U; record B.with_copy.QM r$i 1 B.t1.r$i.out $QM; done
for t in 4 8; do ncbi B.t$t.out $t QM_copy.fa $U; record B.with_copy.QM r1 $t B.t$t.out $QM; done
ncbi B.third.out 1 QM_copy_QS.fa $U; record B.with_copy_and_a_third_query.QM r1 1 B.third.out $QM

# C. A query with a next query: the report follows the next query's length.
ncbi C.qm.out 1 QS_QM.fa $U;   record C.next_is_QM.QS r1 1 C.qm.out $QS
ncbi C.qm2.out 1 QS_QM.fa $U;  record C.next_is_QM.QS r2 1 C.qm2.out $QS
ncbi C.qm8.out 8 QS_QM.fa $U;  record C.next_is_QM.QS r1 8 C.qm8.out $QS
ncbi C.len.out 1 QS_len.fa $U; record C.next_is_another_protein_of_QMs_length.QS r1 1 C.len.out $QS
ncbi C.copy.out 1 QS_copy.fa $U; record C.next_is_a_copy.QS r1 1 C.copy.out $QS

# D. Control without unified P-values (mode 2): the report of a query alone equals its report
#    with a copy after it, so the copy changes nothing else.
ncbi D.alone.out 1 QM.fa -comp_based_stats 2 -evalue 1000;     record D.mode2_alone.QM r1 1 D.alone.out $QM
ncbi D.copy.out 1 QM_copy.fa -comp_based_stats 2 -evalue 1000; record D.mode2_with_copy.QM r1 1 D.copy.out $QM

# F. A next query without a search space: every e-value of the first query is 0.
ncbi F.allx.out 1 QS_allx.fa $U; record F.next_has_no_search_space.QS r1 1 F.allx.out $QS
printf 'F.next_has_no_search_space.QS\te-values\t1\t%s\t-\n' \
    "$(rows_of F.allx.out $QS | cut -f11 | sort -u | tr '\n' ',' )" >> results.tsv
ncbi F.ssp.out 1 QS_allx.fa $U -searchsp 100000000; record F.next_has_no_search_space.searchsp_1e8.QS r1 1 F.ssp.out $QS

# G. A full hit list against tests/db_search/paja: C = BDT63114.1 (86 residues, tests/fasta/x5/c.faa)
#    with -max_target_seqs 1 has one match with an e-value of at most 0.002 and, after it, matches
#    with a preliminary e-value above 0.01, which NCBI's redo ends early once every list of the
#    batch is full (BlastCompo_EarlyTermination). With a copy after it (c_c.faa) the copy's list is
#    empty while C's matches are redone, none is ended early, and C's report is the same in every
#    run; alone, C is the batch's last query and its report changes from run to run.
GD=$R/LOSAT/tests/db_search/paja; QC=BDT63114.1; GU="-max_target_seqs 1 -comp_based_stats 2u"
paja() {  # out query-file: NCBI -num_threads 1 against tests/db_search/paja
    cap "$B/blastp" -db "$GD" -query "$2" -outfmt 6 -num_threads 1 $GU > "$1" 2> "$1.err"
    echo $? > "$1.rc"
}
for i in 1 2 3 4 5 6 7 8 9 10; do
    paja G.copy.r$i.out "$F/x5/c_c.faa";  record G.full_list_with_copy.QC r$i 1 G.copy.r$i.out $QC
    paja G.alone.r$i.out "$F/x5/c.faa";   record G.full_list_last_query.QC r$i 1 G.alone.r$i.out $QC
done

# E. LOSAT: the last query of a batch as NCBI's report with a copy; a query with a next query
#    as NCBI's; in case G's search, C as NCBI's report with the copy.
if [ -n "$L" ]; then
    losat() {  # out threads query-file options...
        local out=$1 t=$2 q=$3; shift 3
        cap "$L" blastp -db db -query "$q" -outfmt 6 -num_threads "$t" "$@" > "$out" 2> "$out.err"
        echo $? > "$out.rc"
    }
    for t in 1 8; do
        losat E.one.t$t.out $t QM.fa $U;      record E.losat.one_query.QM r1 $t E.one.t$t.out $QM
        losat E.two.t$t.out $t QS_QM.fa $U;   record E.losat.QS_then_QM.QS r1 $t E.two.t$t.out $QS
        record E.losat.QS_then_QM.QM r1 $t E.two.t$t.out $QM
    done
    ncbi E.ncbi.out 1 QS_QM_copy.fa $U
    record E.ncbi.QS_QM_copy.QS r1 1 E.ncbi.out $QS
    record E.ncbi.QS_QM_copy.QM r1 1 E.ncbi.out $QM
    losat E.allx.out 1 QS_allx.fa $U
    printf 'E.losat.next_has_no_search_space\texit status\t1\t%s\t%s\n' "$(cat E.allx.out.rc)" "$(wc -l < E.allx.out)" >> results.tsv
    losat E.ssp.out 1 QS_allx.fa $U -searchsp 100000000; record E.losat.next_has_no_search_space.searchsp_1e8.QS r1 1 E.ssp.out $QS
    for t in 1 8; do
        cap "$L" blastp -db "$GD" -query "$F/x5/c.faa" -outfmt 6 -num_threads "$t" $GU > G.losat.t$t.out 2> G.losat.t$t.out.err
        echo $? > G.losat.t$t.out.rc
        record G.losat.full_list.QC r1 $t G.losat.t$t.out $QC
    done
fi
bad=$(cat ./*.rc | grep -vc '^0$' || true)
echo "runs with a non-zero exit status: $bad (LOSAT's refusal in case E is one)"
echo PROBE_END

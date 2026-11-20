#!/usr/bin/env python3
"""Line RG, series Q, Q1c (rows db.pd.q1c.*, db.ps.q1c.*): -parse_deflines query inputs: two queries of one Seq-id, and a query whose id a reported database sequence
holds, each in one query batch and in two (BATCH_SIZE=1). usage (from LOSAT/): python3 tests/db_search/q1c/make_q1c_inputs.py tests/db_search tests/db_search/q1c"""
import os
import random
import sys

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)


def recs(path):
    r, name, seq = [], None, []
    for line in open(os.path.join(src, path)):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                r.append((name, ''.join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line)
    r.append((name, ''.join(seq)))
    return r


def write(name, records):
    with open(os.path.join(out, name), 'w') as f:
        for d, s in records:
            f.write(f'>{d}\n')
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + '\n')


rnd = random.Random(20261113)


def rand_nt(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


# A protein of no database here (WSSV TK-TMK, the first record of Q1's dup.faa).
KIN = ('MQLILSHHLTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEVIEGCTTVSSIKELISVLPEVVDVILIDEGQFFTD'
       'LVLVNRLADKGKRIVIAALDGTSDQQMFSPIHKLLPYTNSIVKLASKCMICKIDTKEAPFTVRFGNDNDNNVICVGG')
pairs = {
    # program key: (id that a reported database sequence holds, unrelated letters, hitting letters)
    'p': ('XP_000001.1 unrelated', KIN[:140], recs('pd_clash_q.faa')[0][1]),
    'n': ('AB000004.1 unrelated', rand_nt(480), recs('pd_clash_n_q.fna')[0][1]),
    'x': ('XP_000001.1 unrelated', rand_nt(420), recs('pd_clash_q.fna')[0][1]),
    't': ('tpe|BN000025.1| unrelated', rand_nt(420), recs('pd_clash_x_q.fna')[0][1]),
    'w': ('gb|AF000082.1| unrelated', rand_nt(660), recs('ids/nweak_qa2.fna')[0][1]),
}
for k, (idline, unrel, hit) in pairs.items():
    ext = 'faa' if k == 'p' else 'fna'
    write(f'{k}_after.{ext}', [(idline, unrel), ('q_hit hitting letters', hit)])
    write(f'{k}_before.{ext}', [('q_hit hitting letters', hit), (idline, unrel)])
# Two queries of one id with other letters, each with hits: 150 letters of two records of the database's own
# query files (blastp/psiblast bslice, tblastn tslice: proteins; blastn nslice, blastx bslice, tblastx xslice: nucleotides).
def two(path, n, a=0, b=1):
    r = recs(path)
    return [('lcl|dup2 first', r[a][1][:n]), ('lcl|dup2 second', r[b][1][:n])]
write('dup2_bp.faa', [('lcl|dup2 first', recs('bslice.faa')[0][1][:150]), ('lcl|dup2 second', recs('bslice.faa')[15][1][:150])])
write('dup2_tn.faa', two('tslice_q.faa', 150))
write('dup2_bn.fna', [('lcl|dup2 first', recs('nslice.fna')[0][1][1000:1400]), ('lcl|dup2 second', recs('nslice.fna')[2][1][1000:1400])])
write('dup2_bx.fna', two('bslice_q.fna', 450))
write('dup2_tx.fna', [('lcl|dup2 first', recs('xslice.fna')[0][1][1000:1450]), ('lcl|dup2 second', recs('xslice.fna')[2][1][500:950])])
# A name with and without an accession (blastp), two queries.
write('p_name.faa', [('gb|AAA00017.1|NAMEQ first', KIN[:140]), ('gb||NAMEQ second', recs('pd_clash_q.faa')[0][1])])
# Line WQA D10 round 4 (2026-11-20): the near inputs of the dup2 files, whose second record's letters are a line that
# NCBI's reader refuses (NCBI stops at that query batch after the reports of the batches before, as LOSAT's refusal of
# the second lcl|dup2 does: the stdout of rows db.pd.q1c.*.dup2stop.*)...
BAD_LINE = 'ACGTAC%%%%'
for name in ('dup2_bp.faa', 'dup2_tn.faa', 'dup2_bn.fna', 'dup2_bx.fna', 'dup2_tx.fna'):
    first, (second, _) = recs(os.path.join(os.path.relpath(out, src), name))
    stem, ext = name.rsplit('.', 1)
    write(f'{stem}_stop.{ext}', [first, (second, BAD_LINE)])
# ...and errors of NCBI's in the search of a batch before the second q1 (rows db.ps.q1c.kastop.*, db.pd.q1c.*.avgstop.*):
# psiblast's first query of only X (no Karlin-Altschul parameters), and subjects of 2, 0 and 0 letters with -searchsp 1000
# (the average subject length is too short), with a first query that fills the first batch of blastp and blastx.
rnd2 = random.Random(20261120)


def rand_aa(n):
    return ''.join(rnd2.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def rand_nt2(n):
    return ''.join(rnd2.choice('ACGT') for _ in range(n))


P = 'WNQQGYFFYHSSPALDGYTNTQRWQFAGQPMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ'
write('ka_ps.faa', [('q1', 'X' * 30), ('q2', P), ('q1', P)])
write('ka_ps_s.faa', [('s1', P)])
write('avg_s.faa', [('s1', 'MK'), ('s2', ''), ('s3', '')])
write('avg_bp.faa', [('q1', rand_aa(10500)), ('q2', rand_aa(300)), ('q1', rand_aa(300))])
write('avg_bx.fna', [('q1', rand_nt2(10500)), ('q2', rand_nt2(900)), ('q1', rand_nt2(900))])
# ...and the batch of the second q1 itself (line WQA D10 round 5, rows db.pd.q1c.*.avgdup.*, .avgdup2.*, .avglong.*,
# .kadup.*): subjects too short for -searchsp 1000 give NCBI's error of the batch's set-up whatever the repeat's letters
# (q1, q2, q1 in one batch; with BATCH_SIZE=30 a first batch of a q1 of only X, then the second q1 with q2), but in BLASTN
# with a repeat longer than its first record (NCBI reads past the first's letters), and with a first q1 of only N NCBI
# sets up no valid context and exits 0, so LOSAT refuses.
rnd3 = random.Random(202611205)


def rand_aa3(n):
    return ''.join(rnd3.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def rand_nt3(n):
    return ''.join(rnd3.choice('ACGT') for _ in range(n))


write('avg_s.fna', [('s1', 'AC'), ('s2', ''), ('s3', '')])
write('avgdup_bp.faa', [('q1', rand_aa3(300)), ('q2', rand_aa3(300)), ('q1', rand_aa3(300))])
write('avgdup_bn.fna', [('q1', rand_nt3(900)), ('q2', rand_nt3(900)), ('q1', rand_nt3(900))])
write('avgdup2_bp.faa', [('q1', 'X' * 30), ('q1', rand_aa3(20)), ('q2', rand_aa3(40))])
write('avglong_bp.faa', [('q1', rand_aa3(20)), ('q2', rand_aa3(100)), ('q1', rand_aa3(300))])
write('avglong_bn.fna', [('q1', rand_nt3(60)), ('q2', rand_nt3(300)), ('q1', rand_nt3(900))])
write('kadup_bx.fna', [('q1', 'N' * 90), ('q2', 'N' * 120), ('q1', rand_nt3(60))])
# NCBI's reader counts a repeated record into its batch by the letters of its first record (line WQA D10 round 5r, rows
# .avgearly.*, .avgmid.*). BATCH_SIZE=600: a first batch of a q1 of only X or N (not searched), then q2 and a second q1 of
# fewer letters that end the input: NCBI's batch is LOSAT's, and its set-up fails for the subjects (blastn, which reads the
# repeat's own interval in letters LOSAT does not hold, refuses). BATCH_SIZE=55: q1 of only X, a second q1 of 20 letters
# that NCBI counts as 30 and so ends its batch there (not searched: q1's warning twice), then q2: refused.
rnd4 = random.Random(202611207)


def rand_aa4(n):
    return ''.join(rnd4.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def rand_nt4(n):
    return ''.join(rnd4.choice('ACGT') for _ in range(n))


write('avgearly_bp.faa', [('q1', 'X' * 600), ('q2', rand_aa4(300)), ('q1', rand_aa4(200))])
write('avgearly_bn.fna', [('q1', 'N' * 600), ('q2', rand_nt4(300)), ('q1', rand_nt4(200))])
write('avgmid_bp.faa', [('q1', 'X' * 30), ('q1', rand_aa4(20)), ('q2', rand_aa4(40))])
# NCBI's reader counts a repeated record by its first record's letters at every cut of the batches (line WQA D10 round 6,
# review 4 of D10, M1; rows .g2a.*, .g2b.*, .avgmid.one.*): g2a, a q1 of 60 letters (180 nucleotides), a second q1 of 10
# (30) and a line that NCBI's reader refuses, with BATCH_SIZE=100 (300): NCBI's first batch is the two q1 (60 + 60),
# searched, then the reader's error, so LOSAT refuses that batch. g2b, the converse (10 then 60; BATCH_SIZE=50, 150): the
# batch counts 10 + 10 and reads the refused line, so the reader's error comes first. avgmid in one batch (q1 of only X or
# N, a second q1, q2; subjects too short for -searchsp 1000): NCBI's set-up of q1 and q2 fails for the subjects. And q1,
# q2, q1 in one batch against a database without its sequence file (review 4 M2, rows .nosq.*: arch_sprot_nopsq,
# arch_tbn_snucl_nonsq): NCBI's set-up passes and its first read of the database fails.
rnd5 = random.Random(202611208)


def rand_aa5(n):
    return ''.join(rnd5.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def rand_nt5(n):
    return ''.join(rnd5.choice('ACGT') for _ in range(n))


g2_aa, g2_nt = rand_aa5(200), rand_nt5(600)
bad = ('q3', 'ACGTAC%%%%')
write('g2_s.faa', [('s1', g2_aa)])
write('g2_s.fna', [('s1', g2_nt)])
write('g2a_bp.faa', [('q1', g2_aa[:60]), ('q1', g2_aa[60:70]), bad])
write('g2b_bp.faa', [('q1', g2_aa[:10]), ('q1', g2_aa[10:70]), bad])
write('g2a_bn.fna', [('q1', g2_nt[:180]), ('q1', g2_nt[180:210]), bad])
write('g2b_bn.fna', [('q1', g2_nt[:30]), ('q1', g2_nt[30:210]), bad])
write('avgmid_bn.fna', [('q1', 'N' * 90), ('q1', rand_nt5(60)), ('q2', rand_nt5(120))])
# A record that -query_loc skips joins NCBI's scope (rows .skipcount.*, .skiprep.*): with -query_loc 20-40 and
# BATCH_SIZE=100, a skipped q1 of 10 letters, q2 of 60 and a second q1 of 60 count 0 + 60 + 10, so the batch reads the
# refused line (NCBI's reader error); with -query_loc 15-50, blastx's second q1 after a skipped q1 of 10 nucleotides is
# set up in the skipped record's letters (NCBI: BLAST engine error, exit 3), so LOSAT refuses it.
write('skip_bp.faa', [('q1', g2_aa[100:110]), ('q2', g2_aa[:60]), ('q1', g2_aa[60:120]), bad])
write('skip_bn.fna', [('q1', g2_nt[300:310]), ('q1', g2_nt[:60])])
# Rows that fail where a cut of the batches loses the skipped record's place (line WQA D10 round 7, review 5 of D10, M1;
# rows .skipcount.*): a q1 that -query_loc skips counts 0 but joins the scope, so with BATCH_SIZE=220 (blastn, blastx,
# tblastx; -query_loc 40-100) skipcount_bn counts 0 + 180 + 30 and the batch reads the refused line (NCBI's reader
# error), where a cut that counted the skipped q1's own 30 nucleotides would end the batch at the second q1 (240) and
# refuse it; tblastn takes skip_bp.faa with BATCH_SIZE=75 (0 + 60 + 10, against 80). And rows that fail where the scope
# is not emptied after a batch against a database (review 5, M2; rows .h4.db.*): q1 of 60 letters, q2 of 40 (the first
# batch, BATCH_SIZE=100; 180, 120 and 300 nucleotides), then q1 of 10 and q1 of 60: in the emptied scope the second
# batch counts 10 + 10 and reads the refused line; a scope kept from the first batch would count 60 + 60 and refuse.
write('skipcount_bn.fna', [('q1', g2_nt[400:430]), ('q2', g2_nt[:180]), ('q1', g2_nt[180:360]), bad])
write('h4_bp.faa', [('q1', g2_aa[:60]), ('q2', g2_aa[60:100]), ('q1', g2_aa[100:110]), ('q1', g2_aa[110:170]), bad])
write('h4_bn.fna', [('q1', g2_nt[:180]), ('q2', g2_nt[180:300]), ('q1', g2_nt[300:330]), ('q1', g2_nt[330:510]), bad])
print(sorted(os.listdir(out)))

#!/usr/bin/env python3
"""Line RG, series Q, Q1c (rows db.pd.q1c.*, db.ps.q1c.*): -parse_deflines query inputs: two queries of one Seq-id, and a query whose id a reported database sequence
holds, each in one query batch and in two (BATCH_SIZE=1). usage (from LOSAT/): python3 tests/db_search/q1c/make_q1c_inputs.py tests/db_search tests/db_search/q1c"""
import os, random, sys

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
rand_nt = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))
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
print(sorted(os.listdir(out)))

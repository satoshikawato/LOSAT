#!/usr/bin/env python3
"""Writes the query files of the db.pd.* rows (-parse_deflines, line DI-7) from the query files of
the other rows: the same letters under deflines of each Seq-id form. Run from LOSAT/."""
import re

def records(path):
    out = []
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            out.append([line[1:], []])
        else:
            out[-1][1].append(line)
    return out

def write(path, deflines, source):
    recs = records(source)
    with open(path, 'w') as f:
        for defline, (_, lines) in zip(deflines, recs):
            f.write('>' + defline + '\n' + '\n'.join(lines) + '\n')

NUC = ['gi|42|ref|NM_000042.1| query gi ref', 'lcl|q2 query local', 'gnl|qdb|q3 query general',
       'gb|AB000004.3|AB4 query gb name', '', 'plain6 query plain']
NUC_NL = ['gi|42|ref|NM_000042.1| query gi ref', 'NC_000913.3 query bare', 'emb|X12345.1|N3 query emb',
          'tpg|BK000004.1| query tpa', 'dbj|AB000005.2| query ddbj']
PROT = ['gi|42|ref|NP_000042.1| query gi ref', 'lcl|q2 query local', 'gnl|qdb|q3 query general',
        'sp|P12345.1|ABC_X query sp', '', 'plain6 query plain', 'pdb|1ABC|A query pdb']
PROT_NL = ['gi|42|ref|NP_000042.1| query gi ref', 'YP_000043.1 query bare', 'sp|P12345.1|ABC_X query sp',
           'pdb|1ABC|A query pdb', 'tpg|DAA00005.1| query tpa']
D = 'tests/db_search/'
for name, source, nuc in (('n', 'nslice_q.fna', True), ('bx', 'bslice_q.fna', True), ('x', 'xslice_q.fna', True),
                          ('tn', 'tslice_q.faa', False), ('bp', 'bpqb_q.faa', False)):
    ext = 'fna' if nuc else 'faa'
    write(f'{D}pd_{name}_q.{ext}', NUC if nuc else PROT, D + source)
    write(f'{D}pd_{name}_nl_q.{ext}', NUC_NL if nuc else PROT_NL, D + source)
# Against the databases with parsed Seq-ids: every id differs from the database's.
write(D + 'pd_ids_q.fna', ['gi|142|ref|NM_000142.1| q one', 'lcl|q2 q two', 'plain3 q three', 'gnl|qdb|q4 q four'], D + 'ids/ids_q.fna')
write(D + 'pd_ids_q.faa', ['gi|142|ref|NP_000142.1| q one', 'lcl|q2 q two', 'plain3 q three', 'gnl|qdb|q4 q four'], D + 'ids/ids_q.faa')
# Refusals and NCBI's reader errors.
letters = ''.join(records(D + 'nslice_q.fna')[0][1])[:300]
amino = ''.join(records(D + 'tslice_q.faa')[0][1])[:120]
def one(path, deflines, seq):
    with open(path, 'w') as f:
        for defline in deflines:
            f.write(('>' + defline + '\n' if defline is not None else '') + seq + '\n')
for ext, seq, acc in (('fna', letters, 'NM_000042.1'), ('faa', amino, 'NP_000042.1')):
    one(f'{D}pd_dup_q.{ext}', ['dup1 a', 'dup1 b'], seq)
    one(f'{D}pd_ver_q.{ext}', [acc[:-2] + ' a', acc + ' b'], seq)
    one(f'{D}pd_ord_q.{ext}', ['gnl|BL_ORD_ID|1 a'], seq)
    one(f'{D}pd_range_q.{ext}', ['rq:10-40 t', 'r2:11-41 u'], seq)
    one(f'{D}pd_bars_q.{ext}', ['|q a'], seq)
    one(f'{D}pd_long_q.{ext}', ['Q' * 60 + ' a'], seq)
    one(f'{D}pd_comma_q.{ext}', ['a,b t'], seq)
    one(f'{D}pd_nodef_q.{ext}', [None], seq)
    one(f'{D}pd_empty_q.{ext}', ['q1 a', 'q2 b'], seq)
one(D + 'pd_other_q.faa', ['NM_000042.1 a'], amino)
# The reader's message names the molecule of the accession.
one(D + 'pd_warn_q.fna', ['NM_000042.1 a', 'plain2 b'], letters[:50] + '!' + letters[50:])
one(D + 'pd_warn_q.faa', ['NP_000042.1 a', 'plain2 b'], amino[:50] + '!' + amino[50:])
one(D + 'pd_nucprot_q.fna', ['NP_000042.1 a', 'AAA12345 b'], letters)
# A query with the id of a database sequence that it hits (ids/prot_v5, ids/nucl_v5).
write(D + 'pd_clash_q.faa', ['XP_000001.1 a'], D + 'ids/ids_q.faa')
write(D + 'pd_clash_q.fna', ['XP_000001.1 a'], D + 'ids/ids_q.fna')
write(D + 'pd_clash_x_q.fna', ['tpe|BN000025.1| a'], D + 'ids/ids_q.fna')
# BLASTN (megablast) against ids/nucl_v5: the query of the db.gis.* rows, which has hits.
write(D + 'pd_ids_n_q.fna', ['gi|142|ref|NM_000142.1| q one', 'lcl|q2 q two'], 'tests/seqid/fasta/qnucl.fa')
write(D + 'pd_clash_n_q.fna', ['AB000004.1 a'], 'tests/seqid/fasta/qnucl.fa')

#!/usr/bin/env python3
"""Writes the inputs of the rows `db.xj.*` (item K4j): BLASTX query and `-subject` titles with bytes that are not
UTF-8 (and the other byte classes of K4h), which NCBI BLAST+ 2.17.0 blastx reads and writes as they are
(`K/SURVEY_K4j.md` of the RG task folder).

Run from LOSAT/: python3 tests/db_search/make_xj_inputs.py. The files are committed. The `-db` rows search the
database `tests/db_search/xb_p` (K4h, `make_xb_inputs.py`: the same 60 residues).

Every record holds the same 60 residues (`PROT`; a nucleotide query their back-translation `NUC`), so every query
hits every subject. A title is `PREFIX` and one of the classes below around the byte(s) `X` (the classes and their
order are K4h's and the K4j probe's, so record N of a multi-record file has the class of the probe's q<N>/s<N>):

- mid: PREFIX ` a` X `b z`; end: PREFIX ` one` X; first: X PREFIX ` ab z`; id: PREFIX X ` rest` (X in the first word).

The files (`xj_<q|s><n|p>_<name>`; `n` = nucleotide `.fna`, `p` = protein `.faa`):

- `xj_qn_mq`: q1 `one` (clean), then one query per class (q2..q44); `xj_sp_ms`: one subject per class (s1..s43).
- `xj_qn_3m_l1x`: q1 and q3 clean, q2 of class `l1x` (Latin-1).
- `xj_qn_long`, `xj_sp_long`: titles of 60..90 (40..80) letters `a` and then an e-acute as UTF-8 or as Latin-1, for
  the byte-counted wrap of `Query=` and of the heading, and the description column's cut by bytes.
- `xj_sp_gd`: ASCII titles that NCBI's definition line generator changes (trailing `.,;~`, runs of spaces) and
  first words that start with `Subject_` (a first word `Subject_N` with a byte NCBI cannot decode is in the rows
  `db.t7.*` of `tests/db_search/t7`).
- `xj_qn_pd`, `xj_sp_pd`: the classes whose byte is in the first word, read with `-parse_deflines`.
- `xj_sp_empty`: two empty subjects (Latin-1 and an undecodable byte in their titles) and a clean one: the warning.
- `xj_rps_l1`: the first four queries of `tests/fasta/rps/nq.fna`, two of them with a Latin-1 byte and an undecodable byte in their
  titles, for rpstblastn (which reads queries with BLASTX's reader).
"""
from pathlib import Path

CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA',
         'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT', 'S': 'TCT', 'T': 'ACC', 'V': 'GTT',
         'W': 'TGG', 'Y': 'TAT'}
PROT = 'MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENTGMFRLGFPKFIG'
NUC = ''.join(CODON[a] for a in PROT)

# (name, bytes X, where), in the K4h / K4j probe order
CLASSES = [
    ('ascii', b'', 'mid'), ('u2', b'\xc3\xa9', 'mid'), ('u3', b'\xe2\x82\xac', 'mid'),
    ('u4', b'\xf0\x9f\x98\x80', 'mid'), ('l1x', b'\xe9', 'mid'), ('cont', b'\xa9', 'mid'), ('c80', b'\x80', 'mid'),
    ('x81', b'\x81', 'mid'), ('ovl2', b'\xc1\x81', 'mid'), ('c080', b'\xc0\x80', 'mid'),
    ('ovl3', b'\xe0\x81\x81', 'mid'), ('ovl3lt', b'\xe0\x80\xbc', 'mid'), ('ovl3c1', b'\xe0\x80\x81', 'mid'),
    ('sur', b'\xed\xa0\x80', 'mid'), ('f4big', b'\xf4\x90\x80\x80', 'mid'), ('f5', b'\xf5\x80\x80\x80', 'mid'),
    ('ff', b'\xff', 'mid'), ('c01', b'\x01', 'mid'), ('c1b', b'\x1b', 'mid'), ('del', b'\x7f', 'mid'),
    ('tab', b'\t', 'mid'), ('cr', b'\r', 'mid'), ('xml', b'&<>\'"', 'mid'), ('ent', b'&eacute;', 'mid'),
    ('ent1', b'&#1;', 'mid'), ('ent0', b'&#0;', 'mid'), ('entbig', b'&#128512;', 'mid'), ('entlt', b'&lt;', 'mid'),
    ('ent9', b'&#9;', 'mid'), ('ent10', b'&#10;', 'mid'), ('ent13', b'&#13;', 'mid'),
    ('u2end', b'\xc3\xa9', 'end'), ('u3end', b'\xe2\x82\xac', 'end'), ('l1end', b'\xe9', 'end'),
    ('contend', b'\xa9', 'end'), ('u4end', b'\xf0\x9f\x98\x80', 'end'), ('c01first', b'\x01', 'first'),
    ('c1bfirst', b'\x1b', 'first'),
    ('idu2', b'\xc3\xa9', 'id'), ('idl1', b'\xe9', 'id'), ('idx81', b'\x81', 'id'), ('idc01', b'\x01', 'id'),
    ('idff', b'\xff', 'id'),
]
CLS = {c[0]: c for c in CLASSES}
ID_CLASSES = ['idu2', 'idl1', 'idx81', 'idc01', 'idff']
GD_TITLES = [b's1 abc z.', b's2 abc z...', b's3 abc z,', b's4 abc z;', b's5 abc z~', b's6 abc z.,;~', b's7 a  b   z',
             b's8 a b  .', b's9 abc z .', b's10', b's11. rest', b'Subject_7 abc z', b'Subject_1 abc z',
             b'lcl|Subject_3 x', b's15 (abc) z', b's16 abc [z]']


def title(prefix: bytes, cls: str) -> bytes:
    _, x, where = CLS[cls]
    if where == 'mid':
        return prefix + b' a' + x + b'b z'
    if where == 'end':
        return prefix + b' one' + x
    if where == 'first':
        return x + prefix + b' ab z'
    return prefix + x + b' rest'


def fa(t: bytes, seq: str) -> bytes:
    return b'>' + t + b'\n' + seq.encode() + b'\n'


def long_title(prefix: bytes, k: int, x: bytes) -> bytes:
    return prefix + b' ' + b'a' * k + x + b' tail of the title words here'


def inputs(rps_query: bytes) -> dict[str, bytes]:
    out = {
        'xj_qn_mq.fna': fa(b'q1 one', NUC) + b''.join(fa(title(b'q%d' % (i + 2), c[0]), NUC)
                                                     for i, c in enumerate(CLASSES)),
        'xj_sp_ms.faa': b''.join(fa(title(b's%d' % (i + 1), c[0]), PROT) for i, c in enumerate(CLASSES)),
        'xj_qn_3m_l1x.fna': fa(b'q1 one', NUC) + fa(title(b'q2', 'l1x'), NUC) + fa(b'q3 three', NUC),
        'xj_sp_gd.faa': b''.join(fa(t, PROT) for t in GD_TITLES),
        'xj_qn_pd.fna': b''.join(fa(title(b'q%d' % (i + 1), c), NUC) for i, c in enumerate(ID_CLASSES)),
        'xj_sp_pd.faa': b''.join(fa(title(b's%d' % (i + 1), c), PROT) for i, c in enumerate(ID_CLASSES)),
        'xj_sp_empty.faa': b'>s1 a\xe9b z\n\n>s2 a\x81b z\n\n' + fa(b's3 subject', PROT),
    }
    xs = [b'\xc3\xa9', b'\xe9']
    out['xj_qn_long.fna'] = fa(b'q1 one', NUC) + b''.join(
        fa(long_title(b'q%d' % (i + 2), k, x), NUC)
        for i, (k, x) in enumerate((k, x) for x in xs for k in range(60, 91)))
    out['xj_sp_long.faa'] = b''.join(
        fa(long_title(b's%d' % (i + 1), k, x), PROT)
        for i, (k, x) in enumerate((k, x) for x in xs for k in range(40, 81)))
    # rpstblastn: the first three records of nq.fna, their titles with a Latin-1 byte and an undecodable byte
    chunks = rps_query.split(b'\n>')
    records = [chunks[0]] + [b'>' + chunk for chunk in chunks[1:4]]
    retitled = []
    for i, record in enumerate(records):
        head, _, rest = record.partition(b'\n')
        if i == 1:
            head = head.replace(b' bt ', b' b\xe9t ', 1)
        elif i == 2:
            head = head.replace(b' bt ', b' b\x81t ', 1)
        retitled.append(head + b'\n' + rest.rstrip(b'\n') + b'\n')
    out['xj_rps_l1.fna'] = b''.join(retitled)
    return out


def main() -> None:
    here = Path(__file__).resolve().parent
    rps_query = (here.parent / 'fasta' / 'rps' / 'nq.fna').read_bytes()
    target = here / 'xj'
    target.mkdir(exist_ok=True)
    for name, data in inputs(rps_query).items():
        (target / name).write_bytes(data)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Writes the inputs of the rows `db.xb.*` (item K4h): titles with bytes >= 0x80 or control bytes in the XML
family (outfmt 5, 13, 14, 15, 16), which NCBI BLAST+ 2.17.0 writes as its serial streams do (`SURVEY_K4h`).

Run from LOSAT/: python3 tests/db_search/make_xb_inputs.py. The files are committed; the databases are not made here:

    makeblastdb -in tests/db_search/xb/xb_sp_clean.faa -dbtype prot -out tests/db_search/xb_p \
        -title "xb proteins" -blastdb_version 4
    makeblastdb -in tests/db_search/xb/xb_sn_clean.fna -dbtype nucl -out tests/db_search/xb_n \
        -title "xb nucleotides" -blastdb_version 4

Every record holds the same 60 residues (`PROT`; the nucleotide records their back-translation `NUC`), so a query
hits the one subject of the database or of the `-subject` file, whatever the program. The title of a record is one
of the classes below, around the byte(s) `X`:

- `u2` `q1 a` X `b z` with X = C3 A9 (U+00E9, two-byte UTF-8): in the middle of the title.
- `u2end` `q1 one` X: the same bytes last. NCBI's `x_CleanAndCompress` drops a last byte >= 0x80, which leaves C3.
- `l1x` X = E9 (Latin-1 e-acute; NCBI converts the title to UTF-8 before it writes it).
- `c01first` 01 `q1 ab z`: a first byte below 0x20 (the reader keeps it; the XML stream writes `#` and posts 803.7).
- `ent1` `q1 a&#1;b z`: a character reference NCBI decodes to 0x01.
- `x81` `q1 a` 81 `b z`: a title whose encoding NCBI cannot guess (`GenerateDefline` throws; no title is written).
- `ovl3` `q1 a` E0 81 81 `b z`: an overlong three-byte sequence of `A`, which NCBI writes as `A`.

`xb_q<p|n>_nohit` is a query of 60 `X` (60 `N`) that hits nothing: a `-subject` title that no hit reaches is never written.

The names: `xb_q<p|n>_<class>` (one query, title `q1 ...`), `xb_q<p|n>_3m_u2` (queries q1 and q3 clean, q2 `u2`) and
`xb_q<p|n>_3f_u2` (q1 `u2`, q2 and q3 clean), `xb_q<p|n>_clean`; `xb_s<p|n>_<class>` (one subject, title `s1 ...`) and
`xb_s<p|n>_clean`. `p` files are `.faa`, `n` files `.fna`.
"""
from pathlib import Path

CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA',
         'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT', 'S': 'TCT', 'T': 'ACC', 'V': 'GTT',
         'W': 'TGG', 'Y': 'TAT'}
PROT = 'MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENTGMFRLGFPKFIG'
NUC = ''.join(CODON[a] for a in PROT)

# class -> (bytes X, where)
CLASSES = {
    'u2': (b'\xc3\xa9', 'mid'),
    'u2end': (b'\xc3\xa9', 'end'),
    'l1x': (b'\xe9', 'mid'),
    'c01first': (b'\x01', 'first'),
    'ent1': (b'&#1;', 'mid'),
    'x81': (b'\x81', 'mid'),
    'ovl3': (b'\xe0\x81\x81', 'mid'),
}


def title(prefix: bytes, cls: str) -> bytes:
    x, where = CLASSES[cls]
    if where == 'mid':
        return prefix + b' a' + x + b'b z'
    if where == 'end':
        return prefix + b' one' + x
    return x + prefix + b' ab z'


def record(t: bytes, seq: str) -> bytes:
    return b'>' + t + b'\n' + seq.encode() + b'\n'


def main() -> None:
    out = Path(__file__).resolve().parent / 'xb'
    out.mkdir(exist_ok=True)
    for mol, seq, ext in (('p', PROT, 'faa'), ('n', NUC, 'fna')):
        def put(name: str, data: bytes) -> None:
            (out / f'{name}.{ext}').write_bytes(data)
        put(f'xb_q{mol}_clean', record(b'q1 one', seq))
        put(f'xb_s{mol}_clean', record(b's1 subject', seq))
        put(f'xb_q{mol}_nohit', record(b'q1 none', ('X' if mol == 'p' else 'N') * 60))
        for cls in CLASSES:
            put(f'xb_q{mol}_{cls}', record(title(b'q1', cls), seq))
            put(f'xb_s{mol}_{cls}', record(title(b's1', cls), seq))
        put(f'xb_q{mol}_3m_u2', record(b'q1 one', seq) + record(title(b'q2', 'u2'), seq) + record(b'q3 three', seq))
        put(f'xb_q{mol}_3f_u2', record(title(b'q1', 'u2'), seq) + record(b'q2 two', seq) + record(b'q3 three', seq))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""The profiles and lists of the rps_odd_* fixture databases (series PR, item PS10e; build.sh): profiles whose
`.psq` letters (their `.pin` length L_pin) differ from their `.rps` columns L_rps, as makeprofiledb writes them
for a scored profile whose seq-data is nucleotide-coded or has more or fewer letters than its length.

Usage: python3 -I make_odd.py SMP_DIR OUT_DIR
SMP_DIR holds the PS10 oracle's c1.smp, c2.smp and c3.smp (LOSAT/tests/fixtures/makeprofiledb/in). OUT_DIR/in gets
the variants below (only the query's seq-data line or its id tag changes; makeprofiledb never reads the letters of
a scored profile and checks only Seq-inst.length against the columns) and one list l_<set>.lst per set, whose
lines name in/<profile>.smp (build.sh runs makeprofiledb -in in/l_<set>.lst -out rps_odd_<set> in OUT_DIR).

  c2_na    c2 with `seq-data iupacna` of 18 bases (its .psq bytes nucleotide-packed: L_pin 5, letters >= 28)
  c1_na    c1 with `seq-data iupacna` of 22 bases (.psq 143 18 159 6 210 210: L_pin 6, four letters >= 28)
  c1_k<n>  c1 under the id tag 91000 + n (n = 1..10): copies of c1 that score as c1
  c1_A     c1 with `seq-data iupacna` of 22 A (.psq 00 00 00 00 00 02: L_pin 6, every letter below 28)
  c1p<n>   c1 (22 columns) with its first n ncbieaa letters (n = 6, 10, 16, 21), or with W, WW after them (23, 24)
  c3p<n>   c3 (20 columns) with W, WW or WWW after its letters (n = 21, 22, 23)
"""
import os
import sys

C1, C2, C3 = 'IVLARIDDRLIHGQVATTWTKE', 'KTPAIRLLDKAKIPYELH', 'RLGVNIDHVATLRNARGTNY'

SETS = {
    'sc': ['c1', 'c2_na', 'c3'],
    'sc2v': ['c1', 'c2_na', 'c3'],  # built with -max_smp_vol 2
    'hl': ['c1_k%d' % k for k in range(1, 11)] + ['c2_na'],
    'aa1': ['c1_A', 'c2', 'c3'],
    'aa1s': ['c1_A'],
    'p6': ['c1p6', 'c2', 'c3'],
    'p10': ['c1p10', 'c2', 'c3'],
    'p16': ['c1p16', 'c2', 'c3'],
    'p21': ['c1p21', 'c2', 'c3'],
    'p23': ['c1p23', 'c2', 'c3'],
    'p24': ['c1p24', 'c2', 'c3'],
    'p23l': ['c1', 'c2', 'c3p21'],
    'p24l': ['c1', 'c2', 'c3p22'],
    'p23x3l': ['c1', 'c2', 'c3p23'],
    'na1': ['c1_na', 'c2', 'c3'],
}


def seq_line(coding: str, letters: str) -> str:
    return '        seq-data %s "%s"' % (coding, letters)


def main() -> None:
    smp_dir, out_dir = sys.argv[1], sys.argv[2]
    into = os.path.join(out_dir, 'in')
    os.makedirs(into, exist_ok=True)
    base = {name: open(os.path.join(smp_dir, name + '.smp')).read() for name in ('c1', 'c2', 'c3')}
    letters = {'c1': C1, 'c2': C2, 'c3': C3}
    smps = dict(base)

    def variant(name: str, of: str, coding: str, text: str) -> None:
        old = seq_line('ncbieaa', letters[of])
        assert base[of].count(old) == 1, (of, old)
        smps[name] = base[of].replace(old, seq_line(coding, text))

    variant('c2_na', 'c2', 'iupacna', 'GATTACAGGCTTAACGTC')
    variant('c1_A', 'c1', 'iupacna', 'A' * 22)
    variant('c1_na', 'c1', 'iupacna', 'GATTACAGGCTTAACGTCAGTC')
    for n in (6, 10, 16, 21):
        variant('c1p%d' % n, 'c1', 'ncbieaa', C1[:n])
    variant('c1p23', 'c1', 'ncbieaa', C1 + 'W')
    variant('c1p24', 'c1', 'ncbieaa', C1 + 'WW')
    for n, tail in ((21, 'W'), (22, 'WW'), (23, 'WWW')):
        variant('c3p%d' % n, 'c3', 'ncbieaa', C3 + tail)
    tag = 'tag id 237975'
    assert base['c1'].count(tag) == 1
    for k in range(1, 11):
        smps['c1_k%d' % k] = base['c1'].replace(tag, 'tag id %d' % (91000 + k))
    for name, text in smps.items():
        with open(os.path.join(into, name + '.smp'), 'w') as f:
            f.write(text)
    for name, members in SETS.items():
        with open(os.path.join(into, 'l_%s.lst' % name), 'w') as f:
            f.write(''.join('in/%s.smp\n' % member for member in members))


if __name__ == '__main__':
    main()

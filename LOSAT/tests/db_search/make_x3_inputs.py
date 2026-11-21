#!/usr/bin/env python3
"""Inputs of the rows db.x3.tn.* and x3.tblastn.* (line RG, X3: tblastn -max_intron_length).

Run from LOSAT/: python3 -I tests/db_search/make_x3_inputs.py
Writes tests/db_search/x3ex_q.faa (the query BDT61866.1 of tslice_q.faa) and tests/db_search/x3ex.fna
(ex_plus and ex_minus, its reverse complement): an exon-like subject, six pieces of the query
(residues [0:45], [60:105], [120:165], [180:225], [240:285], [300:345], each with 25 % substitutions
from a seeded generator, back-translated) separated by 91, 212, 450, 1201 and 4002 random bases
(gaps of about 30, 70, 150, 400 and 1334 letters in the translated subject) between flanks of 300,
so that the largest subject gap of the linking (-max_intron_length) decides which alignments link.
"""
import random

CODON = {
    'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC', 'Q': 'CAG', 'E': 'GAA',
    'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
    'P': 'CCG', 'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG', 'X': 'GCG',
}
AA = 'ACDEFGHIKLMNPQRSTVWY'


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = ''
        elif name:
            seqs[name] += line
    return seqs


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def wrap(s, n=60):
    return '\n'.join(s[i:i + n] for i in range(0, len(s), n))


def main():
    src, out = 'tests/db_search/tslice_q.faa', 'tests/db_search'
    rng = random.Random(20261120)
    q = read_fasta(src)['BDT61866.1']
    exons = [(0, 45), (60, 105), (120, 165), (180, 225), (240, 285), (300, 345)]
    introns = [91, 212, 450, 1201, 4002]

    def rnd(n):
        return ''.join(rng.choice('ACGT') for _ in range(n))

    parts = [rnd(300)]
    for i, (a, b) in enumerate(exons):
        prot = list(q[a:b])
        for j in range(len(prot)):
            if rng.random() < 0.25:
                prot[j] = rng.choice(AA)
        parts.append(''.join(CODON[c] for c in prot))
        if i < len(introns):
            parts.append(rnd(introns[i]))
    parts.append(rnd(300))
    plus = ''.join(parts)
    with open(f'{out}/x3ex_q.faa', 'w') as f:
        f.write(f'>BDT61866.1\n{wrap(q)}\n')
    with open(f'{out}/x3ex.fna', 'w') as f:
        f.write(f'>ex_plus exon-like copy\n{wrap(plus)}\n>ex_minus reverse complement\n{wrap(rc(plus))}\n')
    print(len(plus))


if __name__ == '__main__':
    main()

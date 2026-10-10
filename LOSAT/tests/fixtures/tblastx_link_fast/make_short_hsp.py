#!/usr/bin/env python3
"""Writes short_hsp.query.fna and short_hsp.subject.fna.

A TBLASTX pair in which a 4-residue HSP H is followed, in link_hsps list
order, by a long HSP J whose trimmed start lies beyond H's trimmed end on both
sequences.  NCBI offers H only the HSPs before it in the list, so H and J are
not linked; the large-gap sweep of LOSAT's default linking kernel
(sum_stats_linking/linking_incr.rs) finds J in its prefix-maximum tree and has
to fall back to NCBI's scan.

  query   frame +1: ... E G L D [W C H W] D I L  A V S  D L G E ...
  subject frame +1: ... C I D L [W C H W] L G D  A V S  L D I C ...

Every aligned pair outside WCHW and AVS scores -4 (BLOSUM62).  The word hits
"DWC/LWC" (offset 3 of the block) and "HWD/HWL" (offset 6) are three residues
apart, so they start a two-hit extension; it ends as the HSP WCHW (score 39),
H.  J is the query's frame +2 translation of the same nucleotides, starting
one residue (2 nt) before H, matched by a recoded copy further down the
subject (26 residues after extension).

The expected output is that of NCBI BLAST+ 2.17.0:
  tblastx -query short_hsp.query.fna -subject short_hsp.subject.fna -outfmt 6

usage: make_short_hsp.py OUTDIR
"""
import random
import sys

CODONS = {
    'F': ['TTT', 'TTC'], 'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 'I': ['ATT', 'ATC', 'ATA'],
    'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    '*': ['TAA', 'TAG', 'TGA'],
}
AA = {codon: aa for aa, codons in CODONS.items() for codon in codons}
HD = "WCHW"


def translate(nt):
    return ''.join(AA[nt[i:i + 3]] for i in range(0, len(nt) - 2, 3))


def encode(protein, pick=0):
    return ''.join(CODONS[aa][min(pick, len(CODONS[aa]) - 1)] for aa in protein)


def write(path, name, sequence):
    with open(path, "w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 70):
            handle.write(sequence[start:start + 70] + "\n")


def main():
    outdir = sys.argv[1]
    rng = random.Random(1)

    def background(length):
        return ''.join(rng.choice('ACGT') for _ in range(length))

    # Query: the frame +1 block, then a tail that keeps frame +2 open.
    block = encode("EGLD" + HD + "DIL" + "AVS" + "DLGE")
    insert = block + "ACTCATGCACGTTACCAGAAACGTCATGAACCGTACCAGCAT"
    h_nt = 12                     # H: codon 4 of the block
    j_nt = h_nt - 2               # J: frame +2, one residue earlier
    j_protein = translate(insert[j_nt:j_nt + 3 * 24])
    assert '*' not in j_protein, j_protein
    q_left, q_right = background(300), background(300)
    query = q_left + insert + q_right

    # Subject: the H site in frame +1 and, further down, J recoded between stops.
    h_site = encode("CIDL" + HD + "LGD" + "AVS" + "LDIC")
    j_site = "TAA" + encode(j_protein, pick=1) + "TAA"
    s_left, s_mid, s_right = background(300), background(150), background(300)
    subject = s_left + h_site + s_mid + j_site + s_right

    write(f"{outdir}/short_hsp.query.fna", "shorthsp_query", query)
    write(f"{outdir}/short_hsp.subject.fna", "shorthsp_subject", subject)


main()

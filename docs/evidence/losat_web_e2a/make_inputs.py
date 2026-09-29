#!/usr/bin/env python3
"""Regenerate the constructed outfmt 0 fixture inputs (LOSAT/tests/fasta/outfmt0/).

Every file is derived from fixed random seeds, so the output bytes are always the same;
evidence.sha256 records them. The manifest LOSAT/tests/outfmt0_manifest.tsv says which
fixture uses which file and what it covers.

Usage: make_inputs.py [OUTDIR]   (default: LOSAT/tests/fasta/outfmt0 of this repository)
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
COMPACT = REPO / "LOSAT/tests/fasta/blastn_parity_compact.fasta"
out = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "LOSAT/tests/fasta/outfmt0"
out.mkdir(parents=True, exist_ok=True)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def wrapped(defline, seq, width):
    return ">" + defline + "\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


# strand_*: a Plus/Minus and a Plus/Plus HSP on one subject, gaps in both rows.
random.seed(20260929)


def rnd_g(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate_r(s, rate, rng):
    return "".join(rng.choice([x for x in "ACGT" if x != c.upper()]) if rng.random() < rate else c for c in s)


rng = random.Random(7)
subj = rnd_g(1200)
seg_a, seg_b = subj[100:400], subj[600:900]
a = mutate_r(seg_a, 0.03, rng)
a = a[:120] + a[123:200] + "GT" + a[200:]
b = rc(mutate_r(seg_b, 0.02, rng))
query = rnd_g(40) + a + rnd_g(30) + b + rnd_g(25)
(out / "strand_subject.fasta").write_text(
    wrapped("subj_strand synthetic subject for plus/minus strand, gaps and multiple HSPs", subj, 70))
(out / "strand_query.fasta").write_text(
    wrapped("q_strand synthetic query with plus and minus strand segments", query, 60))

# mask_*: a CA repeat for DUST and lowercase runs for -lcase_masking.
rng = random.Random(424242)


def rnd(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


subj = rnd(220) + "CA" * 24 + rnd(232)
q = list(subj[50:450])
for p in (30, 260, 333):
    q[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[q[p]]
q = "".join(q)
q = q[:100] + q[100:120].lower() + q[120:350] + q[350:360].lower() + q[360:]
(out / "mask_subject.fasta").write_text(wrapped("subj_mask subject with a CA-repeat low-complexity block", subj, 80))
(out / "mask_query.fasta").write_text(wrapped("q_mask query with lowercase runs and a CA-repeat block", q, 80))
s = subj[:300] + subj[300:330].lower() + subj[330:]
(out / "mask_subject_lc.fasta").write_text(wrapped("subj_mask_lc subject with lowercase run at 301-330", s, 80))

# multi_*: 3 queries x 6 subjects.
rng = random.Random(99)


def mutate(s, rate):
    return "".join(rng.choice([x for x in "ACGT" if x != c]) if rng.random() < rate else c for c in s)


base = rnd(600)
other = rnd(400)
subs = [("msA", "close homolog of base (1% divergence)", rnd(50) + mutate(base, 0.01) + rnd(50)),
        ("msB", "medium homolog of base (6% divergence), reverse complemented",
         rc(rnd(30) + mutate(base, 0.06) + rnd(30))),
        ("msC", "distant homolog of base (12% divergence)", mutate(base[100:500], 0.12)),
        ("msD", "two copies of a base fragment, one on each strand",
         rnd(40) + mutate(base[0:200], 0.02) + rnd(60) + rc(mutate(base[300:500], 0.03)) + rnd(40)),
        ("msE", "unrelated random sequence", rnd(500)),
        ("msF", "contains the second query region", rnd(80) + mutate(other, 0.04) + rnd(80))]
(out / "multi_subject.fasta").write_text("".join(wrapped(f"{n} {d}", s, 70) for n, d, s in subs))
qs = [("mq1", "query matching several subjects", base),
      ("mq2", "query with no significant hit", rnd(300)),
      ("mq3", "query matching msF on the minus strand", rc(other[50:350]))]
(out / "multi_query.fasta").write_text("".join(wrapped(f"{n} {d}", s, 70) for n, d, s in qs))

# many_*: 260 subjects (more than the 250 default alignments).
rng = random.Random(2600)
base = rnd(200)
(out / "many_query.fasta").write_text(">many_q query shared by 260 subjects\n" + base + "\n")
with open(out / "many_subject.fasta", "w") as handle:
    for i in range(260):
        rate = 0.005 + 0.10 * i / 260
        handle.write(">ms%03d subject %d of 260 at divergence %.3f\n" % (i + 1, i + 1, rate))
        handle.write(rnd(25) + mutate(base, rate) + rnd(25) + "\n")

# width_*: an HSP ending at subject position 100 (the coordinate-width boundary).
rng = random.Random(100)
s = rnd(100)
(out / "width_subject.fasta").write_text(">width_s 100 bp subject\n" + s + "\n")
(out / "width_query.fasta").write_text(">width_q subject positions 40-100\n" + s[39:100] + "\n")

# longdef_*: deflines that wrap and truncate.
rng = random.Random(31337)
s = rnd(400)
qd = ("longdef_q1 Synthetic query, with commas, hyphen-separated-words and a very-long-token "
      "AAAAAAAAAABBBBBBBBBBCCCCCCCCCCDDDDDDDDDDEEEEEEEEEEFFFFFFFFFFGGGGGGGGGG end of the description "
      "line that keeps going beyond two wrapped lines of sixty-eight columns")
sd = ("longdef_s1 Synthetic subject title whose words keep going past sixty characters and then "
      "contains an unbreakable token XXXXXXXXXXYYYYYYYYYYZZZZZZZZZZXXXXXXXXXXYYYYYYYYYYZZZZZZZZZZ "
      "followed by more words to wrap again at the next whitespace")
(out / "longdef_query.fasta").write_text(">" + qd + "\n" + s[20:380] + "\n")
(out / "longdef_subject.fasta").write_text(">" + sd + "\n" + s + "\n>longdef_s2 short\n" + s[::-1] + "\n")

# edge_*: an invalid (all N) query, and a fully lowercase query next to a valid one.
(out / "edge_allN.fasta").write_text(">allN query of only N\n" + "N" * 40 + "\n")
alpha = "".join(line.strip() for line in COMPACT.read_text().splitlines() if not line.startswith(">"))[:76]
(out / "edge_alllower.fasta").write_text(
    ">alllower fully lowercase copy of alpha\n" + alpha.lower() + "\n>ok_alpha uppercase alpha\n" + alpha + "\n")

# pwidth_*: the same coordinate-width boundary for BLASTP (subject positions 40-100).
rng = random.Random(5)
s = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(100))
(out / "pwidth_s.faa").write_text(">pw_s\n" + s + "\n")
(out / "pwidth_q.faa").write_text(">pw_q\n" + s[39:100] + "\n")

# twidth_*: the same boundary for TBLASTN (frame +2 ends at subject nucleotide 100).
rng = random.Random(11)
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codons: dict[str, list[str]] = {}
for index, amino in enumerate(CODE):
    codons.setdefault(amino, []).append("TCAG"[index // 16] + "TCAG"[index // 4 % 4] + "TCAG"[index % 4])
prot = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(33))
nt = "A" + "".join(rng.choice(codons[amino]) for amino in prot)
(out / "twidth_s.fna").write_text(">tw_s 100 bp subject, frame +2 ends at 100\n" + nt + "\n")
(out / "twidth_q.faa").write_text(">tw_q\n" + prot + "\n")

# edge_mixed_allN / edge_ambig (added in S07): an all-N query between two valid ones, and
# a query with N runs and other IUPAC codes (the ungapped Karlin block depends on the
# query composition).
records = [chunk.split("\n", 1) for chunk in COMPACT.read_text().split(">")[1:]]
alpha_def, alpha_seq = records[0][0], records[0][1].replace("\n", "")
beta_def, beta_seq = records[1][0], records[1][1].replace("\n", "")
(out / "edge_mixed_allN.fasta").write_text(
    ">" + alpha_def + "\n" + alpha_seq + "\n>allN query of only N\n" + "N" * 40 + "\n>" + beta_def + "\n" + beta_seq + "\n")
ambig = list(alpha_seq)
ambig[10:15] = "NNNNN"
ambig[30], ambig[50], ambig[60] = "R", "Y", "K"
(out / "edge_ambig.fasta").write_text(">ambig alpha with N and IUPAC codes\n" + "".join(ambig) + "\n")

# edge_batch_allN (added in S07, after the independent audit): an all-N query that fills
# NCBI's first query batch (5000 residues for these subjects) on its own, then a valid query.
(out / "edge_batch_allN.fasta").write_text(
    ">batchN 5000 N filling the first query batch\n" + "N" * 5000 + "\n>" + alpha_def + "\n" + alpha_seq + "\n")

# lcase_minus (added in S07, after the independent audit): a query with a lowercase run and
# no -lcase_masking, whose only hit is on the minus strand of the subject.
rng = random.Random(4096)
plain = rnd(200)
(out / "lcase_minus_query.fasta").write_text(
    ">lcase_minus query with lowercase 91-110\n" + plain[:90] + plain[90:110].lower() + plain[110:] + "\n")
(out / "lcase_minus_subject.fasta").write_text(">lcase_minus_s reverse complement of the query\n" + rc(plain) + "\n")

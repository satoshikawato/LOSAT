#!/usr/bin/env python3
"""Regenerate the TBLASTX outfmt 0/7 fixture inputs (LOSAT/tests/fasta/outfmt0/tblastx_*).

Every file is cut from genomes of LOSAT/tests/tblastx_v010_parity_manifest.tsv (or built
from fixed letters), so the output bytes are always the same; evidence.sha256 records them.
The manifest LOSAT/tests/outfmt0_manifest.tsv says which fixture uses which file and what
it covers (the `tblastx.*` rows).

Usage: make_inputs.py [OUTDIR]   (default: LOSAT/tests/fasta/outfmt0 of this repository)
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
FASTA = REPO / "LOSAT/tests/fasta"
out = Path(sys.argv[1]) if len(sys.argv) > 1 else FASTA / "outfmt0"
out.mkdir(parents=True, exist_ok=True)


def genome(name: str) -> tuple[str, str]:
    lines = (FASTA / f"{name}.fasta").read_text().splitlines()
    return lines[0][1:], "".join(lines[1:])


def record(defline: str, seq: str, width: int = 70) -> str:
    return ">" + defline + "\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


def window(name: str, title: str, seq: str, start: int, end: int) -> str:
    """A 0-based [start, end) window, named by its 1-based coordinates."""
    return record(f"{name}_{start + 1}_{end} window {start + 1}-{end} of {title}", seq[start:end])


q_title, q_seq = genome("LC738874")      # MelaMJNV (p03, p12 of the parity manifest)
s_title, s_seq = genome("LC738875")      # p12 subject
n_title, n_seq = genome("AvCLPV")        # p11: unrelated to the windows below

# multi: four queries (one without hits) against four subjects, 6 frames, many HSPs.
(out / "tblastx_multi_query.fasta").write_text(
    window("LC738874", q_title, q_seq, 100000, 110000)
    + window("AvCLPV", n_title, n_seq, 120000, 120300)
    + window("LC738874", q_title, q_seq, 120000, 130000)
    + window("LC738874", q_title, q_seq, 30000, 40000))
(out / "tblastx_multi_subject.fasta").write_text(
    window("LC738875", s_title, s_seq, 120000, 130000)
    + window("LC738875", s_title, s_seq, 300000, 310000)
    + window("LC738875", s_title, s_seq, 270000, 280000)
    + window("LC738875", s_title, s_seq, 130000, 140000))

# batch: queries of unequal lengths across NCBI's 10002-nt query batches (the linking
# cutoffs use each batch's average query length and smallest Lambda).
(out / "tblastx_batch_query.fasta").write_text(
    window("LC738874", q_title, q_seq, 100000, 101500)
    + window("LC738874", q_title, q_seq, 120000, 140000)
    + window("LC738874", q_title, q_seq, 30000, 30800)
    + window("LC738874", q_title, q_seq, 180000, 186000))

# invalid: an all-N query and a 2-nt query in a batch with valid queries (searched;
# no warning, no Karlin block, search space 0).
(out / "tblastx_invalid_query.fasta").write_text(
    window("LC738874", q_title, q_seq, 236000, 239000)
    + record("allN 300 N", "N" * 300)
    + record("short2 two bases", "AC")
    + window("LC738874", q_title, q_seq, 120000, 123000))

# unsearched: a batch of only invalid queries after a full batch (not searched: warnings,
# -1 Karlin blocks, no "hits found" line in outfmt 7).
(out / "tblastx_unsearched_query.fasta").write_text(
    window("LC738874", q_title, q_seq, 230000, 240003)
    + record("allN 300 N", "N" * 300)
    + record("short2 two bases", "AC"))

# ambig: IUPAC codons in the displayed rows (B, Z, J, X; display and search translations
# differ), an N island and ambiguity letters in the subject (NCBI's random ncbi2na bases
# in the preliminary search), and a lowercase run in the query (not shown).
query = list(q_seq[234000:242000])
# Frame +2 HSP of the full pair (Query 235430-236149): codon i starts at 1429 + 3i.
for i, codon in zip((10, 20, 30, 40, 50, 60, 70), ("RAY", "SAR", "MTA", "GCN", "NNN", "TAR", "HTA")):
    query[1429 + 3 * i:1432 + 3 * i] = list(codon)
query[3000:3100] = [c.lower() for c in query[3000:3100]]
subject = list(s_seq[264000:270000])
subject[4000:4030] = list("N" * 30)
for offset, letter in zip(range(2000, 2200, 20), "RYKMSWBDHV"):
    subject[offset] = letter
(out / "tblastx_ambig_query.fasta").write_text(
    record(f"LC738874_234001_242000_iupac IUPAC codons in window 234001-242000 of {q_title}", "".join(query)))
(out / "tblastx_ambig_subject.fasta").write_text(
    record(f"LC738875_264001_270000_iupac N island and IUPAC letters in window 264001-270000 of {s_title}", "".join(subject)))

# code4: windows of the genetic code 4 genome pair (d04/p14/d06 of the parity manifest).
c4q_title, c4q_seq = genome("AP027131")
c4s_title, c4s_seq = genome("AP027133")
(out / "tblastx_code4_query.fasta").write_text(window("AP027131", c4q_title, c4q_seq, 100000, 115000))
(out / "tblastx_code4_subject.fasta").write_text(window("AP027133", c4s_title, c4s_seq, 95000, 115000))

# many: 260 subjects that each encode the query's peptide in frame +1 with other
# synonymous codons (and a few substitutions), so that most of them have one HSP: 260
# description rows but 250 alignments without -max_target_seqs.
import random  # noqa: E402

rng = random.Random(20261002)
CODONS: dict[str, list[str]] = {}
for a in "TCAG":
    for b in "TCAG":
        for c in "TCAG":
            aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"["TCAG".index(a) * 16 + "TCAG".index(b) * 4 + "TCAG".index(c)]
            CODONS.setdefault(aa, []).append(a + b + c)
residues = [aa for aa in CODONS if aa != "*"]


def rnd(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


peptide = "".join(rng.choice(residues) for _ in range(40))
query_seq = rnd(30) + "".join(rng.choice(CODONS[aa]) for aa in peptide) + rnd(30)
(out / "tblastx_many_query.fasta").write_text(record("pep40 a 40-residue peptide in frame +1 with random flanks", query_seq))
subjects = []
for index in range(260):
    variant = list(peptide)
    for _ in range(index % 7):
        variant[rng.randrange(len(variant))] = rng.choice(residues)
    coding = "".join(rng.choice(CODONS[aa]) for aa in variant)
    subjects.append(record(f"subject{index:03d} the peptide with {index % 7} substitutions", rnd(30) + coding + rnd(30)))
(out / "tblastx_many_subject.fasta").write_text("".join(subjects))

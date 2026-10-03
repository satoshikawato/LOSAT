#!/usr/bin/env python3
"""Write the E2i (Session SD) fixture inputs: dc-megablast and blastn-short.

Deterministic (slices of the repository genomes and fixed edits); the files are committed
under LOSAT/tests/fasta/outfmt0/, so later runs never regenerate them.

- dc_query.fasta: two windows of LC738874 (MJNV of Melicertus latisulcatus), with a
  lowercase run (for -lcase_masking), a run of N and IUPAC codes in the first one.
- dc_subject.fasta: three windows of LC738875 (PMNV, a distant virus) that the query
  windows hit on the plus strand (dcs1), on the minus strand (dcs2, dcs3), and a random
  record without hits (dcs4).
- dc_small_query.fasta: a 3 kb window (the diagonal array: the query block is at most 8000).
- short_subject.fasta: the first 4000 bases of LC738875 and LC738874, and a window of
  LC738875 with a CA repeat (for -dust).
- short_query.fasta: primer-length queries cut from the subjects (exact, reverse
  complement, mismatches, IUPAC codes, lowercase, a long one, one across the CA repeat)
  and a random one.
- short_nohit_query.fasta: one random primer.

Usage: make_inputs.py  (writes into LOSAT/tests/fasta/outfmt0/)
"""
from __future__ import annotations

import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
FASTA = REPO / "LOSAT/tests/fasta"
OUT = FASTA / "outfmt0"
COMPLEMENT = str.maketrans("ACGTacgtRYKMrykmNn", "TGCAtgcaYRMKyrmkNn")


def read(name: str) -> str:
    lines = (FASTA / name).read_text().splitlines()
    return "".join(line.strip() for line in lines[1:] if not line.startswith(">")).upper()


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def window(seq: str, start: int, end: int) -> str:
    """1-based inclusive window."""
    return seq[start - 1:end]


def write(name: str, records: list[tuple[str, str]]) -> None:
    with open(OUT / name, "w") as handle:
        for title, seq in records:
            handle.write(f">{title}\n")
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")


def edit(seq: str, at: int, text: str) -> str:
    """Replace the letters at 0-based offset `at` with `text`."""
    return seq[:at] + text + seq[at + len(text):]


def main() -> None:
    rng = random.Random(20261003)
    lc74 = read("LC738874.fasta")
    lc75 = read("LC738875.fasta")

    q1 = window(lc74, 178001, 190000)
    q1 = q1[:4000] + q1[4000:4300].lower() + q1[4300:]
    q1 = edit(q1, 7000, "N" * 25)
    q1 = edit(q1, 9500, "R")
    q1 = edit(q1, 9520, "Y")
    q1 = edit(q1, 10500, "K")
    q2 = window(lc74, 232001, 241000)
    write("dc_query.fasta", [
        ("dcq1 LC738874 178001-190000 with a lowercase run, N run and IUPAC codes", q1),
        ("dcq2 LC738874 232001-241000", q2),
    ])
    random_record = "".join(rng.choice("ACGT") for _ in range(2000))
    write("dc_subject.fasta", [
        ("dcs1 LC738875 240001-246500", window(lc75, 240001, 246500)),
        ("dcs2 LC738875 264001-271000", window(lc75, 264001, 271000)),
        ("dcs3 LC738875 251000-249001 reverse complement", revcomp(window(lc75, 249001, 251000))),
        ("dcs4 random 2000", random_record),
    ])
    write("dc_small_query.fasta", [("dcsq LC738874 179001-182000", window(lc74, 179001, 182000))])

    s1 = window(lc75, 1, 4000)
    s2 = window(lc74, 1, 4000)
    # A CA repeat between two windows: without DUST (the blastn-short default) the primer
    # p8 across it hits; -dust yes masks the repeat.
    s3 = window(lc75, 5001, 5300) + "CA" * 30 + window(lc75, 5301, 5600)
    write("short_subject.fasta", [
        ("shs1 LC738875 1-4000", s1),
        ("shs2 LC738874 1-4000", s2),
        ("shs3 LC738875 5001-5300, (CA)30, LC738875 5301-5600", s3),
    ])
    p3 = window(s1, 3001, 3028)
    p3 = edit(p3, 9, "A" if p3[9] != "A" else "C")
    p3 = edit(p3, 18, "G" if p3[18] != "G" else "T")
    p4 = window(s1, 1501, 1522)
    p4 = edit(p4, 5, "N")
    p4 = edit(p4, 14, "R")
    p5 = window(s2, 2001, 2035)
    p5 = p5[:12] + p5[12:22].lower() + p5[22:]
    write("short_query.fasta", [
        ("p1 shs1 1001-1025", window(s1, 1001, 1025)),
        ("p2 shs1 2030-2001 reverse complement", revcomp(window(s1, 2001, 2030))),
        ("p3 shs1 3001-3028 with two mismatches", p3),
        ("p4 shs1 1501-1522 with N and R", p4),
        ("p5 shs2 2001-2035 with a lowercase run", p5),
        ("p6 random 20", "".join(rng.choice("ACGT") for _ in range(20))),
        ("p7 shs2 3045-3001 reverse complement", revcomp(window(s2, 3001, 3045))),
        ("p8 shs3 291-338 across the CA repeat", window(s3, 291, 338)),
    ])
    write("short_nohit_query.fasta", [("nohit random 24", "".join(rng.choice("ACGT") for _ in range(24)))])


if __name__ == "__main__":
    main()

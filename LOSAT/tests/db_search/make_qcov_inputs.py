#!/usr/bin/env python3
"""Writes the inputs of the -qcov_hsp_perc fixtures of the ungapped searches in tests/db_search: qcovp_q.fna (the
query), qcovp.fna (the TBLASTX subjects) and qcovp.faa (the BLASTX subjects, records of the database `qcovp`).

Run from LOSAT/: python3 tests/db_search/make_qcov_inputs.py. The query is 600 random bases. In each subject
file the first record (`s_part`, `p_part`) holds the query's first 330 bases unchanged (55 % of the query, the
best e-value) and the second (`s_full`, `p_full`) the whole query with 25 % substitutions (a weaker e-value).
With -max_target_seqs 1 -qcov_hsp_perc 60 NCBI BLAST+ 2.17.0 reports the second: an ungapped search removes the
HSPs below the coverage before the preliminary hit list (of the final size) is filled (blast_engine.c:1538).
Random letters come from random.Random(20261019).
"""
import random
from pathlib import Path

rng = random.Random(20261019)
OUT = Path("tests/db_search")
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def rnd(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq: str, rate: float) -> str:
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in seq)


def translate(seq: str) -> str:
    return "".join(CODONS[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3))


def wrap(seq: str) -> str:
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


query = rnd(600)
part = rnd(300) + query[:330] + rnd(300)
full = rnd(200) + mutate(query, 0.25) + rnd(200)
(OUT / "qcovp_q.fna").write_text(f">qcovp_q\n{wrap(query)}\n")
(OUT / "qcovp.fna").write_text(f">s_part exact 330 nt of the query\n{wrap(part)}\n"
                               f">s_full the query with 0.25 substitutions\n{wrap(full)}\n")
(OUT / "qcovp.faa").write_text(f">p_part translation of 330 nt of the query\n{wrap(translate(query[:330]))}\n"
                               f">p_full translation of the query with 0.25 substitutions\n"
                               f"{wrap(translate(mutate(query, 0.25)))}\n")

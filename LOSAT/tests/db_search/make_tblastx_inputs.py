#!/usr/bin/env python3
"""Writes the inputs of the TBLASTX -db fixtures (rows `db.tblastx.*`) in tests/db_search: xslice.fna (database
records of `xslice`) and the query files xslice_q.fna, xslice_nohit.fna and
xslice_nohit_masked.fna. TBLASTX is slow, so the database holds
about 14 kb and the queries 5 kb.

Run from LOSAT/: python3 tests/db_search/make_tblastx_inputs.py. Slices of tests/fasta/MelaMJNV.fasta (M) and
LvMJNV.fasta (V), two shrimp viruses; the queries are slices of tests/fasta/PemoMJNVA.fasta (P, record LC738870.1;
its genes BDT61858.1 (plus, 10690-15105) and BDT61866.1 (minus, 31307-33499) are 71 % and 61 % identical to the
MelaMJNV genes at 38595-41240 and 59399-61579, found with NCBI tblastn -subject) and of LvMJNV. Mutations and
ambiguity codes come from random.Random(20261018). Slice indices are 0-based, the coordinates in the comments are
1-based genome positions.
"""
import random
from pathlib import Path

rng = random.Random(20261018)
IUPAC = "RYKMSWBDHVN"


def read(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def mutate(seq: str, rate: float, indels: int = 0) -> str:
    out = []
    for base in seq:
        out.append(rng.choice([b for b in "ACGT" if b != base]) if rng.random() < rate else base)
    for _ in range(indels):
        at = rng.randrange(100, len(out) - 100)
        if rng.random() < 0.5:
            del out[at:at + 3 * rng.randrange(1, 4)]  # whole codons, so the reading frame survives
        else:
            out[at:at] = [rng.choice("ACGT") for _ in range(3 * rng.randrange(1, 4))]
    return "".join(out)


def put(seq: str, at: int, text: str) -> str:
    return seq[:at] + text + seq[at + len(text):]


def scatter(seq: str, start: int, end: int, step: int) -> str:
    """One random ambiguity code every `step` bases in [start, end)."""
    for at in range(start, end, step):
        seq = put(seq, at, rng.choice(IUPAC))
    return seq


def fasta(records: list[tuple[str, str]], width: int) -> str:
    return "".join(f">{name}\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width)) for name, seq in records)


M, V, P = (read(f"tests/fasta/{n}.fasta") for n in ("MelaMJNV", "LvMJNV", "PemoMJNVA"))

# m_plus: MelaMJNV 37001-41500, the plus-strand gene 38595-41240 (BDT61858 of the Mela annotation) and its start.
m_plus = M[37000:41500]
m_plus = put(m_plus, 39000 - 37001, "ACGRYTTGGMKAWSC")        # inside the gene
m_plus = put(m_plus, 39800 - 37001, "N" * 30)                # an N block inside the gene
m_plus = scatter(m_plus, 38700 - 37001, 41200 - 37001, 97)   # single codes inside the gene
# m_minus: MelaMJNV 59301-61700, the minus-strand gene 59399-61579.
m_minus = M[59300:61700]
m_minus = put(m_minus, 60500 - 59301, "RYKMSWBDHVN")
m_minus = scatter(m_minus, 59500 - 59301, 61500 - 59301, 151)
# v_mut / v_far: LvMJNV 28501-30700 with 8 % and 22 % substitutions (the first with whole-codon indels).
v_mut = mutate(V[28500:30700], 0.08, indels=4)
v_far = mutate(V[28500:30700], 0.22)
# v_slice: LvMJNV 37001-39000, genes on both strands, with ambiguity codes.
v_slice = V[37000:39000]
v_slice = put(v_slice, 37900 - 37001, "WSKMRYBDHV")
v_slice = put(v_slice, 38300 - 37001, "N" * 12)
v_slice = scatter(v_slice, 37300 - 37001, 38700 - 37001, 83)

write = Path("tests/db_search")
(write / "xslice.fna").write_text(fasta([
    ("m_plus MelaMJNV 37001-41500 with ambiguity codes and an N block", m_plus),
    ("m_minus MelaMJNV 59301-61700 with ambiguity codes", m_minus),
    ("v_mut_8pct LvMJNV 28501-30700 with substitutions and indels", v_mut),
    ("v_far_22pct LvMJNV 28501-30700 diverged", v_far),
    ("v_slice LvMJNV 37001-39000 with ambiguity codes", v_slice),
    ("tiny shorter than 10 bases", "ACGTACGT"),
], 70))

(write / "xslice_q.fna").write_text(fasta([
    ("q_p_plus PemoMJNVA 10601-12200", P[10600:12200]),
    ("q_p_minus PemoMJNVA 32001-33500", P[32000:33500]),
    ("q_v_a LvMJNV 28501-29500", V[28500:29500]),
    ("q_v_b LvMJNV 37801-38600", V[37800:38600]),
], 70))
# A random query of 45 bases with no hit in xslice at the default E value. Every random query of 60 or more bases has
# weak hits (E value up to 10) in this database; of 20 candidates of 24, 30, 36 and 45 bases from
# random.Random(20261019), found with NCBI tblastx, number 3 (45 bases) has none.
rng = random.Random(20261019)
candidates = ["".join(rng.choice("ACGT") for _ in range(length)) for length in [24, 30, 36, 45] * 5][:20]
(write / "xslice_nohit.fna").write_text(fasta([("q_nohit random 45 bases", candidates[3])], 70))
# A query that SEG masks completely (CAG repeats translate to glutamine, serine and alanine tracts).
(write / "xslice_nohit_masked.fna").write_text(fasta([("q_masked 300 bases of CAG repeats", "CAG" * 100)], 70))

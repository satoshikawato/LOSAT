#!/usr/bin/env python3
"""Writes the inputs of the TBLASTN -db fixtures (rows `db.tblastn.*`) in tests/db_search: tslice.fna (database
records of `tslice`), tgc4.fna (records of `tgc4`) and the query files tslice_q.faa, tslice_q_lc.faa, tgc4_q.faa,
tslice_nohit.faa, tslice_nohit_masked.faa.

Run from LOSAT/: python3 tests/db_search/make_tblastn_inputs.py. Slices of tests/fasta/MelaMJNV.fasta (M) and
LvMJNV.fasta (V), two shrimp viruses, and of AP027133.fasta (H, Candidatus Hepatoplasma scabrum, genetic code 4);
queries are proteins of tests/fasta/PemoMJNVA.faa, LvMJNV.faa and AP027131.faa (found by NCBI tblastn -subject on
the whole genomes). Mutations and ambiguity codes come from random.Random(20261016). The slice borders were
chosen from those hits: slice indices are 0-based, the coordinates in the comments are 1-based genome positions.
"""
import random
from pathlib import Path

rng = random.Random(20261016)
IUPAC = "RYKMSWBDHVN"


def read(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def proteins(path: str) -> dict[str, tuple[str, str]]:
    out, name, title = {}, None, ""
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name, _, title = line[1:].partition(" ")
            out[name] = (title, "")
        elif name:
            out[name] = (out[name][0], out[name][1] + line.strip())
    return out


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


M, V, H = (read(f"tests/fasta/{n}.fasta") for n in ("MelaMJNV", "LvMJNV", "AP027133"))

# m_plus: MelaMJNV 36001-44000, plus-strand genes BDT61857 (37640-38329), BDT61858 (38595-41240), BDT61859 (43175-43576).
m_plus = M[36000:44000]
m_plus = put(m_plus, 39500 - 36001, "ACGRYTTGGMKAWSC")       # inside BDT61858
m_plus = put(m_plus, 40200 - 36001, "N" * 30)               # an N block inside BDT61858
m_plus = scatter(m_plus, 37700 - 36001, 38300 - 36001, 97)  # single codes inside BDT61857
# m_minus: MelaMJNV 59001-62000, minus-strand gene BDT61866 (59399-61579).
m_minus = M[59000:62000]
m_minus = put(m_minus, 60500 - 59001, "RYKMSWBDHVN")
m_minus = scatter(m_minus, 59600 - 59001, 61400 - 59001, 151)
# v_slice: LvMJNV 21001-39000, genes BDT62090 (minus 21385-23397), BDT62091 (plus 24203-24766), BDT62093 (minus
# 30761-31846), BDT62095 (plus 37239-38051), BDT62096 (minus 38073-38708).
v_slice = V[21000:39000]
v_slice = put(v_slice, 22000 - 21001, "WSKMRYBDHV")
v_slice = put(v_slice, 31000 - 21001, "N" * 12)
v_slice = scatter(v_slice, 37300 - 21001, 38000 - 21001, 83)
# v_mut: BDT62092 (plus 28709-30493) region of LvMJNV with 8 % substitutions and whole-codon indels.
v_mut = mutate(V[28500:30700], 0.08, indels=4)
# v_mut_far: the same region with 22 % substitutions, a weak homolog.
v_far = mutate(V[28500:30700], 0.22)

write = Path("tests/db_search")
(write / "tslice.fna").write_text(fasta([
    ("m_plus MelaMJNV 36001-44000 with ambiguity codes and an N block", m_plus),
    ("m_minus MelaMJNV 59001-62000 with ambiguity codes", m_minus),
    ("v_slice LvMJNV 21001-39000 with ambiguity codes", v_slice),
    ("v_mut_8pct LvMJNV 28501-30700 with substitutions and indels", v_mut),
    ("v_far_22pct LvMJNV 28501-30700 diverged", v_far),
], 70))

# tgc4: AP027133 106001-121000 (genetic code 4: TGA is Trp), a mutated copy and a sequence of 8 bases.
h_slice = H[106000:121000]
h_slice = put(h_slice, 119500 - 106001, "RYMKSW")
h_slice = put(h_slice, 112800 - 106001, "N" * 21)
h_slice = scatter(h_slice, 107700 - 106001, 108800 - 106001, 71)
h_mut = mutate(H[118900:120300], 0.1, indels=3)
(write / "tgc4.fna").write_text(fasta([
    ("h_slice AP027133 106001-121000 with ambiguity codes and an N block", h_slice),
    ("h_mut_10pct AP027133 118901-120300 with substitutions and indels", h_mut),
    ("tiny shorter than 10 bases", "ACGTACGT"),
], 70))


def protein_fasta(path: str, names: list[str]) -> str:
    table = proteins(path)
    return fasta([(f"{n} {table[n][0]}", table[n][1]) for n in names], 60)


(write / "tslice_q.faa").write_text(
    protein_fasta("tests/fasta/PemoMJNVA.faa", ["BDT61857.1", "BDT61859.1", "BDT61866.1"])
    + protein_fasta("tests/fasta/LvMJNV.faa", ["BDT62090.1", "BDT62091.1", "BDT62092.1", "BDT62096.1"]))
q95 = proteins("tests/fasta/LvMJNV.faa")["BDT62095.1"]
lc = q95[1]
lc = lc[:40] + lc[40:100].lower() + lc[100:190] + lc[190:230].lower() + lc[230:]
(write / "tslice_q_lc.faa").write_text(fasta([(f"BDT62095.1 {q95[0]}, stretches 41-100 and 191-230 in lower case", lc)], 60))
(write / "tgc4_q.faa").write_text(
    protein_fasta("tests/fasta/AP027131.faa", ["BDV02158.1", "BDV02159.1", "BDV02443.1", "BDV02702.1", "BDV02441.1"]))
low = "GSGSSGGSGSGGSSGSGGSGSSGGSGSGSSGSGGSGSGGSSGSGGSGSSGGSGSGSSGSGGSGSGGSSGSGG"
(write / "tslice_nohit_masked.faa").write_text(fasta([("q_masked low-complexity glycine/serine sequence of 72 residues", low)], 60))
# A random protein with no hit in tslice or tgc4 at the default E value (found with NCBI tblastn: of 20 random
# candidates of 40, 60 and 80 residues, from random.Random(20261017), only number 3 has none in either database).
rng = random.Random(20261017)
candidates = ["".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(length)) for length in [40, 60, 80] * 7][:20]
(write / "tslice_nohit.faa").write_text(fasta([("q_nohit random protein of 40 residues", candidates[3])], 60))

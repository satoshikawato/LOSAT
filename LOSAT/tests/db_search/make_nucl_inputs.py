#!/usr/bin/env python3
"""Writes tests/db_search/nslice.fna (database records) and nslice_q.fna (queries) of the BLASTN -db fixtures.

Run from LOSAT/: python3 tests/db_search/make_nucl_inputs.py. Slices of tests/fasta/LC738875.fasta (L) and
tests/fasta/AP027202.fasta (P), two related viruses; mutations come from random.Random(20261015).
"""
import random
from pathlib import Path

rng = random.Random(20261015)


def read(name: str) -> str:
    return "".join(l.strip() for l in Path(f"tests/fasta/{name}.fasta").read_text().splitlines() if not l.startswith(">"))


def mutate(seq: str, rate: float, indels: int = 0) -> str:
    out = []
    for base in seq:
        out.append(rng.choice([b for b in "ACGT" if b != base]) if rng.random() < rate else base)
    for _ in range(indels):
        at = rng.randrange(100, len(out) - 100)
        if rng.random() < 0.5:
            del out[at:at + rng.randrange(1, 9)]
        else:
            out[at:at] = [rng.choice("ACGT") for _ in range(rng.randrange(1, 9))]
    return "".join(out)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN"))[::-1]


def fasta(records: list[tuple[str, str]]) -> str:
    return "".join(f">{name}\n" + "".join(seq[i:i + 70] + "\n" for i in range(0, len(seq), 70)) for name, seq in records)


def put(seq: str, at: int, text: str) -> str:
    return seq[:at] + text + seq[at + len(text):]


L, P = read("LC738875"), read("AP027202")
a = L[0:6000]
a_mut = mutate(a, 0.08, indels=6)
b = P[0:6000]
c = L[20000:23000]
c = put(c, 300, "RYKMSWBDHVN")
c = put(c, 900, "N" * 60)
c = put(c, 1500, "ACGTNACGRTACGYT")
c = put(c, 2000, "WWWWWSSSSS")
d = mutate(L[40000:44000], 0.15)
e = P[30000:32500]
write = Path("tests/db_search")
(write / "nslice.fna").write_text(fasta([
    ("a_l6k exact slice of LC738875", a),
    ("a_mut_8pct substitutions and indels", a_mut),
    ("b_p6k slice of AP027202", b),
    ("c_amb ambiguity codes and an N block", c),
    ("d_mut_15pct", d),
    ("e_p2k slice of AP027202", e),
]))
rand = "".join(rng.choice("ACGT") for _ in range(900))
(write / "nslice_q.fna").write_text(fasta([
    ("q_a 1500 of a", a[2000:3500]),
    ("q_c_amb 1400 across the ambiguity codes, 4 percent substitutions", mutate(c[100:1500], 0.04)),
    ("q_d_minus reverse complement of d, 12 percent substitutions", revcomp(mutate(d[500:2500], 0.12))),
    ("q_e_p diverged", mutate(e[200:1800], 0.2)),
    ("q_nohit random", rand),
]))
(write / "nslice_nohit.fna").write_text(fasta([("nohit random", rand)]))

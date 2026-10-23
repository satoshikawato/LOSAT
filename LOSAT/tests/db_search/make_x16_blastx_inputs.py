#!/usr/bin/env python3
"""Writes the inputs of the BLASTX outfmt 15 and 16 fixtures (OF5b part 1): tests/fasta/outfmt0/x16_bx_q_unsearched.fna
and tests/db_search/x16_bx_taxid_q.fna.

Run from LOSAT/: python3 tests/db_search/make_x16_blastx_inputs.py.

x16_bx_q_unsearched.fna: a 10,003-base slice of tests/fasta/MelaMJNV.fasta (26001-36003; two lower-case stretches, 301-700
and 1501-1800, for -lcase_masking) closes NCBI's first query batch (10,002 bases, `BATCH_SIZE`), and the next two records,
300 N and two bases, form a second batch with no valid context: NCBI does not search it, so those two queries have no `hits`
member in outfmt 15 and 16 (the same shape as tests/fasta/outfmt0/tblastx_unsearched_query.fasta).

x16_bx_taxid_q2.fna: two queries of 10,003 bases (each closes an NCBI query batch of 10,002 bases) that hit the proteins of
x16_taxid in a minus and a plus frame, so that the taxonomy database is asked in two batches.

x16_bx_taxid_q.fna: three nucleotide queries for the protein database x16_taxid (tests/db_search/x16_taxid.faa): `bq1` the
reverse complement of a back-translation of qp1 (a minus frame), `bq2` the same protein from its second residue with one
extra base in front (frame +2), and `bq3`, 300 random bases without hits. Back-translation: a codon drawn with
random.Random(20261022) from the standard code.
"""
import random
from pathlib import Path


def read(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def fasta(records: list[tuple[str, str]], width: int) -> str:
    return "".join(f">{name}\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width)) for name, seq in records)


M = read("tests/fasta/MelaMJNV.fasta")
big = M[26000:36003]
assert len(big) == 10003
big = big[:300] + big[300:700].lower() + big[700:1500] + big[1500:1800].lower() + big[1800:]
Path("tests/fasta/outfmt0/x16_bx_q_unsearched.fna").write_text(fasta([
    ("m_big MelaMJNV 26001-36003, 301-700 and 1501-1800 in lower case", big),
    ("q_alln 300 N", "N" * 300),
    ("q_short2 two bases", "AC"),
], 70))

CODE = {}
bases = "TCAG"
aminos = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(aminos):
    CODE.setdefault(a, []).append(bases[i // 16] + bases[(i // 4) % 4] + bases[i % 4])
rng = random.Random(20261022)


def back(protein: str) -> str:
    return "".join(rng.choice(CODE[a]) for a in protein)


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


qp1 = read("tests/db_search/x16_taxid_q.faa").split("IRRNWKMVC")[0]
random_query = "".join(rng.choice("ACGT") for _ in range(300))
Path("tests/db_search/x16_bx_taxid_q.fna").write_text(fasta([
    ("bq1 reverse complement of a back-translation of qp1", rc(back(qp1))),
    ("bq2 back-translation of qp1 from residue 2, one base in front", "G" + back(qp1[1:])),
    ("bq3 300 random bases", random_query),
], 70))

# tests/fixtures/outfmt16 (tests/outfmt16_blastx.rs): queries for the proteins of sprot.faa. BLASTX reads its queries in
# batches of 10,002 bases (it does not read the BATCH_SIZE variable), so each valid query is padded with random bases to
# 10,003: `qnx_bad.fna` has two batches that are searched and a third whose line is not a nucleotide sequence
# (the failed-batch test); `qnx_firstbad.fna` is one batch with the bad line.
fixtures = Path("tests/fixtures/outfmt16")
sp1 = "".join(l.strip() for l in (fixtures / "sprot.faa").read_text().split(">")[1].splitlines()[1:])


def padded(core: str) -> str:
    return core + "".join(rng.choice("ACGT") for _ in range(10003 - len(core)))


q1 = padded(back(sp1))
q2 = padded(rc(back(sp1[3:150])))
(fixtures / "qnx_bad.fna").write_text(fasta([("q1 plus frame, padded to 10003", q1), ("q2 minus frame, padded to 10003", q2)], 70) + ">q3\nACGT%%%%\n")
(fixtures / "qnx_firstbad.fna").write_text(">q1 short\n" + back(sp1[:40]) + "\n>q2\nACGT%%%%\n")
# A batch with no hits first, then one with hits (the held-output test: a refusal in the second batch leaves nothing).
(fixtures / "qnx_late.fna").write_text(fasta([("q1 no hits, 10003 bases of a repeat", ("ACGTTGCAAGCTTGCA" * 700)[:10003]), ("q2 plus frame, padded to 10003", padded(back(sp1)))], 70))

# Two batches that both hit x16_taxid: the taxonomy database is opened and the volume asked in the first, the second asks for
# each Hit only (tests/db_search_manifest.tsv: db.x16_bx_taxid.batches.*).
def pad(core: str) -> str:
    return core + "".join(rng.choice("ACGT") for _ in range(10003 - len(core)))


Path("tests/db_search/x16_bx_taxid_q2.fna").write_text(fasta([
    ("bq1b reverse complement of a back-translation of qp1, padded to 10003", pad(rc(back(qp1)))),
    ("bq2b back-translation of qp1 from residue 2, padded to 10003", pad("G" + back(qp1[1:]))),
], 70))

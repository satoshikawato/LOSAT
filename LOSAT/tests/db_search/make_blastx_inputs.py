#!/usr/bin/env python3
"""Writes the inputs of the BLASTX -db fixtures (rows `db.blastx.*`) in tests/db_search: bslice.faa (records of the
protein database `bslice`) and the nucleotide query files bslice_q.fna, bslice_q_lc.fna, bslice_nohit.fna and
bslice_alln.fna.

Run from LOSAT/: python3 tests/db_search/make_blastx_inputs.py. Queries are slices of tests/fasta/MelaMJNV.fasta
(a shrimp virus genome; 1-based coordinates in the comments), subjects are proteins of tests/fasta/PemoMJNVA.faa
(a related virus; NCBI blastx -subject finds them in the genome). The decoy proteins come from random.Random(20261018), the query without hits from random.Random(20261019).
"""
import random
from pathlib import Path

rng = random.Random(20261018)


def read(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def proteins(path: str) -> list[tuple[str, str, str]]:
    out: list[list[str]] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name, _, title = line[1:].partition(" ")
            out.append([name, title, ""])
        else:
            out[-1][2] += line.strip()
    return [tuple(r) for r in out]  # type: ignore[misc]


def fasta(records: list[tuple[str, str]], width: int) -> str:
    return "".join(f">{name}\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width)) for name, seq in records)


M = read("tests/fasta/MelaMJNV.fasta")
write = Path("tests/db_search")

# Database: the PemoMJNVA proteins that hit the two query slices (a minus-strand family of five related proteins
# and its neighbours, a plus- and a minus-strand gene) and decoys: 24 more proteins of at most 1,000 residues,
# chosen with the seeded generator.
hits = ["BDT61857.1", "BDT61858.1", "BDT61861.1", "BDT61866.1", "BDT61876.1", "BDT61959.1", "BDT61960.1",
        "BDT61961.1", "BDT61962.1", "BDT61963.1", "BDT61964.1"]
table = proteins("tests/fasta/PemoMJNVA.faa")
decoys = rng.sample([r for r in table if r[0] not in hits and len(r[2]) <= 1000], 24)
chosen = [r for r in table if r[0] in hits or r in decoys]  # in the order of the source file
records = [(f"{n} {t}", s) for n, t, s in chosen]
(write / "bslice.faa").write_text(fasta(records, 60))

# Queries: m_fam MelaMJNV 26001-31000 (minus-strand genes BDT61960 and BDT61961 and a family of related ones),
# m_mix MelaMJNV 56001-61600 (BDT61861 plus strand, BDT61866 minus strand).
(write / "bslice_q.fna").write_text(fasta([
    ("m_fam MelaMJNV 26001-31000", M[26000:31000]),
    ("m_mix MelaMJNV 56001-61600", M[56000:61600]),
], 70))
# Lower case stretches inside m_fam's genes (for -lcase_masking).
lc = M[26000:29000]
lc = lc[:300] + lc[300:700].lower() + lc[700:1500] + lc[1500:1800].lower() + lc[1800:]
(write / "bslice_q_lc.fna").write_text(fasta([("m_fam_lc MelaMJNV 26001-29000, 301-700 and 1501-1800 in lower case", lc)], 70))
# A random query without hits at the default E value (found with NCBI blastx: of 40 random candidates of 300 bases,
# from random.Random(20261019), only number 25 has none in `bslice`).
rng = random.Random(20261019)
candidates = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(40)]
(write / "bslice_nohit.fna").write_text(fasta([("q_nohit 300 random bases", candidates[25])], 70))
(write / "bslice_alln.fna").write_text(fasta([("q_alln 300 N", "N" * 300)], 70))

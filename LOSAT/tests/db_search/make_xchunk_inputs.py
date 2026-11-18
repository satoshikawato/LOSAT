#!/usr/bin/env python3
"""Writes the inputs of the rows `db.xc.*` (line RG item X-chunk): BLASTN query chunks that NCBI BLAST+ 2.17.0 cannot
make (`X/x2/xchunk/REPORT.md` of the RG task folder).

Run from LOSAT/: python3 tests/db_search/make_xchunk_inputs.py. The files are committed:
`tests/fasta/outfmt0/xch_s.fna` (the subjects `pa`, 4000 random bases, and `pb`, 1500) and the query pairs below, whose
queries are pieces of `pa` with every 37th base changed (`mut`), so that the first query of a pair hits `pa`. The
random bases come from a linear congruential generator (a fixed sequence on every Python). With `BATCH_SIZE=1`
the first query is the first query batch and the second the second:

- `xch_under`: q70 and q90. CHUNK_SIZE=140 OVERLAP_CHUNK_SIZE=100 split q90 into 2 chunks of 96, below the overlap:
  the start of the second wraps past the query (`size_t`), so the first chunk covers q90 and the second holds no query.
- `xch_p95`: q9 and q31. CHUNK_SIZE=100 OVERLAP_CHUNK_SIZE=95 make q31 6 chunks of 85, the first covering it.
- `xch_neg`: q150 and q310. CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200 make q310 3 chunks, the first covering it.
- `xch_gaps`: q150 and q248. The same pair makes q248 2 chunks of 25 with a gap of 200 between them.
- `xch_ctl`: q300 and q1000. CHUNK_SIZE=300 OVERLAP_CHUNK_SIZE=100 split q1000 into 5 chunks of 281, none split again.
- `xch_again`: q70 and q1000. CHUNK_SIZE=150 OVERLAP_CHUNK_SIZE=100 make q1000 20 chunks of 146, which NCBI would
  split again (PD-LOSAT-NCBI-DEFECTS exception 1).

The database `tests/db_search/xch_s` is `xch_s.fna` (`makeblastdb -dbtype nucl -blastdb_version 4 -title "xch
subjects"`, made once with NCBI BLAST+ 2.17.0 on ws).
"""
from pathlib import Path

OUT = Path("tests/fasta/outfmt0")
state = 20261118


def base():
    global state
    state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
    return "ACGT"[(state >> 33) & 3]


def dna(n):
    return "".join(base() for _ in range(n))


def mut(seq, every=37):
    change = {"A": "C", "C": "G", "G": "T", "T": "A"}
    seq = list(seq)
    for i in range(every // 2, len(seq), every):
        seq[i] = change[seq[i]]
    return "".join(seq)


def write(name, records):
    lines = []
    for ident, seq in records:
        lines.append(f">{ident}")
        lines.extend(seq[i : i + 70] for i in range(0, len(seq), 70))
    (OUT / name).write_text("\n".join(lines) + "\n")


PA = dna(4000)
PB = dna(1500)
write("xch_s.fna", [("pa", PA), ("pb", PB)])


def q(ident, offset, n):
    return (ident, mut(PA[offset : offset + n]))


for name, records in {
    "xch_under": [q("q70", 100, 70), q("q90", 500, 90)],
    "xch_p95": [q("q9", 100, 9), q("q31", 500, 31)],
    "xch_neg": [q("q150", 100, 150), q("q310", 500, 310)],
    "xch_gaps": [q("q150", 100, 150), q("q248", 500, 248)],
    "xch_ctl": [q("q300", 50, 300), q("q1000", 300, 1000)],
    "xch_again": [q("q70", 100, 70), q("q1000", 300, 1000)],
}.items():
    write(name + ".fna", records)

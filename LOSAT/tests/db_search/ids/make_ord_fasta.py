#!/usr/bin/env python3
"""Write the FASTA of the line DBID N2a databases (`make_ids.py` builds them) and their queries.

Sequences whose only Seq-id is `gnl|BL_ORD_ID|n`. makeblastdb without -parse_seqids gives each sequence
that id with n its number in the volume; with -parse_seqids it stores such a FASTA id as written, and
NCBI's formatters look the sequence up by n plus the first ordinal id of its volume. Unrelated random
proteins (seed 20261102), each the query `pord_q<k>.faa` of its own record k (`pord_q<k>.fna` a
back-translation); the nucleotide records are back-translations of other ones (queries `nord_q<k>.fna`,
`nord_q<k>.faa` the protein).

`pord.faa`: 0 and 1 `gb|FFF00006.1|` and `gb|FFF00007.1|` (a first record with parsed Seq-ids), 2
`gnl|BL_ORD_ID|2` (its own number, which NCBI finds), 3 `gnl|BL_ORD_ID|6` (no sequence 6 in the database;
through the alias `pord_vx.pal` the second sequence of `pord_x`, whose id is its own number: NCBI writes that
sequence's names), 4 `gnl|BL_ORD_ID|1` (sequence 1, which holds another Seq-id). `pord_x.faa`, the second
database of `pord_vx.pal` (NCBI orders the volumes of an alias by their base names; queries `pord_x_q<k>.faa`):
0 `gb|FFF00008.1|`, 1 `gnl|BL_ORD_ID|1` (its own number in its volume), 2 `gnl|BL_ORD_ID|7` (its ordinal id in
the alias, not in its volume). No record names `gnl|BL_ORD_ID|0`, which makeblastdb stores with a string tag.
`nord.fna`: 0 `gnl|BL_ORD_ID|3` (no such sequence; the first record of a volume tells whether the database
has parsed Seq-ids), 1 `gnl|BL_ORD_ID|1`, 2 `gb|AF000060.1|`; `nord_qall.fna` holds the queries of records 2, 1
and 0 in this order, so that a report stream would write the first two before it reaches record 0.

usage: make_ord_fasta.py [--out LOSAT/tests/db_search/ids]
"""
import argparse
import random
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
         "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT"}
PORD = ["gb|FFF00006.1| control with an accession",
        "gb|FFF00007.1| second control",
        "gnl|BL_ORD_ID|2 ordinal id of its own number",
        "gnl|BL_ORD_ID|6 ordinal id of no sequence here",
        "gnl|BL_ORD_ID|1 ordinal id of the second control"]
PORD_X = ["gb|FFF00008.1| control of the second database",
          "gnl|BL_ORD_ID|1 own number in the second database",
          "gnl|BL_ORD_ID|7 number in the alias, not in its volume"]
NORD = ["gnl|BL_ORD_ID|3 ordinal id of no sequence, first",
        "gnl|BL_ORD_ID|1 ordinal id of its own number",
        "gb|AF000060.1| control with an accession"]


def wrap(seq: str) -> str:
    return "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    out = ap.parse_args().out
    rng = random.Random(20261102)

    def protein() -> str:
        return "M" + "".join(rng.choice(AA) for _ in range(219))

    for name, titles in (("pord", PORD), ("pord_x", PORD_X)):
        records = []
        for k, title in enumerate(titles):
            seq = protein()
            records.append(f">{title}\n{wrap(seq)}")
            (out / f"{name}_q{k}.faa").write_text(f">{name}q{k} query of record {k}\n{wrap(seq)}")
            if name == "pord":
                (out / f"{name}_q{k}.fna").write_text(
                    f">{name}q{k}n back-translation of the query of record {k}\n"
                    f"{wrap(''.join(CODON[a] for a in seq))}")
        (out / f"{name}.faa").write_text("".join(records))
    records = []
    queries = []
    for k, title in enumerate(NORD):
        seq = protein()
        nseq = "".join(CODON[a] for a in seq)
        records.append(f">{title}\n{wrap(nseq)}")
        queries.insert(0, f">nordq{k} query of record {k}\n{wrap(nseq)}")
        (out / f"nord_q{k}.fna").write_text(queries[0])
        (out / f"nord_q{k}.faa").write_text(f">nordq{k}p protein of record {k}\n{wrap(seq)}")
    (out / "nord.fna").write_text("".join(records))
    (out / "nord_qall.fna").write_text("".join(queries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

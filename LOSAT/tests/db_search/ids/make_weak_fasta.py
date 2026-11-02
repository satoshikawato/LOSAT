#!/usr/bin/env python3
"""Write the FASTA of the line DBID N2c databases (`make_ids.py` builds them) and their queries.

NCBI's scope matches Seq-ids weakly (`CSeq_id_Textseq_Tree::FindMatch`): an accession without a version finds the
same accession with a version, an id of accession and name finds the name-only id of that name, a name-only id finds
an id with that name and an accession, and a pir or prf accession-only id finds the name-only id with that text as
its name; not the other way (`gb|X.1|` does not find `gb|X|`, `pir||ABCD` does not find `pir|ABCD|`). Two
consequences, both refused by LOSAT for databases with parsed Seq-ids: a query read with -parse_deflines whose id a
reported sequence's Seq-align id matches weakly (NCBI shows the query for the subject at the default prefetch), and,
with PRE_FETCH_SEQS_LIMIT set to anything but 2147483647 (NCBI reads the reported sequences one by one when the
limit is 0 or below the batch's Seq-align count), a reported sequence S whose Seq-align id matches weakly an id of
another database sequence T that is read first (NCBI writes T's names for S). Unrelated random proteins (seed
20261104), ~220 letters; an "S" record is its "T" record with every tenth letter replaced (a query of T hits T at
100 % and S second; for the nucleotide records every fortieth base).

`pweak.faa`: 0 `sp|P00041.1|NAMEW`, 1 `gb|FFF00040|`, 2 `pir|PIRW1|`, 3 `sp||NAMEX`, 4 `gb|FFF00041.1|`, 5
`sp||NAMEY` (S of 6), 6 `sp|P00042.1|NAMEY` (T), 7 `gb|FFF00043|` (S of 8), 8 `gb|FFF00043.1|` (T), 9
`sp|P00044.1|NAMEZ` (S of 10), 10 `sp||NAMEZ` (T). Queries `pweak_q<k>.faa` (`>q<k>`, record k's letters, no
-parse_deflines) for k = 0, 3, 4, 6, 8, 10 and `pweak_q8.fna` (back-translation, BLASTX). Queries read with
-parse_deflines: `pweak_qa0.faa` `sp||NAMEW` (record 0), `pweak_qa1.faa` `gb|FFF00040.1|` (record 1),
`pweak_qa2.faa` `pir||PIRW1` (record 2), `pweak_qa3.faa` `sp|P00045.1|NAMEX` (record 3), `pweak_qa4.faa`
`gb|FFF00041|` (record 4; the direction NCBI does not match) and `pweak_qa1.fna`, record 1 back-translated, the
same defline (BLASTX).

`nweak.fna`: 0 `gb|AF000081|` (S of 1), 1 `gb|AF000081.1|` (T), 2 `gb|AF000082|`, back-translations of random
proteins. Queries `nweak_q1.fna` (record 1), `nweak_q1.faa` (its protein, TBLASTN) and `nweak_qa2.fna`
(`gb|AF000082.1|`, record 2, read with -parse_deflines).

usage: make_weak_fasta.py [--out LOSAT/tests/db_search/ids]
"""
import argparse
import random
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
         "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT"}
PWEAK = ["sp|P00041.1|NAMEW", "gb|FFF00040|", "pir|PIRW1|", "sp||NAMEX", "gb|FFF00041.1|", "sp||NAMEY",
         "sp|P00042.1|NAMEY", "gb|FFF00043|", "gb|FFF00043.1|", "sp|P00044.1|NAMEZ", "sp||NAMEZ"]
PWEAK_S = (5, 7, 9)  # the records that are the near copy of the next one
PWEAK_QUERIES = (0, 3, 4, 6, 8, 10)
PWEAK_QA = [("qa0", "sp||NAMEW", 0), ("qa1", "gb|FFF00040.1|", 1), ("qa2", "pir||PIRW1", 2),
            ("qa3", "sp|P00045.1|NAMEX", 3), ("qa4", "gb|FFF00041|", 4)]
NWEAK = ["gb|AF000081|", "gb|AF000081.1|", "gb|AF000082|"]


def wrap(seq: str) -> str:
    return "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))


def back(seq: str) -> str:
    return "".join(CODON[a] for a in seq)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    out = ap.parse_args().out
    rng = random.Random(20261104)

    def protein() -> str:
        return "M" + "".join(rng.choice(AA) for _ in range(219))

    def near_copy(seq: str, every: int, alphabet: str) -> str:
        return "".join(rng.choice(alphabet.replace(a, "")) if i % every == every - 1 else a
                       for i, a in enumerate(seq))

    seqs: list[str] = []
    for k in range(len(PWEAK)):
        if k in PWEAK_S:
            target = protein()
            seqs += [near_copy(target, 10, AA), target]
        elif k - 1 not in PWEAK_S:
            seqs.append(protein())
    (out / "pweak.faa").write_text("".join(
        f">{name} w{k}_title\n{wrap(seq)}" for k, (name, seq) in enumerate(zip(PWEAK, seqs))))
    for k in PWEAK_QUERIES:
        (out / f"pweak_q{k}.faa").write_text(f">q{k}\n{wrap(seqs[k])}")
    (out / "pweak_q8.fna").write_text(f">q8\n{wrap(back(seqs[8]))}")
    for name, ident, k in PWEAK_QA:
        (out / f"pweak_{name}.faa").write_text(f">{ident} qa_title\n{wrap(seqs[k])}")
    (out / "pweak_qa1.fna").write_text(f">gb|FFF00040.1| qa_title\n{wrap(back(seqs[1]))}")

    target = protein()
    nseqs = [None, back(target), back(protein())]
    nseqs[0] = near_copy(nseqs[1], 40, "ACGT")
    (out / "nweak.fna").write_text("".join(
        f">{name} n{k}_title\n{wrap(seq)}" for k, (name, seq) in enumerate(zip(NWEAK, nseqs))))
    (out / "nweak_q1.fna").write_text(f">nq1\n{wrap(nseqs[1])}")
    (out / "nweak_q1.faa").write_text(f">nq1p\n{wrap(target)}")
    (out / "nweak_qa2.fna").write_text(f">gb|AF000082.1| qa_title\n{wrap(nseqs[2])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

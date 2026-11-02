#!/usr/bin/env python3
"""Write the FASTA of the line DBID N2b databases (`make_ids.py` builds them) and their queries.

Seq-align ids that LOSAT refused by their kind before line DBID N2b, and that NCBI's formatters look up in the
database like any other (`CSeqDB::SeqidToOid`): version 5 finds `bbs|`, `nat|` and `lcl|<integer>` by their LMDB
keys unless an earlier sequence holds the same key, version 4 finds neither `bbs|` nor `nat|` in its ISAM files,
no version finds `gnl|PIG|` (makeblastdb fills the PIG files from no FASTA id), and `gnl|BL_ORD_ID|n` is n plus
the first ordinal id of a volume, so a sequence finds itself by its own number also beside another id. Unrelated
random proteins (seed 20261103), each the query `plook_q<k>.faa` of its own record k (`plook_q<k>.fna` a
back-translation); the nucleotide records are back-translations of other ones (queries `nlook_q<k>.fna`,
`nlook_q<k>.faa` the protein).

`plook.faa`: 0 `gb|FFF00020.1|` (control), 1 `bbs|123`, 2 `nat|NA000001.1|`, 3 `lcl|777`, 4 `gnl|PIG|124`,
5 `bbs|300` then 6 `lcl|300` (the key `300` gives sequence 5 first), 7 `lcl|301` then 8 `bbs|301` (sequence 7
first), 9 `gnl|BL_ORD_ID|9` with `gb|FFF00021.1|` in a second definition line (its own number), 10
`gnl|BL_ORD_ID|10|gb|FFF00022.1|` (its own number, one line), 11 `gnl|BL_ORD_ID|4` beside `gb|FFF00023.1|` (the
number of sequence 4), 12 `gnl|BL_ORD_ID|40` beside `gb|FFF00024.1|` (no such sequence), 13 `bbs|500|gb|FFF00025.1|`,
14 `lcl|778` with `gb|FFF00026.1|` in a second definition line. `plook_z.faa`, the second database of the alias
`plook_vz.pal` (queries `plook_z_q<k>.faa`): 0 `gb|FFF00027.1|`, 1 `gnl|BL_ORD_ID|1` (its own number in its
volume) with `gb|FFF00028.1|` in a second definition line. `nlook.fna`: 0 `gb|AF000070.1|`, 1 `bbs|123`, 2
`nat|NA000002.1|`, 3 `lcl|777`, 4 `bbs|300` then 5 `lcl|300`, 6 `gnl|BL_ORD_ID|6` with `gb|AF000071.1|` in a second
definition line; `nlook_qall.fna` holds the queries of records 1, 0 and 5 in this order (a report stream would
write the first two before it reaches record 5), `nlook_qfound.fna` those of records 1, 2, 3 and 6.

A Seq-align id that another sequence also holds, which NCBI's scope shows for it where both are reported (review 1
of N2b): `pshare.faa` (versions 5 and 4) 0 `lcl|555`, 1 `gb|FFF00030.1|` with `lcl|555` in a second definition line
(record 0's sequence with every tenth letter replaced, so that its query `pshare_q1.faa` hits both), 2
`gnl|BL_ORD_ID|0` (which makeblastdb stores with a text tag, so NCBI's lookup throws) with `gb|FFF00031.1|` in a
second definition line, 3 `gb|FFF00032.1|`; `pshare_a.faa` `bbs|123` and `pshare_b.faa` `gb|FFF00033.1|` with
`bbs|123` in a second definition line (makeblastdb refuses the two in one database), read through the alias
`pshare_ab.pal` (query `pshare_b_q0.faa`); `nshare.fna` the back-translations of `pshare.faa` 0 and 1 (query
`nshare_q1.fna`).

usage: make_look_fasta.py [--out LOSAT/tests/db_search/ids]
"""
import argparse
import random
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
         "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT"}
PLOOK = ["gb|FFF00020.1| control with an accession",
         "bbs|123 gibbsq id",
         "nat|NA000001.1| named annotation track id",
         "lcl|777 local integer id",
         "gnl|PIG|124 PIG id",
         "bbs|300 gibbsq id before lcl|300",
         "lcl|300 local integer id after bbs|300",
         "lcl|301 local integer id before bbs|301",
         "bbs|301 gibbsq id after lcl|301",
         "gnl|BL_ORD_ID|9 ordinal id of its own number\x01gb|FFF00021.1| second definition line",
         "gnl|BL_ORD_ID|10|gb|FFF00022.1| ordinal id of its own number beside an accession",
         "gnl|BL_ORD_ID|4 ordinal id of the PIG sequence\x01gb|FFF00023.1| second definition line",
         "gnl|BL_ORD_ID|40 ordinal id of no sequence\x01gb|FFF00024.1| second definition line",
         "bbs|500|gb|FFF00025.1| gibbsq id beside an accession",
         "lcl|778 local integer id\x01gb|FFF00026.1| second definition line"]
PLOOK_Z = ["gb|FFF00027.1| control of the second database",
           "gnl|BL_ORD_ID|1 own number in the second database\x01gb|FFF00028.1| second definition line"]
NLOOK = ["gb|AF000070.1| control with an accession",
         "bbs|123 gibbsq id",
         "nat|NA000002.1| named annotation track id",
         "lcl|777 local integer id",
         "bbs|300 gibbsq id before lcl|300",
         "lcl|300 local integer id after bbs|300",
         "gnl|BL_ORD_ID|6 ordinal id of its own number\x01gb|AF000071.1| second definition line"]

PSHARE = ["lcl|555 local integer id",
          "gb|FFF00030.1| accession\x01lcl|555 the same local id in a second definition line",
          "gnl|BL_ORD_ID|0 ordinal id stored with a text tag\x01gb|FFF00031.1| second definition line",
          "gb|FFF00032.1| control with an accession"]


def wrap(seq: str) -> str:
    return "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    out = ap.parse_args().out
    rng = random.Random(20261103)

    def protein() -> str:
        return "M" + "".join(rng.choice(AA) for _ in range(219))

    for name, titles in (("plook", PLOOK), ("plook_z", PLOOK_Z)):
        records = []
        for k, title in enumerate(titles):
            seq = protein()
            records.append(f">{title}\n{wrap(seq)}")
            (out / f"{name}_q{k}.faa").write_text(f">{name}q{k} query of record {k}\n{wrap(seq)}")
            if name == "plook":
                (out / f"{name}_q{k}.fna").write_text(
                    f">{name}q{k}n back-translation of the query of record {k}\n"
                    f"{wrap(''.join(CODON[a] for a in seq))}")
        (out / f"{name}.faa").write_text("".join(records))
    records = []
    queries = {}
    for k, title in enumerate(NLOOK):
        seq = protein()
        nseq = "".join(CODON[a] for a in seq)
        records.append(f">{title}\n{wrap(nseq)}")
        queries[k] = f">nlookq{k} query of record {k}\n{wrap(nseq)}"
        (out / f"nlook_q{k}.fna").write_text(queries[k])
        (out / f"nlook_q{k}.faa").write_text(f">nlookq{k}p protein of record {k}\n{wrap(seq)}")
    (out / "nlook.fna").write_text("".join(records))
    (out / "nlook_qall.fna").write_text("".join(queries[k] for k in (1, 0, 5)))
    (out / "nlook_qfound.fna").write_text("".join(queries[k] for k in (1, 2, 3, 6)))
    first = protein()
    near = "".join(rng.choice(AA.replace(a, "")) if i % 10 == 9 else a for i, a in enumerate(first))
    seqs = [first, near, protein(), protein()]
    (out / "pshare.faa").write_text("".join(f">{t}\n{wrap(q)}" for t, q in zip(PSHARE, seqs)))
    (out / "pshare_q1.faa").write_text(f">pshareq1 query of record 1\n{wrap(near)}")
    (out / "pshare_q2.faa").write_text(f">pshareq2 query of record 2\n{wrap(seqs[2])}")
    (out / "nshare.fna").write_text("".join(
        f">{t}\n{wrap(''.join(CODON[a] for a in q))}" for t, q in zip(PSHARE[:2], seqs[:2])))
    (out / "nshare_q1.fna").write_text(f">nshareq1 query of record 1\n{wrap(''.join(CODON[a] for a in near))}")
    first = protein()
    near = "".join(rng.choice(AA.replace(a, "")) if i % 10 == 9 else a for i, a in enumerate(first))
    (out / "pshare_a.faa").write_text(f">bbs|123 gibbsq id\n{wrap(first)}")
    (out / "pshare_b.faa").write_text(
        f">gb|FFF00033.1| accession\x01bbs|123 the same gibbsq id in a second definition line\n{wrap(near)}")
    (out / "pshare_b_q0.faa").write_text(f">pshare_bq0 query of record 0\n{wrap(near)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

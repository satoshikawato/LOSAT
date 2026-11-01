#!/usr/bin/env python3
"""Write the FASTA of the line DBID N1 databases (`make_ids.py` builds them) and their queries.

Five protein pairs of unrelated random sequences (seed 20261101) in `pfound.faa`; the second record of
each pair is its first's sequence with three substitutions, so the query `pfound_q<k>.faa` of pair k
(the first's sequence) hits both. Pair 1: gb|AAA00001.1| then gb|AAA00001|; pair 2: gb|BBB00002.1| then
emb|BBB00002.1|; pair 3: gb|CCC00003.1| then lcl|CCC00003; pair 4: gb|DDD00004| then gb|DDD00004.1|;
pair 5: pdb|2XYZ| then pdb|2XYZ|AA. `pfound_q1.fna` and `pfound_q4.fna` are back-translations of the
queries of pairs 1 and 4. `pdup.faa` holds one record gb|EEE00005.1| (query `pdup_q.faa`), `pdup_c.faa`
the same record with another title. `nfound.fna`: the sequence of `njoined_a.fna` as gb|AF000050.1| and,
with three substitutions, gb|AF000050|. `../seqidlist/sl_pfound_neg.txt`: AAA00001.1.

usage: make_found_fasta.py [--out LOSAT/tests/db_search/ids]
"""
import argparse
import random
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
         "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT",
         "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT"}
PAIRS = [("gb|AAA00001.1| pair one, versioned", "gb|AAA00001| pair one, unversioned"),
         ("gb|BBB00002.1| pair two, GenBank", "emb|BBB00002.1| pair two, EMBL"),
         ("gb|CCC00003.1| pair three, GenBank", "lcl|CCC00003 pair three, local"),
         ("gb|DDD00004| pair four, unversioned first", "gb|DDD00004.1| pair four, versioned second"),
         ("pdb|2XYZ| pair five, no chain", "pdb|2XYZ|AA pair five, chain AA")]


def wrap(seq: str) -> str:
    return "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))


def mutate(seq: str, alphabet: str, positions) -> str:
    s = list(seq)
    for p in positions:
        s[p] = next(c for c in alphabet if c != s[p])
    return "".join(s)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    out = ap.parse_args().out
    rng = random.Random(20261101)
    records = []
    for k, (first, second) in enumerate(PAIRS, 1):
        seq = "M" + "".join(rng.choice(AA) for _ in range(219))
        records.append(f">{first}\n{wrap(seq)}>{second}\n{wrap(mutate(seq, AA, (40, 110, 180)))}")
        (out / f"pfound_q{k}.faa").write_text(f">pq{k} query of pair {k}\n{wrap(seq)}")
        if k in (1, 4):
            (out / f"pfound_q{k}.fna").write_text(
                f">pq{k}n back-translation of the query of pair {k}\n{wrap(''.join(CODON[a] for a in seq))}")
    (out / "pfound.faa").write_text("".join(records))
    seq = "M" + "".join(rng.choice(AA) for _ in range(219))
    (out / "pdup.faa").write_text(f">gb|EEE00005.1| duplicated record\n{wrap(seq)}")
    (out / "pdup_c.faa").write_text(f">gb|EEE00005.1| duplicated record, another title\n{wrap(seq)}")
    (out / "pdup_q.faa").write_text(f">pdq query of the duplicated record\n{wrap(seq)}")
    nseq = "".join(line.strip() for line in (out / "njoined_a.fna").read_text().splitlines()[1:])
    (out / "nfound.fna").write_text(
        f">gb|AF000050.1| nucleotide pair, versioned\n{wrap(nseq)}"
        f">gb|AF000050| nucleotide pair, unversioned\n{wrap(mutate(nseq, 'ACGT', (100, 250, 400)))}")
    (out.parent / "seqidlist" / "sl_pfound_neg.txt").write_text("AAA00001.1\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

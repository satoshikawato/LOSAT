#!/usr/bin/env python3
"""Write the FASTA of the line DBID N3 databases (`make_ids.py --join-only` builds them) and their queries.

Copies of one record that follow each other in a hit list, which NCBI's XML, XML2, JSON and outfmt 0 reports show
as one hit. A random protein P1 of 200 residues (seed 20261103) and its back-translation N1 (600 nucleotides, the
queries `pjoin_q.faa` and `njoin_q.fna`, and `njoin_q2.fna` with N1 twice for two query batches); the protein records
in `pjoin_a.faa`, `pjoin_b.faa`, `pjoin_c.faa` and the nucleotide records of the same layout in `njoin_a.fna`,
`njoin_b.fna`, `njoin_c.fna`:
  r1 gb|JHA00001.1|  P1 / N1, in a, b and c (three copies)
  r2 gb|JHA00002.1|  in a: about 25 % / 15 % substitutions
  r3 gb|JHA00003.1|  with a second definition line ref|XP_900003.1|, in a and b: two pieces of the query apart (two
                     HSPs in each copy)
  r4 gb|JHA00004.1|  in b: about 40 % / 25 % substitutions
  r6 gb|JHA00006.1|  in b: the letters of r1 under another id; it ties with r1's copies and comes between the copy
                     in c and those in a and b (NCBI orders ties by the ordinal id, the later first)
`pjoin_x.faa` and `njoin_x.fna` hold the definition lines of r1 and r3 with other letters: r1 with four substitutions,
r3 with ten more letters at its end (no copies, though the definition lines are the same). `njoin_t.fna` holds r10,
gb|JHA00010.1| with a second definition line lcl|loc10, with the letters of N1, for two databases with taxid 9606.

usage: make_join_fasta.py [--out LOSAT/tests/db_search/ids]
"""
import argparse
import random
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"], "C": ["TGT", "TGC"], "D": ["GAT", "GAC"], "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"], "G": ["GGT", "GGC", "GGA", "GGG"], "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "M": ["ATG"], "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "Q": ["CAA", "CAG"], "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"], "W": ["TGG"], "Y": ["TAT", "TAC"],
}
TITLES = {
    "r1": "gb|JHA00001.1| joined record one",
    "r2": "gb|JHA00002.1| other record two",
    "r3": "gb|JHA00003.1| two-piece record three\x01ref|XP_900003.1| merged definition line of three",
    "r4": "gb|JHA00004.1| distant record four",
    "r6": "gb|JHA00006.1| twin of record one",
}
LAYOUT = {"a": ["r1", "r2", "r3"], "b": ["r1", "r3", "r4", "r6"], "c": ["r1"]}


def mutate(seq: str, rate: float, alphabet: str, rng: random.Random) -> str:
    return "".join(rng.choice([a for a in alphabet if a != c]) if rng.random() < rate else c for c in seq)


def wrap(seq: str) -> str:
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    out = ap.parse_args().out
    rng = random.Random(20261103)
    p1 = "M" + "".join(rng.choice(AA) for _ in range(199))
    n1 = "".join(rng.choice(CODONS[a]) for a in p1)
    prot = {
        "r1": p1,
        "r2": mutate(p1, 0.25, AA, rng),
        "r3": p1[:90] + "".join(rng.choice(AA) for _ in range(120)) + p1[110:],
        "r4": mutate(p1, 0.40, AA, rng),
        "r6": p1,
    }
    nucl = {
        "r1": n1,
        "r2": mutate(n1, 0.15, "ACGT", rng),
        "r3": n1[:270] + "".join(rng.choice("ACGT") for _ in range(360)) + n1[330:],
        "r4": mutate(n1, 0.25, "ACGT", rng),
        "r6": n1,
    }
    for db, names in LAYOUT.items():
        for kind, seqs, extension in (("pjoin", prot, "faa"), ("njoin", nucl, "fna")):
            (out / f"{kind}_{db}.{extension}").write_text(
                "".join(f">{TITLES[name]}\n{wrap(seqs[name])}" for name in names))
    def substituted(seq: str, alphabet: str) -> str:
        letters = list(seq)
        for position in (20, 60, 100, 140):
            letters[position] = next(a for a in alphabet if a != letters[position])
        return "".join(letters)

    (out / "pjoin_x.faa").write_text(
        f">{TITLES['r1']}\n{wrap(substituted(p1, AA))}>{TITLES['r3']}\n{wrap(prot['r3'] + 'ACDEFGHIKL')}")
    (out / "njoin_x.fna").write_text(
        f">{TITLES['r1']}\n{wrap(substituted(n1, 'ACGT'))}>{TITLES['r3']}\n{wrap(nucl['r3'] + 'ACGTACGTAC')}")
    (out / "njoin_t.fna").write_text(
        f">gb|JHA00010.1| first title of ten\x01lcl|loc10 second title of ten\n{wrap(n1)}")
    (out / "pjoin_q.faa").write_text(f">pjq query of the joined records\n{wrap(p1)}")
    (out / "njoin_q.fna").write_text(f">njq query of the joined records\n{wrap(n1)}")
    (out / "njoin_q2.fna").write_text(
        f">njq1 query of the joined records\n{wrap(n1)}>njq2 the same, a second batch\n{wrap(n1)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Write the TBLASTN tblastn-fast AMR fixture inputs (line DI4; README.md in this directory).

From the ABRicate `ncbi` nucleotide database FASTA (NCBI's AMR reference genes, ABRicate's
db/ncbi/sequences, records named `ncbi~~~SYMBOL~~~ACCESSION~~~...`), write
- OUT/tfast_amr_q.faa: protein queries (translations, standard code) of genes that are in the
  subject, variants of them (other alleles of a family) and genes that are not;
- OUT/tfast_amr_s.fna: two synthetic contigs of random spacer DNA (seeded) with genes inserted
  forward, reverse-complemented, with point mutations (about 5 % of the bases) and truncated,
  and a run of N.

Usage: make_tfast_amr_inputs.py ABRICATE_NCBI_FNA OUT_DIR
"""
import random
import re
import sys
from pathlib import Path

CODE = {}
_bases = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_bases):
    for j, b in enumerate(_bases):
        for k, c in enumerate(_bases):
            CODE[a + b + c] = _aa[16 * i + 4 * j + k]


def read_fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in records.items()}


def by_accession(records, accession):
    hits = [k for k in records if k.split("~~~")[2] == accession]
    assert len(hits) == 1, (accession, hits)
    return hits[0].split("~~~")[1], records[hits[0]]


def translate(nt):
    protein = "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))
    return protein[:-1] if protein.endswith("*") else protein


def revcomp(nt):
    return nt.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def mutate(rng, nt, rate):
    out = list(nt)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def spacer(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    records = read_fasta(sys.argv[1])
    out = Path(sys.argv[2])
    rng = random.Random(20261029)
    gene = lambda acc: by_accession(records, acc)
    # (accession, form) of the genes of each contig: f forward, r reverse complement,
    # m forward with mutations, rm reverse complement with mutations, t forward, first 60 %.
    contigs = [
        ("contig_1", [("NG_050145.1", "f"), ("NG_049253.1", "r"), ("NG_048935.1", "f"),
                      ("NG_049762.1", "m"), ("NG_048102.1", "f"), ("NG_049326.1", "r")]),
        ("contig_2", [("NG_048153.1", "f"), ("NG_050543.1", "r"), ("NG_047797.1", "rm"),
                      ("NG_048323.1", "t")]),
    ]
    with open(out / "tfast_amr_s.fna", "w") as handle:
        for name, genes in contigs:
            parts = [spacer(rng, rng.randrange(300, 800))]
            for accession, form in genes:
                _, nt = gene(accession)
                if "m" in form:
                    nt = mutate(rng, nt, 0.05)
                if form == "t":
                    nt = nt[: len(nt) * 3 // 5]
                if "r" in form:
                    nt = revcomp(nt)
                parts.append(nt)
                parts.append(spacer(rng, rng.randrange(300, 800)))
            if name == "contig_2":
                parts.insert(2, "N" * 40)
            seq = "".join(parts)
            handle.write(f">{name} synthetic contig with AMR genes ({len(genes)} genes)\n{wrap(seq)}\n")
    # Queries: in the subject (TEM-116 and CTX-M-14 are other alleles of TEM-1 and CTX-M-15),
    # and genes that are not in it (SHV-12 is related to TEM; OXA-1 and APH(3')-Ia are not).
    queries = ["NG_050162.1", "NG_048929.1", "NG_049762.1", "NG_049253.1", "NG_049326.1",
               "NG_048102.1", "NG_048153.1", "NG_050543.1", "NG_047797.1", "NG_048323.1",
               "NG_050590.1", "NG_049392.1", "NG_056047.1"]
    with open(out / "tfast_amr_q.faa", "w") as handle:
        for accession in queries:
            symbol, nt = gene(accession)
            protein = translate(nt)
            ident = re.sub(r"[^A-Za-z0-9.-]", "_", symbol)
            handle.write(f">{ident}_{accession} {symbol} (translation of {accession})\n{wrap(protein, 60)}\n")


if __name__ == "__main__":
    main()

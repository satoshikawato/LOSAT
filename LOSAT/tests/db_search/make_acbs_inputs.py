#!/usr/bin/env python3
"""Writes the inputs of the ADAPTIVE_CBS fixtures (rows db.acbs.*, env parity G2+G3) in tests/db_search:
acbs_q.faa (the BLASTP, PSIBLAST and TBLASTN query), acbs_q.fna (the BLASTX query, the same protein
back-translated), and in DIR the subjects acbs.faa and acbs.fna of the databases `acbs` and `acbs_n`, which
are committed as built (not the FASTA files):

    python3 tests/db_search/make_acbs_inputs.py DIR          (run from LOSAT/)
    makeblastdb -in DIR/acbs.faa -dbtype prot -out tests/db_search/acbs -blastdb_version 4 -title acbs
    makeblastdb -in DIR/acbs.fna -dbtype nucl -out tests/db_search/acbs_n -blastdb_version 4 -title acbs_n

The query is 15 background letters, a 160-letter K/E-rich segment, 10 background letters, a 220-letter
unbiased part U and 15 more. The 1,700 subjects (`s0000 bias`, `s0001 ord`, ... in shuffled order): 700 are
random K/E-rich strings of 80-110 letters, which match the query's biased segment by composition only (raw
E-values about 1e-26 to 1e-53, all removed by composition-based statistics); 1,000 are 40-55-letter windows of U
at 40-60 % identity with 2-8 random letters on each side (raw E-values about 1e-22 to 1e-4, all kept). The
back-translation takes one codon per amino acid. With composition-based statistics NCBI BLAST+ 2.17.0 keeps a
preliminary hit list of 1,050 subjects (700 biased, 350 ordinary: 350 reported); with ADAPTIVE_CBS of 1,500,
which its stream close cuts at 1,100 (400 reported; -max_target_seqs 600: 1,250 without, 1,200 with the
variable, 550 and 500 reported). BLASTX and TBLASTN need -seg no (SEG masks the biased segment).
Random letters come from random.Random(20261108 + 7 * 700 + 1000).
"""
import random
import sys
from pathlib import Path

rng = random.Random(20261108 + 7 * 700 + 1000)
OUT = Path("tests/db_search")
BACKGROUND = dict(zip("ARNDCQEGHILKMFPSTWYV",
                      [7.8, 5.1, 4.5, 5.4, 1.9, 4.3, 6.3, 7.4, 2.3, 6.4, 9.0, 5.7, 2.2, 3.9, 5.2, 7.1, 5.8, 1.3, 3.2, 6.6]))
AA = list(BACKGROUND)
WEIGHTS = [BACKGROUND[a] for a in AA]
CODON = dict(zip("ARNDCQEGHILKMFPSTWYV",
                 "GCT CGT AAT GAT TGT CAA GAA GGT CAT ATT CTG AAA ATG TTT CCT TCT ACT TGG TAT GTT".split()))


def background(n: int) -> str:
    return "".join(rng.choices(AA, WEIGHTS, k=n))


def biased(n: int) -> str:
    return "".join(rng.choices("KERQD", [47, 40, 6, 4, 3], k=n))


def mutate(seq: str, identity: float) -> str:
    return "".join(c if rng.random() < identity else rng.choice(AA) for c in seq)


def back_translate(protein: str) -> str:
    return "".join(CODON[c] for c in protein)


def write(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as out:
        for title, seq in records:
            out.write(f">{title}\n")
            for i in range(0, len(seq), 70):
                out.write(seq[i:i + 70] + "\n")


def main() -> None:
    subjects_dir = Path(sys.argv[1])
    unbiased = background(220)
    segment = biased(160)
    query = background(15) + segment + background(10) + unbiased + background(15)
    subjects = [("bias", biased(rng.randint(80, 110))) for _ in range(700)]
    for _ in range(1000):
        length = rng.randint(40, 55)
        start = rng.randint(0, len(unbiased) - length)
        identity = rng.uniform(0.40, 0.60)
        subjects.append(("ord", background(rng.randint(2, 8)) + mutate(unbiased[start:start + length], identity)
                         + background(rng.randint(2, 8))))
    rng.shuffle(subjects)
    write(OUT / "acbs_q.faa", [("q1 query", query)])
    write(OUT / "acbs_q.fna", [("q1 query_backtranslated", back_translate(query))])
    write(subjects_dir / "acbs.faa", [(f"s{i:04d} {kind}", seq) for i, (kind, seq) in enumerate(subjects)])
    write(subjects_dir / "acbs.fna",
          [(f"s{i:04d} {kind}", back_translate(seq)) for i, (kind, seq) in enumerate(subjects)])


if __name__ == "__main__":
    main()

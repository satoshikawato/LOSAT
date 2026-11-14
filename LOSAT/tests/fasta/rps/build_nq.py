#!/usr/bin/env python3
"""Writes the nucleotide queries of the rpstblastn fixtures (line PS, series PR, item PS12a): nq.fna,
nedge.fna and odd.fna, back-translated from q40.faa here with seeded random codons.

No NCBI counterpart: test inputs. Usage: python3 build_nq.py [OUT_DIR] (default: this directory).

- nq.fna: `stopq` (600 nucleotides of stop codons, frame +1 all stops: a first query whose first
  context is not valid), then 15 records of q40.faa (0-13 and 32; 4 and 32 cut to their first 230
  and 330 letters) back-translated with seed 7: the strands alternate, the left flanks are 0-2
  letters plus 0-6 codons, the right flanks 0-20 letters. 12,440 bytes: two query batches of
  10,000 letters.
- nedge.fna: the edges of a back-translation of record 32 (seed 5, no flanks): its first 1-7, 11
  and 12 letters, its first 570-572 letters and their reverse complements, and 600 letters with a
  run of 30 N, a run of 60 lower-case letters, and the ambiguity codes R Y K M S W.
- odd.fna: a query of 4 letters, an empty one and one of 4 N.
"""

import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# The standard genetic code in NCBI's TCAG order; the codons of each amino acid.
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
BASES = "TCAG"
CODONS: dict[str, list[str]] = {}
for index, letter in enumerate(CODE):
    CODONS.setdefault(letter, []).append(
        BASES[index // 16] + BASES[(index // 4) % 4] + BASES[index % 4]
    )


def reverse_complement(text: str) -> str:
    return text[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    title, letters = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if title is not None:
                records.append((title, "".join(letters)))
            title, letters = line[1:], []
        elif line:
            letters.append(line)
    if title is not None:
        records.append((title, "".join(letters)))
    return records


def back_translate(protein: str, rnd: random.Random) -> str:
    return "".join(rnd.choice(CODONS.get(letter, ["NNN"])) for letter in protein)


def fasta(records: list[tuple[str, str]]) -> str:
    text = []
    for title, letters in records:
        text.append(">" + title + "\n")
        text.extend(letters[k : k + 70] + "\n" for k in range(0, len(letters), 70))
    return "".join(text)


def nq(proteins: list[tuple[str, str]]) -> str:
    stops = random.Random(11)
    records = [("stopq", "".join(stops.choice(("TAA", "TGA", "TAG")) for _ in range(200)))]
    sources = [proteins[k] for k in list(range(14)) + [32]]
    sources[4] = (sources[4][0], sources[4][1][:230])
    sources[14] = (sources[14][0], sources[14][1][:330])
    rnd = random.Random(7)
    for i, (title, protein) in enumerate(sources):
        nucleotides = back_translate(protein.upper().replace("*", ""), rnd)
        flank = (i // 2) % 3
        left = "".join(rnd.choice("ACGT") for _ in range(flank + 3 * rnd.randint(0, 6)))
        right = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(0, 20)))
        letters = left + nucleotides + right
        if i % 2:
            letters = reverse_complement(letters)
        strand = "-" if i % 2 else "+"
        records.append((f"{title.split()[0]} bt strand={strand} flank={flank}", letters))
    return fasta(records)


def nedge(proteins: list[tuple[str, str]]) -> str:
    tk = back_translate(proteins[32][1].upper(), random.Random(5))
    records = [(f"len{n}", tk[:n]) for n in (1, 2, 3, 4, 5, 6, 7, 11, 12)]
    for n in (570, 571, 572):
        records += [(f"tkcut{n}", tk[:n]), (f"tkrc{n}", reverse_complement(tk[:n]))]
    head = tk[:600]
    records.append(("nrun", head[:200] + "N" * 30 + head[230:]))
    records.append(("lower", head[:200] + head[200:260].lower() + head[260:]))
    records.append(("amb", head[:100] + "RYKMSW" + head[106:]))
    return fasta(records)


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
    proteins = read_fasta(HERE / "q40.faa")
    (out / "nq.fna").write_text(nq(proteins))
    (out / "nedge.fna").write_text(nedge(proteins))
    (out / "odd.fna").write_text(">e1\nACGT\n>e2\n\n>e3\nNNNN\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

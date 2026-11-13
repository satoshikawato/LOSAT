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
- nwin.fna and nwin0.fna (series PR, item PS12s: the windows of query batches): 130 and 36 records
  of 60-270 letters (seed 12 and 13) for rows with `BATCH_SIZE` 200 and 100, so that the batches
  are many and windows of 4, 8, 16 and 32 batches (1, 2, 4 and 8 threads) cut them in different
  places. The first records are random letters (no hit at a small e-value); later, every third
  record is a back-translated window of 35-70 letters of one of the three profiles' fragments
  (`rps_mini`: q40.faa records 32, 18 and 4), on alternating strands, between random records. A
  record of 250 random letters followed by one of 210 (nwin0: 110) N makes a batch of one record
  without a valid frame, twice; nwin.fna has an empty record early and, in later batches, a record
  with hyphens, one whose title ends with 25 nucleotides and one with the invalid residue J (the
  reader's warnings).
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


def nwin(proteins: list[tuple[str, str]], count: int, seed: int, n_run: int) -> str:
    rnd = random.Random(seed)
    # The profiles' fragments (tests/db_search/rps/build.sh): q40.faa records 32 (1-190), 18
    # (301-480) and 4 (1101-1290).
    regions = [proteins[32][1][0:190], proteins[18][1][300:480], proteins[4][1][1100:1290]]

    def random_letters(n: int) -> str:
        return "".join(rnd.choice("ACGT") for _ in range(n))

    records: list[tuple[str, str]] = []
    first_hit = count // 3
    hit = 0
    for i in range(count):
        if i in (count // 10, (2 * count) // 3):
            records.append((f"w{i} long", random_letters(250)))
            records.append((f"w{i}n nrun", "N" * n_run))
            continue
        if count > 100 and i == count // 5:
            records.append((f"w{i} empty", ""))
            continue
        if i >= first_hit and i % 3 == 0:
            region = regions[hit % 3]
            size = rnd.randint(35, 70)
            start = rnd.randint(0, len(region) - size)
            letters = random_letters(rnd.randint(0, 10))
            letters += back_translate(region[start : start + size].upper(), rnd)
            letters += random_letters(rnd.randint(0, 10))
            if hit % 2:
                letters = reverse_complement(letters)
            records.append((f"w{i} hit {hit % 3} {start + 1}-{start + size}", letters))
            hit += 1
            continue
        letters = random_letters(rnd.randint(60, 240))
        title = f"w{i} rnd"
        if count > 100 and i == (3 * count) // 4:
            letters = letters[:30] + "-" + letters[30:60] + "-" + letters[60:]
            title = f"w{i} hyphens"
        if count > 100 and i == (4 * count) // 5:
            title = f"w{i} title " + random_letters(25)
        if count > 100 and i == (5 * count) // 6:
            letters = letters[:40] + "J" + letters[40:]
            title = f"w{i} invalid"
        records.append((title, letters))
    return fasta(records)


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
    proteins = read_fasta(HERE / "q40.faa")
    (out / "nq.fna").write_text(nq(proteins))
    (out / "nedge.fna").write_text(nedge(proteins))
    (out / "odd.fna").write_text(">e1\nACGT\n>e2\n\n>e3\nNNNN\n")
    (out / "nwin.fna").write_text(nwin(proteins, 130, 12, 210))
    (out / "nwin0.fna").write_text(nwin(proteins, 36, 13, 110))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

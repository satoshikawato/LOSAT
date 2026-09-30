#!/usr/bin/env python3
"""Compare BLASTN on slices of the repository's genomes with NCBI (Session S07+; comparison only).

The seventh S07+ audit found BLASTN differences on slices of real genomes that the
constructed inputs did not show (the preliminary `-subject_besthit` filter, a crash of the
gap reduction at the start of a subject). This sweep cuts, from fixed seeds, a query of 2
to 20 kb and a subject of 5 to 30 kb from two different genomes of `LOSAT/tests/fasta`,
runs each pair in NCBI and LOSAT with the given options and -outfmt 6, and prints each
pair whose stdout, stderr or exit status differs. Exits 1 when any does.

The `ambiguity` pool (the eighth audit round: NCBI resolves a subject's ambiguity codes
with `CRandom` for the preliminary search) cuts windows of EDL933 around its ambiguity
codes instead. The `lcase` pool (the twelfth audit round: the scan ranges of a subject
with lowercase masks start at the last letter of each mask) writes lowercase islands into
the subject; give it `-lcase_masking`. The `iupac` pool (the fourteenth audit round: the
word extension of NCBI's small-query lookup table reads ambiguity codes of the query as
bases, and its word check drops the words that contain them) writes IUPAC codes into a
short query that the subject holds several copies of.

Usage: slice_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--pool viral|ambiguity|lcase|iupac]
                      [--options "..."] [--seed N] [--cases N] [--jobs N]
"""
from __future__ import annotations

import argparse
import random
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
from run_oracle import ENGINE  # noqa: E402

# The distinct genomes (LC738868 = MjeNMV, LC738870 = PemoMJNVA, LC738871 = PemoMJNVB,
# LC738873 = PeseMJNV and LC738874 = MelaMJNV are the same sequences); none has an
# ambiguity code, which the `ambiguity` pool covers with EDL933 (6641 codes).
GENOMES = ["LC738868", "LC738869", "LC738870", "LC738871", "LC738873", "LC738874", "LC738875",
           "MejoMJNV", "TrcuMJNV", "LvMJNV", "MeenMJNV", "MellatMJNV"]
COMPATIBLE = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
              "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def sequence(name: str) -> str:
    lines = (ENGINE / "tests/fasta" / name).read_text().splitlines()
    return "".join(line.strip() for line in lines[1:])


def viral_pair(seqs: dict[str, str], rng: random.Random) -> tuple[str, str]:
    """A query of 2 to 20 kb and a subject of 5 to 30 kb from two distinct genomes."""
    query_genome, subject_genome = rng.sample(GENOMES, 2)
    query_start = rng.randrange(0, len(seqs[query_genome]) - 30000)
    subject_start = rng.randrange(0, len(seqs[subject_genome]) - 40000)
    return (seqs[query_genome][query_start:query_start + rng.randrange(2000, 20000)],
            seqs[subject_genome][subject_start:subject_start + rng.randrange(5000, 30000)])


def ambiguity_pair(edl933: str, positions: list[int], rng: random.Random) -> tuple[str, str]:
    """A subject window of EDL933 around an ambiguity code, and a query from inside it with
    each ambiguity code replaced by a compatible base: the hits depend on how the
    preliminary search resolves the subject's codes."""
    centre = rng.choice(positions)
    subject_start = max(0, centre - rng.randrange(200, 8000))
    subject = edl933[subject_start:centre + rng.randrange(200, 8000)]
    offset = centre - subject_start
    query_start = max(0, offset - rng.randrange(30, 3000))
    query = subject[query_start:offset + rng.randrange(30, 3000)]
    query = "".join(rng.choice(COMPATIBLE[c]) if c in COMPATIBLE else c for c in query)
    return query, subject


def lcase_pair(seqs: dict[str, str], rng: random.Random) -> tuple[str, str]:
    """A query of 200 to 1000 residues from a genome with 0 to 15% of its letters redrawn,
    and a subject window around it with 3 to 30 lowercase islands of 1 to 40 letters."""
    genome = seqs[rng.choice(GENOMES)]
    start = rng.randrange(1000, len(genome) - 20000)
    rate = rng.choice([0.0, 0.05, 0.15])
    query = "".join(rng.choice("ACGT") if rng.random() < rate else c
                    for c in genome[start:start + rng.choice([200, 500, 1000])])
    subject = list(genome[start - 500:start + 2000])
    for _ in range(rng.randint(3, 30)):
        island = rng.randrange(len(subject))
        for index in range(island, min(len(subject), island + rng.choice([1, 1, 2, 3, 5, 10, 40]))):
            subject[index] = subject[index].lower()
    return query, "".join(subject)


def iupac_pair(seqs: dict[str, str], rng: random.Random) -> tuple[str, str]:
    """A query of 30 to 100 residues from a genome with 1 to 3 IUPAC codes (reverse
    complemented for 30% of the pairs), and a subject of 2 to 4 copies of it with 0 to 6%
    of their letters redrawn, between pieces of the genome."""
    genome = seqs[rng.choice(GENOMES)]
    start = rng.randrange(1000, len(genome) - 5000)
    region = genome[start:start + rng.choice([30, 37, 45, 60, 100])]
    rate = lambda: rng.choice([0.0, 0.03, 0.06])  # noqa: E731
    copies = [
        "".join(rng.choice("ACGT") if rng.random() < copy_rate else c for c in region)
        + genome[start + 1000 + index * 50:start + 1000 + index * 50 + rng.randint(0, 40)]
        for index, copy_rate in enumerate(rate() for _ in range(rng.randint(2, 4)))
    ]
    subject = genome[start - rng.randint(0, 40):start] + "".join(copies)
    query = list(region)
    for _ in range(rng.randint(1, 3)):
        query[rng.randrange(len(query))] = rng.choice("NRYSWKMBDHV")
    query = "".join(query)
    if rng.random() < 0.3:
        query = query[::-1].translate(str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB"))
    return query, subject


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--options", default="")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--cases", type=int, default=120)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--pool", choices=("viral", "ambiguity", "lcase", "iupac"), default="viral")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    if args.pool in ("viral", "lcase", "iupac"):
        seqs = {name: sequence(f"{name}.fasta") for name in GENOMES}
        make = {"viral": viral_pair, "lcase": lcase_pair, "iupac": iupac_pair}[args.pool]
        pair = lambda: make(seqs, rng)  # noqa: E731
    else:
        edl933 = sequence("EDL933.fna")
        positions = [index for index, c in enumerate(edl933) if c.upper() not in "ACGT"]
        pair = lambda: ambiguity_pair(edl933, positions, rng)  # noqa: E731
    cases = []
    for index in range(args.cases):
        query_seq, subject_seq = pair()
        query = args.work / f"{args.pool}_q{args.seed}_{index}.fa"
        subject = args.work / f"{args.pool}_s{args.seed}_{index}.fa"
        query.write_text(f">q{index}\n{query_seq}\n")
        subject.write_text(f">s{index}\n{subject_seq}\n")
        cases.append(["-query", str(query), "-subject", str(subject), *shlex.split(args.options), "-outfmt", "6"])

    def run(argv: list[str]) -> tuple[list[str], bool, int, int]:
        ncbi = subprocess.run([str(args.bin_dir / "blastn"), *argv], capture_output=True)
        ours = subprocess.run([str(args.losat.resolve()), "blastn", *argv], capture_output=True)
        same = (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr)
        return argv, same, ncbi.stdout.count(b"\n"), ours.stdout.count(b"\n")

    differing = 0
    with ThreadPoolExecutor(args.jobs) as pool:
        for argv, same, ncbi_lines, our_lines in pool.map(run, cases):
            if not same:
                differing += 1
                print(f"DIFF\t{argv[1]}\t{argv[3]}\tncbi_lines={ncbi_lines}\tlosat_lines={our_lines}")
    print(f"# pool={args.pool} options={args.options!r} seed={args.seed} cases={len(cases)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())

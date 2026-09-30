#!/usr/bin/env python3
"""Compare BLASTN with several queries with NCBI (Session S07++; comparison only).

NCBI searches the queries in adaptive query batches: about 5000 residues first, then sizes
from the successful initial extensions of the batch before (`CBatchSizeMixer`). This sweep
cuts, from fixed seeds, a subject window of a repository genome and 2 to 40 queries whose
total passes the first batch. It runs each case in NCBI and LOSAT with the given options and
-outfmt 6, and prints each case whose stdout, stderr or exit status differs. Exits 1 when
any does. The genomes are those of docs/evidence/losat_web_e2c/slice_sweep.py.

Usage: batch_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--options "..."] [--seed N]
                      [--cases N] [--jobs N]
"""
from __future__ import annotations

import argparse
import importlib.util
import random
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("e2c_slice_sweep", HERE.parent / "losat_web_e2c" / "slice_sweep.py")
e2c = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(e2c)
GENOMES = e2c.GENOMES


def batch_queries(seqs: dict[str, str], rng: random.Random) -> tuple[list[str], str]:
    """A subject window of 5 to 40 kb and 2 to 40 queries of 30 to 3000 residues: pieces of
    the window with 0 to 10% of their letters redrawn (60%), random sequences (20%) and
    pieces of another genome (20%)."""
    name, other = rng.sample(GENOMES, 2)
    genome = seqs[name]
    start = rng.randrange(1000, len(genome) - 60000)
    subject = genome[start:start + rng.randrange(5000, 40000)]
    queries = []
    for _ in range(rng.randint(2, 40)):
        length = rng.choice([30, 60, 150, 500, 1500, 3000])
        kind = rng.random()
        if kind < 0.6:
            offset = rng.randrange(0, max(1, len(subject) - length))
            rate = rng.choice([0.0, 0.03, 0.1])
            queries.append("".join(rng.choice("ACGT") if rng.random() < rate else c
                                   for c in subject[offset:offset + length]))
        elif kind < 0.8:
            queries.append("".join(rng.choice("ACGT") for _ in range(length)))
        else:
            offset = rng.randrange(1000, len(seqs[other]) - length - 1000)
            queries.append(seqs[other][offset:offset + length])
    return queries, subject


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--options", default="")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--cases", type=int, default=120)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    seqs = {name: e2c.sequence(f"{name}.fasta") for name in GENOMES}
    cases = []
    for index in range(args.cases):
        queries, subject_seq = batch_queries(seqs, rng)
        query = args.work / f"batches_q{args.seed}_{index}.fa"
        subject = args.work / f"batches_s{args.seed}_{index}.fa"
        query.write_text("".join(f">q{index}_{k}\n{seq}\n" for k, seq in enumerate(queries)))
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
    print(f"# options={args.options!r} seed={args.seed} cases={len(cases)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())

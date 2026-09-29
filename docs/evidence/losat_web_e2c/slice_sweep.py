#!/usr/bin/env python3
"""Compare BLASTN on slices of the repository's genomes with NCBI (Session S07+; comparison only).

The seventh S07+ audit found BLASTN differences on slices of real genomes that the
constructed inputs did not show (the preliminary `-subject_besthit` filter, a crash of the
gap reduction at the start of a subject). This sweep cuts, from fixed seeds, a query of 2
to 20 kb and a subject of 5 to 30 kb from two different genomes of `LOSAT/tests/fasta`,
runs each pair in NCBI and LOSAT with the given options and -outfmt 6, and prints each
pair whose stdout, stderr or exit status differs. Exits 1 when any does.

Usage: slice_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--options "..."] [--seed N]
                      [--cases N] [--jobs N]
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

GENOMES = ["LC738868", "LC738869", "LC738870", "LC738871", "LC738873", "LC738874", "LC738875",
           "PemoMJNVA", "PemoMJNVB", "MelaMJNV", "MejoMJNV", "PeseMJNV", "TrcuMJNV", "LvMJNV", "MjeNMV"]


def sequence(name: str) -> str:
    lines = (ENGINE / "tests/fasta" / f"{name}.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines[1:])


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
    seqs = {name: sequence(name) for name in GENOMES}
    rng = random.Random(args.seed)
    cases = []
    for index in range(args.cases):
        query_genome, subject_genome = rng.sample(GENOMES, 2)
        query_start = rng.randrange(0, len(seqs[query_genome]) - 30000)
        query_length = rng.randrange(2000, 20000)
        subject_start = rng.randrange(0, len(seqs[subject_genome]) - 40000)
        subject_length = rng.randrange(5000, 30000)
        query = args.work / f"q{args.seed}_{index}.fa"
        subject = args.work / f"s{args.seed}_{index}.fa"
        query.write_text(f">q{index}\n{seqs[query_genome][query_start:query_start + query_length]}\n")
        subject.write_text(f">s{index}\n{seqs[subject_genome][subject_start:subject_start + subject_length]}\n")
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

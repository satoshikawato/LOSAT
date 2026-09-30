#!/usr/bin/env python3
"""Sweep BLASTN word sizes and e-values against NCBI (Session S07+; comparison only).

The sixth S07+ audit found that LOSAT sized the two-hit diagonal table of a single query
from one strand, where NCBI sizes it from the query block of both strands
(`blast_engine.c:1002-1003`), so that the diagonals of the two strands shared entries and
HSPs were lost. Small word sizes and large e-values make many seeds and show it. This
sweep runs the outfmt 0 fixture pairs with both tasks, word sizes 4 to 28 and e-values 10
and 1e5, with -outfmt 6, and prints each combination whose stdout, stderr or exit status
differs from NCBI's. Exits 1 when any does.

Usage: word_size_sweep.py --bin-dir DIR --losat LOSAT [--jobs N]
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
from run_oracle import ENGINE  # noqa: E402

F = "tests/fasta/outfmt0"
PAIRS = [("strand_query", "strand_subject"), ("mask_query", "mask_subject"),
         ("longdef_query", "longdef_subject"), ("many_query", "many_subject"),
         ("multi_query", "multi_subject"), ("width_query", "width_subject"),
         ("lcase_minus_query", "lcase_minus_subject"), ("iupac_query", "multi_subject")]
TASKS = ("blastn", "megablast")
WORD_SIZES = (4, 5, 6, 7, 8, 11, 16, 28)
EVALUES = ("10", "1e5")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    cases = [["-query", f"{F}/{query}.fasta", "-subject", f"{F}/{subject}.fasta", "-task", task,
              "-word_size", str(word_size), "-evalue", evalue, "-outfmt", "6"]
             for (query, subject), task, word_size, evalue in itertools.product(PAIRS, TASKS, WORD_SIZES, EVALUES)]

    def run(argv: list[str]) -> tuple[list[str], bool, int, int]:
        ncbi = subprocess.run([str(args.bin_dir / "blastn"), *argv], cwd=ENGINE, capture_output=True)
        ours = subprocess.run([str(args.losat.resolve()), "blastn", *argv], cwd=ENGINE, capture_output=True)
        same = (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr)
        return argv, same, ncbi.stdout.count(b"\n"), ours.stdout.count(b"\n")

    differing = 0
    with ThreadPoolExecutor(args.jobs) as pool:
        for argv, same, ncbi_lines, our_lines in pool.map(run, cases):
            if not same:
                differing += 1
                print("DIFF\t" + " ".join(argv) + f"\tncbi_lines={ncbi_lines}\tlosat_lines={our_lines}")
    print(f"# cases={len(cases)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())

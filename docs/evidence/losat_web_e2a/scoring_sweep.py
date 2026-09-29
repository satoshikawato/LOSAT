#!/usr/bin/env python3
"""Sweep BLASTN scoring options against NCBI in outfmt 6 (comparison only).

The certified BLASTN profiles use the default scoring of each task (blastn: reward 2,
penalty -3, gaps 5/2; megablast: reward 1, penalty -2, gaps 0/0). This runs the
constructed multi_* input of the outfmt 0 fixtures with other reward/penalty/gap
values, in NCBI and LOSAT, and prints one line per combination: `same`, `DIFF`, or
the exit status of each program when either fails (NCBI rejects many combinations).

Usage: scoring_sweep.py --bin-dir DIR --losat LOSAT
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from run_oracle import ENGINE

INPUT = ["-query", "tests/fasta/outfmt0/multi_query.fasta", "-subject", "tests/fasta/outfmt0/multi_subject.fasta"]
SCORES = [(1, -1), (1, -2), (1, -3), (1, -4), (2, -3), (2, -5), (2, -7), (3, -4), (4, -5), (5, -4)]
GAPS = [[], ["-gapopen", "2", "-gapextend", "2"], ["-gapopen", "5", "-gapextend", "2"],
        ["-gapopen", "0", "-gapextend", "0"]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    args = parser.parse_args()
    print("task\treward\tpenalty\tgaps\tresult\tncbi_stderr")
    for task in ("blastn", "megablast"):
        for reward, penalty in SCORES:
            for gaps in GAPS:
                argv = [*INPUT, "-task", task, "-reward", str(reward), "-penalty", str(penalty), *gaps,
                        "-outfmt", "6"]
                ncbi = subprocess.run([str(args.bin_dir / "blastn"), *argv], cwd=ENGINE, capture_output=True)
                losat = subprocess.run([str(args.losat), "blastn", *argv], cwd=ENGINE, capture_output=True)
                if ncbi.returncode or losat.returncode:
                    result = f"ncbi_exit={ncbi.returncode} losat_exit={losat.returncode}"
                else:
                    result = "same" if ncbi.stdout == losat.stdout else "DIFF"
                message = ncbi.stderr.decode(errors="replace").strip().splitlines()
                message = next((line for line in message if "rror" in line), message[0] if message else "")
                print("\t".join([task, str(reward), str(penalty), " ".join(gaps[1::2]) or "default", result,
                                 message]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

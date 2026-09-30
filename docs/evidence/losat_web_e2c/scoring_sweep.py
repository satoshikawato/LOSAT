#!/usr/bin/env python3
"""Sweep BLASTN scoring options against NCBI (Session S07+; comparison only).

Runs constructed inputs with reward/penalty/gap values in NCBI and LOSAT, with -outfmt 6
by default, and prints one line per combination:

- same: both succeed with the same stdout and stderr;
- same-error: both fail with the same exit status, stdout and stderr (NCBI's rejection);
- losat-rejects: LOSAT fails with a message that names what it does not implement (an
  explicit rejection), where NCBI succeeds or fails otherwise;
- DIFF: anything else, with the first difference.

The pairs cover NCBI's 12 tables, pairs with a common divisor, pairs without a table and
a penalty of 0; the gap costs cover table rows, costs beyond the tables, 0/0 and invalid
values. Exits 1 when any line is DIFF.

Usage: scoring_sweep.py --bin-dir DIR --losat LOSAT [--outfmt 0|6|7] [--jobs N]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
from run_oracle import ENGINE  # noqa: E402

INPUTS = {
    "multi": ["-query", "tests/fasta/outfmt0/multi_query.fasta",
              "-subject", "tests/fasta/outfmt0/multi_subject.fasta"],
    "iupac": ["-query", "tests/fasta/outfmt0/iupac_query.fasta",
              "-subject", "tests/fasta/outfmt0/multi_subject.fasta"],
}
TABLES = [(1, -5), (1, -4), (2, -7), (1, -3), (2, -5), (1, -2), (2, -3), (3, -4), (1, -1), (3, -2),
          (4, -5), (5, -4)]
DIVISIBLE = [(2, -2), (2, -4), (2, -6), (4, -6), (3, -6), (4, -10)]
NO_TABLE = [(1, -6), (3, -5), (2, -1)]
SCORES = [*TABLES, *DIVISIBLE, *NO_TABLE, (1, 0)]
GAPS = [[], ["2", "2"], ["5", "2"], ["0", "0"], ["4", "4"], ["6", "4"], ["10", "10"], ["0", "2"],
        ["3", "0"], ["-1", "2"]]
# A LOSAT failure that NCBI does not have must say that LOSAT does not implement it.
REJECTION_MARKERS = ("which LOSAT does not reproduce", "not supported by LOSAT", "is not implemented")


def first_difference(left: bytes, right: bytes) -> str:
    for number, (a, b) in enumerate(zip(left.split(b"\n"), right.split(b"\n")), 1):
        if a != b:
            return f"line {number}: {a[:80]!r} / {b[:80]!r}"
    return f"{len(left)} / {len(right)} bytes"


def compare(bin_dir: Path, losat: Path, outfmt: str, case: tuple[str, str, int, int, list[str]]) -> list[str]:
    input_name, task, reward, penalty, gaps = case
    gap_args = ["-gapopen", gaps[0], "-gapextend", gaps[1]] if gaps else []
    argv = [*INPUTS[input_name], "-task", task, "-reward", str(reward), "-penalty", str(penalty), *gap_args,
            "-outfmt", outfmt]
    ncbi = subprocess.run([str(bin_dir / "blastn"), *argv], cwd=ENGINE, capture_output=True)
    ours = subprocess.run([str(losat.resolve()), "blastn", *argv], cwd=ENGINE, capture_output=True)
    same_streams = ncbi.stdout == ours.stdout and ncbi.stderr == ours.stderr
    if ncbi.returncode == 0 and ours.returncode == 0:
        result = "same" if same_streams else "DIFF " + (
            "stdout " + first_difference(ncbi.stdout, ours.stdout) if ncbi.stdout != ours.stdout
            else "stderr " + first_difference(ncbi.stderr, ours.stderr))
    elif ncbi.returncode and ncbi.returncode == ours.returncode and same_streams:
        result = "same-error"
    elif ours.returncode and any(marker.encode() in ours.stderr for marker in REJECTION_MARKERS):
        result = "losat-rejects"
    else:
        result = (f"DIFF exit {ncbi.returncode}/{ours.returncode}: "
                  + ("stdout " + first_difference(ncbi.stdout, ours.stdout) if ncbi.stdout != ours.stdout
                     else "stderr " + first_difference(ncbi.stderr, ours.stderr)))
    message = (ncbi.stderr if ncbi.returncode else ours.stderr).decode(errors="replace").strip().splitlines()
    message = next((line for line in message if "rror" in line), message[0] if message else "")
    return [input_name, task, str(reward), str(penalty), " ".join(gaps) or "default", result, message[:160]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--outfmt", default="6")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    cases = [(input_name, task, reward, penalty, gaps) for input_name in INPUTS for task in ("blastn", "megablast")
             for reward, penalty in SCORES for gaps in GAPS]
    with ThreadPoolExecutor(args.jobs) as pool:
        lines = list(pool.map(lambda case: compare(args.bin_dir, args.losat, args.outfmt, case), cases))
    print("input\ttask\treward\tpenalty\tgaps\tresult\tmessage")
    for line in lines:
        print("\t".join(line))
    counts: dict[str, int] = {}
    for line in lines:
        key = line[5].split()[0]
        counts[key] = counts.get(key, 0) + 1
    print("# " + " ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 1 if counts.get("DIFF") else 0


if __name__ == "__main__":
    sys.exit(main())

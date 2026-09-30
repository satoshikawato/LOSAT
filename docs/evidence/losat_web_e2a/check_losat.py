#!/usr/bin/env python3
"""Compare LOSAT with the frozen NCBI outfmt 0 fixtures (LOSAT/tests/outfmt0_manifest.tsv).

Every fixture runs from LOSAT/ as `LOSAT <program> <search_argv>` (the arguments of the
oracle run). Its stdout must equal LOSAT/tests/fixtures/outfmt0/<fixture_id>.out byte for
byte, and its stderr must equal <fixture_id>.err when the manifest records stderr
(otherwise LOSAT must write nothing to stderr). Prints one line per fixture, with the
first differing line of a difference, and exits 1 on any difference.

Usage: check_losat.py --losat LOSAT [--programs blastn,blastp,tblastn] [--threads N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
import subprocess

from run_oracle import ENGINE, read_manifest, search_argv

FIXTURES = ENGINE / "tests/fixtures/outfmt0"


def first_difference(expected: bytes, actual: bytes) -> str:
    for number, (left, right) in enumerate(zip(expected.split(b"\n"), actual.split(b"\n")), 1):
        if left != right:
            return f"line {number}: expected {left[:100]!r}, got {right[:100]!r}"
    return f"lengths differ: expected {len(expected)} bytes, got {len(actual)}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--programs", default="blastn,blastp,tblastn")
    parser.add_argument("--threads", type=int)
    args = parser.parse_args()
    programs = set(args.programs.split(","))
    differing = []
    rows = [row for row in read_manifest()[2] if row["program"] in programs]
    for row in rows:
        threads = ["-num_threads", str(args.threads)] if args.threads else []
        run = subprocess.run([str(args.losat.resolve()), row["program"], *search_argv(row), *threads],
                             cwd=ENGINE, capture_output=True)
        expected_err_path = FIXTURES / f"{row['fixture_id']}.err"
        expected_err = expected_err_path.read_bytes() if row["stderr_sha256"] else b""
        expected_out = (FIXTURES / f"{row['fixture_id']}.out").read_bytes()
        problems = []
        if run.returncode != 0:
            problems.append(f"exit {run.returncode}: {run.stderr.decode(errors='replace').strip()[:200]}")
        elif run.stdout != expected_out:
            problems.append("stdout " + first_difference(expected_out, run.stdout))
        if run.returncode == 0 and run.stderr != expected_err:
            problems.append("stderr " + first_difference(expected_err, run.stderr))
        if problems:
            differing.append(row["fixture_id"])
        print(f"{row['fixture_id']}\t{'same' if not problems else 'DIFF'}\t{'; '.join(problems)}")
    print(f"# fixtures={len(rows)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())

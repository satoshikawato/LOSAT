#!/usr/bin/env python3
"""Compare LOSAT with the frozen NCBI outfmt 0 and 7 fixtures (LOSAT/tests/outfmt0_manifest.tsv).

Every fixture runs from LOSAT/ as `LOSAT <program> <search_argv>` (the arguments of the
oracle run). Its stdout must equal LOSAT/tests/fixtures/outfmt0/<fixture_id>.out byte for
byte, and its stderr must equal <fixture_id>.err when the manifest records stderr
(otherwise LOSAT must write nothing to stderr). Prints one line per fixture, with the
first differing line of a difference, and exits 1 on any difference.

A row with the contract `approved_db_gencode_deviation` (AGENTS.md's approved TBLASTX
exception: LOSAT's local search applies a non-default -db_gencode to the subject) is
compared with <fixture_id>.db.out instead, NCBI's search with the subject as a database,
outside the lines that name the database (`normalize_database_lines`); its line says
`exception` and how many lines differ from NCBI's local search (<fixture_id>.out).

Usage: check_losat.py --losat LOSAT [--programs blastn,blastp,tblastn,tblastx] [--threads N]
       [--only PREFIX]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
import subprocess

from run_oracle import DEVIATION, ENGINE, read_manifest, search_argv

FIXTURES = ENGINE / "tests/fixtures/outfmt0"


def first_difference(expected: bytes, actual: bytes) -> str:
    for number, (left, right) in enumerate(zip(expected.split(b"\n"), actual.split(b"\n")), 1):
        if left != right:
            return f"line {number}: expected {left[:100]!r}, got {right[:100]!r}"
    return f"lengths differ: expected {len(expected)} bytes, got {len(actual)}"


def normalize_database_lines(report: bytes) -> bytes:
    """A report without the lines that name the subjects: the database title of the outfmt 0
    prolog (`Database: ...` up to the counts line) and epilog (`  Database: ...` through the
    `Posted date:` line), the `# Database:` line of outfmt 7, and the space after `>` in a
    subject heading (a local subject is shown as `> title`, a database subject as `>title`)."""
    lines = report.split(b"\n")
    kept: list[bytes] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith(b"Database: "):
            while index < len(lines) and not re.match(rb"^ +\S+ sequences; ", lines[index]):
                index += 1
            kept.append(b"<database>")
            kept.append(lines[index] if index < len(lines) else b"")
        elif line.startswith(b"  Database: "):
            while index < len(lines) and not lines[index].startswith(b"    Posted date:"):
                index += 1
            kept.append(b"<database>")
        elif line.startswith(b"# Database: "):
            kept.append(b"<database>")
        elif line.startswith(b">") and not line.startswith(b"> "):
            kept.append(b"> " + line[1:])
        else:
            kept.append(line)
        index += 1
    return b"\n".join(kept)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--programs", default="blastn,blastp,tblastn,tblastx")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--only", default="")
    args = parser.parse_args()
    programs = set(args.programs.split(","))
    differing = []
    rows = [row for row in read_manifest()[2]
            if row["program"] in programs and row["fixture_id"].startswith(args.only)]
    for row in rows:
        threads = ["-num_threads", str(args.threads)] if args.threads else []
        run = subprocess.run([str(args.losat.resolve()), row["program"], *search_argv(row), *threads],
                             cwd=ENGINE, capture_output=True)
        expected_err_path = FIXTURES / f"{row['fixture_id']}.err"
        expected_err = expected_err_path.read_bytes() if row["stderr_sha256"] else b""
        expected_out = (FIXTURES / f"{row['fixture_id']}.out").read_bytes()
        problems = []
        deviation = row.get("contract") == DEVIATION
        if deviation:
            database_out = (FIXTURES / f"{row['fixture_id']}.db.out").read_bytes()
        if run.returncode != 0:
            problems.append(f"exit {run.returncode}: {run.stderr.decode(errors='replace').strip()[:200]}")
        elif deviation:
            observed, expected = normalize_database_lines(run.stdout), normalize_database_lines(database_out)
            if observed != expected:
                problems.append("stdout (database oracle) " + first_difference(expected, observed))
        elif run.stdout != expected_out:
            problems.append("stdout " + first_difference(expected_out, run.stdout))
        if run.returncode == 0 and run.stderr != expected_err:
            problems.append("stderr " + first_difference(expected_err, run.stderr))
        if problems:
            differing.append(row["fixture_id"])
        status = "DIFF" if problems else "same"
        if deviation and not problems:
            local = sum(1 for a, b in zip(expected_out.split(b"\n"), run.stdout.split(b"\n")) if a != b)
            status = "exception"
            problems.append(f"equal to the database oracle; {local} lines differ from NCBI's local -subject search")
        print(f"{row['fixture_id']}\t{status}\t{'; '.join(problems)}")
    print(f"# fixtures={len(rows)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())

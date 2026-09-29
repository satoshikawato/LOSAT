#!/usr/bin/env python3
"""Run the NCBI oracle over LOSAT/tests/outfmt0_manifest.tsv (comparison only).

Every fixture runs from LOSAT/ as `<program> -query Q -subject S [-task T] EXTRA -outfmt 0`
with the output on stdout. The run writes DIR/<fixture_id>.out and DIR/<fixture_id>.err
(when stderr is not empty), and prints the binaries, then one line per fixture (wall time,
sizes, hashes, command); keep that as the run's log.

- default: compare the outputs with the manifest's hash columns; exit 1 on a difference.
- --freeze: write the hash and size columns of the manifest instead (Session S06 only;
  run with --out LOSAT/tests/fixtures/outfmt0 so the fixtures are the frozen bytes).

NCBI reads these environment variables and .ncbirc, and each one changes the report
(AUTHORITY.md §A.9), so the run refuses to start when one is set or present.

Usage: run_oracle.py --bin-dir DIR --out DIR [--freeze]
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ENGINE = REPO / "LOSAT"
MANIFEST = ENGINE / "tests/outfmt0_manifest.tsv"
HASHES = ("stdout_sha256", "stdout_bytes", "stderr_sha256")
# BATCH_SIZE and CHUNK_SIZE change the query batches (blast_input_aux.cpp:86-90,
# local_blast.cpp:59-62), and with them the report of an invalid query.
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE")


def read_manifest() -> tuple[list[str], list[str], list[dict[str, str]]]:
    lines = MANIFEST.read_text().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    header, *body = [line.split("\t") for line in lines if not line.startswith("#")]
    return comments, header, [dict(zip(header, row)) for row in body]


def search_argv(row: dict[str, str]) -> list[str]:
    """The search arguments after the program name (the same for NCBI and LOSAT)."""
    task = ["-task", row["task"]] if row["task"] else []
    return ["-query", row["query"], "-subject", row["subject"], *task, *shlex.split(row["extra_args"]),
            "-outfmt", "0"]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    found = [name for name in REPORT_ENV if name in os.environ]
    found += [str(path) for path in (ENGINE / ".ncbirc", Path.home() / ".ncbirc") if path.exists()]
    if "NCBI" in os.environ:
        found.append("NCBI")
    if found:
        raise SystemExit(f"report-changing NCBI configuration present: {found}")
    args.out.mkdir(parents=True, exist_ok=True)
    comments, header, rows = read_manifest()
    for program in sorted({row["program"] for row in rows}):
        binary = args.bin_dir / program
        version = subprocess.run([binary, "-version"], capture_output=True, text=True, check=True).stdout
        print(f"# {binary} sha256={digest(binary.read_bytes())} {' / '.join(version.splitlines())}")

    records, differing = [], []
    for row in rows:
        command = [str(args.bin_dir / row["program"]), *search_argv(row)]
        start = time.monotonic()
        run = subprocess.run(command, cwd=ENGINE, capture_output=True)
        wall = time.monotonic() - start
        if run.returncode != 0:
            raise SystemExit(f"{row['fixture_id']}: exit {run.returncode}\n{run.stderr.decode()}")
        (args.out / f"{row['fixture_id']}.out").write_bytes(run.stdout)
        err = args.out / f"{row['fixture_id']}.err"
        if run.stderr:
            err.write_bytes(run.stderr)
        elif err.exists():
            err.unlink()
        observed = {"stdout_sha256": digest(run.stdout), "stdout_bytes": str(len(run.stdout)),
                    "stderr_sha256": digest(run.stderr) if run.stderr else ""}
        if args.freeze:
            row.update(observed)
        elif any(row[key] != observed[key] for key in HASHES):
            differing.append(row["fixture_id"])
        shown = shlex.join([row["program"], *search_argv(row)])
        records.append([row["fixture_id"], f"{wall:.2f}", *(observed[key] for key in HASHES), shown])

    for record in [["fixture_id", "wall_s", *HASHES, "command"], *records]:
        print("\t".join(record))
    if args.freeze:
        MANIFEST.write_text("\n".join([*comments, "\t".join(header),
                                       *("\t".join(row[key] for key in header) for row in rows)]) + "\n")
    print(f"# fixtures={len(rows)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Frozen `blastp -db` fixtures: LOSAT against the bytes NCBI BLAST+ 2.17.0 wrote (comparison only).

LOSAT/tests/db_search_manifest.tsv lists the rows (columns: row_id, kind, db, query, extra_args,
outfmt, rc, stdout_sha256, stdout_bytes, stderr_sha256, covers). Every row runs from LOSAT/ as
`blastp -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT` (an empty `db` cell is the literal empty
argument `-db ""`; `extra_args` is split like a shell line; `outfmt` is one argument, so a custom
`7 qseqid sseqid ...` spec stays whole). The databases are under LOSAT/tests/seqdb/db and
LOSAT/tests/db_search (README.md there says how they were built).

- kind `ncbi`: the expectation is NCBI's stdout, stderr and exit status (`rc`) of a run with
  `-num_threads 1`, every byte, nothing normalised. The frozen bytes are
  LOSAT/tests/fixtures/db_search/<row_id>.out and <row_id>.err (only when stderr is not empty).
- kind `losat-rejection`: NCBI is not run. LOSAT must exit with `rc`, write nothing to stdout, and
  its stderr must contain the text that follows `expect: ` in the `covers` cell (up to the end of
  the cell or a ` | ` separator). The hash columns stay empty.

`freeze --ncbi-bin DIR [--only PREFIX]` runs NCBI for the `ncbi` rows, writes the frozen files and
rewrites the rc, hash and size columns of the manifest (rows outside --only keep theirs). It
refuses to run when ~/.ncbirc exists or when BLASTDB, NCBI, BATCH_SIZE, CHUNK_SIZE, OLD_FSC or
BL2SEQ_LEGACY is set, because each one changes NCBI's result.

`check --losat BIN [--threads N] [--only PREFIX] [--out TSV]` runs LOSAT with `-num_threads N`
(default 1) in a cleaned environment (no NCBI*, BLAST* or LOSAT_* variable), compares stdout,
stderr and exit status with the frozen files (and checks that each frozen file's sha256 is the
manifest's), prints one line per row (`same`, `rejection`, or `differs` with the first differing
line), writes the same lines to --out, and exits 0 only when every row is `same` or `rejection`.
The tool takes no lock; whoever runs it next to other NCBI runs holds the caller's own.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shlex
import subprocess
import sys
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[1]
TESTS = ENGINE / "tests"
MANIFEST = TESTS / "db_search_manifest.tsv"
FROZEN = TESTS / "fixtures" / "db_search"
COLUMNS = ["row_id", "kind", "db", "query", "extra_args", "outfmt", "rc", "stdout_sha256", "stdout_bytes",
           "stderr_sha256", "covers"]
NCBI_ENV = ("BLASTDB", "NCBI", "BATCH_SIZE", "CHUNK_SIZE", "OLD_FSC", "BL2SEQ_LEGACY")
TIMEOUT = 600


def read_manifest() -> tuple[list[str], list[dict[str, str]]]:
    comments, rows = [], []
    header_seen = False
    for line in MANIFEST.read_text().split("\n"):
        if line.startswith("#"):
            comments.append(line)
        elif line and not header_seen:
            if line.split("\t") != COLUMNS:
                raise SystemExit(f"{MANIFEST.name}: unexpected header {line!r}")
            header_seen = True
        elif line:
            cells = line.split("\t")
            if len(cells) != len(COLUMNS):
                raise SystemExit(f"{MANIFEST.name}: {len(cells)} cells in {line[:60]!r}")
            rows.append(dict(zip(COLUMNS, cells)))
    return comments, rows


def write_manifest(comments: list[str], rows: list[dict[str, str]]) -> None:
    lines = comments + ["\t".join(COLUMNS)] + ["\t".join(row[c] for c in COLUMNS) for row in rows]
    MANIFEST.write_text("\n".join(lines) + "\n")


def search_argv(row: dict[str, str]) -> list[str]:
    """The arguments after the program name (the same for NCBI and LOSAT)."""
    return ["-db", row["db"], "-query", row["query"], *shlex.split(row["extra_args"]), "-outfmt", row["outfmt"]]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def first_difference(expected: bytes, actual: bytes) -> str:
    for number, (left, right) in enumerate(zip(expected.split(b"\n"), actual.split(b"\n")), 1):
        if left != right:
            return f"line {number}: expected {left[:100]!r}, got {right[:100]!r}"
    return f"lengths differ: expected {len(expected)} bytes, got {len(actual)}"


def expected_text(row: dict[str, str]) -> str:
    """The text a `losat-rejection` row's stderr must contain."""
    _, marker, rest = row["covers"].partition("expect: ")
    if not marker:
        raise SystemExit(f"{row['row_id']}: the covers cell has no 'expect: ' text")
    return rest.split(" | ")[0]


def freeze(args: argparse.Namespace) -> int:
    if (Path.home() / ".ncbirc").exists():
        raise SystemExit("~/.ncbirc exists: it changes NCBI's result")
    if bad := [name for name in NCBI_ENV if name in os.environ]:
        raise SystemExit(f"environment variables that change NCBI's result are set: {bad}")
    comments, rows = read_manifest()
    FROZEN.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if row["kind"] != "ncbi" or not row["row_id"].startswith(args.only):
            continue
        run = subprocess.run([str(Path(args.ncbi_bin) / "blastp"), *search_argv(row), "-num_threads", "1"],
                             cwd=ENGINE, capture_output=True, timeout=TIMEOUT)
        (FROZEN / f"{row['row_id']}.out").write_bytes(run.stdout)
        err = FROZEN / f"{row['row_id']}.err"
        if run.stderr:
            err.write_bytes(run.stderr)
        else:
            err.unlink(missing_ok=True)
        row.update(rc=str(run.returncode), stdout_sha256=sha256(run.stdout), stdout_bytes=str(len(run.stdout)),
                   stderr_sha256=sha256(run.stderr) if run.stderr else "")
        print(f"{row['row_id']}\trc={run.returncode}\tstdout={len(run.stdout)}\tstderr={len(run.stderr)}")
    write_manifest(comments, rows)
    return 0


def check_row(row: dict[str, str], losat: Path, threads: int, env: dict[str, str]) -> tuple[str, str]:
    run = subprocess.run([str(losat), "blastp", *search_argv(row), "-num_threads", str(threads)],
                         cwd=ENGINE, capture_output=True, env=env, timeout=TIMEOUT)
    if row["kind"] == "losat-rejection":
        problems = []
        if run.returncode != int(row["rc"]):
            problems.append(f"exit {run.returncode}, expected {row['rc']}")
        if run.stdout:
            problems.append(f"stdout is not empty ({len(run.stdout)} bytes)")
        if expected_text(row).encode() not in run.stderr:
            problems.append(f"stderr lacks {expected_text(row)!r}: {run.stderr.decode(errors='replace').strip()[:150]!r}")
        return ("differs", "; ".join(problems)) if problems else ("rejection", "")
    out = (FROZEN / f"{row['row_id']}.out").read_bytes()
    err_path = FROZEN / f"{row['row_id']}.err"
    err = err_path.read_bytes() if row["stderr_sha256"] else b""
    problems = []
    if sha256(out) != row["stdout_sha256"] or len(out) != int(row["stdout_bytes"]):
        problems.append("frozen .out does not match the manifest")
    if (sha256(err) if err else "") != row["stderr_sha256"] or err_path.exists() != bool(row["stderr_sha256"]):
        problems.append("frozen .err does not match the manifest")
    if run.returncode != int(row["rc"]):
        problems.append(f"exit {run.returncode}, expected {row['rc']}: {run.stderr.decode(errors='replace').strip()[:150]!r}")
    if run.stdout != out:
        problems.append("stdout " + first_difference(out, run.stdout))
    if run.stderr != err:
        problems.append("stderr " + first_difference(err, run.stderr))
    return ("differs", "; ".join(problems)) if problems else ("same", "")


def check(args: argparse.Namespace) -> int:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("NCBI", "BLAST", "LOSAT_"))}
    losat = args.losat.resolve()
    lines, bad = [], 0
    for row in read_manifest()[1]:
        if not row["row_id"].startswith(args.only):
            continue
        status, detail = check_row(row, losat, args.threads, env)
        bad += status == "differs"
        lines.append(f"{row['row_id']}\t{status}\t{detail}")
    text = "".join(line + "\n" for line in lines)
    sys.stdout.write(text)
    if args.out:
        args.out.write_text(text)
    print(f"# rows={len(lines)} differing={bad} threads={args.threads}")
    return 1 if bad else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    fz = sub.add_parser("freeze", help="run NCBI and write the frozen bytes and the manifest columns")
    fz.add_argument("--ncbi-bin", required=True)
    fz.add_argument("--only", default="")
    ck = sub.add_parser("check", help="compare LOSAT with the frozen bytes")
    ck.add_argument("--losat", type=Path, required=True)
    ck.add_argument("--threads", type=int, default=1)
    ck.add_argument("--only", default="")
    ck.add_argument("--out", type=Path)
    args = parser.parse_args()
    return freeze(args) if args.command == "freeze" else check(args)


if __name__ == "__main__":
    sys.exit(main())

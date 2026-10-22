#!/usr/bin/env python3
"""Frozen `blastp -db`, `blastn -db`, `tblastn -db`, `blastx -db` and `tblastx -db` fixtures: LOSAT against the bytes NCBI BLAST+ 2.17.0 wrote (comparison only).

LOSAT/tests/db_search_manifest.tsv lists the rows (columns: row_id, program, kind, db, query, extra_args,
outfmt, rc, stdout_sha256, stdout_bytes, stderr_sha256, covers). Every row runs from LOSAT/ as
`PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT`, PROGRAM being the row's `program` cell
(`blastp`, `blastn`, `tblastn`, `blastx` or `tblastx`; NCBI's own binary of that name, or `LOSAT PROGRAM`) (an empty `db` cell is the literal empty
argument `-db ""`; `extra_args` is split like a shell line, its tokens `BLASTDB=DIR` are not arguments but the environment variable of both programs (the directory of a `taxdb.bti`, relative to LOSAT/); `outfmt` is one argument, so a custom
`7 qseqid sseqid ...` spec stays whole). The databases are under LOSAT/tests/seqdb/db and
LOSAT/tests/db_search (README.md there says how they were built).

- kind `ncbi`: the expectation is NCBI's stdout, stderr and exit status (`rc`) of a run with
  `-num_threads 1`, every byte, nothing normalised. The frozen bytes are
  LOSAT/tests/fixtures/db_search/<row_id>.out and <row_id>.err (only when stderr is not empty).
- kind `losat-rejection`: NCBI is not run. LOSAT must exit with `rc`, write nothing to stdout, and
  its stderr must contain the text that follows `expect: ` in the `covers` cell (up to the end of
  the cell or a ` | ` separator). The hash columns stay empty.

`freeze --ncbi-bin DIR [--program NAME] [--only PREFIX]` runs NCBI for the `ncbi` rows, writes the frozen files and
rewrites the rc, hash and size columns of the manifest (rows outside --program and --only keep theirs;
the default is every program). It runs NCBI in a cleaned environment (see `check`), because each
variable and file that is cleaned changes NCBI's result: BLASTDB, NCBI*, BLAST*, BATCH_SIZE, CHUNK_SIZE,
OLD_FSC, BL2SEQ_LEGACY, and `~/.ncbirc` (HOME is an empty directory).

`check --losat BIN [--threads N] [--program NAME] [--only PREFIX] [--out TSV]` runs LOSAT with `-num_threads N`
(default 1) in the same cleaned environment (no NCBI*, BLAST*, LOSAT_*, BATCH_SIZE, CHUNK_SIZE, OLD_FSC or
BL2SEQ_LEGACY variable; HOME and USERPROFILE are an empty directory, so that no `~/.ncbirc` names a
taxonomy database or a BLASTDB), compares stdout,
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
import tempfile
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[1]
TESTS = ENGINE / "tests"
MANIFEST = TESTS / "db_search_manifest.tsv"
FROZEN = TESTS / "fixtures" / "db_search"
COLUMNS = ["row_id", "program", "kind", "db", "query", "extra_args", "outfmt", "rc", "stdout_sha256", "stdout_bytes",
           "stderr_sha256", "covers"]
PROGRAMS = ("blastp", "blastn", "tblastn", "blastx", "tblastx")
NCBI_ENV = ("BATCH_SIZE", "CHUNK_SIZE", "OLD_FSC", "BL2SEQ_LEGACY")
TIMEOUT = 600
# The `NAME=VALUE` tokens of `extra_args` that are environment variables of both programs: BLASTDB
# (the directory of a `taxdb.bti`) and the two that cut a query file into batches.
ENV_NAMES = ("BLASTDB", "BATCH_SIZE", "CHUNK_SIZE")


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
            row = dict(zip(COLUMNS, cells))
            if row["program"] not in PROGRAMS:
                raise SystemExit(f"{MANIFEST.name}: unknown program {row['program']!r} in {row['row_id']}")
            rows.append(row)
    return comments, rows


def write_manifest(comments: list[str], rows: list[dict[str, str]]) -> None:
    lines = comments + ["\t".join(COLUMNS)] + ["\t".join(row[c] for c in COLUMNS) for row in rows]
    MANIFEST.write_text("\n".join(lines) + "\n")


def extra_tokens(row: dict[str, str]) -> list[str]:
    return shlex.split(row["extra_args"])


def is_env_token(token: str) -> bool:
    return token.partition("=")[0] in ENV_NAMES


def search_argv(row: dict[str, str]) -> list[str]:
    """The arguments after the program name (the same for NCBI and LOSAT)."""
    extra = [token for token in extra_tokens(row) if not is_env_token(token)]
    return ["-db", row["db"], "-query", row["query"], *extra, "-outfmt", row["outfmt"]]


def row_env(row: dict[str, str], base: dict[str, str]) -> dict[str, str]:
    """`base` plus the row's `NAME=VALUE` tokens (`BLASTDB=DIR`: the search path of the taxonomy
    database, DIR relative to LOSAT/, the directory both programs run in; `BATCH_SIZE`, `CHUNK_SIZE`)."""
    env = dict(base)
    for token in extra_tokens(row):
        if is_env_token(token):
            name, _, value = token.partition("=")
            env[name] = value
    return env


def clean_env(base: dict[str, str], home: str) -> dict[str, str]:
    """`base` without what changes a result (NCBI*, BLAST*, LOSAT_*, `NCBI_ENV`), HOME and USERPROFILE
    set to the empty directory `home` (`$HOME/.ncbirc` is read by both programs)."""
    env = {k: v for k, v in base.items()
           if not k.startswith(("NCBI", "BLAST", "LOSAT_")) and k not in NCBI_ENV}
    env["HOME"] = env["USERPROFILE"] = home
    return env


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def is_sam_outfmt(outfmt: str) -> bool:
    """Whether `outfmt` (the `-outfmt` argument, a spec included) asks for SAM (17)."""
    return outfmt.split()[:1] == ["17"]


_SAM_PG = b"\n@PG\tID:0\tVN:2.17.0+"


def strip_sam_cl(data: bytes) -> bytes:
    """The SAM document without the CL tag of its @PG line: the bytes from the `\tCL:` that follows
    `\n@PG\tID:0\tVN:2.17.0+` up to, not including, the first following `\tPN:blastn\n`.

    The one rule that the SAM (outfmt 17) rows of every fixture compare with, on the frozen NCBI bytes and on
    LOSAT's (DECISIONS 7 of the outfmt line): LOSAT writes its own command line (its program path, the `blastn`
    sub-command word and the arguments a harness adds, such as `-num_threads N`) where NCBI writes its argv.
    Every other byte must be equal. A document without the tag (no command line) is returned as it is."""
    head = data.find(_SAM_PG + b"\tCL:")
    if head < 0:
        return data
    start = head + len(_SAM_PG)
    end = data.find(b"\tPN:blastn\n", start + len(b"\tCL:"))
    return data if end < 0 else data[:start] + data[end:]


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


def selected(row: dict[str, str], args: argparse.Namespace) -> bool:
    return row["row_id"].startswith(args.only) and args.program in ("", row["program"])


def freeze(args: argparse.Namespace) -> int:
    comments, rows = read_manifest()
    FROZEN.mkdir(parents=True, exist_ok=True)
    home = tempfile.TemporaryDirectory(prefix="db_search_home_")
    env = clean_env(dict(os.environ), home.name)
    for row in rows:
        if row["kind"] != "ncbi" or not selected(row, args):
            continue
        run = subprocess.run([str(Path(args.ncbi_bin) / row["program"]), *search_argv(row), "-num_threads", "1"],
                             cwd=ENGINE, capture_output=True, timeout=TIMEOUT, env=row_env(row, env))
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
    run = subprocess.run([str(losat), row["program"], *search_argv(row), "-num_threads", str(threads)],
                         cwd=ENGINE, capture_output=True, env=row_env(row, env), timeout=TIMEOUT)
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
    # A SAM row differs in the command line of the @PG line only (`strip_sam_cl`).
    expected_out, observed_out = (out, run.stdout)
    if is_sam_outfmt(row["outfmt"]):
        expected_out, observed_out = strip_sam_cl(out), strip_sam_cl(run.stdout)
    if observed_out != expected_out:
        problems.append("stdout " + first_difference(expected_out, observed_out))
    if run.stderr != err:
        problems.append("stderr " + first_difference(err, run.stderr))
    return ("differs", "; ".join(problems)) if problems else ("same", "")


def check(args: argparse.Namespace) -> int:
    home = tempfile.TemporaryDirectory(prefix="db_search_home_")
    env = clean_env(dict(os.environ), home.name)
    losat = args.losat.resolve()
    lines, bad = [], 0
    for row in read_manifest()[1]:
        if not selected(row, args):
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
    fz.add_argument("--program", choices=PROGRAMS, default="", help="only the rows of this program (default: all)")
    fz.add_argument("--only", default="")
    ck = sub.add_parser("check", help="compare LOSAT with the frozen bytes")
    ck.add_argument("--losat", type=Path, required=True)
    ck.add_argument("--threads", type=int, default=1)
    ck.add_argument("--program", choices=PROGRAMS, default="", help="only the rows of this program (default: all)")
    ck.add_argument("--only", default="")
    ck.add_argument("--out", type=Path)
    args = parser.parse_args()
    return freeze(args) if args.command == "freeze" else check(args)


if __name__ == "__main__":
    sys.exit(main())

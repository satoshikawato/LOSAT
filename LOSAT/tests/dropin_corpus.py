#!/usr/bin/env python3
"""Frozen acceptance ledger of the drop-in line: the command lines that downstream tools run (Prokka, ABRicate,
AMRFinderPlus, mlst, BlobToolKit, Galaxy, nf-core ...), each as `losat <command>`, against what NCBI BLAST+ 2.17.0
wrote for the same line (comparison oracle only; `check` runs no NCBI program).

LOSAT/tests/dropin_corpus_manifest.tsv lists the rows (columns: row_id, stdin, check, rc, stdout_sha256, stdout_bytes,
stderr_empty, files, command, note). A row runs in its own temporary directory that holds links to the inputs in
LOSAT/tests/dropin_corpus/; `stdin` names the input it reads on standard input. The command lines are verbatim,
`-num_threads` included. `stderr_empty` is `yes` when NCBI wrote nothing on stderr: LOSAT must be quiet there too
(AMRFinderPlus treats any text on stderr as an error); the stderr bytes are not compared otherwise.

The cell `check` says what is compared besides the exit status:

- `bytes`: stdout and every file the run writes, against the frozen bytes LOSAT/tests/fixtures/dropin_corpus/<row_id>.out
  and <row_id>--<file name>.
- `bytes-dir`: the same after replacing the run directory in stdout by `<RUNDIR>` (`blastdbcmd -info` prints the path).
- `line1`: `-version`: the first line of stdout is NCBI's and a second line follows (it names LOSAT).
- `contains:A|B`: stdout holds each text (a tool greps `-help` for it); NCBI's bytes are not frozen.

A `makeblastdb` row compares its report without the lines that hold the time or the absolute path (`NO_REPORT_LINES`),
and the sha256 of each database file in the `files` cell, without the date of a `.nin`/`.pin` or the keys of a
`.njs`/`.pjs` that follow from it (`database_bytes`); `-logfile` holds the time and the path, so only its name counts.
No bytes of a database are frozen.

`freeze --ncbi-bin DIR [--only PREFIX]` runs NCBI, writes the frozen files and rewrites the rc, hash and size columns.
`check --losat BIN [--only PREFIX] [--out TSV]` prints one line per row (`same`, or `differs` with the first
difference), writes the same lines to --out and exits 0 only when every row is `same`. Both run in the environment of
db_search_fixtures.clean_env.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS))
from db_search_fixtures import clean_env  # noqa: E402

MANIFEST = TESTS / "dropin_corpus_manifest.tsv"
INPUTS = TESTS / "dropin_corpus"
FROZEN = TESTS / "fixtures" / "dropin_corpus"
COLUMNS = ["row_id", "stdin", "check", "rc", "stdout_sha256", "stdout_bytes", "stderr_empty", "files", "command", "note"]
NOT_INPUTS = ("README.md", ".gitattributes")
NO_REPORT_LINES = (b"Building a new DB", b"New DB name:", b"Adding sequences from")
TIMEOUT = 300


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_manifest() -> tuple[list[str], list[dict[str, str]]]:
    comments, rows = [], []
    for line in MANIFEST.read_text().split("\n"):
        if line.startswith("#"):
            comments.append(line)
        elif line and not line.startswith("row_id\t"):
            rows.append(dict(zip(COLUMNS, line.split("\t"), strict=True)))
    return comments, rows


def write_manifest(comments: list[str], rows: list[dict[str, str]]) -> None:
    lines = comments + ["\t".join(COLUMNS)] + ["\t".join(row[name] for name in COLUMNS) for row in rows]
    MANIFEST.write_text("\n".join(lines) + "\n")


def is_database_row(argv: list[str]) -> bool:
    return argv[0] == "makeblastdb" and not {"-version", "-help", "-h"} & set(argv)


def index_without_date(data: bytes) -> bytes:
    """A `.nin`/`.pin` without its date: version, type, (v5) volume number, title, (v5) LMDB name, date, then the rest.
    The strings carry a big-endian length; the date field also holds the padding that aligns what follows."""
    version = int.from_bytes(data[:4], "big")
    at = 8 + (4 if version >= 5 else 0)
    for _ in range(1 if version < 5 else 2):
        at += 4 + int.from_bytes(data[at:at + 4], "big")
    date_end = at + 4 + int.from_bytes(data[at:at + 4], "big")
    return data[:at] + data[date_end:]


def database_bytes(name: str, data: bytes) -> bytes:
    """What is compared of a database file: `.nin`/`.pin` without the date; `.njs`/`.pjs` (the JSON manifest) without
    the time and the two sizes that count the index, whose padding follows the length of the date."""
    if name.endswith((".nin", ".pin")):
        return index_without_date(data)
    if name.endswith((".njs", ".pjs")):
        document = json.loads(data)
        for key in ("last-updated", "bytes-total", "bytes-to-cache"):
            del document[key]
        return json.dumps(document, sort_keys=True).encode()
    return data


class Run:
    def __init__(self, rc, stdout: bytes, stderr: bytes, files: dict[str, bytes]):
        self.rc, self.stdout, self.stderr, self.files = rc, stdout, stderr, files


def run_row(program: Path, ncbi: bool, row: dict[str, str], env: dict[str, str]) -> Run:
    """`program` is LOSAT's binary, or (`ncbi`) the directory of NCBI's programs."""
    argv = shlex.split(row["command"])
    command = [str(program / argv[0])] + argv[1:] if ncbi else [str(program)] + argv
    with tempfile.TemporaryDirectory(prefix="dropin_corpus_") as tmp:
        cwd = Path(tmp).resolve()
        for path in INPUTS.iterdir():
            if path.name not in NOT_INPUTS:
                os.symlink(path, cwd / path.name)
        before = set(os.listdir(cwd))
        stdin = open(INPUTS / row["stdin"], "rb") if row["stdin"] else subprocess.DEVNULL
        try:
            done = subprocess.run(command, cwd=cwd, env=env, stdin=stdin, capture_output=True, timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            return Run("TIMEOUT", b"", b"", {})
        finally:
            if row["stdin"]:
                stdin.close()
        stdout = done.stdout
        if row["check"] == "bytes-dir":
            stdout = stdout.replace(str(cwd).encode(), b"<RUNDIR>")
        if is_database_row(argv):
            stdout = b"\n".join(line for line in stdout.split(b"\n") if not line.startswith(NO_REPORT_LINES))
        files = {name: (cwd / name).read_bytes() for name in sorted(set(os.listdir(cwd)) - before)
                 if (cwd / name).is_file()}
        return Run(done.returncode, stdout, done.stderr, files)


def files_cell(row: dict[str, str], run: Run) -> str:
    """`name:sha256,...`: a database file through `database_bytes`, a `-logfile` by name alone (`-`)."""
    database = is_database_row(shlex.split(row["command"]))
    cells = []
    for name, data in run.files.items():
        if database and name.endswith(".log"):
            cells.append(f"{name}:-")
        else:
            cells.append(f"{name}:{sha256(database_bytes(name, data) if database else data)}")
    return ",".join(cells)


def first_difference(want: bytes, got: bytes) -> str:
    for number, (a, b) in enumerate(zip(want.split(b"\n"), got.split(b"\n")), 1):
        if a != b:
            return f"line {number}: frozen {a[:120]!r}, LOSAT {b[:120]!r}"
    return f"frozen {len(want)} bytes, LOSAT {len(got)} bytes"


def compare(row: dict[str, str], run: Run) -> str:
    """The first difference, or an empty string."""
    if str(run.rc) != row["rc"]:
        return f"exit status {run.rc}, NCBI {row['rc']}: {run.stderr[:200]!r}"
    if row["stderr_empty"] == "yes" and run.stderr.strip():
        return f"stderr where NCBI is quiet: {run.stderr[:200]!r}"
    check = row["check"]
    if check.startswith("contains:"):
        missing = [text for text in check[len("contains:"):].split("|") if text.encode() not in run.stdout]
        return f"stdout lacks {missing}" if missing else ""
    want = (FROZEN / f"{row['row_id']}.out").read_bytes()
    if sha256(want) != row["stdout_sha256"]:
        return "the frozen stdout does not match the manifest"
    if check == "line1":
        lines = run.stdout.split(b"\n")
        if lines[0] != want.split(b"\n")[0] or run.stdout.count(b"\n") != 2:
            return first_difference(want, run.stdout)
        return ""
    if run.stdout != want:
        return "stdout " + first_difference(want, run.stdout)
    got = files_cell(row, run)
    if got == row["files"]:
        return ""
    want_cells, got_cells = dict(c.split(":") for c in row["files"].split(",") if c), \
        dict(c.split(":") for c in got.split(",") if c)
    if set(want_cells) != set(got_cells):
        return f"files {sorted(want_cells)} frozen, LOSAT wrote {sorted(got_cells)}"
    name = next(n for n in want_cells if want_cells[n] != got_cells[n])
    frozen = FROZEN / f"{row['row_id']}--{name}"
    if frozen.exists():
        return f"{name} " + first_difference(frozen.read_bytes(), run.files[name])
    return f"{name}: sha256 {got_cells[name]}, frozen {want_cells[name]}"


def freeze(args) -> int:
    comments, rows = read_manifest()
    FROZEN.mkdir(parents=True, exist_ok=True)
    home = tempfile.TemporaryDirectory(prefix="dropin_corpus_home_")
    env = clean_env(dict(os.environ), home.name)
    failed = 0
    for row in rows:
        if not row["row_id"].startswith(args.only):
            continue
        argv = shlex.split(row["command"])
        run = run_row(args.ncbi_bin, True, row, env)
        if run.rc != 0:
            print(f"{row['row_id']}\tNCBI exit status {run.rc}: {run.stderr[:200]!r}")
            failed += 1
        for old in FROZEN.glob(f"{row['row_id']}.out"):
            old.unlink()
        for old in FROZEN.glob(f"{row['row_id']}--*"):
            old.unlink()
        row["rc"] = str(run.rc)
        row["stderr_empty"] = "no" if run.stderr.strip() else "yes"
        row["files"] = files_cell(row, run)
        if row["check"].startswith("contains:"):
            row["stdout_sha256"], row["stdout_bytes"], row["files"] = "", "", ""
            continue
        (FROZEN / f"{row['row_id']}.out").write_bytes(run.stdout)
        row["stdout_sha256"], row["stdout_bytes"] = sha256(run.stdout), str(len(run.stdout))
        if not is_database_row(argv):
            for name, data in run.files.items():
                (FROZEN / f"{row['row_id']}--{name}").write_bytes(data)
    write_manifest(comments, rows)
    return 1 if failed else 0


def check(args) -> int:
    _, rows = read_manifest()
    losat = args.losat.resolve()
    home = tempfile.TemporaryDirectory(prefix="dropin_corpus_home_")
    env = clean_env(dict(os.environ), home.name)
    lines, failed = [], 0
    for row in rows:
        if not row["row_id"].startswith(args.only):
            continue
        difference = compare(row, run_row(losat, False, row, env))
        failed += bool(difference)
        lines.append(f"{row['row_id']}\tdiffers\t{difference}" if difference else f"{row['row_id']}\tsame")
        print(lines[-1], flush=True)
    if args.out:
        args.out.write_text("\n".join(lines) + "\n")
    print(f"{len(lines) - failed} / {len(lines)} same")
    return 1 if failed or not lines else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("freeze")
    one.add_argument("--ncbi-bin", type=Path, required=True)
    one.add_argument("--only", default="")
    one.set_defaults(handler=freeze)
    two = commands.add_parser("check")
    two.add_argument("--losat", type=Path, required=True)
    two.add_argument("--only", default="")
    two.add_argument("--out", type=Path)
    two.set_defaults(handler=check)
    args = parser.parse_args()
    return args.handler(args)


if __name__ == "__main__":
    sys.exit(main())

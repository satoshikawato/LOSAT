#!/usr/bin/env python3
"""Frozen acceptance ledger of the drop-in line: the command lines that downstream tools run (Prokka, ABRicate,
AMRFinderPlus, mlst, BlobToolKit, Galaxy, nf-core ...) and the common forms of each program, each as `losat <command>`,
against what NCBI BLAST+ 2.17.0 wrote for the same line (comparison oracle only; `check` runs no NCBI program).

LOSAT/tests/dropin_corpus_manifest.tsv lists the rows (columns: row_id, stdin, check, rc, stdout_sha256, stdout_bytes,
stderr_empty, files, command, note). A row runs in its own temporary directory that holds links to the inputs in
LOSAT/tests/dropin_corpus/; `stdin` names the input it reads on standard input. The command lines are verbatim,
`-num_threads` included. `stderr_empty` is `yes` when NCBI wrote nothing on stderr: LOSAT must be quiet there too
(AMRFinderPlus treats any text on stderr as an error); the stderr bytes are not compared otherwise. A run must not
change an input (the rows reach them through links).

The cell `check` says what is compared besides the exit status:

- `bytes`: stdout and every file the run writes, against the frozen bytes LOSAT/tests/fixtures/dropin_corpus/<row_id>.out
  and <row_id>--<file name>.
- `bytes-dir`: the same after replacing the run directory in stdout by `<RUNDIR>` (`blastdbcmd -info` prints the path).
- `line1`: `-version`: the first line of stdout is NCBI's, and one more line names LOSAT.
- `contains:A|B`: stdout holds each text (a tool greps `-help` for it); NCBI's bytes are not frozen.

A `makeblastdb` row compares its report and its `-logfile` without the lines that hold the time or the absolute path
(`NO_REPORT_LINES`), and the sha256 of each database file in the `files` cell, without the date of a `.nin`/`.pin` or
the time of a `.njs`/`.pjs`, whose two sizes are compared with the files they count (`database_bytes`). A date of
another form leaves the file whole, so it differs from the frozen sha256. No bytes of a database are frozen.

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
# NCBI's CTime "b D, Y  H:m P", then the padding that aligns what follows the date.
INDEX_DATE = re.compile(rb"[A-Z][a-z]{2} \d{1,2}, \d{4}  \d{1,2}:\d{2} [AP]M\0*")
MANIFEST_TIME = re.compile(rb'("last-updated": ")\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(")')
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
    The strings carry a big-endian length; the date field also holds the padding that aligns what follows. A date of
    another form or a lost alignment keeps the file whole."""
    version = int.from_bytes(data[:4], "big")
    at = 8 + (4 if version >= 5 else 0)
    for _ in range(1 if version < 5 else 2):
        at += 4 + int.from_bytes(data[at:at + 4], "big")
    date_end = at + 4 + int.from_bytes(data[at:at + 4], "big")
    if date_end % 8 != 0 or INDEX_DATE.fullmatch(data[at + 4:date_end]) is None:
        return data
    return data[:at] + data[date_end:]


def manifest_without_time(data: bytes, files: dict[str, bytes]) -> bytes:
    """A `.njs`/`.pjs` without its time, and with each of its two sizes replaced by its distance from the files it
    counts: `bytes-total` all the files it lists, `bytes-to-cache` the index and sequence files. These sizes follow the
    length of the date (the padding of the index); the rest of the bytes stay as written."""
    try:
        listed = [str(name) for name in json.loads(data)["files"]]
    except (ValueError, KeyError, TypeError):
        return data

    def size(names) -> int:
        return sum(len(files.get(name, b"")) for name in names)

    sizes = {b"bytes-total": size(listed),
             b"bytes-to-cache": size(n for n in listed if n.endswith((".nin", ".pin", ".nsq", ".psq")))}
    data = MANIFEST_TIME.sub(rb"\1\2", data)
    for key, want in sizes.items():
        data = re.sub(rb'("' + key + rb'": )(\d+)', lambda m: m[1] + str(int(m[2]) - want).encode(), data)
    return data


def database_bytes(name: str, files: dict[str, bytes]) -> bytes:
    """What is compared of a file a `makeblastdb` row writes (module docstring)."""
    data = files[name]
    if name.endswith((".nin", ".pin")):
        return index_without_date(data)
    if name.endswith((".njs", ".pjs")):
        return manifest_without_time(data, files)
    if name.endswith(".log"):
        return report_without_time(data)
    return data


def report_without_time(data: bytes) -> bytes:
    return b"\n".join(line for line in data.split(b"\n") if not line.startswith(NO_REPORT_LINES))


class Run:
    def __init__(self, rc, stdout: bytes, stderr: bytes, files: dict[str, bytes], touched: list[str] | None = None):
        self.rc, self.stdout, self.stderr, self.files, self.touched = rc, stdout, stderr, files, touched or []


def input_state(path: Path) -> tuple[int, int]:
    state = path.stat()
    return state.st_mtime_ns, state.st_size


def run_row(program: Path, ncbi: bool, row: dict[str, str], env: dict[str, str]) -> Run:
    """`program` is LOSAT's binary, or (`ncbi`) the directory of NCBI's programs."""
    argv = shlex.split(row["command"])
    command = [str(program / argv[0])] + argv[1:] if ncbi else [str(program)] + argv
    with tempfile.TemporaryDirectory(prefix="dropin_corpus_") as tmp:
        cwd = Path(tmp).resolve()
        inputs = [path for path in INPUTS.iterdir() if path.name not in NOT_INPUTS]
        for path in inputs:
            os.symlink(path, cwd / path.name)
        states = [input_state(path) for path in inputs]
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
            stdout = report_without_time(stdout)
        files = {name: (cwd / name).read_bytes() for name in sorted(set(os.listdir(cwd)) - before)
                 if (cwd / name).is_file()}
        touched = [path.name for path, state in zip(inputs, states)
                   if input_state(path) != state or not (cwd / path.name).is_symlink()]
        return Run(done.returncode, stdout, done.stderr, files, touched)


def files_cell(row: dict[str, str], run: Run) -> str:
    """`name:sha256,...`, a file of a `makeblastdb` row through `database_bytes`."""
    database = is_database_row(shlex.split(row["command"]))
    return ",".join(f"{name}:{sha256(database_bytes(name, run.files) if database else data)}"
                    for name, data in run.files.items())


def first_difference(want: bytes, got: bytes) -> str:
    for number, (a, b) in enumerate(zip(want.split(b"\n"), got.split(b"\n")), 1):
        if a != b:
            return f"line {number}: frozen {a[:120]!r}, LOSAT {b[:120]!r}"
    return f"frozen {len(want)} bytes, LOSAT {len(got)} bytes"


def compare(row: dict[str, str], run: Run) -> str:
    """The first difference, or an empty string."""
    if run.touched:
        return f"the run changed the inputs {run.touched}"
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
        if lines[0] != want.split(b"\n")[0] or len(lines) != 3 or b"LOSAT" not in lines[1] or lines[2]:
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
    ncbi_bin = args.ncbi_bin.resolve()
    home = tempfile.TemporaryDirectory(prefix="dropin_corpus_home_")
    env = clean_env(dict(os.environ), home.name)
    failed = 0
    for row in rows:
        if not row["row_id"].startswith(args.only):
            continue
        argv = shlex.split(row["command"])
        run = run_row(ncbi_bin, True, row, env)
        if run.rc != 0 or run.touched:
            print(f"{row['row_id']}\tNCBI exit status {run.rc}, changed inputs {run.touched}: {run.stderr[:200]!r}")
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
    log = args.out.open("w") if args.out else None
    for row in rows:
        if not row["row_id"].startswith(args.only):
            continue
        try:
            difference = compare(row, run_row(losat, False, row, env))
        except (OSError, ValueError) as error:  # a frozen file that is gone, a file LOSAT wrote in another form
            difference = f"raised {error!r}"
        failed += bool(difference)
        lines.append(f"{row['row_id']}\tdiffers\t{difference}" if difference else f"{row['row_id']}\tsame")
        print(lines[-1], flush=True)
        if log:
            log.write(lines[-1] + "\n")
            log.flush()
    if log:
        log.close()
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

#!/usr/bin/env python3
"""Fixture rows of `-outfmt 13` and `-outfmt 14` (BLAST JSON / XML2 side files and a master file).

An `-outfmt 13` or `14` run needs `-out FILE`: it writes the report of the n-th query to
`<base>_<n>.json` / `<base>_<n>.xml` next to FILE and a master file at FILE that lists them
(DESIGN_OF5 1.7). The manifests `outfmt0_manifest.tsv` and `db_search_manifest.tsv` name such a row
with two more columns after the others (a row without them is a row of the old form):

- `files`: empty (the output is on standard output) or `xinclude` (a row of outfmt 13 or 14). The
  runner adds `-out <tmp>/<fixture_id>.out`, with a fresh `<tmp>`. The base of the side files is then
  `<fixture_id>.out` (its extension is not `.xml`/`.json`), so the side files are
  `<fixture_id>.out_<n>.xml|json` and the master lists those names, whatever `<tmp>` is.
- `side_sha256`: the SHA-256 of each side file in the order of the reports, joined with `,`.

The frozen bytes: the master as `<fixture_id>.out` (its hash and size are the row's `stdout_sha256`
and `stdout_bytes`, the stdout of the run being empty), the side files as `<fixture_id>.out_<n>.<ext>`.
A check compares the master, the number of side files and each side file, the empty standard output,
standard error and the exit status.

This module has no manifest of its own: `run_oracle.py` / `check_losat.py` (outfmt 0 manifest) and
`db_search_fixtures.py` call it.
"""
from __future__ import annotations

import hashlib
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

XINCLUDE = "xinclude"
FILES_COLUMNS = ["files", "side_sha256"]
EXTENSIONS = {"13": "json", "14": "xml"}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def wants_files(row: dict[str, str]) -> bool:
    """Whether the row is a run of outfmt 13 or 14 with side files (the `files` column)."""
    files = row.get("files", "") or ""
    if files not in ("", XINCLUDE):
        raise SystemExit(f"unknown files column {files!r} (empty or {XINCLUDE!r})")
    return files == XINCLUDE


def extension(outfmt: str) -> str:
    """The side files' extension of a row's `-outfmt` value (`13` -> json, `14` -> xml)."""
    try:
        return EXTENSIONS[outfmt.strip()]
    except KeyError:
        raise SystemExit(f"a files=xinclude row needs outfmt 13 or 14, not {outfmt!r}") from None


def side_name(row_id: str, n: int, outfmt: str) -> str:
    return f"{row_id}.out_{n}.{extension(outfmt)}"


@dataclass
class Result:
    """A run with side files: the process, the master file's bytes (None: the file does not exist),
    the side files in the order of their numbers, and the names that do not belong (a file of the
    temporary directory that is neither the master nor `<fixture_id>.out_<n>.<ext>`)."""

    run: subprocess.CompletedProcess
    master: bytes | None
    sides: list[bytes]
    stray: list[str] = field(default_factory=list)


def collect(directory: Path, row_id: str, outfmt: str) -> tuple[bytes | None, list[bytes], list[str]]:
    """The master file and the numbered side files that a run left in `directory`."""
    master = directory / f"{row_id}.out"
    pattern = re.compile(re.escape(f"{row_id}.out_") + r"([1-9][0-9]*)\." + re.escape(extension(outfmt)) + "$")
    numbered: dict[int, bytes] = {}
    stray = []
    for path in sorted(directory.iterdir()):
        if path == master:
            continue
        match = pattern.match(path.name)
        if match:
            numbered[int(match.group(1))] = path.read_bytes()
        else:
            stray.append(path.name)
    sides = [numbered[n] for n in sorted(numbered)]
    if sorted(numbered) != list(range(1, len(sides) + 1)):
        stray.append(f"side file numbers {sorted(numbered)} are not 1..{len(sides)}")
    return (master.read_bytes() if master.exists() else None), sides, stray


def run_with_files(command: list[str], row_id: str, outfmt: str, **kwargs) -> Result:
    """Runs `command` (the search arguments, `-outfmt` included) with `-out <tmp>/<row_id>.out` added
    and collects what it wrote. `kwargs` go to `subprocess.run` (cwd, env, timeout); the output is
    captured."""
    with tempfile.TemporaryDirectory(prefix="losat_xinclude_") as tmp:
        directory = Path(tmp)
        run = subprocess.run([*command, "-out", str(directory / f"{row_id}.out")], capture_output=True, **kwargs)
        master, sides, stray = collect(directory, row_id, outfmt)
    return Result(run, master, sides, stray)


def observed(result: Result) -> dict[str, str]:
    """The hash columns of the row's files: the master's hash and size, and the side files' hashes."""
    master = result.master or b""
    return {"stdout_sha256": sha256(master), "stdout_bytes": str(len(master)),
            "side_sha256": ",".join(sha256(side) for side in result.sides)}


def freeze(result: Result, frozen: Path, row_id: str, outfmt: str) -> None:
    """Writes the master as `<row_id>.out` and the side files; removes the side files of an older freeze."""
    if result.master is None:
        raise SystemExit(f"{row_id}: the run wrote no master file")
    if result.run.stdout:
        raise SystemExit(f"{row_id}: standard output is not empty ({len(result.run.stdout)} bytes)")
    if result.stray:
        raise SystemExit(f"{row_id}: unexpected files {result.stray}")
    for old in frozen.glob(f"{row_id}.out_*"):
        old.unlink()
    (frozen / f"{row_id}.out").write_bytes(result.master)
    for n, side in enumerate(result.sides, 1):
        (frozen / side_name(row_id, n, outfmt)).write_bytes(side)


def problems(result: Result, row: dict[str, str], frozen: Path, first_difference) -> list[str]:
    """What differs between a run and the frozen files of `row` (empty: the same). `first_difference` is
    the caller's `(expected, actual) -> text`. The row's hash columns must also agree with the frozen
    files."""
    row_id = row.get("fixture_id") or row["row_id"]
    outfmt = row.get("outfmt", "")
    out: list[str] = []
    expected = (frozen / f"{row_id}.out").read_bytes()
    names = sorted(frozen.glob(f"{row_id}.out_*"),
                   key=lambda path: int(path.name.removeprefix(f"{row_id}.out_").split(".")[0]))
    expected_sides = [path.read_bytes() for path in names]
    if sha256(expected) != row.get("stdout_sha256") or len(expected) != int(row.get("stdout_bytes") or -1):
        out.append("frozen master does not match the manifest")
    if ",".join(sha256(side) for side in expected_sides) != row.get("side_sha256", ""):
        out.append("frozen side files do not match side_sha256")
    if [path.name for path in names] != [side_name(row_id, n, outfmt) for n in range(1, len(names) + 1)]:
        out.append(f"frozen side files are named {[path.name for path in names]}")
    if result.run.stdout:
        out.append(f"stdout is not empty ({len(result.run.stdout)} bytes)")
    if result.master is None:
        out.append("no master file")
    elif result.master != expected:
        out.append("master " + first_difference(expected, result.master))
    if len(result.sides) != len(expected_sides):
        out.append(f"{len(result.sides)} side files, expected {len(expected_sides)}")
    else:
        for n, (want, got) in enumerate(zip(expected_sides, result.sides), 1):
            if want != got:
                out.append(f"side file {n} " + first_difference(want, got))
    if result.stray:
        out.append(f"unexpected files {result.stray}")
    return out

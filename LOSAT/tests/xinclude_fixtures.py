#!/usr/bin/env python3
"""Fixture rows that write files: `-outfmt 13` and `-outfmt 14` (BLAST JSON / XML2 side files and a
master file), and psiblast's ASCII PSSM (`-out_ascii_pssm`).

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

A row whose `files` column is `pssm` (psiblast's `-out_ascii_pssm`, line PS, PS2) runs with
`-out_ascii_pssm <tmp>/<fixture_id>.pssm` added, with a fresh `<tmp>`. The run writes that file, or
with `-save_each_pssm` the files `<fixture_id>.pssm.<n>`, or none. Its standard output, standard
error and exit status are the row's as for a row without files. The files are frozen under their
names, and `side_sha256` lists them as `<name without the fixture id>=<SHA-256>` (`.pssm=...`,
`.pssm.1=...`), the file `.pssm` first and the numbered ones in order; an empty `side_sha256` is the
expectation that the run writes no file. A check compares the names and the bytes of the files. An
argument of the row may name that `<tmp>` as `{pssm_dir}` (`-out_pssm {pssm_dir}/p.asn`), so that any
other file the run writes there is seen: every name in `<tmp>` other than the PSSM files is reported.

The kinds `pssm-asn` and `pssm-both` (line PS, PS3: psiblast's ASN.1 PSSM, `-out_pssm`) are the same
with `-out_pssm <tmp>/<fixture_id>.asn` added instead of the ASCII file, or with both added. Their
files are `<fixture_id>.asn` (`.asn.<n>`) and `<fixture_id>.pssm` (`.pssm.<n>`), listed in
`side_sha256` in that order. In these two kinds a file is frozen (committed) only when it is 64 KB or
less (`FREEZE_LIMIT`); a larger one is compared by its SHA-256 in `side_sha256` alone, which is a
comparison of all of its bytes, without a first differing line in the report (the ASN.1 PSSM of a
400-letter query is 441 KB, 729 KB with the diagnostics of `-out_ascii_pssm`). The `pssm` kind
freezes every file, as before.

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
PSSM = "pssm"
PSSM_ASN = "pssm-asn"
PSSM_BOTH = "pssm-both"
PSSM_KINDS = (PSSM, PSSM_ASN, PSSM_BOTH)
# The largest PSSM file of the kinds pssm-asn and pssm-both that is frozen (bytes).
FREEZE_LIMIT = 64 * 1024
FILES_COLUMNS = ["files", "side_sha256"]
EXTENSIONS = {"13": "json", "14": "xml"}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def files_kind(row: dict[str, str]) -> str:
    """The row's `files` column: empty, `xinclude`, `pssm`, `pssm-asn` or `pssm-both`."""
    files = row.get("files", "") or ""
    if files not in ("", XINCLUDE, *PSSM_KINDS):
        raise SystemExit(f"unknown files column {files!r} (empty, {XINCLUDE!r} or one of {PSSM_KINDS})")
    return files


def wants_files(row: dict[str, str]) -> bool:
    """Whether the row is a run of outfmt 13 or 14 with side files (the `files` column)."""
    return files_kind(row) == XINCLUDE


def wants_pssm(row: dict[str, str]) -> bool:
    """Whether the row is a psiblast run whose PSSM files are compared (the `files` column)."""
    return files_kind(row) in PSSM_KINDS


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


def defects(result: Result) -> list[str]:
    """What is wrong with a run before any byte is compared: no master file, standard output, or a file
    of the temporary directory that is neither the master nor a numbered side file. A freeze refuses to
    write such a run; a comparison counts it as a difference of that fixture and goes on."""
    found = []
    if result.master is None:
        found.append("the run wrote no master file")
    if result.run.stdout:
        found.append(f"standard output is not empty ({len(result.run.stdout)} bytes)")
    if result.stray:
        found.append(f"unexpected files {result.stray}")
    return found


def freeze(result: Result, frozen: Path, row_id: str, outfmt: str) -> None:
    """Writes the master as `<row_id>.out` and the side files; removes the side files of an older freeze."""
    if found := defects(result):
        raise SystemExit(f"{row_id}: {'; '.join(found)}")
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


# --- psiblast's PSSM files (`files` = pssm, pssm-asn, pssm-both) ----------------------------------


@dataclass
class PssmResult:
    """A run with `-out_ascii_pssm` and/or `-out_pssm`: the process, the files it wrote as (name, bytes)
    in the order of `side_sha256`, and the names that do not belong."""

    run: subprocess.CompletedProcess
    files: list[tuple[str, bytes]]
    stray: list[str] = field(default_factory=list)


def pssm_files(directory: Path, row_id: str) -> tuple[list[tuple[str, bytes]], list[str]]:
    """The PSSM files of `row_id` in `directory` (`<row_id>.asn`, `<row_id>.asn.<n>` by n, then
    `<row_id>.pssm`, `<row_id>.pssm.<n>`) and the other names."""
    found: list[tuple[int, int, str]] = []
    stray = []
    for path in directory.iterdir():
        for group, extension in enumerate(("asn", "pssm")):
            base = f"{row_id}.{extension}"
            match = re.match(re.escape(base) + r"\.([1-9][0-9]*)$", path.name)
            if path.name == base:
                found.append((group, 0, path.name))
                break
            if match:
                found.append((group, int(match.group(1)), path.name))
                break
        else:
            stray.append(path.name)
    return [(name, (directory / name).read_bytes()) for _, _, name in sorted(found)], sorted(stray)


PSSM_DIR = "{pssm_dir}"


def pssm_options(kind: str, directory: Path, row_id: str) -> list[str]:
    """The PSSM file options that a row of `kind` gets."""
    asn = ["-out_pssm", str(directory / f"{row_id}.asn")]
    ascii_pssm = ["-out_ascii_pssm", str(directory / f"{row_id}.pssm")]
    return {PSSM: ascii_pssm, PSSM_ASN: asn, PSSM_BOTH: asn + ascii_pssm}[kind]


def run_with_pssm(command: list[str], row_id: str, kind: str = PSSM, **kwargs) -> PssmResult:
    """Runs `command` with the PSSM file options of `kind` added (`<tmp>/<row_id>.pssm`, `.asn`) and
    collects the files it wrote; `{pssm_dir}` in an argument is `<tmp>`. `kwargs` go to
    `subprocess.run` (cwd, env, timeout); the output is captured."""
    with tempfile.TemporaryDirectory(prefix="losat_pssm_") as tmp:
        directory = Path(tmp)
        command = [word.replace(PSSM_DIR, tmp) for word in command]
        run = subprocess.run([*command, *pssm_options(kind, directory, row_id)],
                             capture_output=True, **kwargs)
        files, stray = pssm_files(directory, row_id)
    return PssmResult(run, files, stray)


def pssm_cell(files: list[tuple[str, bytes]], row_id: str) -> str:
    """The `side_sha256` cell of PSSM files."""
    return ",".join(f"{name.removeprefix(row_id)}={sha256(data)}" for name, data in files)


def observed_pssm(result: PssmResult, row_id: str) -> dict[str, str]:
    """The `side_sha256` column of a PSSM run."""
    return {"side_sha256": pssm_cell(result.files, row_id)}


def frozen(kind: str, data: bytes) -> bool:
    """Whether a PSSM file of a row of `kind` is frozen (`FREEZE_LIMIT`)."""
    return kind == PSSM or len(data) <= FREEZE_LIMIT


def freeze_pssm(result: PssmResult, frozen_dir: Path, row_id: str, kind: str = PSSM) -> None:
    """Writes the PSSM files that are frozen under their names; removes those of an older freeze."""
    if result.stray:
        raise SystemExit(f"{row_id}: unexpected files {result.stray}")
    old, _ = pssm_files(frozen_dir, row_id) if frozen_dir.is_dir() else ([], [])
    for name, _ in old:
        (frozen_dir / name).unlink()
    for name, data in result.files:
        if frozen(kind, data):
            (frozen_dir / name).write_bytes(data)


def parse_pssm_cell(cell: str, row_id: str) -> list[tuple[str, str]]:
    """The (file name, SHA-256) pairs of a `side_sha256` cell."""
    pairs = []
    for item in filter(None, cell.split(",")):
        suffix, _, digest = item.partition("=")
        pairs.append((row_id + suffix, digest))
    return pairs


def frozen_pssm_agrees(row: dict[str, str], frozen_dir: Path) -> bool:
    """Whether the frozen PSSM files of `row` agree with its `side_sha256`: each one is listed with its
    SHA-256, a file of the `pssm` kind is frozen whenever it is listed, and a frozen file of the
    other kinds is not larger than `FREEZE_LIMIT`."""
    row_id = row.get("fixture_id") or row["row_id"]
    kind = files_kind(row)
    digests = dict(parse_pssm_cell(row.get("side_sha256", ""), row_id))
    stored, _ = pssm_files(frozen_dir, row_id) if frozen_dir.is_dir() else ([], [])
    if any(digests.get(name) != sha256(data) or not frozen(kind, data) for name, data in stored):
        return False
    return kind != PSSM or len(stored) == len(digests)


def pssm_problems(result: PssmResult, row: dict[str, str], frozen_dir: Path, first_difference) -> list[str]:
    """What differs between a run's PSSM files and those of `row` (empty: the same): the names and
    SHA-256 of `side_sha256`, with the first differing line where the file is frozen; and whether the
    frozen files agree with `side_sha256` (every file of the `pssm` kind is frozen)."""
    row_id = row.get("fixture_id") or row["row_id"]
    expected = parse_pssm_cell(row.get("side_sha256", ""), row_id)
    digests = dict(expected)
    stored, _ = pssm_files(frozen_dir, row_id) if frozen_dir.is_dir() else ([], [])
    stored_bytes = dict(stored)
    out = []
    if not frozen_pssm_agrees(row, frozen_dir):
        out.append("frozen PSSM files do not match side_sha256")
    names = [name for name, _ in result.files]
    if names != [name for name, _ in expected]:
        out.append(f"PSSM files {names}, expected {[name for name, _ in expected]}")
    else:
        for name, got in result.files:
            if sha256(got) == digests[name]:
                continue
            if name in stored_bytes:
                out.append(f"{name} " + first_difference(stored_bytes[name], got))
            else:
                out.append(f"{name} differs (sha256 {sha256(got)}, expected {digests[name]}; not frozen)")
    if result.stray:
        out.append(f"unexpected files {result.stray}")
    return out

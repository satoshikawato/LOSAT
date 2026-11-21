#!/usr/bin/env python3
"""Frozen `blastp -db`, `blastn -db`, `tblastn -db`, `blastx -db`, `tblastx -db` and `psiblast -db` fixtures: LOSAT against the bytes NCBI BLAST+ 2.17.0 wrote (comparison only).

LOSAT/tests/db_search_manifest/<program>.tsv (one file per program, and <program>.<family>.tsv for a row family
with its own file, `family_files`; read by `read_manifest`, which returns the rows in the order of `PROGRAMS`, each
program's file before its family files, then the order of each file; README.md there) lists the rows (columns: row_id, program, kind, db, query, extra_args,
outfmt, rc, stdout_sha256, stdout_bytes, stderr_sha256, covers, and for the rows of outfmt 13 and 14 two more:
files, side_sha256; xinclude_fixtures.py says what they mean: the run gets `-out <tmp>/<row_id>.out`, the frozen
master is <row_id>.out, the side files <row_id>.out_<n>.json|xml; a row without them has 12 cells; a psiblast row
whose files column is `pssm` gets `-out_ascii_pssm <tmp>/<row_id>.pssm`, and the files it writes are frozen as
<row_id>.pssm and <row_id>.pssm.<n>, listed in side_sha256, none when the cell is empty; `{pssm_dir}` in its
extra_args is that <tmp>, where any other file the run writes is reported; `pssm-asn` gets `-out_pssm
<tmp>/<row_id>.asn` instead and `pssm-both` both, and in these two kinds only a file of 64 KB or less is frozen,
a larger one being compared by its SHA-256 in side_sha256). Every row runs from LOSAT/ as
`PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT`, PROGRAM being the row's `program` cell
(a psiblast or tblastn query cell `pssm-gz:PATH` or `pssm:PATH` is a PSSM query: `-in_pssm FILE` in place of
`-query QUERY`, the gzipped file PATH inflated into a temporary directory, or PATH as it is;
xinclude_fixtures.py `query_argv`; a query cell `stdin:pipe:PATH`, `stdin:file:PATH` or `stdin:empty` is
`-query -` with the file PATH's bytes through a pipe, the file PATH itself, or an empty pipe as the run's
standard input, `query_stdin`)
(`blastp`, `blastn`, `tblastn`, `blastx`, `tblastx` or `psiblast`; NCBI's own binary of that name, or `LOSAT PROGRAM`) (an empty `db` cell is the literal empty
argument `-db ""`, and a cell `subject:PATH` runs `-subject PATH` in place of `-db DB`: psiblast's refusal of a
subject set whose sequences all have no letters, which only a subject set reaches; `extra_args` is split like a shell line, its tokens `BLASTDB=DIR` are not arguments but the environment variable of both programs (the directory of a `taxdb.bti`, relative to LOSAT/); `outfmt` is one argument, so a custom
`7 qseqid sseqid ...` spec stays whole). The databases are under LOSAT/tests/seqdb/db and
LOSAT/tests/db_search (README.md there says how they were built).

- kind `ncbi`: the expectation is NCBI's stdout, stderr and exit status (`rc`) of a run with
  `-num_threads 1`, every byte, nothing normalised. The frozen bytes are
  LOSAT/tests/fixtures/db_search/<row_id>.out and <row_id>.err (only when stderr is not empty).
- a `ncbi` row whose `covers` cell has `ncbi-near: ROW_ID` (up to a space or a ` | `): the expectation is the frozen
  bytes, exit status and hash columns of the `ncbi` row ROW_ID, NCBI's run of a near input (the same search with
  the mask one letter shorter), where NCBI's own run of the row reaches an approved NCBI defect
  (PD-LOSAT-NCBI-DEFECTS exception 13: `-db_hard_mask` with a mask that runs to the end, on which NCBI BLAST+
  2.17.0 crashes, adds rows from the masked letters or drops seeds of the next subject; exception 3: a one-hit gapped
  start whose window reads past the end of a database volume's last memory page, where NCBI faults). The row has no
  frozen files of its own; `freeze` does not run NCBI for it and copies ROW_ID's columns (`near_row_id`).
- a `ncbi` row whose `covers` cell has a part ` | ncbi-valgrind`: the frozen bytes are NCBI's run under valgrind
  memcheck, where NCBI's native run reads past its buffers (PD-LOSAT-NCBI-DEFECTS exception 3: the window of a one-hit
  gapped start of negative length; natively NCBI may crash there or read bytes that differ from run to run, valgrind
  reports the reads and goes on). `freeze` runs such a row only with `--valgrind CMD` (the command line put before
  NCBI's, for example `valgrind -q --log-file=/dev/null`, so that NCBI's stderr stays its own; a stderr with
  valgrind's `==PID==` lines stops the freeze) and otherwise leaves its columns and files as they are
  (`under_valgrind`).
- a `ncbi` row whose `covers` cell has `ncbi-stack: unlimited`: `freeze` runs NCBI with no stack limit
  (`RLIMIT_STACK` infinity). NCBI's reader of a PSSM recurses once per nested value on its main thread's stack,
  and with the default 8 MiB a PSSM whose values nest some thousands of levels deep crashes it (SIGSEGV) at a depth
  that the stack limit, the environment and the stack's random start move; its output with a stack large enough
  is the expectation (PD-LOSAT-NCBI-DEFECTS exception 15). `check` runs LOSAT with the caller's limit.
- kind `losat-rejection`: LOSAT must exit with `rc`, write nothing to stdout, and its stderr must
  contain the text that follows `expect: ` in the `covers` cell (up to the end of the cell or a
  ` | ` separator). The hash columns stay empty. A rejection row with a PSSM files kind must also
  write no PSSM file. A rejection that LOSAT makes after reports that NCBI writes too (psiblast's
  `-out_pssm` file that cannot be written: NCBI stops after that query's report with its own error)
  has a non-empty `stdout_bytes` cell: its stdout must be NCBI's of the same command, frozen as
  <row_id>.out by `freeze` (which runs NCBI for such a row and keeps only its stdout) and checked
  against `stdout_sha256`. Its `rc` and stderr are LOSAT's. Where NCBI crashes on such a row, the frozen
  stdout is the part of the report that LOSAT writes before it refuses; `freeze` keeps it when NCBI's
  stdout stops earlier or later but agrees up to that point (`frozen_stdout`).

`freeze --ncbi-bin DIR [--program NAME] [--only PREFIX] [--valgrind CMD]` runs NCBI for the `ncbi` rows, writes the frozen files and
rewrites the rc, hash and size columns of the manifest (rows outside --program and --only keep theirs;
the default is every program). It runs NCBI in a cleaned environment (see `check`), because each
variable and file that is cleaned changes NCBI's result: BLASTDB, NCBI*, BLAST*, BATCH_SIZE, CHUNK_SIZE,
ARCHIVE_FORMAT, PRE_FETCH_SEQS_LIMIT, OLD_FSC, BL2SEQ_LEGACY, ADAPTIVE_CBS, CTOOLKIT_COMPATIBLE, and `~/.ncbirc`
(HOME is an empty directory); every NCBI run then gets BLAST_USAGE_REPORT=false (`ncbi_env`: NCBI BLAST+
sends no usage report with the run's parameters to NCBI). LOSAT's runs (`check`) keep the cleaned
environment as it is.

`check --losat BIN [--threads N] [--program NAME] [--only PREFIX] [--outfmts N,...] [--out TSV]` runs LOSAT with `-num_threads N`
(default 1) in the same cleaned environment (no NCBI*, BLAST*, LOSAT_*, BATCH_SIZE, CHUNK_SIZE, PRE_FETCH_SEQS_LIMIT, OLD_FSC,
BL2SEQ_LEGACY, ADAPTIVE_CBS or CTOOLKIT_COMPATIBLE variable; HOME and USERPROFILE are an empty directory, so that no `~/.ncbirc` names a
taxonomy database or a BLASTDB), compares stdout,
stderr and exit status with the frozen files (and checks that each frozen file's sha256 is the
manifest's), prints one line per row (`same`, `rejection`, or `differs` with the first differing
line), writes the same lines to --out, and exits 0 only when every row is `same` or `rejection`.
`--outfmts 5,13` keeps only the rows whose output format (the number that starts the `outfmt` cell) is one of
those numbers (`check_wasm_threading.py` runs blastx's XML, XML2 and JSON rows through the WASI commands so).
With `--library` BIN is `examples/library_command_line`, the same search through the library API: a
failure's exit status counts as any other failure's, a rejection row writes nothing to stdout there
too (the library holds its output where the command line holds it), a row whose query the reader refuses is
`library-input-error` when the message is the command line's (the library reads the whole query
before it searches, so the reports of the batches before the error do not exist; so is a rejection
row whose query has a reader error after the batch that the command line refuses). psiblast's rounds
(`-num_iterations` other than 1), its PSSM files (the PSSM kinds) and the PSSM queries of psiblast
and tblastn (`-in_pssm`, the library's `Inputs::pssm`) run through the library as on the command line;
the rows of rpsblast and rpstblastn (the library API has no such program), those that give both
`-query` and `-in_pssm` (the library's query is one or the other) and those that read the query from
standard input (the library's query is records) are `library-skip`.
With `--libm-last-digit NATIVE` (PD-LOSAT-PLATFORM-LIBM; `check_wasm_threading.py` passes it for the WASI
commands, with its native command) a row whose PSSM files differ from the frozen ones only in the last of
the 15 digits of ASN.1 REAL values, by 1 (`xinclude_fixtures.libm_last_digit_differences`), is
`libm-last-digit` (it counts as the same; the files and their number of such values are listed), and the
summary line adds `libm_last_digit=N`. A PSSM file over 64 KB is not frozen: its expected bytes are those of
NATIVE's run of the row, which must have the frozen SHA-256. Standard output, standard error, the exit
status and every other byte stay exact; without the option every byte counts.
The tool takes no lock; whoever runs it next to other NCBI runs holds the caller's own.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import shlex
import subprocess
import sys
import shutil
import tempfile
from collections.abc import Callable
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[1]
TESTS = ENGINE / "tests"
sys.path.insert(0, str(TESTS))
import xinclude_fixtures as xinclude  # noqa: E402  (the rows of outfmt 13 and 14)

MANIFEST = TESTS / "db_search_manifest"  # the directory: <program>.tsv for each of PROGRAMS, <program>.<family>.tsv
FROZEN = TESTS / "fixtures" / "db_search"
COLUMNS = ["row_id", "program", "kind", "db", "query", "extra_args", "outfmt", "rc", "stdout_sha256", "stdout_bytes",
           "stderr_sha256", "covers"]
PROGRAMS = ("blastp", "blastn", "tblastn", "blastx", "tblastx", "psiblast", "rpsblast", "rpstblastn")
NCBI_ENV = ("BATCH_SIZE", "CHUNK_SIZE", "OVERLAP_CHUNK_SIZE", "OLD_FSC", "BL2SEQ_LEGACY", "ARCHIVE_FORMAT",
            "PRE_FETCH_SEQS_LIMIT", "ADAPTIVE_CBS", "CTOOLKIT_COMPATIBLE")
TIMEOUT = 600
# The `NAME=VALUE` tokens of `extra_args` that are environment variables of both programs: BLASTDB
# (the directory of a `taxdb.bti`), the three that cut a query file into batches and chunks, and
# ARCHIVE_FORMAT (how outfmt 11 writes its archive: unset text, `bin...` binary), PRE_FETCH_SEQS_LIMIT
# (the number of Seq-aligns above which NCBI's formatter reads the reported sequences one by one), and
# BL2SEQ_LEGACY (read for -subject only; a -db search ignores it), OLD_FSC (no Gumbel block, any value), and
# ADAPTIVE_CBS (the larger preliminary hit list and its trim with composition-based statistics, any value).
ENV_NAMES = ("BLASTDB", "BATCH_SIZE", "CHUNK_SIZE", "OVERLAP_CHUNK_SIZE", "ARCHIVE_FORMAT", "PRE_FETCH_SEQS_LIMIT",
             "BL2SEQ_LEGACY", "OLD_FSC", "ADAPTIVE_CBS")


def manifest_file(program: str) -> Path:
    """The file of a program's rows, but those of its row families that have their own file."""
    return MANIFEST / f"{program}.tsv"


def row_family(row_id: str) -> str:
    """The family of a row: the second word of its id (`db.ldo.bp.before.0` is a row of `ldo`)."""
    words = row_id.split(".")
    return words[1] if len(words) > 1 else ""


def family_files(program: str) -> dict[str, Path]:
    """The row families of `program` that have their own file `<program>.<family>.tsv` (a file over
    512 KiB moves a family there, `check_manifest_sizes.py`), by family, in the order of their names."""
    return {path.name[len(program) + 1:-len(".tsv")]: path for path in sorted(MANIFEST.glob(f"{program}.*.tsv"))}


def program_files(program: str) -> list[Path]:
    """The files of `program`'s rows in the order `read_manifest` reads them: its file, then its family files."""
    return [manifest_file(program), *family_files(program).values()]


def row_file(row: dict[str, str], families: dict[str, Path]) -> Path:
    """The one file of `row`: its family's own file where its program has one (`families`, of
    `family_files`), else its program's file. `read_manifest` checks it and `write_manifest` writes there."""
    return families.get(row_family(row["row_id"]), manifest_file(row["program"]))


def read_manifest() -> tuple[dict[str, list[str]], list[dict[str, str]]]:
    """The comment lines of each file (by file name) and all rows: the programs in the order of
    `PROGRAMS`, each program's file first, then its family files (`family_files`), each file's rows
    in its order. A row in the file of another program, or in another file of its program than its
    own (`row_file`), is an error."""
    comments: dict[str, list[str]] = {}
    rows: list[dict[str, str]] = []
    for program in PROGRAMS:
        if not manifest_file(program).exists():
            # A missing file would drop its rows from every check without notice.
            raise SystemExit(f"{MANIFEST.name}/{manifest_file(program).name}: missing (every program in PROGRAMS"
                             " has a file)")
        families = family_files(program)
        for path in program_files(program):
            comments[path.name] = []
            header: list[str] = []
            for line in path.read_text().split("\n"):
                if line.startswith("#"):
                    comments[path.name].append(line)
                elif line and not header:
                    header = line.split("\t")
                    if header not in (COLUMNS, COLUMNS + xinclude.FILES_COLUMNS):
                        raise SystemExit(f"{path.name}: unexpected header {line!r}")
                elif line:
                    cells = line.split("\t")
                    if not len(COLUMNS) <= len(cells) <= len(header):
                        raise SystemExit(f"{path.name}: {len(cells)} cells in {line[:60]!r}")
                    row = {column: "" for column in header}
                    row.update(zip(header, cells))
                    if row["program"] != program:
                        raise SystemExit(f"{path.name}: the row {row['row_id']} is a {row['program']!r} row")
                    if row_file(row, families) != path:
                        raise SystemExit(f"{path.name}: the row {row['row_id']} belongs in "
                                         f"{row_file(row, families).name}")
                    rows.append(row)
    return comments, rows


def row_cells(row: dict[str, str]) -> list[str]:
    """The cells of a row in the manifest: the `files` columns are left out when they are empty."""
    cells = [row[column] for column in COLUMNS]
    if row.get("files") or row.get("side_sha256"):
        cells += [row.get(column, "") for column in xinclude.FILES_COLUMNS]
    return cells


def header_cells(rows: list[dict[str, str]]) -> list[str]:
    """The header: the `files` columns when some row has them."""
    files = any(row.get("files") or row.get("side_sha256") for row in rows)
    return COLUMNS + xinclude.FILES_COLUMNS if files else COLUMNS


def write_manifest(comments: dict[str, list[str]], rows: list[dict[str, str]]) -> None:
    """Write each file: its comment lines (`comments`, by file name), the header (the `files` columns
    when some row of any program has them, so the files keep one header) and its rows (`row_file`).
    A file whose text is the same is left alone."""
    header = "\t".join(header_cells(rows))
    MANIFEST.mkdir(exist_ok=True)
    for program in PROGRAMS:
        families = family_files(program)
        for path in program_files(program):
            own = [row for row in rows if row["program"] == program and row_file(row, families) == path]
            if not own and not path.exists():
                continue
            lines = comments.get(path.name, []) + [header] + ["\t".join(row_cells(row)) for row in own]
            text = "\n".join(lines) + "\n"
            if not path.exists() or path.read_text() != text:
                path.write_text(text)


def extra_tokens(row: dict[str, str]) -> list[str]:
    return shlex.split(row["extra_args"])


def is_env_token(token: str) -> bool:
    return token.partition("=")[0] in ENV_NAMES


def rps_cell(cell: str) -> tuple[list[str], set[str], dict[str, str], str | None] | None:
    """The profile databases of an rpsblast `db` cell `rps-gz:NAME[+NAME2...][!EXT...][!EXT=FILE...][@ALIAS]`:
    the names (`-db "NAME NAME2"`), the extensions left out (`!aux`, `!freq`), the files put in their
    place (`!aux=FILE`: LOSAT/tests/db_search/rps/FILE as each NAME.aux, inflated when FILE ends in .gz) and
    the alias file the run searches (`@ALIAS`: LOSAT/tests/db_search/rps/ALIAS.pal beside the inflated files,
    `-db ALIAS`; its DBLIST names them); None for another cell."""
    if not cell.startswith("rps-gz:"):
        return None
    body, at, alias = cell.removeprefix("rps-gz:").partition("@")
    names, *changes = body.split("!")
    dropped = {"." + change for change in changes if "=" not in change}
    replaced = dict(("." + change.split("=")[0], change.split("=")[1]) for change in changes if "=" in change)
    return names.split("+"), dropped, replaced, (alias if at else None)


def inflate_rps(cell: str) -> Path:
    """A fresh directory (its path with links resolved, as NCBI's RPS errors print it) holding the
    files LOSAT/tests/db_search/rps/NAME.*.gz of the cell's databases, inflated, but the extensions it
    leaves out or replaces, and the cell's alias file."""
    names, dropped, replaced, alias = rps_cell(cell)
    rps_dir = TESTS / "db_search" / "rps"
    target = Path(tempfile.mkdtemp(prefix="rps_db_")).resolve()
    for name in names:
        sources = sorted(rps_dir.glob(f"{name}.*.gz"))
        if not sources:
            raise SystemExit(f"no files LOSAT/tests/db_search/rps/{name}.*.gz")
        for source in sources:
            inflated = source.name.removesuffix(".gz")
            if inflated[len(name):] in dropped or inflated[len(name):] in replaced:
                continue
            (target / inflated).write_bytes(gzip.decompress(source.read_bytes()))
        for extension, file in replaced.items():
            data = (rps_dir / file).read_bytes()
            (target / f"{name}{extension}").write_bytes(gzip.decompress(data) if file.endswith(".gz") else data)
    if alias is not None:
        (target / f"{alias}.pal").write_bytes((rps_dir / f"{alias}.pal").read_bytes())
    return target


def rps_dir_token(row: dict[str, str]) -> bool:
    """Whether the row's outputs name the database's directory, frozen as the token `{db_dir}`."""
    return "path: {db_dir}" in row["covers"]


def subjects_argv(cell: str) -> list[str]:
    """The arguments of a `db` cell: `-db DB`, `-subject PATH` for a cell `subject:PATH`, or `-db NAME`
    (`-db ALIAS` with `@ALIAS`) for an rpsblast cell `rps-gz:NAME...` (`rps_cell`; the run's BLASTDB is the
    inflated directory)."""
    if cell.startswith("subject:"):
        return ["-subject", cell.removeprefix("subject:")]
    if (rps := rps_cell(cell)) is not None:
        return ["-db", rps[3] if rps[3] is not None else " ".join(rps[0])]
    return ["-db", cell]


def query_stdin(row: dict[str, str]) -> tuple[str, Path | None] | None:
    """The standard input of a row whose query cell is `stdin:pipe:PATH`, `stdin:file:PATH` or `stdin:empty`
    (`-query -`): ("pipe", the file PATH relative to LOSAT/, whose bytes the run reads through a pipe),
    ("file", PATH, the file itself as the run's standard input) or ("empty", None, an empty pipe). None for
    any other cell (the run keeps the caller's standard input and reads none)."""
    cell = row["query"]
    if not cell.startswith("stdin:"):
        return None
    form = cell.removeprefix("stdin:")
    if form == "empty":
        return "empty", None
    mode, _, path = form.partition(":")
    if mode not in ("pipe", "file") or not path:
        raise SystemExit(f"{row['row_id']}: query cell {cell!r} is not stdin:pipe:PATH, stdin:file:PATH or stdin:empty")
    return mode, ENGINE / path


def run_search(command: list[str], row: dict[str, str], **kw) -> subprocess.CompletedProcess:
    """`subprocess.run` of a search with its output captured and the row's standard input (`query_stdin`)."""
    form = query_stdin(row)
    if form is None:
        return subprocess.run(command, capture_output=True, **kw)
    mode, path = form
    if mode == "file":
        with open(path, "rb") as handle:
            return subprocess.run(command, capture_output=True, stdin=handle, **kw)
    return subprocess.run(command, capture_output=True, input=path.read_bytes() if path else b"", **kw)


def search_argv(row: dict[str, str]) -> list[str]:
    """The arguments after the program name (the same for NCBI and LOSAT)."""
    extra = [token for token in extra_tokens(row) if not is_env_token(token)]
    query = ["-query", "-"] if query_stdin(row) else xinclude.query_argv(row["query"], ENGINE)
    return [*subjects_argv(row["db"]), *query, *extra, "-outfmt", row["outfmt"]]


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
    set to the empty directory `home` (`$HOME/.ncbirc` is read by both programs). LOSAT_X_REACH is kept
    (LOSAT/src/utils/xreach.rs): it records LOSAT's switch paths and changes no output."""
    env = {k: v for k, v in base.items()
           if k == "LOSAT_X_REACH"
           or (not k.startswith(("NCBI", "BLAST", "LOSAT_")) and k not in NCBI_ENV)}
    env["HOME"] = env["USERPROFILE"] = home
    return env


# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_usage_report.cpp:195-220 (the
# report's parameters, sent to NCBI's server when the application ends unless BLAST_USAGE_REPORT, in the
# environment or the registry, is a false value)
NCBI_RUN_ENV = {"BLAST_USAGE_REPORT": "false"}


def ncbi_env(base: dict[str, str]) -> dict[str, str]:
    """The environment of every NCBI run of the harness (`freeze`): `base`, the cleaned environment, with
    BLAST_USAGE_REPORT=false, so that NCBI BLAST+ 2.17.0 does not send the run's parameters to NCBI at exit.
    A row cannot set it (it is not one of `ENV_NAMES`). LOSAT's runs keep `base`."""
    return {**base, **NCBI_RUN_ENV}


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


def frozen_stdout(row: dict[str, str], stdout: bytes, returncode: int = 0) -> bytes:
    """The stdout that `freeze` writes for `row`: NCBI's, except in two cases that keep the frozen bytes.

    1. A SAM document (outfmt 17) whose bytes differ from the frozen file's only in the CL tag of its @PG
       line. NCBI writes its argv[0] in that tag, the path of the NCBI build that ran (the oracle,
       Bioconda's, ws's). Without this rule a freeze with another build of the same NCBI release would
       rewrite the file and its sha256 and byte cells although the search output is the same by the rule
       that `check` (and set D's `run_oracle.py`) compares with, and a freeze of the candidate would leave
       `git status` dirty.
    2. A run that NCBI ends by a signal (`returncode` < 0) whose stdout and the frozen file are the same
       up to where the shorter one stops. Such rows are `losat-rejection` rows of a crash (a segment that
       NCBI frees and still walks): how much of its report NCBI writes before the crash depends on its
       heap and the environment (the same official 2.17.0 binary wrote 17 bytes of
       db.x2g.bx.of.p0127.nd0.1 on ws and 349 on this machine), so a freeze cannot reproduce the cut,
       and the frozen file holds the bytes that LOSAT writes before it refuses. A stdout that differs
       before that point is a new answer and replaces the file."""
    path = FROZEN / f"{row['row_id']}.out"
    if not path.is_file():
        return stdout
    if returncode < 0:
        frozen = path.read_bytes()
        return frozen if frozen.startswith(stdout) or stdout.startswith(frozen) else stdout
    if not is_sam_outfmt(row["outfmt"]):
        return stdout
    frozen = path.read_bytes()
    return frozen if strip_sam_cl(frozen) == strip_sam_cl(stdout) else stdout


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


def near_row_id(row: dict[str, str]) -> str | None:
    """The `ncbi` row whose frozen NCBI run of a near input is the expectation of `row` (`ncbi-near: ROW_ID` in its
    covers cell), or None: the row's own frozen files."""
    _, marker, rest = row["covers"].partition("ncbi-near: ")
    return rest.split(" | ")[0].split()[0] if marker else None


def under_valgrind(row: dict[str, str]) -> bool:
    """Whether the frozen bytes of `row` are NCBI's run under valgrind (a ` | ncbi-valgrind` part of its covers
    cell): rows where NCBI's native run reads past its buffers (PD-LOSAT-NCBI-DEFECTS exception 3)."""
    return any(part.strip() == "ncbi-valgrind" for part in row["covers"].split(" | ")[1:])


def valgrind_messages(stderr: bytes) -> bool:
    """Whether `stderr` holds valgrind's own lines (`==PID== ...`): `freeze` refuses to freeze them as NCBI's
    stderr (the `--valgrind` command needs a `--log-file`)."""
    return re.search(rb"^==\d+==", stderr, re.MULTILINE) is not None


def refuse_valgrind_messages(row: dict[str, str], stderr: bytes) -> None:
    """Stops `freeze` before it writes any file of an `ncbi-valgrind` row whose NCBI stderr holds valgrind's lines."""
    if under_valgrind(row) and valgrind_messages(stderr):
        raise SystemExit(f"{row['row_id']}: valgrind's messages in NCBI's stderr; give --valgrind a --log-file")


def frozen_id(row: dict[str, str]) -> str:
    """The row id that names the frozen files of `row`'s expectation (`near_row_id`, else its own)."""
    return near_row_id(row) or row["row_id"]


def ncbi_stack(row: dict[str, str]) -> Callable[[], None] | None:
    """The `preexec_fn` of `freeze`'s NCBI run of `row`: with `ncbi-stack: unlimited` in its covers cell, the run's
    stack limit is its hard limit (`RLIMIT_STACK` infinity by default, as `prlimit --stack=unlimited`); else None,
    the caller's limit."""
    if "ncbi-stack: unlimited" not in row["covers"]:
        return None
    import resource  # POSIX: NCBI's runs of `freeze` are on Linux

    def unlimited() -> None:
        _, hard = resource.getrlimit(resource.RLIMIT_STACK)
        resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))

    return unlimited


def format_number(row: dict[str, str]) -> str:
    """The number of the row's output format: the first word of its `outfmt` cell (a custom spec follows it)."""
    return (row["outfmt"].split() or [""])[0]


def selected(row: dict[str, str], args: argparse.Namespace) -> bool:
    formats = [number for number in getattr(args, "outfmts", "").split(",") if number]
    return (row["row_id"].startswith(args.only) and args.program in ("", row["program"])
            and (not formats or format_number(row) in formats))


def freeze(args: argparse.Namespace) -> int:
    comments, rows = read_manifest()
    FROZEN.mkdir(parents=True, exist_ok=True)
    home = tempfile.TemporaryDirectory(prefix="db_search_home_")
    # Every NCBI run below (and `freeze_rps`'s) takes this environment, with the usage report off.
    env = ncbi_env(clean_env(dict(os.environ), home.name))
    for row in rows:
        if not selected(row, args):
            continue
        command = [str(Path(args.ncbi_bin) / row["program"]), *search_argv(row), "-num_threads", "1"]
        if row["kind"] == "losat-rejection" and row["stdout_bytes"]:
            # NCBI's stdout of the command that LOSAT refuses after the reports NCBI writes too.
            if xinclude.wants_pssm(row):
                run = xinclude.run_with_pssm(command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE,
                                             timeout=TIMEOUT, env=row_env(row, env), preexec_fn=ncbi_stack(row)).run
            else:
                run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env=row_env(row, env),
                                 preexec_fn=ncbi_stack(row))
            out = frozen_stdout(row, run.stdout, run.returncode)
            (FROZEN / f"{row['row_id']}.out").write_bytes(out)
            row.update(stdout_sha256=sha256(out), stdout_bytes=str(len(out)))
            print(f"{row['row_id']}\tncbi_rc={run.returncode}\tstdout={len(out)} (rejection row: stdout only)")
            continue
        if row["kind"] != "ncbi" or near_row_id(row):
            continue
        if under_valgrind(row):
            if not args.valgrind:
                print(f"{row['row_id']}\tkept: frozen from NCBI under valgrind (freeze it with --valgrind)")
                continue
            command = [*shlex.split(args.valgrind), *command]
        files = xinclude.wants_files(row)
        pssm = xinclude.wants_pssm(row)
        if (files or pssm) and query_stdin(row):
            raise SystemExit(f"{row['row_id']}: a row with files reads no query from standard input")
        if rps_cell(row["db"]) is not None:
            freeze_rps(row, command, env)
            continue
        if files:
            written = xinclude.run_with_files(command, row["row_id"], row["outfmt"], cwd=ENGINE,
                                              timeout=TIMEOUT, env=row_env(row, env), preexec_fn=ncbi_stack(row))
            run = written.run
            refuse_valgrind_messages(row, run.stderr)
            out = run.stdout
            xinclude.freeze(written, FROZEN, row["row_id"], row["outfmt"])
        elif pssm:
            pssms = xinclude.run_with_pssm(command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE,
                                           timeout=TIMEOUT, env=row_env(row, env), preexec_fn=ncbi_stack(row))
            run = pssms.run
            refuse_valgrind_messages(row, run.stderr)
            xinclude.freeze_pssm(pssms, FROZEN, row["row_id"], xinclude.files_kind(row))
            out = frozen_stdout(row, run.stdout, run.returncode)
            (FROZEN / f"{row['row_id']}.out").write_bytes(out)
        else:
            run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env=row_env(row, env),
                             preexec_fn=ncbi_stack(row))
            refuse_valgrind_messages(row, run.stderr)
            out = frozen_stdout(row, run.stdout, run.returncode)
            (FROZEN / f"{row['row_id']}.out").write_bytes(out)
        err = FROZEN / f"{row['row_id']}.err"
        if run.stderr:
            err.write_bytes(run.stderr)
        else:
            err.unlink(missing_ok=True)
        row.update(rc=str(run.returncode), stdout_sha256=sha256(out), stdout_bytes=str(len(out)),
                   stderr_sha256=sha256(run.stderr) if run.stderr else "")
        if files:
            row.update(xinclude.observed(written))
        if pssm:
            row.update(xinclude.observed_pssm(pssms, row["row_id"], xinclude.files_kind(row)))
        print(f"{row['row_id']}\trc={run.returncode}\tstdout={len(out)}\tstderr={len(run.stderr)}"
              + (f"\tsides={len(written.sides)}" if files else "")
              + (f"\tpssm_files={len(pssms.files)}" if pssm else ""))
    # The rows whose expectation is a near input's NCBI run (`near_row_id`) take that row's columns, after it
    # is frozen (it may come later in the file).
    by_id = {row["row_id"]: row for row in rows}
    for row in rows:
        near = near_row_id(row)
        if not selected(row, args) or row["kind"] != "ncbi" or not near:
            continue
        source = by_id.get(near)
        if source is None or source["kind"] != "ncbi" or near_row_id(source) or source["outfmt"] != row["outfmt"]:
            raise SystemExit(f"{row['row_id']}: ncbi-near {near} is not an ncbi row of its own run with the same outfmt")
        row.update({key: source[key] for key in ("rc", "stdout_sha256", "stdout_bytes", "stderr_sha256")})
        print(f"{row['row_id']}\tncbi-near={near}\trc={row['rc']}\tstdout={row['stdout_bytes']}")
    write_manifest(comments, rows)
    return 0


def freeze_rps(row: dict[str, str], command: list[str], env: dict[str, str]) -> None:
    """`freeze` of an rpsblast row: NCBI's run with BLASTDB the inflated database, its directory
    written as `{db_dir}` where the row says so."""
    db_dir = inflate_rps(row["db"])
    try:
        run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT,
                         env={**row_env(row, env), "BLASTDB": str(db_dir)})
    finally:
        shutil.rmtree(db_dir, ignore_errors=True)
    stdout, stderr = run.stdout, run.stderr
    if rps_dir_token(row):
        stdout = stdout.replace(str(db_dir).encode(), b"{db_dir}")
        stderr = stderr.replace(str(db_dir).encode(), b"{db_dir}")
    (FROZEN / f"{row['row_id']}.out").write_bytes(stdout)
    err = FROZEN / f"{row['row_id']}.err"
    if stderr:
        err.write_bytes(stderr)
    else:
        err.unlink(missing_ok=True)
    row.update(rc=str(run.returncode), stdout_sha256=sha256(stdout), stdout_bytes=str(len(stdout)),
               stderr_sha256=sha256(stderr) if stderr else "")
    print(f"{row['row_id']}\trc={run.returncode}\tstdout={len(stdout)}\tstderr={len(stderr)}")


def status_differs(observed: int, expected: int, library: bool) -> bool:
    """Whether a run's exit status is not the frozen one. `examples/library_command_line.rs` ends every
    failure with status 1: the library API's errors carry the command line's message, not its status."""
    return (observed != 0) != (expected != 0) if library else observed != expected


def check_row(row: dict[str, str], losat: Path, threads: int, env: dict[str, str], library: bool = False,
              libm_native: Path | None = None) -> tuple[str, str]:
    """The row's status (`same`, `rejection`, `differs`, `libm-last-digit`, `library-skip`,
    `library-input-error`) and its detail. With `libm_native` (`--libm-last-digit`), the PSSM files
    follow PD-LOSAT-PLATFORM-LIBM's rule (`xinclude.pssm_problems`): a file over 64 KB, which is not
    frozen, is compared with the file of `libm_native`'s run of the row, which must have the frozen
    SHA-256. Standard output, standard error and the exit status are compared byte for byte."""
    files = xinclude.wants_files(row)
    pssm = xinclude.wants_pssm(row)
    # The library API refuses rpsblast and rpstblastn (`api::Program::parse`, line PS items PS11a, PS12a).
    # Its query is records or a PSSM file (`Inputs::new`, `Inputs::pssm`), never both: a row that gives
    # the command line `-query` and `-in_pssm` (the parser's exclusion) has no library call.
    argv = search_argv(row)
    if library and (row["program"] in ("rpsblast", "rpstblastn")
                    or ("-query" in argv and "-in_pssm" in argv) or query_stdin(row)):
        return "library-skip", ""
    if (files or pssm) and query_stdin(row):
        raise SystemExit(f"{row['row_id']}: a row with files reads no query from standard input")
    command = [str(losat), row["program"], *search_argv(row), "-num_threads", str(threads)]
    db_dir = inflate_rps(row["db"]) if rps_cell(row["db"]) is not None else None
    if db_dir is not None:
        try:
            run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT,
                             env={**row_env(row, env), "BLASTDB": str(db_dir)})
        finally:
            shutil.rmtree(db_dir, ignore_errors=True)
        if rps_dir_token(row):
            run.stdout = run.stdout.replace(str(db_dir).encode(), b"{db_dir}")
            run.stderr = run.stderr.replace(str(db_dir).encode(), b"{db_dir}")
    elif files:
        written = xinclude.run_with_files(command, row["row_id"], row["outfmt"], cwd=ENGINE,
                                          env=row_env(row, env), timeout=TIMEOUT)
        run = written.run
    elif pssm:
        pssms = xinclude.run_with_pssm(command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE,
                                       env=row_env(row, env), timeout=TIMEOUT)
        run = pssms.run
    else:
        run = run_search(command, row, cwd=ENGINE, env=row_env(row, env), timeout=TIMEOUT)
    if row["kind"] == "losat-rejection":
        # The library reads the whole query before it searches: a reader error ends the call before
        # the batch that the command line refuses, with the reader's message.
        if library and run.returncode != 0 and run.stdout == b"" and run.stderr.startswith(b"BLAST query error:") \
                and expected_text(row).encode() not in run.stderr:
            return "library-input-error", ""
        problems = []
        if status_differs(run.returncode, int(row["rc"]), library):
            problems.append(f"exit {run.returncode}, expected {row['rc']}")
        # A refusal follows no byte in either path (the library holds its output where the command
        # line holds it, BLASTX's XML documents through `blast_xml::documents_may_refuse`); the rows
        # whose refusal comes after NCBI's reports are psiblast's PSSM file and rounds rows, which the
        # library skips.
        if row["stdout_bytes"]:
            out_path = FROZEN / f"{row['row_id']}.out"
            out = out_path.read_bytes() if out_path.exists() else b""
            if sha256(out) != row["stdout_sha256"] or len(out) != int(row["stdout_bytes"]):
                problems.append("frozen .out does not match the manifest")
            if run.stdout != out:
                problems.append("stdout " + first_difference(out, run.stdout))
        elif run.stdout:
            problems.append(f"stdout is not empty ({len(run.stdout)} bytes)")
        if expected_text(row).encode() not in run.stderr:
            problems.append(f"stderr lacks {expected_text(row)!r}: {run.stderr.decode(errors='replace').strip()[:150]!r}")
        if pssm and (pssms.files or pssms.stray):
            problems.append(f"a refused run wrote {[name for name, _ in pssms.files] + pssms.stray}")
        return ("differs", "; ".join(problems)) if problems else ("rejection", "")
    out = (FROZEN / f"{frozen_id(row)}.out").read_bytes()
    err_path = FROZEN / f"{frozen_id(row)}.err"
    err = err_path.read_bytes() if row["stderr_sha256"] else b""
    expected_out_losat = out.replace(
        b"Please refer to the BLAST+ user manual.",
        b"Please refer to the LOSAT user manual.",
    )
    expected_err_losat = err.replace(
        b"Please refer to the BLAST+ user manual.",
        b"Please refer to the LOSAT user manual.",
    )
    # The library reads the whole query before it searches: a reader error ends the call there,
    # with the command line's message of it, before the reports of the batches before it.
    if library and run.returncode != 0 and int(row["rc"]) != 0 and run.stdout == b"" \
            and run.stderr.startswith(b"BLAST query error:") and run.stderr in expected_err_losat:
        return "library-input-error", ""
    problems = []
    if files:
        problems += xinclude.problems(written, row, FROZEN, first_difference)
    elif sha256(out) != row["stdout_sha256"] or len(out) != int(row["stdout_bytes"]):
        problems.append("frozen .out does not match the manifest")
    if (sha256(err) if err else "") != row["stderr_sha256"] or err_path.exists() != bool(row["stderr_sha256"]):
        problems.append("frozen .err does not match the manifest")
    if status_differs(run.returncode, int(row["rc"]), library):
        problems.append(f"exit {run.returncode}, expected {row['rc']}: {run.stderr.decode(errors='replace').strip()[:150]!r}")
    # A SAM row differs in the command line of the @PG line only (`strip_sam_cl`).
    if not files:
        expected_out, observed_out = (expected_out_losat, run.stdout)
        if is_sam_outfmt(row["outfmt"]):
            expected_out, observed_out = strip_sam_cl(expected_out_losat), strip_sam_cl(run.stdout)
        if observed_out != expected_out:
            problems.append("stdout " + first_difference(expected_out, observed_out))
    if run.stderr != expected_err_losat:
        problems.append("stderr " + first_difference(expected_err_losat, run.stderr))
    accepted: list[str] = []
    if pssm:
        reference = None
        if libm_native is not None:
            native_files: list[dict[str, bytes]] = []

            def reference(name: str) -> bytes | None:
                # `libm_native`'s run of the row, once, for the files that are not frozen.
                if not native_files:
                    command = [str(libm_native), row["program"], *search_argv(row), "-num_threads", str(threads)]
                    native_files.append(dict(xinclude.run_with_pssm(
                        command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE, env=row_env(row, env),
                        timeout=TIMEOUT).files))
                return native_files[0].get(name)
        problems += xinclude.pssm_problems(pssms, row, FROZEN, first_difference, reference, accepted)
    if problems:
        return "differs", "; ".join(problems)
    return ("libm-last-digit", "; ".join(accepted)) if accepted else ("same", "")


def check(args: argparse.Namespace) -> int:
    home = tempfile.TemporaryDirectory(prefix="db_search_home_")
    env = clean_env(dict(os.environ), home.name)
    losat = args.losat.resolve()
    libm_native = args.libm_last_digit.resolve() if args.libm_last_digit else None
    lines, bad, libm = [], 0, 0
    for row in read_manifest()[1]:
        if not selected(row, args):
            continue
        status, detail = check_row(row, losat, args.threads, env, args.library, libm_native)
        bad += status == "differs"
        libm += status == "libm-last-digit"
        lines.append(f"{row['row_id']}\t{status}\t{detail}")
    text = "".join(line + "\n" for line in lines)
    sys.stdout.write(text)
    if args.out:
        args.out.write_text(text)
    print(f"# rows={len(lines)} differing={bad} threads={args.threads}"
          + (f" libm_last_digit={libm}" if libm_native else ""))
    return 1 if bad else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    fz = sub.add_parser("freeze", help="run NCBI and write the frozen bytes and the manifest columns")
    fz.add_argument("--ncbi-bin", required=True)
    fz.add_argument("--program", choices=PROGRAMS, default="", help="only the rows of this program (default: all)")
    fz.add_argument("--only", default="")
    fz.add_argument("--valgrind", default="",
                    help="the command before NCBI's for the rows frozen under valgrind (`ncbi-valgrind` in covers)")
    ck = sub.add_parser("check", help="compare LOSAT with the frozen bytes")
    ck.add_argument("--losat", type=Path, required=True)
    ck.add_argument("--threads", type=int, default=1)
    ck.add_argument("--program", choices=PROGRAMS, default="", help="only the rows of this program (default: all)")
    ck.add_argument("--only", default="")
    ck.add_argument("--outfmts", default="",
                    help="only the rows of these output formats, numbers separated by commas (default: all)")
    ck.add_argument("--out", type=Path)
    ck.add_argument("--library", action="store_true",
                    help="--losat is examples/library_command_line (the search through the library API)")
    ck.add_argument("--libm-last-digit", type=Path, metavar="NATIVE",
                    help="PD-LOSAT-PLATFORM-LIBM (the WASI commands): a PSSM file may differ in the last of the 15 "
                         "digits of ASN.1 REAL values by 1; NATIVE is a LOSAT command whose run gives the frozen "
                         "bytes of a PSSM file that is not frozen")
    args = parser.parse_args()
    return freeze(args) if args.command == "freeze" else check(args)


if __name__ == "__main__":
    sys.exit(main())

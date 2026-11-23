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
`-query QUERY`, the gzipped file PATH inflated into a temporary directory, or PATH as it is; a psiblast query
cell `msa:PATH` is an alignment query, `-in_msa PATH`;
xinclude_fixtures.py `query_argv`; a query cell `stdin:pipe:PATH`, `stdin:file:PATH` or `stdin:empty` (for psiblast also `stdin-msa:...`, `-in_msa -`) is
`-query -` with the file PATH's bytes through a pipe, the file PATH itself, or an empty pipe as the run's
standard input, `query_stdin`)
(`blastp`, `blastn`, `tblastn`, `blastx`, `tblastx` or `psiblast`; NCBI's own binary of that name, or `LOSAT PROGRAM`) (an empty `db` cell is the literal empty
argument `-db ""`, and a cell `subject:PATH` runs `-subject PATH` in place of `-db DB`: psiblast's refusal of a
subject set whose sequences all have no letters, which only a subject set reaches; `extra_args` is split like a shell line, its tokens `BLASTDB=DIR` are not arguments but the environment variable of both programs (the directory of a `taxdb.bti`, relative to LOSAT/; an RPS row's search path is its inflated database's directory, then DIR, `rps_env`); `outfmt` is one argument, so a custom
`7 qseqid sseqid ...` spec stays whole). The databases are under LOSAT/tests/seqdb/db and
LOSAT/tests/db_search (README.md there says how they were built).

- kind `ncbi`: the expectation is NCBI's stdout, stderr and exit status (`rc`) of a run with
  `-num_threads 1`, every byte, nothing normalised. The frozen bytes are
  LOSAT/tests/fixtures/db_search/<row_id>.out and <row_id>.err (only when stderr is not empty).
- a `ncbi` row whose `covers` cell has `ncbi-near: ROW_ID` (up to a ` | ` or the end of the cell): the expectation
  is the frozen bytes, exit status and hash columns of the `ncbi` row ROW_ID, NCBI's run of a near input (the same
  search with the mask one letter shorter, in TBLASTN and TBLASTX the fewest bases after which every reading frame
  keeps a letter after its part; for exception 21 the freed letters stored as `X`, or `N`, so that they stay masked
  in the alignments and the statistics), where NCBI's own run of the row reaches an approved NCBI defect
  (PD-LOSAT-NCBI-DEFECTS exceptions 13 and 21: `-db_hard_mask` with a mask that runs to the end, on which NCBI BLAST+
  2.17.0 crashes, adds rows from the masked letters, drops seeds of the next subject or writes an HSP many times;
  exception 3: a one-hit gapped start whose window reads past the end of a database volume's last memory page, where
  NCBI faults; exception 19: a hit list size whose preliminary size passes NCBI's Int4, with the near size
  1073741798; exception 23: BLASTX `-use_sw_tback`'s heap test in e-value mode past NCBI's
  array of a batch's heaps, with the near input that follows each query batch by queries without a hit, where
  they add no byte: outfmt 6, or a batch that fails and writes nothing; exception 26: `-db_hard_mask` with
  overlapping masks after which NCBI's word hits reach diagonal entries of another subject, with the near database
  that holds a spacer of `X` between the two subjects under an alias with the database's statistics). The row has no
  frozen files of its own; `freeze` does not run NCBI for it and copies ROW_ID's columns (`near_row`). ROW_ID is
  an `ncbi` row of the same outfmt without a near row of its own, and an `ncbi` row takes all of its bytes.
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
- a `ncbi` row whose `covers` cell has `ncbi-copies: PATH [NAME=VALUE ...]` (up to a ` | ` or the end; a "copy row",
  `copies_cell`): the expectation is made from NCBI's run on the oracle query file PATH (relative to LOSAT/) in place
  of the row's query file. In PATH each query that ends a query batch of LOSAT's run (`query_batches`) is followed by
  a copy of itself (its residues under an id that starts with `COPY`: `>COPY`, `>COPY2`), and only the copies end the
  oracle run's batches (`copy_oracle_problems` checks both files). With unified P-values (`blastp -comp_based_stats
  1u`, `2u`, `3u`) NCBI BLAST+ 2.17.0 adjusts a query's E-values with the values of the next query of its batch and,
  for a batch's last query, with bytes past its array that change from run to run; in the oracle run each original
  query's next query is itself or its copy, so its report is NCBI's with its own values where LOSAT's run has no next
  query. The `NAME=VALUE` words (BATCH_SIZE, CHUNK_SIZE or OVERLAP_CHUNK_SIZE: `COPY_ENV_NAMES`) are the oracle
  run's environment in place of the row's own words of those names. Each batch of both runs is searched in one
  chunk (`query_chunks`: the copy would move the split of a longer batch, and the split can change the preliminary
  hits). `freeze` runs NCBI
  `-num_threads 1` on PATH, which must exit 0 with a stderr that names no copy (`COPY`: a warning of the arguments,
  such as NCBI's advice to examine 5 or more matches under a smaller hit list, is the row's own and is frozen as
  <row_id>.err), keeps its whole stdout unchanged as <row_id>.copies.out (the evidence) and writes as <row_id>.out
  (with the row's hash and size cells; rc 0) the cut of it by a fixed rule per output format
  (`cut_copies`): outfmt 6, the lines whose first field starts with `COPY`; outfmt 7, each copy's block (its `# BLASTP
  2.17.0+` line up to the next such line or the final line), and the number of the final line `# BLAST processed N
  queries` lowered by the copies (the one NCBI byte that the rule rewrites); outfmt 0, each copy's report (its `Query=
  COPY` line up to the next `Query= ` line or, for the last query, up to the final `  Database: ` block); outfmt 5, the
  copy's `<Iteration>` element, a copy that is the oracle file's last query only (the iteration numbers of later
  queries would change). The copies' own reports are a batch's last query's, which change from run to run: a later
  freeze keeps the frozen whole output when its cut is the new run's cut. A copy row is a BLASTP row of outfmt 0, 5, 6
  or 7 with the default fields; the loader refuses any other. A hit list that fills does not limit a copy row
  (`-max_target_seqs`, `-num_descriptions` and `-num_alignments` are allowed): NCBI's composition-based redo takes the
  matches in the order of their queries and ends the remaining ones early (`BlastCompo_EarlyTermination`) only while
  the list of every query of the batch is full, and the lists of the queries after the one in work are still empty.
  So in the oracle run no match of an original query is ended early (its copy, or another query, follows it); only
  the copy's matches can be, and the copy's report is cut out. LOSAT ends no match early with unified P-values
  (`BlastCompoHeap::without_early_end`), so the cut is its expectation also where a list fills. `freeze` prints the
  number of subjects of each query of the oracle output (`subjects_per_query`). `check` cuts the frozen whole output
  again, fails the row when that differs from <row_id>.out, and compares LOSAT's run of the row's own query file with
  <row_id>.out as for any `ncbi` row.
- kind `losat-rejection`: LOSAT must exit with `rc`, write nothing to stdout, and its stderr must
  contain the text that follows `expect: ` in the `covers` cell (up to the end of the cell or a
  ` | ` separator). The hash columns stay empty. A rejection row with a PSSM files kind must also
  write no PSSM file. A rejection that LOSAT makes after reports that NCBI writes too (psiblast's
  `-out_pssm` file that cannot be written: NCBI stops after that query's report with its own error)
  has a non-empty `stdout_bytes` cell: its stdout must be NCBI's of the same command, frozen as
  <row_id>.out by `freeze` (which runs NCBI for such a row and keeps only its stdout) and checked
  against `stdout_sha256`. Its `rc` and stderr are LOSAT's. A rejection row whose `covers` cell has
  `ncbi-near: ROW_ID` takes its stdout from the `ncbi` row ROW_ID of the same outfmt, NCBI's run of a near
  input that writes the bytes LOSAT writes before it refuses, and one with `ncbi-near: ROW_ID first N bytes`
  takes the first N bytes of that run: an input that fails at the same step where NCBI crashes on the row's
  own (line RG item X4: the failed batch of a Smith-Waterman redo), the query batches before the batch that
  LOSAT refuses and NCBI searches on (X4: the negative CHUNK_SIZE pair whose chunks leave gaps; X5: a copy
  row whose query file holds the batches that LOSAT writes before it refuses, where NCBI's own run of the
  row reads past its array there; line WQA D10: two queries of one ID under `-parse_deflines`, the refused
  one's letters in the near input a line that NCBI's reader refuses), the rounds before the one that LOSAT
  refuses (psiblast with
  `-num_iterations 1`), or a search of the same database and
  output format whose first N bytes are the prolog that LOSAT writes before it refuses the report (line RG
  item T3c: the same search with a description table, which NCBI ends with an exception at an HSP of length
  0, or without a hit). `freeze` runs nothing for such a row and writes the sha256 and size of those bytes
  into its stdout columns, `check` reads them from ROW_ID's frozen stdout (`frozen_id`, `expected_stdout`),
  and the row keeps its own rc and stderr (LOSAT's). Where NCBI ends a row by a signal (a crash), its stdout
  stops where the crash cuts it, which its heap and the environment move (the same official 2.17.0 binary
  wrote 349 bytes of db.x2g.bx.of.p0127.nd0.1 when it was frozen and 17 in later runs): no expectation.
  `freeze` stops on such a run before it writes a file (`refuse_signal`); the row takes its stdout from a near
  row, or has empty stdout cells where LOSAT writes nothing before it refuses.

`freeze --ncbi-bin DIR [--program NAME] [--only PREFIX] [--valgrind CMD]` runs NCBI for the `ncbi` rows, writes the frozen files and
rewrites the rc, hash and size columns of the manifest (rows outside --program and --only keep theirs;
the default is every program), and stops with an error where NCBI ends a run by a signal, before it writes
that row's files or the manifest (`refuse_signal`). It runs NCBI in a cleaned environment (see `check`), because each
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
`-query` and `-in_pssm` (the library's query is one or the other), those that give `-in_msa` with either of them,
those with `-msa_master_idx` or `-ignore_msa_master` and no `-in_msa` (the library takes the options with the
alignment, `Inputs::msa`), those that write the search strategy file (`-export_search_strategy`, which the
library refuses) and those that read the query from
standard input (the library's query is records) are `library-skip`. psiblast's alignment queries
(`-in_msa`, a query cell `msa:PATH`, the library's `Inputs::msa`) run through the library as on the command line.
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
import losat_wording  # noqa: E402  (LOSAT's words in NCBI's frozen bytes, AGENTS.md rule 8)

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
    the names (`-db "NAME NAME2"`; `DIR/NAME`: NAME's files inflated into the subdirectory DIR, a copy of a
    volume under its own base name), the extensions left out (`!aux`, `!freq`), the files put in their
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


DB_DIR_LEN_FILE = TESTS / "db_search" / "rps" / "db_dir_len.txt"


def db_dir_len() -> int:
    """The length of a `{db_dir}` directory's path in the rows that keep it (`inflate_rps`); the
    file is read by `wording::tests::archive_messages_are_wrapped_again_as_ncbi_wraps_them` too."""
    return int(DB_DIR_LEN_FILE.read_text())


def temp_dir_of_length(prefix: str, length: int) -> Path:
    """A fresh directory under the temporary directory (links resolved; `/tmp` or `/var/tmp` when that
    path is too long) whose path is exactly `length` characters: the name's prefix is lengthened. The
    wraps of an ASN.1 message that holds the path (NCBI's RPS errors in an archive) depend on its
    length, so a frozen `{db_dir}` stands for one."""
    random_part = 8  # `tempfile.mkdtemp` appends eight random characters.
    parents = [Path(tempfile.gettempdir()).resolve(), Path("/tmp").resolve(), Path("/var/tmp").resolve()]
    for parent in parents:
        pad = length - len(str(parent)) - len(os.sep) - len(prefix) - random_part
        if pad >= 0 and parent.is_dir():
            target = Path(tempfile.mkdtemp(prefix=prefix + "p" * pad, dir=parent))
            if len(str(target)) != length:
                raise SystemExit(f"{target} has {len(str(target))} characters, not {length}")
            return target
    raise SystemExit(f"no temporary directory of {[str(p) for p in parents]} is short enough for a directory path "
                     f"of {length} characters")


def inflate_rps(cell: str, token: bool = False) -> Path:
    """A fresh directory (its path with links resolved, as NCBI's RPS errors print it) holding the
    files LOSAT/tests/db_search/rps/NAME.*.gz of the cell's databases, inflated, but the extensions it
    leaves out or replaces, and the cell's alias file. With `token` (the row keeps the directory as
    `{db_dir}`, `rps_dir_token`) its path has `db_dir_len()` characters."""
    names, dropped, replaced, alias = rps_cell(cell)
    rps_dir = TESTS / "db_search" / "rps"
    target = temp_dir_of_length("rps_db_", db_dir_len()) if token else Path(tempfile.mkdtemp(prefix="rps_db_")).resolve()
    for name in names:
        sub, _, name = name.rpartition("/")
        into = target / sub
        into.mkdir(parents=True, exist_ok=True)
        sources = sorted(rps_dir.glob(f"{name}.*.gz"))
        if not sources:
            raise SystemExit(f"no files LOSAT/tests/db_search/rps/{name}.*.gz")
        for source in sources:
            inflated = source.name.removesuffix(".gz")
            if inflated[len(name):] in dropped or inflated[len(name):] in replaced:
                continue
            (into / inflated).write_bytes(gzip.decompress(source.read_bytes()))
        for extension, file in replaced.items():
            data = (rps_dir / file).read_bytes()
            (into / f"{name}{extension}").write_bytes(gzip.decompress(data) if file.endswith(".gz") else data)
    if alias is not None:
        (target / f"{alias}.pal").write_bytes((rps_dir / f"{alias}.pal").read_bytes())
    return target


def rps_dir_token(row: dict[str, str]) -> bool:
    """Whether the row's outputs name the database's directory, frozen as the token `{db_dir}`."""
    return "path: {db_dir}" in row["covers"]


def rps_env(row: dict[str, str], env: dict[str, str], db_dir: Path) -> dict[str, str]:
    """`row_env` of an RPS row whose database is inflated in `db_dir` (`inflate_rps`): BLASTDB is that
    directory, then the row's own `BLASTDB=DIR` (the directory of a `taxdb.bti`) after it in the search
    path."""
    run_env = row_env(row, env)
    own = run_env.get("BLASTDB") if any(token.startswith("BLASTDB=") for token in extra_tokens(row)) else None
    run_env["BLASTDB"] = str(db_dir) if own is None else f"{db_dir}{os.pathsep}{own}"
    return run_env


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
    if not cell.startswith(STDIN_PREFIXES):
        return None
    form = cell.partition(":")[2]
    if form == "empty":
        return "empty", None
    mode, _, path = form.partition(":")
    if mode not in ("pipe", "file") or not path:
        raise SystemExit(f"{row['row_id']}: query cell {cell!r} is not stdin:pipe:PATH, stdin:file:PATH or stdin:empty "
                         "(or the same with stdin-msa:)")
    return mode, ENGINE / path


# A query cell that reads standard input is `stdin:...` (`-query -`) or, for psiblast's alignment,
# `stdin-msa:...` (`-in_msa -`, line PS item PS6); the part after the first colon is the same.
STDIN_PREFIXES = ("stdin:", "stdin-msa:")


def stdin_query_argv(row: dict[str, str]) -> list[str]:
    """The query arguments of a row that reads standard input: `-query -`, or `-in_msa -`."""
    return ["-in_msa", "-"] if row["query"].startswith("stdin-msa:") else ["-query", "-"]


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
    query = stdin_query_argv(row) if query_stdin(row) else xinclude.query_argv(row["query"], ENGINE)
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


# NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:195-220 (the
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


# The @PG line's start before its CL tag: NCBI's version (frozen bytes) or LOSAT's (LOSAT's output and
# NCBI's bytes in LOSAT's words, `losat_wording`).
_SAM_PGS = tuple(b"\n@PG\tID:0\tVN:" + version for version in (losat_wording.NCBI_VERSION, losat_wording.LOSAT_VERSION))


def strip_sam_cl(data: bytes) -> bytes:
    """The SAM document without the CL tag of its @PG line: the bytes from the `\tCL:` that follows
    `\n@PG\tID:0\tVN:2.17.0+` (NCBI's) or `\n@PG\tID:0\tVN:LOSAT 0.2.0` (LOSAT's) up to, not including,
    the first following `\tPN:blastn\n`.

    The one rule that the SAM (outfmt 17) rows of every fixture compare with, on the frozen NCBI bytes (as
    frozen, or in LOSAT's words) and on LOSAT's (DECISIONS 7 of the outfmt line): LOSAT writes its own command
    line (its program path, the `blastn` sub-command word and the arguments a harness adds, such as
    `-num_threads N`) where NCBI writes its argv. Every other byte must be equal. A document without the tag
    (no command line) is returned as it is."""
    for pg in _SAM_PGS:
        head = data.find(pg + b"\tCL:")
        if head >= 0:
            break
    else:
        return data
    start = head + len(pg)
    end = data.find(b"\tPN:blastn\n", start + len(b"\tCL:"))
    return data if end < 0 else data[:start] + data[end:]


def frozen_stdout(row: dict[str, str], stdout: bytes) -> bytes:
    """The stdout that `freeze` writes for `row`: NCBI's, except a SAM document (outfmt 17) whose bytes differ
    from the frozen file's only in the CL tag of its @PG line, which keeps the frozen bytes. NCBI writes its
    argv[0] in that tag, the path of the NCBI build that ran (the oracle, Bioconda's, ws's). Without this rule
    a freeze with another build of the same NCBI release would rewrite the file and its sha256 and byte cells
    although the search output is the same by the rule that `check` (and set D's `run_oracle.py`) compares
    with, and a freeze of the candidate would leave `git status` dirty."""
    path = FROZEN / f"{row['row_id']}.out"
    if not path.is_file() or not is_sam_outfmt(row["outfmt"]):
        return stdout
    frozen = path.read_bytes()
    return frozen if strip_sam_cl(frozen) == strip_sam_cl(stdout) else stdout


def refuse_signal(row: dict[str, str], returncode: int) -> None:
    """Stops `freeze` before it writes any file of `row` when NCBI ended its run by a signal (`returncode` < 0, a
    crash). How much of its report NCBI writes before the crash depends on its heap and the environment (the same
    official 2.17.0 binary wrote 349 bytes of db.x2g.bx.of.p0127.nd0.1 when it was frozen and 17 in later runs),
    so the run's stdout is no expectation: such a row takes its stdout from a near row that NCBI runs to its end
    (`near_row`), or has empty stdout cells where LOSAT writes nothing before it refuses."""
    if returncode < 0:
        raise SystemExit(f"{row['row_id']}: NCBI ended its run by signal {-returncode}; a row whose NCBI run "
                         "crashes takes its stdout from a near row that NCBI runs to its end ('ncbi-near: ROW_ID' "
                         "or 'ncbi-near: ROW_ID first N bytes' in covers)")


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


def near_row(row: dict[str, str]) -> tuple[str, int | None] | None:
    """The `ncbi` row whose frozen NCBI run of a near input is the expectation of `row` and the number of its first
    stdout bytes that are the expectation: (ROW_ID, None) for `ncbi-near: ROW_ID` in the covers cell (all of them),
    (ROW_ID, N) for `ncbi-near: ROW_ID first N bytes` (each up to a ` | ` or the end of the cell); None without the
    marker (the row's own frozen files)."""
    _, marker, rest = row["covers"].partition("ncbi-near: ")
    if not marker:
        return None
    words = rest.split(" | ")[0].split()
    if len(words) == 1:
        return words[0], None
    if len(words) == 4 and words[1] == "first" and words[2].isdecimal() and words[3] == "bytes":
        return words[0], int(words[2])
    raise SystemExit(f"{row.get('row_id')}: 'ncbi-near: {rest.split(' | ')[0]}' is neither "
                     "'ncbi-near: ROW_ID' nor 'ncbi-near: ROW_ID first N bytes'")


def near_row_id(row: dict[str, str]) -> str | None:
    """The ROW_ID of `near_row`, or None."""
    near = near_row(row)
    return near[0] if near else None


def near_cut(row: dict[str, str]) -> int | None:
    """The N of `ncbi-near: ROW_ID first N bytes` (`near_row`), or None: all the bytes, or no near row."""
    near = near_row(row)
    return near[1] if near else None


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


def expected_stdout(row: dict[str, str]) -> bytes:
    """The frozen stdout that a `losat-rejection` row with stdout cells expects: the `.out` of `frozen_id`, its first N
    bytes for `ncbi-near: ROW_ID first N bytes` (`near_cut`); empty where the file does not exist."""
    path = FROZEN / f"{frozen_id(row)}.out"
    return path.read_bytes()[:near_cut(row)] if path.exists() else b""


# Copy rows (`ncbi-copies:`): NCBI's run on a query file in which each query that ends a query batch of LOSAT's run
# is followed by a copy of itself, cut to the reports of the original queries.
COPY_PREFIX = "COPY"
COPY_OUTFMTS = ("0", "5", "6", "7")
# NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:86-88, 135-137
BLASTP_BATCH_SIZE = 10000
# A full hit list does not limit a copy row. The composition-based redo ends the remaining matches early only while
# the heap of every query of the batch is full (below), and it takes the matches in the order of their queries
# (blast_hspstream.c:289-319), so the heaps of the queries after the one in work are still empty: in the oracle run no
# match of an original query is ended early, only the copy's can be (the batch's last query), and its report is cut
# out. LOSAT ends no match early with unified P-values (`BlastCompoHeap::without_early_end`).
# NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1565-1566, 1575-1577
# The query chunks of a BLASTP batch: CHUNK_SIZE (else 10000) and OVERLAP_CHUNK_SIZE (else 100) residues.
# NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:59-61, 91-93
# NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:53-57
BLASTP_CHUNK_SIZE = 10000
QUERY_CHUNK_OVERLAP = 100
# The variables that a copy row's oracle run may set in place of the row's own: only those that place the queries in
# batches and chunks (another variable would give the oracle conditions that the row's LOSAT run does not have).
COPY_ENV_NAMES = ("BATCH_SIZE", "CHUNK_SIZE", "OVERLAP_CHUNK_SIZE")


def copies_cell(row: dict[str, str]) -> tuple[str, dict[str, str]] | None:
    """The oracle query file of a copy row (`ncbi-copies: PATH [NAME=VALUE ...]` in its covers cell, up to a ` | `
    or the end; PATH relative to LOSAT/) and the variables of the oracle run that replace the row's own; None for
    another row."""
    _, marker, rest = row["covers"].partition("ncbi-copies: ")
    if not marker:
        return None
    words = rest.split(" | ")[0].split()
    if not words or not all(is_env_token(word) and word.split("=", 1)[0] in COPY_ENV_NAMES for word in words[1:]):
        raise SystemExit(f"{row['row_id']}: ncbi-copies takes a path and NAME=VALUE words of {', '.join(COPY_ENV_NAMES)}")
    return words[0], dict(word.split("=", 1) for word in words[1:])


def fasta_records(path: Path) -> list[tuple[str, str]]:
    """The records of a FASTA file: each defline without its `>`, and its residues without white space."""
    records: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line[1:], ""))
        elif line.strip():
            if not records:
                raise SystemExit(f"{path}: letters before the first defline")
            records[-1] = (records[-1][0], records[-1][1] + "".join(line.split()))
    return records


def record_id(defline: str) -> str:
    """The id of a record: the first word of its defline."""
    return (defline.split() or [""])[0]


def is_copy(name: str | bytes) -> bool:
    """Whether a query id (or a line that starts with it) is a copy's: it starts with `COPY`."""
    return name.startswith(COPY_PREFIX.encode() if isinstance(name, bytes) else COPY_PREFIX)


def query_batches(lengths: list[int], batch_size: int) -> list[list[int]]:
    """NCBI's query batches of a query file (the indexes of its queries): a batch takes queries while the sum of
    their lengths is below `batch_size`; the query that brings it to `batch_size` or more is the batch's last.

    NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140, 163, 170
            """
    batches: list[list[int]] = []
    current: list[int] = []
    size = 0
    for index, length in enumerate(lengths):
        current.append(index)
        size += length
        if size >= batch_size:
            batches.append(current)
            current, size = [], 0
    return batches + [current] if current else batches


def batch_size(env: dict[str, str]) -> int:
    """BLASTP's query batch size in an environment: BATCH_SIZE, else 10000."""
    return int(env.get("BATCH_SIZE", BLASTP_BATCH_SIZE))


def query_chunks(length: int, env: dict[str, str]) -> int:
    """The number of chunks in which NCBI searches a BLASTP query batch of `length` residues in an environment:
    `length` // (CHUNK_SIZE - OVERLAP_CHUNK_SIZE) when CHUNK_SIZE is the larger (each a blank value or none: 10000,
    100), and 1 when that is 1 or less (19,799 residues and fewer by default).

    NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:124-131

           """
    chunk, overlap = (int(env.get(name, "").strip() or default) for name, default in
                      (("CHUNK_SIZE", BLASTP_CHUNK_SIZE), ("OVERLAP_CHUNK_SIZE", QUERY_CHUNK_OVERLAP)))
    return max(length // (chunk - overlap) if chunk > overlap else 0, 1)


def copy_oracle_problems(row: dict[str, str]) -> list[str]:
    """Why the oracle file of a copy row cannot stand for its query file (an empty list when it can): the row is a
    BLASTP `ncbi` row of outfmt 0, 5, 6 or 7 with the default fields and a query file of its own; no query id of the
    row starts with `COPY`; the oracle's records without the copies are the row's records (deflines and residues);
    each copy has a new id and the residues of the original query before it; each batch of the oracle run (its
    BATCH_SIZE: the row's, replaced by the cell's words) ends with a copy and holds no other, and without its copy is
    a batch of LOSAT's run; each batch of both runs is searched in one chunk (`query_chunks`, with the run's
    CHUNK_SIZE and OVERLAP_CHUNK_SIZE); with outfmt 5 the one copy is the last query. The lengths are the residues of
    the FASTA records (`fasta_records`), which are NCBI's lengths for plain FASTA files of letters only."""
    oracle_path, words = copies_cell(row)
    problems = []
    if row["program"] != "blastp" or row["kind"] != "ncbi":
        problems.append("a copy row is a blastp ncbi row")
    if row["outfmt"] not in COPY_OUTFMTS:
        problems.append(f"outfmt {row['outfmt']!r} has no cut rule (0, 5, 6 or 7 with the default fields)")
    if near_row_id(row) or under_valgrind(row) or xinclude.wants_files(row) or xinclude.wants_pssm(row):
        problems.append("a copy row is not ncbi-near, ncbi-valgrind or a row with files")
    if row["query"].startswith(("stdin:", "pssm:", "pssm-gz:")):
        return problems + ["a copy row reads its query from a FASTA file"]
    originals = fasta_records(ENGINE / row["query"])
    oracle = fasta_records(ENGINE / oracle_path)
    if any(is_copy(record_id(defline)) for defline, _ in originals):
        problems.append(f"a query id of {row['query']} starts with {COPY_PREFIX}")
    copies = [index for index, (defline, _) in enumerate(oracle) if is_copy(record_id(defline))]
    if [record for record in oracle if not is_copy(record_id(record[0]))] != originals:
        problems.append(f"the queries of {oracle_path} without the copies are not those of {row['query']}")
    if len({record_id(oracle[index][0]) for index in copies}) != len(copies):
        problems.append("two copies have one id")
    for index in copies:
        if index == 0 or index - 1 in copies or oracle[index][1] != oracle[index - 1][1]:
            problems.append(f"the copy {record_id(oracle[index][0])} does not follow a query with its residues")
    original_index = {index: number for number, index in
                      enumerate(index for index in range(len(oracle)) if index not in copies)}
    env = row_env(row, {})
    own = query_batches([len(residues) for _, residues in originals], batch_size(env))
    batches = query_batches([len(residues) for _, residues in oracle], batch_size({**env, **words}))
    for batch in batches:
        if batch[-1] not in copies or any(index in copies for index in batch[:-1]):
            problems.append(f"the oracle batch of queries {batch} does not end with its one copy")
    if [[original_index[index] for index in batch if index not in copies] for batch in batches] != own:
        problems.append(f"the oracle batches without the copies are not LOSAT's batches {own}")
    for run, records, run_batches, run_env in (("oracle", oracle, batches, {**env, **words}),
                                               ("LOSAT", originals, own, env)):
        for batch in run_batches:
            length = sum(len(records[index][1]) for index in batch)
            if query_chunks(length, run_env) > 1:
                problems.append(f"the {run} batch of queries {batch} ({length} residues) is searched in "
                                f"{query_chunks(length, run_env)} chunks; a copy row's batches are searched in one")
    if row["outfmt"] == "5" and copies != [len(oracle) - 1]:
        problems.append("with outfmt 5 the one copy is the oracle's last query")
    return problems


def copy_count(row: dict[str, str]) -> int:
    """The number of copies in a copy row's oracle file."""
    return sum(is_copy(record_id(defline)) for defline, _ in fasta_records(ENGINE / copies_cell(row)[0]))


def output_lines(data: bytes) -> list[bytes]:
    """The lines of an output, each with its newline (the last without one when the output does not end with it)."""
    parts = data.split(b"\n")
    return [part + b"\n" for part in parts[:-1]] + ([parts[-1]] if parts[-1] else [])


_PROCESSED = re.compile(rb"# BLAST processed (\d+) queries\n")
_QUERY_DEF = b"  <Iteration_query-def>"


def outfmt7_blocks(lines: list[bytes]) -> list[list[bytes]]:
    """The query blocks of an outfmt 7 document (each from its program line, the document's first line, to the next
    one), without the final line `# BLAST processed N queries`, whose N must be their number."""
    final = _PROCESSED.fullmatch(lines[-1]) if lines else None
    if final is None:
        raise ValueError("outfmt 7 without its final line '# BLAST processed N queries'")
    starts = [index for index, line in enumerate(lines[:-1]) if line == lines[0]]
    blocks = [lines[start:end] for start, end in zip(starts, starts[1:] + [len(lines) - 1])]
    if int(final.group(1)) != len(blocks) or any(not block[1:2] or not block[1].startswith(b"# Query: ")
                                                 for block in blocks):
        raise ValueError(f"outfmt 7: {len(blocks)} blocks with a '# Query: ' line, final line {lines[-1]!r}")
    return blocks


def outfmt0_reports(lines: list[bytes]) -> list[tuple[int, int]]:
    """The line ranges of the query reports of an outfmt 0 document: each from its `Query= ` line to the next one
    or, for the last, to the document's final `  Database: ` block."""
    starts = [index for index, line in enumerate(lines) if line.startswith(b"Query= ")]
    footers = [index for index, line in enumerate(lines) if line.startswith(b"  Database: ")]
    if not starts or not footers or footers[-1] < starts[-1]:
        raise ValueError("outfmt 0 without 'Query= ' lines before a final '  Database: ' block")
    return list(zip(starts, starts[1:] + [footers[-1]]))


def outfmt5_iterations(lines: list[bytes]) -> list[tuple[int, int]]:
    """The line ranges of the `<Iteration>` elements of an outfmt 5 document, each with its end line."""
    starts = [index for index, line in enumerate(lines) if line == b"<Iteration>\n"]
    ends = [index + 1 for index, line in enumerate(lines) if line == b"</Iteration>\n"]
    if not starts or len(starts) != len(ends) or any(not start < end for start, end in zip(starts, ends)):
        raise ValueError("outfmt 5 without whole <Iteration> elements")
    return list(zip(starts, ends))


def iteration_is_copy(lines: list[bytes]) -> bool:
    """Whether an outfmt 5 `<Iteration>` element (its lines) is a copy's: its `<Iteration_query-def>` starts with
    `COPY` (an element without one is refused)."""
    names = [line[len(_QUERY_DEF):] for line in lines if line.startswith(_QUERY_DEF)]
    if len(names) != 1:
        raise ValueError("an outfmt 5 <Iteration> without one <Iteration_query-def>")
    return is_copy(names[0])


def cut_copies(data: bytes, outfmt: str, copies: int) -> bytes:
    """NCBI's whole output `data` of a copy row's oracle run without the reports of its `copies` copies (queries
    whose id starts with `COPY`), by the rule of the output format (`copy_oracle_problems` allows 0, 5, 6 and 7):

    - 6: the lines whose first field starts with `COPY` are left out.
    - 7: each copy's block (`outfmt7_blocks`) is left out, and the final line `# BLAST processed N queries` gets
      the number of the blocks kept: the one byte of NCBI's that a rule writes.
    - 0: each copy's report (`outfmt0_reports`: its `Query= COPY` line up to the next `Query= ` line or, for the
      last query, up to the final `  Database: ` block) is left out.
    - 5: the last `<Iteration>` element, which must be the copy's and the only copy's, is left out.

    ValueError where the output does not have that shape or another number of copies."""
    lines = output_lines(data)
    if outfmt == "6":
        return b"".join(line for line in lines if not is_copy(line.split(b"\t", 1)[0]))
    if outfmt == "7":
        blocks = outfmt7_blocks(lines)
        kept = [block for block in blocks if not is_copy(block[1][len(b"# Query: "):])]
        if len(blocks) - len(kept) != copies:
            raise ValueError(f"outfmt 7: {len(blocks) - len(kept)} copy blocks for {copies} copies")
        return b"".join(line for block in kept for line in block) + b"# BLAST processed %d queries\n" % len(kept)
    if outfmt == "0":
        dropped = [(start, end) for start, end in outfmt0_reports(lines) if is_copy(lines[start][len(b"Query= "):])]
        if len(dropped) != copies:
            raise ValueError(f"outfmt 0: {len(dropped)} copy reports for {copies} copies")
        return b"".join(line for index, line in enumerate(lines)
                        if not any(start <= index < end for start, end in dropped))
    if outfmt == "5":
        iterations = outfmt5_iterations(lines)
        flags = [iteration_is_copy(lines[start:end]) for start, end in iterations]
        if copies != 1 or flags != [False] * (len(flags) - 1) + [True]:
            raise ValueError(f"outfmt 5: the copy must be the last of {len(flags)} iterations and the only one")
        start, end = iterations[-1]
        return b"".join(lines[:start] + lines[end:])
    raise ValueError(f"no cut rule for outfmt {outfmt!r}")


def subjects_per_query(data: bytes, outfmt: str) -> list[int]:
    """The number of subjects in each query's report of an output of outfmt 0, 5, 6 or 7 (the copies' included):
    the distinct second fields of a query's rows (6, 7), the lines of its description table (0), its `<Hit>`
    elements (5)."""
    lines = output_lines(data)
    if outfmt in ("6", "7"):
        blocks = ([[line for line in block if not line.startswith(b"#")] for block in outfmt7_blocks(lines)]
                  if outfmt == "7" else [lines])
        counts: dict[tuple[int, bytes], set[bytes]] = {}
        for number, block in enumerate(blocks):
            for line in block:
                fields = line.split(b"\t")
                counts.setdefault((number, fields[0]), set()).add(fields[1] if len(fields) > 1 else b"")
        return [len(subjects) for subjects in counts.values()]
    if outfmt == "0":
        result = []
        for start, end in outfmt0_reports(lines):
            table = [index for index in range(start, end)
                     if lines[index].startswith(b"Sequences producing significant alignments:")]
            count = 0
            if table:
                index = table[0] + 2
                while index < end and lines[index].strip():
                    count, index = count + 1, index + 1
            result.append(count)
        return result
    if outfmt == "5":
        return [sum(line == b"<Hit>\n" for line in lines[start:end]) for start, end in outfmt5_iterations(lines)]
    raise ValueError(f"no subject count for outfmt {outfmt!r}")


def copy_row_problems(row: dict[str, str]) -> list[str]:
    """The problems of a copy row's frozen files (`check` and the manifest's test): those of its oracle file
    (`copy_oracle_problems`), and a cut of the frozen whole output <row_id>.copies.out that is not <row_id>.out."""
    problems = copy_oracle_problems(row)
    whole = FROZEN / f"{row['row_id']}.copies.out"
    if not whole.is_file():
        return problems + [f"no frozen {whole.name}"]
    try:
        cut = cut_copies(whole.read_bytes(), row["outfmt"], copy_count(row))
    except ValueError as error:
        return problems + [f"{whole.name}: {error}"]
    expected = FROZEN / f"{row['row_id']}.out"
    if not expected.is_file() or cut != expected.read_bytes():
        problems.append(f"the cut of {whole.name} is not {expected.name}")
    return problems


def freeze_copies(row: dict[str, str], ncbi_bin: str, env: dict[str, str]) -> None:
    """`freeze` of a copy row: NCBI `-num_threads 1` on the oracle file with the oracle's variables (exit 0 and a
    stderr that names no copy, else the freeze stops), its whole stdout as <row_id>.copies.out, its cut
    (`cut_copies`) as <row_id>.out and its stderr, a warning of the arguments that NCBI writes for the row's own
    query file too, as <row_id>.err, with the row's cells (rc 0); it prints the number of subjects of each query of
    the whole output (`subjects_per_query`, the copies' included), which shows a full hit list."""
    problems = copy_oracle_problems(row)
    if problems:
        raise SystemExit(f"{row['row_id']}: " + "; ".join(problems))
    oracle_path, words = copies_cell(row)
    command = [str(Path(ncbi_bin) / row["program"]), *search_argv({**row, "query": oracle_path}), "-num_threads", "1"]
    run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env={**row_env(row, env), **words})
    if run.returncode != 0 or COPY_PREFIX.encode() in run.stderr:
        raise SystemExit(f"{row['row_id']}: NCBI's oracle run exited {run.returncode} with stderr "
                         f"{run.stderr.decode(errors='replace')[:200]!r}; a copy row needs exit 0 and a stderr that "
                         f"names no copy")
    try:
        out = cut_copies(run.stdout, row["outfmt"], copy_count(row))
    except ValueError as error:
        raise SystemExit(f"{row['row_id']}: {error}") from None
    # The copies' own reports are those of a batch's last query, which NCBI computes with bytes past its array
    # that change from run to run: a frozen whole output whose cut is this run's cut is kept as it is.
    whole = FROZEN / f"{row['row_id']}.copies.out"
    try:
        kept = whole.is_file() and cut_copies(whole.read_bytes(), row["outfmt"], copy_count(row)) == out
    except ValueError:
        kept = False
    if not kept:
        whole.write_bytes(run.stdout)
    (FROZEN / f"{row['row_id']}.out").write_bytes(out)
    if run.stderr:
        (FROZEN / f"{row['row_id']}.err").write_bytes(run.stderr)
    else:
        (FROZEN / f"{row['row_id']}.err").unlink(missing_ok=True)
    row.update(rc="0", stdout_sha256=sha256(out), stdout_bytes=str(len(out)),
               stderr_sha256=sha256(run.stderr) if run.stderr else "")
    print(f"{row['row_id']}\tcopies={oracle_path}\twhole={len(run.stdout)}{' (frozen whole kept)' if kept else ''}"
          f"\tstdout={len(out)}\tsubjects={subjects_per_query(run.stdout, row['outfmt'])}")


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
        if row["kind"] == "losat-rejection" and row["stdout_bytes"] and not near_row_id(row):
            # NCBI's stdout of the command that LOSAT refuses after the reports NCBI writes too.
            if xinclude.wants_pssm(row):
                run = xinclude.run_with_pssm(command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE,
                                             timeout=TIMEOUT, env=row_env(row, env), preexec_fn=ncbi_stack(row)).run
            elif rps_cell(row["db"]) is not None:
                # An RPS row: NCBI's run with BLASTDB the inflated database (`freeze_rps`).
                if rps_dir_token(row):
                    raise SystemExit(f"{row['row_id']}: a rejection row with NCBI's stdout writes no {{db_dir}}")
                db_dir = inflate_rps(row["db"], False)
                try:
                    run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env=rps_env(row, env, db_dir))
                finally:
                    shutil.rmtree(db_dir, ignore_errors=True)
            else:
                run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env=row_env(row, env),
                                 preexec_fn=ncbi_stack(row))
            refuse_signal(row, run.returncode)
            out = frozen_stdout(row, run.stdout)
            (FROZEN / f"{row['row_id']}.out").write_bytes(out)
            row.update(stdout_sha256=sha256(out), stdout_bytes=str(len(out)))
            print(f"{row['row_id']}\tncbi_rc={run.returncode}\tstdout={len(out)} (rejection row: stdout only)")
            continue
        if row["kind"] != "ncbi" or near_row_id(row):
            continue
        if copies_cell(row) is not None:
            freeze_copies(row, args.ncbi_bin, env)
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
            refuse_signal(row, run.returncode)
            out = run.stdout
            xinclude.freeze(written, FROZEN, row["row_id"], row["outfmt"])
        elif pssm:
            pssms = xinclude.run_with_pssm(command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE,
                                           timeout=TIMEOUT, env=row_env(row, env), preexec_fn=ncbi_stack(row))
            run = pssms.run
            refuse_valgrind_messages(row, run.stderr)
            refuse_signal(row, run.returncode)
            xinclude.freeze_pssm(pssms, FROZEN, row["row_id"], xinclude.files_kind(row))
            out = frozen_stdout(row, run.stdout)
            (FROZEN / f"{row['row_id']}.out").write_bytes(out)
        else:
            run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env=row_env(row, env),
                             preexec_fn=ncbi_stack(row))
            refuse_valgrind_messages(row, run.stderr)
            refuse_signal(row, run.returncode)
            out = frozen_stdout(row, run.stdout)
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
    # The rows whose expectation is a near input's NCBI run (`near_row`) take that row's columns, after it is
    # frozen (it may come later in the file); a rejection row takes its stdout columns only.
    by_id = {row["row_id"]: row for row in rows}
    for row in rows:
        near = near_row_id(row)
        if not selected(row, args) or row["kind"] not in ("ncbi", "losat-rejection") or not near:
            continue
        source = by_id.get(near)
        if source is None or source["kind"] != "ncbi" or near_row_id(source) or source["outfmt"] != row["outfmt"]:
            raise SystemExit(f"{row['row_id']}: ncbi-near {near} is not an ncbi row of its own run with the same outfmt")
        cut = near_cut(row)
        if row["kind"] == "ncbi":
            if cut is not None:
                raise SystemExit(f"{row['row_id']}: an ncbi row takes all of its near row's bytes (no 'first N bytes')")
            row.update({key: source[key] for key in ("rc", "stdout_sha256", "stdout_bytes", "stderr_sha256")})
        else:
            # A rejection row keeps its own rc and stderr (LOSAT's) and takes the near run's stdout, or its first
            # N bytes, which the near run must still write.
            out = (FROZEN / f"{near}.out").read_bytes()
            if cut is not None and cut > len(out):
                raise SystemExit(f"{row['row_id']}: ncbi-near {near} wrote {len(out)} bytes, fewer than its first {cut}")
            out = out[:cut]
            row.update(stdout_sha256=sha256(out), stdout_bytes=str(len(out)))
        print(f"{row['row_id']}\tncbi-near={near}" + (f" first {cut} bytes" if cut is not None else "")
              + f"\trc={row['rc']}\tstdout={row['stdout_bytes']}")
    write_manifest(comments, rows)
    return 0


def freeze_rps(row: dict[str, str], command: list[str], env: dict[str, str]) -> None:
    """`freeze` of an rpsblast or rpstblastn row: NCBI's run with BLASTDB the inflated database, its
    directory written as `{db_dir}` where the row says so; a row of outfmt 13 or 14 with its files
    (`xinclude.wants_files`) freezes the master and the side files as other programs' rows."""
    files = xinclude.wants_files(row)
    if files and rps_dir_token(row):
        raise SystemExit(f"{row['row_id']}: a row with files writes no {{db_dir}}")
    db_dir = inflate_rps(row["db"], rps_dir_token(row))
    run_env = rps_env(row, env, db_dir)
    try:
        if files:
            written = xinclude.run_with_files(command, row["row_id"], row["outfmt"], cwd=ENGINE,
                                              timeout=TIMEOUT, env=run_env)
            run = written.run
        else:
            run = run_search(command, row, cwd=ENGINE, timeout=TIMEOUT, env=run_env)
    finally:
        shutil.rmtree(db_dir, ignore_errors=True)
    refuse_signal(row, run.returncode)
    stdout, stderr = run.stdout, run.stderr
    if rps_dir_token(row):
        stdout = stdout.replace(str(db_dir).encode(), b"{db_dir}")
        stderr = stderr.replace(str(db_dir).encode(), b"{db_dir}")
    if files:
        xinclude.freeze(written, FROZEN, row["row_id"], row["outfmt"])
    else:
        (FROZEN / f"{row['row_id']}.out").write_bytes(stdout)
    err = FROZEN / f"{row['row_id']}.err"
    if stderr:
        err.write_bytes(stderr)
    else:
        err.unlink(missing_ok=True)
    row.update(rc=str(run.returncode), stdout_sha256=sha256(stdout), stdout_bytes=str(len(stdout)),
               stderr_sha256=sha256(stderr) if stderr else "")
    if files:
        row.update(xinclude.observed(written))
    print(f"{row['row_id']}\trc={run.returncode}\tstdout={len(stdout)}\tstderr={len(stderr)}"
          + (f"\tsides={len(written.sides)}" if files else ""))


def status_differs(observed: int, expected: int, library: bool) -> bool:
    """Whether a run's exit status is not the frozen one. `examples/library_command_line.rs` ends every
    failure with status 1: the library API's errors carry the command line's message, not its status."""
    return (observed != 0) != (expected != 0) if library else observed != expected


def check_row(row: dict[str, str], losat: Path, threads: int, env: dict[str, str],
              library: bool = False) -> tuple[str, str]:
    """The row's status (`same`, `rejection`, `differs`, `library-skip`, `library-input-error`) and its
    detail. Standard output, standard error, the exit status and the PSSM files are compared byte for
    byte."""
    files = xinclude.wants_files(row)
    pssm = xinclude.wants_pssm(row)
    # The library API refuses rpsblast and rpstblastn (`api::Program::parse`, line PS items PS11a, PS12a).
    # Its query is records or a PSSM file (`Inputs::new`, `Inputs::pssm`), never both: a row that gives
    # the command line `-query` and `-in_pssm` (the parser's exclusion) has no library call, nor has
    # a row that gives `-in_msa` with either of them, or the alignment's master options
    # (`xinclude.ALIGNMENT_OPTIONS`) without `-in_msa` (`Inputs::msa` is the library's alignment
    # query, line LIB item 1.6), nor a row that writes the search strategy file
    # (`-export_search_strategy`: the library writes none and refuses it, `api::search::STRATEGY_FILES`,
    # line RG item F1b-1).
    argv = search_argv(row)
    if library and (row["program"] in ("rpsblast", "rpstblastn")
                    or ("-query" in argv and "-in_pssm" in argv)
                    or ("-in_msa" in argv and ("-query" in argv or "-in_pssm" in argv))
                    or ("-in_msa" not in argv
                        and any(word in xinclude.ALIGNMENT_OPTIONS for word in argv))
                    or "-export_search_strategy" in argv
                    or query_stdin(row)):
        return "library-skip", ""
    if (files or pssm) and query_stdin(row):
        raise SystemExit(f"{row['row_id']}: a row with files reads no query from standard input")
    command = [str(losat), row["program"], *search_argv(row), "-num_threads", str(threads)]
    # An RPS row runs with BLASTDB the inflated database (`inflate_rps`), its files rows as any other.
    db_dir = inflate_rps(row["db"], rps_dir_token(row)) if rps_cell(row["db"]) is not None else None

    def in_losat_words(ncbi: bytes, output: bool = True) -> bytes:
        """NCBI's frozen bytes in LOSAT's words (AGENTS.md rule 8, `losat_wording`). A row that keeps the
        database's directory as `{db_dir}` is mapped with this run's directory in the token's place, as
        NCBI's run wrote it: the mapping re-wraps an ASN.1 text string by its length."""
        words = losat_wording.losat_output_bytes if output else losat_wording.losat_bytes
        if db_dir is None or not rps_dir_token(row):
            return words(ncbi)
        path = str(db_dir).encode()
        return words(ncbi.replace(b"{db_dir}", path)).replace(path, b"{db_dir}")

    run_env = row_env(row, env) if db_dir is None else rps_env(row, env, db_dir)
    try:
        if files:
            written = xinclude.run_with_files(command, row["row_id"], row["outfmt"], cwd=ENGINE,
                                              env=run_env, timeout=TIMEOUT)
            run = written.run
        elif pssm:
            pssms = xinclude.run_with_pssm(command, row["row_id"], xinclude.files_kind(row), cwd=ENGINE,
                                           env=run_env, timeout=TIMEOUT)
            run = pssms.run
        else:
            run = run_search(command, row, cwd=ENGINE, env=run_env, timeout=TIMEOUT)
    finally:
        if db_dir is not None:
            shutil.rmtree(db_dir, ignore_errors=True)
    if db_dir is not None and rps_dir_token(row):
        run.stdout = run.stdout.replace(str(db_dir).encode(), b"{db_dir}")
        run.stderr = run.stderr.replace(str(db_dir).encode(), b"{db_dir}")
    if row["kind"] == "losat-rejection":
        # The library reads the whole query before it searches: a reader error ends the call before
        # the batch that the command line refuses, with the reader's message.
        if library and run.returncode != 0 and run.stdout == b"" \
                and run.stderr.startswith(losat_wording.QUERY_ERROR) and expected_text(row).encode() not in run.stderr:
            return "library-input-error", ""
        problems = []
        if status_differs(run.returncode, int(row["rc"]), library):
            problems.append(f"exit {run.returncode}, expected {row['rc']}")
        # A refusal follows no byte in either path (the library holds its output where the command
        # line holds it, BLASTX's XML documents through `blast_xml::documents_may_refuse`); the rows
        # whose refusal comes after NCBI's reports are psiblast's PSSM file and rounds rows, which the
        # library skips, and the refusals at a later query batch (`ncbi-near:`), which the library
        # writes as the command line does.
        if row["stdout_bytes"]:
            out = expected_stdout(row)
            if sha256(out) != row["stdout_sha256"] or len(out) != int(row["stdout_bytes"]):
                problems.append("frozen .out does not match the manifest")
            if run.stdout != in_losat_words(out):
                problems.append("stdout " + first_difference(in_losat_words(out), run.stdout))
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
    expected_out_losat = in_losat_words(out)
    expected_err_losat = in_losat_words(err, output=False)
    # The library reads the whole query before it searches: a reader error ends the call there,
    # with the command line's message of it, before the reports of the batches before it.
    if library and run.returncode != 0 and int(row["rc"]) != 0 and run.stdout == b"" \
            and run.stderr.startswith(losat_wording.QUERY_ERROR) and run.stderr in expected_err_losat:
        return "library-input-error", ""
    # A copy row's expectation is the cut of its frozen whole NCBI output (`copy_row_problems`).
    problems = copy_row_problems(row) if copies_cell(row) is not None else []
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
    if pssm:
        problems += xinclude.pssm_problems(pssms, row, FROZEN, first_difference)
    if problems:
        return "differs", "; ".join(problems)
    return "same", ""


def check(args: argparse.Namespace) -> int:
    home = tempfile.TemporaryDirectory(prefix="db_search_home_")
    env = clean_env(dict(os.environ), home.name)
    losat = args.losat.resolve()
    lines, bad = [], 0
    for row in read_manifest()[1]:
        if not selected(row, args):
            continue
        status, detail = check_row(row, losat, args.threads, env, args.library)
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
    args = parser.parse_args()
    return freeze(args) if args.command == "freeze" else check(args)


if __name__ == "__main__":
    sys.exit(main())

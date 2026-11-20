# LOSAT -db fixture rows, one file per program and row family files

`<program>.tsv` holds the rows of one program (`blastp blastn tblastn blastx tblastx psiblast rpsblast rpstblastn`; the `program`
column decides the file), but those of a row family with its own file: `<program>.<family>.tsv` holds the program's rows
`db.<family>.*` (the second word of the row id; for example `blastp.ldo.tsv`, the BLASTP rows `db.ldo.*`). Each row has
one file, decided by its `program` cell and its id. `LOSAT/tests/db_search_fixtures.py` is the one loader: `read_manifest`
reads the files in that program order, each program's file before its family files (by name), each in its own row order,
and refuses a row in another file than its own; `freeze` writes a row back to its file. Every file starts with the same
header line, after its `# Rows ...` comment lines (a comment that describes rows of several programs is repeated in each
of their files). `LOSAT/tests/check_manifest_sizes.py` fails CI when a file passes 512 KiB: move a row family of it, its
rows and their `# Rows` comment lines, to its own file (the loader needs no change).

Columns: row_id, program, kind, db, query, extra_args, outfmt, rc, stdout_sha256, stdout_bytes, stderr_sha256, covers,
and for the rows of outfmt 13 and 14 and psiblast's PSSM files two more, files and side_sha256 (a row without them has
12 cells; `xinclude_fixtures.py` says what they mean).

- `kind ncbi`: stdout, stderr and exit status (rc) equal the frozen bytes of NCBI BLAST+ 2.17.0
  `PROGRAM -db ... -num_threads 1` (`LOSAT/tests/fixtures/db_search/<row_id>.out`, `.err` when stderr is not empty).
- `kind ncbi` with `ncbi-near: ROW_ID` in covers (after a ` | `): the expectation is the frozen bytes and columns of the
  `ncbi` row ROW_ID, NCBI's run of a near input (the mask one letter shorter; for exception 3 the same volume with zero
  bytes after its last memory page), where NCBI's own run of the row reaches an approved NCBI defect
  (PD-LOSAT-NCBI-DEFECTS exceptions 13 and 3); `freeze` copies ROW_ID's columns and runs nothing.
- `kind ncbi` with a part ` | ncbi-valgrind` of covers: the frozen bytes are NCBI's run under valgrind memcheck, where
  NCBI's native run reads past its buffers (PD-LOSAT-NCBI-DEFECTS exception 3) and may crash or vary; `freeze` runs
  such a row only with `--valgrind CMD` (for example `--valgrind "valgrind -q --log-file=/dev/null"`) and otherwise
  keeps its columns and files.
- `kind ncbi` with `ncbi-stack: unlimited` in covers: `freeze` runs NCBI with no stack limit, where its default 8 MiB
  stack overflows on a PSSM's values nested thousands of levels deep (PD-LOSAT-NCBI-DEFECTS exception 15); `check`
  runs LOSAT as for any row.
- `kind losat-rejection`: LOSAT exits with rc, stdout is empty, stderr contains the text after `expect: ` in covers
  (to a ` | ` or the end); the hash columns are empty.
- `db`: path relative to LOSAT/ (an empty cell is `-db ""`; names separated by spaces are one `-db` argument;
  `subject:PATH` runs `-subject PATH` in place of `-db`, in any program: the rows of a subject set, such as psiblast's
  losat-rejection row of a subject set whose sequences all have no letters); `extra_args` is split like a shell line;
  `outfmt` is one argument; `stdout_sha256`, `stdout_bytes` and `stderr_sha256` are written by `freeze`.
- A row runs `LOSAT PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT` from `LOSAT/` (from `LOSAT/DIR` with
  `cwd=DIR`, below).
- `query`: a path relative to LOSAT/; `pssm-gz:PATH` or `pssm:PATH` is a PSSM query (`-in_pssm`); `stdin:pipe:PATH`,
  `stdin:file:PATH` or `stdin:empty` is `-query -` with PATH's bytes through a pipe, the file PATH itself, or an empty pipe
  as standard input (`db_search_fixtures.query_stdin`; such a row has no PSSM or XInclude files and the library skips it).

## Tokens of extra_args that are not arguments

- `NAME=VALUE` with NAME one of `db_search_fixtures.ENV_NAMES`: an environment variable of both programs (NCBI's run
  of `freeze` and LOSAT's of `check`). `BLASTDB=DIR` (the directory of a `taxdb.bti`) and `NCBI=DIR` (the directory of
  `$NCBI/.ncbirc` and of the formatter's `$NCBI/<name>`) are paths relative to LOSAT/; `BATCH_SIZE`, `CHUNK_SIZE`,
  `OVERLAP_CHUNK_SIZE`, `ARCHIVE_FORMAT`, `PRE_FETCH_SEQS_LIMIT`, `BL2SEQ_LEGACY`, `OLD_FSC`, `ADAPTIVE_CBS`, and the
  `-html` formatter's `FMTCFG` (the name of its configuration file in the run directory, default `.ncbirc`),
  `GETURL_DEBUG` (a `REGISTRY: <file>` line on stderr) and `NCBI_DONT_USE_NCBIRC` are values as they are (NCBI
  reference: `align_format_util.cpp:3383-3415`, `CAlignFormatUtil::InitConfig`). The caller's own variables of these
  names (and every NCBI* and BLAST*) are dropped first (`clean_env`).
- `cwd=DIR` (lower case: neither an argument nor a variable): the row runs in `LOSAT/DIR` in place of `LOSAT/`, in
  `freeze` and in `check`. DIR, relative to LOSAT/, holds the row's `.ncbirc` and the other files NCBI reads from its
  current directory (the file that FMTCFG names). The row's cells keep every path relative to LOSAT/; the programs get
  them relative to DIR (`os.path.relpath`): each name of the `db` cell, the PATH of `subject:PATH`, the query, and the
  values of `BLASTDB=` and `NCBI=`. So the argv's text is the same on every machine: NCBI writes the `-db` text into
  the `[LOCAL_ID] TOOL_URL` link (`db_name`), and an absolute path would put the checkout's path into the frozen bytes.
  Not taken: absolute paths (bytes of one machine), paths relative to DIR in the cells (every reader counts on paths
  relative to LOSAT/), a temporary directory with copies of the files (absolute paths again), a new column (a new
  header for every file and every reader). A `cwd=` row with a PSSM or standard-input query, an `rps-gz:` db cell, a
  `files` column or a path under `tests/` in extra_args (an argument that would not be rebased) is refused. `freeze`
  and `check` refuse to start when `LOSAT/.ncbirc` exists (it would change every row).

## Who runs which rows

- `db_search_fixtures.py check` (and `freeze`): every row, the `cwd=` and `NAME=VALUE` rows included;
  `ci_fast_regressions.py` runs it for every program at 1 and 4 threads, and the `losat-gates` skill uses the same
  checkers.
- `db_search_fixtures.py check --library` (`examples/library_command_line`, the library API): `library-skip` for the
  rows of rpsblast and rpstblastn, those with both `-query` and `-in_pssm`, those with a standard-input query, and those
  with `cwd=` or a token of `PROCESS_SETTINGS` (FMTCFG, NCBI, GETURL_DEBUG, NCBI_DONT_USE_NCBIRC): the API's contract
  (`LOSAT/src/api/mod.rs:78-86`) names what a call reads of its process, and these formatter settings are not among
  them; `-html` reaches the library's host, LOSAT Web, in T4b (line RG, decision T4-1).
- `check_wasm_threading.py`: a few rows by id, every psiblast row and blastx's rows of outfmt 5 and 13 to 16; no row of
  outfmt 0 or 18 (the `-html` rows) and no psiblast `-html` row (psiblast refuses `-html`).
- The Rust tests that read this directory read named program files (`blastp.tsv`, `psiblast.tsv`, `blastx.tsv`,
  `tblastn.tsv`, `tblastx.tsv`) and select their own rows, never a family file such as `<program>.html.tsv`;
  `check_manifest_sizes.py` checks the size of every file.

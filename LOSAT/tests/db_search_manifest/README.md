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
  `ncbi` row ROW_ID, NCBI's run of a near input (the mask one letter shorter, in TBLASTN and TBLASTX the fewest bases
  after which every reading frame keeps a letter after its part, for exception 21 with the freed letters stored as `X`
  or `N`; for exception 3 the same volume with zero bytes after its last memory page; for exception 19 the hit list
  size 1073741798; for exception 27 the database with the short mask dropped and its bases stored as `N`; for
  exception 31 a PSSM file's gapped Lambda of 8e-9, in outfmt 7, which prints no Lambda), where NCBI's own run of the
  row reaches an approved NCBI defect (PD-LOSAT-NCBI-DEFECTS exceptions 13, 21, 3, 19, 27 and 31); `freeze` copies
  ROW_ID's columns and runs nothing.
- `kind ncbi` with a part ` | ncbi-valgrind` of covers: the frozen bytes are NCBI's run under valgrind memcheck, where
  NCBI's native run reads past its buffers (PD-LOSAT-NCBI-DEFECTS exception 3) and may crash or vary; `freeze` runs
  such a row only with `--valgrind CMD` (for example `--valgrind "valgrind -q --log-file=/dev/null"`) and otherwise
  keeps its columns and files.
- `kind ncbi` with `ncbi-stack: unlimited` in covers: `freeze` runs NCBI with no stack limit, where its default 8 MiB
  stack overflows on a PSSM's values nested thousands of levels deep (PD-LOSAT-NCBI-DEFECTS exception 15); `check`
  runs LOSAT as for any row.
- `kind ncbi` with `device: MODE` in covers (Linux only; line RG item T6a fix 1; the rows `db.device.*` of
  `<program>.device.tsv`): the report goes to a device that fails, and the expectation is NCBI's exit status and
  stderr and the bytes that reached the device (`<row_id>.out`). MODE `stdout-full`: standard output is /dev/full;
  `out-full`: `-out /dev/full`; `fsize=N`: `-out` a temporary file under `RLIMIT_FSIZE` N, SIGXFSZ ignored (the
  write past N bytes fails with EFBIG); `check` gives LOSAT's file N plus the bytes that LOSAT's words add to NCBI's
  frozen bytes (`device_growth`: LOSAT's version statement moves the same point of the report that far). Inputs: each program's base search (blastn nslice, blastp prot_v5, tblastn
  tslice, blastx bslice, tblastx xslice; formats 0, 5-12, 15-18 where LOSAT dev and NCBI agree without a device; SAM
  only for blastn) with stdout-full, out-full and `fsize` of half the report (at most 100000) and of its size minus 1
  (where the report is at most 100001 bytes), and for outfmt 0 inside the last query footer, and just inside and 20
  bytes after the epilog's flushed line feeds; `dev_nw.fna` (a warning for an all-N query between two queries, also
  with BATCH_SIZE=1; outfmt 6, 7, 10); `dev_pblk.faa` (1000 copies of p3.faa's first record, 61,000 letters) against
  `-subject dev_psubj.faa` (prot_v5 as FASTA; outfmt 0, 5, 8, 9, 12 with `fsize` 100000); and outfmt 16 with
  `-subject` (`dev_nsubj.fna`, `dev_psubj.faa`; out-full and half). Frozen files stay at or below 100 KB.
  Fix 3 (2026-11-22) adds the modes `stdout-fsize=N` (standard output a file limited to N bytes, no `-out`),
  `stdout-stderr-fsize=N` (standard output and standard error one such file, `> file 2>&1`), `stderr-full`
  (standard error /dev/full), and `stdout-closed` / `stderr-closed` (a pipe whose reader closed: LOSAT only, NCBI is
  ended by SIGPIPE; the expectation is the `ncbi-near` row of the full device, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES
  exception 5), and rows (`<program>.device.tsv`, psiblast, rpsblast and rpstblastn files new): psiblast
  `Search has CONVERGED!` limits (`paja`, `tk.faa` 2 rounds and `tk_o.faa` 3 rounds; outfmt 0, 6, 7) and its warnings
  on a failed standard error; rpsblast and rpstblastn on each device; tblastx and rpstblastn outfmt 6 at 4,096,
  8,191 and 8,192 bytes on a file and on standard output; the review's outfmt 0 exit transitions on standard output
  (`dev_nm.fna` against `of18/ns_ids` with `-task blastn`, `dev_pm.faa` against `of18/ps_ids`), each boundary on
  both sides; `dev_nw.fna`'s warning on a failed standard error and in one limited file with the report; a strategy
  file that LOSAT refuses on a device that fails at NCBI's prolog, in every program; and closed pipes for formats 0, 5, 6, 11 (and
  blastn's 17).
  Line RG item T8f (2026-11-28) adds the rows `db.device.*.apost.*` (out-full, outfmt 11, every search program;
  `tests/db_search/apost`, README.md there): a query or subject title with a byte that the archive's stream replaces,
  whose 803.7 post the error archive holds where the archive's write throws (with a strategy, over two query batches,
  and over psiblast's queries and rounds, where an archive that was written drops it).
- `kind ncbi` with `ncbi-copies: PATH [NAME=VALUE ...]` in covers (a copy row; up to a ` | ` or the end): the
  expectation is made from NCBI's run on the oracle query file PATH (relative to LOSAT/) in place of the row's query
  file. With unified P-values (`blastp -comp_based_stats 1u`, `2u`, `3u`) NCBI BLAST+ 2.17.0 adjusts a query's E-values
  with the values of the next query of its batch and, for a batch's last query, with bytes past its array that change
  from run to run. In PATH each query that ends a query batch of LOSAT's run is followed by a copy of itself (its
  residues under an id that starts with `COPY`: `>COPY`, `>COPY2`), and only the copies end the oracle run's batches
  (a batch takes queries while the sum of their lengths is below BATCH_SIZE, 10000 for BLASTP; the query that brings it
  there ends it); so each original's report is NCBI's with its own values. The `NAME=VALUE` words set the oracle run's
  environment in place of the row's words of those names (BATCH_SIZE). Each batch of both runs is searched in one
  chunk: NCBI splits a batch of L residues into L / (CHUNK_SIZE - OVERLAP_CHUNK_SIZE) chunks when that is 2 or more
  (10000 and 100 by default, so from 19,800 residues), the copy would move the split, and the split can change the
  preliminary hits. The loader checks both files (`copy_oracle_problems`). `freeze` runs NCBI `-num_threads 1` on PATH
  (exit 0 and a stderr that names no copy, or it stops; such a stderr is a warning of the arguments that NCBI writes for
  the row's own query file too, for example the advice to examine 5 or more matches under `-max_target_seqs 1`, and
  is frozen as `<row_id>.err`), keeps the whole stdout unchanged as `<row_id>.copies.out` and writes as `<row_id>.out`
  (and the row's cells; rc 0) the cut of it by a fixed rule (`cut_copies`):
  - outfmt 6: the lines whose first field (the query id) starts with `COPY` are left out;
  - outfmt 7: each copy's block (its `# BLASTP 2.17.0+` line up to the next such line or the final line) is left out,
    and the final line `# BLAST processed N queries` gets N lowered by the number of copies: **the one place where a
    byte of NCBI's output is rewritten**;
  - outfmt 0: each copy's report (its `Query= COPY` line up to the next `Query= ` line or, for the last query, up to
    the final `  Database: ` block) is left out;
  - outfmt 5: the copy's `<Iteration>` element is left out; the one copy must be the oracle file's last query (the
    iteration numbers of later queries would change);
  - any other outfmt, a custom field list, another program or kind: the loader refuses the row.

  The copies' own reports are those of a batch's last query and change from run to run; a later `freeze` keeps the
  frozen `<row_id>.copies.out` when its cut is the new run's cut (so a freeze of an unchanged expectation leaves the
  files as they are).

  A hit list that fills does not limit a copy row (`-max_target_seqs`, `-num_descriptions` and `-num_alignments` are
  allowed). NCBI's composition-based redo takes the matches in the order of their queries and ends the remaining ones
  early (`BlastCompo_EarlyTermination`) only while the list of every query of the batch is full; the lists of the
  queries after the one in work are still empty. So in the oracle run no match of an original query is ended early,
  only the copy's matches can be, and the copy's report is cut out. LOSAT ends no match early with unified P-values,
  so the cut is its expectation also where a list fills. `freeze` prints the number of subjects of each query of the
  oracle output. `check` cuts `<row_id>.copies.out` again and fails the row when that is not `<row_id>.out`, then compares LOSAT's
  run of the row's own query file with `<row_id>.out` as for any `ncbi` row.
- `kind losat-rejection`: LOSAT exits with rc, stdout is empty, stderr contains the text after `expect: ` in covers
  (to a ` | ` or the end), in the library's run (`check --library`) the text after `library-expect: ` where covers has
  one: the library checks an argv before it knows the inputs (`api::validate`), so a refusal that comes from an input's
  word on the command line (NCBI's exclusion of `-phi_pattern` by `-in_msa`) is the library's message of the argv
  alone; the hash columns are empty. A refusal after reports that NCBI writes too has stdout
  cells: NCBI's stdout of the same command, or with `ncbi-near: ROW_ID` that of the `ncbi` row ROW_ID of the same
  outfmt (a near input that NCBI runs to its end, for example the copy row of the batches that LOSAT writes before it
  refuses, or line WQA D10's near input, whose refused batch NCBI's reader stops at), and with
  `ncbi-near: ROW_ID first N bytes` its first N bytes (line RG item T4b-f: where LOSAT refuses after the prolog
  that BLASTP, PSIBLAST, RPSBLAST and TBLASTN write before their first batch, the prolog of a search of the same
  database and outfmt; rows `db.pro.*` and `db.id.pro.*` are such searches with a query without hits). Where NCBI
  crashes on the row's own command,
  its stdout depends on where the crash cuts it: `freeze` stops on such a run, and the row takes a near row (or has
  empty stdout cells where LOSAT writes nothing before it refuses).
- `db`: path relative to LOSAT/ (an empty cell is `-db ""`; names separated by spaces are one `-db` argument;
  `subject:PATH` runs `-subject PATH` in place of `-db`, in any program: the rows of a subject set, such as psiblast's
  losat-rejection row of a subject set whose sequences all have no letters); `extra_args` is split like a shell line;
  `outfmt` is one argument; `stdout_sha256`, `stdout_bytes` and `stderr_sha256` are written by `freeze`.
- A row runs `LOSAT PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT` from `LOSAT/` (from `LOSAT/DIR` with
  `cwd=DIR`, below).
- `query`: a path relative to LOSAT/; `pssm-gz:PATH` or `pssm:PATH` is a PSSM query (`-in_pssm`), and a psiblast query `msa:PATH` is an alignment query (`-in_msa PATH`, line PS item PS6; the library API's `Inputs::msa`, line LIB item 1.6); `stdin:pipe:PATH`,
  `stdin:file:PATH` or `stdin:empty` is `-query -` with PATH's bytes through a pipe, the file PATH itself, or an empty pipe
  as standard input (`db_search_fixtures.query_stdin`; such a row has no PSSM or XInclude files and the library skips it).
  psiblast's `stdin-msa:pipe:PATH`, `stdin-msa:file:PATH` and `stdin-msa:empty` are the same with `-in_msa -`; `none:` gives no query argument
  at all (a row that tests an option that needs `-in_msa`, `-msa_master_idx` or `-ignore_msa_master` without it). A row with `-msa_master_idx` or `-ignore_msa_master` and no `-in_msa`, or with `-in_msa` and
  `-query` or `-in_pssm`, is `library-skip` (`xinclude_fixtures.ALIGNMENT_OPTIONS`); a row of `-in_msa PATH` runs through the library (`Inputs::msa`).

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
  header for every file and every reader). A `cwd=` row with a PSSM, alignment, standard-input or `none:` query, a
  `files` column or a path under `tests/` in extra_args (an argument that would not be rebased) is refused; an `rps-gz:`
  db cell is taken (its names are databases of the inflated directory, not paths). `freeze`
  and `check` refuse to start when `LOSAT/.ncbirc` exists (it would change every row).

## Who runs which rows

- `db_search_fixtures.py check` (and `freeze`): every row, the `cwd=` and `NAME=VALUE` rows included;
  `ci_fast_regressions.py` runs it for every program at 1 and 4 threads, and the `losat-gates` skill uses the same
  checkers.
- `db_search_fixtures.py check --library` (`examples/library_command_line`, the library API): `library-skip` for the
  rows of rpsblast and rpstblastn, those with both `-query` and `-in_pssm`, those with `-in_msa` and `-query` or
  `-in_pssm`, those with the alignment's master options and no `-in_msa`, those with `-export_search_strategy`, those
  with a standard-input query, and those with `cwd=` or a token of `PROCESS_SETTINGS` (FMTCFG, NCBI, GETURL_DEBUG, NCBI_DONT_USE_NCBIRC): the API's contract
  (`LOSAT/src/api/mod.rs:78-86`) names what a call reads of its process, and these formatter settings are not among
  them; the library's host, LOSAT Web, runs `-html` with an empty environment and a run directory of the files that
  the argv names (`docs/web/abi_v2.md` §3), so no Web run has them.
- `check_wasm_threading.py`: a few rows by id, every psiblast row (`psiblast.html.tsv` included) and blastx's rows of outfmt 5 and 13 to 16; no row of
  outfmt 0 to 4 or 18 (the `-html` rows) of the other blast programs (rpsblast and rpstblastn have none there).
- The Rust tests that read this directory read named program files (`blastp.tsv`, `psiblast.tsv`, `blastx.tsv`,
  `tblastn.tsv`, `tblastx.tsv`) and select their own rows, never a family file such as `<program>.html.tsv`;
  `check_manifest_sizes.py` checks the size of every file.

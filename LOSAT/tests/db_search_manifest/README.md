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
  bytes after its last memory page; for exception 19 the hit list size 1073741798), where NCBI's own run of the
  row reaches an approved NCBI defect (PD-LOSAT-NCBI-DEFECTS exceptions 13, 3 and 19); `freeze` copies ROW_ID's
  columns and runs nothing.
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
  write past N bytes fails with EFBIG). Inputs: each program's base search (blastn nslice, blastp prot_v5, tblastn
  tslice, blastx bslice, tblastx xslice; formats 0, 5-12, 15-18 where LOSAT dev and NCBI agree without a device; SAM
  only for blastn) with stdout-full, out-full and `fsize` of half the report (at most 100000) and of its size minus 1
  (where the report is at most 100001 bytes), and for outfmt 0 inside the last query footer, and just inside and 20
  bytes after the epilog's flushed line feeds; `dev_nw.fna` (a warning for an all-N query between two queries, also
  with BATCH_SIZE=1; outfmt 6, 7, 10); `dev_pblk.faa` (1000 copies of p3.faa's first record, 61,000 letters) against
  `-subject dev_psubj.faa` (prot_v5 as FASTA; outfmt 0, 5, 8, 9, 12 with `fsize` 100000); and outfmt 16 with
  `-subject` (`dev_nsubj.fna`, `dev_psubj.faa`; out-full and half). Frozen files stay at or below 100 KB.
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
  (to a ` | ` or the end); the hash columns are empty. A refusal after reports that NCBI writes too has stdout
  cells: NCBI's stdout of the same command, or with `ncbi-near: ROW_ID` that of the `ncbi` row ROW_ID of the same
  outfmt (a near input that NCBI runs to its end, for example the copy row of the batches that LOSAT writes before it
  refuses), and with `ncbi-near: ROW_ID first N bytes` its first N bytes. Where NCBI crashes on the row's own command,
  its stdout depends on where the crash cuts it: `freeze` stops on such a run, and the row takes a near row (or has
  empty stdout cells where LOSAT writes nothing before it refuses).
- `db`: path relative to LOSAT/ (an empty cell is `-db ""`; `subject:PATH` runs `-subject PATH` in place of `-db`, for a
  psiblast losat-rejection row of a subject set); `extra_args` is split like a shell line; `outfmt` is one argument;
  `stdout_sha256`, `stdout_bytes` and `stderr_sha256` are written by `freeze`.
- A row runs `LOSAT PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT` from `LOSAT/`.
- `query`: a path relative to LOSAT/; `pssm-gz:PATH` or `pssm:PATH` is a PSSM query (`-in_pssm`), and a psiblast query `msa:PATH` is an alignment query (`-in_msa PATH`, line PS item PS6; the library API's `Inputs::msa`, line LIB item 1.6); `stdin:pipe:PATH`,
  `stdin:file:PATH` or `stdin:empty` is `-query -` with PATH's bytes through a pipe, the file PATH itself, or an empty pipe
  as standard input (`db_search_fixtures.query_stdin`; such a row has no PSSM or XInclude files and the library skips it).
  psiblast's `stdin-msa:pipe:PATH`, `stdin-msa:file:PATH` and `stdin-msa:empty` are the same with `-in_msa -`; `none:` gives no query argument
  at all (a row that tests an option that needs `-in_msa`, `-msa_master_idx` or `-ignore_msa_master` without it). A row with `-msa_master_idx` or `-ignore_msa_master` and no `-in_msa`, or with `-in_msa` and
  `-query` or `-in_pssm`, is `library-skip` (`xinclude_fixtures.ALIGNMENT_OPTIONS`); a row of `-in_msa PATH` runs through the library (`Inputs::msa`).

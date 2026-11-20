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
- `db`: path relative to LOSAT/ (an empty cell is `-db ""`; `subject:PATH` runs `-subject PATH` in place of `-db`, for a
  psiblast losat-rejection row of a subject set); `extra_args` is split like a shell line; `outfmt` is one argument;
  `stdout_sha256`, `stdout_bytes` and `stderr_sha256` are written by `freeze`.
- A row runs `LOSAT PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT` from `LOSAT/`.
- `query`: a path relative to LOSAT/; `pssm-gz:PATH` or `pssm:PATH` is a PSSM query (`-in_pssm`), and a psiblast query `msa:PATH` is an alignment query (`-in_msa PATH`, line PS item PS6; the library API takes none, so such a row is `library-skip`); `stdin:pipe:PATH`,
  `stdin:file:PATH` or `stdin:empty` is `-query -` with PATH's bytes through a pipe, the file PATH itself, or an empty pipe
  as standard input (`db_search_fixtures.query_stdin`; such a row has no PSSM or XInclude files and the library skips it).

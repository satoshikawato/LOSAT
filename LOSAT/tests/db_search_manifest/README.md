# LOSAT -db fixture rows, one file per program

`<program>.tsv` holds the rows of one program (`blastp blastn tblastn blastx tblastx psiblast rpsblast`; the `program`
column decides the file). `LOSAT/tests/db_search_fixtures.py` is the one loader: `read_manifest` reads the files in that
program order (each in its own row order), `freeze` writes a row back to its own program's file. Every file starts with
the same header line, after its `# Rows ...` comment lines (a comment that describes rows of several programs is
repeated in each of their files). `LOSAT/tests/check_manifest_sizes.py` fails CI when a file passes 512 KiB: split
that program's rows by family into more files and teach `read_manifest` and `write_manifest` the new names.

Columns: row_id, program, kind, db, query, extra_args, outfmt, rc, stdout_sha256, stdout_bytes, stderr_sha256, covers,
and for the rows of outfmt 13 and 14 and psiblast's PSSM files two more, files and side_sha256 (a row without them has
12 cells; `xinclude_fixtures.py` says what they mean).

- `kind ncbi`: stdout, stderr and exit status (rc) equal the frozen bytes of NCBI BLAST+ 2.17.0
  `PROGRAM -db ... -num_threads 1` (`LOSAT/tests/fixtures/db_search/<row_id>.out`, `.err` when stderr is not empty).
- `kind losat-rejection`: LOSAT exits with rc, stdout is empty, stderr contains the text after `expect: ` in covers
  (to a ` | ` or the end); the hash columns are empty.
- `db`: path relative to LOSAT/ (an empty cell is `-db ""`; `subject:PATH` runs `-subject PATH` in place of `-db`, for a
  psiblast losat-rejection row of a subject set); `extra_args` is split like a shell line; `outfmt` is one argument;
  `stdout_sha256`, `stdout_bytes` and `stderr_sha256` are written by `freeze`.
- A row runs `LOSAT PROGRAM -db DB -query QUERY EXTRA_ARGS -outfmt OUTFMT` from `LOSAT/`.

- `blastdbcheck -cdd_delta`: validates DELTA-BLAST profile database headers (`.aux`, `.loo`,
  `.rps`, `.freq`, `.wcounts`, `.obsr`) including magic numbers, profile counts, offset tables,
  carried state across volumes, and exact NCBI exception reporting; empty (0-byte) profile files
  where NCBI crashes remain refused.
- Shared database reader (`seqdb`): reads zero-filled bytes past EOF within the last 4 KiB memory
  page on truncated index, sequence, and header files to match NCBI BLAST+ 2.17.0 across
  `blastdbcheck`, `blastdbcmd`, and the search programs; handles non-ascending index offsets
  (including multithreaded DB-split scan stop at the first negative-length OID — whether included or
  excluded by a list — across `blastp`, `blastx`, `psiblast`, and `blastn`, and in-mapping
  out-of-order nucleotide ambiguity offsets in `blastn`, `blastdbcmd`, and `blastdbcheck`, while
  `tblastn` and `tblastx` continue to refuse out-of-order ambiguity offsets, and offsets reading
  before the mapping or multithreaded `blastn` traceback on `amb > next` remain refused),
  `SEQDB_FILE_ASSERT` exception messages in `blastdbcheck`, and `Error: (803.7)` diagnostics for
  `#`-replaced definition-line bytes in `blastdbcheck` and failed `blastdbcmd` header decodes.
- `blastdbcheck`: widens alias count checks (`i128`) and `-ends` OID bounds (`i64`) so
  signed integer overflow does not cause false alias errors or wrapped OIDs
  (`PD-LOSAT-NCBI-DEFECTS` version 1.120, approved exception 35).

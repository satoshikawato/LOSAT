- `blastdbcheck -cdd_delta`: validates DELTA-BLAST profile database headers (`.aux`, `.loo`,
  `.rps`, `.freq`, `.wcounts`, `.obsr`) including magic numbers, profile counts, offset tables,
  carried state across volumes, and exact NCBI exception reporting; empty (0-byte) profile files
  where NCBI crashes remain refused.
- Shared database reader (`seqdb`): reads zero-filled bytes past EOF within the last 4 KiB memory
  page on truncated index, sequence, and header files to match NCBI BLAST+ 2.17.0 across
  `blastdbcheck`, `blastdbcmd`, and the search programs; handles non-ascending index offsets
  (including multithreaded DB-split scan stop at the first negative-length OID and in-mapping
  out-of-order nucleotide ambiguity offsets, refusing only offsets before the mapping and
  multithreaded `blastn` traceback where NCBI crashes or aborts),
  `SEQDB_FILE_ASSERT` exception messages in `blastdbcheck`, and `Error: (803.7)` diagnostics for
  `#`-replaced definition-line bytes in `blastdbcheck` and failed `blastdbcmd` header decodes.
- `blastdbcheck`: widens alias count checks (`i128`) and `-ends` OID bounds (`i64`) so
  signed integer overflow does not cause false alias errors or wrapped OIDs
  (`PD-LOSAT-NCBI-DEFECTS` version 1.117, approved exception 33).

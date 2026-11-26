- Non-UTF-8 file paths, database names, mask keys, option values, and `WINDOW_MASKER_PATH`
  entries now work across all search programs (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`,
  `psiblast`, `rpsblast`, `rpstblastn`), `makeblastdb`, and `rpsbproc`, preserving the raw bytes
  of paths, keys, delimiters (`-outfmt "6 delim=..."`), `-export_search_strategy`, `-use_index`,
  reports, and diagnostics to match NCBI BLAST+ 2.17.0 and RpsbProc 0.5.1.
- `makeprofiledb` now stops where NCBI BLAST+ 2.17.0 stops when standard output or `-logfile` is a
  closed pipe, with exit status 0 where NCBI is ended by SIGPIPE (exit status 141): fresh builds
  and `-dbtype delta` exclusion lines within NCBI's stream buffer complete with all database files,
  exclusion lines that overflow the stream buffer (`cout` or `-logfile`) stop before creating
  volumes, diagnostics written to `cerr` flush buffered `cout` bytes first (`cerr` tied to `cout`,
  preserving NCBI's `2>&1` line order and closed-pipe stop points), overwriting an existing
  database deletes the old files and stops at the log line, and `CWriteDB` stdout lines stop inside
  `AddSequence` without flushing buffered `-logfile` bytes or running destructor clean-up (recorded
  as exception 10 in `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` v1.35).

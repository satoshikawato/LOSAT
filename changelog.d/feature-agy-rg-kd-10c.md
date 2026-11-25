- Non-UTF-8 file paths, database names, mask keys, option values, and `WINDOW_MASKER_PATH`
  entries now work across all search programs (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`,
  `psiblast`, `rpsblast`, `rpstblastn`), `makeblastdb`, and `rpsbproc`, preserving the raw bytes
  of paths, keys, delimiters (`-outfmt "6 delim=..."`), `-export_search_strategy`, `-use_index`,
  reports, and diagnostics to match NCBI BLAST+ 2.17.0 and RpsbProc 0.5.1.
- `makeprofiledb` now matches NCBI BLAST+ 2.17.0 when standard output or `-logfile` is a closed
  pipe: fresh builds complete with exit status 0, overwriting an existing database deletes the old
  files and stops at the log line with exit status 0 and no message, and `CWriteDB` stdout lines
  stop inside `AddSequence` with exit status 0 and no message (recorded as exception 10 in
  `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` v1.33).

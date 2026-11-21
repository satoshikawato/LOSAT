- Changed: outputs name LOSAT and its version where NCBI BLAST+ names its own. The first line of
  `-outfmt 0`, the `# BLASTP (LOSAT 0.2.0)` line of `-outfmt 7`, `BlastOutput_version` of `-outfmt 5`
  and the `version` of `-outfmt 13` to `16` read `BLASTP (LOSAT 0.2.0)` in place of `BLASTP 2.17.0+`
  (the program's name stays first); an archive's `ident`, SAM's `@PG VN:` and `-version` read
  `LOSAT 0.2.0` (`blastn: LOSAT 0.2.0`), and the program's `<version>` of `-xmlhelp` reads `0.2.0`
  (NCBI's `2.17.0`). Error messages start with `LOSAT engine error:`,
  `LOSAT query/options error:`, `LOSAT database error:` and the other `LOSAT ... error:` prefixes in
  place of NCBI's `BLAST ... error:`. Every other byte of the search results, the message bodies and
  the exit statuses are unchanged. To compare a LOSAT output with NCBI's, map NCBI's file with
  `python3 -I LOSAT/tests/losat_wording.py map-output <file>`.

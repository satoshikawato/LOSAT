- Documentation: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 1.22 records `blastn`'s progress bar as exception 11
  (pending the Owner's review): on a terminal, `losat blastn` with `-num_threads` 2 or more and two subjects
  or more draws it on the standard error, where NCBI BLAST+ 2.17.0 writes nothing; a standard error sent to
  a file or a pipe gets none of it, and the library API, pylosat and the Web draw none. The manual's
  differences page and AGENTS.md name it. No output changes.

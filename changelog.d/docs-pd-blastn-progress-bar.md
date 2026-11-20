- `blastn`: the progress bar that the command line draws on a terminal (`-num_threads` 2 or more, two subjects or
  more) is now cleared when the search of each query batch ends, so nothing of it stays on the screen; before, the
  finished bar stayed, one line per query batch. A standard error sent to a file or a pipe still gets none of it,
  and the standard output, the `-out` file and the exit status do not change.
- Documentation: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` 1.22 records the bar as exception 11 (the Owner approved its
  text on 2026-11-20 by pressing "Ready for review" on #1355): NCBI BLAST+ 2.17.0 writes nothing to the standard
  error while it searches; the library API, pylosat and the Web draw no bar. The manual's differences page and
  AGENTS.md name it.

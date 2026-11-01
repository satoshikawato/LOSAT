- The variables `LOSAT_X_BXDBSTREAM`, `LOSAT_X_BXEVALEMPTY` and `LOSAT_X_BXLUTPOOL` are gone: the blastx paths
  they chose were already the default and now run unconditionally, so `=0` no longer reaches the previous code.
  What went is that previous code: the copy of the whole `-db` database made before the first query batch
  (a pass now reads blocks of 256 KiB, which `LOSAT_X_BXDBSTREAM=KiB` set), the full statistics evaluation of a
  subject without HSPs, and the lookup table build without shared neighbour lists and kept vectors for a search of
  two threads or more (a search on one thread still builds its tables without the pool, `LOSAT_X_BXPOOLN1`). The
  shadow modes `LOSAT_X_BXDBSTREAMSHADOW` and `LOSAT_X_BXEVALEMPTYSHADOW` and their `[X_SHADOW]` lines went with
  them. Outputs do not change.
- `blastx -db` with a `-seqidlist`, `-taxids` or other list reads the scanned shortest length in its hit-saving cutoff as NCBI BLAST+ does (the default path answered 10, and could keep rows NCBI drops).

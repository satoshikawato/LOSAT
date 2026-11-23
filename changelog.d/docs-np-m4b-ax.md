- Documentation: `PD-LOSAT-NCBI-DEFECTS` 1.83 extends approved exception 3 to NCBI BLAST+ 2.17.0's traceback read of one
  letter past a hard-masked copy of a subject in `blastp` and `blastx` with `-db_hard_mask`: where that byte happens to
  equal the fence byte twice, NCBI crashes (SIGSEGV), sometimes after part of its rows, and LOSAT writes the rows of
  NCBI's runs that do not meet it (the Owner's decision of 2026-11-23). LOSAT's code and output do not change.

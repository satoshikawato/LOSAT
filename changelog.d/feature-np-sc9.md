- `-word_size` 2, 4, 6 and 7 in BLASTX (NCBI BLAST+ 2.17.0 takes 2 to 7): the lookup table of 2 and 4 letters, and above 4 the
  compressed lookup table with NCBI's thresholds 21 (6) and 20.25 (7); with `-task blastx-fast`, every matrix, `-ungapped`,
  `-window_size 0` and the masking options. IDENTITY with the compressed lookup table stays refused (NCBI crashes). LOSAT Web's
  BLASTX Word size list is NCBI's (2, 3, 5, 6).
- `blastx -word_size 4` searches one query unit at a time, or `ceil(threads / subjects)` with few subjects, each unit's lookup table built on the threads (`LOSAT_X_BXW4CAP[=K]`, on unless `=0`; the Owner approved it as a deadline fallback until the table is made smaller). 50 cDNAs x yeast proteins at 8 threads hold 189 MB instead of 820 MB (NCBI 241 MB) in 0.95 of the time; 800 cDNAs 463 MB instead of 1,222 MB in 1.16 of the time. Reports do not change.

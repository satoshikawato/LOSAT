- `-word_size` 2, 4, 6 and 7 in BLASTX (NCBI BLAST+ 2.17.0 takes 2 to 7): the lookup table of 2 and 4 letters, and above 4 the
  compressed lookup table with NCBI's thresholds 21 (6) and 20.25 (7); with `-task blastx-fast`, every matrix, `-ungapped`,
  `-window_size 0` and the masking options. IDENTITY with the compressed lookup table stays refused (NCBI crashes). LOSAT Web's
  BLASTX Word size list is NCBI's (2, 3, 5, 6).

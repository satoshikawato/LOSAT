- BLASTX and TBLASTN use less memory for the compressed lookup table of `-word_size 5` (also
  `-task blastx-fast` and `-task tblastn-fast`) at low `-threshold` values, at most what NCBI BLAST+
  2.17.0 uses, and build it faster. One 8,714-base BLASTX query, `-word_size 5 -threshold 2`: 3.82 GB of
  maximum RSS instead of 5.07 GB (NCBI 3.85 GB) in about half the time; the largest table LOSAT
  builds (8,693 bases, `-threshold 1`): 4.98 GB instead of 6.59 GB (NCBI 5.00 GB). TBLASTN
  `-task tblastn-fast` on an 18,960-letter query batch at `-threshold 5`: 1.35 GB instead of 2.62 GB
  (NCBI 2.04 GB), 0.63 of the time; a BLASTX query batch of more than 65,536 letters at `-threshold 8`:
  1.73 GB instead of 3.37 GB (NCBI 2.59 GB). Other matrices than BLOSUM62 (`-matrix BLOSUM45`,
  `PAM70`): 0.75 to 0.81 of the memory; BLASTP's `-task blastp-fast` at low thresholds: 0.54 to 0.65.
  The default thresholds' tables, and every output, are unchanged.

- `tblastn -task tblastn-fast` and `tblastn -word_size 5`, `6` and `7` run NCBI's compressed-alphabet lookup table
  (BLOSUM62) with `-subject` and `-db`, and TBLASTN with BLOSUM62 11/1 takes any `-threshold` and any `-window_size`
  above 0 (word size 3, 5, 6 or 7): AMRFinderPlus's two tblastn calls (`-task tblastn-fast -word_size 5
  -window_size 15 -threshold 100` and `-word_size 3 -window_size 15 -threshold 100`) give NCBI's output in the
  standard tabular columns. The compressed table with another matrix is refused.

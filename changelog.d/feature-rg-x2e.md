- `tblastx` with `-window_size 0` (or `-lcase_masking`) and `-word_size 4`, `-word_size 2` or a matrix
  other than BLOSUM62, and `tblastn` with `-word_size 4`, use less memory: their word lookup table is
  built the way `blastx` builds it. At `-word_size 4 -threshold 1`, `tblastx -window_size 0` used
  1.31 times NCBI BLAST+'s memory and `tblastn` 1.12 times. Output is unchanged.

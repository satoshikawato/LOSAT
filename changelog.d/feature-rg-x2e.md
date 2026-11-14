- `tblastx` with `-word_size 4`, `-word_size 2` or a matrix other than BLOSUM62, where it builds one word table
  for a whole query batch (`-window_size 0`, `-lcase_masking`, or a short query on one thread), and `tblastn`
  with `-word_size 4` use less memory and time. At `-word_size 4 -threshold 1`, `tblastx -window_size 0` now uses
  0.57 times NCBI BLAST+ 2.17.0's memory (was 1.31 times) and runs about 6 times faster than before; `tblastn` uses
  0.54 times NCBI's memory (was 1.12 times) and runs about 2.7 times faster. Output is unchanged.

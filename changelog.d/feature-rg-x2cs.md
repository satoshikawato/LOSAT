- `tblastx -window_size 0` (the one-hit extension) and `tblastx` with an `-xdrop_ungap` above about
  58 bits (BLOSUM62) take about 12 % and 4 % fewer instructions (one thread), and about 31 % fewer
  with `-xdrop_ungap 20000` on three test inputs: their ungapped extensions read each score in two
  loads, as NCBI BLAST+ does, instead of testing both letters first. Reports do not change.

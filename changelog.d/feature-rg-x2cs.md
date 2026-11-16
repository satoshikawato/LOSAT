- `tblastx -window_size 0` (the one-hit extension) and `tblastx` with an `-xdrop_ungap` above about
  58 bits (BLOSUM62) read each score of their ungapped extensions in two loads, as NCBI BLAST+ does,
  instead of testing both letters first: they take about 12 % and 4 % fewer instructions (one
  thread), and with an `-xdrop_ungap` of about 15,000 bits or more the one-hit extension takes 1.14
  times NCBI BLAST+'s instructions where it took 1.68 times. Reports do not change.

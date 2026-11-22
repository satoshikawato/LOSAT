- `makeblastdb -input_type blastdb -randomize` copies the sequences in a random order, as NCBI
  BLAST+ 2.17.0 does: the order NCBI's Linux builds give for the same second of the clock, the
  OIDs of alias files included. LOSAT refused the option.

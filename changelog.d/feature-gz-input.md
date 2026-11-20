- Input files named `.gz` are read through gzip: `-query` and `-subject` of every search program
  (NCBI BLAST+ reads such a file as text and stops at its first line) and `makeblastdb -in`, which
  NCBI's makeblastdb also inflates; concatenated members (bgzip) are read whole, and the name alone
  decides, as in makeblastdb (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 8). `makeblastdb -in`
  names ending in `.bz2` or `.zst` stay refused.

- BLASTX `-word_size 4` searches one query batch at a time also when there are few subjects: at 16 or more threads with a handful of
  `-subject` proteins it no longer holds two or four lookup tables at once, so its memory stays below NCBI's (800 cDNAs x 8
  proteins at 16 / 32 threads: 309 / 316 MB instead of 520 / 795 MB, NCBI 386 MB) in 1.4-1.5 times the time (4.5 / 3.4 s,
  NCBI 141 s). Reports do not change.

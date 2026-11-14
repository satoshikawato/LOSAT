- BLASTX `-word_size 6` and `7` use less memory with few queries at many threads: each query batch's lookup table is
  built a run of its cells at a time, so the build holds one run's list of word positions instead of the whole
  table's. 50 cDNAs x yeast proteins at 8 threads: word size 7 about 410 MB instead of 800 MB, word size 6 about
  305 MB instead of 455 MB. Reports do not change.

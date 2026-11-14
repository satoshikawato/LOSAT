- BLASTX `-word_size 5` and `-task blastx-fast` at 4 threads or more build each query batch's lookup table on the
  threads the batch is given, as word sizes 6 and 7 do, and the tables of word sizes 5 to 7 are built into memory kept
  from the tables before. Reports do not change.

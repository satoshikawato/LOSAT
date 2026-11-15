- BLASTX `-word_size 6` and `7` use less memory and run faster at 8 to 32 threads: each query batch's lookup table is
  built a run of its cells at a time, so the build holds one run's list of word positions instead of the whole
  table's. 50 cDNAs x yeast proteins: word size 7 about 380-430 MB instead of 680-770 MB, word size 6 about 290 MB
  instead of 430-440 MB (NCBI BLAST+ 435 and 367 MB), in 0.79-0.93 of the time; 800 cDNAs in 0.87-0.96 of the time
  with 20-55 % less memory. Reports do not change.

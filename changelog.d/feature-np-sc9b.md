- BLASTX and TBLASTX `-word_size 4` are faster and use less memory: their lookup tables of 4-letter words are built by
  stripe of cells on all the threads. BLASTX 800 cDNAs x yeast proteins at 8 / 16 / 32 threads takes 16.5 / 9.9 / 7.0 s
  instead of 21.4 / 18.5 / 16.4 s (NCBI 48.2 / 27.8 / 33.8 s) and 325 / 421 / 495 MB instead of 494 / 451 / 565 MB
  (NCBI 1,563 / 2,338 / 4,404 MB); 50 cDNAs 1.30 / 1.87 / 2.18 times as fast, 148 / 154 / 174 MB instead of 178 / 222 /
  232 MB (NCBI 237 / 234 / 230 MB). TBLASTX of two genomes (AP027131 x AP027133, genetic code 4) is 1.07 times as fast
  at 1 thread with 7 % less memory, 50 cDNAs x a genome 1.36 / 1.48 times as fast at 1 / 8 threads with 13 / 22 % less
  memory. Reports do not change.

- BLASTX `-word_size 5` and `-task blastx-fast` are 1.3 to 3 times as fast at 8 to 32 threads, and word sizes 6 and 7
  1.1 to 1.4 times as fast with less memory: each query batch's lookup table is built on the threads the batch is
  given. 800 cDNAs x yeast proteins with blastx-fast at 8 / 16 / 32 threads take 2.05 / 1.19 / 0.90 s instead of
  2.76 / 2.48 / 2.46 s (NCBI 4.08 / 2.47 / 1.85 s); word size 7 at 32 threads 11.5 s and 965 MB instead of 15.9 s
  and 1,289 MB (NCBI 31.4 s and 9,640 MB). BLOSUM45 at word size 5 uses 25-50 % less
  memory. Reports do not change.

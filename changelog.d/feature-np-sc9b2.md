- BLASTX `-word_size 6` and `7` (also with `-task blastx-fast` and other matrices) are 2.4 to 5.7 times as fast at 8 to
  32 threads and use less memory: the lookup table of each query batch is built on the threads the batch is given
  instead of one. 800 cDNAs x yeast proteins, word size 6, at 8 / 16 / 32 threads takes 18.4 / 12.7 / 9.6 s instead
  of 45.0 / 43.7 / 46.6 s (NCBI 47.4 / 28.2 / 18.7 s) and 911 / 995 / 1,052 MB instead of 1,174 / 1,215 / 1,420 MB
  (NCBI 2,495 / 4,701 / 9,177 MB); word size 7 29.4 / 20.6 / 14.7 s instead of 81.0 / 81.2 / 83.7 s (NCBI 85.5 /
  49.4 / 31.6 s); BLOSUM45 word size 6 at 8 threads 76 s instead of 224 s (NCBI 224 s). Reports do not change.

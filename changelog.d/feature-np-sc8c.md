- BLASTP with the compressed lookup table (blastp-fast, `-word_size` 5 to 7) and a lowered `-threshold` keeps
  its lookup table near the size of NCBI BLAST+'s: a block of queries is also closed on the words of its table,
  at most four per backbone cell plus the words of the largest NCBI query batch. `-word_size 6 -threshold 11` on
  89 proteins took 14.7 GB at one thread and now takes about 3.3 GB (NCBI BLAST+ 2.6 GB); outputs are unchanged.

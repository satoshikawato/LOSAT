- BLASTP with the compressed lookup table (blastp-fast, `-word_size` 5 to 7) keeps the table of each block of
  queries within about the size of NCBI BLAST+'s largest query-batch table: a block of LOSAT's block search (more
  than 50,000 query letters) also ends before a batch whose words would pass four per backbone cell plus the
  words of the largest batch. With a lowered `-threshold` the largest RSS falls: `-word_size 6 -threshold 11` on
  89 proteins (64,069 letters) at one thread from 14,662 MiB to 3,231 MiB (NCBI BLAST+ 2,579 MiB). At the default
  thresholds the blocks of blastp-fast and of word size 6 are unchanged; those of word sizes 5 and 7 can end one
  batch earlier. Outputs are unchanged. Query sets of at most 50,000 letters still build one table.

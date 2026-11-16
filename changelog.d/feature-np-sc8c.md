- BLASTP with the compressed lookup table (blastp-fast, `-word_size` 5 to 7) keeps the table of each block of
  queries near the size of NCBI BLAST+'s largest query-batch table: a block of LOSAT's block search (more than
  50,000 query letters) also ends before a batch whose words would pass four per backbone cell plus the words of
  the largest batch. With a lowered `-threshold` the largest RSS falls by about three quarters. On 89 proteins
  (64,069 letters) against WSSV proteins: `-word_size 6 -threshold 11` from 14,663 MiB to 3,231 MiB at one thread
  (NCBI BLAST+ 2,582 MiB) and from 15,066 MiB to 3,585 MiB at 8 threads; `-word_size 6 -threshold 15` from
  3,415 MiB to 810 MiB; `-task blastp-fast -threshold 13` from 454 MiB to 110 MiB. At the default thresholds
  `-word_size 7` falls by about 110 MiB (532 MiB to 425 MiB) and word sizes 5 and 6 and blastp-fast fall by 0 to
  8 MiB. Some searches with a lowered `-threshold` take longer: `-word_size 6 -threshold 11` at 8 threads 3.5 %
  (8.38 s to 8.66 s), `-word_size 6 -threshold 15` and `-task blastp-fast -threshold 13` at one thread about
  1 %. Others are faster: `-word_size 6 -threshold 11` at one thread 74.4 s to 68.1 s, `-task blastp-fast
  -threshold 11` at 8 threads 10 %. Outputs are unchanged. Query sets of at most 50,000 letters still build one
  table.

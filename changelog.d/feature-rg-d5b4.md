- BLASTN: `-use_index true` (a search with a MegaBLAST index made by `makembindex`) is faster at
  every `-num_threads` and uses less memory, with the same reports: each index volume is searched
  on all the threads, and a query whose words hit many places of one sequence at once (repeat
  islands) no longer slows the search down. On the workstation (16 cores of one socket),
  EDL933 against an index of 2,000-letter pieces of two E. coli genomes took 0.48 s instead of
  0.64 s at one thread and 0.16 s instead of 0.54 s at 16 threads (NCBI BLAST+ 2.17.0: 1.12 s and
  0.98 s), with 117 MB instead of 140 MB; a query with repeat islands took 7.8 s instead of 17.4 s
  at one thread and 2.9 s instead of 12.5 s at 16 threads (NCBI: 13.0 s and 3.9 s).

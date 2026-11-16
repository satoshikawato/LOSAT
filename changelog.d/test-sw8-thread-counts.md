- Tests: 136 loops in 48 files of LOSAT's tests that compared the output at 1 and 4 threads now
  also compare it at 8 and 16 (blastn, blastp, blastx, tblastn and tblastx; pairwise, tabular,
  XML, JSON, archive and side-file reports; database and `-subject` searches), where the test took
  3 s or less. Slower tests, and one whose fixture meets the TBLASTX cost rule only up to four
  threads, keep 1 and 4. Eight tests named for one and four threads are renamed
  `..._at_1_4_8_and_16_threads`. Two TBLASTX linking tests now compare every thread count, not only
  1 and 4, and `tblastn_newton_drop` compares the Newton results dropped at one thread only (on a
  pool, results left after a thread's last prefetch are not counted, and which matches a thread
  walks last depends on the schedule). Also: NCBI's `ParseSequenceRange_InvalidDelimiter` case
  (`3,4`) is ported and the six `ParseSequenceRange` cases are cited above their test, three
  `-task` citations now point at `blast_args.cpp:105`, and the `writedb` module documentation says
  what makeblastdb writes today.

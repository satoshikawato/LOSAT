- Tests: 130 loops of LOSAT's tests that compared the output at 1 and 4 threads now also compare it
  at 8 and 16 (blastn, blastp, blastx, tblastn and tblastx; pairwise, tabular, XML, JSON, archive
  and side-file reports; database and `-subject` searches), where the test takes 3 s or less.
  Slower tests and two that compare counters a thread-count rule sets keep 1 and 4. Eight tests
  named for one and four threads are renamed `..._at_1_4_8_and_16_threads`. Also: NCBI's
  `ParseSequenceRange_InvalidDelimiter` case (`3,4`) is ported, three `-task` citations now point
  at `blast_args.cpp:105`, and the `writedb` module documentation says what makeblastdb writes today.

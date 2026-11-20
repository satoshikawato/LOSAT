- Megablast's greedy alignment runs each distance in two passes per block of diagonals instead of three on x86_64
  CPUs with AVX2 (`LOSAT_X_GREEDYROWFUSE`, on by default; `=0` runs the three passes): the eight-byte match test
  runs as scalar code inside the pass over the diagonals. On the workstation megablast of Sakai against MG1655 takes
  0.90-0.91 of the cycles at one thread (0.623 -> 0.567 s) and 0.92 at 8 and 16 threads; 16S rRNAs against a 16S
  database 0.99. Outputs do not change; other CPUs and targets are unchanged.

- `blastp -use_sw_tback` and `rpsblast -use_sw_tback` take about 30 % less time on x86-64 processors
  with AVX2 (less on others): the Smith-Waterman alignment that finds each match's score and start does
  less work per cell and, with AVX2, takes the processor's packed maximum. At `-num_threads 1` on a Xeon
  Gold 6338, `blastp -use_sw_tback` took 4.17 s and 20.92 s on two queries where it took 5.87 s and
  29.34 s before (NCBI BLAST+ 2.17.0: 5.68 s and 29.68 s), and `rpsblast -use_sw_tback` 4.88 s and
  23.19 s where it took 7.09 s and 34.97 s (NCBI: 23.16 s and 131.45 s). `rpstblastn -use_sw_tback` runs
  the same alignment. Outputs and peak memory do not change.

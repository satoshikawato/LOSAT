- `blastp -use_sw_tback` and `rpsblast -use_sw_tback` take about 30 % less time on x86-64 processors
  with AVX2 (measured on a Xeon Gold 6338). At `-num_threads 1`, `blastp -use_sw_tback` took 4.17 s and
  20.92 s on two queries where it took 5.87 s and 29.34 s before (NCBI BLAST+ 2.17.0: 5.68 s and 29.68 s),
  and `rpsblast -use_sw_tback` took 4.88 s and 23.19 s where it took 7.09 s and 34.97 s (NCBI: 23.16 s and
  131.45 s). Outputs and peak memory do not change.

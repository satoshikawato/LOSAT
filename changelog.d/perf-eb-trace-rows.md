- The X-drop traceback stores its rows in less memory: each row's record takes two bytes instead of eight (the row's cells and
  the step of its first column; rows of 255 cells or more take nine), the rows go into 64 KiB chunks instead of 1 MiB, and the
  walk over the rows runs in a function of its own. Reports do not change. Sakai x EDL933 (5.5 Mb genomes) `-subject` at
  16 threads: peak memory 0.86 of the previous version's for blastn and 0.98 for dc-megablast; 0.97 at 1 thread for both.
  Time is unchanged at 16 threads; at 1 thread the dc-megablast search takes 0.6 % more CPU time, and tblastn `-db` 0.3 %.

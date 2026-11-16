- The flat compressed lookup table is built a little faster where many query offsets share a cell (BLASTP
  `-task blastp-fast` and `-word_size` 6 or 7, TBLASTN `-word_size` 6 or 7, BLASTX-fast): each cell of more than
  five offsets is put into NCBI's scan order in place instead of through a copy. On 89 proteins (64,069 letters)
  against WSSV proteins at one thread: BLASTP `-task blastp-fast -threshold 13` 1.67 s to 1.59 s, `-word_size 6
  -threshold 15` 17.04 s to 16.47 s, `-word_size 7` 2.81 s to 2.75 s; at 8 threads `-word_size 6 -threshold 11`
  9.41 s to 9.13 s. Tables and outputs are unchanged.

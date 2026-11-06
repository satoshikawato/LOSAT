- The flat compressed lookup table of `-word_size` 6 and 7 (BLASTP, TBLASTN) is built faster at one thread: the
  builder that counts and places each cell's query offsets reads the cells a few appends ahead (a cache hint) and
  no longer moves the cells' starts after the fill. BLASTP at one thread on 89 proteins (64,069 letters) against
  WSSV proteins: `-word_size 7` 3.58 s to 2.81 s, `-word_size 6` 2.09 s to 1.71 s, `-word_size 6 -threshold 15`
  22.37 s to 16.96 s. Tables, outputs and memory are unchanged.

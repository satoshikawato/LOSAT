- The flat compressed lookup table of `-word_size` 6 and 7 (BLASTP, TBLASTN) is built faster: the builder that
  counts and places each cell's query offsets reads the cells a few appends ahead on x86_64 (a cache hint) and no
  longer moves the cells' starts after the fill. On 89 proteins (64,069 letters): BLASTP at one thread against WSSV
  proteins `-word_size 7` 3.58 s to 2.83 s, `-word_size 6` 2.10 s to 1.72 s, `-word_size 6 -threshold 15`
  22.44 s to 17.06 s; TBLASTN against 5 nucleotide subjects (31 kb) `-word_size 6` 1.57 s to 1.31 s (0.84 s to
  0.74 s at 8 threads), `-word_size 7` 2.52 s to 1.97 s. Tables, outputs and memory are unchanged.

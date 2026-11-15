- Speed: BLASTN no longer spends time on a subject chunk without an initial hit after its word
  finder (NCBI goes on to the next chunk there), and drops other fixed work of every chunk. A
  megablast database search of many short subjects runs about 10 % faster at one thread (about
  20 % with `-db_hard_mask`), and a 16S amplicon against SILVA about 4 % faster (6 % at 16
  threads). Output is unchanged.

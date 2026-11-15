- Speed: BLASTN no longer copies the query batch's search values for every database sequence and
  every chunk of a sequence (NCBI also hands them over by reference), and a database read into a
  reused buffer no longer clears it first. A megablast search of 40,000 short sequences runs about
  4 % faster at one thread (7 % with `-db_hard_mask`), a 16S amplicon against SILVA about 1 %
  faster (2 % at 16 threads), and 20 cDNAs against the zebrafish genome about 1.6 % faster.
  Output is unchanged.

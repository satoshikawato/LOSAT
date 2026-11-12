- blastx and tblastn: `-subject` sequences whose average length is below one letter, in a search without
  NCBI's Gumbel block (Spouge's finite-size correction, which `-searchsp`, `OLD_FSC` and blastx
  `-ungapped` turn off), end with NCBI BLAST+'s `BLAST engine error: The average subject length is too
  short` (exit status 3) after the start of the report, as NCBI does. tblastn refused these runs before
  the report (exit status 1), blastx wrote its own text in that place, and blastx `-ungapped` searched
  them (exit status 0).

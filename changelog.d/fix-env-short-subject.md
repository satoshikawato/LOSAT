- blastx and tblastn: `-subject` sequences whose average length is below one letter, in a search without
  the finite-size correction (`-searchsp`, `OLD_FSC`, or blastx `-ungapped`), end with NCBI BLAST+'s
  `BLAST engine error: The average subject length is too short` (exit status 3) after the start of the
  report, as NCBI does. tblastn refused these runs before the report (exit status 1), and blastx
  `-ungapped` searched them (exit status 0).

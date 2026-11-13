- `tblastn` and `tblastx` accept `-window_size 0`, NCBI's one-hit word finder (`s_BlastAaWordFinder_OneHit`),
  with `-subject` and `-db`, every matrix and word size, `-in_pssm` (`tblastn`) and on the Web, giving NCBI
  BLAST+ 2.17.0's output. Two rare cases of a `-threshold` below the default stop with an error instead: a
  `tblastx` one-hit alignment whose length is not positive, and a `tblastn` one whose gapped start reads past
  the query batch or the subject's last frame (NCBI BLAST+ reads past its buffer there).

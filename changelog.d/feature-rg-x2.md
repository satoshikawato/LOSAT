- `tblastn` and `tblastx` accept `-window_size 0`, NCBI's one-hit word finder (`s_BlastAaWordFinder_OneHit`),
  with `-subject` and `-db`, every matrix and word size, `-in_pssm` (`tblastn`) and on the Web, giving NCBI
  BLAST+ 2.17.0's output. Three rare cases of a low `-threshold` stop with an error instead: a `tblastx`
  one-hit alignment whose length is not positive (LOSAT does not reproduce NCBI BLAST+'s result there yet),
  and a `tblastn` one of negative length whose gapped start reads past NCBI BLAST+'s buffer of the query
  batch or of the subject's frames, or with `-in_pssm` past the PSSM (NCBI BLAST+ reads outside its memory
  there).

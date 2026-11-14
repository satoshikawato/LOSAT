- `tblastn` and `tblastx` accept `-window_size 0`, NCBI's one-hit word finder (`s_BlastAaWordFinder_OneHit`),
  with `-subject` and `-db`, every matrix and word size, `-in_pssm` (`tblastn`) and on the Web, giving NCBI
  BLAST+ 2.17.0's output. Five rare cases of a low `-threshold` stop with an error instead: a `tblastx`
  one-hit alignment whose length is not positive, and with `-in_pssm` a `tblastn` one of negative length
  whose gapped start lies after the last letter of its frame (LOSAT does not reproduce NCBI BLAST+'s results
  there yet); a `tblastn` one of negative length whose gapped start reads past NCBI BLAST+'s buffer of the
  query batch or of the subject's frames, or with `-in_pssm` past the PSSM (NCBI BLAST+ reads outside its
  memory there); and two whose gapped starts lie at the end of their 5,000,000-letter part of a frame
  or past it (NCBI BLAST+ never finishes there).

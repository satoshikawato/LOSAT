- `tblastx` with a very large `-xdrop_ungap` (about 15,000 bits or more with BLOSUM62) now lets an
  ungapped extension that reaches the end of a query frame go on past it into the next frame of the
  query batch, as NCBI BLAST+ 2.17.0 does, so where it does not score higher there the alignments
  that NCBI BLAST+ leaves out (later hits that the extension covers) are left out too, at every
  `-num_threads`. With `-window_size 0` the one-hit extension of a word at the end of a frame whose
  score never rises above 0 (such as `XXX`) now starts its right part on the next frame's first
  letter, as NCBI BLAST+'s code does. Where such an extension scores higher in the next frame,
  NCBI BLAST+ reports an alignment of two frames with query coordinates outside the query's frame
  (past the query's length) or none, and `tblastx` now stops with an error instead.

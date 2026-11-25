- `tblastn` and `tblastx` with `-db_soft_mask ID` on a database sequence of more than about
  15,000,000 bases whose stored masks are out of order give NCBI BLAST+ 2.17.0's `-num_threads 1`
  rows at every `-num_threads` also where NCBI restarts its table of diagonals inside that sequence
  (more than about 270 Mbp searched before it), and a gapped `tblastn` search now gives NCBI's rows
  for an alignment of at most 11 letters that starts before one of NCBI's 5,000,000-letter parts of
  a reading frame with its middle inside the part, which LOSAT refused. Only the alignments on
  which NCBI BLAST+ crashes (longer than 11 letters, or with the middle before the part) are still
  refused (line NP, M4b-T3).

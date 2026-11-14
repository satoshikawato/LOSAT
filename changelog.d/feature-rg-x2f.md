- `blastx -window_size 0` with a low `-threshold` no longer crashes (a panic, exit status 134) where
  NCBI's one-hit word finder finds an ungapped alignment of negative length, and gives NCBI BLAST+
  2.17.0's output there, with `-subject` and `-db`. Where the gapped start of such an alignment reads
  past NCBI BLAST+'s buffer of the query batch or of the `-subject` protein, or past the last memory
  page of a database volume's sequence file, NCBI BLAST+ reads outside its memory (it can crash, or
  its report can change from run to run), and LOSAT stops with an error instead.

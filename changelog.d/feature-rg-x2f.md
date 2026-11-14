- `blastx -window_size 0` with a low `-threshold` no longer crashes (a panic, exit status 134) where
  NCBI's one-hit word finder finds an ungapped alignment of negative length, and gives NCBI BLAST+
  2.17.0's output there, with `-subject` and `-db`. In the rare searches where the gapped start of
  such an alignment reads past NCBI BLAST+'s buffer of the query batch, of the `-subject` protein or
  of the database volume, NCBI BLAST+ reads outside its memory (it can crash, or its report can
  change from run to run), and LOSAT stops with an error instead.

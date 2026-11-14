- `blastx -ungapped -window_size 0` with a low `-threshold` no longer crashes (a panic) where NCBI's
  one-hit word finder finds ungapped alignments of negative length, and gives NCBI BLAST+ 2.17.0's
  output there, with `-subject` and `-db`: such alignments take part in the linking of the others and
  are then deleted, as in NCBI BLAST+. Where an alignment of length 0 is among the hits that the
  reports show, which NCBI BLAST+ reports in forms of its own for each format, LOSAT stops with an
  error for now.

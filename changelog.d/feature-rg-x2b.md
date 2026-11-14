- `tblastx -window_size 0` with a low `-threshold` now runs the searches it stopped on before, where
  NCBI's one-hit word finder finds an ungapped alignment of length 0 or less, and gives NCBI BLAST+
  2.17.0's output: such an alignment is never reported, but NCBI BLAST+ can join it to a chain of
  other alignments, which changes their E-values, which alignments are reported and their order.
  One case still stops with an error: on a subject frame of more than 5,000,000 letters, two such
  alignments in the overlap of two of its parts whose lengths add up to 0 (NCBI BLAST+'s result is
  undefined there).

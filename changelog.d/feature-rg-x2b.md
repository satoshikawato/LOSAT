- `tblastx -window_size 0` with a low `-threshold` now runs the searches it stopped on before, where
  NCBI's one-hit word finder finds an ungapped alignment of length 0 or less, and gives NCBI BLAST+
  2.17.0's output: such an alignment is never reported, but NCBI BLAST+ can join it to a chain of
  other alignments, which changes their E-values, which alignments are reported and their order.
  One case still stops with an error: on a subject frame of more than 5,000,000 letters, such an
  alignment that meets another in the overlap of two of its parts where their lengths add up to 0
  (the same alignment of length 0 found in both parts, or one of length -1 or -2 that starts or
  ends inside one of length 1 or 2; NCBI BLAST+'s result is undefined there).
- `tblastx` on a subject frame of more than 5,000,000 letters (a subject of more than 15,000,002
  letters; it is searched in overlapping parts) now merges the alignment lists of the parts in
  NCBI BLAST+ 2.17.0's order, so alignments anywhere on that frame with the same score and
  coordinates in different query frames can no longer come out in another order or join another
  chain than NCBI BLAST+'s.

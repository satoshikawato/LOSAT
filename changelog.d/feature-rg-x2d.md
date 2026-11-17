- `tblastn -in_pssm -window_size 0` with a PSSM whose `-` column scores above NCBI BLAST+'s minimum (a
  hand-made PSSM) now gives NCBI BLAST+ 2.17.0's output where a one-hit alignment of negative length has
  its gapped start on the end of its frame's translation or after it and NCBI BLAST+ drops it before its
  traceback (below the score cutoff, for example); before, every such search stopped with an error. Where
  such an alignment reaches the traceback, `tblastn` still stops with an error: after the frame's end
  NCBI BLAST+ reads a byte its translation does not set, and on it NCBI BLAST+ reports a column past the
  subject's end, which LOSAT does not reproduce yet.

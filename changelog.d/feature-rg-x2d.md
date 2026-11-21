- `tblastn -in_pssm -window_size 0` with a PSSM whose `-` column scores above NCBI BLAST+'s minimum (a
  hand-made PSSM) now gives NCBI BLAST+ 2.17.0's output (apart from the version line, which names LOSAT) where a one-hit alignment of negative length has
  its gapped start on the end of its frame's translation or after it and NCBI BLAST+ drops it before its
  reports (below the score cutoff, or past `-max_target_seqs`, for example); before, every such search
  stopped with an error. `tblastn` still stops with an error where NCBI BLAST+'s traceback starts after
  the frame's end (it reads a byte its translation does not set) and where the reports would show an
  alignment whose traceback starts on it (NCBI BLAST+ places one of its columns past the subject's end,
  which LOSAT does not reproduce yet).

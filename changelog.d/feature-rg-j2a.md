- `blastx` and `tblastn` with `-window_size 0` and a low `-threshold` now run the searches they stopped
  on before, where the gapped start of a one-hit alignment of negative length scores a window of 11
  letters that passes NCBI BLAST+'s buffer of the query batch or of the subject (`blastx`: the
  `-subject` protein or, with `-db`, the last memory page of the database volume's sequence file;
  `tblastn`: the subject's six frames). NCBI BLAST+ 2.17.0 reads outside its memory there and can
  crash; LOSAT reads the letters past the sequence as its end marker, as `blastp` does. That gives
  NCBI BLAST+'s own output where those bytes are zero: under valgrind, which reports the reads, NCBI
  BLAST+ wrote LOSAT's output in every such search of ours that LOSAT runs and valgrind finished. `tblastn
  -in_pssm` still stops with an error where such a window passes the PSSM's last position, and
  where a PSSM whose `-` column scores above NCBI BLAST+'s minimum puts the gapped start past the
  end of the subject's translation.

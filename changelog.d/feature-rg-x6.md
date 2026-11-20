- BLASTN: `-task rmblastn` and `-reward 0 -penalty 0` run (matrix-only scoring), as NCBI BLAST+ 2.17.0
  does: the search scores with the default nucleotide matrix (1/-3), reads the X-drops as raw scores,
  keeps preliminary alignments by raw score (`-min_raw_gapped_score`), and format 0 opens with the
  RMBLASTN lines and shows raw scores only, without the list of subjects or the Karlin-Altschul lines.
  It runs with every task NCBI runs it with; `-task megablast` with its greedy extension, where NCBI
  BLAST+ divides by zero or does not finish, is refused, as is a reward of 0 with a penalty other
  than 0 (no query has Karlin-Altschul values).

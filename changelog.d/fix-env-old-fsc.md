- blastp, blastx and tblastn read the environment variable `OLD_FSC` (any value, an empty one too) as NCBI BLAST+
  2.17.0 does: E-values and cutoffs without the finite-size correction, the statistics of `-searchsp 0`. blastx
  ignored it and gave the corrected E-values; blastp and tblastn refused it. psiblast and rpsblast, which take no
  `-searchsp`, still refuse it. Without the variable nothing changes.

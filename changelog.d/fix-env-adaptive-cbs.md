- The environment variable `ADAPTIVE_CBS` (any value, an empty one too) is read as NCBI BLAST+ 2.17.0 reads it
  in blastp, psiblast, rpsblast, blastx and tblastn: with composition-based statistics the preliminary hit list
  holds 1,500 subjects (twice the hit list size and 50 more from 1,000), and with a hit list size below 1,000 it
  is cut by E-value before the final alignments (not in psiblast). blastp kept the larger list but did not cut
  it, blastx ignored the variable, and tblastn and rpsblast refused it. On an input built to show it, the
  reported subjects now match NCBI's (for example 400 where blastp reported 500 and blastx 350). Without the
  variable nothing changes.
- blastx refuses a `-max_target_seqs` whose preliminary hit list size overflows NCBI's 32-bit int
  (1073741799 to 2147483622 with composition-based statistics), with which NCBI BLAST+ crashes, as blastp and
  tblastn already did. Larger sizes, and sizes from 2^30 without composition-based statistics, wrap to NCBI's
  small preliminary size (2147483647: 48 subjects; without the statistics 10) instead of keeping every subject.

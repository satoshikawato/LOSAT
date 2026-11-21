- A hit list size whose preliminary size overflows NCBI BLAST+'s 32-bit int now keeps every hit in blastn,
  blastp, psiblast, rpsblast, blastx, rpstblastn and tblastn: a `-max_target_seqs`, `-num_alignments`, or
  `-num_descriptions` in `-outfmt` 0 to 4, of 2^30 or more, or 1073741799 or more with composition-based
  statistics. NCBI BLAST+ keeps 10 or 2 to 48 subjects per query at these sizes, or crashes; LOSAT reproduced
  the small sizes and refused the others. For example `-max_target_seqs 2147483647`, often used to ask for every
  hit, kept 48 subjects per query in blastp, psiblast, blastx and tblastn and was refused in blastn; it now keeps
  all of them, as NCBI does with 1073741798. With `-best_hit_overhang` or `-best_hit_score_edge`, sizes from
  1073741774 (536870887 with composition-based statistics, 2^30 ungapped) give the result of the size just below,
  where NCBI's overflowing count keeps some alignments the best-hit filter removes.

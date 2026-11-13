- psiblast takes `-sum_stats true` in every round of a FASTA query and with `-in_pssm`, as NCBI
  BLAST+ 2.17.0: the HSPs of each subject are linked (an HSP's E-value is its linked set's, and the
  HSPs are kept from the sum statistics' lower score), outfmt 0 shows a set's size as `Expect(2)`,
  and the alignments a round includes and its PSSM follow the linked E-values (line RG, Q3b).

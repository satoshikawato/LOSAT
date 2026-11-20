- psiblast takes `-sum_stats true` in every round of a FASTA query and with `-in_pssm`, as NCBI
  BLAST+ 2.17.0: the HSPs of each subject are linked (an HSP's E-value is its linked set's, and the
  HSPs are kept from the sum statistics' lower score), outfmt 0 shows a set's size as `Expect(2)`,
  and the alignments a round includes and its PSSM follow the linked E-values (line RG, Q3b).
- blastp and psiblast: a query that NCBI splits into chunks (19,800 letters and more) searches each
  chunk with the effective search spaces of its whole batch, as NCBI does. Without the finite-size
  correction (the environment variable `OLD_FSC`) the chunks' cutoffs read them, so such searches
  now give NCBI's hits (before, more HSPs were kept); with it nothing changes (line RG, Q3b).

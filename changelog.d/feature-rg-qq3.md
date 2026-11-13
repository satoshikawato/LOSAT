- psiblast takes `-culling_limit`, `-gap_trigger`, `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`
  and `-sum_stats false` in every round of a FASTA query and with `-in_pssm`, as NCBI BLAST+ 2.17.0:
  the later rounds cull and convert the X-drops and the gap trigger with the PSSM's statistics, and
  the PSSMs are made of the culled alignments (line RG, Q3). `-sum_stats true`, with which NCBI
  links the HSPs of each subject, is still refused.

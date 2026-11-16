- Fixed `psiblast`'s first round of a FASTA query: it aligned with BLASTP's approximate (restricted)
  gapped alignment, which NCBI BLAST+ uses for `blastp` only, and could miss a hit whose exact
  alignment passes the E-value threshold (and so build another PSSM for the next round). It now
  aligns exactly, as NCBI's `psiblast` does. 3,000 human proteins against the yeast proteome
  (`-evalue 0.001`) reported 12,330 rows, BLASTP's; NCBI BLAST+ 2.17.0 and LOSAT now report 12,356,
  with 0.3 % fewer instructions than before. On generated inputs (distant copies of 4-8 queries among
  500-3,000 random proteins), 16 of 152 command lines (all with `-evalue 1e-3`, `-subject` and `-db`,
  one and three rounds) differed from NCBI BLAST+ 2.17.0 and now give its output at 1, 4 and 8
  threads. `blastp` does not change.

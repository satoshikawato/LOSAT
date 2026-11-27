- Fixed: `tblastn -in_pssm` with `-searchsp` (no Gumbel block) and a PSSM whose gapped H is negative
  now gives NCBI's output: every e-value is -1 (divided by the gap-decay divisor with sum statistics,
  printed 0.0) and the hit-saving cutoff stays 1, so the low-scoring alignments NCBI keeps are no
  longer dropped (LOSAT used to compute both from the block).
- Internal: one function computes the e-value of a single alignment without a Gumbel block (NCBI's
  `BLAST_KarlinStoE_simple`) with NCBI's check of the block (a negative Lambda, K or H gives -1);
  TBLASTN, BLASTX, BLASTP, PSI-BLAST and BLASTN call it in place of their own copies of the formula
  (TBLASTX already did). The port of `BlastKarlinEtoS_simple` that TBLASTN, BLASTX, BLASTP,
  PSI-BLAST and TBLASTX share checks the block as NCBI does; BLASTN keeps its own two copies of
  that code for now (they have the check, and another rule for a quotient past the double). BLASTN
  no longer computes an e-value that NCBI does not compute while it scores a trimmed alignment
  again. Other searches give the same output (line RG, X-evalue-owner).

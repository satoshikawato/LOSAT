- Fixed: `tblastn -in_pssm` with `-searchsp` (no Gumbel block) and a PSSM whose gapped H is negative
  now gives NCBI's output: every e-value is -1 (printed 0.0) and the hit-saving cutoff stays 1, so the
  low-scoring alignments NCBI keeps are no longer dropped (LOSAT used to compute both from the block).
- Internal: one function computes the e-value of a single alignment without a Gumbel block (NCBI's
  `BLAST_KarlinStoE_simple`) with NCBI's check of the block (a negative Lambda, K or H gives -1), and
  the port of `BlastKarlinEtoS_simple` checks the block as NCBI does; TBLASTN, BLASTX, BLASTP,
  PSI-BLAST, BLASTN and TBLASTX call them in place of their own copies. BLASTN no longer computes an
  e-value that NCBI does not compute while it scores a trimmed alignment again. Other searches give
  the same output (line RG, X-evalue-owner).

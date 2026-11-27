- Internal: one function computes the cutoff score of an e-value (NCBI's `BLAST_Cutoffs` with its
  `BlastKarlinEtoS_simple`) for every program: the hit-saving cutoffs, the ungapped word cutoffs
  and the sum statistics' cutoffs of BLASTN, BLASTP, PSI-BLAST, RPS-BLAST, BLASTX, TBLASTN and
  TBLASTX call it in place of their own copies, and it owns the cutoff's start at 1. Each caller
  names its rule for a cutoff past the range of a double (BLASTN NCBI's x86-64 value, the others
  the finite value of NCBI-defect exception 8), as before. The sum statistics' cutoffs carry their
  e-value from one query context to the next as NCBI does, also across the parts of a BLASTX
  query batch of more than 10002 letters (`BATCH_SIZE`) that LOSAT searches apart (line RG,
  X-etos-owner).
- Fixed: blastp and psiblast with an `-evalue` of NaN (`-nan`, which NCBI BLAST+ accepts) and
  `-searchsp` write NCBI's hits: the score cutoff of such an e-value is 1, so every HSP is kept,
  where LOSAT wrote none; with the composition-based redo no match is kept, as in NCBI (rows
  `db.etos.*`; line RG, X-etos-owner).

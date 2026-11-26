- Internal: unified all single-alignment Karlin-Altschul E-value computations (`BLAST_KarlinStoE_simple`)
  across TBLASTN, BLASTX, BLASTN, BLASTP, PSI-BLAST, RPS-BLAST, and TBLASTX onto
  `stats::karlin::evalue_from_raw_score`, adding NCBI's `H < 0.0` guard alongside `Lambda < 0.0`
  and `K < 0.0` so invalid Karlin blocks (including `tblastn -in_pssm` with negative gapped
  `Lambda`, `K`, or `H` and `-searchsp`) return `-1.0` as in NCBI BLAST+ 2.17.0. Output on all
  valid searches is unchanged.

- `psiblast -in_pssm` and `-in_msa` now accept every supported `-matrix`, `-gapopen` and
  `-gapextend` combination (lifting the BLOSUM62 11/1 refusal), write the search gap costs into
  `params.rpsdbparams` for frequency-ratios-only PSSMs, and match NCBI BLAST+ 2.17.0's
  ungapped-before-gapped score-block error order and `tblastn -in_pssm` `PSITBLASTN` first-batch
  error prolog for split-environment errors (`PS7b`).

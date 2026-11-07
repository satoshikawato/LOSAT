- Engine: blastn, blastp, blastx, tblastn, tblastx and `common` under `algorithm/` no longer use
  each other in cycles (psiblast and pssm_search still form one with blastp): the translated-mask
  functions and the translated-chunk schedule move to `common`, `display_codon` to
  `tblastx/report.rs`, and the branch-free scan kernel to `tblastx/x_aascan_bf.rs`. No output
  changes.

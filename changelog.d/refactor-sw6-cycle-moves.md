- Engine: the programs under `algorithm/` no longer use each other in cycles: the translated-mask
  functions and the translated-chunk schedule move to `common`, `display_codon` to
  `tblastx/report.rs`, and the branch-free scan kernel to `tblastx/x_aascan_bf.rs`. No output
  changes.

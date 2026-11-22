# TBLASTX regression fixtures

`tests/tblastx_regression_fixtures.py` lists the cases and freezes NCBI BLAST+ 2.17.0's outputs.

`inputs/tie_nisland_query.fa`, `tie_nisland_subject.fa`, `tie_trim_subject.fa` and `sume_subject.fa`: inputs of the cases
`tie.nisland_seg_no`, `tie.trim_seg_no_fmt0` and `sume.large_gap_seg_no`, from the S08 investigations (commit 2049103b0,
2026-10-02). `tblastx_regression_fixtures.py generate` does not write them; no generator or source of the letters is
recorded.

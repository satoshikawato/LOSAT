- Tests: the three placeholder tests that were always `#[ignore]`d "until reference data is populated"
  are removed (two in `tests/unit/common/evalue_ncbi.rs`; `tests/unit/helpers/ncbi_blast_comparison_test.rs`,
  which no test binary compiled, is deleted). No test that ran is lost.

- Tests: the unit tests of `LOSAT/tests/unit/` (178 tests in 19 files, and one in `tests/unit_tests.rs`) are in the
  crate, each file next to the module it tests in a `#[cfg(test)]` module, so they no longer reach the engine's hidden
  modules from outside the crate (line LIB, P-V2). The four tblastx tests that run the `losat` binary are in the tblastx
  suite (`tests/tblastx_run_engine.rs`). The test binary `unit_tests` is gone; `cargo test --lib` runs these tests.
  `test_diagnostics_enabled` gives its six values of `LOSAT_DIAGNOSTICS` to the value check instead of setting the
  variable, since the lib's tests share one process. The command line's output is unchanged.

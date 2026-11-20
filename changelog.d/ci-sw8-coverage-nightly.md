- CI: Nightly gains the coverage of the test set of ci.yml's `rust-tests`: four `coverage (k/4)`
  jobs run `cargo llvm-cov nextest --all-features` on a quarter of the tests each and export lcov,
  and `coverage report` adds them (`LOSAT/tests/coverage_summary.py`) into the line and function
  coverage of LOSAT's sources, as the job summary with a table by module and the artifact
  `coverage` (lcov, a table per file). Test files are left out; inline `#[cfg(test)]` modules are
  counted. The jobs have no threshold and do not change the run's result; they run in the nightly
  run of dev and in dispatched runs of other branches, not on `main`. `rust lint` runs the script's
  tests.

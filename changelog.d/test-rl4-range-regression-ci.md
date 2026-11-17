- Test tooling: the pull request and nightly regressions (`ci_fast_regressions.py`) run the 75 frozen NCBI range-error
  fixtures of `range_regression_fixtures.py` (`-query_loc` / `-subject_loc`, all five classic programs) whenever one of
  those programs is selected. No workflow ran them before. No product change.

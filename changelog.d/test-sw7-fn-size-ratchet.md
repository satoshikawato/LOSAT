- Tests: `LOSAT/tests/function_size_limit.rs` caps the length of the functions over 1,000 lines at
  their length today and fails when one of them grows past its cap or another function passes
  1,000 lines (test code is not counted); a pull request that grows one on purpose raises its
  number in the test.

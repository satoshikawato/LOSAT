- Tests: the two trace-row tests in `xdrop_simd_tests.rs` take the row pointer through `NonNull::new` before
  writing to it, so CodeQL no longer reports a possible null pointer there (alerts #47 and #48, the code of
  #43 and #44 after the tests moved to their own file). Outputs do not change.

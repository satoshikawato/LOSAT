- The integration tests build into eight test binaries instead of 113: six suites (`blastn`, `blastp`,
  `blastx`, `tblastn`, `tblastx`, `common`) and the two files that change the process environment.
  `cargo test --test blastn_strand` becomes `cargo test --test blastn blastn_strand::`, and nextest
  names a test `<file>::<test>` in the binary `LOSAT::<suite>`. `tests/suite_layout.rs` keeps every
  file in exactly one binary. Outputs do not change.

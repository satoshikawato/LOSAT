- The published `losat` crate leaves out the 118 source files that only unit tests read (`exclude` in
  `LOSAT/Cargo.toml`): 7.84 -> 7.27 MiB, 78.4 % -> 72.7 % of the 10 MiB that crates.io takes (dev of
  2026-11-20). What `cargo build`, the verification of `cargo publish` and docs.rs read is unchanged; the CI
  check `LOSAT/tests/check_crate_test_files.py` keeps the list exact. No change to results, speed or memory.

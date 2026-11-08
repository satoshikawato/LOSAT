- Library: the body of the `losat` binary (the argument handling, the program dispatch and the diagnostic summaries at
  exit) is in the crate as one hidden function without arguments, `main_entry` (`src/cli/main_entry.rs`), which
  `src/main.rs` calls; the binary no longer names the engine's hidden modules (line LIB, P-V3). The command line's
  output, error messages and exit statuses are unchanged.

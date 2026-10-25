- CI: `rust` runs `cargo test --all-features` in two shards on two runners with cargo-nextest (`cargo test (1/2)`,
  `cargo test (2/2)`, `--partition count:k/2`) and its other checks in `rust lint`; the check `rust` passes when all three
  pass. The test profile in CI is opt-level 3 with line tables only (was opt-level 1; debug assertions and overflow
  checks stay on). The platform checks run their tests with cargo-nextest at the same level.

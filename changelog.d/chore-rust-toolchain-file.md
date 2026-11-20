- `rust-toolchain.toml` at the repository root selects Rust 1.92.0 (the CI toolchain) with clippy,
  rustfmt and the two WASI targets, and `LOSAT/Cargo.toml` declares `rust-version = "1.92"`. Builds with
  CI's toolchain are unchanged.

//! The examples of the library API (`examples/`) as the tests of the `common` binary run them.
//! `cargo test` builds every example; a run of one test binary may not have, and then the
//! example is built here, with the profile of the binary under test.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The example `name`, built (with the profile of the binary) when a single test target was
/// asked for.
pub(crate) fn example_binary(name: &str) -> PathBuf {
    let bin_dir = Path::new(env!("CARGO_BIN_EXE_losat"))
        .parent()
        .unwrap()
        .to_path_buf();
    let example = bin_dir
        .join("examples")
        .join(format!("{name}{}", std::env::consts::EXE_SUFFIX));
    if !example.exists() {
        let profile: &[&str] = if bin_dir.ends_with("release") {
            &["--release"]
        } else {
            &["--profile", "test"]
        };
        let status = Command::new(env!("CARGO"))
            .args(["build", "--example", name])
            .args(profile)
            .arg("--target-dir")
            .arg(bin_dir.parent().unwrap())
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .status()
            .expect("build the example");
        assert!(status.success(), "cargo build --example {name}");
    }
    example
}

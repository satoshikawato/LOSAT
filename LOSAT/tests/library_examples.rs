//! The examples of the library API (`examples/`) as the tests of the `common` binary run them.
//! `cargo test` builds every example; a run of one test binary may not have, and an example left
//! by an earlier run does not follow a later source change. So every call builds the example
//! here (a no-op for cargo when it is fresh), with the profile and the features of the binary
//! under test.
//!
//! No NCBI counterpart: test tooling.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The features of `LOSAT/Cargo.toml` `[features]` with whether the test crate, which is
/// compiled with the features of the package, has each on. `default` is one of them: cargo gives
/// the compiler `feature="default"` when it is on (`cargo test --all-features` too), so a build
/// without it is another variant of the library. A new feature goes in this list:
/// `every_feature_of_the_package_is_listed` fails otherwise.
const FEATURES: &[(&str, bool)] = &[
    ("default", cfg!(feature = "default")),
    ("parallel", cfg!(feature = "parallel")),
    ("xstats", cfg!(feature = "xstats")),
    ("diagnostics", cfg!(feature = "diagnostics")),
    ("wasm-threads", cfg!(feature = "wasm-threads")),
    ("blastx-worker-probe", cfg!(feature = "blastx-worker-probe")),
    ("tblastx-wasm-scalar", cfg!(feature = "tblastx-wasm-scalar")),
];

/// The example `name`, built by `cargo build --example` on every call, with the profile of the
/// binary under test and the features that binary has on (`--no-default-features --features
/// <those>`, so that a `cargo test --all-features` run neither rebuilds the library with the
/// default features nor replaces the example by another variant). A file that already exists
/// is not enough: it can be older than the source (a run of one test target after a change
/// would run it, and a stale example can pass as well as fail). `--locked` as CI.
pub(crate) fn example_binary(name: &str) -> PathBuf {
    let bin_dir = Path::new(env!("CARGO_BIN_EXE_losat"))
        .parent()
        .unwrap()
        .to_path_buf();
    let example = bin_dir
        .join("examples")
        .join(format!("{name}{}", std::env::consts::EXE_SUFFIX));
    let profile: &[&str] = if bin_dir.ends_with("release") {
        &["--release"]
    } else {
        &["--profile", "test"]
    };
    let on: Vec<&str> = FEATURES
        .iter()
        .filter(|(_, on)| *on)
        .map(|(feature, _)| *feature)
        .collect();
    let mut command = Command::new(env!("CARGO"));
    command
        .args(["build", "--example", name])
        .args(profile)
        .arg("--no-default-features");
    if !on.is_empty() {
        command.arg("--features").arg(on.join(","));
    }
    let status = command
        .arg("--locked")
        .arg("--target-dir")
        .arg(bin_dir.parent().unwrap())
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .status()
        .expect("build the example");
    assert!(status.success(), "cargo build --example {name}");
    example
}

/// The names of the features in the `[features]` table of `Cargo.toml`.
fn package_features(manifest: &str) -> Vec<&str> {
    let mut in_table = false;
    let mut names = Vec::new();
    for line in manifest.lines() {
        let line = line.trim();
        if line.starts_with('[') {
            in_table = line == "[features]";
        } else if in_table && !line.starts_with('#') {
            if let Some((key, _)) = line.split_once('=') {
                names.push(key.trim());
            }
        }
    }
    names
}

#[test]
fn every_feature_of_the_package_is_listed() {
    let manifest =
        std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml"))
            .expect("read Cargo.toml");
    let in_toml = package_features(&manifest);
    assert!(
        !in_toml.is_empty(),
        "no [features] table found in Cargo.toml"
    );
    for feature in &in_toml {
        assert!(
            FEATURES.iter().any(|(name, _)| name == feature),
            "feature `{feature}` of Cargo.toml is missing from FEATURES in tests/library_examples.rs; \
             example_binary would build the example without it"
        );
    }
    for (name, _) in FEATURES {
        assert!(
            in_toml.contains(name),
            "FEATURES lists `{name}`, which Cargo.toml no longer has"
        );
    }
}

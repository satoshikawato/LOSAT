//! The examples of the library API (`examples/`) as the tests of the `common` binary run them.
//! `cargo test` builds every example; a run of one test binary may not have, and an example left
//! by an earlier run does not follow a later source change. So the first call of a run of the
//! tests builds the example here (a no-op for cargo when it is fresh), with the profile and the
//! features of the binary under test, and the later calls of that run use it.
//!
//! Once per run, not once per call: cargo puts the example in place again at every build, fresh
//! or not, and on macOS it does so by removing the file and copying it back (elsewhere the file
//! is a hard link that is already in place, and nothing is done), so a build in one test removed
//! the example under a test that was starting it (os error 2, Platform checks 45315823856 on
//! macOS arm64 and x64).
//! cargo 0.93.0: src/cargo/core/compiler/mod.rs:233-234 (the fresh unit's link step)
//! cargo 0.93.0: crates/cargo-util/src/paths.rs:610-675 (`_link_or_copy`)
//!
//! No NCBI counterpart: test tooling.

use std::fs::File;
use std::io::{Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;
use std::time::{SystemTime, UNIX_EPOCH};

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

/// The run of the tests that this process belongs to: nextest's (`NEXTEST_RUN_ID`, the same in
/// the process of each test of one `cargo nextest run`), else this process's (`cargo test` runs
/// the tests of a binary as threads of one process), told from an earlier process of the same id
/// by the time of its first call.
fn test_run() -> &'static str {
    static RUN: OnceLock<String> = OnceLock::new();
    RUN.get_or_init(|| match std::env::var("NEXTEST_RUN_ID") {
        Ok(id) if !id.is_empty() => format!("nextest run {id}"),
        _ => {
            let since = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default();
            format!(
                "process {} from {} ns",
                std::process::id(),
                since.as_nanos()
            )
        }
    })
}

/// The example `name`, built by `cargo build --example` at the first call of a run of the tests
/// (`test_run`), with the profile of the binary under test and the features that binary has on
/// (`--no-default-features --features <those>`, so that a `cargo test --all-features` run neither
/// rebuilds the library with the default features nor replaces the example by another variant).
/// A file that already exists is not enough: it can be older than the source (a run of one test
/// target after a change would run it, and a stale example can pass as well as fail).
/// `--locked` as CI.
///
/// The calls of all processes take turns on a lock file in `CARGO_TARGET_TMPDIR`, which records
/// the run and the build that made the example; a call that finds its own there builds nothing,
/// so no build of the run moves the example while a test starts it.
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
    let built = format!("{}\n{}\n{}\n", test_run(), bin_dir.display(), on.join(","));
    let tmp = Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(tmp).expect("make CARGO_TARGET_TMPDIR");
    let profile_dir = bin_dir.file_name().unwrap().to_string_lossy();
    let mut record = File::options()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(tmp.join(format!("example-{name}-{profile_dir}.built")))
        .expect("open the record of the example's build");
    record
        .lock()
        .expect("lock the record of the example's build");
    let mut last = String::new();
    // A record that is not UTF-8 is no build of this run.
    let _ = record.read_to_string(&mut last);
    if last == built {
        return example;
    }
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
    record
        .set_len(0)
        .expect("empty the record of the example's build");
    record
        .rewind()
        .expect("rewind the record of the example's build");
    record
        .write_all(built.as_bytes())
        .expect("write the record of the example's build");
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

//! `LOSAT_X_REACH` takes only an absolute directory (`src/utils/xreach.rs`, `reach_directory`).
//! A relative value would be a directory under the run's current directory: on 2026-10-31 runs
//! with the values `0` and `1` left `LOSAT/0/` and `LOSAT/1/` in a checkout (12 files, removed by
//! #715). A run given such a value writes no reach file, says why once on standard error, and
//! reports what it reports without the variable.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A blastn search of LC738874 against LC738868 (`-outfmt 6`, two threads, so the decision
/// points also run on pool threads) in `cwd`, with `reach` as its only `LOSAT_*` variable besides
/// `LOSAT_BLASTX_WORKER_LOG` (which a `--all-features` binary needs and which selects nothing).
fn search(cwd: &Path, reach: Option<&str>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    if let Some(value) = reach {
        command.env("LOSAT_X_REACH", value);
    }
    command
        .current_dir(cwd)
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(["-outfmt", "6", "-num_threads", "2"])
        .output()
        .expect("run LOSAT")
}

#[test]
fn a_relative_reach_directory_is_refused_once_and_writes_nothing() {
    let cwd = std::env::temp_dir().join(format!("losat-x-reach-relative-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&cwd);
    std::fs::create_dir_all(&cwd).expect("working directory");
    let plain = search(&cwd, None);
    assert!(
        plain.status.success(),
        "{}",
        String::from_utf8_lossy(&plain.stderr)
    );
    assert!(!plain.stdout.is_empty(), "no hits");
    for value in ["0", "1", "reach/out"] {
        let run = search(&cwd, Some(value));
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert!(run.status.success(), "LOSAT_X_REACH={value}: {stderr}");
        assert!(
            run.stdout == plain.stdout,
            "LOSAT_X_REACH={value}: the report differs"
        );
        let line = format!("[X_REACH] ignored: LOSAT_X_REACH={value} is not an absolute path\n");
        assert_eq!(
            stderr.matches(&line).count(),
            1,
            "LOSAT_X_REACH={value}: {stderr}"
        );
        assert_eq!(
            stderr.replace(&line, ""),
            String::from_utf8_lossy(&plain.stderr),
            "LOSAT_X_REACH={value}: other standard error"
        );
    }
    let left: Vec<_> = std::fs::read_dir(&cwd)
        .expect("working directory")
        .map(|entry| entry.expect("entry").file_name())
        .collect();
    assert!(
        left.is_empty(),
        "the runs left {left:?} in {}",
        cwd.display()
    );
    let _ = std::fs::remove_dir_all(&cwd);
}

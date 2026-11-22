//! dc-megablast's CSR table builds with `LOSAT_X_DCQSAHEAD` (`src/algorithm/blastn/x_dc_csr.rs`
//! `build_with`, `src/algorithm/blastn/x_disc_lookup_par.rs` `csr_on_pool_steps`) on whole
//! searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes when the place pass of a CSR build brings in the lines of `qs`, not what it writes
//! where, so the property under test is that the report does not depend on it, and that the
//! default's builds took its path (`[TIMING] disc_lookup_csr: ... qs_ahead=`, which counts them).
//! MG1655 (4.6 Mb, one table of about 9.2 million words, above the CSR's size rule) against
//! AP027131, built by the serial build at one thread and on the pool at four. That those bytes
//! equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DCQSAHEAD", "0");
const TIMING: (&str, &str) = ("LOSAT_TIMING", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One dc-megablast search (`-outfmt 6` on standard output) with the given variables as its only
/// `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. Returns the report and the
/// standard error.
fn search(threads: u32, vars: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .args(["-task", "dc-megablast"])
        .arg("-query")
        .arg(fasta("MG1655.fna"))
        .arg("-subject")
        .arg(fasta("AP027131.fasta"))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// `name=<n>` of the `[TIMING] disc_lookup_csr: ` line of `stderr`.
fn csr_field(stderr: &str, name: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] disc_lookup_csr: "))
        .unwrap_or_else(|| panic!("no CSR timing line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|word| word.strip_prefix(name))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
}

// At one thread (the serial build) and four (the build on the pool): the default and `=0` give the
// same report and standard error; with LOSAT_TIMING every CSR build of the default took the
// switch's path and none of `=0`'s did.
#[test]
fn dc_megablast_csr_builds_with_qs_ahead() {
    for threads in [1, 4] {
        let (expected, expected_stderr) = search(threads, &[SWITCH_OFF]);
        assert!(!expected.is_empty(), "threads {threads}: no hits");
        let (report, stderr) = search(threads, &[]);
        assert!(report == expected, "threads {threads}: default differs");
        assert_eq!(stderr, expected_stderr, "threads {threads}");
        let (report, stderr) = search(threads, &[TIMING]);
        assert!(
            report == expected,
            "threads {threads}: default with LOSAT_TIMING differs"
        );
        let builds = csr_field(&stderr, "builds=");
        assert!(builds > 0, "threads {threads}: no CSR build:\n{stderr}");
        let parallel = if threads > 1 { builds } else { 0 };
        assert_eq!(
            csr_field(&stderr, "parallel="),
            parallel,
            "threads {threads}"
        );
        assert_eq!(csr_field(&stderr, "qs_ahead="), builds, "threads {threads}");
        let (report, stderr) = search(threads, &[SWITCH_OFF, TIMING]);
        assert!(
            report == expected,
            "threads {threads}: =0 with LOSAT_TIMING differs"
        );
        assert_eq!(csr_field(&stderr, "builds="), builds, "threads {threads}");
        assert_eq!(csr_field(&stderr, "qs_ahead="), 0, "threads {threads}");
    }
}

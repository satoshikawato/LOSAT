//! The NEON build of the vector X-drop kernels (`LOSAT_X_DPNEON`, `src/utils/xdrop_simd.rs`) on
//! whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. On aarch64 it lets `LOSAT_X_DPFAST`
//! run the vector kernels; on the other targets it does nothing, so these tests mean most on the
//! arm64 runners of the platform checks, which run `cargo test`. It chooses between two LOSAT
//! code paths with the same result, so the property under test is that the report does not
//! depend on it: the default run (the NEON kernels on aarch64), `LOSAT_X_DPNEON=0` (the scalar
//! X-drop with the other switches on), `LOSAT_X_DEFAULT=0` and the shadow run (which runs the
//! kernels and the scalar port on every call and stops on the first difference) give the same
//! bytes, and no run prints anything on standard error. That those bytes equal NCBI's is the
//! business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DPNEON", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_DPSHADOW", "1");

/// A blastp search (gapped DP in 8- and 16-bit lanes, with traceback) and a blastn search with
/// `-task blastn` (the X-drop DP of nucleotide matches), both from the platform-check cases.
const SEARCHES: [(&str, &str, &str, &[&str]); 2] = [
    ("blastp", "WSSV.faa", "CoBV.faa", &[]),
    (
        "blastn",
        "LC738874.fasta",
        "LC738870.fasta",
        &["-task", "blastn"],
    ),
];

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One search (`-outfmt 6` on standard output) with the given variables as its only `LOSAT_*`
/// variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. Returns the report and the
/// standard error.
fn search(
    (program, query, subject, extra): (&str, &str, &str, &[&str]),
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg(program)
        .arg("-query")
        .arg(fasta(query))
        .arg("-subject")
        .arg(fasta(subject))
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{program} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

// The default search runs the NEON kernels on aarch64; `LOSAT_X_DPNEON=0` and `LOSAT_X_DEFAULT=0`
// run the scalar port there. The reports are equal, and no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_the_neon_kernels() {
    for case in SEARCHES {
        for threads in [1, 4] {
            let (default, default_stderr) = search(case, threads, &[]);
            assert!(!default.is_empty(), "{} threads {threads}: no hits", case.0);
            assert_eq!(default_stderr, "", "{} threads {threads}", case.0);
            for off in [SWITCH_OFF, DEFAULT_OFF] {
                let (report, stderr) = search(case, threads, &[off]);
                assert_eq!(report, default, "{} {off:?} threads {threads}", case.0);
                assert_eq!(stderr, "", "{} {off:?} threads {threads}", case.0);
            }
        }
    }
}

// `LOSAT_X_DEFAULT=0` with `LOSAT_X_DPFAST` and the variable set to a value other than 0 runs the
// kernels, as the variable did before the switch joined the set. The shadow (set to any value,
// not in the set) runs and compares the kernels on every call, whatever `LOSAT_X_DPNEON` says,
// and does not change the report.
#[test]
fn the_shadow_and_the_explicit_switch_leave_the_report_unchanged() {
    for case in SEARCHES {
        let (expected, _) = search(case, 1, &[]);
        let dpfast = ("LOSAT_X_DPFAST", "1");
        let neon = ("LOSAT_X_DPNEON", "1");
        for switches in [
            &[DEFAULT_OFF, dpfast, neon][..],
            &[SHADOW][..],
            &[SWITCH_OFF, SHADOW][..],
        ] {
            let (report, _) = search(case, 1, switches);
            assert_eq!(report, expected, "{} {switches:?}", case.0);
        }
    }
}

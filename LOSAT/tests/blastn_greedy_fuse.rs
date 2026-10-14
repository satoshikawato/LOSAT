//! The fused loop of megablast's score-only greedy alignment (`LOSAT_X_GREEDYFUSE`,
//! `src/algorithm/blastn/alignment/x_greedy_fused.rs`) on a whole search.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. It chooses between two LOSAT code
//! paths with the same result, so the property under test is that the report does not depend
//! on it: the default run, `LOSAT_X_GREEDYFUSE=0`, `LOSAT_X_DEFAULT=0` and the shadow run (which
//! runs both loops on every call and stops on the first difference) give the same bytes, and
//! the default run prints nothing on standard error (the `[X_SHADOW] GREEDYFUSE` line is
//! printed only with `LOSAT_X_GREEDYFUSESHADOW`). That those bytes equal NCBI's is the
//! business of `megablast_xdrop.rs` and of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_GREEDYFUSE", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_GREEDYFUSESHADOW", "1");

fn input() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/NZ_CP006932.fasta")
}

/// One megablast search of the fixture against itself (`-outfmt 6` on standard output) with the
/// given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more
/// than one thread, and it selects nothing. Returns the report and the standard error.
fn megablast(threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(input())
        .arg("-subject")
        .arg(input())
        .args(["-task", "megablast"])
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

// The default search runs the fused loop; `LOSAT_X_GREEDYFUSE=0` and `LOSAT_X_DEFAULT=0` run the
// ported loop. The reports are equal, and no run prints a line on standard error.
#[test]
fn the_report_does_not_depend_on_the_fused_loop() {
    for threads in [1, 4] {
        let (default, default_stderr) = megablast(threads, &[]);
        assert!(!default.is_empty(), "threads {threads}: no hits");
        assert_eq!(default_stderr, "", "threads {threads}");
        for off in [SWITCH_OFF, DEFAULT_OFF] {
            let (report, stderr) = megablast(threads, &[off]);
            assert_eq!(report, default, "{off:?} threads {threads}");
            assert_eq!(stderr, "", "{off:?} threads {threads}");
        }
    }
}

// `LOSAT_X_DEFAULT=0` with the variable set to a value other than 0 runs the fused loop, as the
// variable did before the switch joined the set. The shadow (set to any value, not in the set)
// compares the two loops on every score-only call, prints its summary line, and does not change
// the report. The calls it compares are those of the preliminary stage with
// `LOSAT_X_GREEDYPRELIM` (in the set), so under `LOSAT_X_DEFAULT=0` that variable is set too.
#[test]
fn the_shadow_compares_the_two_loops_and_the_report_is_unchanged() {
    let (expected, _) = megablast(1, &[]);
    let (report, _) = megablast(1, &[DEFAULT_OFF, ("LOSAT_X_GREEDYFUSE", "1")]);
    assert_eq!(report, expected);
    let prelim = ("LOSAT_X_GREEDYPRELIM", "1");
    for switches in [&[SHADOW][..], &[DEFAULT_OFF, prelim, SHADOW][..]] {
        let (report, stderr) = megablast(1, switches);
        assert_eq!(report, expected, "{switches:?}");
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[X_SHADOW] GREEDYFUSE compared="))
            .unwrap_or_else(|| panic!("{switches:?}: no summary line:\n{stderr}"));
        assert!(
            !line.starts_with("[X_SHADOW] GREEDYFUSE compared=0 "),
            "{switches:?}: {line}"
        );
    }
}

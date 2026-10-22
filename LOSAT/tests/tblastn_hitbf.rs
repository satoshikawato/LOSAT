//! The branch-free two-hit loop of tblastn's unit path (`LOSAT_X_TBNHITBF`,
//! `src/algorithm/tblastn/x_tbn_hitbf.rs`, module doc) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. It acts on every offset array the unit path's scan fills;
//! its shadow (`LOSAT_X_TBNHITBFSHADOW`) turns it on and compares it with the loop before the switch
//! on every array.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-551
//! ```c
//!                 if (diff >= window) {
//!                     /* We are beyond the window for this diagonal; start a
//!                        new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//!                 }
//!
//!                 /* If the difference is less than the wordsize (i.e. last
//!                    hit and this hit overlap), give up */
//!
//!                 if (diff < wordsize) {
//!                     continue;
//!                 }
//! ```
//! The loop takes these two outcomes without a branch and hands the others to the per-hit body,
//! with the same cell states and the same saved init HSPs (the module's argument), so the property
//! under test is that the report does not depend on it. That reports equal NCBI's is the business of
//! the fixture comparisons. The input (the first 40 proteins of AvCLPV against its genome) gives
//! extensions, flag-1 cells and HSPs. The tests compare the default, cleared tables per unit
//! (`LOSAT_X_TBNSTATE=0`), the switch alone in the unit path (`LOSAT_X_DEFAULT=0` with
//! `LOSAT_X_TBNPAR=1`) and the shadow (which also shows that the loop ran on every offset array and
//! that none fell outside its range) with `=0`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const OFF: (&str, &str) = ("LOSAT_X_TBNHITBF", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBNHITBFSHADOW", "1");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
/// A cleared table for every unit instead of the kept one (`XTwoHit::reset`).
const CLEARED: (&str, &str) = ("LOSAT_X_TBNSTATE", "0");
/// The unit path (`LOSAT_X_DEFAULT=0` turns it off with the rest of the strict set).
const UNIT_PATH: (&str, &str) = ("LOSAT_X_TBNPAR", "1");
const ON: (&str, &str) = ("LOSAT_X_TBNHITBF", "1");
/// The `[X_STATS] TBNSCANBF` and `[X_STATS] TBNSCANBYTE` lines: the scan calls of the unit path,
/// one offset array each.
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
/// The proteins of the query.
const PROTEINS: usize = 40;

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(switches.iter().copied());
    command
}

/// The query file (the first `PROTEINS` records of AvCLPV.faa), in a directory of the test's own.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-tblastn-hitbf-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.faa");
        let text = std::fs::read_to_string(fasta("AvCLPV.faa")).unwrap();
        let end = text
            .match_indices('>')
            .nth(PROTEINS)
            .map_or(text.len(), |(at, _)| at);
        std::fs::write(&query, &text[..end]).unwrap();
        Self { dir, query }
    }

    /// One tblastn search of the AvCLPV genome on standard output (`-outfmt 6`); returns the report
    /// and the standard error.
    fn tblastn(&self, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
        let output = losat(switches)
            .arg("tblastn")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(fasta("AvCLPV.fasta"))
            .args(["-outfmt", "6"])
            .args(["-num_threads", &threads.to_string()])
            .output()
            .expect("run LOSAT");
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(
            output.status.success(),
            "{switches:?} threads {threads}: {stderr}"
        );
        (output.stdout, stderr)
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The lines of `stderr` that start with `prefix`.
fn lines<'a>(stderr: &'a str, prefix: &str) -> Vec<&'a str> {
    stderr
        .lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

/// The scan calls of the unit path: those of the `[X_STATS] TBNSCANBF` line (both loops) and of
/// the `[X_STATS] TBNSCANBYTE` line (LOSAT_X_TBNSCANBYTE), each line counted where it is printed.
fn scan_calls(stderr: &str) -> u64 {
    let mut calls = 0;
    if !lines(stderr, "[X_STATS] TBNSCANBF").is_empty() {
        calls += field(stderr, "[X_STATS] TBNSCANBF", "calls_bf")
            + field(stderr, "[X_STATS] TBNSCANBF", "calls_ref");
    }
    if !lines(stderr, "[X_STATS] TBNSCANBYTE").is_empty() {
        calls += field(stderr, "[X_STATS] TBNSCANBYTE", "calls");
    }
    assert!(calls > 0, "no scan line in {stderr}");
    calls
}

/// The number after `name=` in the first line of `stderr` that starts with `prefix`.
fn field(stderr: &str, prefix: &str, name: &str) -> u64 {
    let line = lines(stderr, prefix)
        .first()
        .copied()
        .unwrap_or_else(|| panic!("no {prefix} line in {stderr}"));
    let key = format!("{name}=");
    let at = line
        .find(&key)
        .unwrap_or_else(|| panic!("no {name} in {line}"))
        + key.len();
    line[at..]
        .split(|c: char| !c.is_ascii_digit())
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

#[test]
fn the_branch_free_loop_gives_the_same_report() {
    let inputs = Inputs::new("same");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastn(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        for switches in [&[][..], &[CLEARED][..]] {
            let (report, _) = inputs.tblastn(threads, switches);
            assert_eq!(report, off, "threads {threads} {switches:?}: against =0");
        }
        let (cleared_off, _) = inputs.tblastn(threads, &[CLEARED, OFF]);
        assert_eq!(cleared_off, off, "threads {threads}: cleared =0 against =0");
    }
}

#[test]
fn the_shadow_compares_every_offset_array() {
    let inputs = Inputs::new("shadow");
    for threads in [1, 2] {
        let (off, _) = inputs.tblastn(threads, &[OFF]);
        for switches in [&[SHADOW][..], &[SHADOW, OFF][..], &[SHADOW, CLEARED][..]] {
            let switches = &[switches, &[STATS][..]].concat();
            let (shadow, stderr) = inputs.tblastn(threads, switches);
            assert_eq!(
                shadow, off,
                "threads {threads} {switches:?}: shadow against =0"
            );
            assert!(
                field(&stderr, "[X_SHADOW] TBNHITBF", "arrays_compared") > 0,
                "threads {threads} {switches:?}: no array compared: {stderr}"
            );
            assert!(
                field(&stderr, "[X_SHADOW] TBNHITBF", "hits_compared") > 0,
                "threads {threads} {switches:?}: no hit compared: {stderr}"
            );
            assert_eq!(
                field(&stderr, "[X_SHADOW] TBNHITBF", "outside"),
                0,
                "threads {threads} {switches:?}: an array outside the range: {stderr}"
            );
            assert_eq!(
                field(&stderr, "[X_SHADOW] TBNHITBF", "arrays_compared"),
                scan_calls(&stderr),
                "threads {threads} {switches:?}: one offset array per scan call: {stderr}"
            );
        }
    }
}

// No NCBI counterpart: where the switch acts. With `LOSAT_X_DEFAULT=0` the unit path is off and so
// is the switch: no loop of the unit path runs (no shadow line). With `LOSAT_X_TBNPAR=1` as well the
// unit path runs again with the loop before the switch, unless LOSAT_X_TBNHITBF=1 or the shadow
// turns it on.
#[test]
fn the_switch_acts_in_the_unit_path_only() {
    let inputs = Inputs::new("unit");
    let (off, _) = inputs.tblastn(1, &[OFF]);
    for switches in [
        &[DEFAULT_OFF][..],
        &[DEFAULT_OFF, SHADOW][..],
        &[DEFAULT_OFF, UNIT_PATH][..],
        &[DEFAULT_OFF, UNIT_PATH, ON][..],
    ] {
        let (report, stderr) = inputs.tblastn(1, switches);
        assert_eq!(report, off, "{switches:?}: against =0");
        assert!(
            lines(&stderr, "[X_SHADOW] TBNHITBF").is_empty(),
            "{switches:?}: {stderr}"
        );
    }
    let (report, stderr) = inputs.tblastn(1, &[DEFAULT_OFF, UNIT_PATH, SHADOW]);
    assert_eq!(report, off, "the unit path with the shadow");
    assert!(field(&stderr, "[X_SHADOW] TBNHITBF", "arrays_compared") > 0);
    assert!(field(&stderr, "[X_SHADOW] TBNHITBF", "hits_compared") > 0);
}

#[test]
fn the_switch_off_prints_no_shadow_line() {
    let inputs = Inputs::new("quiet");
    for switches in [&[][..], &[OFF][..]] {
        let (_, stderr) = inputs.tblastn(1, switches);
        assert!(
            lines(&stderr, "[X_SHADOW] TBNHITBF").is_empty(),
            "{switches:?}: {stderr}"
        );
    }
}

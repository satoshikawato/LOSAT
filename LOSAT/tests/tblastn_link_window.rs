//! tblastn's linking candidates chosen by subject position (`src/algorithm/tblastn/stage_d_linking/x_link_window.rs`;
//! the former switch `LOSAT_X_TBNLINKWIN` graduated on 2026-10-26) on whole searches.
//!
//! NCBI reference: c++/src/algo/blast/core/link_hsps.c:1713-1723
//! ```c
//!        for (index1 = hsp_index_left; index1 < hsp_index_right; ++index1) {
//!            BlastLinkedHSPSet* lhsp = offset_hsp_array[index1];
//!            ...
//!            if (s_LinkedHSPSetsAdmissible(head_hsp, lhsp,
//!                                          link_hsp_params, program)) {
//! ```
//! The search tests the candidates of NCBI's window less those that cannot be admissible, in
//! NCBI's order (the unit test `x_chosen_candidates_give_the_links_of_the_whole_window` compares
//! the links with NCBI's whole window), so its report is that of NCBI's loop (the frozen NCBI
//! hashes of the fast regressions check that). These tests compare the reports at one and four
//! threads, and with `LOSAT_X_TBNLINKSETS` (the sets of two or more come from the subject lookup)
//! on and off at one and two threads, and read the
//! `[X_STATS] TBNLINKWIN heads= window= chosen= sets=` line (`LOSAT_X_STATS=1`, printed when a
//! list was linked).
//!
//! The queries are the first protein of `tests/fasta/PemoMJNVA.faa` and all 115 of them; the
//! subject is `tests/fasta/MelaMJNV.fasta`.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The two query files: the first protein of `PemoMJNVA.faa` (written once per test process),
/// and the whole file.
fn queries() -> &'static [PathBuf; 2] {
    static FILES: std::sync::OnceLock<[PathBuf; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| {
        let all = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PemoMJNVA.faa");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_link_window");
        std::fs::create_dir_all(&dir).expect("work directory");
        let one = dir.join(format!("one.{}.faa", std::process::id()));
        let text = std::fs::read_to_string(&all).expect("query file");
        let second = text[1..].find('>').map_or(text.len(), |at| at + 1);
        std::fs::write(&one, &text[..second]).expect("query file");
        [one, all]
    })
}

/// One TBLASTN run against `MelaMJNV.fasta` (`-outfmt 6` on standard output) with `variables`
/// and `LOSAT_X_STATS` as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build
/// with the worker probe (`--all-features`) writes to that file as soon as a pool starts, so the
/// caller's value is kept, or a file of the tests is named. Returns the report and the standard
/// error.
fn tblastn(query: &Path, threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/MelaMJNV.fasta"))
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .env("LOSAT_X_STATS", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// The fields of the one `[X_STATS] TBNLINKWIN` line; `None` without the line.
fn stats(stderr: &str) -> Option<[u64; 4]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNLINKWIN "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(["heads", "window", "chosen", "sets"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

// NCBI reference: c++/src/algo/blast/core/link_hsps.c:1713-1723 (above). The report does not
// depend on the threads; the line is there, and fewer candidates were tested than the windows held.
#[test]
fn the_report_is_the_same_at_every_thread_count_and_the_lookup_ran() {
    for query in queries() {
        let (expected, stderr) = tblastn(query, 1, &[]);
        assert!(
            !expected.is_empty(),
            "{query:?}: the search must report hits"
        );
        let [heads, window, chosen, _] = stats(&stderr).expect("a line");
        assert!(heads > 0 && chosen <= window, "{query:?}");
        let (report, _) = tblastn(query, 4, &[]);
        assert!(report == expected, "{query:?}: four threads");
    }
}

// NCBI reference: c++/src/algo/blast/core/link_hsps.c:1713-1723 (above). With the switch on, the
// sets of two or more come from the subject lookup (the whole file of queries forms sets that the
// lookup meets); with `=0` every set in the window is a candidate, as before. The reports are the
// same at one and two threads.
#[test]
fn the_sets_from_the_lookup_give_the_report_of_every_set() {
    let query = &queries()[1];
    for threads in [1, 2] {
        let (expected, stderr) = tblastn(query, threads, &[("LOSAT_X_TBNLINKSETS", "0")]);
        assert!(!expected.is_empty(), "the search must report hits");
        let [_, _, chosen_off, sets_off] = stats(&stderr).expect("a line");
        assert_eq!(sets_off, 0, "{threads} threads, =0: no set from the lookup");
        let (report, stderr) = tblastn(query, threads, &[]);
        assert!(report == expected, "{threads} threads: the switch on");
        let [_, _, chosen, sets] = stats(&stderr).expect("a line");
        assert!(
            sets > 0 && chosen <= chosen_off,
            "{threads} threads: the lookup must bring sets in ({sets}), no more candidates \
             ({chosen} vs {chosen_off})"
        );
    }
}

//! tblastn's linking candidates chosen by subject position (`LOSAT_X_TBNLINKWIN`, in the strict
//! set; `src/algorithm/tblastn/stage_d_linking/x_link_window.rs`) on whole searches.
//!
//! NCBI reference: c++/src/algo/blast/core/link_hsps.c:1713-1723
//! ```c
//!        for (index1 = hsp_index_left; index1 < hsp_index_right; ++index1) {
//!            BlastLinkedHSPSet* lhsp = offset_hsp_array[index1];
//!            ...
//!            if (s_LinkedHSPSetsAdmissible(head_hsp, lhsp,
//!                                          link_hsp_params, program)) {
//! ```
//! The switch tests the candidates of NCBI's window less those that cannot be admissible, in
//! NCBI's order (the unit test `x_chosen_candidates_give_the_links_of_the_whole_window` compares
//! the links), so a search with it gives the report of NCBI's loop (`LOSAT_X_TBNLINKWIN=0`). It is
//! on by default and off with `=0`. These
//! tests compare the reports and read the `[X_STATS] TBNLINKWIN heads= window= chosen=
//! shadow_heads=` line (`LOSAT_X_STATS=1`, printed when a list was linked with the switch).
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
    Some(["heads", "window", "chosen", "shadow_heads"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

const OFF: (&str, &str) = ("LOSAT_X_TBNLINKWIN", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBNLINKWINSHADOW", "1");

// NCBI reference: c++/src/algo/blast/core/link_hsps.c:1713-1723 (above). The default, the shadow
// (every head's whole window tested as well) and four threads give the report of NCBI's loop
// (`=0`); the line is there with the switch and not without it, and fewer candidates were tested
// than the windows held.
#[test]
fn the_report_is_that_of_ncbis_loop() {
    for query in queries() {
        let (expected, stderr) = tblastn(query, 1, &[OFF]);
        assert!(
            !expected.is_empty(),
            "{query:?}: the search must report hits"
        );
        assert_eq!(stats(&stderr), None, "{query:?}: no line with =0");
        let (report, stderr) = tblastn(query, 1, &[]);
        assert!(
            report == expected,
            "{query:?}: the default changes the report"
        );
        let [heads, window, chosen, shadow] = stats(&stderr).expect("a line by default");
        assert!(heads > 0 && chosen <= window, "{query:?}");
        assert_eq!(shadow, 0, "{query:?}");
        let (report, stderr) = tblastn(query, 1, &[SHADOW]);
        assert!(report == expected, "{query:?}: shadow");
        let [heads, .., shadow] = stats(&stderr).expect("a line with the shadow");
        assert!(heads > 0 && shadow > 0, "{query:?}");
        for variables in [&[][..], &[OFF]] {
            let (report, _) = tblastn(query, 4, variables);
            assert!(report == expected, "{query:?} {variables:?}: four threads");
        }
    }
}

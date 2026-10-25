//! tblastn's lookup scan without a branch on the presence vector (`src/algorithm/tblastn/x_scan_bf.rs`;
//! the former switch `LOSAT_X_TBNSCANBF` graduated on 2026-10-26) on whole searches.
//!
//! NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127
//! ```c
//! for (s = s_first; s <= s_last; s++) {
//!     index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
//!     if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
//!         numhits = bbc[index].num_used;
//!         if (numhits <= (array_size - totalhits)) {
//!             ...
//!         } else {
//!             s_range[1] = (Int4)(s - subject->sequence);
//!             return totalhits;
//!         }
//!     }
//! }
//! ```
//! The loop returns NCBI's offset pairs and scan range on every call (the unit test
//! `x_scan_bf_equals_ncbi_scan_on_random_subjects` compares the calls), so a search gives the
//! report of NCBI's loop (the frozen NCBI hashes of the fast regressions check that). It acts in
//! the unit path, which every search takes. These tests compare the reports at one and four
//! threads and read the `[X_STATS] TBNSCANBF units= calls=` line (`LOSAT_X_STATS=1`) to check that
//! the loop scanned and that the offset array filled.
//!
//! The queries are the first protein of `tests/fasta/PemoMJNVA.faa` (the presence bit is set at
//! few subject words) and all 115 of them (72,314 residues, four query batches; the offset array
//! fills, so the scan stops and resumes); the subject is `tests/fasta/MelaMJNV.fasta`. Every run
//! sets `LOSAT_X_TBNSCANBYTE=0`: that switch (on by default, `tests/tblastn_scan_byte.rs`) scans
//! tblastn's table in place of the loop tested here.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The two query files: the first protein of `PemoMJNVA.faa` (written once per test process),
/// and the whole file.
fn queries() -> &'static [PathBuf; 2] {
    static FILES: std::sync::OnceLock<[PathBuf; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| {
        let all = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PemoMJNVA.faa");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_scan_bf");
        std::fs::create_dir_all(&dir).expect("work directory");
        let one = dir.join(format!("one.{}.faa", std::process::id()));
        let text = std::fs::read_to_string(&all).expect("query file");
        let second = text[1..].find('>').map_or(text.len(), |at| at + 1);
        std::fs::write(&one, &text[..second]).expect("query file");
        [one, all]
    })
}

/// One TBLASTN run against `MelaMJNV.fasta` (`-outfmt 6` on standard output) with `variables`,
/// `LOSAT_X_STATS` and `LOSAT_X_TBNSCANBYTE=0` as its only `LOSAT_*` variables, except
/// `LOSAT_BLASTX_WORKER_LOG`: a build
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
        .env("LOSAT_X_TBNSCANBYTE", "0")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// The scan units and calls of the one `[X_STATS] TBNSCANBF` line: `[units, calls]`; `None`
/// without the line (no unit of the unit path was scanned).
fn scans(stderr: &str) -> Option<[u64; 2]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNSCANBF "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(["units", "calls"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (above). The report does not
// depend on the threads, and the line says the loop scanned at one and at four threads. The search
// of all proteins took more calls than units: the offset array filled.
#[test]
fn the_loop_scans_and_the_report_is_the_same_at_every_thread_count() {
    for (query, many) in queries().iter().zip([false, true]) {
        let (expected, stderr) = tblastn(query, 1, &[]);
        assert!(
            !expected.is_empty(),
            "{query:?}: the search must report hits"
        );
        let [units, calls] = scans(&stderr).expect("a line");
        assert!(units > 0, "{query:?}");
        if many {
            assert!(calls > units, "{query:?}: the offset array never filled");
        }

        let (report, stderr) = tblastn(query, 4, &[]);
        assert!(report == expected, "{query:?}: four threads");
        let [units, _] = scans(&stderr).expect("a line at four threads");
        assert!(units > 0, "{query:?}: four threads");
    }
}

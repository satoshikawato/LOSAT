//! tblastn's lookup scan without a branch on the presence vector (`LOSAT_X_TBNSCANBF`, in the
//! strict set; `src/algorithm/tblastn/x_scan_bf.rs`) on whole searches.
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
//! The switch's loop returns NCBI's offset pairs and scan range on every call (the unit test
//! `x_scan_bf_equals_ncbi_scan_on_random_subjects` compares the calls), so a search with it gives
//! the report of NCBI's loop (`LOSAT_X_TBNSCANBF=0`) and takes as many scan calls. The switch is
//! on by default, off with `LOSAT_X_TBNSCANBF=0`, and off with `LOSAT_X_DEFAULT=0` unless its
//! variable is set; it acts in the unit path (`LOSAT_X_TBNPAR`). These tests compare the reports
//! and read the `[X_STATS] TBNSCANBF units_bf= calls_bf= units_ref= calls_ref=` line
//! (`LOSAT_X_STATS=1`) to check which loop scanned.
//!
//! The queries are the first protein of `tests/fasta/PemoMJNVA.faa` (the presence bit is set at
//! few subject words) and all 115 of them (72,314 residues, four query batches; the offset array
//! fills, so the scan stops and resumes); the subject is `tests/fasta/MelaMJNV.fasta`.

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

/// The scan units and calls of the one `[X_STATS] TBNSCANBF` line, by the switch's loop and by
/// NCBI's: `[units_bf, calls_bf, units_ref, calls_ref]`; `None` without the line (no unit of the
/// unit path was scanned).
fn scans(stderr: &str) -> Option<[u64; 4]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNSCANBF "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(
        ["units_bf", "calls_bf", "units_ref", "calls_ref"].map(|field| {
            line.split_whitespace()
                .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
                .unwrap_or_else(|| panic!("no {field} in {line}"))
                .parse()
                .expect("a number")
        }),
    )
}

const OFF: (&str, &str) = ("LOSAT_X_TBNSCANBF", "0");
const ON: (&str, &str) = ("LOSAT_X_TBNSCANBF", "1");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const UNIT_PATH: (&str, &str) = ("LOSAT_X_TBNPAR", "1");

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (above). Each search gives the
// report of NCBI's loop with the switch on (the default) and off, at one and four threads; the
// line says which loop scanned, and at one thread both loops took the same number of scan calls.
// The search of all proteins took more calls than units: the offset array filled.
#[test]
fn the_report_and_the_scan_calls_are_those_of_ncbis_loop() {
    for (query, many) in queries().iter().zip([false, true]) {
        let (expected, stderr) = tblastn(query, 1, &[OFF]);
        assert!(
            !expected.is_empty(),
            "{query:?}: the search must report hits"
        );
        let [units_bf, _, units_ref, calls_ref] = scans(&stderr).expect("a line with =0");
        assert_eq!(units_bf, 0, "{query:?}: =0 scans with NCBI's loop");
        assert!(units_ref > 0, "{query:?}");
        if many {
            assert!(
                calls_ref > units_ref,
                "{query:?}: the offset array never filled"
            );
        }

        let (report, stderr) = tblastn(query, 1, &[]);
        assert!(
            report == expected,
            "{query:?}: the default changes the report"
        );
        let [units_bf, calls_bf, units, _] = scans(&stderr).expect("a line by default");
        assert_eq!(
            units, 0,
            "{query:?}: the default scans with the switch's loop"
        );
        assert_eq!([units_bf, calls_bf], [units_ref, calls_ref], "{query:?}");

        for variables in [&[][..], &[OFF]] {
            let (report, stderr) = tblastn(query, 4, variables);
            assert!(report == expected, "{query:?} {variables:?}: four threads");
            let [units_bf, _, units_ref, _] = scans(&stderr).expect("a line at four threads");
            assert_eq!(
                [units_bf > 0, units_ref > 0],
                [variables.is_empty(), !variables.is_empty()],
                "{query:?} {variables:?}"
            );
        }
    }
}

// No NCBI counterpart: where the switch acts. With `LOSAT_X_DEFAULT=0` the reference loop searches
// each subject, no unit of the unit path is scanned (no line), and the variable changes nothing.
// With `LOSAT_X_TBNPAR` set as well the unit path runs again, and the switch is off unless its
// variable is set, as before it joined the set.
#[test]
fn the_switch_acts_in_the_unit_path_only() {
    let [one, _] = queries();
    let (expected, _) = tblastn(one, 1, &[]);
    for variables in [&[PREVIOUS][..], &[PREVIOUS, ON]] {
        let (report, stderr) = tblastn(one, 1, variables);
        assert!(report == expected, "{variables:?}");
        assert_eq!(scans(&stderr), None, "{variables:?}");
    }
    let (report, stderr) = tblastn(one, 1, &[PREVIOUS, UNIT_PATH]);
    assert!(report == expected);
    let [units_bf, _, units_ref, _] = scans(&stderr).expect("a line in the unit path");
    assert_eq!(units_bf, 0);
    assert!(units_ref > 0);
    let (report, stderr) = tblastn(one, 1, &[PREVIOUS, UNIT_PATH, ON]);
    assert!(report == expected);
    let [units_bf, _, units_ref, _] = scans(&stderr).expect("a line in the unit path");
    assert!(units_bf > 0);
    assert_eq!(units_ref, 0);
}

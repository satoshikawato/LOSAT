//! tblastn's lookup scan on presence bytes (`LOSAT_X_TBNSCANBYTE`, in the strict set;
//! `src/algorithm/tblastn/x_scan_byte.rs`) on whole searches.
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
//! The switch's scan returns NCBI's offset pairs and scan range on every call (the unit test
//! `x_scan_byte_equals_ncbi_scan_on_random_subjects` compares the calls of both forms), so a search
//! with it gives the report of the scan before it (`LOSAT_X_TBNSCANBYTE=0`) and takes as many scan
//! calls. The switch is on by default: the AVX2 forms where the CPU has AVX2, the plain loops for
//! sparse tables (at most three tenths of the cells with hits) on x86-64 without AVX2 and in
//! WebAssembly, as before on aarch64; `=2` takes the plain loops on any CPU for every table and
//! `=3` for sparse tables only. It is off with `=0`, and off with `LOSAT_X_DEFAULT=0` unless its
//! variable is set;
//! it acts in the unit path, which every search takes. These tests compare the reports and read the
//! `[X_STATS] TBNSCANBYTE units= calls= avx2=` line (`LOSAT_X_STATS=1`), the
//! `[X_STATS] TBNSCANBF` line of the scans before, and the `[X_SHADOW] TBNSCANBYTE` line of
//! `LOSAT_X_TBNSCANBYTESHADOW`.
//!
//! The queries are the first protein of `tests/fasta/PemoMJNVA.faa` (few candidate words, a sparse
//! table) and all 115 of them (four query batches, no sparse table; the offset array fills, so the
//! scan stops and resumes); the
//! subject is `tests/fasta/MelaMJNV.fasta`.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The two query files: the first protein of `PemoMJNVA.faa` (written once per test process),
/// and the whole file.
fn queries() -> &'static [PathBuf; 2] {
    static FILES: std::sync::OnceLock<[PathBuf; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| {
        let all = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PemoMJNVA.faa");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_scan_byte");
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
        command.env(
            WORKER_LOG,
            query.with_file_name(format!("worker.{}.log", std::process::id())),
        );
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

/// The numbers of `fields` in the one line of `stderr` that starts with `prefix`; `None` without
/// the line.
fn line_fields<const N: usize>(stderr: &str, prefix: &str, fields: [&str; N]) -> Option<[u64; N]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with(prefix))
        .collect();
    assert!(lines.len() <= 1, "one {prefix} line per process:\n{stderr}");
    let line = lines.first()?;
    Some(fields.map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

/// `[units, calls]` of the `[X_STATS] TBNSCANBYTE` line, and whether it says the AVX2 forms ran.
fn byte_scans(stderr: &str) -> Option<([u64; 2], bool)> {
    let counts = line_fields(stderr, "[X_STATS] TBNSCANBYTE ", ["units", "calls"])?;
    let avx2 = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] TBNSCANBYTE "))
        .and_then(|line| {
            line.split_whitespace()
                .find_map(|w| w.strip_prefix("avx2="))
        })
        .expect("avx2= in the line");
    Some((counts, avx2 == "true"))
}

/// Scan units and calls of the `[X_STATS] TBNSCANBF` line, both loops together.
fn earlier_scans(stderr: &str) -> Option<[u64; 2]> {
    let [units_bf, calls_bf, units_ref, calls_ref] = line_fields(
        stderr,
        "[X_STATS] TBNSCANBF ",
        ["units_bf", "calls_bf", "units_ref", "calls_ref"],
    )?;
    Some([units_bf + units_ref, calls_bf + calls_ref])
}

/// Whether the default (and `=1`) takes the plain loops for sparse tables: no AVX2, on the targets
/// the rule was timed on.
fn sparse_default() -> bool {
    !avx2() && cfg!(any(target_arch = "x86_64", target_arch = "wasm32"))
}

fn avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

const OFF: (&str, &str) = ("LOSAT_X_TBNSCANBYTE", "0");
const ON: (&str, &str) = ("LOSAT_X_TBNSCANBYTE", "1");
const PLAIN: (&str, &str) = ("LOSAT_X_TBNSCANBYTE", "2");
const SPARSE: (&str, &str) = ("LOSAT_X_TBNSCANBYTE", "3");
const SHADOW: (&str, &str) = ("LOSAT_X_TBNSCANBYTESHADOW", "1");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (above). Each search gives the
// report of the scan before the switch with the switch on (the default; =2, the plain loops) and in the
// shadow, at one and four threads; the lines say which scan ran, and at one thread every scan took
// the same number of calls. The search of all proteins took more calls than units: the offset
// array filled.
#[test]
fn the_report_and_the_scan_calls_are_those_of_the_scan_before() {
    for (query, many) in queries().iter().zip([false, true]) {
        let (expected, stderr) = tblastn(query, 1, &[OFF]);
        assert!(
            !expected.is_empty(),
            "{query:?}: the search must report hits"
        );
        assert_eq!(byte_scans(&stderr), None, "{query:?}: =0 scans as before");
        let [units, calls] = earlier_scans(&stderr).expect("a TBNSCANBF line with =0");
        assert!(units > 0, "{query:?}");
        if many {
            assert!(calls > units, "{query:?}: the offset array never filled");
        }

        // The default scans this way where the CPU has AVX2 (with the AVX2 forms), with the
        // plain loops for the sparse table of one protein without AVX2 (`=3` on any CPU), as before
        // elsewhere; `=2` takes the plain loops on any CPU.
        let (report, stderr) = tblastn(query, 1, &[]);
        assert!(
            report == expected,
            "{query:?}: the default changed the report"
        );
        let plain_sparse = |stderr: &str, what: &str| {
            if many {
                assert_eq!(
                    byte_scans(stderr),
                    None,
                    "{query:?} {what}: no sparse table"
                );
                assert_eq!(
                    earlier_scans(stderr),
                    Some([units, calls]),
                    "{query:?} {what}"
                );
            } else {
                let (counts, ran_avx2) = byte_scans(stderr).expect("a TBNSCANBYTE line");
                assert_eq!(counts, [units, calls], "{query:?} {what}");
                assert!(!ran_avx2, "{query:?} {what}: the plain loops");
                assert_eq!(earlier_scans(stderr), None, "{query:?} {what}");
            }
        };
        if avx2() {
            let (counts, ran_avx2) = byte_scans(&stderr).expect("a TBNSCANBYTE line");
            assert_eq!(counts, [units, calls], "{query:?}");
            assert!(ran_avx2, "{query:?}: the default takes the AVX2 forms");
            assert_eq!(earlier_scans(&stderr), None, "{query:?}");
        } else if sparse_default() {
            plain_sparse(&stderr, "default without AVX2");
        } else {
            assert_eq!(
                byte_scans(&stderr),
                None,
                "{query:?}: aarch64, the scan before"
            );
            assert_eq!(earlier_scans(&stderr), Some([units, calls]), "{query:?}");
        }
        let (report, stderr) = tblastn(query, 1, &[SPARSE]);
        assert!(report == expected, "{query:?}: =3 changed the report");
        plain_sparse(&stderr, "=3");
        let (report, stderr) = tblastn(query, 1, &[PLAIN]);
        assert!(report == expected, "{query:?}: =2 changed the report");
        let (counts, ran_avx2) = byte_scans(&stderr).expect("a TBNSCANBYTE line with =2");
        assert_eq!(counts, [units, calls], "{query:?}: =2");
        assert!(!ran_avx2, "{query:?}: =2 takes the plain loops");
        assert_eq!(earlier_scans(&stderr), None, "{query:?}: =2");

        let (report, stderr) = tblastn(query, 1, &[SHADOW]);
        assert!(
            report == expected,
            "{query:?}: the shadow changed the report"
        );
        let [compared, pairs] = line_fields(
            &stderr,
            "[X_SHADOW] TBNSCANBYTE ",
            ["calls_compared", "pairs_compared"],
        )
        .expect("a shadow line");
        assert_eq!(compared, calls, "{query:?}: one comparison per scan call");
        assert!(pairs > 0, "{query:?}: no pair compared");

        for variables in [&[][..], &[OFF]] {
            let (report, stderr) = tblastn(query, 4, variables);
            assert!(report == expected, "{query:?} {variables:?}: four threads");
            let byte = variables.is_empty() && (avx2() || (sparse_default() && !many));
            assert_eq!(
                [
                    byte_scans(&stderr).is_some(),
                    earlier_scans(&stderr).is_some()
                ],
                [byte, !byte],
                "{query:?} {variables:?}"
            );
        }
    }
}

// No NCBI counterpart: where the switch acts. With `LOSAT_X_DEFAULT=0` the switch is off unless its
// variable (or the shadow's) is set, as before it joined the set; `=1` then scans this way where
// the CPU has AVX2 (and the sparse table of one protein without AVX2 off aarch64), `=2`, `=3` and
// the shadow on any CPU.
#[test]
fn the_default_off_turns_the_switch_off_unless_its_variable_is_set() {
    let [one, _] = queries();
    let (expected, _) = tblastn(one, 1, &[]);
    let (report, stderr) = tblastn(one, 1, &[PREVIOUS]);
    assert!(report == expected);
    assert_eq!(byte_scans(&stderr), None);
    assert!(earlier_scans(&stderr).is_some());
    for (variables, byte) in [
        (&[PREVIOUS, ON][..], avx2() || sparse_default()),
        (&[PREVIOUS, PLAIN], true),
        (&[PREVIOUS, SPARSE], true),
        (&[PREVIOUS, SHADOW], true),
    ] {
        let (report, stderr) = tblastn(one, 1, variables);
        assert!(report == expected, "{variables:?}");
        assert_eq!(byte_scans(&stderr).is_some(), byte, "{variables:?}");
        assert_eq!(earlier_scans(&stderr).is_some(), !byte, "{variables:?}");
    }
}

//! tblastn's frames of a packed subject read from its packed bytes, 32 codons at a time with AVX2
//! (`LOSAT_X_TBNXLATESIMD`, in the strict set; `src/algorithm/tblastn/x_subject_translate.rs`),
//! on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596
//! The vector loop and the word loop give the frames of `x_generate_frames` (the unit test
//! `x_packed_frames_read_from_the_bytes_are_the_frames_of_the_codes`), so the report does not
//! depend on the switch. It is on by default and off with `LOSAT_X_TBNXLATESIMD=0`;
//! it acts on packed subjects (`LOSAT_X_TBNSTORE`, here `-subject`) and only where the CPU has
//! AVX2. Where the AVX2 loop does not run, `LOSAT_X_TBNXLATE128` reads 16 codons at a time with
//! a 128-bit loop (wasm32 simd128, x86-64 SSSE3 and SSE4.1; `=2` takes it in place of AVX2, as
//! the last test does). These tests compare the reports and read the `[X_STATS] TBNXLATESIMD
//! codons_vector= codons_scalar=` line (`LOSAT_X_STATS=1`), which counts the codons of either
//! vector loop.
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
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_xlate_simd");
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

/// The fields of the one `[X_STATS] TBNXLATESIMD` line; `None` without the line.
fn stats(stderr: &str) -> Option<[u64; 2]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNXLATESIMD "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(["codons_vector", "codons_scalar"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

const OFF: (&str, &str) = ("LOSAT_X_TBNXLATESIMD", "0");
const NO_STORE: (&str, &str) = ("LOSAT_X_TBNSTORE", "0");

/// Whether a vector loop can run here: AVX2, or the 128-bit loop (`LOSAT_X_TBNXLATE128`, SSSE3
/// and SSE4.1).
fn vector_cpu() -> bool {
    #[cfg(target_arch = "x86_64")]
    return std::arch::is_x86_feature_detected!("avx2") || w128_cpu();
    #[cfg(not(target_arch = "x86_64"))]
    return false;
}

/// Whether the 128-bit loop can run here (`LOSAT_X_TBNXLATE128=2` takes it on any such CPU).
fn w128_cpu() -> bool {
    #[cfg(target_arch = "x86_64")]
    return std::arch::is_x86_feature_detected!("ssse3")
        && std::arch::is_x86_feature_detected!("sse4.1");
    #[cfg(not(target_arch = "x86_64"))]
    return false;
}

const W128: (&str, &str) = ("LOSAT_X_TBNXLATE128", "2");

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (above). The default and `=0`
// give the same report at one and four threads; by default the vector loop read codons where
// the CPU has AVX2, with `=0` none.
#[test]
fn the_report_does_not_depend_on_the_loop() {
    for query in queries() {
        let (expected, stderr) = tblastn(query, 1, &[OFF]);
        assert!(
            !expected.is_empty(),
            "{query:?}: the search must report hits"
        );
        let [vector, scalar] = stats(&stderr).expect("a line with =0");
        assert!(vector == 0 && scalar > 0, "{query:?}");
        let (report, stderr) = tblastn(query, 1, &[]);
        assert!(
            report == expected,
            "{query:?}: the default changes the report"
        );
        let [vector, _] = stats(&stderr).expect("a line by default");
        assert_eq!(vector > 0, vector_cpu(), "{query:?}");
        for variables in [&[][..], &[OFF]] {
            let (report, _) = tblastn(query, 4, variables);
            assert!(report == expected, "{query:?} {variables:?}: four threads");
        }
    }
}

// No NCBI counterpart: where the switch acts. With `LOSAT_X_TBNSTORE=0` the subjects are not held
// packed: no packed subject, no line.
#[test]
fn the_switch_acts_on_packed_subjects_only() {
    let [one, _] = queries();
    let (expected, _) = tblastn(one, 1, &[]);
    let (report, stderr) = tblastn(one, 1, &[NO_STORE]);
    assert!(report == expected);
    assert_eq!(stats(&stderr), None);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (above). `LOSAT_X_TBNXLATE128`'s
// 128-bit loop (the loop of the wasm32 build, here its x86-64 SSSE3 / SSE4.1 build, taken with `=2`
// in place of AVX2) gives the report of the word loop (`LOSAT_X_TBNXLATESIMD=0`) at one and four
// threads, and reads codons where the CPU has SSSE3 and SSE4.1.
#[test]
fn the_128_bit_loop_gives_the_same_report() {
    for query in queries() {
        let (expected, _) = tblastn(query, 1, &[OFF]);
        let (report, stderr) = tblastn(query, 1, &[W128]);
        assert!(
            report == expected,
            "{query:?}: the 128-bit loop changes the report"
        );
        let [vector, _] = stats(&stderr).expect("a line with =2");
        assert_eq!(vector > 0, w128_cpu(), "{query:?}");
        let (report, _) = tblastn(query, 4, &[W128]);
        assert!(report == expected, "{query:?}: four threads");
    }
}

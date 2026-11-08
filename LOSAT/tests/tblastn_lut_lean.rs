//! tblastn's lookup tables built without chains (`LOSAT_X_TBNLUTLEAN`, in the strict set;
//! `src/algorithm/tblastx/lookup/x_lut_direct.rs` `x_build_lean`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
//! ```c
//!     for (i = 0; i < lookup->backbone_size; i++) {
//!         if (exact_backbone[i] != NULL) {
//!             s_AddWordHits(lookup, matrix, query->sequence,
//!                           exact_backbone[i], query_bias, row_max);
//!             sfree(exact_backbone[i]);
//!         }
//!     }
//! ```
//! The switch builds each query batch's table with the count, layout and fill passes of the other
//! protein programs' in-place builder, every temporary freed by the build; the table is the chain
//! builder's (the unit test `lean_table_equals_reference_for_thresholds_masks_matrices_and_word_lengths`
//! compares them cell by cell), so a search with it gives the report of the chain builder
//! (`LOSAT_X_TBNLUTLEAN=0`). These tests compare the reports at the two word-finder settings
//! LOSAT's TBLASTN takes (BLOSUM62 word size 3 threshold 13; BLOSUM45 word size 2 threshold 16) and
//! read the `[X_STATS] TBNLUTLEAN lean= chain=` line (`LOSAT_X_STATS=1`): every table lean by
//! default, every table with chains with `=0`.
//!
//! The queries are the first protein of `tests/fasta/PemoMJNVA.faa` and all 115 of them (four
//! query batches); the subject is `tests/fasta/MelaMJNV.fasta`.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The two query files: the first protein of `PemoMJNVA.faa` (written once per test process),
/// and the whole file.
fn queries() -> &'static [PathBuf; 2] {
    static FILES: std::sync::OnceLock<[PathBuf; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| {
        let all = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PemoMJNVA.faa");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_lut_lean");
        std::fs::create_dir_all(&dir).expect("work directory");
        let one = dir.join(format!("one.{}.faa", std::process::id()));
        let text = std::fs::read_to_string(&all).expect("query file");
        let second = text[1..].find('>').map_or(text.len(), |at| at + 1);
        std::fs::write(&one, &text[..second]).expect("query file");
        [one, all]
    })
}

/// One TBLASTN run against `MelaMJNV.fasta` (`-outfmt 6` on standard output) with `options` and
/// `variables` and `LOSAT_X_STATS` as its only `LOSAT_*` variables, except
/// `LOSAT_BLASTX_WORKER_LOG`: a build with the worker probe (`--all-features`) writes to that file
/// as soon as a pool starts, so the caller's value is kept, or a file of the tests is named.
/// Returns the report and the standard error.
fn tblastn(
    query: &Path,
    options: &[&str],
    threads: usize,
    variables: &[(&str, &str)],
) -> (Vec<u8>, String) {
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
        .args(["-outfmt", "6"])
        .args(options)
        .args(["-num_threads", &threads.to_string()])
        .env("LOSAT_X_STATS", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// `[lean, chain]` of the one `[X_STATS] TBNLUTLEAN` line of `stderr`.
fn builds(stderr: &str) -> [u64; 2] {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNLUTLEAN "))
        .collect();
    assert_eq!(lines.len(), 1, "one TBNLUTLEAN line per process:\n{stderr}");
    ["lean", "chain"].map(|field| {
        lines[0]
            .split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {}", lines[0]))
            .parse()
            .expect("a number")
    })
}

/// The default's report is the chain builder's (`=0`) at both settings, at one and four threads,
/// and the default builds every table lean, `=0` every table with chains.
#[test]
fn the_lean_tables_give_the_chain_builders_report() {
    let [one, all] = queries();
    let blosum45: &[&str] = &[
        "-matrix",
        "BLOSUM45",
        "-gapopen",
        "14",
        "-gapextend",
        "2",
        "-word_size",
        "2",
        "-threshold",
        "16",
        "-window_size",
        "60",
        "-comp_based_stats",
        "0",
    ];
    for (query, options) in [(all, &[][..]), (one, &[][..]), (all, blosum45)] {
        for threads in [1, 4] {
            let (lean, lean_err) = tblastn(query, options, threads, &[]);
            let (chain, chain_err) =
                tblastn(query, options, threads, &[("LOSAT_X_TBNLUTLEAN", "0")]);
            assert!(!lean.is_empty(), "{options:?}: the search finds hits");
            assert_eq!(
                lean, chain,
                "{options:?}, {threads} threads: the reports differ"
            );
            let [built, chains] = builds(&lean_err);
            assert!(built > 0 && chains == 0, "default: {lean_err}");
            assert_eq!(
                builds(&chain_err),
                [0, built],
                "=0 builds as many tables, with chains: {chain_err}"
            );
        }
    }
}

//! The megablast greedy rows in two passes per block (`LOSAT_X_GREEDYROWFUSE`, in the strict set;
//! `src/algorithm/blastn/alignment/x_greedy_row_fuse.rs`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537,571-573
//! ```c
//!         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
//! ...
//!             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
//! ```
//! The switch chooses between two builds of one distance of `BLAST_GreedyAlign` on x86-64 CPUs
//! with AVX2, BMI1, BMI2 and LZCNT: the three passes of `x_greedy_row_avx2` (`=0`) or the two passes
//! of the fused build (the unit tests compare the two on single distances and whole alignments).
//! Other CPUs and targets run the plain three passes either way. So the report must not depend on
//! it, and with `LOSAT_X_STATS` the `[X_STATS] GREEDYROWFUSE kernel= fused_rows=` line says which
//! build ran and how many distances the fused build ran. Every traceback runs the rows; the
//! preliminary stage runs its own fused loop (GREEDYFUSE).
//!
//! Searches: NZ_CP006932 against itself (the frozen NCBI 2.17.0 report of `megablast_xdrop.rs`),
//! and Sakai against MG1655 (`-outfmt 0`, the alignments letter by letter).

use std::path::{Path, PathBuf};
use std::process::Command;

const OFF: (&str, &str) = ("LOSAT_X_GREEDYROWFUSE", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const LINE: &str = "[X_STATS] GREEDYROWFUSE ";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One megablast search with `variables` as its only `LOSAT_*` variables, except
/// `LOSAT_BLASTX_WORKER_LOG`: a build with the worker probe (`--all-features`, as in CI) writes
/// to that file as soon as a pool starts, so the caller's value is kept, or a file of the temporary
/// directory is named. Returns the report and the standard error.
fn megablast(
    query: &str,
    subject: &str,
    outfmt: &str,
    threads: u32,
    variables: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-greedyrowfuse-worker.log"),
        );
    }
    let output = command
        .arg("blastn")
        .args(["-task", "megablast"])
        .arg("-query")
        .arg(fasta(query))
        .arg("-subject")
        .arg(fasta(subject))
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{query} x {subject} {variables:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// `kernel=` and `fused_rows=` of the one `[X_STATS] GREEDYROWFUSE` line; `None` without it.
fn row_line(stderr: &str) -> Option<(String, u64)> {
    let lines: Vec<_> = stderr.lines().filter(|l| l.starts_with(LINE)).collect();
    assert!(lines.len() <= 1, "one {LINE}line per process:\n{stderr}");
    let line = lines.first()?;
    let value = |field: &str| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .to_string()
    };
    Some((
        value("kernel"),
        value("fused_rows").parse().expect("a number"),
    ))
}

/// The build the default runs: `Some("fused")` on x86-64 with AVX2, BMI1, BMI2 and LZCNT,
/// `Some("plain")` on other x86-64 CPUs, `None` (no line) on other targets.
fn default_kernel() -> Option<&'static str> {
    #[cfg(target_arch = "x86_64")]
    {
        let avx2 = std::is_x86_feature_detected!("avx2")
            && std::is_x86_feature_detected!("bmi1")
            && std::is_x86_feature_detected!("bmi2")
            && std::is_x86_feature_detected!("lzcnt");
        Some(if avx2 { "fused" } else { "plain" })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        None
    }
}

/// The line of a run with `LOSAT_X_STATS`: the default's build with its distances counted, and
/// with `=0` the three passes (`avx2`, or `plain` without AVX2) and no fused distance.
fn check_line(stderr: &str, off: bool, what: &str) {
    let got = row_line(stderr);
    match default_kernel() {
        None => assert_eq!(got, None, "{what}: {stderr}"),
        Some(kernel) => {
            let (got_kernel, fused_rows) =
                got.unwrap_or_else(|| panic!("{what}: no line\n{stderr}"));
            if kernel == "fused" && !off {
                assert_eq!(got_kernel, "fused", "{what}");
                assert!(fused_rows > 0, "{what}: no fused row ran");
            } else {
                let three = if kernel == "fused" { "avx2" } else { "plain" };
                assert_eq!((&*got_kernel, fused_rows), (three, 0), "{what}");
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ```
/// NZ_CP006932 against itself: the default's report (the fused build, and with `LOSAT_X_STATS` the
/// build that counts its distances) and `=0`'s are the frozen NCBI report, and the line says which
/// build ran.
#[test]
fn megablast_with_fused_rows_gives_the_ncbi_report() {
    let ncbi: &[u8] = include_bytes!("fixtures/megablast_automatic_xdrop.outfmt6");
    let (report, stderr) = megablast("NZ_CP006932.fasta", "NZ_CP006932.fasta", "6", 1, &[]);
    assert!(
        report == ncbi,
        "NZ_CP006932: the report differs from NCBI's"
    );
    assert_eq!(stderr, "", "NZ_CP006932: no line without LOSAT_X_STATS");
    for off in [false, true] {
        let mut variables = vec![STATS];
        if off {
            variables.push(OFF);
        }
        let (report, stderr) =
            megablast("NZ_CP006932.fasta", "NZ_CP006932.fasta", "6", 1, &variables);
        let what = format!("NZ_CP006932 off={off}");
        assert!(report == ncbi, "{what}: the report differs from NCBI's");
        check_line(&stderr, off, &what);
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ```
/// Sakai against MG1655 with the alignments (`-outfmt 0`): the default's report (the fused build at
/// one thread, the counting build at four) is `=0`'s, and no `[X_STATS]` line without
/// `LOSAT_X_STATS`.
#[test]
fn megablast_with_fused_rows_gives_the_report_without_them() {
    const OUTFMT: &str = "0";
    let (reference, stderr) = megablast("Sakai.fna", "MG1655.fna", OUTFMT, 1, &[OFF]);
    assert!(!reference.is_empty(), "no hits");
    assert_eq!(stderr, "");
    let (report, stderr) = megablast("Sakai.fna", "MG1655.fna", OUTFMT, 1, &[]);
    assert!(
        report == reference,
        "Sakai x MG1655: the report differs from =0's"
    );
    assert_eq!(stderr, "", "Sakai x MG1655");
    let (report, stderr) = megablast("Sakai.fna", "MG1655.fna", OUTFMT, 4, &[STATS]);
    assert!(
        report == reference,
        "Sakai x MG1655 threads 4: the report differs from =0's"
    );
    check_line(&stderr, false, "Sakai x MG1655 threads 4");
}

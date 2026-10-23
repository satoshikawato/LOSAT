//! BLASTX's two-hit seed stage with the subject scan and the diagonal test fused, and the hits of
//! one lookup cell tested in batches of eight (`LOSAT_X_BXDIAG`, `seed/x_seed_diag.rs`), on whole
//! searches.
//!
//! The switch is in the strict set (on unless `LOSAT_X_BXDIAG=0`; not set, it acts only with the
//! AVX2 cell test), and the reference runs here set it to `0`; its value `per-hit` or `scalar`
//! forces a cell test, and `LOSAT_X_BXDIAGSHADOW` runs the reference beside it for
//! every subject and stops at a difference. `LOSAT_X_BXDIAGMIN` sends the query batches with few
//! hits per lookup cell to the other loop (`0`: none). `LOSAT_X_BXDEFER` (strict set) makes the
//! AVX2 cell test the deferred calls with a masked last batch on the tables with at least
//! `LOSAT_X_BXDEFERMIN` hits per used cell (the others count in `tables_avx2=`); `0` gives the
//! AVX2 batch. `LOSAT_X_BXCROSS` (an opt-in prototype, off unless set to a value other than `0`)
//! makes the AVX2 cell test the batches cut across the cells' borders on every fused table, run at
//! the staged hits `LOSAT_X_BXCROSSCHUNK` gives (an experiment knob: a multiple of eight from 64 to
//! 4096, else 1024); the timing line then ends with `cross_chunk=` and `cross_splits=` (the
//! batches cut short before a hit that repeats an earlier one's entry). So the
//! properties under test are that the report, standard error without the diagnostic lines and the exit status equal those of the
//! search with the switch off, with the default and every cell test, with `-subject` and `-db`, outfmt 6 and 7,
//! one and three threads and a two-hit window other than the default, and with every batch, some
//! or none fused; that the shadow compared subjects; that the timing line names the cell test that
//! ran; and that it counts the batches each loop took. That the bytes equal NCBI's is the
//! business of the BLASTX fixtures.
//!
//! The queries are 12 pieces of AP027131 of 1,500 to 3,900 letters each; the subjects are the
//! first 150 proteins of AP027131, so the pieces have many hits and the lookup cells many
//! entries.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The switch off: the reference seed loop (the switch is in the strict set).
const OFF: (&str, &str) = ("LOSAT_X_BXDIAG", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_diag_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The records of a FASTA file as (title, letters).
fn records(name: &str) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(title) = line.strip_prefix('>') {
            records.push((title.to_string(), Vec::new()));
        } else if let Some(record) = records.last_mut() {
            record.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    records
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) {
    let mut text = Vec::new();
    for (title, letters) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in letters.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(path, text).expect("write FASTA");
}

/// The inputs of the searches, removed when dropped.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
    db: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let genome = &records("AP027131.fasta")[0].1;
        let queries: Vec<(String, Vec<u8>)> = (0..12)
            .map(|index| {
                let start = index * 40_000;
                let length = 1_500 + (index * 373) % 2_400;
                (
                    format!("q{index} piece {start}"),
                    genome[start..start + length].to_vec(),
                )
            })
            .collect();
        let mut proteins = records("AP027131.faa");
        proteins.truncate(150);
        let query = dir.join("query.fna");
        let subject = dir.join("subject.faa");
        write_fasta(&query, &queries);
        write_fasta(&subject, &proteins);
        let db = dir.join("db");
        let made = losat(&[])
            .arg("makedb")
            .arg("-in")
            .arg(&subject)
            .args(["-dbtype", "prot", "-out"])
            .arg(&db)
            .output()
            .expect("run losat makedb");
        assert!(
            made.status.success(),
            "makedb: {}",
            String::from_utf8_lossy(&made.stderr)
        );
        Inputs {
            dir,
            query,
            subject,
            db,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The binary with `env` as its only variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// names a file of its own: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one thread,
/// and it selects nothing.
fn losat(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
}

impl Run {
    /// Standard error without the diagnostic lines (`[TIMING]`, `[X_SHADOW]`).
    fn messages(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with("[TIMING]") && !line.starts_with("[X_SHADOW]"))
            .map(|line| format!("{line}\n"))
            .collect()
    }
}

fn blastx(args: &[&str], env: &[(&str, &str)]) -> Run {
    let output = losat(env)
        .arg("blastx")
        .args(args)
        .output()
        .expect("run LOSAT");
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

/// The searches of the test: (`-subject` or `-db`, outfmt, threads, window).
fn searches(inputs: &Inputs) -> Vec<Vec<String>> {
    let mut searches = Vec::new();
    for target in [
        ["-subject", inputs.subject.to_str().unwrap()],
        ["-db", inputs.db.to_str().unwrap()],
    ] {
        for (outfmt, threads, window) in [("6", "1", "40"), ("7", "3", "40"), ("6", "1", "200")] {
            searches.push(
                [
                    "-query",
                    inputs.query.to_str().unwrap(),
                    target[0],
                    target[1],
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    threads,
                    "-window_size",
                    window,
                ]
                .map(String::from)
                .to_vec(),
            );
        }
    }
    searches
}

#[test]
fn the_fused_seed_loop_gives_the_reports_of_the_reference() {
    let inputs = Inputs::new();
    for args in searches(&inputs) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let reference = blastx(&args, &[OFF]);
        assert_eq!(reference.code, Some(0), "{args:?}: {}", reference.stderr);
        assert!(
            reference.stdout.split(|&b| b == b'\n').count() > 50,
            "{args:?}: too few rows"
        );
        // LOSAT_X_BXDIAGMIN=0 makes every batch take the fused loop; the built-in value and one
        // above every table's hits per cell (no batch fused) give the same report. An empty value
        // is the default (the switch on, with the AVX2 cell test only). LOSAT_X_BXDEFER=0 gives
        // the AVX2 batch without the deferred calls and the masked last batch, and
        // LOSAT_X_BXDEFERMIN=0 takes LOSAT_X_BXDEFER's cell test for every fused table.
        for (value, min, defer, defer_min) in [
            ("", "", "", ""),
            ("", "0", "", ""),
            ("", "0", "", "0"),
            ("", "0", "0", ""),
            ("1", "0", "", ""),
            ("per-hit", "0", "", ""),
            ("scalar", "0", "", ""),
            ("1", "", "", ""),
            ("1", "1000000", "", ""),
        ] {
            let run = blastx(
                &args,
                &[
                    ("LOSAT_X_BXDIAG", value),
                    ("LOSAT_X_BXDIAGMIN", min),
                    ("LOSAT_X_BXDEFER", defer),
                    ("LOSAT_X_BXDEFERMIN", defer_min),
                ],
            );
            let what = format!(
                "{args:?} LOSAT_X_BXDIAG={value} LOSAT_X_BXDIAGMIN={min} LOSAT_X_BXDEFER={defer} \
                 LOSAT_X_BXDEFERMIN={defer_min}"
            );
            assert_eq!(run.code, reference.code, "{what}: {}", run.stderr);
            assert!(run.stdout == reference.stdout, "{what}: report differs");
            assert_eq!(run.messages(), reference.messages(), "{what}");
        }
    }
}

#[test]
fn the_shadow_compares_every_subject_and_the_timing_line_names_the_cell_test() {
    let inputs = Inputs::new();
    let args = searches(&inputs).remove(0);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[OFF]);
    // LOSAT_X_BXDEFERMIN=0: LOSAT_X_BXDEFER's cell test on every table the shadow compares.
    for (value, defer) in [("1", ""), ("1", "0"), ("per-hit", ""), ("scalar", "")] {
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXDIAG", value),
                ("LOSAT_X_BXDEFER", defer),
                ("LOSAT_X_BXDEFERMIN", "0"),
                ("LOSAT_X_BXDIAGSHADOW", "1"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(
            run.code,
            Some(0),
            "{value} LOSAT_X_BXDEFER={defer}: {}",
            run.stderr
        );
        assert!(
            run.stdout == reference.stdout,
            "{value} LOSAT_X_BXDEFER={defer}: report differs"
        );
        // Every query batch scans every subject; the shadow compares each scan.
        let compared = run
            .stderr
            .lines()
            .find_map(|line| line.strip_prefix("[X_SHADOW] BXDIAG subjects_compared="))
            .and_then(|rest| rest.split(' ').next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| {
                panic!(
                    "{value} LOSAT_X_BXDEFER={defer}: no shadow line: {}",
                    run.stderr
                )
            });
        assert!(
            compared >= 150,
            "{value} LOSAT_X_BXDEFER={defer}: {compared} subjects compared"
        );
        let line = run
            .stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] blastx_diag: "))
            .unwrap_or_else(|| {
                panic!(
                    "{value} LOSAT_X_BXDEFER={defer}: no timing line: {}",
                    run.stderr
                )
            });
        let mode = line
            .split_whitespace()
            .find_map(|field| field.strip_prefix("mode="))
            .unwrap();
        match value {
            "per-hit" => assert_eq!(mode, "per-hit", "{line}"),
            "scalar" => assert_eq!(mode, "scalar", "{line}"),
            _ => {
                #[cfg(target_arch = "x86_64")]
                let expected = if !std::arch::is_x86_feature_detected!("avx2") {
                    "per-hit"
                } else if defer == "0" {
                    "avx2"
                } else {
                    "avx2-defer-masked"
                };
                #[cfg(not(target_arch = "x86_64"))]
                let expected = "per-hit";
                assert_eq!(mode, expected, "{line}");
            }
        }
        if mode != "per-hit" {
            assert!(!line.contains(" cells=0 "), "{line}");
        }
    }
}

/// The number after `name=` in a timing line.
fn field(line: &str, name: &str) -> u64 {
    line.split_whitespace()
        .find_map(|f| f.strip_prefix(name)?.strip_prefix('='))
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| panic!("no {name}= in {line}"))
}

#[test]
fn the_batches_with_few_hits_per_cell_take_the_other_loop() {
    let inputs = Inputs::new();
    let args = searches(&inputs).remove(0);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[OFF]);
    let timing = |min: &str| {
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXDIAG", "1"),
                ("LOSAT_X_BXDIAGMIN", min),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{min}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{min}: report differs");
        run.stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] blastx_diag: "))
            .unwrap_or_else(|| panic!("{min}: no timing line: {}", run.stderr))
            .to_string()
    };
    // Every table fused: the fused loop scanned every subject.
    let all = timing("0");
    assert_eq!(field(&all, "min_hits_per_cell"), 0, "{all}");
    assert!(field(&all, "tables_fused") > 0, "{all}");
    assert_eq!(field(&all, "tables_other"), 0, "{all}");
    assert!(field(&all, "subjects") >= 150, "{all}");
    let range = all
        .split_whitespace()
        .find_map(|f| f.strip_prefix("hits_per_cell="))
        .unwrap();
    let (least, most) = range.split_once('-').unwrap();
    let (least, most): (u64, u64) = (least.parse().unwrap(), most.parse().unwrap());
    assert!(0 < least && least <= most, "{all}");
    // EXPERIMENT (LOSAT_X_BXDEFERMIN): every fused table takes LOSAT_X_BXDEFER's cell test with 0
    // and the AVX2 batch with a minimum above every table's (the per-hit loop without AVX2).
    #[cfg(target_arch = "x86_64")]
    let avx2 = std::arch::is_x86_feature_detected!("avx2");
    #[cfg(not(target_arch = "x86_64"))]
    let avx2 = false;
    for (defer_min, unmasked) in [("0".to_string(), 0), ((most + 1).to_string(), 1)] {
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXDIAG", "1"),
                ("LOSAT_X_BXDIAGMIN", "0"),
                ("LOSAT_X_BXDEFERMIN", &defer_min),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{defer_min}: {}", run.stderr);
        assert!(
            run.stdout == reference.stdout,
            "{defer_min}: report differs"
        );
        let line = run
            .stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] blastx_diag: "))
            .unwrap_or_else(|| panic!("{defer_min}: no timing line: {}", run.stderr));
        let fused = field(line, "tables_fused");
        assert!(fused > 0, "{line}");
        let expected = if avx2 { unmasked * fused } else { 0 };
        assert_eq!(field(line, "tables_avx2"), expected, "{line}");
    }
    // The most as the minimum: the tables below it take the AVX2 batch and the others
    // LOSAT_X_BXDEFER's cell test, in one search.
    assert!(
        least < most,
        "the batches need different hits per cell: {all}"
    );
    let run = blastx(
        &args,
        &[
            ("LOSAT_X_BXDIAG", "1"),
            ("LOSAT_X_BXDIAGMIN", "0"),
            ("LOSAT_X_BXDEFERMIN", &most.to_string()),
            ("LOSAT_TIMING", "1"),
        ],
    );
    assert_eq!(run.code, Some(0), "both: {}", run.stderr);
    assert!(run.stdout == reference.stdout, "both: report differs");
    let line = run
        .stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] blastx_diag: "))
        .unwrap_or_else(|| panic!("both: no timing line: {}", run.stderr));
    let unmasked = field(line, "tables_avx2");
    if avx2 {
        assert!(
            0 < unmasked && unmasked < field(line, "tables_fused"),
            "{line}"
        );
    } else {
        assert_eq!(unmasked, 0, "{line}");
    }
    // A minimum above every table's: no subject went through the fused loop.
    let none = timing(&(most + 1).to_string());
    assert_eq!(field(&none, "tables_fused"), 0, "{none}");
    assert!(field(&none, "tables_other") > 0, "{none}");
    assert_eq!(field(&none, "subjects"), 0, "{none}");
    // The default (the switch not set) acts only with the AVX2 cell test.
    let run = blastx(&args, &[("LOSAT_X_BXDIAGMIN", "0"), ("LOSAT_TIMING", "1")]);
    assert_eq!(run.code, Some(0), "default: {}", run.stderr);
    assert!(run.stdout == reference.stdout, "default: report differs");
    let line = run
        .stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] blastx_diag: "));
    match line {
        Some(line) => {
            assert!(avx2, "default without AVX2: {line}");
            assert!(line.contains("mode=avx2-defer-masked "), "{line}");
            assert!(field(line, "subjects") >= 150, "{line}");
        }
        None => assert!(!avx2, "default with AVX2: no timing line: {}", run.stderr),
    }
    // LOSAT_X_BXDEFER=0 (LOSAT_X_BXDIAG not set), and LOSAT_X_DEFAULT=0 with LOSAT_X_BXDIAG set:
    // the fused loop with the AVX2 batch where the CPU has AVX2.
    for (name, value, named) in [
        ("LOSAT_X_BXDEFER", "0", false),
        ("LOSAT_X_DEFAULT", "0", true),
    ] {
        let what = format!("{name}={value}");
        let diag = if named { "1" } else { "" };
        let run = blastx(
            &args,
            &[
                (name, value),
                ("LOSAT_X_BXDIAG", diag),
                ("LOSAT_X_BXDIAGMIN", "0"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        let line = run
            .stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] blastx_diag: "));
        match line {
            Some(line) => {
                assert!(avx2 || named, "{what} without AVX2: {line}");
                let mode = if avx2 { "mode=avx2 " } else { "mode=per-hit " };
                assert!(line.contains(mode), "{what}: {line}");
                assert!(field(line, "subjects") >= 150, "{what}: {line}");
            }
            None => assert!(!avx2 && !named, "{what}: no timing line: {}", run.stderr),
        }
    }
}

/// The AVX2 cell test this machine runs with LOSAT_X_BXCROSS (the per-hit loop without AVX2; it
/// also needs POPCNT, else LOSAT_X_BXDEFER's cell test runs).
fn cross_mode() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        return if std::arch::is_x86_feature_detected!("popcnt") {
            "avx2-cross"
        } else {
            "avx2-defer-masked"
        };
    }
    "per-hit"
}

/// The `[TIMING] blastx_diag` line of a run, if it printed one.
fn diag_line(run: &Run) -> Option<&str> {
    run.stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] blastx_diag: "))
}

// EXPERIMENT (LOSAT_X_BXCROSS): the batches across the cells' borders give the reports of the
// reference with every search, with the built-in LOSAT_X_BXDIAGMIN and with every table fused, and
// LOSAT_X_BXDEFERMIN does not reach them (`tables_avx2=0`); the timing line names `avx2-cross`.
// With every table fused, the batches are cut short before repeated entries (`cross_splits=`: the
// hits of a diagonal in consecutive cells, the two-hit pairs of these related sequences), with the
// built-in chunk and with LOSAT_X_BXCROSSCHUNK=64 (a drain every few cells).
#[test]
fn the_cross_cell_batches_give_the_reports_of_the_reference() {
    let inputs = Inputs::new();
    let mode = cross_mode();
    for args in searches(&inputs) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let reference = blastx(&args, &[OFF]);
        assert_eq!(reference.code, Some(0), "{args:?}: {}", reference.stderr);
        for (min, defer_min, chunk) in [
            ("", "", ""),
            ("0", "", ""),
            ("0", "1000000", ""),
            ("0", "", "64"),
        ] {
            let run = blastx(
                &args,
                &[
                    ("LOSAT_X_BXCROSS", "1"),
                    ("LOSAT_X_BXCROSSCHUNK", chunk),
                    ("LOSAT_X_BXDIAGMIN", min),
                    ("LOSAT_X_BXDEFERMIN", defer_min),
                    ("LOSAT_TIMING", "1"),
                ],
            );
            let what = format!(
                "{args:?} LOSAT_X_BXCROSS=1 LOSAT_X_BXCROSSCHUNK={chunk} LOSAT_X_BXDIAGMIN={min} \
                 LOSAT_X_BXDEFERMIN={defer_min}"
            );
            assert_eq!(run.code, reference.code, "{what}: {}", run.stderr);
            assert!(run.stdout == reference.stdout, "{what}: report differs");
            assert_eq!(run.messages(), reference.messages(), "{what}");
            // LOSAT_X_BXDIAG not set: the fused loop only with an AVX2 cell test.
            match diag_line(&run) {
                Some(line) if mode == "avx2-cross" => {
                    assert!(line.contains("mode=avx2-cross "), "{what}: {line}");
                    assert_eq!(field(line, "tables_avx2"), 0, "{what}: {line}");
                    let expected = if chunk.is_empty() { 1024 } else { 64 };
                    assert_eq!(field(line, "cross_chunk"), expected, "{what}: {line}");
                    if min == "0" {
                        assert!(field(line, "subjects") >= 150, "{what}: {line}");
                        assert!(field(line, "batched_hits") > field(line, "cells"), "{line}");
                        assert!(field(line, "candidates") > 0, "{what}: {line}");
                        assert!(field(line, "cross_splits") > 0, "{what}: {line}");
                    }
                }
                Some(line) => assert_eq!(mode, "avx2-defer-masked", "{what}: {line}"),
                None => assert_eq!(mode, "per-hit", "{what}: no timing line: {}", run.stderr),
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_BXCROSSCHUNK): the knob takes a multiple of eight from 64 to 4096 (white
// space around it ignored); any other value, or none, gives the built-in 1024. The reports equal
// the reference's with each.
#[test]
fn the_cross_chunk_takes_the_multiples_of_eight_from_64_to_4096() {
    let inputs = Inputs::new();
    let args = searches(&inputs).remove(0);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[OFF]);
    for (value, chunk) in [
        ("64", 64),
        ("72", 72),
        ("4096", 4096),
        (" 128 ", 128),
        ("65", 1024),
        ("56", 1024),
        ("4104", 1024),
        ("x", 1024),
        ("", 1024),
    ] {
        let what = format!("LOSAT_X_BXCROSSCHUNK={value:?}");
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXCROSS", "1"),
                ("LOSAT_X_BXCROSSCHUNK", value),
                ("LOSAT_X_BXDIAGMIN", "0"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        let line = diag_line(&run);
        if cross_mode() == "avx2-cross" {
            let line = line.unwrap_or_else(|| panic!("{what}: no timing line: {}", run.stderr));
            assert_eq!(field(line, "cross_chunk"), chunk, "{what}: {line}");
        } else if let Some(line) = line {
            assert!(!line.contains("cross_chunk="), "{what}: {line}");
        }
    }
}

// EXPERIMENT (LOSAT_X_BXCROSS): the shadow compares the batches across the cells' borders with the
// reference for every subject, and `0` or an empty value leaves LOSAT_X_BXDEFER's cell test.
#[test]
fn the_shadow_compares_the_cross_cell_batches_and_zero_turns_them_off() {
    let inputs = Inputs::new();
    let args = searches(&inputs).remove(0);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[OFF]);
    // The built-in chunk and the smallest.
    for chunk in ["", "64"] {
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXCROSS", "1"),
                ("LOSAT_X_BXCROSSCHUNK", chunk),
                ("LOSAT_X_BXDIAG", "1"),
                ("LOSAT_X_BXDIAGSHADOW", "1"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "shadow {chunk}: {}", run.stderr);
        assert!(
            run.stdout == reference.stdout,
            "shadow {chunk}: report differs"
        );
        let compared = run
            .stderr
            .lines()
            .find_map(|line| line.strip_prefix("[X_SHADOW] BXDIAG subjects_compared="))
            .and_then(|rest| rest.split(' ').next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("shadow {chunk}: no shadow line: {}", run.stderr));
        assert!(
            compared >= 150,
            "shadow {chunk}: {compared} subjects compared"
        );
        let line = diag_line(&run)
            .unwrap_or_else(|| panic!("shadow {chunk}: no timing line: {}", run.stderr));
        assert!(line.contains(&format!("mode={} ", cross_mode())), "{line}");
    }
    // Off: LOSAT_X_BXDEFER's cell test (on every table with LOSAT_X_BXDEFERMIN=0).
    #[cfg(target_arch = "x86_64")]
    let off = if std::arch::is_x86_feature_detected!("avx2") {
        "avx2-defer-masked"
    } else {
        "per-hit"
    };
    #[cfg(not(target_arch = "x86_64"))]
    let off = "per-hit";
    for value in ["0", " 0 ", ""] {
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXCROSS", value),
                ("LOSAT_X_BXDIAG", "1"),
                ("LOSAT_X_BXDIAGMIN", "0"),
                ("LOSAT_X_BXDEFERMIN", "0"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(
            run.code,
            Some(0),
            "LOSAT_X_BXCROSS={value:?}: {}",
            run.stderr
        );
        assert!(
            run.stdout == reference.stdout,
            "LOSAT_X_BXCROSS={value:?}: report differs"
        );
        let line = diag_line(&run)
            .unwrap_or_else(|| panic!("LOSAT_X_BXCROSS={value:?}: no timing line: {}", run.stderr));
        assert!(
            line.contains(&format!("mode={off} ")),
            "LOSAT_X_BXCROSS={value:?}: {line}"
        );
    }
}

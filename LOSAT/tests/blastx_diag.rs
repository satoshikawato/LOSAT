//! BLASTX's two-hit seed stage with the subject scan and the diagonal test fused, and the hits of
//! one lookup cell tested in batches of eight (`LOSAT_X_BXDIAG`, `seed/x_seed_diag.rs`), on whole
//! searches.
//!
//! The switch is in the strict set (on unless `LOSAT_X_BXDIAG=0`; not set, it acts only with
//! `LOSAT_X_BXCROSS`'s cell test), and the reference runs here set it to `0`; its value `per-hit`
//! or `scalar` forces a cell test, and `LOSAT_X_BXDIAGSHADOW` runs the reference beside it for
//! every subject and stops at a difference. `LOSAT_X_BXDIAGMIN` sends the query batches with few
//! hits per lookup cell to the other loop (`0`: none; built-in 670 for every cell test).
//! `LOSAT_X_BXCROSS` (strict set) makes the AVX2 cell test, where the CPU also has POPCNT, the
//! batches cut across the cells' borders on every fused table; the timing line then ends with
//! `cross_splits=` (the windows of eight staged hits with a repeated entry, tested in parts); `0`
//! leaves the direct loop on every table when `LOSAT_X_BXDIAG` is not set, and the masked AVX2
//! batch (`mode=avx2-defer-masked`) when it has a value.
//! `LOSAT_X_BXCROSSDEFER` (strict set) defers those batches' candidate calls (`mode=avx2-cross-defer`);
//! `0` calls them at once (`mode=avx2-cross`).
//! So the
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
            .arg("makeblastdb")
            .arg("-in")
            .arg(&subject)
            .args(["-dbtype", "prot", "-out"])
            .arg(&db)
            .output()
            .expect("run losat makeblastdb");
        assert!(
            made.status.success(),
            "makeblastdb: {}",
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
        // is the default (the switch on, with LOSAT_X_BXCROSS's cell test only, where the CPU has
        // POPCNT). LOSAT_X_BXCROSS=0 gives the masked AVX2 batch with a value (with the default it
        // is the direct loop, the reference's own path: the other test checks that no fused subject
        // runs).
        for (value, min, cross) in [
            ("", "", ""),
            ("", "0", ""),
            ("1", "0", "0"),
            ("1", "0", ""),
            ("per-hit", "0", ""),
            ("scalar", "0", ""),
            ("1", "", ""),
            ("1", "1000000", ""),
        ] {
            let run = blastx(
                &args,
                &[
                    ("LOSAT_X_BXDIAG", value),
                    ("LOSAT_X_BXDIAGMIN", min),
                    ("LOSAT_X_BXCROSS", cross),
                ],
            );
            let what = format!(
                "{args:?} LOSAT_X_BXDIAG={value} LOSAT_X_BXDIAGMIN={min} LOSAT_X_BXCROSS={cross}"
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
    // The shadow compares every table: with LOSAT_X_BXCROSS=0 the masked AVX2 batch (the default
    // takes LOSAT_X_BXCROSS's cell test).
    // (an empty LOSAT_X_BXDIAG is not set: the shadow alone, with LOSAT_X_BXCROSS=0, also runs the
    // masked batch).
    for (value, cross) in [
        ("1", ""),
        ("1", "0"),
        ("", "0"),
        ("per-hit", ""),
        ("scalar", ""),
    ] {
        let what = format!("{value} LOSAT_X_BXCROSS={cross}");
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXDIAG", value),
                ("LOSAT_X_BXCROSS", cross),
                ("LOSAT_X_BXDIAGSHADOW", "1"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        // Every query batch scans every subject; the shadow compares each scan.
        let compared = run
            .stderr
            .lines()
            .find_map(|line| line.strip_prefix("[X_SHADOW] BXDIAG subjects_compared="))
            .and_then(|rest| rest.split(' ').next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("{what}: no shadow line: {}", run.stderr));
        assert!(compared >= 150, "{what}: {compared} subjects compared");
        let line =
            diag_line(&run).unwrap_or_else(|| panic!("{what}: no timing line: {}", run.stderr));
        let mode = line
            .split_whitespace()
            .find_map(|field| field.strip_prefix("mode="))
            .unwrap();
        let expected = match value {
            "per-hit" => "per-hit",
            "scalar" => "scalar",
            _ if cross != "0" => cross_mode(""),
            _ => masked_mode(),
        };
        assert_eq!(mode, expected, "{what}: {line}");
        if mode != "per-hit" {
            assert!(!line.contains(" cells=0 "), "{what}: {line}");
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
    let avx2 = masked_mode() != "per-hit";
    // The most hits per cell as the minimum (the cutoff is `>=`): the tables with fewer take the
    // direct loop and the others the fused one, in one search.
    assert!(
        least < most,
        "the tables need different hits per cell: {all}"
    );
    let both = timing(&most.to_string());
    assert_eq!(field(&both, "min_hits_per_cell"), most, "{both}");
    assert!(field(&both, "tables_fused") > 0, "{both}");
    assert!(field(&both, "tables_other") > 0, "{both}");
    assert!(field(&both, "subjects") > 0, "{both}");
    // A minimum above every table's: no subject went through the fused loop.
    let none = timing(&(most + 1).to_string());
    assert_eq!(field(&none, "tables_fused"), 0, "{none}");
    assert!(field(&none, "tables_other") > 0, "{none}");
    assert_eq!(field(&none, "subjects"), 0, "{none}");
    // The default (the switch not set) acts only with LOSAT_X_BXCROSS's cell test, on every fused
    // table: it needs AVX2 and POPCNT.
    let cross = cross_mode("").starts_with("avx2-cross");
    let run = blastx(&args, &[("LOSAT_X_BXDIAGMIN", "0"), ("LOSAT_TIMING", "1")]);
    assert_eq!(run.code, Some(0), "default: {}", run.stderr);
    assert!(run.stdout == reference.stdout, "default: report differs");
    match diag_line(&run) {
        Some(line) => {
            assert!(cross, "default without the cross cell test: {line}");
            assert!(
                line.contains(&format!("mode={} ", cross_mode(""))),
                "{line}"
            );
            assert!(field(line, "subjects") >= 150, "{line}");
            assert!(field(line, "cross_splits") > 0, "{line}");
        }
        None => assert!(
            !cross,
            "default with the cross cell test: no timing line: {}",
            run.stderr
        ),
    }
    // LOSAT_X_BXCROSS=0, with LOSAT_X_BXDIAG not set and set: not set, no subject goes through the
    // fused loop (the timing line is not printed, as the switch is off) and the report is the
    // reference's; set, the fused loop takes the masked AVX2 batch where the CPU has AVX2 (the
    // per-hit test elsewhere).
    for named in [false, true] {
        let diag = if named { "1" } else { "" };
        let what = format!("LOSAT_X_BXCROSS=0 LOSAT_X_BXDIAG={diag:?}");
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXCROSS", "0"),
                ("LOSAT_X_BXDIAG", diag),
                ("LOSAT_X_BXDIAGMIN", "0"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        match (named, diag_line(&run)) {
            (true, Some(line)) => {
                let mode = if avx2 {
                    "mode=avx2-defer-masked "
                } else {
                    "mode=per-hit "
                };
                assert!(line.contains(mode), "{what}: {line}");
                assert!(field(line, "subjects") >= 150, "{what}: {line}");
            }
            (true, None) => panic!("{what}: no timing line: {}", run.stderr),
            (false, line) => assert!(line.is_none(), "{what}: a fused subject: {line:?}"),
        }
    }
}

/// The AVX2 cell test this machine runs with LOSAT_X_BXCROSS (the default; the per-hit loop
/// without AVX2; it also needs POPCNT, else the masked AVX2 batch runs where `LOSAT_X_BXDIAG` has a
/// value) and LOSAT_X_BXCROSSDEFER=`defer` (`0`: the calls at once).
fn cross_mode(defer: &str) -> &'static str {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        return if !std::arch::is_x86_feature_detected!("popcnt") {
            "avx2-defer-masked"
        } else if defer == "0" {
            "avx2-cross"
        } else {
            "avx2-cross-defer"
        };
    }
    let _ = defer;
    "per-hit"
}

/// The cell test this machine runs with LOSAT_X_BXCROSS=0 and a value of `LOSAT_X_BXDIAG`: the
/// masked AVX2 batch, the per-hit loop without AVX2.
fn masked_mode() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        return "avx2-defer-masked";
    }
    "per-hit"
}

/// The `[TIMING] blastx_diag` line of a run, if it printed one.
fn diag_line(run: &Run) -> Option<&str> {
    run.stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] blastx_diag: "))
}

// EXPERIMENT (LOSAT_X_BXCROSS): the batches across the cells' borders (the default) give the
// reports of the reference with every search, with the built-in LOSAT_X_BXDIAGMIN (670) and with
// every table fused; the timing line names `avx2-cross-defer` (LOSAT_X_BXCROSSDEFER, the default)
// or, with `LOSAT_X_BXCROSSDEFER=0`, `avx2-cross`. With every table fused, the batches are cut
// short before repeated entries (`cross_splits=`: the hits of a diagonal in consecutive cells, the
// two-hit pairs of these related sequences).
#[test]
fn the_cross_cell_batches_give_the_reports_of_the_reference() {
    let inputs = Inputs::new();
    for args in searches(&inputs) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let reference = blastx(&args, &[OFF]);
        assert_eq!(reference.code, Some(0), "{args:?}: {}", reference.stderr);
        for (cross, min, cross_defer) in [
            ("", "", ""),
            ("", "0", ""),
            ("1", "0", ""),
            ("", "0", "0"),
            ("", "", "0"),
        ] {
            let run = blastx(
                &args,
                &[
                    ("LOSAT_X_BXCROSS", cross),
                    ("LOSAT_X_BXDIAGMIN", min),
                    ("LOSAT_X_BXCROSSDEFER", cross_defer),
                    ("LOSAT_TIMING", "1"),
                ],
            );
            let what = format!(
                "{args:?} LOSAT_X_BXCROSS={cross} LOSAT_X_BXDIAGMIN={min} \
                 LOSAT_X_BXCROSSDEFER={cross_defer}"
            );
            assert_eq!(run.code, reference.code, "{what}: {}", run.stderr);
            assert!(run.stdout == reference.stdout, "{what}: report differs");
            assert_eq!(run.messages(), reference.messages(), "{what}");
            // LOSAT_X_BXDIAG not set: the fused loop only with LOSAT_X_BXCROSS's cell test.
            let expected = cross_mode(cross_defer);
            match diag_line(&run) {
                Some(line) => {
                    assert!(expected.starts_with("avx2-cross"), "{what}: {line}");
                    assert!(
                        line.contains(&format!("mode={expected} ")),
                        "{what}: {line}"
                    );
                    if min == "0" {
                        assert!(field(line, "subjects") >= 150, "{what}: {line}");
                        assert!(field(line, "batched_hits") > field(line, "cells"), "{line}");
                        assert!(field(line, "candidates") > 0, "{what}: {line}");
                        assert!(field(line, "cross_splits") > 0, "{what}: {line}");
                    }
                    // The built-in LOSAT_X_BXDIAGMIN.
                    if min.is_empty() {
                        assert_eq!(field(line, "min_hits_per_cell"), 670, "{what}: {line}");
                    }
                }
                None => assert!(
                    !expected.starts_with("avx2-cross"),
                    "{what}: no timing line: {}",
                    run.stderr
                ),
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_BXCROSS): the shadow compares the batches across the cells' borders with the
// reference for every subject, and `0` (white space around it ignored) leaves the masked AVX2
// batch with a value of `LOSAT_X_BXDIAG`; an empty value is the default.
#[test]
fn the_shadow_compares_the_cross_cell_batches_and_zero_turns_them_off() {
    let inputs = Inputs::new();
    let args = searches(&inputs).remove(0);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[OFF]);
    for value in ["", "1"] {
        let what = format!("shadow LOSAT_X_BXCROSS={value:?}");
        let run = blastx(
            &args,
            &[
                ("LOSAT_X_BXCROSS", value),
                ("LOSAT_X_BXDIAG", "1"),
                ("LOSAT_X_BXDIAGSHADOW", "1"),
                ("LOSAT_TIMING", "1"),
            ],
        );
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        let compared = run
            .stderr
            .lines()
            .find_map(|line| line.strip_prefix("[X_SHADOW] BXDIAG subjects_compared="))
            .and_then(|rest| rest.split(' ').next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("{what}: no shadow line: {}", run.stderr));
        assert!(compared >= 150, "{what}: {compared} subjects compared");
        let line =
            diag_line(&run).unwrap_or_else(|| panic!("{what}: no timing line: {}", run.stderr));
        assert!(
            line.contains(&format!("mode={} ", cross_mode(""))),
            "{line}"
        );
    }
    // Off: the masked AVX2 batch (on every table with LOSAT_X_BXDIAGMIN=0).
    for env in [
        &[("LOSAT_X_BXCROSS", "0")][..],
        &[("LOSAT_X_BXCROSS", " 0 ")][..],
    ] {
        let what = format!("{env:?}");
        let mut env = env.to_vec();
        env.extend([
            ("LOSAT_X_BXDIAG", "1"),
            ("LOSAT_X_BXDIAGMIN", "0"),
            ("LOSAT_TIMING", "1"),
        ]);
        let run = blastx(&args, &env);
        assert_eq!(run.code, Some(0), "{what}: {}", run.stderr);
        assert!(run.stdout == reference.stdout, "{what}: report differs");
        let line =
            diag_line(&run).unwrap_or_else(|| panic!("{what}: no timing line: {}", run.stderr));
        assert!(
            line.contains(&format!("mode={} ", masked_mode())),
            "{what}: {line}"
        );
        assert!(!line.contains("cross_splits="), "{what}: {line}");
    }
}

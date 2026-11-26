//! dc-megablast's strips stage with `LOSAT_X_DCLATEST` (`src/algorithm/blastn/blast_engine/x_dc_latest.rs`,
//! `x_pair_seed_disc.rs` `seed_stage_latest`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes where the two-hit fold keeps its diagonals' inserts (a map of each diagonal's latest
//! insert instead of NCBI's diagonal hash, per partition of `LOSAT_X_PAIRPARFOLD` or the owner's
//! without partitions), not the decisions any seed makes, so the property under test is that the
//! report does not depend on it, and that the default took the switch's path (`LOSAT_X_REACH`).
//! Against `=0` at one and four threads: 9,800 letters of LC738874 against LC738868 (a query
//! block above the 8,000 letters of the diagonal hash) at the default window (40), above the word
//! length (24), and below it (12, and 4 in the band 2w + 2 < wl where NCBI's chains and LOSAT's
//! differ: below the word length no seed is extended, so no hit), and with the subject's lower case
//! masked
//! (`-lcase_masking`, the stage in the subject's ranges), and `-ungapped`.
//! On a database (`tests/db_search/dbmask/edge1`, five subjects) at one thread, gapped and
//! `-ungapped`, with a query of repeats of s5's letters (22,000 letters): unmasked, every subject
//! takes the maps; with `-db_soft_mask 11` s5's ranges keep its stored masks' order (600-699 before
//! 100-199), step back, and its chunk goes to the chains (`x_pair_seed_disc.rs`
//! `latest_ranges_ascend`, counted as `back=` on the `[TIMING] x_dc_latest:` line), the other four
//! take the maps.
//! EXPERIMENT (LOSAT_X_DCFOLDREG, `x_dc_fold_reg.rs`, also in the strict set): the fold on the maps
//! with the map in registers. The tests named `dc_fold_reg_*` compare the default with
//! `LOSAT_X_DCFOLDREG=0` on the same searches (and with `LOSAT_X_FOLDLEAN=0`, whose fold calls
//! `type_of_word` for every seed) and check `reg_calls=` on `[TIMING] x_dc_latest:`.
//! This test compares the default with `=0` only and runs no NCBI BLAST+. That the `=0` bytes
//! equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_DCLATEST";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_DCLATEST", "0");
const REG_OFF: (&str, &str) = ("LOSAT_X_DCFOLDREG", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-dclatest-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

/// The letters of the first record of `name`.
fn letters(name: &str) -> String {
    let text = std::fs::read_to_string(fasta(name)).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// Letters 20,000 to 29,800 of LC738874 as one record in `dir`.
fn lc874_9800(dir: &TestDir) -> PathBuf {
    let query = dir.path("lc874_9800.fna");
    let text = letters("LC738874.fasta");
    std::fs::write(&query, format!(">lc874_9800\n{}\n", &text[20_000..29_800]))
        .expect("write the query");
    query
}

/// LC738868 with every third piece of 2,000 letters in lower case.
fn lc868_lower(dir: &TestDir) -> PathBuf {
    let subject = dir.path("lc868_lower.fna");
    let text: String = letters("LC738868.fasta")
        .as_bytes()
        .chunks(2_000)
        .enumerate()
        .map(|(k, piece)| {
            let piece = String::from_utf8_lossy(piece).into_owned();
            if k % 3 == 1 {
                piece.to_lowercase()
            } else {
                piece
            }
        })
        .collect();
    std::fs::write(&subject, format!(">lc868_lower\n{text}\n")).expect("write the subject");
    subject
}

/// The database of `tests/db_search/dbmask` with s5's masks stored out of order (algorithm 11).
fn edge1() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/dbmask/edge1")
}

/// A query of s5's letters (`tests/db_search/dbmask/subj.fna`) in `dir`: 200 copies of letters 200
/// to 259, each followed by 20 letters of a fixed generator, then 30 copies of letters 400 to 599
/// (22,000 letters). With `-db_soft_mask 11` s5's ranges are 0-599, an inverted one and 200-1199
/// (the complement of its masks 600-699 and 100-199 in their stored order): the early copies'
/// diagonals get seeds at 200-259 in the scan of 0-599 and again when the scan steps back to 200.
fn s5_repeats(dir: &TestDir) -> PathBuf {
    let text = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/dbmask/subj.fna"),
    )
    .expect("read subj.fna");
    let s5: String = text
        .lines()
        .skip_while(|line| !line.starts_with(">s5 "))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect();
    assert_eq!(s5.len(), 1_200, "s5");
    let mut state = 0x5EED_0005_u64;
    let mut query = String::with_capacity(22_000);
    for _ in 0..200 {
        query.push_str(&s5[200..260]);
        for _ in 0..20 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            query.push(b"ACGT"[(state >> 62) as usize] as char);
        }
    }
    for _ in 0..30 {
        query.push_str(&s5[400..600]);
    }
    assert_eq!(query.len(), 22_000);
    let path = dir.path("s5_repeats.fna");
    std::fs::write(&path, format!(">s5_repeats\n{query}\n")).expect("write the query");
    path
}

/// One dc-megablast search of `query` against `target` (`("-subject", file)` or `("-db", name)`;
/// `-outfmt 6` on standard output) with `args` and the given variables as its only `LOSAT_*`
/// variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing. Returns the report and the standard error.
fn search(
    query: &Path,
    target: (&str, &Path),
    args: &[&str],
    threads: u32,
    vars: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .args(["-task", "dc-megablast"])
        .arg("-query")
        .arg(query)
        .arg(target.0)
        .arg(target.1)
        .args(args)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{args:?} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// For `query` against `subject` with `args`, at one and four threads: the report with the switch
/// off and by default, the same bytes and standard error; the default took the switch's path and
/// never the reference's (in the stage without a pool at one thread, `LOSAT_X_BNDCONE`, and in a
/// partitioned stage at four, `LOSAT_X_PAIRPARFOLD`), and `=0` records nothing for the switch.
/// `hits`: whether the report has hits.
fn check(query: &Path, subject: &Path, args: &[&str], hits: bool) {
    for threads in [1, 4] {
        let at = format!("{args:?} threads {threads}");
        let off_reach = reach_support::Reach::new("dclatest-off");
        let (expected, expected_stderr) = search(
            query,
            ("-subject", subject),
            args,
            threads,
            &[SWITCH_OFF, off_reach.var()],
        );
        assert_eq!(!expected.is_empty(), hits, "{at}: hits");
        assert_eq!(off_reach.sides(SWITCH), (false, false), "{at}: =0");
        let reach = reach_support::Reach::new("dclatest");
        let (report, stderr) = search(query, ("-subject", subject), args, threads, &[reach.var()]);
        assert!(report == expected, "{at}: default differs");
        assert_eq!(stderr, expected_stderr, "{at}");
        assert_eq!(
            reach.sides(SWITCH),
            (true, false),
            "{at}: sides of {SWITCH}"
        );
        assert!(
            reach.sides("LOSAT_X_PAIRPARDISC").0,
            "{at}: no strips stage"
        );
        if threads == 1 {
            assert!(
                reach.sides("LOSAT_X_BNDCONE").0,
                "{at}: no stage without a pool"
            );
        } else {
            assert!(
                reach.sides("LOSAT_X_PAIRPARFOLD").0,
                "{at}: no partitioned stage"
            );
        }
    }
}

// The default window and one above the word length (hits), one below it and one in the band (no
// hit).
#[test]
fn dc_megablast_windows() {
    let dir = TestDir::new("windows");
    let query = lc874_9800(&dir);
    let subject = fasta("LC738868.fasta");
    check(&query, &subject, &[], true);
    check(&query, &subject, &["-window_size", "24"], true);
    check(&query, &subject, &["-window_size", "12"], false);
    check(&query, &subject, &["-window_size", "4"], false);
}

// A subject with lower case masked: the stage scans its ranges.
#[test]
fn dc_megablast_masked_subject() {
    let dir = TestDir::new("masked");
    let query = lc874_9800(&dir);
    let subject = lc868_lower(&dir);
    check(&query, &subject, &["-lcase_masking"], true);
}

// `-ungapped`: the stage's hits are the report.
#[test]
fn dc_megablast_ungapped() {
    let dir = TestDir::new("ungapped");
    let query = lc874_9800(&dir);
    let subject = fasta("LC738868.fasta");
    check(&query, &subject, &["-ungapped"], true);
}

/// The counts `stages=`, `back=` and `reg_calls=` of the `[TIMING] x_dc_latest:` line in `stderr`.
fn latest_counts_all(stderr: &str) -> (u64, u64, u64) {
    let line = stderr
        .lines()
        .find_map(|line| line.strip_prefix("[TIMING] x_dc_latest:"))
        .expect("the x_dc_latest line");
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(name))
            .and_then(|n| n.parse().ok())
            .unwrap_or_else(|| panic!("{name} in {line}"))
    };
    (field("stages="), field("back="), field("reg_calls="))
}

/// The counts `stages=` and `back=` of the `[TIMING] x_dc_latest:` line in `stderr`.
fn latest_counts(stderr: &str) -> (u64, u64) {
    let (stages, back, _) = latest_counts_all(stderr);
    (stages, back)
}

/// `stderr` without the `[TIMING]` lines.
fn without_timing(stderr: &str) -> Vec<&str> {
    stderr
        .lines()
        .filter(|line| !line.starts_with("[TIMING]"))
        .collect()
}

/// For `query` against the database edge1 with `args`, at one thread, gapped and `-ungapped`: the
/// report with the switch off and by default (with `LOSAT_TIMING`), the same bytes and standard
/// error outside the timing lines; the default's sides of the switch, the stages on the maps and
/// the chunks sent to the chains because their ranges step back (`back`) are `sides`, `stages`
/// and `back`, all in the stage without a pool (`LOSAT_X_BNDCONE`).
fn check_db(query: &Path, args: &[&str], sides: (bool, bool), stages: u64, back: u64) {
    let db = edge1();
    for ungapped in [false, true] {
        let mut args = args.to_vec();
        if ungapped {
            args.push("-ungapped");
        }
        let at = format!("{args:?}");
        let off_reach = reach_support::Reach::new("dclatest-db-off");
        let (expected, expected_stderr) = search(
            query,
            ("-db", db.as_path()),
            &args,
            1,
            &[SWITCH_OFF, off_reach.var()],
        );
        assert!(!expected.is_empty(), "{at}: no hits");
        assert_eq!(off_reach.sides(SWITCH), (false, false), "{at}: =0");
        let reach = reach_support::Reach::new("dclatest-db");
        let (report, stderr) = search(
            query,
            ("-db", db.as_path()),
            &args,
            1,
            &[reach.var(), ("LOSAT_TIMING", "1")],
        );
        assert!(report == expected, "{at}: default differs");
        assert_eq!(
            without_timing(&stderr),
            without_timing(&expected_stderr),
            "{at}"
        );
        assert_eq!(reach.sides(SWITCH), sides, "{at}: sides of {SWITCH}");
        assert_eq!(
            latest_counts(&stderr),
            (stages, back),
            "{at}: stages and back"
        );
        assert!(
            reach.sides("LOSAT_X_BNDCONE").0,
            "{at}: no stage without a pool"
        );
    }
}

// `-db` with `-db_soft_mask 11`: s5's ranges step back and its chunk goes to the chains (in a debug
// build the maps' order assertion would stop the search there); s1 to s4 take the maps.
#[test]
fn dc_megablast_db_soft_mask_out_of_order() {
    let dir = TestDir::new("db-soft-mask");
    let query = s5_repeats(&dir);
    check_db(
        &query,
        &["-db_soft_mask", "11", "-dust", "no"],
        (true, true),
        4,
        1,
    );
}

// `-db` over five subjects with a query block above 8,000 letters: every subject takes the maps.
#[test]
fn dc_megablast_db_several_subjects() {
    let dir = TestDir::new("db");
    let query = s5_repeats(&dir);
    check_db(&query, &["-dust", "no"], (true, false), 5, 0);
}

/// `check` for `LOSAT_X_DCFOLDREG` with extra `vars` in both runs: byte-identical output and
/// non-timing stderr, `reg_calls > 0` by default and `reg_calls == 0` with `=0`.
fn check_reg(vars: &[(&str, &str)], query: &Path, subject: &Path, args: &[&str], hits: bool) {
    for threads in [1, 4] {
        let at = format!("DCFOLDREG {vars:?} {args:?} threads {threads}");
        let mut off_vars = vars.to_vec();
        off_vars.extend([REG_OFF, ("LOSAT_TIMING", "1")]);
        let (expected, expected_stderr) =
            search(query, ("-subject", subject), args, threads, &off_vars);
        assert_eq!(!expected.is_empty(), hits, "{at}: hits");
        assert_eq!(
            latest_counts_all(&expected_stderr).2,
            0,
            "{at}: =0 reg_calls"
        );
        let reach = reach_support::Reach::new("dcfoldreg");
        let mut on_vars = vars.to_vec();
        on_vars.extend([reach.var(), ("LOSAT_TIMING", "1")]);
        let (report, stderr) = search(query, ("-subject", subject), args, threads, &on_vars);
        assert!(report == expected, "{at}: default differs");
        assert_eq!(
            without_timing(&stderr),
            without_timing(&expected_stderr),
            "{at}"
        );
        assert!(latest_counts_all(&stderr).2 > 0, "{at}: reg_calls == 0");
        assert_eq!(
            reach.sides(SWITCH),
            (true, false),
            "{at}: sides of {SWITCH}"
        );
        assert!(
            reach.sides("LOSAT_X_PAIRPARDISC").0,
            "{at}: strips stage did not run"
        );
        let stage = if threads == 1 {
            "LOSAT_X_BNDCONE"
        } else {
            "LOSAT_X_PAIRPARFOLD"
        };
        assert!(reach.sides(stage).0, "{at}: {stage} did not run");
    }
}

// EXPERIMENT (LOSAT_X_DCFOLDREG): the fold with the map in registers against `=0` at the default
// window and one above the word length (hits, so extended seeds), in the band (no hit), with
// LOSAT_X_FOLDLEAN=0 (`type_of_word` for every seed), on a masked subject and `-ungapped`.
#[test]
fn dc_fold_reg_subjects() {
    let dir = TestDir::new("fold-reg");
    let query = lc874_9800(&dir);
    let subject = fasta("LC738868.fasta");
    check_reg(&[], &query, &subject, &[], true);
    check_reg(&[], &query, &subject, &["-window_size", "24"], true);
    check_reg(&[], &query, &subject, &["-window_size", "4"], false);
    check_reg(&[("LOSAT_X_FOLDLEAN", "0")], &query, &subject, &[], true);
    check_reg(&[], &query, &subject, &["-ungapped"], true);
    let lower = lc868_lower(&dir);
    check_reg(&[], &query, &lower, &["-lcase_masking"], true);
}

// EXPERIMENT (LOSAT_X_DCFOLDREG): on the database edge1 with `-db_soft_mask 11`, the four chunks
// on the maps take the fold in registers (`reg_calls > 0`) and s5's chunk, whose ranges step back,
// takes the chains (`back=1`).
#[test]
fn dc_fold_reg_db_soft_mask_out_of_order() {
    let dir = TestDir::new("fold-reg-db");
    let query = s5_repeats(&dir);
    let db = edge1();
    for ungapped in [false, true] {
        let mut args = vec!["-db_soft_mask", "11", "-dust", "no"];
        if ungapped {
            args.push("-ungapped");
        }
        let at = format!("DCFOLDREG {args:?}");
        let (expected, expected_stderr) = search(
            &query,
            ("-db", db.as_path()),
            &args,
            1,
            &[REG_OFF, ("LOSAT_TIMING", "1")],
        );
        assert!(!expected.is_empty(), "{at}: no hits");
        assert_eq!(
            latest_counts_all(&expected_stderr).2,
            0,
            "{at}: =0 reg_calls"
        );
        let reach = reach_support::Reach::new("dcfoldreg-db");
        let (report, stderr) = search(
            &query,
            ("-db", db.as_path()),
            &args,
            1,
            &[reach.var(), ("LOSAT_TIMING", "1")],
        );
        assert!(report == expected, "{at}: default differs");
        assert_eq!(
            without_timing(&stderr),
            without_timing(&expected_stderr),
            "{at}"
        );
        assert_eq!(reach.sides(SWITCH), (true, true), "{at}: sides of {SWITCH}");
        let (stages, back, reg_calls) = latest_counts_all(&stderr);
        assert_eq!((stages, back), (4, 1), "{at}: x_dc_latest counts");
        assert!(reg_calls > 0, "{at}: reg_calls == 0");
        assert!(reach.sides("LOSAT_X_BNDCONE").0, "{at}: no serial stage");
    }
}

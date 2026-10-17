//! The early release of the TBLASTX seed-stage tables (`LOSAT_X_SEEDRELEASE`,
//! `src/algorithm/tblastx/x_seed_release.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is
//! `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. With it on, the unit path of
//! `LOSAT_X_PAIRPAR` with its group tables is freed when the units of the batch's last subject
//! have run, before their lists are merged; the batch's lookup table is freed when the last
//! seed stage has ended; and the worker's diagonal table is either not created (the unit path
//! takes the batch) or dropped at that point (the reference frame loop ran). Every test reads
//! the `[X_STATS] SEEDRELEASE` line, which is printed only while the switch is on, to check
//! which of these happened, not only that the output is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
//! ```c
//! word_params = BlastInitialWordParametersFree(word_params);
//! s_BlastCoreAuxStructFree(aux_struct);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:110-113
//! ```c
//! static BlastCoreAuxStruct*
//! s_BlastCoreAuxStructFree(BlastCoreAuxStruct* aux_struct)
//! {
//!     BlastExtendWordFree(aux_struct->ewp);
//! ```
//! NCBI keeps its lookup table and its diagonal table until every subject has been searched and
//! linked. Freeing LOSAT's copies earlier must give the same report. For the fixture pair the
//! expected output (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was
//! produced by NCBI BLAST+ 2.17.0.
//!
//! The searches with several subject records are compared with the same LOSAT search with the
//! switch off (`LOSAT_X_SEEDRELEASE=0`), not with an NCBI file: the switch computes nothing, it
//! only frees memory, so the property under test is that the report does not depend on it. That
//! those reports equal NCBI's is the business of the fixture suites, which run with the switch
//! on by default and with it off under `LOSAT_X_DEFAULT=0`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

const SWITCH: (&str, &str) = ("LOSAT_X_SEEDRELEASE", "1");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_SEEDRELEASE", "0");
const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 2] = [1, 4];

/// One TBLASTX run of the fixture query against `subject` (`-outfmt 6` on standard output)
/// with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_X_DEFAULT` of the test's own environment is dropped as well: every test states the
/// default it runs with. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name
/// as soon as a search has more than one thread, and it selects nothing. `LOSAT_X_PAIRPARGATE=0`
/// comes first (a test's own variables may set it again): the fixture query is short, so the
/// gate, in the strict set, would send its batch to the reference frame loop, and these tests are
/// about the unit path. Returns the report and the standard error.
fn tblastx(subject: &Path, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.env("LOSAT_X_PAIRPARGATE", "0");
    let output = command
        .arg("tblastx")
        .arg("-query")
        .arg(fixture("short_hsp.query.fna"))
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

fn ncbi_report() -> Vec<u8> {
    std::fs::read(fixture("short_hsp.ncbi_2.17.0.fmt6.out")).expect("NCBI output")
}

fn stats_line(stderr: &str) -> Option<&str> {
    stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] SEEDRELEASE "))
}

/// The five counters of the `[X_STATS] SEEDRELEASE` line: batches, tables_released,
/// unit_tables_released, batches_without_worker_table, worker_tables_released.
fn counters(stderr: &str) -> [u64; 5] {
    let line =
        stats_line(stderr).unwrap_or_else(|| panic!("no [X_STATS] SEEDRELEASE line:\n{stderr}"));
    [
        "batches",
        "tables_released",
        "unit_tables_released",
        "batches_without_worker_table",
        "worker_tables_released",
    ]
    .map(|name| {
        let key = format!("{name}=");
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .expect("counter")
    })
}

fn fasta_sequence(path: &Path) -> String {
    std::fs::read_to_string(path)
        .expect("FASTA")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

struct TempFasta(PathBuf);

impl TempFasta {
    fn new(stem: &str, records: &[(&str, &str)]) -> Self {
        let path = std::env::temp_dir().join(format!(
            "losat_seed_release_{}_{stem}.fna",
            std::process::id()
        ));
        let text: String = records
            .iter()
            .map(|(name, sequence)| format!(">{name}\n{sequence}\n"))
            .collect();
        std::fs::write(&path, text).expect("write FASTA");
        Self(path)
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
// ```c
// word_params = BlastInitialWordParametersFree(word_params);
// s_BlastCoreAuxStructFree(aux_struct);
// ```
// The default search of the fixture pair (no variable for any switch but
// `LOSAT_X_PAIRPARGATE=0`) runs on the unit path with the switch on: one batch, its unit path released when the units of the only subject have
// run, its lookup table (a stub) by the only seed stage, and no worker ever holds a diagonal
// table. The report is NCBI's.
#[test]
fn the_unit_path_releases_the_tables_and_holds_no_worker_table() {
    let expected = ncbi_report();
    let subject = fixture("short_hsp.subject.fna");
    for threads in THREADS {
        let (report, stderr) = tblastx(&subject, threads, &[STATS]);
        assert_eq!(report, expected, "threads {threads}");
        assert_eq!(
            counters(&stderr),
            [1, 1, 1, 1, 0],
            "threads {threads}:\n{stderr}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:110-113
// ```c
// static BlastCoreAuxStruct*
// s_BlastCoreAuxStructFree(BlastCoreAuxStruct* aux_struct)
// {
//     BlastExtendWordFree(aux_struct->ewp);
// ```
// `LOSAT_X_PAIRPAR=0`, and `LOSAT_X_DEFAULT=0` with the switch set (the run that the variable
// gave before the switch joined the set), leave the search to the reference frame loop, which
// needs the worker's diagonal table (NCBI's `ewp`): the worker has it and drops it when the
// batch's lookup table is released. There is no unit path to release. The report is NCBI's.
#[test]
fn the_frame_loop_releases_the_tables_and_its_worker_table() {
    let expected = ncbi_report();
    let subject = fixture("short_hsp.subject.fna");
    for frame_loop in [&[("LOSAT_X_PAIRPAR", "0")][..], &[DEFAULT_OFF, SWITCH][..]] {
        for threads in THREADS {
            let switches = [&[STATS][..], frame_loop].concat();
            let (report, stderr) = tblastx(&subject, threads, &switches);
            assert_eq!(report, expected, "{frame_loop:?} threads {threads}");
            assert_eq!(
                counters(&stderr),
                [1, 1, 0, 0, 1],
                "{frame_loop:?} threads {threads}:\n{stderr}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
// ```c
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ```
// The shadow of the unit path runs the unit path and the frame loop for the subject, so the
// worker keeps its diagonal table until the release. The shadow of the direct group tables
// builds the full table as well and searches on the unit path alone, so no worker holds one.
// In both the unit path is released before the merge of its lists (and before the frame loop
// of the first shadow, which scans with the lookup table). Both shadows end without a
// difference and the report is NCBI's. No variable names the switch: it is on by default.
#[test]
fn the_shadows_of_the_unit_path_run_with_the_switch() {
    let expected = ncbi_report();
    let subject = fixture("short_hsp.subject.fna");
    for (shadow, expected_counters) in [
        (("LOSAT_X_PAIRPARSHADOW", "1"), [1, 1, 1, 0, 1]),
        (("LOSAT_X_PAIRLUTSHADOW", "1"), [1, 1, 1, 1, 0]),
    ] {
        for threads in THREADS {
            let (report, stderr) = tblastx(&subject, threads, &[STATS, shadow]);
            assert_eq!(report, expected, "{shadow:?} threads {threads}");
            assert_eq!(
                counters(&stderr),
                expected_counters,
                "{shadow:?} threads {threads}:\n{stderr}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
// ```c
// word_params = BlastInitialWordParametersFree(word_params);
// s_BlastCoreAuxStructFree(aux_struct);
// ```
// With `LOSAT_X_SEEDRELEASE=0`, or with `LOSAT_X_DEFAULT=0` and no variable for the switch, the
// tables live until the batch has been searched, as in NCBI: no SEEDRELEASE line is printed,
// and the report is NCBI's. `LOSAT_X_SEEDRELEASE=0` wins over the default, and
// `LOSAT_X_SEEDRELEASE=1` over `LOSAT_X_DEFAULT=0` (the frame-loop test above).
#[test]
fn the_switch_follows_its_variable_and_the_global_one() {
    let expected = ncbi_report();
    let subject = fixture("short_hsp.subject.fna");
    for off in [SWITCH_OFF, DEFAULT_OFF] {
        for threads in THREADS {
            let (report, stderr) = tblastx(&subject, threads, &[STATS, off]);
            assert_eq!(report, expected, "{off:?} threads {threads}");
            assert_eq!(
                stats_line(&stderr),
                None,
                "{off:?} threads {threads}:\n{stderr}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1429-1432
// ```c
// if (seq_arg.seq->length < min_subj_seq_length) {
//     BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
//     continue;
// }
// ```
// Several subject records, one of them shorter than a word of three codons plus two bases (not
// searched: its seed stage ends without a scan), as the last record and between the others.
// With one thread the subjects are searched in order by one worker; with four, by the pool.
// The lookup table and the unit path are each released once per search, by the stage of the
// last record that ends, and the report equals the report of the same search with the switch
// off (see the file header for why this comparison is against LOSAT).
#[test]
fn several_subjects_with_a_short_one_give_the_report_of_the_search_without_the_switch() {
    let first = fasta_sequence(&fixture("short_hsp.subject.fna"));
    let second = fasta_sequence(&fixture("short_hsp.query.fna"));
    let third: String = first.chars().rev().collect();
    let short = "ACGTACGTAC";
    let layouts = [
        (
            "short_last",
            vec![
                ("first", first.as_str()),
                ("second", second.as_str()),
                ("third", third.as_str()),
                ("short", short),
            ],
            // The last stage is the short record's, which is dropped: no worker table to drop.
            0,
        ),
        (
            "short_between",
            vec![
                ("first", first.as_str()),
                ("short", short),
                ("second", second.as_str()),
                ("third", third.as_str()),
            ],
            1,
        ),
    ];
    for (stem, records, serial_worker_tables) in &layouts {
        let subjects = TempFasta::new(stem, records);
        for path in [&[][..], &[("LOSAT_X_PAIRPAR", "0")][..]] {
            let frame_loop = !path.is_empty();
            for threads in THREADS {
                let what = format!("{stem} {path:?} threads {threads}");
                let off = [&[SWITCH_OFF, STATS][..], path].concat();
                let (expected, off_stderr) = tblastx(&subjects.0, threads, &off);
                assert!(!expected.is_empty(), "{what}: no hits");
                assert_eq!(stats_line(&off_stderr), None, "{what}:\n{off_stderr}");
                let switches = [&[STATS][..], path].concat();
                let (report, stderr) = tblastx(&subjects.0, threads, &switches);
                assert_eq!(report, expected, "{what}");
                let [batches, tables_released, unit_tables_released, without_worker_table, worker_tables_released] =
                    counters(&stderr);
                assert_eq!((batches, tables_released), (1, 1), "{what}:\n{stderr}");
                assert_eq!(
                    unit_tables_released,
                    u64::from(!frame_loop),
                    "{what}:\n{stderr}"
                );
                assert_eq!(
                    without_worker_table,
                    u64::from(!frame_loop),
                    "{what}:\n{stderr}"
                );
                if !frame_loop {
                    assert_eq!(worker_tables_released, 0, "{what}:\n{stderr}");
                } else if threads == 1 {
                    assert_eq!(
                        worker_tables_released, *serial_worker_tables,
                        "{what}:\n{stderr}"
                    );
                } else {
                    assert!(worker_tables_released <= 1, "{what}:\n{stderr}");
                }
            }
        }
    }
}

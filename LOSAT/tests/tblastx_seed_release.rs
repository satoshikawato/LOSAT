//! The early release of the TBLASTX seed-stage tables
//! (`src/algorithm/tblastx/x_seed_release.rs`) on whole searches.
//!
//! The unit path (`x_pair_par.rs`) with its group tables is freed when the units of the
//! batch's last subject have run, before their lists are merged; the batch's lookup table is
//! freed when the last seed stage has ended; and the worker's diagonal table is either not
//! created (the unit path takes the batch) or dropped at that point (the reference frame loop
//! ran). Every test reads the `[X_STATS] SEEDRELEASE` line to check which of these happened, not
//! only that the output is equal.
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
//! The searches with several subject records are compared with the same LOSAT search on the
//! reference frame loop at one thread, not with an NCBI file: the release computes nothing, it
//! only frees memory, so the property under test is that the report does not depend on the path
//! or the thread count. That those reports equal NCBI's is the business of the fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

/// The route to the reference frame loop: the cost rule of `LOSAT_X_PAIRPARGATE`, which the helper
/// switches off first. For the fixture query (696 bases, a diagonal table below 2^18 cells) it
/// sends a batch to the frame loop at one thread and, at four, a batch whose largest subject is
/// at most 1.5 threads' share of the bases (`2 * longest * threads <= 3 * total`).
const FRAME_LOOP: (&str, &str) = ("LOSAT_X_PAIRPARGATE", "1");
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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

/// The counters `units` and `fallback_cost` of the `[X_STATS] PAIRPAR` line.
fn pair_par(stderr: &str) -> [u64; 2] {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] PAIRPAR "))
        .unwrap_or_else(|| panic!("no [X_STATS] PAIRPAR line:\n{stderr}"));
    ["units", "fallback_cost"].map(|name| {
        let key = format!("{name}=");
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {line}"))
            .parse()
            .expect("counter")
    })
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
// The search of the fixture pair with `LOSAT_X_PAIRPARGATE=0` (the helper's default) runs on the
// unit path: one batch, its unit path released when the units of the only subject have run, its
// lookup table (a stub) by the only seed stage, and no worker ever holds a diagonal table. The
// report is NCBI's.
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
// The cost rule (`FRAME_LOOP`, at one thread) leaves the search to the reference frame loop, which
// needs the worker's diagonal table (NCBI's `ewp`): the worker has it and drops it when the
// batch's lookup table is released. There is no unit path to release (`units=0`). The report is
// NCBI's. At four threads the fixture subject holds all the bases, so the units run there.
#[test]
fn the_frame_loop_releases_the_tables_and_its_worker_table() {
    let expected = ncbi_report();
    let subject = fixture("short_hsp.subject.fna");
    let (report, stderr) = tblastx(&subject, 1, &[STATS, FRAME_LOOP]);
    assert_eq!(report, expected);
    assert_eq!(pair_par(&stderr), [0, 1], "{stderr}");
    assert_eq!(counters(&stderr), [1, 1, 0, 0, 1], "{stderr}");
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
// last record that ends, and every report equals the report of the same layout on the frame
// loop at one thread (see the file header for why this comparison is against LOSAT). The frame-loop
// runs take no unit (`units=0`) at either thread count: the total is 2,470 bases and the longest
// subject 882, so the cost rule holds up to four threads (2 * 882 * 4 = 7,056 <= 3 * 2,470 =
// 7,410).
#[test]
fn several_subjects_with_a_short_one_give_one_report_on_every_path() {
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
        let (expected, _) = tblastx(&subjects.0, 1, &[STATS, FRAME_LOOP]);
        assert!(!expected.is_empty(), "{stem}: no hits");
        for path in [&[][..], &[FRAME_LOOP][..]] {
            let frame_loop = !path.is_empty();
            for threads in THREADS {
                let what = format!("{stem} {path:?} threads {threads}");
                let switches = [&[STATS][..], path].concat();
                let (report, stderr) = tblastx(&subjects.0, threads, &switches);
                assert_eq!(report, expected, "{what}");
                let [units, _] = pair_par(&stderr);
                assert_eq!(units == 0, frame_loop, "{what}:\n{stderr}");
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

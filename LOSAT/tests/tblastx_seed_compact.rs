//! The compact, exactly sized seed-stage lists of the TBLASTX unit path
//! (`LOSAT_X_SEEDCOMPACT`, `src/algorithm/tblastx/x_pair_par.rs`, the section above `XUnitHsp`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_SEEDCOMPACT`
//! is `0`, or `LOSAT_X_DEFAULT=0` and it is not set. With it on, the units of
//! `LOSAT_X_PAIRPAR` keep their saved init HSPs as 28-byte records, copied out of a per-thread
//! vector at their exact length; the merge expands them into the same init lists, converts each
//! chunk's list into a list of its exact length and moves the first chunk's list into the empty
//! frame list. Every test reads the `[X_STATS] SEEDCOMPACT` line, which is printed only while the
//! switch is on, to check that the compact lists were made (or, where no unit path runs, that
//! none was), not only that the report is equal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
//! ```c
//!     ungapped_data->q_start = q_start;
//!     ungapped_data->s_start = s_start;
//!     ungapped_data->length = len;
//!     ungapped_data->score = score;
//!
//!     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
//! ```
//! NCBI saves an initial HSP as these values and its two seed offsets; the compact record is that
//! record and the context index, and the switch changes how LOSAT holds the lists, not a value.
//! For the fixture pair the expected output
//! (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was produced by NCBI
//! BLAST+ 2.17.0. The database search is compared with the same LOSAT search with the switch off
//! (`LOSAT_X_SEEDCOMPACT=0`): the property under test is that the report does not depend on it;
//! that `-db` reports equal NCBI's is the business of the database fixtures, which run on the
//! default paths.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_link_fast")
        .join(name)
}

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const SWITCH: (&str, &str) = ("LOSAT_X_SEEDCOMPACT", "1");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_SEEDCOMPACT", "0");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 2] = [1, 4];

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_X_DEFAULT` of the test's own environment is dropped as well: every test states the
/// default it runs with. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. `LOSAT_X_PAIRPARGATE=0`
/// comes first (a test's own variables may set it again): the fixtures' query batches are short,
/// so the gate, in the strict set, would send them to the reference frame loop, and these tests
/// are about the unit path.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.env("LOSAT_X_PAIRPARGATE", "0");
    command.envs(switches.iter().copied());
    command
}

/// One TBLASTX search on standard output (`-outfmt 6`); returns the report and the standard
/// error.
fn tblastx(
    query: &Path,
    subjects: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat(switches)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .args(subjects)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{subjects:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The fixture pair.
fn pair(threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    let subject = fixture("short_hsp.subject.fna");
    tblastx(
        &fixture("short_hsp.query.fna"),
        &["-subject", subject.to_str().unwrap()],
        threads,
        switches,
    )
}

fn ncbi_report() -> Vec<u8> {
    std::fs::read(fixture("short_hsp.ncbi_2.17.0.fmt6.out")).expect("NCBI output")
}

fn line<'a>(stderr: &'a str, prefix: &str) -> Option<&'a str> {
    stderr.lines().find(|line| line.starts_with(prefix))
}

/// The counters `names` of the line that starts with `prefix`.
fn fields<const N: usize>(stderr: &str, prefix: &str, names: [&str; N]) -> [u64; N] {
    let found = line(stderr, prefix).unwrap_or_else(|| panic!("no {prefix} line:\n{stderr}"));
    names.map(|name| {
        let key = format!("{name}=");
        found
            .split_whitespace()
            .find_map(|field| field.strip_prefix(key.as_str()))
            .unwrap_or_else(|| panic!("no {name} in {found}"))
            .parse()
            .expect("counter")
    })
}

/// The four counters of `[X_STATS] SEEDCOMPACT`: subjects, units, records, bytes_saved.
fn compact(stderr: &str) -> [u64; 4] {
    fields(
        stderr,
        "[X_STATS] SEEDCOMPACT ",
        ["subjects", "units", "records", "bytes_saved"],
    )
}

/// The bytes a record saves: the size of an `InitHSP` minus the 28 bytes of the record.
fn saved_per_record() -> u64 {
    (std::mem::size_of::<LOSAT::algorithm::tblastx::blast_gapalign::InitHSP>() - 28) as u64
}

/// Checks the counters of a search whose units took `subjects` subjects on the compact path, and
/// that they are those of `[X_STATS] PAIRPAR` (one chunk per frame: every init HSP of a chunk
/// list was one record).
fn assert_compact_ran(stderr: &str, subjects: u64, what: &str) {
    let [compact_subjects, units, records, bytes_saved] = compact(stderr);
    let [pairpar_subjects, pairpar_units, init_hsps] = fields(
        stderr,
        "[X_STATS] PAIRPAR ",
        ["subjects", "units", "init_hsps"],
    );
    assert_eq!(compact_subjects, subjects, "{what}:\n{stderr}");
    assert_eq!(
        (compact_subjects, units, records),
        (pairpar_subjects, pairpar_units, init_hsps),
        "{what}:\n{stderr}"
    );
    assert!(units > 0 && records > 0, "{what}:\n{stderr}");
    assert_eq!(
        bytes_saved,
        records * saved_per_record(),
        "{what}:\n{stderr}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
// ```c
//     ungapped_data->q_start = q_start;
//     ungapped_data->s_start = s_start;
//     ungapped_data->length = len;
//     ungapped_data->score = score;
//
//     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
// ```
// With no variable for the switch (it is on by default) the fixture pair runs on the compact lists
// at one and four threads (one subject, all its units and init HSPs as records) and the report is
// NCBI's.
#[test]
fn the_switch_keeps_compact_lists_and_gives_ncbis_report() {
    let expected = ncbi_report();
    for threads in THREADS {
        let what = format!("threads {threads}");
        let (report, stderr) = pair(threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        assert_compact_ran(&stderr, 1, &what);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:367-372
// ```c
//     BLAST_SaveInitialHit(ungapped_hsps, q_off, s_off, ungapped_data);
// ```
// `LOSAT_X_SEEDCOMPACT=0` and `LOSAT_X_DEFAULT=0` turn the switch off: no SEEDCOMPACT line, the
// same report. `LOSAT_X_DEFAULT=0` with `LOSAT_X_PAIRPAR=1` and `LOSAT_X_SEEDCOMPACT=1` is the run
// the variable gave before the switch joined the set: the compact lists are made.
#[test]
fn the_switch_follows_its_variable_and_the_global_one() {
    let expected = ncbi_report();
    for threads in THREADS {
        for off in [SWITCH_OFF, PREVIOUS] {
            let what = format!("{off:?} threads {threads}");
            let (report, stderr) = pair(threads, &[STATS, off]);
            assert_eq!(report, expected, "{what}");
            assert_eq!(
                line(&stderr, "[X_STATS] SEEDCOMPACT "),
                None,
                "{what}:\n{stderr}"
            );
        }
        let what = format!("previous paths with the switch set, threads {threads}");
        let (report, stderr) = pair(
            threads,
            &[STATS, PREVIOUS, ("LOSAT_X_PAIRPAR", "1"), SWITCH],
        );
        assert_eq!(report, expected, "{what}");
        assert_compact_ran(&stderr, 1, &what);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
// ```c
// if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
// ```
// The switch with the variables of the paths around it: `LOSAT_X_PAIRPARSHADOW=1` runs the unit
// path (on the compact lists) and the reference frame loop and compares the subject's combined
// list, every field; `LOSAT_X_SEEDRELEASE=0` keeps the group tables through the merge;
// `LOSAT_X_PAIRLUT=0` builds the full table and filters the groups out of it. The report is NCBI's
// in every run and the compact path ran.
#[test]
fn the_switch_with_the_shadow_and_the_neighbouring_switches() {
    let expected = ncbi_report();
    for other in [
        ("LOSAT_X_PAIRPARSHADOW", "1"),
        ("LOSAT_X_SEEDRELEASE", "0"),
        ("LOSAT_X_PAIRLUT", "0"),
    ] {
        for threads in THREADS {
            let what = format!("{other:?} threads {threads}");
            let (report, stderr) = pair(threads, &[STATS, other]);
            assert_eq!(report, expected, "{what}");
            assert_compact_ran(&stderr, 1, &what);
            if other.0 == "LOSAT_X_PAIRPARSHADOW" {
                let [subjects_compared, hsps_compared] = fields(
                    &stderr,
                    "[X_SHADOW] PAIRPAR ",
                    ["subjects_compared", "hsps_compared"],
                );
                assert_eq!(subjects_compared, 1, "{what}:\n{stderr}");
                assert!(hsps_compared > 0, "{what}:\n{stderr}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
// ```c
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ```
// `LOSAT_X_PAIRPAR=0` (the switch is on by default), and `LOSAT_X_DEFAULT=0` with the switch set,
// leave the search to the reference frame loop: no unit runs, so no compact list is made. The
// line is still printed (the switch is on) with zero counters, and the report is NCBI's.
#[test]
fn without_the_unit_path_no_compact_list_is_made() {
    let expected = ncbi_report();
    for (no_units, switches) in [
        (("LOSAT_X_PAIRPAR", "0"), &[STATS][..]),
        (PREVIOUS, &[STATS, SWITCH][..]),
    ] {
        for threads in THREADS {
            let what = format!("{no_units:?} threads {threads}");
            let mut switches = switches.to_vec();
            switches.push(no_units);
            let (report, stderr) = pair(threads, &switches);
            assert_eq!(report, expected, "{what}");
            assert_eq!(compact(&stderr), [0, 0, 0, 0], "{what}:\n{stderr}");
        }
    }
}

/// The letters of a one-record FASTA file.
fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

/// A FASTA record of `letters` in lines of 70.
fn record(defline: &str, letters: &str) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.as_bytes().chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:479-491
// ```c
// status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
//                                dbseq_chunk_overlap);
// ...
// BlastInitHitListReset(init_hitlist);
// ```
// A database search (the inputs of tests/tblastx_db_search.rs: the first 6 kb of MjeNMV against
// a database of MelaMJNV and PemoMJNVA made by `LOSAT makedb`) at one and four threads: with the
// switch (on by default) both subjects run on the compact lists, and the report is the report of
// the same search with `LOSAT_X_SEEDCOMPACT=0`.
#[test]
fn a_database_search_gives_the_report_of_the_search_without_the_switch() {
    let dir = std::env::temp_dir().join(format!(
        "losat-tblastx-seed-compact-db-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let query = dir.join("q.fna");
    std::fs::write(
        &query,
        record(
            "mje6k first 6 kb of MjeNMV",
            &letters("MjeNMV.fasta")[..6000],
        ),
    )
    .unwrap();
    let set = dir.join("set.fna");
    std::fs::write(
        &set,
        record("mela MelaMJNV", &letters("MelaMJNV.fasta"))
            + &record("pemo PemoMJNVA", &letters("PemoMJNVA.fasta")),
    )
    .unwrap();
    let db = dir.join("set").to_str().unwrap().to_string();
    let made = losat(&[])
        .arg("makedb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        made.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&made.stderr)
    );
    for threads in THREADS {
        let what = format!("-db threads {threads}");
        let (expected, off_stderr) = tblastx(&query, &["-db", &db], threads, &[STATS, SWITCH_OFF]);
        assert!(!expected.is_empty(), "{what}: no hits");
        assert_eq!(
            line(&off_stderr, "[X_STATS] SEEDCOMPACT "),
            None,
            "{what}:\n{off_stderr}"
        );
        let (report, stderr) = tblastx(&query, &["-db", &db], threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        assert_compact_ran(&stderr, 2, &what);
    }
    let _ = std::fs::remove_dir_all(&dir);
}

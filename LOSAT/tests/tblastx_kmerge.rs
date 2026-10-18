//! The merge of a subject's unit lists without full sorts (`LOSAT_X_TBXKMERGE`,
//! `src/algorithm/tblastx/x_kmerge.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_TBXKMERGE`
//! is `0`, or `LOSAT_X_DEFAULT=0` and it is not set. With it on, each unit of `LOSAT_X_PAIRPAR`
//! sorts its own list, a chunk's init list is the k-way merge of its units' lists, the conversion
//! sorts only the runs of HSPs equal in score and subject range, and each frame list is merged into
//! the subject's list. `[X_STATS] TBXKMERGE` (printed while the switch or its shadow is on) counts
//! the subjects, chunks, sorted tie runs and the two fallbacks, which must stay zero;
//! `LOSAT_X_TBXKMERGESHADOW=1` compares every field of every HSP with the merge of the switch off.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766
//! ```c
//! for (index=combined_hsp_list->hspcnt, index1=0;
//!      index1<hsp_list->hspcnt; index1++) {
//!    combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
//! }
//! combined_hsp_list->hspcnt = new_hspcnt;
//! Blast_HSPListSortByScore(combined_hsp_list);
//! ```
//! The switch changes how LOSAT orders the same lists, not a value. For the fixture pair the
//! expected output (`tests/fixtures/tblastx_link_fast/short_hsp.ncbi_2.17.0.fmt6.out`) was produced
//! by NCBI BLAST+ 2.17.0. The search with tied HSPs is compared with the same LOSAT search with the
//! switch off: the property under test is that the report does not depend on it.

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

const SWITCH: (&str, &str) = ("LOSAT_X_TBXKMERGE", "1");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_TBXKMERGE", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_TBXKMERGESHADOW", "1");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 5] = [1, 2, 4, 8, 16];

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
    subject: &Path,
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let output = losat(switches)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The fixture pair.
fn pair(threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
    tblastx(
        &fixture("short_hsp.query.fna"),
        &fixture("short_hsp.subject.fna"),
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

/// The five counters of `[X_STATS] TBXKMERGE`: subjects, chunks, tie_runs, k1_fallback,
/// frame_unsorted.
fn kmerge(stderr: &str) -> [u64; 5] {
    fields(
        stderr,
        "[X_STATS] TBXKMERGE ",
        [
            "subjects",
            "chunks",
            "tie_runs",
            "k1_fallback",
            "frame_unsorted",
        ],
    )
}

/// Checks that the switch merged `subjects` subjects, the subjects of `[X_STATS] PAIRPAR`, with
/// no fallback; returns the sorted tie runs.
fn assert_kmerge_ran(stderr: &str, subjects: u64, what: &str) -> u64 {
    let [merged, chunks, tie_runs, k1_fallback, frame_unsorted] = kmerge(stderr);
    let [pairpar_subjects] = fields(stderr, "[X_STATS] PAIRPAR ", ["subjects"]);
    assert_eq!(
        (merged, pairpar_subjects),
        (subjects, subjects),
        "{what}:\n{stderr}"
    );
    assert!(chunks >= 6 * subjects, "{what}:\n{stderr}");
    assert_eq!((k1_fallback, frame_unsorted), (0, 0), "{what}:\n{stderr}");
    tie_runs
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
// ```c
// void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
// {
//     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// With no variable for the switch (it is on by default) the fixture pair is merged by the switch
// at 1, 2, 4, 8 and 16 threads, and the report is NCBI's.
#[test]
fn the_switch_merges_the_lists_and_gives_ncbis_report() {
    let expected = ncbi_report();
    for threads in THREADS {
        let what = format!("threads {threads}");
        let (report, stderr) = pair(threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        assert_kmerge_ran(&stderr, 1, &what);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2766
// ```c
// Blast_HSPListSortByScore(combined_hsp_list);
// ```
// `LOSAT_X_TBXKMERGE=0` and `LOSAT_X_DEFAULT=0` turn the switch off: no TBXKMERGE line, the same
// report. `LOSAT_X_DEFAULT=0` with `LOSAT_X_PAIRPAR=1` and `LOSAT_X_TBXKMERGE=1` merges with the
// switch on the paths of the previous default.
#[test]
fn the_switch_follows_its_variable_and_the_global_one() {
    let expected = ncbi_report();
    for threads in [1, 4] {
        for off in [SWITCH_OFF, PREVIOUS] {
            let what = format!("{off:?} threads {threads}");
            let (report, stderr) = pair(threads, &[STATS, off]);
            assert_eq!(report, expected, "{what}");
            assert_eq!(
                line(&stderr, "[X_STATS] TBXKMERGE "),
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
        assert_kmerge_ran(&stderr, 1, &what);
    }
}

/// The `[X_SHADOW] TBXKMERGE` counters: subjects and HSPs compared.
fn shadow_counts(stderr: &str) -> [u64; 2] {
    fields(
        stderr,
        "[X_SHADOW] TBXKMERGE ",
        ["subjects_compared", "hsps_compared"],
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:844
// ```c
// if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
// ```
// The shadow on both forms of the units' lists (`LOSAT_X_SEEDCOMPACT` on by default, and off),
// and with `LOSAT_X_PAIRPARSHADOW=1`, which also compares the subject's list with the reference
// frame loop's: every HSP compared, the report NCBI's.
#[test]
fn the_shadow_finds_the_lists_of_the_switch_off() {
    let expected = ncbi_report();
    for others in [
        &[][..],
        &[("LOSAT_X_SEEDCOMPACT", "0")][..],
        &[("LOSAT_X_PAIRPARSHADOW", "1")][..],
    ] {
        for threads in [1, 4, 16] {
            let what = format!("{others:?} threads {threads}");
            let mut switches = vec![STATS, SHADOW];
            switches.extend_from_slice(others);
            let (report, stderr) = pair(threads, &switches);
            assert_eq!(report, expected, "{what}");
            let [subjects, hsps] = shadow_counts(&stderr);
            assert_eq!(subjects, 1, "{what}:\n{stderr}");
            assert!(hsps > 0, "{what}:\n{stderr}");
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1355
// ```c
// if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
//     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
//     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
// ```
// Ties that the conversion's sort reorders: query 1 is 900 bases of MjeNMV followed by the 600
// bases X (bases 2,000-2,600), query 2 is X alone, and the subject is the first 6 kb of MjeNMV
// (short enough for the unit path's exactness guard: the queries' 4.2 k residues leave a slack of
// 4 k cells in the 8,192-cell diagonal table, more than a context's 500 residues). Every
// subject hit of X is found in both queries with the same score and subject range; in the
// concatenated query (the init order) query 1 comes first, in the context (the HSP order) query
// 2's offset is 300 residues (900 bases) lower. The switch sorts those runs (`tie_runs` > 0) and the
// report, at 1, 2, 4, 8 and 16 threads, without the statistics (the units' scan without a probe)
// and with the shadow, is the report of the search with the switch off. With
// `LOSAT_TBLASTX_TEST_CHUNK_SIZE=700` (debug builds) each 2,000-residue subject frame is searched as
// three chunks, merged by `merge_tblastx_subject_chunk_hits`: the report is that of the switch off
// with the same chunks.
#[test]
fn tied_hsps_of_two_queries_give_the_report_of_the_switch_off() {
    let dir = std::env::temp_dir().join(format!("losat-tblastx-kmerge-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let genome = letters("MjeNMV.fasta");
    let x = &genome[2000..2600];
    let query = dir.join("q.fna");
    std::fs::write(
        &query,
        record("q1 900 bases then X", &(genome[9000..9900].to_string() + x)) + &record("q2 X", x),
    )
    .unwrap();
    let subject = dir.join("s.fna");
    std::fs::write(&subject, record("s first 6 kb of MjeNMV", &genome[..6000])).unwrap();
    let (expected, off_stderr) = tblastx(&query, &subject, 1, &[STATS, SWITCH_OFF]);
    assert!(!expected.is_empty(), "no hits");
    assert_eq!(line(&off_stderr, "[X_STATS] TBXKMERGE "), None);
    for threads in THREADS {
        let what = format!("threads {threads}");
        let (report, stderr) = tblastx(&query, &subject, threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        let tie_runs = assert_kmerge_ran(&stderr, 1, &what);
        assert!(tie_runs > 0, "{what}:\n{stderr}");
    }
    for threads in [1, 4, 16] {
        let (report, _) = tblastx(&query, &subject, threads, &[]);
        assert_eq!(report, expected, "no statistics, threads {threads}");
    }
    let (report, stderr) = tblastx(&query, &subject, 4, &[STATS, SHADOW]);
    assert_eq!(report, expected, "shadow");
    let [subjects, hsps] = shadow_counts(&stderr);
    assert!(subjects == 1 && hsps > 0, "shadow:\n{stderr}");
    let chunks = ("LOSAT_TBLASTX_TEST_CHUNK_SIZE", "700");
    let (chunked, _) = tblastx(&query, &subject, 1, &[chunks, SWITCH_OFF]);
    assert!(!chunked.is_empty(), "no hits with chunks");
    for threads in [1, 4] {
        let what = format!("chunks, threads {threads}");
        let (report, stderr) = tblastx(&query, &subject, threads, &[STATS, chunks]);
        assert_eq!(report, chunked, "{what}");
        assert_kmerge_ran(&stderr, 1, &what);
        let [_, chunks_merged, _, _, _] = kmerge(&stderr);
        assert!(chunks_merged >= 18, "{what}:\n{stderr}");
    }
    let (report, stderr) = tblastx(&query, &subject, 4, &[STATS, SHADOW, chunks]);
    assert_eq!(report, chunked, "chunks, shadow");
    let [subjects, hsps] = shadow_counts(&stderr);
    assert!(subjects == 1 && hsps > 0, "chunks, shadow:\n{stderr}");
    let _ = std::fs::remove_dir_all(&dir);
}

//! The large per-split blocks of the TBLASTX unit path mapped by glibc (`LOSAT_X_SEEDUNMAP`,
//! `src/algorithm/tblastx/x_pair_par.rs`, the section above `X_MAPPED_BYTES`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_SEEDUNMAP`
//! is `0`, or `LOSAT_X_DEFAULT=0` and it is not set. It acts on Linux with glibc only, on the units
//! of `LOSAT_X_PAIRPAR` whose diagonal table has at least 128 KiB of cells and whose subject chunk
//! is at least a quarter of the table: their split's diagonal table and record vector, and a
//! unit's copy of its records of at least 128 KiB, get a capacity glibc always maps. It changes
//! capacities, never a value. Every test compares the report with the report of the same search
//! with `LOSAT_X_SEEDUNMAP=0` and reads the `[X_STATS] SEEDUNMAP` line, which is printed only while
//! the switch is on, to check where blocks were given the mapped capacity and where none was.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
//! ```c
//!         diag_table->hit_level_array = (DiagStruct *)
//!             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
//! ```
//! NCBI allocates one diagonal table per search thread; the switch chooses only the capacity of
//! LOSAT's per-split tables and lists, so the property under test is that the report does not
//! depend on it. That reports equal NCBI's is the business of the fixture comparisons, which run
//! on the default paths. The inputs are pieces of two fixture genomes (MjeNMV, MelaMJNV): 1,500
//! bases against 100,000 (contexts of 500 letters against subject frames of 33,333: tables of about
//! 33,800 cells, the chunk almost all of each; the rule holds), and the reverse (the same width,
//! a chunk of 500 letters: less than a quarter, so the arenas keep the blocks).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const SWITCH: (&str, &str) = ("LOSAT_X_SEEDUNMAP", "1");
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_SEEDUNMAP", "0");
const PREVIOUS: (&str, &str) = ("LOSAT_X_DEFAULT", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const THREADS: [u32; 2] = [1, 4];
const GLIBC: bool = cfg!(all(target_os = "linux", target_env = "gnu"));

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_X_DEFAULT` of the test's own environment is dropped as well: every test states the
/// default it runs with. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(switches.iter().copied());
    command
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

/// The query and subject files of one test: the first `query_bases` of MjeNMV and the first
/// `subject_bases` of MelaMJNV, in a directory of the test's own.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new(test: &str, query_bases: usize, subject_bases: usize) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-seed-unmap-{test}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.fna");
        let subject = dir.join("s.fna");
        std::fs::write(
            &query,
            record("mje MjeNMV piece", &letters("MjeNMV.fasta")[..query_bases]),
        )
        .unwrap();
        std::fs::write(
            &subject,
            record(
                "mela MelaMJNV piece",
                &letters("MelaMJNV.fasta")[..subject_bases],
            ),
        )
        .unwrap();
        Self {
            dir,
            query,
            subject,
        }
    }

    /// One TBLASTX search on standard output (`-outfmt 6`); returns the report and the standard
    /// error.
    fn tblastx(&self, threads: u32, switches: &[(&str, &str)]) -> (Vec<u8>, String) {
        let output = losat(switches)
            .arg("tblastx")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(&self.subject)
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

    /// The report with the switch off: no SEEDUNMAP line, at least one hit.
    fn report_without_the_switch(&self, threads: u32) -> Vec<u8> {
        let (report, stderr) = self.tblastx(threads, &[STATS, SWITCH_OFF]);
        assert!(!report.is_empty(), "threads {threads}: no hits");
        assert_eq!(
            line(&stderr, "[X_STATS] SEEDUNMAP "),
            None,
            "threads {threads}:\n{stderr}"
        );
        report
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// Contexts of 500 letters against subject frames of 33,333: the rule holds for every unit.
fn rule_holds(test: &str) -> Inputs {
    Inputs::new(test, 1_500, 100_000)
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

/// The units the unit path ran (`[X_STATS] PAIRPAR units=`).
fn units(stderr: &str) -> u64 {
    fields(stderr, "[X_STATS] PAIRPAR ", ["units"])[0]
}

/// The three counters of `[X_STATS] SEEDUNMAP`: diag, records, copies.
fn unmap(stderr: &str) -> [u64; 3] {
    fields(
        stderr,
        "[X_STATS] SEEDUNMAP ",
        ["diag", "records", "copies"],
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
// ```c
//         diag_table->hit_level_array = (DiagStruct *)
//             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
// ```
// Where the rule holds, with no variable for the switch (it is on by default), at one and four
// threads: the report is that of `LOSAT_X_SEEDUNMAP=0`; on Linux with glibc diagonal tables and
// record vectors got the mapped capacity (the units save too few records here for a copy of
// 128 KiB, which the unit tests and the larger searches of the entry in diagnostics.md reach),
// elsewhere none. `LOSAT_X_DEFAULT=0` turns it off (no line, the same report), and
// `LOSAT_X_DEFAULT=0` with `LOSAT_X_PAIRPAR=1` and the switch set maps the tables of the unit path
// (the record vector of `LOSAT_X_SEEDCOMPACT`, which is then off, does not exist).
#[test]
fn where_the_rule_holds_the_blocks_are_mapped_and_the_report_is_the_same() {
    let inputs = rule_holds("holds");
    for threads in THREADS {
        let what = format!("threads {threads}");
        let expected = inputs.report_without_the_switch(threads);
        let (report, stderr) = inputs.tblastx(threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        assert!(units(&stderr) > 0, "{what}:\n{stderr}");
        let [diag, records, copies] = unmap(&stderr);
        if GLIBC {
            assert!(diag > 0 && records > 0, "{what}:\n{stderr}");
        } else {
            assert_eq!([diag, records, copies], [0, 0, 0], "{what}:\n{stderr}");
        }

        let what = format!("LOSAT_X_DEFAULT=0 threads {threads}");
        let (report, stderr) = inputs.tblastx(threads, &[STATS, PREVIOUS]);
        assert_eq!(report, expected, "{what}");
        assert_eq!(
            line(&stderr, "[X_STATS] SEEDUNMAP "),
            None,
            "{what}:\n{stderr}"
        );

        let what = format!("previous paths with the unit path and the switch, threads {threads}");
        let (report, stderr) = inputs.tblastx(
            threads,
            &[STATS, PREVIOUS, ("LOSAT_X_PAIRPAR", "1"), SWITCH],
        );
        assert_eq!(report, expected, "{what}");
        let [diag, records, copies] = unmap(&stderr);
        assert_eq!(diag > 0, GLIBC, "{what}:\n{stderr}");
        assert_eq!([records, copies], [0, 0], "{what}:\n{stderr}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
// ```c
//         diag_table->hit_level_array = (DiagStruct *)
//             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
// ```
// With `LOSAT_X_SEEDCOMPACT=0` the units keep today's lists: only the diagonal tables are given
// the mapped capacity, and the report is the same.
#[test]
fn without_the_compact_lists_only_the_tables_are_mapped() {
    let inputs = rule_holds("nocompact");
    let threads = 4;
    let expected = inputs.report_without_the_switch(threads);
    let (report, stderr) = inputs.tblastx(threads, &[STATS, ("LOSAT_X_SEEDCOMPACT", "0")]);
    assert_eq!(report, expected);
    let [diag, records, copies] = unmap(&stderr);
    assert_eq!(diag > 0, GLIBC, "{stderr}");
    assert_eq!([records, copies], [0, 0], "{stderr}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
// ```c
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ```
// A short subject against a long query (subject frames of 500 letters, contexts of 33,333: tables
// of the same size as above, but the chunk is less than a quarter of them), and the search of the
// test above with `LOSAT_X_PAIRPAR=0` (no unit runs): the line is printed (the switch is on) with
// zero counters, and the report is that of `LOSAT_X_SEEDUNMAP=0`.
#[test]
fn where_the_rule_does_not_hold_nothing_is_mapped() {
    let short_subject = Inputs::new("short", 100_000, 1_500);
    for threads in THREADS {
        let what = format!("short subject, threads {threads}");
        let expected = short_subject.report_without_the_switch(threads);
        let (report, stderr) = short_subject.tblastx(threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        assert!(units(&stderr) > 0, "{what}:\n{stderr}");
        assert_eq!(unmap(&stderr), [0, 0, 0], "{what}:\n{stderr}");
    }
    let inputs = rule_holds("nounits");
    let threads = 4;
    let expected = inputs.report_without_the_switch(threads);
    let (report, stderr) = inputs.tblastx(threads, &[STATS, ("LOSAT_X_PAIRPAR", "0")]);
    assert_eq!(report, expected, "LOSAT_X_PAIRPAR=0");
    assert_eq!(unmap(&stderr), [0, 0, 0], "LOSAT_X_PAIRPAR=0:\n{stderr}");
}

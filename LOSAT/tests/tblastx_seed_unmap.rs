//! The large per-split blocks of the TBLASTX unit path mapped by glibc (`LOSAT_X_SEEDUNMAP`,
//! `src/algorithm/tblastx/x_pair_par.rs`, the section above `X_MAPPED_BYTES`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_SEEDUNMAP`
//! is `0`. It acts on Linux with glibc only, where transparent huge pages do not round mapped
//! blocks up (THP "madvise" or "never", or disabled for the process), at two or more threads (the search pool; at one thread nothing is retained), on
//! the units of the unit path whose query-group span and subject chunk both have at least
//! 65,536 letters: their split's diagonal table and record vector, and a unit's copy of its records
//! of at least 128 KiB, get a capacity glibc always maps. It changes capacities, never a value.
//! Every test compares the report with the report of the same search with `LOSAT_X_SEEDUNMAP=0`
//! and reads the `[X_STATS] SEEDUNMAP` line, which is printed only while the switch is on: its
//! `active` field says whether the switch acts in that process (it is 0 on a machine whose THP
//! setting is "always"), and the counters say where blocks were given the mapped capacity.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
//! ```c
//!         diag_table->hit_level_array = (DiagStruct *)
//!             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
//! ```
//! NCBI allocates one diagonal table per search thread; the switch chooses only the capacity of
//! LOSAT's per-split tables and lists, so the property under test is that the report does not
//! depend on it. That reports equal NCBI's is the business of the fixture comparisons, which run
//! on the default paths. The inputs are pieces of two fixture genomes (MelaMJNV, PemoMJNVA): 210,000
//! bases against 210,000 (contexts and subject frames of about 70,000 letters: the rule holds), and
//! 3,000 bases against 210,000 and the reverse (one side of about 1,000 letters: the arenas keep
//! the blocks).

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

const SWITCH_OFF: (&str, &str) = ("LOSAT_X_SEEDUNMAP", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const GLIBC: bool = cfg!(all(target_os = "linux", target_env = "gnu"));
/// The bases of a long piece: frames of about 70,000 letters, above the rule's 65,536.
const LONG: usize = 210_000;
/// The bases of a short piece.
const SHORT: usize = 3_000;

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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

/// The query and subject files of one test: the first `query_bases` of MelaMJNV and the first
/// `subject_bases` of PemoMJNVA, in a directory of the test's own.
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
            record(
                "mela MelaMJNV piece",
                &letters("MelaMJNV.fasta")[..query_bases],
            ),
        )
        .unwrap();
        std::fs::write(
            &subject,
            record(
                "pemo PemoMJNVA piece",
                &letters("PemoMJNVA.fasta")[..subject_bases],
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

/// Both sides long: the rule holds for every unit.
fn rule_holds(test: &str) -> Inputs {
    Inputs::new(test, LONG, LONG)
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

/// Whether the switch acts in the process (`active`), and the three counters of
/// `[X_STATS] SEEDUNMAP`: diag, records, copies. `active` is 0 or 1, and 0 off Linux with glibc.
fn unmap(stderr: &str) -> (bool, [u64; 3]) {
    let [active, diag, records, copies] = fields(
        stderr,
        "[X_STATS] SEEDUNMAP ",
        ["active", "diag", "records", "copies"],
    );
    assert!(active <= 1, "{stderr}");
    if !GLIBC {
        assert_eq!(active, 0, "{stderr}");
    }
    (active == 1, [diag, records, copies])
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
// ```c
//         diag_table->hit_level_array = (DiagStruct *)
//             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
// ```
// Where the rule holds, with no variable for the switch (it is on by default): the report is that
// of `LOSAT_X_SEEDUNMAP=0` at one and four threads; at one thread (the serial path keeps today's
// capacities) the counters are 0; at four threads, where the switch acts (`active=1`), diagonal
// tables and record vectors got the mapped capacity, and where it does not, none. The units save
// too few records here for a copy of 128 KiB, which the unit tests and the larger searches of the
// entry in diagnostics.md reach.
#[test]
fn where_the_rule_holds_the_blocks_are_mapped_and_the_report_is_the_same() {
    let inputs = rule_holds("holds");

    let expected = inputs.report_without_the_switch(1);
    let (report, stderr) = inputs.tblastx(1, &[STATS]);
    assert_eq!(report, expected, "threads 1");
    assert!(units(&stderr) > 0, "threads 1:\n{stderr}");
    assert_eq!(unmap(&stderr).1, [0, 0, 0], "threads 1:\n{stderr}");

    let threads = 4;
    let expected = inputs.report_without_the_switch(threads);
    let (report, stderr) = inputs.tblastx(threads, &[STATS]);
    assert_eq!(report, expected, "threads 4");
    assert!(units(&stderr) > 0, "threads 4:\n{stderr}");
    let (active, [diag, records, copies]) = unmap(&stderr);
    if active {
        assert!(diag > 0 && records > 0, "threads 4:\n{stderr}");
    } else {
        assert_eq!([diag, records, copies], [0, 0, 0], "threads 4:\n{stderr}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
// ```c
//         diag_table->hit_level_array = (DiagStruct *)
//             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
// ```
// With `LOSAT_X_SEEDCOMPACT=0` the units keep today's lists: only the diagonal tables are given
// the mapped capacity (where the switch acts), and the report is the same.
#[test]
fn without_the_compact_lists_only_the_tables_are_mapped() {
    let inputs = rule_holds("nocompact");
    let threads = 4;
    let expected = inputs.report_without_the_switch(threads);
    let (report, stderr) = inputs.tblastx(threads, &[STATS, ("LOSAT_X_SEEDCOMPACT", "0")]);
    assert_eq!(report, expected);
    let (active, [diag, records, copies]) = unmap(&stderr);
    assert_eq!(diag > 0, active, "{stderr}");
    assert_eq!([records, copies], [0, 0], "{stderr}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
// ```c
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ```
// A short query against a long subject (query-group spans of about 1,000 letters, chunks of about
// 70,000) and the reverse, at four threads: the line is printed (the switch is on) with zero
// counters, and the report is that of `LOSAT_X_SEEDUNMAP=0`.
#[test]
fn where_the_rule_does_not_hold_nothing_is_mapped() {
    let threads = 4;
    for (test, query_bases, subject_bases) in [("shortq", SHORT, LONG), ("shorts", LONG, SHORT)] {
        let inputs = Inputs::new(test, query_bases, subject_bases);
        let what = format!("{query_bases} bases against {subject_bases}");
        let expected = inputs.report_without_the_switch(threads);
        let (report, stderr) = inputs.tblastx(threads, &[STATS]);
        assert_eq!(report, expected, "{what}");
        assert!(units(&stderr) > 0, "{what}:\n{stderr}");
        assert_eq!(unmap(&stderr).1, [0, 0, 0], "{what}:\n{stderr}");
    }
}

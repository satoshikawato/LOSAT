//! The SEG of a TBLASTX query batch's frames on the search pool (`LOSAT_X_TBXQFRAMES`, in the
//! strict set: on unless `=0`; `src/algorithm/tblastx/x_query_frames.rs`) on whole runs of the
//! binary.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1281
//! ```c
//!     for (context = query_info->first_context;
//!          context <= query_info->last_context; ++context) {
//! ```
//! The switch chooses between two LOSAT code paths with the same result: the frames masked one
//! after the other, or one task per frame on the search pool. So the property under test is that
//! the report does not depend on it: for several queries with low-complexity runs (pieces of
//! MelaMJNV, a query of four letters, whose four frames hold one residue each, and one of two
//! letters, which has no frame) at 1 and 4 threads, the searches with the switch off and by
//! default give the same bytes; with several query batches (`BATCH_SIZE`) too; and `-seg no` is not
//! acted on. That the frames are masked on the pool's threads is checked by the unit test
//! `frames_are_masked_on_the_pool_threads`; here the reach rows (`LOSAT_X_REACH`,
//! `tests/reach_support`) show that the pool path ran at 4 threads and not at 1. That those bytes
//! equal NCBI's is the business of the TBLASTX fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

#[allow(clippy::duplicate_mod)]
mod reach_support;

const OFF: (&str, &str) = ("LOSAT_X_TBXQFRAMES", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of the first record of a FASTA file.
fn letters(name: &str) -> Vec<u8> {
    std::fs::read(fasta(name))
        .expect("read FASTA")
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A file name in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_tbxqframes_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The query and subject files, removed when dropped. Queries: 8 pieces of the first 60,000
/// letters of MelaMJNV, of 1,200 to 3,300 letters, every second with a run of codons that SEG
/// masks in some frames (`CAG`, `GCA`, `AAAAAG`, `CCTCCA` repeats), then a query of four letters
/// and one of two letters. Subject: those 60,000 letters.
struct Inputs {
    queries: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let genome = letters("MelaMJNV.fasta");
        let genome = &genome[..60_000];
        let mut text = Vec::new();
        let mut put = |title: String, seq: &[u8]| {
            text.extend_from_slice(format!(">{title}\n").as_bytes());
            for line in seq.chunks(60) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
        };
        for i in 0..8usize {
            let start = i * 7_100 + 211 * i;
            let len = 1_200 + (i * 977) % 2_100;
            let mut seq = genome[start..start + len].to_vec();
            if i % 2 == 0 {
                let unit: &[u8] = [b"CAG".as_slice(), b"GCA", b"AAAAAG", b"CCTCCA"][(i / 2) % 4];
                let run: Vec<u8> = unit.iter().copied().cycle().take(150).collect();
                let at = len / 3;
                seq.splice(at..at, run);
            }
            put(format!("piece{i}"), &seq);
        }
        put("four".into(), b"ACGT");
        put("two".into(), b"AC");
        let queries = temp_path("queries.fa");
        std::fs::write(&queries, text).expect("write the queries");
        let mut text = b">subject\n".to_vec();
        for line in genome.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
        let subject = temp_path("subject.fa");
        std::fs::write(&subject, text).expect("write the subject");
        Inputs { queries, subject }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.queries);
        let _ = std::fs::remove_file(&self.subject);
    }
}

/// One search, `-outfmt 6` on standard output, with `vars` as its only `LOSAT_*` variables (and
/// `BATCH_SIZE`). `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe`
/// feature (`cargo test --all-features`, as in CI) needs that file name as soon as a search has
/// more than one thread, and it selects nothing. Returns the report and the standard error.
fn tblastx(
    inputs: &Inputs,
    extra: &[&str],
    threads: u32,
    vars: &[(&str, &str)],
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
    let output = command
        .arg("tblastx")
        .arg("-query")
        .arg(&inputs.queries)
        .arg("-subject")
        .arg(&inputs.subject)
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{extra:?} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// `env` and `more` as one list of variables for `tblastx`.
fn with<'a>(
    env: &[(&'static str, &'static str)],
    more: &[(&'a str, &'a str)],
) -> Vec<(&'a str, &'a str)> {
    let mut vars: Vec<(&'a str, &'a str)> = env.to_vec();
    vars.extend_from_slice(more);
    vars
}

/// The search with the switch off at one thread, then off and by default at 1 and 4 threads, all
/// with the same bytes and no standard error. The default runs write the reach rows: the frames
/// are masked on the pool at 4 threads and not at 1.
fn check(env: &[(&'static str, &'static str)], extra: &[&str]) {
    let inputs = Inputs::new();
    let (expected, stderr) = tblastx(&inputs, extra, 1, &with(env, &[OFF]));
    assert!(!expected.is_empty(), "{env:?} {extra:?}: no hits");
    assert_eq!(stderr, "", "{env:?} {extra:?}");
    for threads in [1, 4] {
        let (report, stderr) = tblastx(&inputs, extra, threads, &with(env, &[OFF]));
        assert!(
            report == expected,
            "{env:?} {extra:?} {OFF:?} threads {threads}: the report differs"
        );
        assert_eq!(stderr, "", "{env:?} {extra:?} {OFF:?} threads {threads}");
        let reach = reach_support::Reach::new("qframes");
        let (report, stderr) = tblastx(&inputs, extra, threads, &with(env, &[reach.var()]));
        assert!(
            report == expected,
            "{env:?} {extra:?} threads {threads}: the report differs"
        );
        assert_eq!(stderr, "", "{env:?} {extra:?} threads {threads}");
        let sides = reach.sides("LOSAT_X_TBXQFRAMES");
        assert_eq!(
            sides.0,
            threads > 1,
            "{env:?} {extra:?} threads {threads}: sides {sides:?}"
        );
    }
}

// The default query batches (TBLASTX's 10,002 letters: two batches of the queries' 20 kb): the
// 10 queries' frames (8 x 6 for the pieces, 4 of one residue for the four letters, none for the
// two) are masked on the pool, some of them with intervals.
#[test]
fn query_frames_masked_on_the_pool() {
    check(&[], &[]);
}

// More query batches (`BATCH_SIZE` of 5,000 letters): each batch's frames are masked on the
// pool.
#[test]
fn query_frames_of_several_batches() {
    check(&[("BATCH_SIZE", "5000")], &[]);
}

// Other SEG parameters (`-seg 10 1.8 2.1`): other intervals, the same report with and without
// the switch.
#[test]
fn other_seg_parameters() {
    check(&[], &["-seg", "10 1.8 2.1"]);
}

// `-seg no` masks nothing: the pool path does not act, and the report is the same.
#[test]
fn seg_no_is_not_acted_on() {
    let inputs = Inputs::new();
    let extra = ["-seg", "no"];
    let (expected, _) = tblastx(&inputs, &extra, 1, &[OFF]);
    assert!(!expected.is_empty());
    for threads in [1, 4] {
        let reach = reach_support::Reach::new("qframes_seg_no");
        let (report, _) = tblastx(&inputs, &extra, threads, &[reach.var()]);
        assert!(report == expected, "threads {threads}: default differs");
        let sides = reach.sides("LOSAT_X_TBXQFRAMES");
        assert!(!sides.0, "threads {threads}: sides {sides:?}");
    }
}

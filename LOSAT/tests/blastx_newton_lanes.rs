//! BLASTX's composition-based redo with the Newton problems of a match's frame calls solved ahead
//! (`LOSAT_X_BXNEWTONLANES`, `x_bx_newton_lanes_prefetch_frames` in
//! `src/core/composition_adjustment/redo_alignment.rs`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1232-1244
//! ```c
//!                 if (compo_adjust_mode != eNoCompositionBasedStats &&
//!                         (subject_is_translated || hsp_index == 0
//!                                 || (nearIdenticalStatus != oldNearIdenticalStatus))) {
//!                     Blast_AminoAcidComposition subject_composition;
//!                     s_GetComposition(&subject_composition,
//!                             alphsize, &subject,
//!                             &window->subject_range,
//!                             in_align, FALSE, subject_is_translated);
//!                     adjust_search_failed =
//!                             Blast_AdjustScores(matrix, query_composition,
//!                                     query.length,
//!                                     &subject_composition,
//!                                     subject.length,
//! ```
//! The switch is in the strict set (on unless it is `0`). A `Blast_AdjustScores` call uses a result
//! solved ahead only when its Newton input has the same bits, so the property under test is that
//! the report does not depend on the switch, at one thread and at four, and that the default run
//! solved problems ahead and used them (`[X_STATS] NEWTONLANES prefetched= hits=`). That reports equal NCBI's is the business of the
//! BLASTX fixtures, which run the default.
//!
//! The query is the first 60,000 bases of AP027131 (genetic code 4), the subjects the first 200
//! proteins of the related AP027133, so that matches hold several frames and windows.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const OFF: (&str, &str) = ("LOSAT_X_BXNEWTONLANES", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_newton_lanes_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The first `n` records of a FASTA file, as text.
fn head(name: &str, n: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut out = String::new();
    let mut records = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            records += 1;
            if records > n {
                break;
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    assert_eq!(records, n + 1, "{name} has more than {n} records");
    out
}

/// The first `n` bases of the one record of a FASTA file, as FASTA text.
fn bases(name: &str, n: usize) -> String {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut lines = text.lines();
    let header = lines.next().expect("a header");
    let sequence: String = lines.flat_map(str::chars).take(n).collect();
    assert_eq!(sequence.len(), n, "{name} is shorter than {n} bases");
    format!("{header}\n{sequence}\n")
}

/// The inputs of the searches, removed when dropped.
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("query.fa");
        let subject = dir.join("subject.faa");
        std::fs::write(&query, bases("AP027131.fasta", 60_000)).unwrap();
        std::fs::write(&subject, head("AP027133.faa", 200)).unwrap();
        Inputs {
            dir,
            query,
            subject,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// One search (`-outfmt 6`) with `env` as its only variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` names a file of its own: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. The run must exit with 0.
/// Returns the report and the standard error.
fn blastx(inputs: &Inputs, threads: u32, env: &[(&str, &str)]) -> (Vec<u8>, String) {
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
    let output = command
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .envs(env.iter().copied())
        .arg("blastx")
        .arg("-query")
        .arg(&inputs.query)
        .arg("-subject")
        .arg(&inputs.subject)
        .args(["-query_gencode", "4", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "threads {threads} {env:?}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The number after `key=` on the `[X_STATS] NEWTONLANES ` line.
fn lanes(stderr: &str, key: &str) -> u64 {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] NEWTONLANES "))
        .unwrap_or_else(|| panic!("no NEWTONLANES line:\n{stderr}"));
    line.split_whitespace()
        .find_map(|field| field.strip_prefix(key)?.strip_prefix('='))
        .unwrap_or_else(|| panic!("no {key}= in {line}"))
        .parse()
        .expect("a counter")
}

// The default report equals `LOSAT_X_BXNEWTONLANES=0`'s at one thread and at four, its standard
// error stays empty, and with `LOSAT_X_STATS=1` the default run shows problems solved ahead and
// used, `=0` none. At one thread every Newton problem of `=0` (all solved where they arise, the
// misses) is a problem the default either took from the store or solved itself.
#[test]
fn the_report_does_not_depend_on_the_switch_and_the_lanes_run() {
    let inputs = Inputs::new();
    for threads in [1, 4] {
        let (expected, expected_stderr) = blastx(&inputs, threads, &[OFF]);
        assert!(!expected.is_empty(), "threads {threads}: no hits");
        assert_eq!(expected_stderr, "", "threads {threads}");
        let (report, stderr) = blastx(&inputs, threads, &[]);
        assert_eq!(report, expected, "threads {threads}");
        assert_eq!(stderr, "", "threads {threads}");
        let (report, stats) = blastx(&inputs, threads, &[STATS]);
        assert_eq!(report, expected, "threads {threads}");
        let (report, stats_off) = blastx(&inputs, threads, &[STATS, OFF]);
        assert_eq!(report, expected, "threads {threads}");
        let hits = lanes(&stats, "hits");
        assert!(lanes(&stats, "prefetched") >= hits, "{stats}");
        assert!(
            hits > 0,
            "threads {threads}: no problem used from the store\n{stats}"
        );
        assert_eq!(lanes(&stats_off, "prefetched"), 0, "{stats_off}");
        assert_eq!(lanes(&stats_off, "hits"), 0, "{stats_off}");
        if threads == 1 {
            assert_eq!(
                hits + lanes(&stats, "misses"),
                lanes(&stats_off, "misses"),
                "the same Newton problems\n{stats}\n{stats_off}"
            );
        }
    }
}

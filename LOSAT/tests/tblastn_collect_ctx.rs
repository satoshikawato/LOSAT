//! The preliminary collector of tblastn split by query context (`LOSAT_X_TBNCOLLECTCTX`,
//! `src/algorithm/tblastn/stage_d_pipeline.rs` `PreliminaryCollector::x_add_by_context`) on a
//! whole search.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. The hit lists are
//! those of one filter per context (the unit test in `stage_d_pipeline.rs` compares every field on
//! random lists), so the property under test is that the report is `=0`'s at 1 and 2 threads,
//! with and without the counters. The `[X_STATS] TBNCOLLECTCTX` line (`LOSAT_X_STATS=1`) shows
//! that the collector split lists that touch several contexts. Input: the 115 proteins of
//! PemoMJNVA x MelaMJNV + LvMJNV (two subjects of 291 and 285 kb, so each collector takes two
//! lists and 2 threads search the subjects side by side).

use std::path::{Path, PathBuf};
use std::process::Command;

const OFF: (&str, &str) = ("LOSAT_X_TBNCOLLECTCTX", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

fn work_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("tblastn_collect_ctx");
    std::fs::create_dir_all(&dir).expect("work directory");
    dir
}

/// MelaMJNV and LvMJNV in one FASTA file of this process.
fn two_subjects() -> PathBuf {
    let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let mut text = String::new();
    for name in ["MelaMJNV.fasta", "LvMJNV.fasta"] {
        text += &std::fs::read_to_string(fasta.join(name)).expect("subject FASTA");
        if !text.ends_with('\n') {
            text.push('\n');
        }
    }
    let path = work_dir().join(format!("subjects.{}.fasta", std::process::id()));
    std::fs::write(&path, text).expect("write the subjects");
    path
}

fn tblastn(subjects: &Path, threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(
            WORKER_LOG,
            work_dir().join(format!("worker.{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(fasta.join("PemoMJNVA.faa"))
        .arg("-subject")
        .arg(subjects)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// lists, contexts, hsps of the `[X_STATS] TBNCOLLECTCTX` line, if there is one.
fn stats(stderr: &str) -> Option<[u64; 3]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNCOLLECTCTX "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(["lists", "contexts", "hsps"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

// No NCBI counterpart (LOSAT_X_TBNCOLLECTCTX): the default report equals `=0`'s at 1 and 2
// threads, with and without LOSAT_X_STATS, and 1 thread's equals 2 threads'; the default splits
// lists that touch more contexts than there are lists, `=0` splits none.
#[test]
fn the_collector_split_by_context_gives_the_report_of_the_filters() {
    let subjects = two_subjects();
    let (off, _) = tblastn(&subjects, 1, &[OFF]);
    assert!(!off.is_empty(), "no hits");
    for threads in [1, 2] {
        let (off_stats, stderr) = tblastn(&subjects, threads, &[OFF, STATS]);
        assert!(
            off_stats == off,
            "threads {threads}: =0 against =0 at 1 thread"
        );
        assert_eq!(stats(&stderr), None, "threads {threads}: =0 splits none");
        let (report, _) = tblastn(&subjects, threads, &[]);
        assert!(report == off, "threads {threads}: default against =0");
        let (report, stderr) = tblastn(&subjects, threads, &[STATS]);
        assert!(
            report == off,
            "threads {threads}: default with stats against =0"
        );
        let [lists, contexts, hsps] = stats(&stderr).expect("a TBNCOLLECTCTX line");
        assert!(
            lists >= 2 && contexts > lists && hsps >= contexts,
            "threads {threads}: lists {lists}, contexts {contexts}, hsps {hsps}"
        );
    }
    std::fs::remove_file(&subjects).ok();
}

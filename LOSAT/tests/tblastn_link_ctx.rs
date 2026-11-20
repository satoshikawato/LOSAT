//! The preliminary linking of tblastn one query context at a time (`LOSAT_X_TBNLINKCTX`,
//! `src/algorithm/tblastn/stage_d_linking/x_link_ctx.rs`) on a whole search.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. A subject list
//! with HSPs of two or more query contexts is linked one context at a time; a subject longer than a
//! chunk's budget links them on the pool its units ran on. The linked values are those of the one
//! pass (the unit test of `x_link_ctx.rs` compares every field on random lists), so the property
//! under test is that the report does not depend on the switch or the threads. The
//! `[X_STATS] TBNLINKCTX` line (`LOSAT_X_STATS=1`) shows that lists were linked apart, and on the
//! pool where the subject's units ran. Input: PemoMJNVA x MelaMJNV (one subject of 291 kb): its first
//! 20 proteins (one query batch, whose units run on the pool from two threads) and all 115 (query
//! batches side by side, each on one thread: the contexts are linked apart one after another).

use std::path::{Path, PathBuf};
use std::process::Command;

const OFF: (&str, &str) = ("LOSAT_X_TBNLINKCTX", "0");

/// The first 20 proteins of PemoMJNVA (a file of the test's own) and all 115.
fn queries() -> &'static [PathBuf; 2] {
    static FILES: std::sync::OnceLock<[PathBuf; 2]> = std::sync::OnceLock::new();
    FILES.get_or_init(|| {
        let all = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PemoMJNVA.faa");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_link_ctx");
        std::fs::create_dir_all(&dir).expect("work directory");
        let first = dir.join(format!("first20.{}.faa", std::process::id()));
        let text = std::fs::read_to_string(&all).expect("query file");
        let end = text
            .match_indices('>')
            .nth(20)
            .map_or(text.len(), |(at, _)| at);
        std::fs::write(&first, &text[..end]).expect("query file");
        [first, all]
    })
}

fn tblastn(query: &Path, threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("tblastn_link_ctx");
        std::fs::create_dir_all(&dir).expect("work directory");
        command.env(
            WORKER_LOG,
            dir.join(format!("worker.{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(fasta.join("MelaMJNV.fasta"))
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .env("LOSAT_X_STATS", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// lists, contexts, pooled of the `[X_STATS] TBNLINKCTX` line, if there is one.
fn stats(stderr: &str) -> Option<[u64; 3]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNLINKCTX "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(["lists", "contexts", "pooled"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

// No NCBI counterpart (LOSAT_X_TBNLINKCTX): the default report equals `=0`'s at 1 and 2 threads;
// the default links lists apart (on the pool at 2 threads when the subject's units run there, never
// at 1), `=0` links none apart.
#[test]
fn linking_one_context_at_a_time_gives_the_report_of_one_pass() {
    for (query, pool_at_two) in queries().iter().zip([true, false]) {
        for threads in [1, 2] {
            let label = format!("{query:?}, threads {threads}");
            let (off, stderr) = tblastn(query, threads, &[OFF]);
            assert!(!off.is_empty(), "{label}: no hits");
            assert_eq!(stats(&stderr), None, "{label}: =0 links none apart");
            let (report, stderr) = tblastn(query, threads, &[]);
            assert!(report == off, "{label}: against =0");
            let [lists, contexts, pooled] = stats(&stderr).expect("a TBNLINKCTX line");
            assert!(lists > 0, "{label}: lists {lists}, contexts {contexts}");
            let expected = if threads == 2 && pool_at_two {
                lists
            } else {
                0
            };
            assert_eq!(pooled, expected, "{label}: lists linked on the pool");
        }
    }
}

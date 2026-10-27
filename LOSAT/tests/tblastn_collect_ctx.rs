//! The preliminary collector of tblastn split by query context in one pass
//! (`LOSAT_X_TBNCOLLECTCTX`, `src/algorithm/tblastn/stage_d_pipeline.rs`
//! `PreliminaryCollector::x_add_by_context`) on a whole search.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. The hit lists are
//! those of one filter per context (the unit test in `stage_d_pipeline.rs` compares every field on
//! random lists), so the property under test is that the report does not depend on the switch or
//! the threads. The `[X_STATS] TBNCOLLECTCTX` line (`LOSAT_X_STATS=1`) shows that the collector
//! split lists in one pass. Input: the 115 proteins of PemoMJNVA x MelaMJNV (one subject of 291 kb).

use std::path::{Path, PathBuf};
use std::process::Command;

const OFF: (&str, &str) = ("LOSAT_X_TBNCOLLECTCTX", "0");

fn tblastn(threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
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
        let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("tblastn_collect_ctx");
        std::fs::create_dir_all(&dir).expect("work directory");
        command.env(
            WORKER_LOG,
            dir.join(format!("worker.{}.log", std::process::id())),
        );
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(fasta.join("PemoMJNVA.faa"))
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

/// lists, hsps of the `[X_STATS] TBNCOLLECTCTX` line, if there is one.
fn stats(stderr: &str) -> Option<[u64; 2]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with("[X_STATS] TBNCOLLECTCTX "))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(["lists", "hsps"].map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

// No NCBI counterpart (LOSAT_X_TBNCOLLECTCTX): the default report equals `=0`'s at 1 and 2
// threads; the default splits lists in one pass, `=0` splits none.
#[test]
fn the_collector_split_in_one_pass_gives_the_report_of_the_filters() {
    for threads in [1, 2] {
        let (off, stderr) = tblastn(threads, &[OFF]);
        assert!(!off.is_empty(), "threads {threads}: no hits");
        assert_eq!(stats(&stderr), None, "threads {threads}: =0 splits none");
        let (report, stderr) = tblastn(threads, &[]);
        assert!(report == off, "threads {threads}: against =0");
        let [lists, hsps] = stats(&stderr).expect("a TBNCOLLECTCTX line");
        assert!(
            lists > 0 && hsps > 0,
            "threads {threads}: lists {lists}, hsps {hsps}"
        );
    }
}

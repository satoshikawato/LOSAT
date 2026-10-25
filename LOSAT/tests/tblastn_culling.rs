//! TBLASTN `-culling_limit` (NCBI's culling writer of the preliminary stage and its culling
//! pipe after the traceback) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the option is read as NCBI reads it, and a
//! search with a limit keeps a subset of the HSPs of the search without it, at every thread
//! count the same. That the bytes equal NCBI BLAST+ `tblastn -culling_limit` is the business
//! of the NCBI parity probe (`ncbi-parity-20261019/O3b`).

use std::path::{Path, PathBuf};
use std::process::Command;

use LOSAT::cli::{try_parse_from, Cli, Commands};

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta/outfmt0")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path (`LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature needs that file name as soon
/// as a search has more than one thread, and it selects nothing).
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

fn limit(words: &[&str]) -> Result<i32, String> {
    let mut argv = vec!["losat", "tblastn", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend(words);
    match try_parse_from::<Cli, _, _>(argv) {
        Ok(cli) => match cli.command {
            Commands::Tblastn(args) => Ok(args.culling_limit_value()),
            _ => unreachable!("tblastn"),
        },
        Err(error) => Err(LOSAT::cli::render_message(&error)),
    }
}

// NCBI blast_args.cpp (CHspFilteringArgs): `-culling_limit` is an integer >= 0, default 0,
// hexadecimal accepted as NCBI's argument parser does.
#[test]
fn the_option_is_read_as_ncbi_reads_it() {
    assert_eq!(limit(&[]), Ok(0));
    assert_eq!(limit(&["-culling_limit", "0"]), Ok(0));
    assert_eq!(limit(&["-culling_limit", "1"]), Ok(1));
    assert_eq!(limit(&["-culling_limit", "0x2"]), Ok(2));
    assert_eq!(limit(&["-culling_limit", "2147483647"]), Ok(2147483647));
    for bad in ["-1", "1.5", "abc", "2147483648"] {
        let error = limit(&["-culling_limit", bad]).unwrap_err();
        assert!(!error.contains("is not supported by"), "{bad}: {error}");
    }
}

// The best-hit options stay refused (they are not ported), with or without a limit.
#[test]
fn the_best_hit_options_stay_refused() {
    for words in [
        &["-best_hit_overhang", "0.1"][..],
        &["-culling_limit", "1", "-best_hit_score_edge", "0.1"],
    ] {
        let error = limit(words).unwrap_err();
        assert!(
            error.contains("is not supported by LOSAT's TBLASTN"),
            "{words:?}: {error}"
        );
    }
}

fn search(cbs: &str, extra: &[&str], threads: u32) -> Vec<String> {
    let output = losat()
        .arg("tblastn")
        .arg("-query")
        .arg(fasta("e2e_many_query.faa"))
        .arg("-subject")
        .arg(fasta("e2e_many_subject.fna"))
        .args(["-outfmt", "6"])
        .args(["-evalue", "1000", "-comp_based_stats", cbs])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

#[test]
fn a_limit_keeps_a_subset_of_the_hsps_at_every_thread_count() {
    let all = search("0", &[], 1);
    assert!(!all.is_empty());
    let zero = search("0", &["-culling_limit", "0"], 1);
    assert_eq!(zero, all, "a limit of 0 is no culling");
    let culled = search("0", &["-culling_limit", "1"], 1);
    assert!(!culled.is_empty());
    assert!(
        culled.len() < all.len(),
        "limit 1 keeps {} of {} HSPs",
        culled.len(),
        all.len()
    );
    for row in &culled {
        assert!(all.contains(row), "culling made a row of its own: {row}");
    }
    // A limit above the number of overlapping HSPs culls nothing.
    assert_eq!(search("0", &["-culling_limit", "100000"], 1), all);
    // The same rows at other thread counts.
    for threads in [2, 4] {
        assert_eq!(
            search("0", &["-culling_limit", "1"], threads),
            culled,
            "{threads} threads"
        );
    }
    // CBS 2 runs the pipe after the redo.
    let cbs2_all = search("2", &[], 1);
    let cbs2 = search("2", &["-culling_limit", "1"], 1);
    assert!(cbs2.len() <= cbs2_all.len());
    for row in &cbs2 {
        assert!(
            cbs2_all.contains(row),
            "culling made a row of its own: {row}"
        );
    }
}

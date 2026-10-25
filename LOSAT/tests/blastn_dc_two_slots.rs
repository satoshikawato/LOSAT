//! dc-megablast's two LOSAT_X_QBATCHPAR slots of half the threads where the slot model takes
//! more (`LOSAT_X_DCTWOSLOTS`, `src/algorithm/blastn/blast_engine/x_qbatch_par.rs`) on whole
//! searches.
//!
//! The switch is in the strict set: on unless it is `0`. It chooses only how many slots search
//! the planned query batches ahead of the batch loop, whose results the loop takes under the keys
//! it computes itself, so the property under test is that the report does not depend on it:
//! pieces of LC738874 against LC738874, one query per batch (`BATCH_SIZE=1`), with the
//! slot budget that NCBI-sized query files get by default (`LOSAT_X_QBATCHPARMEM`, 28 MiB per
//! thread beyond the first), give the same bytes on standard output and standard error (without
//! the `[X_STATS]` lines) with the default and with `LOSAT_X_DCTWOSLOTS=0`. The QBATCHPAR line
//! of `LOSAT_X_STATS` and `LOSAT_X_REACH` show where the cap acted: for 60 queries at 16 threads
//! the model takes three slots or more for dc-megablast and the default runs 2 x 8 (the switch's
//! own side); for 3 queries (two batches after the first) at 8 threads it takes 2 x 4 itself (the
//! reference side); megablast has no cap. That the rows
//! equal NCBI's is the business of the blastn fixtures; the unit tests in `x_qbatch_par.rs`
//! compare the capped choice with the model's on random inputs.

#![cfg(feature = "parallel")]

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

mod reach_support;

const OFF: (&str, &str) = ("LOSAT_X_DCTWOSLOTS", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const ONE_QUERY_PER_BATCH: (&str, &str) = ("BATCH_SIZE", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A query file removed when dropped: `count` pieces of LC738874 of 300 to 2,700 letters, every
/// fifth reverse-complemented.
struct Queries(PathBuf, usize);

impl Queries {
    fn new(count: usize) -> Queries {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let letters: Vec<u8> = std::fs::read(fasta("LC738874.fasta"))
            .expect("read FASTA")
            .split(|&byte| byte == b'\n')
            .skip(1)
            .flat_map(|line| line.iter().copied().filter(|b| !b.is_ascii_whitespace()))
            .collect();
        let mut text = String::new();
        for i in 0..count {
            let start = (i * 4_513) % (letters.len() - 3_000);
            let len = 300 + (i * 977) % 2_400;
            let mut piece = letters[start..start + len].to_vec();
            if i % 5 == 4 {
                piece.reverse();
                for base in piece.iter_mut() {
                    *base = match base.to_ascii_uppercase() {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        other => other,
                    };
                }
            }
            text.push_str(&format!(">piece_{i} LC738874 {start}+{len}\n"));
            for line in piece.chunks(70) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
        }
        let path = std::env::temp_dir().join(format!(
            "losat_dctwoslots_{}_{}.fa",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, text).expect("write queries");
        Queries(path, count)
    }
}

impl Drop for Queries {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// One search's standard output, its standard error without the `[X_STATS]` lines, and the
/// fields of its `[X_STATS] QBATCHPAR` line.
struct Run {
    stdout: Vec<u8>,
    stderr: String,
    qbatchpar: HashMap<String, String>,
}

/// `blastn -task <task> -outfmt 6` of the queries against LC738874 at `threads` threads, one
/// query per batch, the default slot budget of a query file that NCBI splits, `LOSAT_X_STATS`,
/// and no `LOSAT_*` variable of the caller but `switches` (`LOSAT_BLASTX_WORKER_LOG` stays: a
/// binary built with the `blastx-worker-probe` feature needs it at two threads or more).
fn blastn(queries: &Queries, task: &str, threads: u64, switches: &[(&str, &str)]) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    let budget = (28 * (threads - 1)).to_string();
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(&queries.0)
        .arg("-subject")
        .arg(fasta("LC738874.fasta"))
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .env("LOSAT_X_QBATCHPARMEM", &budget)
        .envs([ONE_QUERY_PER_BATCH, STATS])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let what = format!("{task} threads {threads} {switches:?}");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{what}: {stderr}");
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] QBATCHPAR "))
        .unwrap_or_else(|| panic!("{what}: no QBATCHPAR line: {stderr}"));
    let qbatchpar = line
        .split_whitespace()
        .filter_map(|field| field.split_once('='))
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect();
    let stderr = stderr
        .lines()
        .filter(|line| !line.starts_with("[X_STATS]"))
        .map(|line| format!("{line}\n"))
        .collect();
    Run {
        stdout: output.stdout,
        stderr,
        qbatchpar,
    }
}

impl Run {
    fn stat(&self, name: &str) -> usize {
        self.qbatchpar[name].parse().expect("a count")
    }

    /// (k, m, rule_k, rule_m) of the QBATCHPAR line.
    fn layouts(&self) -> (usize, usize, usize, usize) {
        (
            self.stat("k"),
            self.stat("m"),
            self.stat("rule_k"),
            self.stat("rule_m"),
        )
    }
}

#[test]
fn the_report_does_not_depend_on_the_cap() {
    let (many, few) = (Queries::new(60), Queries::new(3));
    // (task, queries, threads, cap acts)
    for (task, queries, threads, capped) in [
        ("dc-megablast", &many, 16, true),
        ("dc-megablast", &few, 8, false),
        ("megablast", &many, 16, false),
    ] {
        let what = format!("{task} {} queries threads {threads}", queries.1);
        let reach = reach_support::Reach::new(&format!("dctwoslots-{task}-{threads}"));
        let default = blastn(queries, task, threads, &[reach.var()]);
        let off = blastn(queries, task, threads, &[OFF]);
        let rows = default
            .stdout
            .split(|&b| b == b'\n')
            .filter(|l| !l.is_empty())
            .count();
        assert!(rows >= queries.1, "{what}: {rows} rows");
        assert!(
            default.stdout == off.stdout,
            "{what}: the default's report differs from LOSAT_X_DCTWOSLOTS=0's"
        );
        assert_eq!(default.stderr, off.stderr, "{what}");
        for run in [&default, &off] {
            assert_eq!(run.qbatchpar["off"], "none", "{what}: {:?}", run.qbatchpar);
            assert!(run.stat("precomputed") > 0, "{what}: {:?}", run.qbatchpar);
        }
        let (k, m, rule_k, rule_m) = off.layouts();
        assert_eq!(
            (k, m),
            (rule_k, rule_m),
            "{what}: =0 runs the model's layout"
        );
        assert_eq!(
            off.qbatchpar["predicted"], default.qbatchpar["predicted"],
            "{what}"
        );
        if capped {
            assert!(rule_k > 2, "{what}: {:?}", off.qbatchpar);
            assert_eq!(
                default.layouts(),
                (2, threads as usize / 2, rule_k, rule_m),
                "{what}: {:?}",
                default.qbatchpar
            );
        } else {
            assert_eq!(default.layouts(), off.layouts(), "{what}");
        }
        let sides = reach.sides("LOSAT_X_DCTWOSLOTS");
        let expected = match (task, capped) {
            ("dc-megablast", true) => (true, false),
            ("dc-megablast", false) => (false, true),
            _ => (false, false),
        };
        assert_eq!(sides, expected, "{what}: sides of LOSAT_X_DCTWOSLOTS");
    }
}

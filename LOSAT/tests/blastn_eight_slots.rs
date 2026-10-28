//! blastn's best layout of at most eight LOSAT_X_QBATCHPAR slots where the slot model takes
//! more (`LOSAT_X_BNEIGHTSLOTS`, `src/algorithm/blastn/blast_engine/x_qbatch_par.rs`) on whole
//! searches.
//!
//! The switch is in the strict set: on unless it is `0` (the slots themselves, `LOSAT_X_QBATCHPAR`,
//! are off by default since 2026-10-28, so every run here sets `LOSAT_X_QBATCHPAR=1`). It chooses only how many slots search
//! the planned query batches ahead of the batch loop, whose results the loop takes under the keys
//! it computes itself, so the property under test is that the report does not depend on it:
//! pieces of LC738874 against LC738874, one query per batch (`BATCH_SIZE=1`), with the
//! slot budget that NCBI-sized query files get by default (`LOSAT_X_QBATCHPARMEM`, 96 MiB per
//! thread beyond the first), give the same bytes on standard output and standard error (without
//! the `[X_STATS]` lines) with the default and with `LOSAT_X_BNEIGHTSLOTS=0`. The QBATCHPAR line
//! of `LOSAT_X_STATS` and `LOSAT_X_REACH` show where the eight-slot rule acted: for 40 queries at
//! 16 threads (28 batch units) the model takes 14 one-thread slots for `-task blastn` and the
//! default runs 7 x 2 (the switch's own side); for 60 queries (43 units) it takes 16 x 1 and the
//! default runs 8 x 2; for 40 queries at 11 threads it takes 11 x 1 and the rule declines (fewer
//! than 12 threads); for 3 queries (two batches after the first) at 8 threads it takes 2 x 4
//! itself (the reference side in these two cases); blastn-short has no rule; with `=0` no side
//! is recorded, nor where the model does not run (2 queries: one batch after the first) or
//! `LOSAT_X_QBATCHPARSLOTS` sets the slots (12 x 1, more than eight, with and without the
//! switch). That the rows equal NCBI's is the business of the blastn fixtures; the unit tests in
//! `x_qbatch_par.rs` compare the ruled choice with the model's on random inputs and on the
//! measured cells.

#![cfg(feature = "parallel")]

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

mod reach_support;

const OFF: (&str, &str) = ("LOSAT_X_BNEIGHTSLOTS", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");
const ONE_QUERY_PER_BATCH: (&str, &str) = ("BATCH_SIZE", "1");
/// The slots are off by default since 2026-10-28; every run here is meant to use them.
const SLOTS_ON: (&str, &str) = ("LOSAT_X_QBATCHPAR", "1");

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
            "losat_bneightslots_{}_{}.fa",
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
            || matches!(
                &*text,
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    let budget = (96 * (threads - 1)).to_string();
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
        .envs([ONE_QUERY_PER_BATCH, STATS, SLOTS_ON])
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
fn the_report_does_not_depend_on_the_rule() {
    let (forty, sixty, three) = (Queries::new(40), Queries::new(60), Queries::new(3));
    // (task, queries, threads, the rule's side (None: no rule), the model's (k, units), the
    // slots that run with the rule)
    for (task, queries, threads, side, model, ruled) in [
        ("blastn", &forty, 16, Some(true), (14, 28), 7),
        ("blastn", &sixty, 16, Some(true), (16, 43), 8),
        ("blastn", &forty, 11, Some(false), (11, 28), 11),
        ("blastn", &three, 8, Some(false), (2, 2), 2),
        ("blastn-short", &forty, 16, None, (14, 28), 14),
    ] {
        let what = format!("{task} {} queries threads {threads}", queries.1);
        let reach = reach_support::Reach::new(&format!("bneightslots-{task}-{threads}"));
        let default = blastn(queries, task, threads, &[reach.var()]);
        let off_reach = reach_support::Reach::new(&format!("bneightslots-off-{task}-{threads}"));
        let off = blastn(queries, task, threads, &[OFF, off_reach.var()]);
        let rows = default
            .stdout
            .split(|&b| b == b'\n')
            .filter(|l| !l.is_empty())
            .count();
        assert!(rows >= queries.1, "{what}: {rows} rows");
        assert!(
            default.stdout == off.stdout,
            "{what}: the default's report differs from LOSAT_X_BNEIGHTSLOTS=0's"
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
            (rule_k, off.stat("units")),
            model,
            "{what}: the model's slots and the batch units"
        );
        assert_eq!(
            off.qbatchpar["predicted"], default.qbatchpar["predicted"],
            "{what}"
        );
        assert_eq!(
            default.layouts(),
            (ruled, threads as usize / ruled, rule_k, rule_m),
            "{what}: {:?}",
            default.qbatchpar
        );
        let sides = reach.sides("LOSAT_X_BNEIGHTSLOTS");
        let expected = match side {
            Some(true) => (true, false),
            Some(false) => (false, true),
            None => (false, false),
        };
        assert_eq!(sides, expected, "{what}: sides of LOSAT_X_BNEIGHTSLOTS");
        // A switch that is off records no side.
        assert_eq!(
            off_reach.sides("LOSAT_X_BNEIGHTSLOTS"),
            (false, false),
            "{what}: sides of LOSAT_X_BNEIGHTSLOTS=0"
        );
    }
}

// No side where the rule is not consulted: the model does not run (2 queries at 16 threads: one
// batch after the first, `off=single`), or `LOSAT_X_QBATCHPARSLOTS=12` sets the slots (12 x 1
// with and without the switch).
#[test]
fn no_side_without_the_model_s_slots() {
    let (forty, two) = (Queries::new(40), Queries::new(2));
    // (queries, extra variables, the slots that run, off=)
    for (queries, extra, layout, off) in [
        (&two, &[][..], (0, 0), "single"),
        (
            &forty,
            &[("LOSAT_X_QBATCHPARSLOTS", "12")][..],
            (12, 1),
            "none",
        ),
    ] {
        let what = format!("{} queries {extra:?}", queries.1);
        let reach = reach_support::Reach::new(&format!("bneightslots-none-{}", queries.1));
        let mut on: Vec<(&str, &str)> = extra.to_vec();
        on.push(reach.var());
        let default = blastn(queries, "blastn", 16, &on);
        let mut zero: Vec<(&str, &str)> = extra.to_vec();
        zero.push(OFF);
        let off_run = blastn(queries, "blastn", 16, &zero);
        assert!(
            default.stdout == off_run.stdout,
            "{what}: the default's report differs from LOSAT_X_BNEIGHTSLOTS=0's"
        );
        assert_eq!(default.stderr, off_run.stderr, "{what}");
        for run in [&default, &off_run] {
            assert_eq!(
                (run.stat("k"), run.stat("m")),
                layout,
                "{what}: {:?}",
                run.qbatchpar
            );
            assert_eq!(run.qbatchpar["off"], off, "{what}: {:?}", run.qbatchpar);
        }
        assert_eq!(
            reach.sides("LOSAT_X_BNEIGHTSLOTS"),
            (false, false),
            "{what}: sides of LOSAT_X_BNEIGHTSLOTS"
        );
    }
}

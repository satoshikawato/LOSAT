//! The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle): blastp's
//! `LOSAT_X_BPREDOPAR`, `LOSAT_X_BPREDOLPT`, `LOSAT_X_BPDBCUT`, `LOSAT_X_BPSUBJMETA` and
//! `LOSAT_X_BPLKPAHEAD`, and blastn's `LOSAT_X_PAIRPARDISC`, `LOSAT_X_PAIRPARFOLD`,
//! `LOSAT_X_QDUSTPAR`, `LOSAT_X_DIAGLV`, `LOSAT_X_QBATCHPAR`, `LOSAT_X_QPREPPAR` and
//! `LOSAT_X_RPTSTREAM`. Each has its own test file, whose searches run with the others off (the
//! default before). Here the searches run with the default (no variable) and
//! give the report of each switch off (`=0`) and of
//! all of them off, with the same standard error. And the default runs them: the
//! `[TIMING] x_bpredopar:` line of BPREDOPAR, with BPREDOLPT's `lpt_windows=`, comes without a
//! variable and not with `LOSAT_X_BPREDOPAR=0`. QBATCHPAR left the default set on 2026-10-28
//! (Owner decision: the slots' max RSS rise exceeds their time saving): it is off without a
//! variable (no `[X_STATS]` line) and runs only with `LOSAT_X_QBATCHPAR=1`, where its off reason
//! for these small query files is the budget, which is 0 by default for a query file that NCBI's
//! blastn would not split by queries (`x_qbatch_par.rs`, `Profile::of_search`), and
//! `LOSAT_X_QBATCHPARMEM` still sets one. `LOSAT_X_REACH`, set in the
//! default run of a search, shows that BPLKPAHEAD, BPSUBJMETA, BPDBCUT, DIAGLV, QDUSTPAR and
//! QPREPPAR take their own paths on these searches.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

/// The blastp switches of the bundle.
const BLASTP: &[&str] = &[
    "LOSAT_X_BPDBCUT",
    "LOSAT_X_BPLKPAHEAD",
    "LOSAT_X_BPREDOLPT",
    "LOSAT_X_BPREDOPAR",
    "LOSAT_X_BPSUBJMETA",
];

/// The blastn switches of the bundle.
const BLASTN: &[&str] = &[
    "LOSAT_X_DIAGLV",
    "LOSAT_X_PAIRPARDISC",
    "LOSAT_X_PAIRPARFOLD",
    "LOSAT_X_QBATCHPAR",
    "LOSAT_X_QDUSTPAR",
    "LOSAT_X_QPREPPAR",
    "LOSAT_X_RPTSTREAM",
];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of this test file (`label`: one per test).
fn work_dir(label: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("default_bundle")
        .join(label);
    std::fs::create_dir_all(&dir).expect("work directory");
    dir
}

/// One LOSAT run with `args` and `vars` as its only variables that choose a code path or NCBI's
/// batches (`BATCH_SIZE` too). `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. Must exit with 0; returns
/// the report and standard error.
fn losat(args: &[&str], vars: &[(&str, &str)]) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-default-bundle-worker.log"),
        );
    }
    let output = command
        .args(args)
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{args:?} {vars:?}: {stderr}");
    (output.stdout, stderr)
}

/// The search `args` (with `base` set in every run) gives one report and one standard error with
/// no variable, with each switch of `switches` at `0`, and with all of
/// them at `0`; returns the report. `first` is set in the first (default) run only, for example
/// `LOSAT_X_REACH` (which writes nothing to standard error): the `=0` arms must not get it, or
/// their rows would be ORed into the sides read afterwards.
fn one_report(
    args: &[&str],
    base: &[(&str, &str)],
    first: &[(&str, &str)],
    switches: &[&'static str],
) -> Vec<u8> {
    let with = |extra: Vec<(&'static str, &'static str)>| -> Vec<(&str, &str)> {
        base.iter().copied().chain(extra).collect()
    };
    let default_vars: Vec<(&str, &str)> = base.iter().chain(first).copied().collect();
    let (report, stderr) = losat(args, &default_vars);
    assert!(!report.is_empty(), "{args:?}: the search has hits");
    let mut arms = Vec::new();
    for &name in switches {
        arms.push(with(vec![(name, "0")]));
    }
    arms.push(with(
        BLASTP
            .iter()
            .chain(BLASTN)
            .map(|&name| (name, "0"))
            .collect(),
    ));
    for vars in &arms {
        let (other, other_stderr) = losat(args, vars);
        assert!(
            other == report,
            "{args:?}: the report differs with {vars:?}"
        );
        assert_eq!(
            other_stderr, stderr,
            "{args:?}: standard error with {vars:?}"
        );
    }
    report
}

/// The numbers of the first line of `stderr` that starts with `prefix` (`name=value` fields whose
/// value is a whole number).
fn counts(stderr: &str, prefix: &str) -> std::collections::HashMap<String, u64> {
    let line = stderr
        .lines()
        .find(|line| line.starts_with(prefix))
        .unwrap_or_else(|| panic!("no {prefix} line: {stderr}"));
    line.split_whitespace()
        .filter_map(|word| word.split_once('='))
        .filter_map(|(name, value)| Some((name.to_string(), value.parse().ok()?)))
        .collect()
}

/// A protein database of the subject proteome, made with `losat makeblastdb`; returns its name.
fn protein_database(label: &str) -> PathBuf {
    let name = work_dir(label).join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// Twelve pieces of LC738874 (1.5-3.7 kb) as queries, in a file of this test.
fn pieces(label: &str) -> PathBuf {
    let text = std::fs::read_to_string(fixture("LC738874.fasta")).expect("fixture");
    let genome: String = text.lines().filter(|line| !line.starts_with('>')).collect();
    let mut fasta = String::new();
    for i in 0..12usize {
        let start = i * 13_000;
        fasta.push_str(&format!(
            ">piece{i}\n{}\n",
            &genome[start..start + 1_500 + 200 * i]
        ));
    }
    let file = work_dir(label).join("pieces.fna");
    std::fs::write(&file, fasta).expect("query file");
    file
}

/// blastp-fast with `-subject` and exact blastp with `-db`, each in blocks (two blocks of the default
/// `LOSAT_X_BPQBLOCK` for the 77,448 letters of the query proteome) at four threads.
#[test]
fn blastp_gives_one_report_with_the_switches_on_and_off() {
    let query = fixture("MeenMJNV.faa");
    let subject = fixture("MejoMJNV.faa");
    let (query, subject) = (query.to_str().unwrap(), subject.to_str().unwrap());
    let fast = [
        "blastp",
        "-task",
        "blastp-fast",
        "-query",
        query,
        "-subject",
        subject,
        "-outfmt",
        "6",
        "-evalue",
        "1e-3",
        "-num_threads",
        "4",
    ];
    // The first (default) run of each search also writes the reach rows (S2: no extra run).
    let reach_fast = reach_support::Reach::new("blastp_fast");
    one_report(&fast, &[], &[reach_fast.var()], BLASTP);
    let db = protein_database("blastp");
    let db = db.to_str().unwrap();
    let exact = [
        "blastp",
        "-query",
        query,
        "-db",
        db,
        "-outfmt",
        "6",
        "-evalue",
        "1e-3",
        "-num_threads",
        "4",
    ];
    let reach_db = reach_support::Reach::new("blastp_db");
    one_report(&exact, &[], &[reach_db.var()], BLASTP);
    // The rows come from the default runs whose reports every `=0` arm above matched, so the
    // equalities compare these paths (a later change of scope or input would otherwise make
    // them vacuous). BPSUBJMETA: subject records were read by the metadata path. BPLKPAHEAD: the
    // look-ahead build was chosen for exact blastp's blocks.
    let sides = reach_fast.sides("LOSAT_X_BPSUBJMETA");
    assert!(
        sides.0,
        "LOSAT_X_BPSUBJMETA: sides {sides:?}: its path did not run"
    );
    let sides = reach_db.sides("LOSAT_X_BPLKPAHEAD");
    assert!(
        sides.0,
        "LOSAT_X_BPLKPAHEAD: sides {sides:?}: its path did not run"
    );
    // BPDBCUT: the exact `-db` search was cut into blocks of the switch's own size.
    let sides = reach_db.sides("LOSAT_X_BPDBCUT");
    assert!(
        sides.0,
        "LOSAT_X_BPDBCUT: sides {sides:?}: its path did not run"
    );
    // The default runs BPREDOPAR's path (and BPREDOLPT's windows); its switch at 0 does not.
    let timing = [("LOSAT_TIMING", "1")];
    let (_, stderr) = losat(&fast, &timing);
    let lines: Vec<&str> = stderr
        .lines()
        .filter(|line| line.starts_with("[TIMING] x_bpredopar: "))
        .collect();
    assert!(
        !lines.is_empty(),
        "no x_bpredopar line by default: {stderr}"
    );
    assert!(
        lines.iter().all(|line| line.contains(" lpt_windows=")),
        "{lines:?}"
    );
    let (_, stderr) = losat(&fast, &[("LOSAT_TIMING", "1"), ("LOSAT_X_BPREDOPAR", "0")]);
    assert!(!stderr.contains("[TIMING] x_bpredopar: "), "{stderr}");
}

/// dc-megablast and blastn with `-subject`, several query batches (`BATCH_SIZE=4000`), at four
/// threads.
#[test]
fn blastn_gives_one_report_with_the_switches_on_and_off() {
    let query = pieces("blastn");
    let subject = fixture("LC738868.fasta");
    let (query, subject) = (query.to_str().unwrap(), subject.to_str().unwrap());
    let batches = [("BATCH_SIZE", "4000")];
    for task in ["dc-megablast", "blastn", "megablast"] {
        let args = [
            "blastn",
            "-task",
            task,
            "-query",
            query,
            "-subject",
            subject,
            "-outfmt",
            "6",
            "-num_threads",
            "4",
        ];
        // The default run also writes the reach rows (no extra run).
        let reach = reach_support::Reach::new(task);
        one_report(&args, &batches, &[reach.var()], BLASTN);
        // The default runs the strips stages: PAIRPARDISC's for dc-megablast, PAIRPAR's
        // (megablast, blastn) otherwise, each fold split by PAIRPARFOLD (partitions), and
        // RPTSTREAM writes batch by batch; their switches at 0 do not.
        let stage = if task == "dc-megablast" {
            "[TIMING] x_pairpardisc_seed:"
        } else {
            "[TIMING] x_pairpar_seed:"
        };
        let observed = [
            ("BATCH_SIZE", "4000"),
            ("LOSAT_TIMING", "1"),
            ("LOSAT_X_STATS", "1"),
        ];
        let (_, stderr) = losat(&args, &observed);
        let seeds = counts(&stderr, stage);
        assert!(
            seeds["chunks"] > 0 && seeds["partitions"] > 0,
            "{task}: {stderr}"
        );
        let streamed = counts(&stderr, "[X_STATS] RPTSTREAM");
        assert!(streamed["batches"] > 1, "{task}: {stderr}");
        assert!(
            stderr
                .lines()
                .any(|line| line.starts_with("[X_STATS] RPTSTREAM ") && line.contains(" off=none")),
            "{task}: {stderr}"
        );
        let folded_off = [
            ("BATCH_SIZE", "4000"),
            ("LOSAT_TIMING", "1"),
            ("LOSAT_X_PAIRPARFOLD", "0"),
        ];
        let (_, stderr) = losat(&args, &folded_off);
        let unfolded = counts(&stderr, stage);
        assert_eq!(
            unfolded.get("partitions").copied().unwrap_or(0),
            0,
            "{task}: {stderr}"
        );
        // The rows of the default run above. DIAGLV: the reference fold took the diag-hash path
        // with level bounds (its mark is in the reference fold only, not in the strips folds).
        // QDUSTPAR: the query dusting ran on the pool. QPREPPAR: the pool ran the batch (the
        // mark fires for any batch with more than one context at 2+ threads).
        for switch in ["LOSAT_X_DIAGLV", "LOSAT_X_QDUSTPAR", "LOSAT_X_QPREPPAR"] {
            let sides = reach.sides(switch);
            assert!(
                sides.0,
                "{task} {switch}: sides {sides:?}: its path did not run"
            );
        }
    }
    // QBATCHPAR is off by default (since 2026-10-28) and gives no `[X_STATS]` line; with
    // `LOSAT_X_QBATCHPAR=1` these query files get no default budget, and the variable still sets
    // one.
    let args = [
        "blastn",
        "-task",
        "dc-megablast",
        "-query",
        query,
        "-subject",
        subject,
        "-outfmt",
        "6",
        "-num_threads",
        "4",
    ];
    let (_, stderr) = losat(&args, &[("BATCH_SIZE", "4000"), ("LOSAT_X_STATS", "1")]);
    assert!(!stderr.contains("[X_STATS] QBATCHPAR "), "{stderr}");
    let stats = |vars: &[(&str, &str)]| -> String {
        let (_, stderr) = losat(&args, vars);
        stderr
            .lines()
            .find(|line| line.starts_with("[X_STATS] QBATCHPAR "))
            .unwrap_or_else(|| panic!("no QBATCHPAR line with {vars:?}: {stderr}"))
            .to_string()
    };
    let line = stats(&[
        ("BATCH_SIZE", "4000"),
        ("LOSAT_X_STATS", "1"),
        ("LOSAT_X_QBATCHPAR", "1"),
    ]);
    assert!(
        line.contains(" budget_mib=0 ") && line.ends_with(" off=budget"),
        "{line}"
    );
    let line = stats(&[
        ("BATCH_SIZE", "4000"),
        ("LOSAT_X_STATS", "1"),
        ("LOSAT_X_QBATCHPAR", "1"),
        ("LOSAT_X_QBATCHPARMEM", "84"),
    ]);
    assert!(line.contains(" budget_mib=84 "), "{line}");
    let (_, stderr) = losat(
        &args,
        &[
            ("BATCH_SIZE", "4000"),
            ("LOSAT_X_STATS", "1"),
            ("LOSAT_X_QBATCHPAR", "0"),
        ],
    );
    assert!(!stderr.contains("[X_STATS] QBATCHPAR "), "{stderr}");
    // With a budget and two slots forced, the slots search batches ahead of the loop, and the
    // report is that of the default before (every switch of the bundle at 0).
    let slots = [
        ("BATCH_SIZE", "4000"),
        ("LOSAT_X_STATS", "1"),
        ("LOSAT_X_QBATCHPAR", "1"),
        ("LOSAT_X_QBATCHPARMEM", "4096"),
        ("LOSAT_X_QBATCHPARSLOTS", "2"),
    ];
    let (report, stderr) = losat(&args, &slots);
    let speculated = counts(&stderr, "[X_STATS] QBATCHPAR ");
    assert_eq!((speculated["k"], speculated["m"]), (2, 2), "{stderr}");
    let before: Vec<(&str, &str)> = [("BATCH_SIZE", "4000")]
        .into_iter()
        .chain(BLASTP.iter().chain(BLASTN).map(|&name| (name, "0")))
        .collect();
    let (expected, _) = losat(&args, &before);
    assert!(report == expected, "the report differs with two slots");
}

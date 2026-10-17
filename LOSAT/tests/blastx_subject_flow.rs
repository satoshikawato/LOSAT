//! BLASTX's subject loop searched as a flow of tasks (`LOSAT_X_BXSUBJFLOW`,
//! `src/algorithm/blastx/x_subject_flow.rs`) on whole searches.
//!
//! The switch is off unless set (it is not in the strict set). With two or more threads the
//! default path searches a query chunk's subjects in waves of `threads x 64` and hands each wave
//! to the sink after its slowest worker; with the switch the other threads search tasks of `k`
//! consecutive subjects ahead and the loop's thread hands them to the sink in subject order. So
//! the property under test is that the report does not depend on it: standard output equals that
//! of the same search without the switch, for two, three and four threads, several `k` (64, the
//! default, and 2, 7 and 16), outfmt 0, 6 and 7, a database (`-db`) and a
//! FASTA subject set (`-subject`), and an ungapped search (whose subjects carry link cutoffs).
//! The shadow (`LOSAT_X_BXSUBJFLOWSHADOW`) searches every task again on another table and must
//! find the same subjects. `LOSAT_TIMING` shows that the flow ran (and the waves did not). The
//! fallback to the waves where a diagonal offset could pass 2^30 (`-window_size` near 2^28, a
//! diagonal table of 1 GiB or more) is a unit test of the module. That the bytes equal NCBI's
//! is the business of the BLASTX fixtures and the campaign's comparison runs.
//!
//! The query is 9,000 letters of AP027131 (one batch, one chunk, so the pool searches its
//! subjects); the subjects are the proteins of AP027131, AP027132 and AP027133 (1,711 records).

#![cfg(feature = "parallel")]

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// The temporary directory of this test process.
fn temp_dir() -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-blastx-subject-flow-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The query (letters 100,000-109,000 of AP027131) and the subject set (its proteins and those
/// of AP027132 and AP027133) as files, and the subject set's `LOSAT makedb` database.
fn inputs(dir: &Path) -> (PathBuf, PathBuf, String) {
    let genome = std::fs::read_to_string(fixture("AP027131.fasta")).unwrap();
    let letters: String = genome
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(str::chars)
        .filter(|c| !c.is_ascii_whitespace())
        .skip(100_000)
        .take(9_000)
        .collect();
    let query = dir.join("query.fna");
    std::fs::write(
        &query,
        format!(">piece AP027131 100000-109000\n{letters}\n"),
    )
    .unwrap();
    let set = dir.join("set.faa");
    let mut bytes = Vec::new();
    for name in ["AP027131.faa", "AP027132.faa", "AP027133.faa"] {
        bytes.extend(std::fs::read(fixture(name)).unwrap());
    }
    std::fs::write(&set, bytes).unwrap();
    let db = dir.join("set").to_str().unwrap().to_string();
    let output = losat()
        .arg("makedb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "prot", "-out", &db])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, set, db)
}

/// One search; returns standard output and standard error.
fn blastx(query: &Path, args: &[&str], threads: u32, env: &[(&str, &str)]) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(query)
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{args:?} {env:?}: {stderr}");
    (output.stdout, stderr)
}

/// The numbers of a `[TIMING] blastx_subject_flow` or `blastx_subject_waves` line: `name=value`.
fn timing(stderr: &str, line: &str, name: &str) -> u64 {
    let text = stderr
        .lines()
        .find(|l| l.starts_with(&format!("[TIMING] {line}: ")))
        .unwrap_or_else(|| panic!("no {line} line: {stderr}"));
    text.split_whitespace()
        .find_map(|word| word.strip_prefix(&format!("{name}=")))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in {text}"))
}

#[test]
fn the_report_does_not_depend_on_the_subject_flow() {
    let dir = temp_dir();
    let (query, set, db) = inputs(&dir);
    let set = set.to_str().unwrap();
    for (subjects, outfmts) in [
        (["-db", db.as_str()], &["6", "0", "7"][..]),
        (["-subject", set], &["6"][..]),
    ] {
        for &outfmt in outfmts {
            let args = [subjects[0], subjects[1], "-outfmt", outfmt];
            let (reference, _) = blastx(&query, &args, 3, &[]);
            assert!(
                reference.iter().filter(|&&b| b == b'\n').count() > 20,
                "{args:?}: too few hits to test"
            );
            for threads in [2, 3, 4] {
                for k in ["1", "2", "7", "16"] {
                    let (report, _) = blastx(&query, &args, threads, &[("LOSAT_X_BXSUBJFLOW", k)]);
                    assert!(
                        report == reference,
                        "{args:?} threads {threads} LOSAT_X_BXSUBJFLOW={k}: report differs"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn an_ungapped_search_and_the_shadow_give_the_same_report() {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-subject-flow-ungapped-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let (query, _, db) = inputs(&dir);
    for extra in [&[][..], &["-ungapped", "-comp_based_stats", "0"][..]] {
        let mut args = vec!["-db", db.as_str(), "-outfmt", "6"];
        args.extend_from_slice(extra);
        let (reference, _) = blastx(&query, &args, 4, &[]);
        let (report, stderr) = blastx(
            &query,
            &args,
            4,
            &[
                ("LOSAT_X_BXSUBJFLOW", "2"),
                ("LOSAT_X_BXSUBJFLOWSHADOW", "1"),
            ],
        );
        assert!(report == reference, "{extra:?}: report differs");
        let tasks = 1_711_u64.div_ceil(2);
        assert!(
            stderr.contains(&format!(
                "[X_SHADOW] BXSUBJFLOW tasks_compared={tasks} (subject results identical)"
            )),
            "{extra:?}: {stderr}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_flow_takes_the_subject_loop() {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-subject-flow-timing-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let (query, _, db) = inputs(&dir);
    let args = ["-db", db.as_str(), "-outfmt", "6"];
    let (reference, stderr) = blastx(&query, &args, 2, &[("LOSAT_TIMING", "1")]);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 0);
    assert!(timing(&stderr, "blastx_subject_waves", "waves") > 0);
    // The 1,711 subjects as 27 tasks of 64, and no wave.
    let (report, stderr) = blastx(
        &query,
        &args,
        2,
        &[("LOSAT_X_BXSUBJFLOW", "1"), ("LOSAT_TIMING", "1")],
    );
    assert!(report == reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 1);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "tasks"), 27);
    assert_eq!(timing(&stderr, "blastx_subject_waves", "waves"), 0);
    // One thread: no pool, so the subject loop of the waves runs (on this thread).
    let (report, stderr) = blastx(
        &query,
        &args,
        1,
        &[("LOSAT_X_BXSUBJFLOW", "1"), ("LOSAT_TIMING", "1")],
    );
    assert!(report == reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 0);
    std::fs::remove_dir_all(&dir).unwrap();
}

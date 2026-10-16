//! The DUST masks of a query batch on the search pool (`LOSAT_X_QDUSTPAR`,
//! `src/algorithm/blastn/x_query_dust_par.rs`) on whole searches.
//!
//! The switch is off unless set to any value (it is not in the strict set). It chooses between
//! two LOSAT code paths with the same result: the queries of a batch masked one after the other,
//! or on the search pool, one query per job, collected in query order. So the property under
//! test is that the report does not depend on it: for megablast, blastn and dc-megablast with
//! several queries (pieces of LC738874 with microsatellites, runs of `N`, a poly-A query and a
//! query shorter than a triplet) at 1 and 4 threads, the searches without the switch, with it,
//! and with the shadow (`LOSAT_X_QDUSTPARSHADOW`, which makes the masks serially as well and
//! compares them) give the same bytes, and `-dust no` is not acted on (the shadow compares
//! nothing). That those bytes equal NCBI's is the business of the blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: (&str, &str) = ("LOSAT_X_QDUSTPAR", "1");
const SHADOW: (&str, &str) = ("LOSAT_X_QDUSTPARSHADOW", "1");

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
        "losat_qdustpar_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A query file removed when dropped: 24 pieces of LC738874 of 900 to 4,500 letters, every third
/// with a microsatellite, every fourth with a run of `N`s and lower case letters, then a poly-A
/// query, a query of two letters, an empty-after-the-title-free short one and one with only `N`.
struct Queries(PathBuf);

impl Queries {
    fn new() -> Self {
        let genome = letters("LC738874.fasta");
        let mut text = Vec::new();
        let mut put = |title: String, mut seq: Vec<u8>| {
            text.extend_from_slice(format!(">{title}\n").as_bytes());
            seq.push(b'\n');
            for line in seq.chunks(61) {
                text.extend_from_slice(line);
            }
            if !text.ends_with(b"\n") {
                text.push(b'\n');
            }
        };
        for i in 0..24usize {
            let start = i * 11_000 + 137 * i;
            let len = 900 + (i * 1_531) % 3_600;
            let mut seq = genome[start..start + len].to_vec();
            if i % 3 == 0 {
                let at = len / 3;
                let unit: &[u8] = [b"CA".as_slice(), b"GAA", b"T", b"ACGGT"][(i / 3) % 4];
                let run: Vec<u8> = unit.iter().copied().cycle().take(180).collect();
                seq.splice(at..at, run);
            }
            if i % 4 == 0 {
                let at = len / 2;
                seq.splice(at..at, std::iter::repeat_n(b'N', 25));
                for byte in &mut seq[at + 40..at + 140] {
                    *byte = byte.to_ascii_lowercase();
                }
            }
            put(format!("piece{i}"), seq);
        }
        put("polyA".into(), vec![b'A'; 700]);
        put("two".into(), b"AC".to_vec());
        put("short".into(), b"ACGTACGTACGTTTGACAGT".to_vec());
        let mut mostly_n = vec![b'N'; 300];
        mostly_n.extend_from_slice(&genome[5_000..5_400]);
        put("mostlyN".into(), mostly_n);
        let path = temp_path("queries.fa");
        std::fs::write(&path, text).expect("write the queries");
        Queries(path)
    }
}

impl Drop for Queries {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// One search of the queries against LC738868, `-outfmt 6` on standard output, with the given
/// variables as its only `LOSAT_*` variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as
/// in CI) needs that file name as soon as a search has more than one thread, and it selects
/// nothing. Returns the report and the standard error.
fn blastn(
    queries: &Queries,
    extra: &[&str],
    threads: u32,
    switches: &[(&str, &str)],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(&queries.0)
        .arg("-subject")
        .arg(fasta("LC738868.fasta"))
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{extra:?} {switches:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The numbers of `[X_SHADOW] QDUSTPAR batches=B queries=Q intervals=I identical` in `stderr`.
fn shadow_counts(stderr: &str) -> (u64, u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] QDUSTPAR "))
        .unwrap_or_else(|| panic!("no shadow line:\n{stderr}"));
    assert!(line.ends_with(" identical"), "{line}");
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
    };
    assert_eq!(stderr.lines().count(), 1, "{stderr}");
    (field("batches="), field("queries="), field("intervals="))
}

/// For one task: the search without the switch at one thread, then off, on and the shadow at 1
/// and 4 threads, all with the same bytes; the switch writes nothing, the shadow one line, which
/// counts masks made on the pool only at 4 threads.
fn check_task(task: &str, extra: &[&str]) {
    let queries = Queries::new();
    let mut args = vec!["-task", task];
    args.extend_from_slice(extra);
    let (expected, expected_stderr) = blastn(&queries, &args, 1, &[]);
    assert!(!expected.is_empty(), "{task}: no hits");
    assert_eq!(expected_stderr, "", "{task}");
    for threads in [1, 4] {
        let (report, stderr) = blastn(&queries, &args, threads, &[]);
        assert!(report == expected, "{task} threads {threads}: off differs");
        assert_eq!(stderr, "", "{task} threads {threads}");
        let (report, stderr) = blastn(&queries, &args, threads, &[SWITCH]);
        assert!(report == expected, "{task} threads {threads}: on differs");
        assert_eq!(stderr, "", "{task} threads {threads}");
        let (report, stderr) = blastn(&queries, &args, threads, &[SHADOW]);
        assert!(
            report == expected,
            "{task} threads {threads}: shadow differs"
        );
        let (batches, count, intervals) = shadow_counts(&stderr);
        if threads == 1 {
            assert_eq!((batches, count, intervals), (0, 0, 0), "{task}");
        } else {
            assert!(
                batches > 0 && count > 1 && intervals > 0,
                "{task}: {stderr}"
            );
        }
    }
}

#[test]
fn megablast_queries_masked_on_the_pool() {
    check_task("megablast", &[]);
}

#[test]
fn blastn_queries_masked_on_the_pool() {
    check_task("blastn", &[]);
}

#[test]
fn dc_megablast_queries_masked_on_the_pool() {
    check_task("dc-megablast", &[]);
}

// Another window and level, and the lower-case masking that merges with the DUST masks.
#[test]
fn dust_parameters_and_lowercase_masking() {
    check_task("megablast", &["-dust", "10 32 2", "-lcase_masking"]);
}

// `-dust no` makes no masks: the shadow has nothing to compare, and the report is the same.
#[test]
fn dust_no_is_not_acted_on() {
    let queries = Queries::new();
    let args = ["-task", "megablast", "-dust", "no"];
    let (expected, _) = blastn(&queries, &args, 1, &[]);
    assert!(!expected.is_empty());
    for threads in [1, 4] {
        let (report, _) = blastn(&queries, &args, threads, &[SWITCH]);
        assert!(report == expected, "threads {threads}: on differs");
        let (report, stderr) = blastn(&queries, &args, threads, &[SHADOW]);
        assert!(report == expected, "threads {threads}: shadow differs");
        assert_eq!(shadow_counts(&stderr), (0, 0, 0));
    }
}

// With LOSAT_TIMING (or LOSAT_X_STATS) the DUST step prints its wall time at exit, with or without
// the switch.
#[test]
fn timing_line_is_printed() {
    let queries = Queries::new();
    let args = ["-task", "megablast"];
    let (expected, _) = blastn(&queries, &args, 1, &[]);
    for switches in [
        &[("LOSAT_TIMING", "1")][..],
        &[("LOSAT_TIMING", "1"), SWITCH][..],
    ] {
        let (report, stderr) = blastn(&queries, &args, 4, switches);
        assert!(report == expected, "{switches:?}: the report differs");
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] query_dust: "))
            .unwrap_or_else(|| panic!("{switches:?}: no timing line:\n{stderr}"));
        assert!(line.ends_with('s'), "{line}");
    }
    let (report, stderr) = blastn(&queries, &args, 4, &[("LOSAT_X_STATS", "1"), SWITCH]);
    assert!(report == expected, "the report differs with LOSAT_X_STATS");
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] QDUSTPAR dust_s="))
        .unwrap_or_else(|| panic!("no stats line:\n{stderr}"));
    assert!(
        line["[X_STATS] QDUSTPAR dust_s=".len()..]
            .parse::<f64>()
            .is_ok(),
        "{line}"
    );
    assert!(!stderr.contains("[TIMING] query_dust"), "{stderr}");
}

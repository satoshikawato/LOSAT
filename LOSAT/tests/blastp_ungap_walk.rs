//! blastp's ungapped X-drop walks with `LOSAT_X_BPWALK` (`src/algorithm/blastp/x_ungapped_fast.rs`)
//! on whole searches.
//!
//! The switch is not in the strict set (`src/utils/xswitch.rs`): it is off unless its variable is
//! set, and it acts only where `LOSAT_X_BPUNGAP` (or its shadow `LOSAT_X_BPUNGAPSHADOW`) runs the
//! fast ungapped extension; there the walks keep the drop below the best score instead of the
//! best score. These tests compare the reports with the reference (no variable) and with
//! `LOSAT_X_BPUNGAP` alone, at one and four threads, for two-hit and one-hit (`-window_size 0`)
//! searches, in the tabular and the pairwise report, on a generated pair (related proteins with
//! long ungapped stretches, repeats, unrelated proteins) and on the fixture pair; the switch
//! without `LOSAT_X_BPUNGAP` must change nothing; and they run the shadow, which compares the
//! drop walks' extension with the reference for every call and names the walk in its line.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:903-916
//! ```c
//!     for (i = n; i >= 0; i--) {
//!         score += matrix[q[i]][s[i]];
//!
//!         if (score > maxscore) {
//!             maxscore = score;
//!             best_i = i;
//!         }
//!         ...
//!         if ((maxscore - score) >= dropoff)
//!             break;
//!     }
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// One blastp run of `query` against `subject` with the given variables as its only `LOSAT_*`
/// variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing; without the caller's value
/// a file of the tests is named. `BATCH_SIZE` and the chunk sizes are removed. The run must exit
/// with 0. Returns the standard output and the standard error.
fn blastp(
    query: &Path,
    subject: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    extra: &[&str],
) -> (Vec<u8>, String) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpwalk-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

fn next(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 33) as u32
}

/// A generated query and subject set in the test's own directory: random proteins, a protein
/// with long repeats, and subjects that are copies of the queries with every seventh letter
/// changed (long ungapped extensions with drops and recoveries) or with a block of unrelated
/// letters inside (extensions that stop on the drop), plus unrelated proteins.
fn generated_pair() -> (PathBuf, PathBuf) {
    const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("blastp_ungap_walk");
    std::fs::create_dir_all(&dir).expect("test directory");
    let mut seed = 0x5eed_u64;
    let mut protein = |n: usize| -> Vec<u8> {
        (0..n)
            .map(|_| AMINO[next(&mut seed) as usize % AMINO.len()])
            .collect()
    };
    let mut queries: Vec<Vec<u8>> = (0..12).map(|i| protein(150 + 41 * i)).collect();
    let mut repeats = Vec::new();
    for _ in 0..10 {
        repeats.extend_from_slice(b"WKLSEELLKKAEEC");
    }
    repeats.extend(protein(80));
    queries.push(repeats);
    let mut subjects: Vec<Vec<u8>> = Vec::new();
    for (k, q) in queries.iter().enumerate() {
        let mut s: Vec<u8> = q
            .iter()
            .enumerate()
            .map(|(i, &a)| if i % 7 == 6 { AMINO[(i * 5) % 20] } else { a })
            .collect();
        if k % 2 == 1 {
            let at = s.len() / 2;
            let block = protein(12);
            s[at..at + 12].copy_from_slice(&block);
        }
        subjects.push(s);
    }
    subjects.extend((0..6).map(|i| protein(100 + 53 * i)));
    let write = |name: &str, seqs: &[Vec<u8>], tag: &str| {
        let mut text = String::new();
        for (i, seq) in seqs.iter().enumerate() {
            text.push_str(&format!(">{tag}{i}\n{}\n", String::from_utf8_lossy(seq)));
        }
        let path = dir.join(name);
        std::fs::write(&path, text).expect("write FASTA");
        path
    };
    (
        write("query.faa", &queries, "q"),
        write("subject.faa", &subjects, "s"),
    )
}

const UNGAP: (&str, &str) = ("LOSAT_X_BPUNGAP", "1");
const WALK: (&str, &str) = ("LOSAT_X_BPWALK", "1");

#[test]
fn drop_walks_give_the_reports_of_the_reference() {
    let (query, subject) = generated_pair();
    let pairs = [
        (query, subject),
        (fixture("MeenMJNV.faa"), fixture("MejoMJNV.faa")),
    ];
    for (query, subject) in &pairs {
        // two-hit, then one-hit (`-window_size 0`)
        for window in [&[][..], &["-window_size", "0"][..]] {
            for threads in [1, 4] {
                for outfmt in ["6", "0"] {
                    let extra = [&["-outfmt", outfmt][..], window].concat();
                    let what = format!(
                        "{} threads {threads} outfmt {outfmt} {window:?}",
                        query.display()
                    );
                    let (reference, _) = blastp(query, subject, threads, &[], &extra);
                    assert!(!reference.is_empty(), "{what}: the pair has a report");
                    let (walk, _) = blastp(query, subject, threads, &[UNGAP, WALK], &extra);
                    assert!(walk == reference, "{what}: LOSAT_X_BPWALK differs");
                    if outfmt == "0" {
                        continue;
                    }
                    // the walk of LOSAT_X_BPUNGAP alone
                    let (fast, _) = blastp(query, subject, threads, &[UNGAP], &extra);
                    assert!(fast == reference, "{what}: LOSAT_X_BPUNGAP differs");
                    // the switch alone does not turn the fast path on
                    if threads == 1 {
                        let (alone, stderr) = blastp(query, subject, threads, &[WALK], &extra);
                        assert!(alone == reference, "{what}: LOSAT_X_BPWALK alone differs");
                        assert!(!stderr.contains("[X_SHADOW]"), "{what}: {stderr}");
                    }
                }
            }
        }
    }
}

/// The count of the `[X_SHADOW] BPUNGAP compared=<n> walk=drop` line, which only the shadow
/// with the drop walks prints.
fn compared_with_drop_walk(stderr: &str) -> Option<u64> {
    stderr.lines().find_map(|line| {
        let rest = line.strip_prefix("[X_SHADOW] BPUNGAP compared=")?;
        let (count, rest) = rest.split_once(' ')?;
        rest.starts_with("walk=drop (ungapped extensions identical)")
            .then(|| count.parse::<u64>().ok())?
    })
}

#[test]
fn the_shadow_compares_the_drop_walks_with_the_reference() {
    let (query, subject) = generated_pair();
    let pairs = [
        (query, subject),
        (fixture("MeenMJNV.faa"), fixture("MejoMJNV.faa")),
    ];
    // the shadow compares every extension; one thread is enough
    let threads = 1;
    for (query, subject) in &pairs {
        for window in [&[][..], &["-window_size", "0"][..]] {
            let extra = [&["-outfmt", "6"][..], window].concat();
            let what = format!("{} {window:?}", query.display());
            let (reference, _) = blastp(query, subject, threads, &[], &extra);
            let (rows, stderr) = blastp(
                query,
                subject,
                threads,
                &[("LOSAT_X_BPUNGAPSHADOW", "1"), WALK],
                &extra,
            );
            assert!(rows == reference, "{what}: rows differ");
            let compared = compared_with_drop_walk(&stderr);
            assert!(
                compared.is_some_and(|n| n > 0),
                "{what}: no shadow line of the drop walks: {stderr}"
            );
            // without the switch the line does not name the drop walk
            let (rows, stderr) = blastp(
                query,
                subject,
                threads,
                &[("LOSAT_X_BPUNGAPSHADOW", "1")],
                &extra,
            );
            assert!(
                rows == reference,
                "{what}: rows differ (shadow without the switch)"
            );
            assert!(
                stderr
                    .lines()
                    .any(|line| line.starts_with("[X_SHADOW] BPUNGAP compared=")
                        && line.ends_with(" (ungapped extensions identical)")
                        && !line.contains("walk=drop")),
                "{what}: {stderr}"
            );
        }
    }
}

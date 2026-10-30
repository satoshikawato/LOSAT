//! BLASTX's compressed lookup table (`-word_size 5`, the table of `-task blastx-fast`) built flat
//! (`LOSAT_X_BPCLUT`, `src/algorithm/blastx/seed.rs` `x_compressed_lookup`) on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
//! ```c
//!     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//!                                     query, locations);
//!     s_CompressedLookupFinalize(lookup);
//! ```
//! The flat table lists each cell's hits in the order of the reference table's scan (the lib tests
//! of `x_clut_direct.rs`), so the property under test is that the report does not depend on the
//! switch, for subjects read from FASTA (`-subject`) and from a database (`-db`), on one thread
//! and on several. That reports equal NCBI's is the business of the BLASTX fixtures. The same for
//! the units searched at a time with this table (`LOSAT_X_BXCLUTCAP`, `blastx/x_bx_mem.rs`), which
//! is scheduling only, and the side of its rule each run takes.
//!
//! The query is long enough to be split into chunks, and a batch of short records follows it, so
//! that several tables are built in one search.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const REFERENCE_TABLE: (&str, &str) = ("LOSAT_X_BPCLUT", "0");
const NO_UNIT_RULE: (&str, &str) = ("LOSAT_X_BXCLUTCAP", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory that no other test of this process uses.
fn temp_dir(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!(
        "losat_blastx_fast_flat_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The binary with `env` as its only variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// names a file in `dir`: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(env: &[(&str, &str)], dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"))
        .envs(env.iter().copied());
    command
}

/// The search's stdout, after checking that it succeeded.
fn search(
    env: &[(&str, &str)],
    dir: &Path,
    query: &Path,
    target: &[&str],
    threads: &str,
) -> Vec<u8> {
    let output = losat(env, dir)
        .args([
            "blastx",
            "-word_size",
            "5",
            "-outfmt",
            "6",
            "-num_threads",
            threads,
            "-query",
        ])
        .arg(query)
        .args(target)
        .output()
        .expect("run LOSAT blastx");
    assert!(
        output.status.success(),
        "blastx: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

/// The query: a long record (split into chunks) and short ones (several batches).
fn write_query(dir: &Path) -> PathBuf {
    let genome = std::fs::read_to_string(fixture("AP027133.fasta")).unwrap();
    let residues: String = genome
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let mut text = format!(">long\n{}\n", &residues[..60_000]);
    for k in 0..16 {
        let at = 60_000 + 6_000 * k;
        text.push_str(&format!(">short{k}\n{}\n", &residues[at..at + 5_500]));
    }
    let query = dir.join("query.fna");
    std::fs::write(&query, text).unwrap();
    query
}

/// The database of `subject`, in `dir`.
fn make_db(dir: &Path, subject: &Path) -> PathBuf {
    let db = dir.join("db");
    let made = losat(&[], dir)
        .arg("makeblastdb")
        .arg("-in")
        .arg(subject)
        .args(["-dbtype", "prot", "-out"])
        .arg(&db)
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        made.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&made.stderr)
    );
    db
}

/// The sides of `LOSAT_X_BXCLUTCAP` a search recorded: (its own K, K raised for few subjects).
fn unit_rule_sides(dir: &Path, query: &Path, target: &[&str], threads: &str) -> (bool, bool) {
    let reach = dir.join(format!("reach-{threads}"));
    let _ = std::fs::remove_dir_all(&reach);
    let reach_text = reach.to_str().unwrap().to_string();
    search(
        &[("LOSAT_X_REACH", &reach_text)],
        dir,
        query,
        target,
        threads,
    );
    let mut sides = (false, false);
    for file in std::fs::read_dir(&reach).unwrap() {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields[0] == "LOSAT_X_BXCLUTCAP" {
                sides.0 |= fields[1] == "1";
                sides.1 |= fields[2] == "1";
            }
        }
    }
    sides
}

#[test]
fn the_flat_table_gives_the_reference_tables_report() {
    let dir = temp_dir("report");
    let subject = fixture("AP027133.faa");
    let db = make_db(&dir, &subject);
    let query = write_query(&dir);
    let subject = subject.to_str().unwrap();
    let db = db.to_str().unwrap();
    for target in [["-subject", subject], ["-db", db]] {
        for threads in ["1", "4"] {
            let flat = search(&[], &dir, &query, &target, threads);
            let reference = search(&[REFERENCE_TABLE], &dir, &query, &target, threads);
            assert!(!flat.is_empty(), "{target:?} n{threads}: no hits");
            assert_eq!(flat, reference, "{target:?} n{threads}");
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_units_at_a_time_do_not_change_the_report() {
    let dir = temp_dir("units");
    let subject = fixture("AP027133.faa");
    let db = make_db(&dir, &subject);
    let query = write_query(&dir);
    let subject = subject.to_str().unwrap();
    let db = db.to_str().unwrap();
    for target in [["-subject", subject], ["-db", db]] {
        let reference = search(&[NO_UNIT_RULE], &dir, &query, &target, "1");
        assert!(!reference.is_empty(), "{target:?}: no hits");
        for threads in ["3", "4", "8"] {
            for env in [
                &[][..],
                &[REFERENCE_TABLE][..],
                &[NO_UNIT_RULE][..],
                &[("LOSAT_X_BXCLUTCAP", "3")][..],
            ] {
                let got = search(env, &dir, &query, &target, threads);
                assert_eq!(got, reference, "{target:?} n{threads} {env:?}");
            }
        }
    }
    // Many subjects: the rule's own K. One subject: K rises to the threads.
    assert_eq!(
        unit_rule_sides(&dir, &query, &["-db", db], "4"),
        (true, false)
    );
    let one = dir.join("one.faa");
    let proteins = std::fs::read_to_string(&subject).unwrap();
    let second = proteins[1..].find('>').unwrap() + 1;
    std::fs::write(&one, &proteins[..second]).unwrap();
    assert_eq!(
        unit_rule_sides(&dir, &query, &["-subject", one.to_str().unwrap()], "4"),
        (false, true)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

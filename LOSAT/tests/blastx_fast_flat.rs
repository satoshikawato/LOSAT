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
//! and on several. That reports equal NCBI's is the business of the BLASTX fixtures.
//!
//! The query is long enough to be split into chunks, and a batch of short records follows it, so
//! that several tables are built in one search.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const REFERENCE_TABLE: (&str, &str) = ("LOSAT_X_BPCLUT", "0");

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

#[test]
fn the_flat_table_gives_the_reference_tables_report() {
    let dir = temp_dir("report");
    let subject = fixture("AP027133.faa");
    let db = dir.join("db");
    let made = losat(&[], &dir)
        .arg("makeblastdb")
        .arg("-in")
        .arg(&subject)
        .args(["-dbtype", "prot", "-out"])
        .arg(&db)
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        made.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&made.stderr)
    );
    let genome = std::fs::read_to_string(fixture("AP027133.fasta")).unwrap();
    let residues: String = genome
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let mut text = format!(">long\n{}\n", &residues[..60_000]);
    for k in 0..4 {
        let at = 60_000 + 3_000 * k;
        text.push_str(&format!(">short{k}\n{}\n", &residues[at..at + 2_500]));
    }
    let query = dir.join("query.fna");
    std::fs::write(&query, text).unwrap();

    let subject = subject.to_str().unwrap();
    let db = db.to_str().unwrap();
    let targets = [["-subject", subject], ["-db", db]];
    for target in &targets {
        for threads in ["1", "4"] {
            let flat = search(&[], &dir, &query, target, threads);
            let reference = search(&[REFERENCE_TABLE], &dir, &query, target, threads);
            assert!(!flat.is_empty(), "{target:?} n{threads}: no hits");
            assert_eq!(flat, reference, "{target:?} n{threads}");
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

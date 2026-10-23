//! BLASTX's memory at many threads on whole searches: the workers' diagonal tables kept one per
//! thread (`LOSAT_X_BXDIAGKEEP`), the lookup tables' neighbour lists shared and their vectors kept
//! between builds (`LOSAT_X_BXLUTPOOL`), and at most K query batches or chunks searched at a time
//! (`LOSAT_X_BXUNITCAP[=K]`; `src/algorithm/blastx/x_bx_mem.rs`,
//! `src/algorithm/tblastx/lookup/x_lut_direct.rs`).
//!
//! The three are in the strict set (on unless `=0` or `LOSAT_X_DEFAULT=0`). They change which
//! thread computes what, and the vectors the tables are built on, not what is computed, so the
//! property under test is that the report does not depend on them: standard output equals that of
//! the same search with all three `=0`, for several query batches (the batch flow,
//! `LOSAT_X_BXFLOW`) and for one query split into chunks (the chunk waves, `LOSAT_X_BXPAR`), with
//! the cap acting (`LOSAT_X_BXUNITCAP=2` at 4 threads; its built-in K = 8 needs 9 threads or more),
//! not acting, and with each part off. `LOSAT_X_REACH` shows that the cap acted. That a pooled
//! table equals the reference table is a unit test of `x_lut_direct.rs`; that the bytes equal
//! NCBI's is the business of the BLASTX fixtures and the campaign's comparison runs.
//!
//! The queries are pieces of AP027131; the subjects are the proteins of AP027131, AP027132 and
//! AP027133 (1,711 records, fewer than `LOSAT_X_BXCHUNKFLOW` needs, so a split query's chunks are
//! searched in waves).

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
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// A directory of this test process and test (nextest runs each test in a process of its own).
fn temp_dir(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-unit-cap-{test}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The letters of AP027131.
fn genome() -> String {
    std::fs::read_to_string(fixture("AP027131.fasta"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(str::chars)
        .filter(|c| !c.is_ascii_whitespace())
        .collect()
}

/// `pieces` queries of `length` letters of AP027131, 20,000 letters apart from position 100,000,
/// and the three protein sets (1,711 proteins), as files.
fn inputs(dir: &Path, pieces: usize, length: usize) -> (PathBuf, PathBuf) {
    let genome = genome();
    let mut text = String::new();
    for piece in 0..pieces {
        let start = 100_000 + piece * 20_000;
        text.push_str(&format!(
            ">piece{piece} AP027131 {start}-{}\n{}\n",
            start + length,
            &genome[start..start + length]
        ));
    }
    let query = dir.join("query.fna");
    std::fs::write(&query, text).unwrap();
    let set = dir.join("set.faa");
    let mut proteins = String::new();
    for name in ["AP027131.faa", "AP027132.faa", "AP027133.faa"] {
        proteins.push_str(&std::fs::read_to_string(fixture(name)).unwrap());
    }
    std::fs::write(&set, proteins).unwrap();
    (query, set)
}

/// The set's `losat makedb` database.
fn database(dir: &Path, set: &Path) -> String {
    let db = dir.join("set").to_str().unwrap().to_string();
    let output = losat()
        .arg("makedb")
        .arg("-in")
        .arg(set)
        .args(["-dbtype", "prot", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

/// One search at `threads` threads; returns standard output.
fn blastx(query: &Path, args: &[&str], threads: u32, env: &[(&str, &str)]) -> Vec<u8> {
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(query)
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{args:?} {env:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

/// The sides of `LOSAT_X_BXUNITCAP` one search reached (`LOSAT_X_REACH`): (own path, reference).
fn cap_sides(dir: &Path, query: &Path, args: &[&str], env: &[(&str, &str)]) -> (bool, bool) {
    let reach = dir.join(format!("reach-{}", env.len()));
    let _ = std::fs::remove_dir_all(&reach);
    let mut env = env.to_vec();
    let reach_text = reach.to_str().unwrap().to_string();
    env.push(("LOSAT_X_REACH", &reach_text));
    blastx(query, args, 4, &env);
    let mut sides = (false, false);
    for file in std::fs::read_dir(&reach).unwrap() {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields[0] == "LOSAT_X_BXUNITCAP" {
                sides.0 |= fields[1] == "1";
                sides.1 |= fields[2] == "1";
            }
        }
    }
    sides
}

const ALL_OFF: &[(&str, &str)] = &[
    ("LOSAT_X_BXDIAGKEEP", "0"),
    ("LOSAT_X_BXLUTPOOL", "0"),
    ("LOSAT_X_BXUNITCAP", "0"),
];

#[test]
fn several_query_batches_give_the_report_without_the_switches() {
    let dir = temp_dir("batches");
    // Eight queries of 4,000 letters: four query batches of 10,002 bases or more for the batch
    // flow at 4 threads.
    let (query, set) = inputs(&dir, 8, 4_000);
    let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
    let reference = blastx(&query, &args, 4, ALL_OFF);
    assert!(
        reference.iter().filter(|&&b| b == b'\n').count() > 20,
        "too few hits to test"
    );
    for env in [
        &[][..],
        &[("LOSAT_X_BXUNITCAP", "2")][..],
        &[("LOSAT_X_BXUNITCAP", "3"), ("LOSAT_X_BXLUTPOOL", "0")][..],
        &[("LOSAT_X_BXUNITCAP", "2"), ("LOSAT_X_BXDIAGKEEP", "0")][..],
        &[("LOSAT_X_DEFAULT", "0")][..],
    ] {
        assert!(
            blastx(&query, &args, 4, env) == reference,
            "{env:?}: report differs"
        );
    }
    // The cap acts with K below the thread count and only with the workers' tables kept.
    assert_eq!(
        cap_sides(&dir, &query, &args, &[("LOSAT_X_BXUNITCAP", "2")]),
        (true, false)
    );
    assert_eq!(cap_sides(&dir, &query, &args, &[]), (false, true));
    assert_eq!(
        cap_sides(
            &dir,
            &query,
            &args,
            &[("LOSAT_X_BXUNITCAP", "2"), ("LOSAT_X_BXDIAGKEEP", "0")]
        ),
        (false, true)
    );
}

#[test]
fn a_split_query_gives_the_report_without_the_switches() {
    let dir = temp_dir("split");
    // 45,000 letters: a query split into chunks (blastx/split.rs), searched in waves of two with
    // LOSAT_X_BXUNITCAP=2, each chunk's 1,711 subjects as tasks on the pool.
    let (query, set) = inputs(&dir, 1, 45_000);
    let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
    let reference = blastx(&query, &args, 4, ALL_OFF);
    assert!(
        reference.iter().filter(|&&b| b == b'\n').count() > 20,
        "too few hits to test"
    );
    for env in [
        &[][..],
        &[("LOSAT_X_BXUNITCAP", "2")][..],
        &[("LOSAT_X_BXUNITCAP", "2"), ("LOSAT_X_BXLUTPOOL", "0")][..],
    ] {
        assert!(
            blastx(&query, &args, 4, env) == reference,
            "{env:?}: report differs"
        );
    }
    assert_eq!(
        cap_sides(&dir, &query, &args, &[("LOSAT_X_BXUNITCAP", "2")]),
        (true, false)
    );
    // A database: the chunks' subject tasks read their subjects as blocks (LOSAT_X_BXDBSTREAM).
    let db = database(&dir, &set);
    let args = ["-db", db.as_str(), "-outfmt", "6"];
    assert!(
        blastx(&query, &args, 4, &[("LOSAT_X_BXUNITCAP", "2")])
            == blastx(&query, &args, 4, ALL_OFF),
        "-db: report differs"
    );
}

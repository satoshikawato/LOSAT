//! BLASTX's memory at many threads on whole searches: the workers' diagonal tables kept one per
//! thread (`LOSAT_X_BXDIAGKEEP`), the lookup tables' neighbour lists shared and their vectors kept
//! between builds (`LOSAT_X_BXLUTPOOL`; not in a search on one thread, `LOSAT_X_BXPOOLN1`), and at
//! most K query batches or chunks searched at a time
//! (`LOSAT_X_BXUNITCAP[=K]`; `src/algorithm/blastx/x_bx_mem.rs`,
//! `src/algorithm/tblastx/lookup/x_lut_direct.rs`).
//!
//! The four are in the strict set (on unless `=0`). They change which
//! thread computes what, and the vectors the tables are built on, not what is computed, so the
//! property under test is that the report does not depend on them: standard output equals that of
//! the same search with all four `=0`, for several query batches (the batch flow,
//! `LOSAT_X_BXFLOW`) and for one query split into chunks (the chunk waves, `LOSAT_X_BXPAR`), with
//! the cap acting (`LOSAT_X_BXUNITCAP=2` at 4 threads; its built-in K = 8 needs 9 threads or more),
//! not acting, and with each part off. `LOSAT_X_REACH` shows that the cap acted, and that it does
//! not act where K units would leave threads without a task (too few subjects). That a pooled
//! table equals the reference table is a unit test of `x_lut_direct.rs`; that the bytes equal
//! NCBI's is the business of the BLASTX fixtures and the campaign's comparison runs.
//!
//! The queries are pieces of AP027131; the subjects are its 595 proteins (more than the subject
//! flow's 64, fewer than `LOSAT_X_BXCHUNKFLOW` needs, so a split query's chunks are searched in
//! waves), or the first few of them.

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

/// `pieces` queries of `length` letters of AP027131 from position 100,000, 20,000 letters apart
/// (closer where that would pass the end), and the first `proteins` proteins of AP027131, as files.
fn inputs(dir: &Path, pieces: usize, length: usize, proteins: usize) -> (PathBuf, PathBuf) {
    let genome = genome();
    let step = ((genome.len() - 100_000 - length) / pieces).min(20_000);
    let mut text = String::new();
    for piece in 0..pieces {
        let start = 100_000 + piece * step;
        text.push_str(&format!(
            ">piece{piece} AP027131 {start}-{}\n{}\n",
            start + length,
            &genome[start..start + length]
        ));
    }
    let query = dir.join("query.fna");
    std::fs::write(&query, text).unwrap();
    let set = dir.join(format!("set{proteins}.faa"));
    let text = std::fs::read_to_string(fixture("AP027131.faa")).unwrap();
    let mut kept = String::new();
    let mut records = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            records += 1;
        }
        if records > proteins {
            break;
        }
        kept.push_str(line);
        kept.push('\n');
    }
    std::fs::write(&set, kept).unwrap();
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
    sides(dir, query, args, 4, env, "LOSAT_X_BXUNITCAP")
}

/// The sides of `switch` one search at `threads` threads reached: (own path, reference).
fn sides(
    dir: &Path,
    query: &Path,
    args: &[&str],
    threads: u32,
    env: &[(&str, &str)],
    switch: &str,
) -> (bool, bool) {
    let reach = dir.join(format!("reach-{}-{threads}", env.len()));
    let _ = std::fs::remove_dir_all(&reach);
    let mut env = env.to_vec();
    let reach_text = reach.to_str().unwrap().to_string();
    env.push(("LOSAT_X_REACH", &reach_text));
    blastx(query, args, threads, &env);
    let mut sides = (false, false);
    for file in std::fs::read_dir(&reach).unwrap() {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields[0] == switch {
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
    ("LOSAT_X_BXPOOLN1", "0"),
    ("LOSAT_X_BXUNITCAP", "0"),
];

#[test]
fn several_query_batches_give_the_report_without_the_switches() {
    let dir = temp_dir("batches");
    // Eight queries of 4,000 letters: four query batches of 10,002 bases or more for the batch
    // flow at 4 threads.
    let (query, set) = inputs(&dir, 8, 4_000, 595);
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
    // 40 proteins: the subject flow does not run (64 subjects or fewer), so a capped batch would
    // search on one thread; the cap does not act.
    let (query, few) = inputs(&dir, 8, 4_000, 40);
    let args = ["-subject", few.to_str().unwrap(), "-outfmt", "6"];
    assert_eq!(
        cap_sides(&dir, &query, &args, &[("LOSAT_X_BXUNITCAP", "2")]),
        (false, true)
    );
    // Forty queries: sixteen batches, more than the capped flow reads ahead of the writing with
    // K = 2 (`flow_bound`), so the reader waits and batches are spawned by searches that end.
    let (query, set) = inputs(&dir, 40, 4_000, 595);
    let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
    let reference = blastx(&query, &args, 4, ALL_OFF);
    assert!(
        blastx(&query, &args, 4, &[("LOSAT_X_BXUNITCAP", "2")]) == reference,
        "forty queries: report differs"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_split_query_gives_the_report_without_the_switches() {
    let dir = temp_dir("split");
    // 45,000 letters: a query split into chunks (blastx/split.rs), searched in waves of two with
    // LOSAT_X_BXUNITCAP=2, each chunk's 1,711 subjects as tasks on the pool.
    let (query, set) = inputs(&dir, 1, 45_000, 595);
    let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
    let reference = blastx(&query, &args, 4, ALL_OFF);
    assert!(
        reference.iter().filter(|&&b| b == b'\n').count() > 20,
        "too few hits to test"
    );
    for env in [&[][..], &[("LOSAT_X_BXUNITCAP", "2")][..]] {
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
    // 20 proteins make two tasks of 16 per chunk: two chunks at a time give four threads a task,
    // and the cap acts on the chunk waves (not on the batch flow, whose batch has no subject flow
    // with 64 subjects or fewer); 10 proteins make one task, and it acts nowhere.
    let (query, twenty) = inputs(&dir, 1, 45_000, 20);
    let args = ["-subject", twenty.to_str().unwrap(), "-outfmt", "6"];
    assert_eq!(
        cap_sides(&dir, &query, &args, &[("LOSAT_X_BXUNITCAP", "2")]),
        (true, true)
    );
    let (query, ten) = inputs(&dir, 1, 45_000, 10);
    let args = ["-subject", ten.to_str().unwrap(), "-outfmt", "6"];
    assert_eq!(
        cap_sides(&dir, &query, &args, &[("LOSAT_X_BXUNITCAP", "2")]),
        (false, true)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `LOSAT_X_BXKAONCE` (a query batch's Karlin blocks computed once, `blastx/search.rs`): the report
/// equals that without it. The reuse runs in the serial chunk loop for a batch that is not split,
/// so the search is run at one thread (the loop) and at four (the batch flow).
#[test]
fn several_queries_give_the_report_without_ka_once() {
    let dir = temp_dir("ka-once");
    let (query, set) = inputs(&dir, 4, 4_000, 595);
    let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
    let off = [("LOSAT_X_BXKAONCE", "0")];
    for threads in [1, 4] {
        let reference = blastx(&query, &args, threads, &off);
        assert!(
            reference.iter().filter(|&&b| b == b'\n').count() > 20,
            "too few hits to test"
        );
        assert!(
            blastx(&query, &args, threads, &[]) == reference,
            "{threads} threads: report differs"
        );
    }
    assert_eq!(
        sides(&dir, &query, &args, 1, &[], "LOSAT_X_BXKAONCE"),
        (true, false)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `LOSAT_X_BXPOOLN1` (a search on one thread builds and frees its lookup tables without the pool,
/// `tblastx/lookup/x_lut_direct.rs`): at one thread the report equals that without it, for several
/// query batches (one table each) and for a split query (a table per chunk); `LOSAT_X_REACH` shows
/// that it acts at one thread and not at four (where both build pooled).
#[test]
fn queries_give_the_report_without_pool_n1() {
    let dir = temp_dir("pool-n1");
    let off = [("LOSAT_X_BXPOOLN1", "0")];
    for (pieces, length) in [(4, 4_000), (1, 45_000)] {
        let (query, set) = inputs(&dir, pieces, length, 595);
        let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
        let reference = blastx(&query, &args, 1, &off);
        assert!(
            reference.iter().filter(|&&b| b == b'\n').count() > 20,
            "too few hits to test"
        );
        assert!(
            blastx(&query, &args, 1, &[]) == reference,
            "{pieces} x {length} letters: report differs"
        );
        for threads in [1, 4] {
            assert_eq!(
                sides(&dir, &query, &args, threads, &[], "LOSAT_X_BXPOOLN1"),
                (threads == 1, threads > 1),
                "{pieces} x {length} letters, {threads} threads: sides"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

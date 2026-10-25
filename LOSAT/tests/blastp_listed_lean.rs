//! blastp `-db` reports with `LOSAT_X_BPLISTEDLEAN`
//! (`src/algorithm/blastp/blast_engine/database/x_listed_lean.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless its variable is `0`.
//! With it on, the reports of a `-db` search hold the titles
//! of the reported subjects only and read the letters of the listed subjects from their stored
//! residues, instead of one record (ASCII letters and title) per listed subject. These tests
//! compare the reports with the switch off and on (one search of all queries and the block search
//! of `LOSAT_X_BPQBLOCK`; pairwise, tabular with `stitle` and the rendered rows, Seq-align text and
//! XML); the reach rows of the run with the switch on (`LOSAT_X_REACH`, `tests/reach_support`)
//! show that the lean path ran. The switch reads the titles search by search
//! (`x_read_titles_lean` once per search), the same for a block of the block search and for the
//! search of all queries, so three queries of the fixture are enough; the block search runs them
//! in two blocks (`BLOCKS`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1582
//! ```c++
//!     for (int index = 0; index < hit_list->hsplist_count; index++) {
//!         BlastHSPList* hsp_list = hit_list->hsplist_array[index];
//!         if (!hsp_list)
//!             continue;
//! ...
//!         const Uint4 kOid = hsp_list->oid;
//!         TSeqPos subj_length = 0;
//!         CRef<CSeq_id> subject_id;
//!         GetSequenceLengthAndId(seqinfo_src, kOid, subject_id, &subj_length);
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// What one run printed: the report and standard error.
struct Run {
    stdout: Vec<u8>,
    stderr: String,
}

/// A protein database of the fixture's subject proteome (`MejoMJNV.faa`, 124 sequences), made with
/// `losat makedb` in a directory of this test; returns the database name.
fn subject_database() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bplistedlean-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// Three records of the fixture's query proteome (`MeenMJNV.faa`, a nimavirus related to the
/// subjects') from BDT62583.1 on, written to `dir`: BDT62583.1 (629 letters), whose preliminary
/// list's subjects all leave the final one, so that it reports no hits; BDT62584.1 (944) and
/// BDT62585.1 (918), which report the same five subjects, read once by a search of both. A search
/// of all three lists 53 subjects and reports 5.
fn queries(dir: &Path) -> PathBuf {
    let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let start = text
        .find(">BDT62583.1 ")
        .expect("the record is in the fixture");
    let mut end = start;
    for _ in 0..3 {
        end = text[end + 1..]
            .find("\n>")
            .map_or(text.len(), |at| end + 1 + at + 1);
    }
    let path = dir.join("queries.faa");
    std::fs::write(&path, &text[start..end]).unwrap();
    path
}

/// One blastp run of `query` against the database `db`, with `vars` as its only `LOSAT_*`
/// variables that choose a code path (and its `BATCH_SIZE`), at most 5 subjects per query (so that
/// subjects of the preliminary lists leave the final ones) and `-outfmt outfmt`; it must exit with
/// 0. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn blastp(query: &Path, db: &Path, threads: u32, vars: &[(&str, &str)], outfmt: &str) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
            std::env::temp_dir().join("losat-bplistedlean-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-db")
        .arg(db)
        .args([
            "-evalue",
            "0.001",
            "-max_target_seqs",
            "5",
            "-outfmt",
            outfmt,
        ])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let run = Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    };
    assert!(
        output.status.success(),
        "{vars:?} threads {threads} outfmt {outfmt}: {:?} {}",
        output.status.code(),
        run.stderr
    );
    run
}

/// The one search of all queries.
const ALL: &[(&str, &str)] = &[("LOSAT_X_BPQBLOCK", "0")];
/// The block search in two blocks of one query batch each: BDT62583.1 and BDT62584.1 (the batch
/// ends at the query that reaches 900 letters), then BDT62585.1. A search runs in blocks only when
/// it has more than one; at the default `BATCH_SIZE` the three queries are one batch.
const BLOCKS: &[(&str, &str)] = &[("BATCH_SIZE", "900"), ("LOSAT_X_BPQBLOCK", "1")];

#[test]
fn lean_reports_are_the_records_reports() {
    let db = subject_database();
    let query = queries(db.parent().unwrap());
    // Each output format once, at one or two threads, in the search of all queries (`ALL`) or in
    // two blocks (`BLOCKS`).
    for (outfmt, threads, base) in [
        ("0", 2, ALL),
        ("7 qseqid sacc stitle qseq sseq btop ppos", 1, BLOCKS),
        (
            "6 qseqid sseqid stitle pident nident positive mismatch evalue bitscore",
            2,
            BLOCKS,
        ),
        ("8", 1, ALL),
        ("5", 2, BLOCKS),
    ] {
        let with = |extra: (&'static str, &'static str)| {
            base.iter().copied().chain([extra]).collect::<Vec<_>>()
        };
        let what = format!("outfmt {outfmt} threads {threads} {base:?}");
        let off = blastp(
            &query,
            &db,
            threads,
            &with(("LOSAT_X_BPLISTEDLEAN", "0")),
            outfmt,
        );
        assert!(!off.stdout.is_empty(), "{what}: the queries have rows");
        let reach = reach_support::Reach::new("listed_lean");
        let mut on_vars: Vec<(&str, &str)> = with(("LOSAT_X_BPLISTEDLEAN", "1"));
        on_vars.push(reach.var());
        let on = blastp(&query, &db, threads, &on_vars, outfmt);
        assert!(on.stdout == off.stdout, "{what}: report");
        assert_eq!(on.stderr, off.stderr, "{what}: standard error");
        let sides = reach.sides("LOSAT_X_BPLISTEDLEAN");
        assert!(sides.0, "{what}: the lean path did not run: {sides:?}");
        if base == BLOCKS {
            // Outfmt 5 and 7 have rules of their own for blocks: this case runs in two blocks, a
            // `[TIMING] blastp_total` line each.
            let timing = blastp(&query, &db, threads, &with(("LOSAT_TIMING", "1")), outfmt);
            let blocks = timing
                .stderr
                .lines()
                .filter(|line| line.starts_with("[TIMING] blastp_total:"))
                .count();
            assert_eq!(blocks, 2, "{what}: {}", timing.stderr);
        }
    }
}

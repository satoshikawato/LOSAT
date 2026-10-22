// CI reference step: skipped (each LOSAT run here sets its own LOSAT_* variables; tests/reference_step_targets.sh).
//! blastp `-db` reports with `LOSAT_X_BPLISTEDLEAN`
//! (`src/algorithm/blastp/blast_engine/database/x_listed_lean.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless its variable is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. With it on, the reports of a `-db` search hold the titles
//! of the reported subjects only and read the letters of the listed subjects from their stored
//! residues, instead of one record (ASCII letters and title) per listed subject. These tests
//! compare the reports with the switch off and on (one search of all queries and the block search
//! of `LOSAT_X_BPQBLOCK`; pairwise, tabular with `stitle` and the rendered rows, Seq-align text and
//! XML)
//! and run its shadow mode (`LOSAT_X_BPLISTEDLEANSHADOW`), which makes the records too and compares
//! them with the lean titles rank by rank; its line at exit shows that the path ran. A difference
//! stops the run.
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

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

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

const SHADOW_PREFIX: &str = "[X_SHADOW] BPLISTEDLEAN ";

impl Run {
    /// The counters of the `[X_SHADOW] BPLISTEDLEAN` line, which must end with `identical`.
    fn shadow(&self) -> HashMap<String, u64> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with(SHADOW_PREFIX))
            .unwrap_or_else(|| panic!("no {SHADOW_PREFIX}line in {}", self.stderr));
        assert!(line.ends_with(" identical"), "{line}");
        line[SHADOW_PREFIX.len()..]
            .split(' ')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.parse().expect("a count")))
            .collect()
    }

    /// Standard error without the shadow's line.
    fn stderr_without_shadow(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with(SHADOW_PREFIX))
            .map(|line| format!("{line}\n"))
            .collect()
    }
}

/// A protein database of the fixture's subject proteome (`MejoMJNV.faa`, 124 sequences), made with
/// `LOSAT makedb` in a directory of this test; returns the database name.
fn subject_database() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bplistedlean-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let name = dir.join("MejoMJNV");
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
        .arg("makedb")
        .arg("-in")
        .arg(fixture("MejoMJNV.faa"))
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// One blastp run of the fixture pair (two related nimavirus proteomes) against the database `db`,
/// with `vars` as its only `LOSAT_*` variables that choose a code path, at most 5 subjects per
/// query (so that subjects of the preliminary lists leave the final ones) and `-outfmt outfmt`;
/// it must exit with 0. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn blastp(db: &Path, threads: u32, vars: &[(&str, &str)], outfmt: &str) -> Run {
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
            std::env::temp_dir().join("losat-bplistedlean-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(fixture("MeenMJNV.faa"))
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

#[test]
fn lean_reports_are_the_records_reports() {
    let db = subject_database();
    // Each output format once; one and two threads, one search of all queries and the block
    // search of `LOSAT_X_BPQBLOCK` (its default) alternate.
    for (outfmt, threads, block) in [
        ("0", 2, Some("0")),
        ("7 qseqid sacc stitle qseq sseq btop ppos", 1, None),
        (
            "6 qseqid sseqid stitle pident nident positive mismatch evalue bitscore",
            2,
            None,
        ),
        ("8", 1, Some("0")),
        ("5", 2, None),
    ] {
        let base: Vec<(&str, &str)> = block.map(|b| ("LOSAT_X_BPQBLOCK", b)).into_iter().collect();
        let with = |extra: (&'static str, &'static str)| {
            base.iter().copied().chain([extra]).collect::<Vec<_>>()
        };
        let what = format!("outfmt {outfmt} threads {threads} block {block:?}");
        let off = blastp(&db, threads, &with(("LOSAT_X_BPLISTEDLEAN", "0")), outfmt);
        assert!(!off.stdout.is_empty(), "{what}: the fixture pair has rows");
        let on = blastp(&db, threads, &with(("LOSAT_X_BPLISTEDLEAN", "1")), outfmt);
        assert!(on.stdout == off.stdout, "{what}: report");
        assert_eq!(on.stderr, off.stderr, "{what}: standard error");
        let shadow = blastp(
            &db,
            threads,
            &with(("LOSAT_X_BPLISTEDLEANSHADOW", "1")),
            outfmt,
        );
        assert!(shadow.stdout == off.stdout, "{what}: shadow report");
        assert_eq!(shadow.stderr_without_shadow(), off.stderr, "{what}");
        let counts = shadow.shadow();
        assert!(counts["searches"] > 0, "{what}: {counts:?}");
        assert!(
            counts["ranks"] > counts["reported"] && counts["reported"] > 0,
            "{what}: subjects of the preliminary lists leave the final ones: {counts:?}"
        );
        assert!(counts["letters"] > 0, "{what}: {counts:?}");
    }
}

//! BLASTP's `-db` window of `LOSAT_X_BPDBWIDE`
//! (`src/algorithm/blastp/blast_engine/database/x_db_keep.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. It acts where a
//! `-db` search runs its preliminary stage on the pool (more than one thread): the window of the
//! blocks is wider, and narrows while many bytes of results wait for the merge. The width decides
//! only when a block is taken, so these tests compare the reports with the switch on and off with
//! the report at one thread, for a database of more blocks than the narrow window holds, and read
//! the sides of the decision (`LOSAT_X_REACH`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413,1553-1554
//! ```c
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ...
//!          /* Save the results. */
//!          status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const WIDE_OFF: (&str, &str) = ("LOSAT_X_BPDBWIDE", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// What one run printed: the report and standard error.
struct Run {
    code: Option<i32>,
    stdout: Vec<u8>,
    stderr: String,
}

/// A protein database of `copies` numbered copies of the fixture's proteome (`MejoMJNV.faa`, about
/// 90 KB of sequence), with 40 copies about 3.6 MB: about 38 blocks of 96 KiB, more than the 8 blocks
/// per thread of the narrow window at 2 and 4 threads, and the query file `query.faa` beside it, the first 10 records
/// of the fixture's other proteome. Returns the database name.
fn database(label: &str, copies: usize) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpdbwide-{label}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let query = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let end = query
        .match_indices('>')
        .nth(10)
        .map_or(query.len(), |(at, _)| at);
    std::fs::write(dir.join("query.faa"), &query[..end]).unwrap();
    let text = std::fs::read_to_string(fixture("MejoMJNV.faa")).unwrap();
    let mut records = String::with_capacity(copies * text.len());
    for copy in 0..copies {
        let mut record = 0;
        for line in text.lines() {
            if line.starts_with('>') {
                records.push_str(&format!(">c{copy}_{record}\n"));
                record += 1;
            } else {
                records.push_str(line);
                records.push('\n');
            }
        }
    }
    let fasta = dir.join("copies.faa");
    std::fs::write(&fasta, records).unwrap();
    let name = dir.join("copies");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(&fasta)
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

/// One blastp run of the query file beside `db` against it, with `vars` as its only `LOSAT_*`
/// variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn blastp(task: &str, db: &Path, threads: u32, vars: &[(&str, &str)]) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpdbwide-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .args(["-task", task])
        .arg("-query")
        .arg(db.with_file_name("query.faa"))
        .arg("-db")
        .arg(db)
        .args(["-evalue", "0.001", "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let run = Run {
        code: output.status.code(),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    };
    assert!(
        run.code == Some(0),
        "{task} {vars:?} threads {threads}: {:?} {}",
        run.code,
        run.stderr
    );
    run
}

#[test]
fn reports_with_the_wide_window_are_the_reports_with_the_narrow_one() {
    let db = database("reports", 40);
    for (task, thread_counts) in [("blastp", &[4u32][..]), ("blastp-fast", &[2, 4][..])] {
        let one_thread = blastp(task, &db, 1, &[]);
        assert!(
            !one_thread.stdout.is_empty(),
            "{task}: the search found nothing"
        );
        for &threads in thread_counts {
            let what = format!("{task} threads {threads}");
            let off = blastp(task, &db, threads, &[WIDE_OFF]);
            assert!(off.stdout == one_thread.stdout, "{what}: off differs");
            let reach = reach_support::Reach::new("bpdbwide");
            let on = blastp(task, &db, threads, &[reach.var()]);
            assert!(on.stdout == one_thread.stdout, "{what}: on differs");
            assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
            assert_eq!(reach.sides("LOSAT_X_BPDBWIDE"), (true, false), "{what}");
        }
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// At one thread the blocks are searched one after another, and a database of fewer blocks than the
// narrow window holds behaves as without the switch: the switch's decision takes its reference
// side only.
#[test]
fn the_window_does_not_act_at_one_thread_or_on_few_blocks() {
    let db = database("one", 1);
    for threads in [1u32, 2] {
        let reach = reach_support::Reach::new("bpdbwide-one");
        let run = blastp("blastp-fast", &db, threads, &[reach.var()]);
        assert_eq!(
            reach.sides("LOSAT_X_BPDBWIDE"),
            (false, true),
            "threads {threads}: {}",
            run.stderr
        );
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

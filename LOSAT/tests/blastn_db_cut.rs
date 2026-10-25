//! BLASTN's `-db` blocks cut by the database size and the thread count with `LOSAT_X_BNDBCUT`
//! (`src/algorithm/blastn/blast_engine/x_db_cut.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`, and it acts only
//! where a `-db` search runs its preliminary stage in rounds on the pool (more than one thread).
//! With it on, the blocks of the preliminary stage hold about a quarter of the database's sequence
//! bytes per thread (at least 4 KiB, at most the 1 MiB of the search without it), and the results
//! are merged in ordinal order. These tests compare the reports with the switch off and on
//! (megablast, blastn and dc-megablast, tabular and pairwise) with the report at one thread, and
//! read the sides of the switch's decision (`LOSAT_X_REACH`) to show that the cut ran where it
//! should: blocks other than the 1 MiB of the search without it at 2 threads and more, and not at
//! one thread.
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

use LOSAT::seqdb::{SeqDb, SeqType};

mod reach_support;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of the one record of the FASTA fixture `name`.
fn genome(name: &str) -> Vec<u8> {
    let text = std::fs::read(fixture(name)).unwrap();
    text.split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// What one run printed: the report and standard error.
struct Run {
    code: Option<i32>,
    stdout: Vec<u8>,
    stderr: String,
}

/// A directory of the test (`label`, which no other test uses), empty.
fn temp_dir(label: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bndbcut-{label}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Pieces of MG1655 of 1,000 bases every 5,000 bases (928 sequences, about 233 KB of sequence
/// file: one block of 1 MiB, 8 blocks of the cut at 2 threads, 16 at 4), `count` of them from
/// piece `first` on; piece 7 has a run of ten `N` and piece 300 an `R` and a `Y` (ambiguity data
/// in the sequence file).
fn pieces(first: usize, count: usize) -> Vec<u8> {
    let letters = genome("MG1655.fna");
    let mut fasta = Vec::new();
    for (index, start) in (0..letters.len() - 1_000)
        .step_by(5_000)
        .enumerate()
        .skip(first)
        .take(count)
    {
        let mut piece = letters[start..start + 1_000].to_vec();
        if index == 7 {
            piece[100..110].copy_from_slice(b"NNNNNNNNNN");
        }
        if index == 300 {
            piece[40] = b'R';
            piece[700] = b'Y';
        }
        fasta.extend(format!(">piece{index} MG1655 {}-{}\n", start + 1, start + 1_000).bytes());
        for line in piece.chunks(70) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
    }
    fasta
}

/// A nucleotide database of `pieces(first, count)` made with `losat makedb` in the directory
/// `label`, with `extra` arguments (`-max_file_sz` for volumes); returns the database name.
fn database(label: &str, first: usize, count: usize, extra: &[&str]) -> PathBuf {
    let dir = temp_dir(label);
    let input = dir.join("pieces.fna");
    std::fs::write(&input, pieces(first, count)).unwrap();
    let name = dir.join("pieces");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(&input)
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&name)
        .args(extra)
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// The queries: MG1655 223,001-226,000 (rrnH: hits in the pieces of the other rRNA operons),
/// MG1655 1,000,001-1,010,000 (hits across two pieces) and Sakai 1,000,001-1,020,000 (a related
/// strain), written to the directory of `db`.
fn queries(db: &Path) -> PathBuf {
    let (k12, sakai) = (genome("MG1655.fna"), genome("Sakai.fna"));
    let mut fasta = Vec::new();
    for (name, letters) in [
        ("rrnH", &k12[223_000..226_000]),
        ("k12_1M", &k12[1_000_000..1_010_000]),
        ("sakai_1M", &sakai[1_000_000..1_020_000]),
    ] {
        fasta.extend(format!(">{name}\n").bytes());
        for line in letters.chunks(70) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
    }
    let path = db.with_file_name("queries.fna");
    std::fs::write(&path, fasta).unwrap();
    path
}

/// One blastn run of `queries(db)` with `-task task` against the database `db`, with `vars` as
/// its only `LOSAT_*` variables, and `-outfmt outfmt`. `LOSAT_BLASTX_WORKER_LOG` stays: a binary
/// built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that
/// file name as soon as a search has more than one thread, and it selects nothing; without the
/// caller's value a file of the tests is named. The run is returned whatever its exit code.
fn blastn_any_exit(
    task: &str,
    db: &Path,
    threads: u32,
    vars: &[(&str, &str)],
    outfmt: &str,
) -> Run {
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
            std::env::temp_dir().join("losat-bndbcut-worker.log"),
        );
    }
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(queries(db))
        .arg("-db")
        .arg(db)
        .args(["-evalue", "1e-5", "-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    Run {
        code: output.status.code(),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// `blastn_any_exit`, which must exit with 0.
fn blastn(task: &str, db: &Path, threads: u32, vars: &[(&str, &str)], outfmt: &str) -> Run {
    let run = blastn_any_exit(task, db, threads, vars, outfmt);
    assert!(
        run.code == Some(0),
        "{task} {vars:?} threads {threads}: {:?} {}",
        run.code,
        run.stderr
    );
    run
}

const OFF: (&str, &str) = ("LOSAT_X_BNDBCUT", "0");

/// The switch off and on (the default) give the same report and standard error, the report of
/// `one_thread` (the same search at one thread); returns the sides of the switch's decision in the
/// default run, `(new, ref)`: `new` when the search cut its blocks.
fn same_reports(
    task: &str,
    db: &Path,
    threads: u32,
    outfmt: &str,
    one_thread: &Run,
) -> (bool, bool) {
    let what = format!("{task} threads {threads} outfmt {outfmt}");
    assert!(
        !one_thread.stdout.is_empty(),
        "{what}: the queries have hits"
    );
    let off = blastn(task, db, threads, &[OFF], outfmt);
    assert!(
        off.stdout == one_thread.stdout,
        "{what}: the report differs from one thread's without the switch"
    );
    let reach = reach_support::Reach::new("bndbcut");
    let on = blastn(task, db, threads, &[reach.var()], outfmt);
    assert!(
        on.stdout == one_thread.stdout,
        "{what}: the report differs with the switch"
    );
    assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
    reach.sides("LOSAT_X_BNDBCUT")
}

/// The blocks of the cut of the database `db` at `threads` threads: its 2-bit bytes and one byte
/// per sequence over 4 x threads, at least 4 KiB, as `SeqDb::next_block` cuts the sequences.
fn blocks_of_the_cut(db: &Path, threads: u64) -> u64 {
    let seqdb = SeqDb::open(db.to_str().unwrap(), SeqType::Nucleotide).unwrap();
    let bytes = seqdb.total_length() / 4 + u64::from(seqdb.num_seqs());
    let cut = (bytes / (4 * threads)).clamp(4096, 1 << 20);
    let (mut blocks, mut next) = (0, 0);
    while let Some(range) = seqdb.next_block(next, cut) {
        blocks += 1;
        next = range.end_ordinal();
    }
    blocks
}

// megablast, blastn and dc-megablast at 2 and 4 threads, tabular; megablast pairwise at 4.
#[test]
fn reports_with_the_cut_are_the_reports_without_it() {
    let db = database("tasks", 0, 10_000, &[]);
    for threads in [2u64, 4] {
        let blocks = blocks_of_the_cut(&db, threads);
        assert!(blocks >= 2 * threads, "{threads} threads: {blocks} blocks");
    }
    for task in ["megablast", "blastn", "dc-megablast"] {
        let one_thread = blastn(task, &db, 1, &[], "6");
        for threads in [2u64, 4] {
            let sides = same_reports(task, &db, threads as u32, "6", &one_thread);
            assert!(sides.0, "{task} threads {threads}: sides {sides:?}");
        }
    }
    let one_thread = blastn("megablast", &db, 1, &[], "0");
    let sides = same_reports("megablast", &db, 4, "0", &one_thread);
    assert!(sides.0, "pairwise report: sides {sides:?}");
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// At one thread the blocks keep 1 MiB, and a database below 4 KiB x 4 x threads of sequence bytes
// is one block of the cut too.
#[test]
fn the_cut_does_not_act_at_one_thread_or_on_a_database_of_one_block() {
    let db = database("one", 0, 10_000, &[]);
    let reach = reach_support::Reach::new("bndbcut-one");
    let run = blastn("megablast", &db, 1, &[reach.var()], "6");
    let sides = reach.sides("LOSAT_X_BNDBCUT");
    assert_eq!(sides, (false, true), "one thread: the blocks keep 1 MiB");
    assert!(
        run.stdout == blastn("megablast", &db, 1, &[OFF], "6").stdout,
        "the report differs"
    );
    let _ = std::fs::remove_dir_all(db.parent().unwrap());

    // Pieces 196-207 (about 3 KB of sequence file, with the hits of MG1655 1,000,001-1,010,000 in
    // pieces 200 and 201) are below the 4 KiB floor: one block at every thread count.
    let small = database("small", 196, 12, &[]);
    assert_eq!(blocks_of_the_cut(&small, 4), 1);
    let one_thread = blastn("megablast", &small, 1, &[], "6");
    let sides = same_reports("megablast", &small, 4, "6", &one_thread);
    assert!(sides.0, "one block of the 4 KiB floor: sides {sides:?}");
    let _ = std::fs::remove_dir_all(small.parent().unwrap());
}

// A damaged database (its `.nsq` cut short by 500 bytes) through the command line: the search
// fails with the error of the one-thread order, so standard error and the exit code are the same at
// one thread and, with the switch off and on, at 4 threads.
#[test]
fn a_damaged_database_fails_with_the_same_error_with_the_cut() {
    let db = database("damaged", 0, 10_000, &[]);
    let nsq = db.with_extension("nsq");
    let bytes = std::fs::read(&nsq).unwrap();
    std::fs::write(&nsq, &bytes[..bytes.len() - 500]).unwrap();
    let one_thread = blastn_any_exit("megablast", &db, 1, &[], "6");
    assert!(
        one_thread.code.is_some_and(|code| code != 0),
        "the damaged database must fail: {:?} {}",
        one_thread.code,
        one_thread.stderr
    );
    assert!(!one_thread.stderr.is_empty());
    for (what, vars) in [("off", vec![OFF]), ("on", vec![])] {
        let run = blastn_any_exit("megablast", &db, 4, &vars, "6");
        assert_eq!(run.code, one_thread.code, "{what}: exit code differs");
        assert_eq!(
            run.stderr, one_thread.stderr,
            "{what}: standard error differs from one thread's"
        );
        assert!(run.stdout == one_thread.stdout, "{what}: report differs");
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// A database of several volumes (`-max_file_sz 60000`): the cut ends blocks at volume ends, and 4
// threads with the switch on give one thread's report.
#[test]
fn a_database_of_volumes_reports_with_the_cut_are_the_reports_without_it() {
    let db = database("volumes", 0, 10_000, &["-max_file_sz", "60000"]);
    let volumes = (0..100)
        .take_while(|volume| {
            db.with_file_name(format!("pieces.{volume:02}.nsq"))
                .exists()
        })
        .count();
    assert!(volumes >= 3, "{volumes} volumes");
    let one_thread = blastn("megablast", &db, 1, &[], "6");
    assert!(!one_thread.stdout.is_empty());
    let reach = reach_support::Reach::new("bndbcut-volumes");
    let on = blastn("megablast", &db, 4, &[reach.var()], "6");
    assert!(
        on.stdout == one_thread.stdout,
        "the report differs with the switch"
    );
    let sides = reach.sides("LOSAT_X_BNDBCUT");
    assert!(sides.0, "volumes: sides {sides:?}");
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

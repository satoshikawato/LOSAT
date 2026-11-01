//! BLASTN's `-db` preliminary stage without rounds with `LOSAT_X_BNDBSTREAM`
//! (`src/algorithm/blastn/blast_engine/x_db_stream.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`, and it acts only
//! where a `-db` search runs its preliminary stage on the pool (more than one thread). With it on,
//! the threads claim the blocks in ordinal order and one of them hands their results to the hit
//! lists in that order while the others search the next blocks. These tests compare the reports
//! with the switch off (rounds) and on (megablast, blastn and dc-megablast, tabular and pairwise,
//! with the blocks of LOSAT_X_BNDBCUT and of 1 MiB) with the report at one thread, and read the
//! sides of the switch's decision (`LOSAT_X_REACH`) to show that the window ran at 2 threads and
//! more and not at one.
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
    let dir = std::env::temp_dir().join(format!("losat-bndbstream-{label}-{}", std::process::id()));
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

/// A nucleotide database of `pieces(first, count)` made with `LOSAT makeblastdb` in the directory
/// `label`, with `extra` arguments (`-max_file_sz` for volumes); returns the database name.
fn database(label: &str, first: usize, count: usize, extra: &[&str]) -> PathBuf {
    let dir = temp_dir(label);
    let input = dir.join("pieces.fna");
    std::fs::write(&input, pieces(first, count)).unwrap();
    let name = dir.join("pieces");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(&input)
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&name)
        .args(extra)
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
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
            std::env::temp_dir().join("losat-bndbstream-worker.log"),
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

const OFF: (&str, &str) = ("LOSAT_X_BNDBSTREAM", "0");
const CUT_OFF: (&str, &str) = ("LOSAT_X_BNDBCUT", "0");

/// The switch off (rounds) and on (the default), each with `extra`, give the same report and
/// standard error, the report of `one_thread` (the same search at one thread); returns the sides
/// of the switch's decision in the default run, `(new, ref)`: `new` when the window ran.
fn same_reports(
    task: &str,
    db: &Path,
    threads: u32,
    outfmt: &str,
    extra: &[(&'static str, &'static str)],
    one_thread: &Run,
) -> (bool, bool) {
    let what = format!("{task} threads {threads} outfmt {outfmt} {extra:?}");
    assert!(
        !one_thread.stdout.is_empty(),
        "{what}: the queries have hits"
    );
    let off_vars: Vec<(&str, &str)> = extra.iter().copied().chain([OFF]).collect();
    let off = blastn(task, db, threads, &off_vars, outfmt);
    assert!(
        off.stdout == one_thread.stdout,
        "{what}: the report differs from one thread's without the switch"
    );
    let reach = reach_support::Reach::new("bndbstream");
    let mut on_vars: Vec<(&str, &str)> = extra.to_vec();
    on_vars.push(reach.var());
    let on = blastn(task, db, threads, &on_vars, outfmt);
    assert!(
        on.stdout == one_thread.stdout,
        "{what}: the report differs with the switch"
    );
    assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
    reach.sides("LOSAT_X_BNDBSTREAM")
}

// megablast, blastn and dc-megablast at 2, 4 and 8 threads, tabular, with the blocks of the cut;
// megablast at 4 threads with blocks of 1 MiB (the database is one block: a window of one block);
// megablast pairwise at 4 threads.
#[test]
fn reports_with_the_window_are_the_reports_with_rounds() {
    let db = database("tasks", 0, 10_000, &[]);
    for task in ["megablast", "blastn", "dc-megablast"] {
        let one_thread = blastn(task, &db, 1, &[], "6");
        for threads in [2u64, 4, 8] {
            let sides = same_reports(task, &db, threads as u32, "6", &[], &one_thread);
            assert!(sides.0, "{task} threads {threads}: sides {sides:?}");
        }
    }
    let one_thread = blastn("megablast", &db, 1, &[], "6");
    let sides = same_reports("megablast", &db, 4, "6", &[CUT_OFF], &one_thread);
    assert!(sides.0, "blocks of 1 MiB: sides {sides:?}");
    let one_thread = blastn("megablast", &db, 1, &[], "0");
    let sides = same_reports("megablast", &db, 4, "0", &[], &one_thread);
    assert!(sides.0, "pairwise report: sides {sides:?}");
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// `LOSAT_X_BNDBWIDE` (the window of 16 blocks per thread while few bytes of results wait) changes no
// report: on and off at 2, 4 and 8 threads, for each task. Its own side runs where the list of
// blocks is longer than 2 blocks per thread (4 per thread of the cut), its reference side where
// the database is one block (blocks of 1 MiB).
#[test]
fn reports_with_the_wide_window_are_the_reports_with_the_narrow_one() {
    const WIDE_OFF: (&str, &str) = ("LOSAT_X_BNDBWIDE", "0");
    let db = database("wide", 0, 10_000, &[]);
    for task in ["megablast", "blastn", "dc-megablast"] {
        let one_thread = blastn(task, &db, 1, &[], "6");
        for threads in [2u32, 4, 8] {
            let what = format!("{task} threads {threads}");
            let off = blastn(task, &db, threads, &[WIDE_OFF], "6");
            assert!(off.stdout == one_thread.stdout, "{what}: off differs");
            let reach = reach_support::Reach::new("bndbwide");
            let on = blastn(task, &db, threads, &[reach.var()], "6");
            assert!(on.stdout == one_thread.stdout, "{what}: on differs");
            assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
            assert_eq!(reach.sides("LOSAT_X_BNDBWIDE"), (true, false), "{what}");
        }
    }
    let reach = reach_support::Reach::new("bndbwide-one-block");
    let on = blastn("megablast", &db, 4, &[CUT_OFF, reach.var()], "6");
    assert_eq!(
        reach.sides("LOSAT_X_BNDBWIDE"),
        (false, true),
        "{}",
        on.stderr
    );
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// `LOSAT_X_BNDBSMALL` (a largest block of 256 KiB) changes no report. Pieces of two genomes, about
// 2.5 MiB of sequence file: above 256 KiB per block at 2 threads (its own side), below it at 4 (its
// reference side), and at one thread the blocks are those of 1 MiB.
#[test]
fn reports_with_small_blocks_are_the_reports_with_large_ones() {
    const SMALL_OFF: (&str, &str) = ("LOSAT_X_BNDBSMALL", "0");
    let dir = temp_dir("small");
    let mut fasta = Vec::new();
    for (name, letters) in [
        ("k12", genome("MG1655.fna")),
        ("sakai", genome("Sakai.fna")),
    ] {
        for (index, piece) in letters.chunks(5_000).enumerate() {
            fasta.extend(format!(">{name}{index}\n").bytes());
            fasta.extend_from_slice(piece);
            fasta.push(b'\n');
        }
    }
    let input = dir.join("two.fna");
    std::fs::write(&input, fasta).unwrap();
    let db = dir.join("two");
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(&input)
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&db)
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(output.status.success(), "makeblastdb failed");
    let one_thread = blastn("megablast", &db, 1, &[], "6");
    assert!(!one_thread.stdout.is_empty());
    for (threads, sides) in [(2u32, (true, false)), (4, (false, true))] {
        let what = format!("threads {threads}");
        let off = blastn("megablast", &db, threads, &[SMALL_OFF], "6");
        assert!(off.stdout == one_thread.stdout, "{what}: off differs");
        let reach = reach_support::Reach::new("bndbsmall");
        let on = blastn("megablast", &db, threads, &[reach.var()], "6");
        assert!(on.stdout == one_thread.stdout, "{what}: on differs");
        assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
        assert_eq!(reach.sides("LOSAT_X_BNDBSMALL"), sides, "{what}");
    }
    let _ = std::fs::remove_dir_all(dir);
}

// At one thread the blocks are searched one after another: the switch's decision takes its
// reference side only.
#[test]
fn the_window_does_not_act_at_one_thread() {
    let db = database("one", 196, 12, &[]);
    let reach = reach_support::Reach::new("bndbstream-one");
    let run = blastn("megablast", &db, 1, &[reach.var()], "6");
    assert_eq!(
        reach.sides("LOSAT_X_BNDBSTREAM"),
        (false, true),
        "{}",
        run.stderr
    );
    assert!(
        run.stdout == blastn("megablast", &db, 1, &[OFF], "6").stdout,
        "the report differs"
    );
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// A damaged database (its `.nsq` cut short by 500 bytes) through the command line: the search
// fails with the error of the one-thread order, so standard error and the exit code are the same at
// one thread and, with the switch off and on, at 4 threads.
#[test]
fn a_damaged_database_fails_with_the_same_error_with_the_window() {
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

// A database of several volumes (`-max_file_sz 60000`): at 4 threads the switch off and on give
// one thread's report, with the blocks of the cut and of 1 MiB.
#[test]
fn a_database_of_volumes_reports_with_the_window_are_the_reports_with_rounds() {
    let db = database("volumes", 0, 10_000, &["-max_file_sz", "60000"]);
    let volumes = (0..100)
        .take_while(|volume| {
            db.with_file_name(format!("pieces.{volume:02}.nsq"))
                .exists()
        })
        .count();
    assert!(volumes >= 3, "{volumes} volumes");
    let one_thread = blastn("megablast", &db, 1, &[], "6");
    // The blocks of the cut and, without it, of 1 MiB: a volume ends a block, so the window has a
    // block per volume at least.
    for extra in [&[][..], &[CUT_OFF][..]] {
        let sides = same_reports("megablast", &db, 4, "6", extra, &one_thread);
        assert!(sides.0, "{extra:?}: sides {sides:?}");
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

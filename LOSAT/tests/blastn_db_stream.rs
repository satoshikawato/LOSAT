//! BLASTN's `-db` preliminary stage without rounds with `LOSAT_X_BNDBSTREAM`
//! (`src/algorithm/blastn/blast_engine/x_db_stream.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`, and it acts only
//! where a `-db` search runs its preliminary stage on the pool (more than one thread). With it on,
//! the threads claim the blocks in ordinal order and one of them hands their results to the hit
//! lists in that order while the others search the next blocks. These tests compare the reports
//! with the switch off (rounds) and on (megablast, blastn and dc-megablast, tabular and pairwise,
//! with the blocks of LOSAT_X_BNDBCUT and of 1 MiB) with the report at one thread, and run its
//! shadow mode (`LOSAT_X_BNDBSTREAMSHADOW`), which searches the database again in rounds and
//! compares every result handed to the hit lists; its line at exit shows that the path ran. A
//! difference stops the run.
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

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use LOSAT::seqdb::{SeqDb, SeqType};

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

const SHADOW_PREFIX: &str = "[X_SHADOW] BNDBSTREAM ";

impl Run {
    /// The counters of the `[X_SHADOW] BNDBSTREAM` line, which must end with `identical`.
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

/// A nucleotide database of `pieces(first, count)` made with `LOSAT makedb` in the directory
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
        .expect("run LOSAT makedb");
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

const OFF: (&str, &str) = ("LOSAT_X_BNDBSTREAM", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_BNDBSTREAMSHADOW", "1");
const CUT_OFF: (&str, &str) = ("LOSAT_X_BNDBCUT", "0");

/// The switch off (rounds), on (the default) and in its shadow mode, each with `extra`, give the
/// same report and standard error, the report of `one_thread` (the same search at one thread);
/// returns the shadow's counters.
fn same_reports(
    task: &str,
    db: &Path,
    threads: u32,
    outfmt: &str,
    extra: &[(&'static str, &'static str)],
    one_thread: &Run,
) -> HashMap<String, u64> {
    let what = format!("{task} threads {threads} outfmt {outfmt} {extra:?}");
    assert!(
        !one_thread.stdout.is_empty(),
        "{what}: the queries have hits"
    );
    let with = |var: Option<(&'static str, &'static str)>| -> Vec<(&'static str, &'static str)> {
        extra.iter().copied().chain(var).collect()
    };
    let off = blastn(task, db, threads, &with(Some(OFF)), outfmt);
    assert!(
        off.stdout == one_thread.stdout,
        "{what}: the report differs from one thread's without the switch"
    );
    let on = blastn(task, db, threads, &with(None), outfmt);
    assert!(
        on.stdout == one_thread.stdout,
        "{what}: the report differs with the switch"
    );
    assert_eq!(on.stderr, off.stderr, "{what}: standard error differs");
    let shadow = blastn(task, db, threads, &with(Some(SHADOW)), outfmt);
    assert!(
        shadow.stdout == one_thread.stdout,
        "{what}: the report differs in the shadow mode"
    );
    assert_eq!(
        shadow.stderr_without_shadow(),
        off.stderr,
        "{what}: standard error differs in the shadow mode"
    );
    shadow.shadow()
}

/// The counters of a shadow run where the window acted: `blocks_per_search` blocks per search, and
/// results compared.
fn assert_engaged(counts: &HashMap<String, u64>, blocks_per_search: u64, what: &str) {
    let searches = counts["searches"];
    assert!(searches >= 1, "{what}: {counts:?}");
    assert_eq!(
        counts["blocks"],
        blocks_per_search * searches,
        "{what}: {counts:?}"
    );
    assert!(counts["merged"] > 0, "{what}: {counts:?}");
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
            let counts = same_reports(task, &db, threads as u32, "6", &[], &one_thread);
            assert_engaged(
                &counts,
                blocks_of_the_cut(&db, threads),
                &format!("{task} threads {threads}"),
            );
        }
    }
    let one_thread = blastn("megablast", &db, 1, &[], "6");
    let counts = same_reports("megablast", &db, 4, "6", &[CUT_OFF], &one_thread);
    assert_engaged(&counts, 1, "blocks of 1 MiB");
    let one_thread = blastn("megablast", &db, 1, &[], "0");
    let counts = same_reports("megablast", &db, 4, "0", &[], &one_thread);
    assert_engaged(&counts, blocks_of_the_cut(&db, 4), "pairwise report");
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// At one thread the blocks are searched one after another: the shadow line counts no search.
#[test]
fn the_window_does_not_act_at_one_thread() {
    let db = database("one", 196, 12, &[]);
    let run = blastn("megablast", &db, 1, &[SHADOW], "6");
    assert_eq!(run.shadow()["searches"], 0, "{}", run.stderr);
    assert!(
        run.stdout == blastn("megablast", &db, 1, &[], "6").stdout,
        "the report differs"
    );
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// A damaged database (its `.nsq` cut short by 500 bytes) through the command line: the search
// fails with the error of the one-thread order, so standard error and the exit code are the same at
// one thread and, with the switch off, on and in its shadow mode, at 4 threads.
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
    for (what, vars) in [("off", vec![OFF]), ("on", vec![]), ("shadow", vec![SHADOW])] {
        let run = blastn_any_exit("megablast", &db, 4, &vars, "6");
        assert_eq!(run.code, one_thread.code, "{what}: exit code differs");
        assert_eq!(
            run.stderr_without_shadow(),
            one_thread.stderr_without_shadow(),
            "{what}: standard error differs from one thread's"
        );
        assert!(run.stdout == one_thread.stdout, "{what}: report differs");
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

// A database of several volumes (`-max_file_sz 60000`): at 4 threads the switch off, on and in its
// shadow mode give one thread's report, with the blocks of the cut and of 1 MiB.
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
        let counts = same_reports("megablast", &db, 4, "6", extra, &one_thread);
        assert!(counts["searches"] >= 1, "{extra:?} {counts:?}");
        assert!(
            counts["blocks"] >= volumes as u64 * counts["searches"],
            "{extra:?} {counts:?}"
        );
        assert!(counts["merged"] > 0, "{extra:?} {counts:?}");
    }
    let _ = std::fs::remove_dir_all(db.parent().unwrap());
}

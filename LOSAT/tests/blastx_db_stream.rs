//! BLASTX's `-db` search with the database read in blocks where it is searched
//! (`blastx/x_db_stream.rs`; graduated 2026-10-27, the switch was `LOSAT_X_BXDBSTREAM`) on a whole
//! run of the binary.
//!
//! The properties under test are LOSAT's own: the reports of one thread and of four threads are
//! the same for the gapped search with and without composition-based statistics and for the
//! ungapped search, with every scheduling of the passes over the database, on a database of
//! several volumes (so that each pass reads several blocks); and a database that cannot be read
//! fails the same batch with the same error. That a block holds the bytes and lengths of the
//! files is the unit test of `x_db_stream.rs`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path but those of `env`.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing.
fn losat(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(env.iter().copied());
    command
}

/// A temporary directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-db-stream-{tag}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The protein database of `fasta_file` made by `LOSAT makedb` in `dir` under `name`, in volumes
/// of at most `max_file_sz` bytes.
fn makedb(dir: &Path, fasta_file: &Path, name: &str, max_file_sz: &str) -> String {
    let db = dir.join(name).to_str().unwrap().to_string();
    let output = losat(&[])
        .arg("makedb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "prot", "-out", &db, "-max_file_sz", max_file_sz])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

/// One search of `query` on `db` with `args` at `threads` threads; returns the exit status, the
/// report and the standard error.
fn blastx(
    query: &Path,
    db: &str,
    args: &[&str],
    threads: u32,
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let output = losat(env)
        .arg("blastx")
        .arg("-query")
        .arg(query)
        .args(["-db", db])
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    (
        output.status.code(),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// The first `bases` bases of the MelaMJNV genome as one query.
fn genome_slice(dir: &Path, bases: usize) -> PathBuf {
    let genome = std::fs::read_to_string(fasta("MelaMJNV.fasta")).unwrap();
    let letters: String = genome
        .lines()
        .skip(1)
        .flat_map(str::chars)
        .take(bases)
        .collect();
    let query = dir.join("slice.fna");
    std::fs::write(&query, format!(">slice MelaMJNV slice\n{letters}\n")).unwrap();
    query
}

/// `count` pieces of `bases` bases of the MelaMJNV genome, one query each.
fn genome_pieces(dir: &Path, count: usize, bases: usize) -> PathBuf {
    let genome = std::fs::read_to_string(fasta("MelaMJNV.fasta")).unwrap();
    let letters: Vec<char> = genome.lines().skip(1).flat_map(str::chars).collect();
    let mut fasta = String::new();
    for i in 0..count {
        let piece: String = letters[i * bases..(i + 1) * bases].iter().collect();
        fasta.push_str(&format!(">piece{i} MelaMJNV piece\n{piece}\n"));
    }
    let query = dir.join("pieces.fna");
    std::fs::write(&query, fasta).unwrap();
    query
}

/// The scheduling switches whose paths read the database differently at more than one thread:
/// the default (the subject flow, the query chunks side by side, the query batches side by side),
/// the waves of the pool (`LOSAT_X_BXSUBJFLOW=0`), and the waves with the query batches as waves
/// and the index held (`LOSAT_X_BXSUBJFLOW=0 LOSAT_X_BXFLOW=0 LOSAT_X_DBIDXREAD=0`).
const SCHEDULES: [&[(&str, &str)]; 3] = [
    &[],
    &[("LOSAT_X_BXSUBJFLOW", "0")],
    &[
        ("LOSAT_X_BXSUBJFLOW", "0"),
        ("LOSAT_X_BXFLOW", "0"),
        ("LOSAT_X_DBIDXREAD", "0"),
    ],
];

#[test]
fn the_blocks_give_the_same_reports_at_every_thread_count_and_schedule() {
    let dir = temp_dir("reports");
    // A 30,000-base slice (one query batch of three query chunks) and 30 queries of 1,000 bases
    // (three query batches) against the PemoMJNVA proteins in volumes of at most 20,000 bytes.
    let slice = genome_slice(&dir, 30_000);
    let pieces = genome_pieces(&dir, 30, 1_000);
    let db = makedb(&dir, &fasta("PemoMJNVA.faa"), "pemo", "20000");
    let volumes = std::fs::read_dir(&dir)
        .unwrap()
        .filter(|entry| {
            let name = entry.as_ref().unwrap().file_name();
            let name = name.to_string_lossy();
            name.starts_with("pemo.") && name.ends_with(".psq")
        })
        .count();
    assert!(volumes > 2, "{volumes} volumes");
    let cases: [(&Path, &[&str]); 5] = [
        (&slice, &["-outfmt", "6"]),
        (&slice, &["-outfmt", "0"]),
        (&slice, &["-outfmt", "6", "-comp_based_stats", "0"]),
        (
            &slice,
            &["-outfmt", "6", "-ungapped", "-comp_based_stats", "0"],
        ),
        (&pieces, &["-outfmt", "6"]),
    ];
    for (query, args) in cases {
        let (status, expected, stderr) = blastx(query, &db, args, 1, &[]);
        assert_eq!((status, stderr.as_str()), (Some(0), ""), "{args:?}");
        assert!(
            expected.iter().filter(|&&b| b == b'\n').count() > 10,
            "{args:?}"
        );
        for schedule in SCHEDULES {
            let (status, report, stderr) = blastx(query, &db, args, 4, schedule);
            assert_eq!(
                (status, stderr.as_str()),
                (Some(0), ""),
                "{args:?} {schedule:?}"
            );
            assert!(
                report == expected,
                "{query:?} {args:?} at 4 threads with {schedule:?}: the report differs from one thread's"
            );
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// local_blast.cpp:177-224 (`CheckInternalData`): NCBI does not search a query batch without a
// valid context, so the sequence file is read by the first batch with one (as
// `blastx_db_search.rs` `the_sequence_file_is_read_by_the_first_batch_with_a_valid_context`).
#[test]
fn an_unreadable_database_fails_the_same_batch_with_the_same_error() {
    let dir = temp_dir("unreadable");
    let db = makedb(&dir, &fasta("PemoMJNVA.faa"), "pemo", "20000");
    // The last volume's sequence file cut short: the index names bytes the file does not hold.
    let last = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.to_string_lossy().ends_with(".psq"))
        .max()
        .unwrap();
    let bytes = std::fs::read(&last).unwrap();
    std::fs::write(&last, &bytes[..bytes.len() / 2]).unwrap();
    let all_n = format!(">nn all N\n{}\n", "N".repeat(10_005));
    let mut both = all_n.into_bytes();
    let slice = genome_slice(&dir, 3_000);
    both.extend_from_slice(&std::fs::read(&slice).unwrap());
    let query = dir.join("nn.fna");
    std::fs::write(&query, both).unwrap();
    let one = blastx(&query, &db, &[], 1, &[]);
    assert_ne!(one.0, Some(0), "{}", one.2);
    assert!(
        String::from_utf8_lossy(&one.1).matches("Query= ").count() == 1,
        "{}",
        String::from_utf8_lossy(&one.1)
    );
    assert!(one.2.contains("passes the end of the file"), "{}", one.2);
    assert!(blastx(&query, &db, &[], 4, &[]) == one, "4 threads");
    std::fs::remove_dir_all(&dir).unwrap();
}

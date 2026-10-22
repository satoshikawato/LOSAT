//! EXPERIMENT (LOSAT_X_BXDBSTREAM): BLASTX's `-db` search with the database read in blocks where
//! it is searched (`blastx/x_db_stream.rs`) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the reports are those of the search on the copy of
//! the whole database (`LOSAT_X_BXDBSTREAM=0`) at one and four threads, for the gapped search with
//! and without composition-based statistics and for the ungapped search, with the default blocks
//! and with blocks of 2 KiB (a few sequences each, many blocks in every pass), on a database of
//! several volumes; the shadow mode (`LOSAT_X_BXDBSTREAMSHADOW=1`) finds every read identical to
//! the copy; and a database that cannot be read fails the same batch with the same error.

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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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

#[test]
fn the_blocks_give_the_reports_of_the_copy() {
    let dir = temp_dir("reports");
    // A 30,000-base slice (three query batches) against the PemoMJNVA proteins in volumes of at
    // most 20,000 bytes.
    let query = genome_slice(&dir, 30_000);
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
    let copy = [("LOSAT_X_BXDBSTREAM", "0")];
    let small = [("LOSAT_X_BXDBSTREAM", "2")];
    for (args, outfmt) in [
        (&["-outfmt", "6"][..], "6"),
        (&["-outfmt", "0"][..], "0"),
        (&["-outfmt", "6", "-comp_based_stats", "0"][..], "6"),
        (
            &["-outfmt", "6", "-ungapped", "-comp_based_stats", "0"][..],
            "6",
        ),
    ] {
        for threads in [1, 4] {
            let (status, expected, stderr) = blastx(&query, &db, args, threads, &copy);
            assert_eq!((status, stderr.as_str()), (Some(0), ""), "{args:?}");
            assert!(!expected.is_empty(), "{args:?}");
            if outfmt == "6" {
                assert!(expected.iter().filter(|&&b| b == b'\n').count() > 10);
            }
            for env in [&[][..], &small[..]] {
                let (status, report, stderr) = blastx(&query, &db, args, threads, env);
                assert_eq!((status, stderr.as_str()), (Some(0), ""), "{args:?} {env:?}");
                assert!(
                    report == expected,
                    "{args:?} at {threads} threads with {env:?}: the report differs from the copy's"
                );
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_shadow_mode_finds_every_read_identical() {
    let dir = temp_dir("shadow");
    let query = genome_slice(&dir, 30_000);
    let db = makedb(&dir, &fasta("PemoMJNVA.faa"), "pemo", "20000");
    let shadow = [
        ("LOSAT_X_BXDBSTREAM", "2"),
        ("LOSAT_X_BXDBSTREAMSHADOW", "1"),
    ];
    // (arguments, sequences read by the passes, lengths read by the sink): the gapped search
    // splits the query into three chunks, each a pass of the preliminary search and of the sink
    // over all 115 sequences; the ungapped search searches one chunk, and its sink reads the
    // residues in a pass of its own.
    for (args, sequences, lengths) in [
        (&["-outfmt", "0"][..], 3 * 115, 3 * 115),
        (
            &["-outfmt", "6", "-ungapped", "-comp_based_stats", "0"][..],
            2 * 115,
            115,
        ),
    ] {
        let (_, expected, _) = blastx(&query, &db, args, 1, &[("LOSAT_X_BXDBSTREAM", "0")]);
        for threads in [1, 4] {
            let (status, report, stderr) = blastx(&query, &db, args, threads, &shadow);
            assert_eq!(status, Some(0), "{args:?}: {stderr}");
            assert!(report == expected, "{args:?} at {threads} threads");
            let line = stderr
                .lines()
                .find(|line| line.starts_with("[X_SHADOW] BXDBSTREAM "))
                .unwrap_or_else(|| panic!("{args:?}: no shadow line in {stderr}"));
            let count = |key: &str| -> u64 {
                line.split(' ')
                    .find_map(|field| field.strip_prefix(key))
                    .and_then(|value| value.parse().ok())
                    .unwrap_or_else(|| panic!("{key} in {line}"))
            };
            // Every pass reads each of the 115 sequences once; the composition redo and the
            // reports read the subjects they need one by one.
            assert_eq!(count("sequences="), sequences, "{line}");
            assert!(count("blocks=") > 2, "{line}");
            assert!(count("residues=") > 0, "{line}");
            assert_eq!(count("lengths="), lengths, "{line}");
            assert!(count("reads=") > 0, "{line}");
            assert!(line.ends_with(" identical"), "{line}");
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
    let copy = blastx(&query, &db, &[], 1, &[("LOSAT_X_BXDBSTREAM", "0")]);
    assert_ne!(copy.0, Some(0), "{}", copy.2);
    assert!(
        String::from_utf8_lossy(&copy.1).matches("Query= ").count() == 1,
        "{}",
        String::from_utf8_lossy(&copy.1)
    );
    for threads in [1, 4] {
        assert!(
            blastx(&query, &db, &[], threads, &[]) == copy,
            "{threads} threads"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

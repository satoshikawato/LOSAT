//! The blocks of the cut (`x_cut_block_bytes`, cut by `SeqDb::next_block`) of the fixture
//! databases of `tests/blastn_db_cut.rs`. That integration test searches these databases with
//! the command line; this test makes them in process and counts the blocks the engine cuts, the
//! fact the searches there rest on (a cut of several blocks at 2 and 4 threads for the large
//! database, one block for the small one). The generator below is the one of that file; both
//! assert the size of the sequence file, so that a change of one fails in the other too.

use std::path::{Path, PathBuf};

use super::*;
use crate::seqdb::{SeqDb, SeqType};

/// The letters of the one record of the FASTA fixture `name`.
fn genome(name: &str) -> Vec<u8> {
    let text = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fasta")
            .join(name),
    )
    .unwrap();
    text.split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// Pieces of MG1655 of 1,000 bases every 5,000 bases, `count` of them from piece `first` on;
/// piece 7 has a run of ten `N` and piece 300 an `R` and a `Y` (ambiguity data in the sequence
/// file).
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

/// A nucleotide database of `pieces(first, count)` made by `run_makeblastdb`; returns the
/// directory (to remove) and the database name.
fn database(label: &str, first: usize, count: usize) -> (PathBuf, String) {
    let dir = std::env::temp_dir().join(format!(
        "losat-bndbcut-fixture-{label}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("pieces.fna");
    std::fs::write(&input, pieces(first, count)).unwrap();
    let name = dir.join("pieces").to_str().unwrap().to_string();
    let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
    crate::api::run_makeblastdb(
        &[
            "-in",
            input.to_str().unwrap(),
            "-dbtype",
            "nucl",
            "-out",
            &name,
        ],
        &mut std::io::empty(),
        &mut output,
        &mut diagnostics,
    )
    .unwrap_or_else(|error| panic!("makeblastdb: {error}"));
    (dir, name)
}

/// The blocks that the cut of the database at `threads` threads gives.
fn blocks_of_the_cut(db: &SeqDb, threads: usize) -> u64 {
    let cut = x_cut_block_bytes(x_sequence_bytes(db), threads, BLOCK_BYTES);
    let (mut blocks, mut next) = (0, 0);
    while let Some(range) = db.next_block(next, cut) {
        blocks += 1;
        next = range.end_ordinal();
    }
    blocks
}

// The database of `reports_with_the_cut_are_the_reports_without_it`: 928 pieces, about 233 KB of
// sequence file, so one block of 1 MiB without the cut.
#[test]
fn the_large_fixture_database_is_cut_into_several_blocks() {
    let (dir, name) = database("large", 0, 10_000);
    let db = SeqDb::open(&name, SeqType::Nucleotide).unwrap();
    assert_eq!(
        std::fs::metadata(format!("{name}.nsq")).unwrap().len(),
        LARGE_NSQ_BYTES
    );
    assert_eq!(blocks_of_the_cut(&db, 1), 1);
    for threads in [2usize, 4] {
        let blocks = blocks_of_the_cut(&db, threads);
        assert!(
            blocks >= 2 * threads as u64,
            "{threads} threads: {blocks} blocks"
        );
    }
    assert_eq!(blocks_of_the_cut(&db, 2), 9);
    assert_eq!(blocks_of_the_cut(&db, 4), 17);
    drop(db);
    std::fs::remove_dir_all(dir).unwrap();
}

// The database of `the_cut_does_not_act_at_one_thread_or_on_a_database_of_one_block`: pieces
// 196-207, below the 4 KiB floor of the cut at every thread count.
#[test]
fn the_small_fixture_database_is_one_block_of_the_cut() {
    let (dir, name) = database("small", 196, 12);
    let db = SeqDb::open(&name, SeqType::Nucleotide).unwrap();
    assert_eq!(
        std::fs::metadata(format!("{name}.nsq")).unwrap().len(),
        SMALL_NSQ_BYTES
    );
    for threads in [1usize, 2, 4, 16] {
        assert_eq!(blocks_of_the_cut(&db, threads), 1, "{threads} threads");
    }
    drop(db);
    std::fs::remove_dir_all(dir).unwrap();
}

/// The size of the `.nsq` file of the large fixture database; `tests/blastn_db_cut.rs` asserts
/// the same number.
const LARGE_NSQ_BYTES: u64 = 233_200;
/// The size of the `.nsq` file of the small fixture database; `tests/blastn_db_cut.rs` asserts
/// the same number.
const SMALL_NSQ_BYTES: u64 = 3_013;

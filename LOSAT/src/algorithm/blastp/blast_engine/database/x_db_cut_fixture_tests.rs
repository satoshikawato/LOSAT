//! The blocks of the cut (`x_cut_block_bytes`, cut by `SeqDb::next_block`) of the fixture database
//! of `tests/blastp_db_cut.rs`. That integration test searches the database with the command line;
//! this test makes it in process and counts the blocks the engine cuts, the fact the searches
//! there rest on (several blocks of the cut at 2 and 4 threads). Both assert the size of the
//! sequence file, so that a change of the fixture fails in the other too.

use super::*;
use crate::blastinput::ncbi_environment::ApplicationSettings;

// The database of `search_of_all_queries_reports_with_the_cut_are_the_reports_without_it`: the
// subject proteome `MejoMJNV.faa`, 124 sequences and about 90 KB of sequence file, one block of
// 1 MiB without the cut.
#[test]
fn the_fixture_database_is_cut_into_several_blocks() {
    let dir = std::env::temp_dir().join(format!("losat-bpdbcut-fixture-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let input = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/MejoMJNV.faa");
    let name = dir.join("MejoMJNV").to_str().unwrap().to_string();
    let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
    crate::api::run_makeblastdb(
        &[
            "-in",
            input.to_str().unwrap(),
            "-dbtype",
            "prot",
            "-out",
            &name,
        ],
        &mut std::io::empty(),
        &mut output,
        &mut diagnostics,
    )
    .unwrap_or_else(|error| panic!("makeblastdb: {error}"));
    let database = BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
    assert_eq!(database.num_sequences(), 124);
    let bytes = x_sequence_bytes(&database);
    assert_eq!(
        bytes,
        std::fs::metadata(format!("{name}.psq")).unwrap().len()
    );
    assert_eq!(bytes, FIXTURE_PSQ_BYTES);
    let blocks = |threads| x_count_blocks(&database.db, x_cut_block_bytes(bytes, threads));
    assert_eq!(blocks(1), 1);
    for threads in [2usize, 4] {
        assert!(
            blocks(threads) > 2 * threads as u64,
            "{threads} threads: {} blocks",
            blocks(threads)
        );
    }
    assert_eq!(blocks(2), 9);
    assert_eq!(blocks(4), 19);
    drop(database);
    std::fs::remove_dir_all(dir).unwrap();
}

/// The size of the `.psq` file of the fixture database; `tests/blastp_db_cut.rs` asserts the
/// same number.
const FIXTURE_PSQ_BYTES: u64 = 89_564;

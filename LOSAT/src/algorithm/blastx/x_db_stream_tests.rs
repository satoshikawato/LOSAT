use super::*;
use crate::algorithm::blastx::database::BlastxDatabase;
use crate::blastinput::ncbi_environment::ApplicationSettings;

fn open(name: &str) -> BlastxDatabase {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/db")
        .join(name);
    BlastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
}

// No NCBI counterpart: LOSAT's own reads of the stored bytes (`SeqDb::read_sequence`, one
// sequence at a time, is the reference, as the copy of the whole database was before it was
// deleted). LOSAT/tests/seqdb/expected/prot_mv.info.tsv: 10 sequences, 801 residues, in six
// volumes; prot_v5 holds them in one.
#[test]
fn every_read_gives_the_bytes_and_lengths_of_the_files() {
    for name in ["prot_mv", "prot_v5"] {
        let database = open(name);
        let db = database.x_seqdb();
        let n = db.num_seqs() as usize;
        // The residues of every sequence, one read each.
        let residues: Vec<Vec<u8>> = (0..n)
            .map(|oid| db.read_sequence(oid as u32).unwrap().protein(0).to_vec())
            .collect();
        let total: usize = residues.iter().map(Vec::len).sum();
        if name == "prot_mv" {
            assert_eq!((n, total), (10, 801));
        }
        let longest = database.x_check_stores().unwrap();
        assert_eq!(
            longest,
            residues.iter().map(Vec::len).max().unwrap(),
            "{name}"
        );
        assert_eq!(database.x_open_stream().unwrap(), (longest, n), "{name}");
        let with_sentinel = |oid: usize| [&residues[oid][..], &[0u8][..]].concat();
        // The same bytes with the sentinel before them.
        let with_sentinels = |oid: usize| [&[0u8][..], &with_sentinel(oid)[..]].concat();
        // One byte: a block of one sequence; 200: a few; 1 MiB: the rest of the volume.
        for block_bytes in [1, 200, 1 << 20] {
            let stream = XStreamDb::with_block_bytes(db, longest, block_bytes);
            assert_eq!(
                (stream.len(), stream.total_length()),
                (n, db.total_length())
            );
            let mut pass = XPass::new(&stream);
            for oid in 0..n {
                let view = pass.view(oid).unwrap();
                assert_eq!(
                    &*view.encoded_with_sentinel(oid),
                    &with_sentinel(oid)[..],
                    "{name} {block_bytes} {oid}"
                );
                assert_eq!(
                    &*view.x_encoded_with_sentinels(oid),
                    &with_sentinels(oid)[..],
                    "{name} {block_bytes} {oid}"
                );
                assert_eq!(&*view.encoded(oid), &residues[oid][..]);
                assert_eq!(view.length(oid), residues[oid].len());
                assert_eq!(&*view.read_encoded(oid).unwrap(), &residues[oid][..]);
            }
            for first in 0..n {
                for end in first + 1..=n {
                    let block = stream.read_range(first, end).unwrap();
                    for oid in first..end {
                        assert_eq!(
                            block.residues_with_sentinel(oid),
                            &with_sentinel(oid)[..],
                            "{name} {first}..{end} {oid}"
                        );
                        assert_eq!(
                            block.x_residues_with_sentinels(oid),
                            &with_sentinels(oid)[..],
                            "{name} {first}..{end} {oid}"
                        );
                    }
                }
            }
            // The sink's order restarts at 0 for the next pass of a chunk.
            let mut cursor = XCursor::new(&stream);
            for oid in (0..n).chain(0..n) {
                assert_eq!(cursor.length(oid).unwrap(), residues[oid].len(), "{oid}");
            }
            let outside = XStream {
                db: &stream,
                block: None,
            };
            for oid in 0..n {
                assert_eq!(&*outside.read_encoded(oid).unwrap(), &residues[oid][..]);
            }
        }
    }
}

#[test]
fn a_read_that_fails_during_the_search_gives_the_error_of_the_check() {
    // A copy of prot_v5 whose sequence file is cut short after the check.
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-bxdbstream-read-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    for extension in ["pdb", "phr", "pin", "psq"] {
        let file = format!("prot_v5.{extension}");
        std::fs::copy(fixtures.join(&file), dir.join(&file)).unwrap();
    }
    let name = dir.join("prot_v5");
    let database =
        BlastxDatabase::open(name.to_str().unwrap(), &ApplicationSettings::default()).unwrap();
    let longest = database.x_check_stores().unwrap();
    let stream = XStreamDb::new(database.x_seqdb(), longest);
    let psq = dir.join("prot_v5.psq");
    let bytes = std::fs::read(&psq).unwrap();
    std::fs::write(&psq, &bytes[..bytes.len() / 2]).unwrap();
    let error = XPass::new(&stream).view(0).map(|_| ()).unwrap_err();
    let XReadError(kept) = error.downcast::<XReadError>().unwrap();
    assert_eq!(
        kept.to_string(),
        database.x_check_stores().unwrap_err().to_string()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

/// A copy of fixture `name` (`tests/seqdb/db`) in a directory of its own, and its path.
fn copy_of(name: &str, tag: &str) -> (std::path::PathBuf, std::path::PathBuf) {
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-bxdbstream-{tag}-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    for extension in ["pdb", "phr", "pin", "psq"] {
        let file = format!("{name}.{extension}");
        std::fs::copy(fixtures.join(&file), dir.join(&file)).unwrap();
    }
    let path = dir.join(name);
    (dir, path)
}

/// Writes `byte` over the first residue of sequence `ordinal` in the sequence file `psq` (a
/// volume of one sequence file: the residues follow the `ordinal + 1`-th NUL byte).
fn damage(psq: &std::path::Path, ordinal: usize, byte: u8) {
    let mut bytes = std::fs::read(psq).unwrap();
    let at = bytes
        .iter()
        .enumerate()
        .filter(|&(_, &b)| b == 0)
        .nth(ordinal)
        .unwrap()
        .0
        + 1;
    assert_ne!(bytes[at], 0, "sequence {ordinal} has residues");
    bytes[at] = byte;
    std::fs::write(psq, bytes).unwrap();
}

// No NCBI counterpart: LOSAT's check of the stored bytes. The head of the run reads no
// sequence: without its sequence file the database opens and gives the longest that the
// check counts, and the first read fails.
#[test]
fn the_head_of_the_stream_reads_no_sequence() {
    let (dir, path) = copy_of("prot_v5", "head");
    let longest = open("prot_v5").x_check_stores().unwrap();
    std::fs::remove_file(dir.join("prot_v5.psq")).unwrap();
    // The index header's longest, big-endian after the total length, set to 1: the head
    // counts the offsets' lengths, not the header's.
    let pin = dir.join("prot_v5.pin");
    let mut index = std::fs::read(&pin).unwrap();
    let header = (longest as u32).to_be_bytes();
    let at: Vec<usize> = (0..index.len() - 3)
        .filter(|&i| index[i..i + 4] == header)
        .collect();
    assert_eq!(at.len(), 1, "one longest length in the index");
    index[at[0]..at[0] + 4].copy_from_slice(&1u32.to_be_bytes());
    std::fs::write(&pin, index).unwrap();
    let database =
        BlastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default()).unwrap();
    assert_eq!(
        database.x_open_stream().unwrap(),
        (longest, database.x_seqdb().num_seqs() as usize)
    );
    let stream = XStreamDb::new(database.x_seqdb(), longest);
    let error = XPass::new(&stream).view(0).map(|_| ()).unwrap_err();
    let XReadError(kept) = error.downcast::<XReadError>().unwrap();
    assert_eq!(
        kept.to_string(),
        database.x_check_stores().unwrap_err().to_string()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// No NCBI counterpart: LOSAT's check of the stored bytes. A sequence with a byte that is no
// NCBIstdaa residue fails the read of every block that holds it, in a pass, a range or a
// read of its own, with the check's error, until the reads have checked every sequence; a
// read after that checks nothing (the bytes of a sequence file changed after the reads are
// not checked again, as they were not after the check of the whole database).
#[test]
fn the_reads_check_the_sequences_until_every_one_is_checked() {
    let (dir, path) = copy_of("prot_v5", "check");
    let database =
        BlastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default()).unwrap();
    let (longest, _) = database.x_open_stream().unwrap();
    let n = database.x_seqdb().num_seqs() as usize;
    damage(&dir.join("prot_v5.psq"), 6, 200);
    let want = database.x_check_stores().unwrap_err().to_string();
    assert!(want.contains("database sequence 6 holds a byte"), "{want}");
    let read_error = |error: anyhow::Error| error.downcast::<XReadError>().unwrap().0.to_string();
    for block_bytes in [1, 200, 1 << 20] {
        let stream = XStreamDb::with_block_bytes(database.x_seqdb(), longest, block_bytes);
        let mut pass = XPass::new(&stream);
        let failed = (0..n).find_map(|oid| pass.view(oid).err().map(|e| (oid, e)));
        let (oid, error) = failed.expect("the pass reads sequence 6");
        assert!(oid <= 6, "{block_bytes}: {oid}");
        assert_eq!(read_error(error), want, "{block_bytes}");
        assert_eq!(
            read_error(stream.read_range(5, 8).map(|_| ()).unwrap_err()),
            want
        );
        assert_eq!(
            read_error(stream.read_one(6).map(|_| ()).unwrap_err()),
            want
        );
        assert!(stream.read_range(7, n).is_ok() && stream.read_one(5).is_ok());
        assert!(!stream.x_checked.load(std::sync::atomic::Ordering::Acquire));
    }
    // Sequence 6 a residue again: the reads check 7 to the end in a pass, then 0..7 in a
    // range, and the database is checked; sequence 6 damaged after that is read as it is.
    damage(&dir.join("prot_v5.psq"), 6, 20);
    let stream = XStreamDb::with_block_bytes(database.x_seqdb(), longest, 200);
    let mut pass = XPass::new(&stream);
    for oid in 7..n {
        pass.view(oid).unwrap();
    }
    assert!(!stream.x_checked.load(std::sync::atomic::Ordering::Acquire));
    stream.read_range(0, 7).unwrap();
    assert!(stream.x_checked.load(std::sync::atomic::Ordering::Acquire));
    damage(&dir.join("prot_v5.psq"), 6, 200);
    assert!(stream.read_one(6).is_ok());
    std::fs::remove_dir_all(&dir).unwrap();
}

// The error of the first batch with a valid context for a database without sequences is NCBI's
// (`blastx_db_search.rs`): the reader's, not LOSAT's.
#[test]
fn the_check_of_an_empty_database_fails_as_ncbis_preliminary_search() {
    let error = open("empty").x_open_stream().unwrap_err();
    let error = error
        .downcast::<crate::algorithm::common::database::PreliminaryReadError>()
        .expect("NCBI's error of the preliminary search");
    assert_eq!(error.as_bytes(), b"Database memory map file error");
}

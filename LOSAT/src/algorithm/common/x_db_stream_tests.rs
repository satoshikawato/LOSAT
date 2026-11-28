use super::*;
use crate::seqdb::{DbList, SeqType, SeqidList, TaxidList};
use std::path::Path;
use std::sync::Barrier;

fn fixture(path: &str, seqtype: SeqType) -> SeqDb {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(path);
    SeqDb::open(path.to_str().unwrap(), seqtype).unwrap()
}

fn nucleotides(name: &str) -> SeqDb {
    fixture(&format!("seqdb/db/{name}"), SeqType::Nucleotide)
}

/// The nucleotide fixture of parsed ids without the sequence `mynuc6`, in the middle of it: a
/// block is cut there as a volume ends (`SeqDb::included_end`).
fn listed() -> SeqDb {
    let mut db = fixture("db_search/ids/nucl_v5", SeqType::Nucleotide);
    let mut warnings = Vec::new();
    let list = SeqidList {
        negative: true,
        file: crate::seqdb::seqidlist::read_si_list(b"mynuc6\n", &mut warnings).unwrap(),
    };
    db.restrict(&DbList::Seqids(list), &mut Vec::new()).unwrap();
    db
}

/// The fixture with NCBI's OIDs of a taxid list (`tests/seqdb/lmdb/expected/tx_oracle.txt`:
/// 0 1 2, 4 5 6, 9, 11 12 13), which the grid reads as it reads any other list (proteins here:
/// the grid does not read the sequences).
fn taxids() -> SeqDb {
    let mut db = fixture("seqdb/lmdb/tx", SeqType::Protein);
    db.restrict_taxids(&TaxidList::Positive([100, 200].into_iter().collect()))
        .unwrap();
    db
}

/// The database of several volumes without every third OID: runs of two OIDs, some across
/// the end of a volume.
fn volumes_and_list() -> SeqDb {
    let mut db = nucleotides("nucl_mv");
    db.x_keep_oids(|oid| oid % 3 != 2);
    db
}

/// `db`'s stream with the blocks of `max_bytes`: the fixtures are smaller than a block of the
/// search.
fn stream_of(db: &SeqDb, threads: usize, max_bytes: u64) -> XDbStream<'_> {
    let mut stream = XDbStream::new(db, threads, "TBLASTN").unwrap();
    stream.grid = x_block_grid(db, max_bytes);
    stream
}

/// The bytes and the length of sequence `oid` that a store gives.
fn stored(stream: &XDbStream<'_>, oid: usize) -> (Vec<u8>, usize) {
    let store = stream.read(oid).unwrap();
    (store.packed().to_vec(), store.len())
}

/// The blocks of more than one sequence, which `with_block` reads as blocks.
fn read_blocks(stream: &XDbStream<'_>) -> usize {
    stream.grid.windows(2).filter(|b| b[1] - b[0] > 1).count()
}

/// The grid is the chain of `next_block` from the first sequence, whatever cuts its blocks,
/// and each of its blocks is found again from its two ends, without the index: both give one
/// read of the same bytes. The check's blocks are such a chain too.
fn check_grid(db: &SeqDb, max_bytes: u64, what: &str) -> Vec<u32> {
    let grid = x_block_grid(db, max_bytes);
    assert_eq!(grid[0], 0, "{what}");
    assert_eq!(*grid.last().unwrap(), db.num_seqs(), "{what}");
    for block in grid.windows(2) {
        let range = db.next_block(block[0], max_bytes).unwrap();
        assert_eq!(range.first_ordinal, block[0], "{what}");
        assert_eq!(range.end_ordinal(), block[1], "{what}");
        assert_eq!(db.x_block_range(block[0], block[1]), Some(range), "{what}");
    }
    grid
}

#[test]
fn the_grid_is_the_chain_of_the_blocks_from_the_first_sequence() {
    for max_bytes in [0, 16, 64, 128, 256, 1024, PLAN_BLOCK_BYTES] {
        for name in ["nucl_v4", "nucl_v5", "nucl_mv"] {
            check_grid(&nucleotides(name), max_bytes, name);
        }
        check_grid(&listed(), max_bytes, "listed");
        check_grid(&taxids(), max_bytes, "taxids");
        check_grid(&volumes_and_list(), max_bytes, "volumes and list");
        // A volume without sequences, first in the order.
        check_grid(
            &fixture("seqdb/db/alias_emptyvol", SeqType::Protein),
            max_bytes,
            "empty volume",
        );
    }
    // Every volume of several starts a block, and with no limit on the bytes its last sequence
    // is in the block of its volume; with a tight one that sequence is alone.
    let db = nucleotides("nucl_mv");
    let grid = check_grid(&db, u64::MAX, "nucl_mv");
    assert!(grid.len() > 2, "several volumes: {grid:?}");
    assert!(grid.windows(2).any(|block| block[1] - block[0] > 1));
    let tight = check_grid(&db, 0, "nucl_mv");
    assert_eq!(tight, (0..=db.num_seqs()).collect::<Vec<u32>>());
    // A sequence that a list excludes cuts the block.
    let whole = check_grid(
        &fixture("db_search/ids/nucl_v5", SeqType::Nucleotide),
        u64::MAX,
        "whole",
    );
    assert_eq!(whole.len(), 2);
    assert_eq!(check_grid(&listed(), u64::MAX, "listed").len(), 3);
    // With several volumes and a list, a block ends where its run of OIDs or its volume does.
    let db = volumes_and_list();
    let mut starts = Vec::new();
    let mut volume_start = 0;
    for volume in db.volumes() {
        starts.push(volume_start);
        volume_start += volume.num_oids();
    }
    let oids: Vec<u32> = (0..db.num_seqs())
        .map(|ordinal| db.oid(ordinal).get())
        .collect();
    let want: Vec<u32> = (0..db.num_seqs())
        .filter(|&ordinal| {
            let oid = oids[ordinal as usize];
            ordinal == 0 || oids[ordinal as usize - 1] + 1 != oid || starts.contains(&oid)
        })
        .chain([db.num_seqs()])
        .collect();
    assert!(want.len() > 2);
    assert_eq!(check_grid(&db, u64::MAX, "volumes and list"), want);
    // The empty volume adds no block.
    let empty = fixture("seqdb/db/alias_emptyvol", SeqType::Protein);
    assert_eq!(
        check_grid(&empty, u64::MAX, "empty volume"),
        [0, empty.num_seqs()]
    );
    // The stream's blocks are the search's.
    let db = nucleotides("nucl_mv");
    assert_eq!(
        XDbStream::new(&db, 1, "TBLASTN").unwrap().grid,
        x_block_grid(&db, PLAN_BLOCK_BYTES)
    );
    assert_eq!(check_grid(&taxids(), u64::MAX, "taxids"), [0, 3, 6, 7, 10]);
}

/// The bytes and the length of `with_block` are those of the sequence's store, for a thread
/// that goes up, comes down, or starts anywhere; a long sequence (a block of its own) too.
#[test]
fn a_block_gives_the_bytes_of_the_stores() {
    for db in [
        nucleotides("nucl_v4"),
        nucleotides("nucl_v5"),
        nucleotides("nucl_mv"),
        listed(),
        volumes_and_list(),
    ] {
        for max_bytes in [0, 64, 256, u64::MAX] {
            let stream = stream_of(&db, 1, max_bytes);
            let count = stream.count();
            let oids = (0..count)
                .chain((0..count).rev())
                .chain((0..count).step_by(3));
            for oid in oids {
                let (packed, length) = stored(&stream, oid);
                assert_eq!(length, stream.length(oid).unwrap(), "{oid}");
                assert_eq!(length, stream.read_length(oid).unwrap(), "{oid}");
                stream
                    .with_block(oid, |sequence| {
                        assert_eq!(sequence.len(), length, "{oid}");
                        assert_eq!(sequence.packed(), packed, "{oid}");
                        Ok(())
                    })
                    .unwrap();
            }
        }
    }
}

/// A nested call replaces the thread's block; the outer call's bytes are its own, and what
/// the outer call puts back is the block of its sequence.
#[test]
fn a_nested_call_does_not_disturb_the_block_of_the_outer_one() {
    let db = nucleotides("nucl_mv");
    let stream = stream_of(&db, 1, u64::MAX);
    let blocks: Vec<&[u32]> = stream.grid.windows(2).filter(|b| b[1] - b[0] > 1).collect();
    assert!(blocks.len() > 1, "two blocks of several sequences");
    let (first, other) = (blocks[0][0] as usize, blocks[1][0] as usize);
    let (packed, _) = stored(&stream, first);
    let (inner, _) = stored(&stream, other);
    stream
        .with_block(first, |sequence| {
            stream.with_block(other, |nested| {
                assert_eq!(nested.packed(), inner);
                Ok(())
            })?;
            assert_eq!(sequence.packed(), packed);
            Ok(())
        })
        .unwrap();
    let reads = stream.block_reads.load(Ordering::Relaxed);
    stream.with_block(first, |_| Ok(())).unwrap();
    assert_eq!(stream.block_reads.load(Ordering::Relaxed), reads);
}

/// A thread that searches a chunk of subjects lets its block go with the chunk, and so does an
/// inner chunk; the block is the grid's whole block wherever the chunk ends, and while another
/// thread holds it a chunk takes it again without a read.
#[test]
fn a_block_goes_with_the_chunk_of_its_thread() {
    let db = nucleotides("nucl_v4");
    let stream = stream_of(&db, 2, u64::MAX);
    let block = stream
        .grid
        .windows(2)
        .find(|block| block[1] - block[0] > 2)
        .expect("a block of three sequences or more");
    let (first, end) = (block[0], block[1]);
    let ask = |ordinal: u32| stream.with_block(ordinal as usize, |_| Ok(())).unwrap();
    let held = || {
        let kept = X_PLAN_BLOCK.with(std::cell::Cell::take);
        let seen = kept.as_ref().map(|kept| {
            assert_eq!(
                (kept.first, kept.block.len() as u32),
                (first, end - first),
                "the grid's block"
            );
            Arc::clone(&kept.block)
        });
        X_PLAN_BLOCK.with(|cell| cell.set(kept));
        seen
    };
    let reads = || stream.block_reads.load(Ordering::Relaxed);
    // Outside a chunk the block stays.
    ask(first);
    assert!(held().is_some());
    X_PLAN_BLOCK.with(std::cell::Cell::take);
    // `None` with LOSAT_X_TBNPLANCUT=0, as in the gate's arm of the switches off.
    let Some(chunk) = x_plan_chunk() else {
        return;
    };
    let before = reads();
    ask(first);
    let other = held().expect("the chunk's block");
    {
        let _inner = x_plan_chunk();
        ask(first + 1);
        assert!(held().is_some());
    }
    assert!(held().is_none(), "the inner chunk's end lets the block go");
    ask(first + 2);
    assert!(Arc::ptr_eq(&held().unwrap(), &other));
    drop(chunk);
    assert!(held().is_none());
    assert_eq!(X_PLAN_CHUNKS.with(std::cell::Cell::get), 0);
    assert_eq!(reads(), before + 1, "taken from the thread that holds it");
    drop(other);
    ask(first);
    assert_eq!(reads(), before + 2, "read again once nothing holds it");
    X_PLAN_BLOCK.with(std::cell::Cell::take);
}

/// Threads that ask for the same sequences at once read each block once: the one that comes
/// first reads it, and the others wait for that read or take what it holds.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn threads_that_ask_together_read_each_block_once() {
    let db = nucleotides("nucl_v5");
    let stream = stream_of(&db, 4, 256);
    assert!(read_blocks(&stream) > 1);
    let together = Barrier::new(4);
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for oid in 0..stream.count() {
                    together.wait();
                    stream.with_block(oid, |_| Ok(())).unwrap();
                }
            });
        }
    });
    assert_eq!(
        stream.block_reads.load(Ordering::Relaxed),
        read_blocks(&stream)
    );
}

/// Threads that take the next sequence of one counter (TBLASTX's hand-out) get the same bytes
/// and read about one block each per pass. It is a limit, not the number: a thread held up
/// after it took its sequence reads a block again when the others have let it go.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn threads_that_take_the_next_sequence_read_each_block_about_once() {
    let db = nucleotides("nucl_v5");
    let stream = stream_of(&db, 4, 256);
    let stores: Vec<_> = (0..stream.count())
        .map(|oid| stored(&stream, oid))
        .collect();
    const PASSES: usize = 16;
    for _ in 0..PASSES {
        let next = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..4 {
                scope.spawn(|| {
                    loop {
                        let oid = next.fetch_add(1, Ordering::Relaxed);
                        if oid >= stream.count() {
                            break;
                        }
                        stream
                            .with_block(oid, |sequence| {
                                assert_eq!(
                                    (sequence.packed(), sequence.len()),
                                    (&*stores[oid].0, stores[oid].1)
                                );
                                Ok(())
                            })
                            .unwrap();
                    }
                    // `scope` can return before the destructors of its threads' locals run: the
                    // block goes now, so that no thread's block serves the next pass.
                    drop(X_PLAN_BLOCK.with(std::cell::Cell::take));
                });
            }
        });
    }
    let reads = stream.block_reads.load(Ordering::Relaxed);
    assert!(reads >= PASSES * read_blocks(&stream), "{reads}");
    assert!(reads <= 2 * PASSES * read_blocks(&stream), "{reads}");
}

#[test]
fn the_stream_holds_the_longest_length_and_their_sum() {
    for db in [
        nucleotides("nucl_v4"),
        nucleotides("nucl_v5"),
        nucleotides("nucl_mv"),
        listed(),
        volumes_and_list(),
    ] {
        let stream = XDbStream::new(&db, 1, "TBLASTN").unwrap();
        let lengths: Vec<usize> = shared::read_nucleotide_stores(&db, "TBLASTN")
            .unwrap()
            .iter()
            .map(PackedNucleotide::len)
            .collect();
        assert_eq!(stream.count(), lengths.len());
        assert_eq!(stream.first_length(), lengths[0]);
        assert_eq!(
            stream.longest_and_total(),
            (
                lengths.iter().copied().max().unwrap(),
                lengths.iter().map(|&length| length as u64).sum()
            )
        );
    }
}

/// The index's lengths (`x_approx_lengths`) are those of NCBI's scan: one per sequence in
/// ordinal order (a list's included ones: a seqid list and OIDs kept as a taxid list keeps
/// them, `SeqDb::x_keep_oids`), at most 3 bases off the stored length, exact in the bits
/// above the last byte's. When the index cannot be read, the sequences from there on
/// are given by the blocks of the grid.
#[test]
fn the_index_gives_a_length_near_the_stored_one_for_every_sequence() {
    for db in [
        nucleotides("nucl_v4"),
        nucleotides("nucl_v5"),
        nucleotides("nucl_mv"),
        listed(),
        volumes_and_list(),
    ] {
        let stream = XDbStream::new(&db, 1, "TBLASTN").unwrap();
        let mut approx = Vec::new();
        stream.x_approx_lengths(|length| approx.push(length as usize));
        assert_eq!(approx.len(), stream.count());
        for (oid, &length) in approx.iter().enumerate() {
            let stored = stream.read_length(oid).unwrap();
            assert_eq!(length / 4, stored / 4, "{oid}: {length} for {stored}");
        }
        let mut given = Vec::new();
        db.x_approx_lengths(|length| given.push(length as usize))
            .unwrap();
        assert_eq!(approx, given);
    }
    // An index read that fails after `visited` sequences: the grid's blocks give the rest,
    // a whole block's bases at the first sequence of each.
    let db = nucleotides("nucl_mv");
    let whole = (PLAN_BLOCK_BYTES * 4) as u32;
    for max_bytes in [0, 64, u64::MAX] {
        let stream = stream_of(&db, 1, max_bytes);
        let count = stream.count() as u32;
        for visited in 0..=count {
            let mut lengths = Vec::new();
            stream.x_grid_lengths(visited, |length| lengths.push(length));
            assert_eq!(lengths.len() as u32, count - visited, "{visited}");
            assert!(lengths.iter().all(|&length| length == whole || length == 0));
            let starts: Vec<u32> = stream.grid[..stream.grid.len() - 1]
                .iter()
                .copied()
                .filter(|&first| first >= visited)
                .collect();
            let given: Vec<u32> = (visited..count)
                .filter(|ordinal| lengths[(ordinal - visited) as usize] == whole)
                .collect();
            assert_eq!(given, starts, "{max_bytes} {visited}");
        }
    }
}

/// A cut is made once for its key and kept; another key makes its own.
#[test]
fn a_cut_is_made_once_for_its_key() {
    let db = nucleotides("nucl_v4");
    let stream = XDbStream::new(&db, 1, "TBLASTN").unwrap();
    let first = stream.x_kept_cut(7, || vec![1, 2]);
    let again = stream.x_kept_cut(7, || panic!("made again"));
    assert!(Arc::ptr_eq(&first, &again));
    assert_eq!(&*stream.x_kept_cut(8, || vec![3]), [3]);
    assert_eq!(&*stream.x_kept_cut(7, || panic!("made again")), [1, 2]);
}

/// A short sequence's length is read with its block, which the `with_block` that follows
/// uses: the lengths add no read of a block; a long one's is read from the index.
#[test]
fn a_length_reads_the_block_that_its_sequence_reads() {
    for db in [
        nucleotides("nucl_v4"),
        nucleotides("nucl_mv"),
        volumes_and_list(),
    ] {
        for max_bytes in [64, 256, u64::MAX] {
            let stream = stream_of(&db, 1, max_bytes);
            for oid in 0..stream.count() {
                stream.length(oid).unwrap();
                stream.with_block(oid, |_| Ok(())).unwrap();
            }
            assert_eq!(
                stream.block_reads.load(Ordering::Relaxed),
                read_blocks(&stream),
                "{max_bytes}"
            );
        }
    }
}

/// A value is read once while something holds it and again when nothing does.
#[test]
fn a_shared_read_lives_as_long_as_what_holds_it() {
    let reads = XSharedReads::<u32>::new();
    let held = reads.read_once(7, || Ok(1)).unwrap();
    let again = reads
        .read_once(7, || panic!("read while it is held"))
        .unwrap();
    assert!(Arc::ptr_eq(&held, &again));
    assert_eq!(*reads.read_once(8, || Ok(2)).unwrap(), 2);
    drop((held, again));
    assert_eq!(*reads.read_once(7, || Ok(3)).unwrap(), 3);
}

/// A thread that asks while another reads waits for that read, and takes its value.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_shared_read_is_waited_for() {
    let reads = XSharedReads::<u32>::new();
    let (started, asked) = (Barrier::new(2), Barrier::new(2));
    std::thread::scope(|scope| {
        scope.spawn(|| {
            let value = reads
                .read_once(1, || {
                    started.wait();
                    Ok(5)
                })
                .unwrap();
            asked.wait();
            drop(value);
        });
        started.wait();
        let value = reads.read_once(1, || panic!("read twice")).unwrap();
        assert_eq!(*value, 5);
        asked.wait();
    });
}

/// A read that fails or unwinds frees its key: the next thread reads, and nobody is left
/// waiting.
#[test]
#[cfg_attr(panic = "abort", ignore = "catches a panic: the build aborts on panic")]
fn a_failed_read_leaves_the_next_thread_to_read() {
    let reads = XSharedReads::<u32>::new();
    let error = reads.read_once(1, || anyhow::bail!("no")).unwrap_err();
    assert_eq!(error.to_string(), "no");
    assert_eq!(*reads.read_once(1, || Ok(1)).unwrap(), 1);
    let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = reads.read_once(2, || panic!("unwound"));
    }));
    assert!(unwound.is_err());
    assert_eq!(*reads.read_once(2, || Ok(2)).unwrap(), 2);
}

/// A thread that waits for a read that fails makes the read itself, whenever it comes.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_failed_read_wakes_the_thread_that_waits_for_it() {
    let reads = XSharedReads::<u32>::new();
    let started = Barrier::new(2);
    std::thread::scope(|scope| {
        scope.spawn(|| {
            let error = reads
                .read_once(1, || {
                    started.wait();
                    std::thread::sleep(std::time::Duration::from_millis(20));
                    anyhow::bail!("no")
                })
                .unwrap_err();
            assert_eq!(error.to_string(), "no");
        });
        started.wait();
        assert_eq!(*reads.read_once(1, || Ok(9)).unwrap(), 9);
    });
}

/// Dead entries do not pile up in the table of a long pass over a database.
#[test]
fn a_table_lets_go_of_what_nothing_holds() {
    let reads = XSharedReads::<u32>::new();
    for key in 0..10_000 {
        reads.read_once(key, || Ok(0)).unwrap();
    }
    assert!(reads.table.lock().unwrap().held.len() <= 64);
}

/// The files of the fixture `seqdb/db/{name}` copied into a directory of their own (`tag`), for
/// a test that damages them.
fn copied(name: &str, tag: &str) -> std::path::PathBuf {
    let from = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-x-db-stream-{tag}-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for entry in std::fs::read_dir(&from).unwrap() {
        let path = entry.unwrap().path();
        let file = path.file_name().unwrap().to_str().unwrap();
        if file.starts_with(&format!("{name}.")) {
            std::fs::copy(&path, dir.join(file)).unwrap();
        }
    }
    dir
}

/// The check of a stream fails where the copy of the database
/// (`common::database::read_nucleotide_stores`) fails, with its error: a sequence file
/// missing, cut at every length or with any byte changed, and an index with any byte changed,
/// in one volume and in several.
#[test]
fn the_check_fails_where_the_copy_fails_with_its_error() {
    let mut failed = 0;
    for (name, files) in [
        ("nucl_v4", &["nucl_v4.nsq", "nucl_v4.nin"][..]),
        ("nucl_mv", &["nucl_mv.01.nsq", "nucl_mv.05.nsq"][..]),
    ] {
        let dir = copied(name, &format!("check-{name}"));
        let path = dir.join(name);
        for file in files {
            let file = dir.join(file);
            let original = std::fs::read(&file).unwrap();
            let mut damaged = vec![None];
            damaged.extend((0..original.len()).map(|cut| Some(original[..cut].to_vec())));
            damaged.extend((0..original.len()).map(|at| {
                let mut bytes = original.clone();
                bytes[at] ^= 0xFF;
                Some(bytes)
            }));
            for (number, bytes) in damaged.iter().enumerate() {
                match bytes {
                    None => std::fs::remove_file(&file).unwrap(),
                    Some(bytes) => std::fs::write(&file, bytes).unwrap(),
                }
                let Ok(db) = SeqDb::open(path.to_str().unwrap(), SeqType::Nucleotide) else {
                    continue;
                };
                let what = format!("{} damage {number}", file.display());
                let check = XDbStream::new(&db, 1, "TBLASTX").map(|_| ());
                let copy = shared::read_nucleotide_stores(&db, "TBLASTX").map(|_| ());
                match (check, copy) {
                    (Ok(()), Ok(())) => {}
                    (Err(check), Err(copy)) => {
                        assert_eq!(format!("{check:#}"), format!("{copy:#}"), "{what}");
                        assert_eq!(
                            check.is::<shared::PreliminaryReadError>(),
                            copy.is::<shared::PreliminaryReadError>(),
                            "{what}"
                        );
                        failed += 1;
                    }
                    (check, copy) => panic!("{what}: check {check:?}, copy {copy:?}"),
                }
            }
            std::fs::write(&file, &original).unwrap();
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
    assert!(failed > 100, "{failed} damaged databases fail");
}

/// NCBI's scan of a database with a sequence of negative length (its ambiguity offset at its
/// start): with threads (`SeqDb::set_mt_stop_at_negative`, NCBI's multithreaded scan) the stream
/// ends before it (`count`); at one thread it is a sequence of no bases. The other sequences read
/// as without the damage, in blocks of several sequences and of one. An ambiguity offset past
/// the next sequence stays the check's error.
#[test]
fn a_negative_length_ends_the_stream_with_threads_and_is_empty_without() {
    let intact = nucleotides("nucl_v4");
    let count = intact.num_seqs() as usize;
    let dir = copied("nucl_v4", "negative");
    let name = dir.join("nucl_v4");
    let nin = dir.join("nucl_v4.nin");
    let original = std::fs::read(&nin).unwrap();
    // Nucleotide `.nin` ends with 3 arrays of (num_oids + 1) u32s: header, sequence, ambiguity.
    let stride = 4 * (count + 1);
    let seq_at = |ordinal: usize| original.len() - 2 * stride + 4 * ordinal;
    let amb_at = |ordinal: usize| original.len() - stride + 4 * ordinal;
    for ordinal in [3, count - 1] {
        let mut bytes = original.clone();
        let start = original[seq_at(ordinal)..seq_at(ordinal) + 4].to_vec();
        bytes[amb_at(ordinal)..amb_at(ordinal) + 4].copy_from_slice(&start);
        std::fs::write(&nin, &bytes).unwrap();
        for stop in [false, true] {
            let mut db = SeqDb::open(name.to_str().unwrap(), SeqType::Nucleotide).unwrap();
            db.set_mt_stop_at_negative(stop);
            for max_bytes in [0, 64] {
                let what = format!("negative at {ordinal}, stop {stop}, blocks of {max_bytes}");
                let stream = stream_of(&db, 2, max_bytes);
                let whole = stream_of(&intact, 2, max_bytes);
                assert_eq!(stream.count(), if stop { ordinal } else { count }, "{what}");
                for oid in 0..stream.count() {
                    let length = stream.with_block(oid, |sequence| Ok(sequence.len()));
                    if oid == ordinal {
                        assert_eq!(length.unwrap(), 0, "{what}");
                        assert_eq!(stream.length(oid).unwrap(), 0, "{what}");
                        assert_eq!(stream.read_length(oid).unwrap(), 0, "{what}");
                        assert_eq!(stored(&stream, oid).1, 0, "{what}");
                    } else {
                        assert_eq!(length.unwrap(), whole.length(oid).unwrap(), "{what} {oid}");
                        assert_eq!(stored(&stream, oid), stored(&whole, oid), "{what} {oid}");
                    }
                }
            }
        }
        let next = u32::from_be_bytes(original[seq_at(ordinal + 1)..][..4].try_into().unwrap());
        bytes[amb_at(ordinal)..amb_at(ordinal) + 4].copy_from_slice(&(next + 4).to_be_bytes());
        std::fs::write(&nin, &bytes).unwrap();
        for stop in [false, true] {
            let mut db = SeqDb::open(name.to_str().unwrap(), SeqType::Nucleotide).unwrap();
            db.set_mt_stop_at_negative(stop);
            let error = XDbStream::new(&db, 2, "TBLASTX").map(|_| ()).unwrap_err();
            assert!(
                format!("{error:#}").contains("could not get ambiguity data"),
                "amb > next at {ordinal}, stop {stop}: {error:#}"
            );
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// A read that fails after the check (`x_failing`) is the error of every thread that asks for
/// it, the one that reads and those that wait for its read, for a block of several sequences
/// and for a sequence alone, with NCBI's read error and with LOSAT's refusal; the other blocks
/// are read.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_read_that_fails_fails_the_reader_and_the_threads_that_wait() {
    let db = nucleotides("nucl_v4");
    for max_bytes in [0, 64] {
        let grid = x_block_grid(&db, max_bytes);
        let block = grid
            .windows(2)
            .find(|block| block[1] - block[0] > 1 || max_bytes == 0)
            .expect("a block of several sequences");
        let (first, failing) = (block[0] as usize, block[1] as usize - 1);
        for refusal in [false, true] {
            let what = format!("blocks of {max_bytes} refusal {refusal}");
            let stream = stream_of(&db, 4, max_bytes).x_failing(&[failing], refusal);
            let together = Barrier::new(4);
            let errors: Vec<String> = std::thread::scope(|scope| {
                let threads: Vec<_> = (0..4)
                    .map(|_| {
                        scope.spawn(|| {
                            together.wait();
                            let error = stream.with_block(first, |_| Ok(())).unwrap_err();
                            let XReadError(error) = error.downcast::<XReadError>().unwrap();
                            assert_eq!(
                                error.is::<shared::PreliminaryReadError>(),
                                !refusal,
                                "{what}"
                            );
                            error.to_string()
                        })
                    })
                    .collect();
                threads.into_iter().map(|t| t.join().unwrap()).collect()
            });
            assert!(errors.iter().all(|error| *error == errors[0]), "{what}");
            assert!(stream.read_shared(failing).is_err(), "{what}");
            for oid in block[1] as usize..stream.count() {
                stream.with_block(oid, |_| Ok(())).unwrap();
            }
        }
    }
}

/// `XDbStream::lead_start` returns an internal error on an unmasked or hard-masked stream, and
/// when asked for an ordinal that the one-thread walk did not record as a lead sequence.
#[test]
fn lead_start_rejects_unmasked_and_non_lead_ordinals() {
    use crate::algorithm::common::db_mask_search::SearchMask;
    use crate::blastinput::db_mask::MaskKind;

    let db = nucleotides("nucl_v4");
    for mask in [
        None,
        Some(SearchMask {
            algo_id: 11,
            kind: MaskKind::Hard,
        }),
    ] {
        let stream = XDbStream::new(&db, 1, "TBLASTN").unwrap().with_mask(mask);
        let err = stream.lead_start(0, 3, 0).unwrap_err();
        assert!(
            err.to_string().contains("LOSAT internal error"),
            "unexpected error for {mask:?}: {err}"
        );
    }

    let soft = XDbStream::new(&db, 1, "TBLASTN")
        .unwrap()
        .with_mask(Some(SearchMask {
            algo_id: 11,
            kind: MaskKind::Soft,
        }));
    let err = soft.lead_start(0, 3, 0).unwrap_err();
    assert!(
        err.to_string()
            .contains("not recorded as a soft-mask lead sequence"),
        "unexpected error: {err}"
    );
}

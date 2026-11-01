//! The database of a `-db` search, read where it is searched, instead of a copy of the whole
//! database made before the first query batch: the core of TBLASTN's stream
//! (`tblastn/x_db_stream.rs` holds what only TBLASTN's stages do with it).
//!
//! A copy packs every sequence of the database into anonymous memory before the first query
//! batch, and it stays until the reports, where NCBI maps the sequence file.
//!
//! The stream:
//!
//! - `XDbStream::new` reads the database once in blocks of `STORE_BLOCK_BYTES`, with the checks
//!   of the copy in its order (`common::database::x_check_nucleotide_stores`), and keeps the
//!   length of every sequence. A block of one sequence (a long one) is checked without its 2-bit
//!   data: its last 2-bit byte and its ambiguity data only (`SeqDb::x_nucleotide_check`). Its
//!   error is the error that the copy gave, which fails the same query batch with the same text.
//!   The check keeps no more than a `u32` per sequence, with the longest sequence and the sum of
//!   the lengths, so that a search does not walk the lengths for them.
//! - The blocks of 256 KiB of the sequence file are fixed once, from the index alone
//!   (`x_block_grid`), since `SeqDb::next_block` starts a block at the sequence asked for and
//!   threads that ask for different sequences would cut different blocks. The search of a short
//!   subject makes what it needs from the bytes of its block, without a store
//!   (`XDbStream::with_block`): each thread keeps the last block it used, and a block that another
//!   thread holds or reads is not read again. The threads that hand out the sequences one at a
//!   time read each block about once, and the text of a failed read, which holds the block's
//!   bytes, does not depend on the threads or their order.
//! - A long sequence (a block of its own) is read once while anything holds it: the stages and
//!   the query batches that run at the same time share it, and a thread that asks for it while
//!   another reads it waits for that read (`XDbStream::read_shared`). With more than one thread
//!   it is read into a buffer that glibc maps (`utils::x_mapped`), taken from the pool
//!   `utils::x_mapped::X_BYTES` when it holds one, and its store is that buffer cut to the 2-bit
//!   data, not a copy; dropped, its pages go back instead of staying in the worker's arena. At one
//!   thread it is read into the heap, where a dropped store's pages serve the next read.
//! - While a stream is open (`x_frames_pooled`), a search may give its idle large buffers to
//!   `X_BYTES`. When the last stream of the process closes, the idle buffers are freed.
//!
//! Every read packs the stored bytes as the copy did (`common::database::x_read_nucleotide_store`:
//! the packing of `PackedNucleotide::from_stored`), so every stage computes the same values on the
//! same subjects in the same order. A read that fails during the search (a database changed or
//! damaged after the check, an error of the file system) is an `XReadError` holding the error that
//! the check gives for the same read; the program decides how the batch that read it fails.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:374-376 (the sequence of
//! the search is the mapped volume's bytes, not a copy)
//! ```c++
//!     else {
//!         len = seqdb.GetSequence(oid, &buf);
//!     }
//! ```

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, Weak};

use anyhow::Result;

use crate::algorithm::common::database as shared;
use crate::seqdb::{SeqDb, SequenceBlock};
use crate::sequence::packed_subject::PackedNucleotide;

/// A read of the database that failed during the search, with
/// the error that the check (`XDbStream::new`) gives for the same read.
#[derive(Debug)]
pub(crate) struct XReadError(pub(crate) anyhow::Error);

impl std::fmt::Display for XReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for XReadError {}

/// What the threads of a stream read once while anything holds it, by a number of its own: the
/// long sequences by ordinal id (`XDbStream::read_shared`), the blocks by their number in the grid
/// (`XDbStream::with_block`).
struct XSharedReads<T> {
    table: Mutex<XShared<T>>,
    /// Woken when a read ends.
    read_done: Condvar,
}

struct XShared<T> {
    /// What something holds.
    held: HashMap<usize, Weak<T>>,
    /// What a thread reads now.
    reading: HashSet<usize>,
}

/// A read of `XSharedReads::read_once`; dropped (the read made, failed or unwound), its key is no
/// longer read and the waiters wake.
struct XReading<'s, T> {
    reads: &'s XSharedReads<T>,
    key: usize,
}

impl<T> Drop for XReading<'_, T> {
    fn drop(&mut self) {
        let mut table = self.reads.table.lock().unwrap_or_else(|e| e.into_inner());
        table.reading.remove(&self.key);
        drop(table);
        self.reads.read_done.notify_all();
    }
}

impl<T> XSharedReads<T> {
    fn new() -> Self {
        Self {
            table: Mutex::new(XShared {
                held: HashMap::new(),
                reading: HashSet::new(),
            }),
            read_done: Condvar::new(),
        }
    }

    /// What `key` names: the value that something holds, or the value of the read that a thread
    /// makes now (this thread waits for it), or else `read` of this thread, whose value the next
    /// threads share. It is read again once nothing holds it, and after a read that failed. `read`
    /// runs no rayon work, so a waiting worker depends on no task that could be stolen below it.
    fn read_once(&self, key: usize, read: impl FnOnce() -> Result<T>) -> Result<Arc<T>> {
        let mut table = self.table.lock().unwrap_or_else(|e| e.into_inner());
        loop {
            if let Some(value) = table.held.get(&key).and_then(Weak::upgrade) {
                return Ok(value);
            }
            if !table.reading.contains(&key) {
                break;
            }
            table = self
                .read_done
                .wait(table)
                .unwrap_or_else(|e| e.into_inner());
        }
        table.reading.insert(key);
        drop(table);
        let reading = XReading { reads: self, key };
        let value = Arc::new(read()?);
        let mut table = self.table.lock().unwrap_or_else(|e| e.into_inner());
        // A dead `Weak` still holds its allocation.
        if table.held.len() >= 64 && table.held.len().is_power_of_two() {
            table.held.retain(|_, value| value.strong_count() > 0);
        }
        table.held.insert(key, Arc::downgrade(&value));
        drop(table);
        drop(reading);
        Ok(value)
    }
}

/// The database of a `-db` search, read where it is searched.
pub(crate) struct XDbStream<'a> {
    db: &'a SeqDb,
    /// The program that the errors of the reads name (`common::database`).
    program: &'static str,
    /// A number of its own in the process (`XPlanBlock::stream`).
    id: u64,
    /// The length of every sequence, by ordinal id (`x_check_nucleotide_stores`).
    lengths: Vec<u32>,
    /// The largest of `lengths`.
    longest: u32,
    /// The sum of `lengths`.
    total: u64,
    /// The ordinal ids that the blocks start at, and the number of sequences (`x_block_grid`).
    grid: Vec<u32>,
    /// Whether a sequence is read into a buffer that glibc maps (`x_read_nucleotide_store`): with
    /// more than one thread. The freed stores of the workers stayed resident in their arenas; at
    /// one thread the heap's freed store serves the next read without faults (TBLASTN, C. elegans
    /// genome: 19 reads of its 7 sequences, against one read of the copy).
    mapped: bool,
    /// The long sequences that something holds or a thread reads, by ordinal id (`read_shared`).
    stores: XSharedReads<PackedNucleotide>,
    /// The blocks that something holds or a thread reads, by their number in `grid` (`with_block`).
    blocks: XSharedReads<SequenceBlock>,
    /// The blocks read (`with_block`).
    #[cfg(test)]
    block_reads: AtomicUsize,
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as the copy of the database was read (`common::database::x_check_nucleotide_stores`),
    /// keeps the lengths and fixes the blocks (`x_block_grid`). The errors name `program`.
    pub(crate) fn new(db: &'a SeqDb, threads: usize, program: &'static str) -> Result<Self> {
        static STREAMS: AtomicU64 = AtomicU64::new(0);
        let lengths = shared::x_check_nucleotide_stores(db, program)?;
        let (longest, total) = lengths.iter().fold((0, 0), |(longest, total), &length| {
            (u32::max(longest, length), total + u64::from(length))
        });
        let grid = x_block_grid(db, PLAN_BLOCK_BYTES);
        X_OPEN_STREAMS.fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            db,
            program,
            id: STREAMS.fetch_add(1, Ordering::Relaxed),
            lengths,
            longest,
            total,
            grid,
            mapped: threads > 1,
            stores: XSharedReads::new(),
            blocks: XSharedReads::new(),
            #[cfg(test)]
            block_reads: AtomicUsize::new(0),
        })
    }

    /// The number of sequences.
    pub(crate) fn count(&self) -> usize {
        self.lengths.len()
    }

    /// The number of bases of sequence `oid`.
    pub(crate) fn length(&self, oid: usize) -> usize {
        self.lengths[oid] as usize
    }

    /// The number of bases of the longest sequence and of all the sequences.
    pub(crate) fn longest_and_total(&self) -> (usize, u64) {
        (self.longest as usize, self.total)
    }

    /// Sequence `oid`, packed as the copy of the database packed it, with one read (`XReadError` if it fails).
    pub(crate) fn read(&self, oid: usize) -> Result<PackedNucleotide> {
        let ordinal = u32::try_from(oid).expect("a database ordinal id fits Uint4");
        shared::x_read_nucleotide_store(self.db, ordinal, self.program, self.mapped)
            .map_err(|error| XReadError(error).into())
    }

    /// Sequence `oid` (`read`), shared with whatever holds it already: a stage or a query batch
    /// that runs at the same time reads it once, and a thread that asks for it while another
    /// reads it waits for that read (the threads of a stage that first align the same subject at
    /// once made one copy each). It is read again once nothing holds it.
    pub(crate) fn read_shared(&self, oid: usize) -> Result<Arc<PackedNucleotide>> {
        self.stores.read_once(oid, || self.read(oid))
    }

    /// `make` of the packed bytes and the length of sequence `oid`, from the bytes of its block
    /// (`x_block_grid`), with the errors of its store but without the store. The thread keeps the
    /// last block it used; a block that another thread holds or reads is shared
    /// (`XSharedReads::read_once`). A block of one sequence (a long one) is read alone and shared
    /// (`read_shared`). The thread's block is out of its cell while `make` runs: a nested call
    /// (a pool thread that takes another subject while it waits) does not borrow it.
    pub(crate) fn with_block<R>(
        &self,
        oid: usize,
        make: impl FnOnce(&[u8], usize) -> Result<R>,
    ) -> Result<R> {
        let ordinal = u32::try_from(oid).expect("a database ordinal id fits Uint4");
        let mut kept = X_PLAN_BLOCK.with(std::cell::Cell::take);
        let held = |kept: &Option<XPlanBlock>| {
            kept.as_ref().is_some_and(|kept| {
                kept.stream == self.id
                    && (kept.first..kept.first + kept.block.len() as u32).contains(&ordinal)
            })
        };
        if !held(&kept) {
            kept = None;
            let number = self.grid.partition_point(|&first| first <= ordinal) - 1;
            let (first, end) = (self.grid[number], self.grid[number + 1]);
            if end - first == 1 {
                let store = self.read_shared(oid)?;
                return make(store.packed(), store.len());
            }
            kept = Some(XPlanBlock {
                stream: self.id,
                first,
                block: self
                    .blocks
                    .read_once(number, || self.read_grid_block(first, end))?,
            });
        }
        let block = kept.as_ref().expect("a block holds the sequence");
        let i = (ordinal - block.first) as usize;
        let made = self
            .length_in(&block.block, i, oid)
            .and_then(|length| make(block.block.nucleotide_packed(i), length));
        X_PLAN_BLOCK.with(|cell| cell.set(kept));
        made
    }

    /// The block of the sequences `first..end` of the grid (`XReadError` if it cannot be read).
    fn read_grid_block(&self, first: u32, end: u32) -> Result<SequenceBlock> {
        #[cfg(test)]
        self.block_reads.fetch_add(1, Ordering::Relaxed);
        let range = self
            .db
            .x_block_range(first, end)
            .expect("a block of the grid is a range of the database");
        self.db
            .read_block(range)
            .map_err(|error| XReadError(shared::search_read_error(error, self.program)).into())
    }

    /// The length of sequence `oid`, sequence `i` of `block`, with the errors of its store
    /// (`XReadError`).
    fn length_in(&self, block: &SequenceBlock, i: usize, oid: usize) -> Result<usize> {
        shared::x_nucleotide_length_in(block, i, oid as u32, self.program)
            .map_err(|error| XReadError(error).into())
    }
}

impl Drop for XDbStream<'_> {
    /// The calling thread (the pool's first worker) outlives the run: its block goes with the
    /// stream. The pool's other threads end with the pool.
    fn drop(&mut self) {
        let _ = X_PLAN_BLOCK.try_with(std::cell::Cell::take);
        if X_OPEN_STREAMS.fetch_sub(1, Ordering::Relaxed) == 1 {
            crate::utils::x_mapped::X_BYTES.x_clear_idle();
        }
    }
}

/// The streams open in the process (`XDbStream::new` to its drop).
static X_OPEN_STREAMS: AtomicUsize = AtomicUsize::new(0);

/// Whether a `-db` stream is open in the process: then TBLASTN's large frames of a dropped subject
/// plan go to `utils::x_mapped::X_BYTES`, from which the next frames and reads take them
/// (`tblastn/x_db_stream.rs`). A `-subject` search that runs at the same time (library use) pools
/// its frames too; a frame taken from the pool is zero-filled before it is written, so no report
/// changes. No NCBI counterpart: only where the buffers go.
pub(crate) fn x_frames_pooled() -> bool {
    X_OPEN_STREAMS.load(Ordering::Relaxed) > 0
}

/// No NCBI counterpart: the bytes of the sequence file in a block of the grid
/// (`XDbStream::with_block`; the blocks of BLASTX's `-db` stream are 256 KiB too).
const PLAN_BLOCK_BYTES: u64 = 256 << 10;

/// No NCBI counterpart: where the reads of a stream start, and no value depends on them. The
/// ordinal ids that the blocks of at most `max_bytes` of the sequence file start at, chained from
/// the first sequence (`SeqDb::next_block`, which cuts them at the end of a volume and at the
/// sequences a list excludes), then the number of sequences. It reads the index, not the sequences.
fn x_block_grid(db: &SeqDb, max_bytes: u64) -> Vec<u32> {
    let mut grid = vec![0];
    let mut first = 0;
    while let Some(range) = db.next_block(first, max_bytes) {
        first = range.end_ordinal();
        grid.push(first);
    }
    grid.shrink_to_fit();
    grid
}

/// A block of sequences that a thread keeps (`XDbStream::with_block`), and its stream.
struct XPlanBlock {
    /// `XDbStream::id` of the stream that read it.
    stream: u64,
    first: u32,
    block: Arc<SequenceBlock>,
}

thread_local! {
    /// This thread's block.
    static X_PLAN_BLOCK: std::cell::Cell<Option<XPlanBlock>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
mod tests {
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
        ] {
            for max_bytes in [0, 64, 256, u64::MAX] {
                let stream = stream_of(&db, 1, max_bytes);
                let count = stream.count();
                let oids = (0..count)
                    .chain((0..count).rev())
                    .chain((0..count).step_by(3));
                for oid in oids {
                    let (packed, length) = stored(&stream, oid);
                    assert_eq!(length, stream.length(oid));
                    stream
                        .with_block(oid, |bytes, made| {
                            assert_eq!(made, length, "{oid}");
                            assert_eq!(bytes, packed, "{oid}");
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
        let (first, last) = (0, stream.count() - 1);
        assert!(stream.grid[1] as usize <= last, "two blocks");
        let (packed, _) = stored(&stream, first);
        let (inner, _) = stored(&stream, last);
        stream
            .with_block(first, |bytes, _| {
                stream.with_block(last, |nested, _| {
                    assert_eq!(nested, inner);
                    Ok(())
                })?;
                assert_eq!(bytes, packed);
                Ok(())
            })
            .unwrap();
        let reads = stream.block_reads.load(Ordering::Relaxed);
        stream.with_block(first, |_, _| Ok(())).unwrap();
        assert_eq!(stream.block_reads.load(Ordering::Relaxed), reads);
    }

    /// Threads that ask for the same sequences at once read each block once: the one that comes
    /// first reads it, and the others wait for that read or take what it holds.
    #[test]
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
                        stream.with_block(oid, |_, _| Ok(())).unwrap();
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
                    scope.spawn(|| loop {
                        let oid = next.fetch_add(1, Ordering::Relaxed);
                        if oid >= stream.count() {
                            break;
                        }
                        stream
                            .with_block(oid, |bytes, length| {
                                assert_eq!((bytes, length), (&*stores[oid].0, stores[oid].1));
                                Ok(())
                            })
                            .unwrap();
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
        for db in [nucleotides("nucl_mv"), listed()] {
            let stream = XDbStream::new(&db, 1, "TBLASTN").unwrap();
            let lengths: Vec<usize> = (0..stream.count())
                .map(|oid| stored(&stream, oid).1)
                .collect();
            assert_eq!(
                stream.longest_and_total(),
                (
                    lengths.iter().copied().max().unwrap(),
                    lengths.iter().map(|&length| length as u64).sum()
                )
            );
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

    /// Dead entries do not pile up in the table of a long pass over a database.
    #[test]
    fn a_table_lets_go_of_what_nothing_holds() {
        let reads = XSharedReads::<u32>::new();
        for key in 0..10_000 {
            reads.read_once(key, || Ok(0)).unwrap();
        }
        assert!(reads.table.lock().unwrap().held.len() <= 64);
    }
}

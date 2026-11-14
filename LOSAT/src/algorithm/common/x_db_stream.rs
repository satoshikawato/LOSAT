//! The database of a `-db` search, read where it is searched, instead of a copy of the whole
//! database made before the first query batch: the core of TBLASTN's stream
//! (`tblastn/x_db_stream.rs` holds what only TBLASTN's stages do with it) and TBLASTX's
//! (`tblastx::subjects::Subjects::Stream`).
//!
//! A copy packs every sequence of the database into anonymous memory before the first query
//! batch, and it stays until the reports, where NCBI maps the sequence file.
//!
//! The stream:
//!
//! - `XDbStream::new` reads the database once in blocks of `STORE_BLOCK_BYTES`, with the checks
//!   of the copy in its order (`common::database::x_check_nucleotide_stores`). A block of one
//!   sequence (a long one) is checked without its 2-bit data: its last 2-bit byte and its
//!   ambiguity data only (`SeqDb::x_nucleotide_check`). Its error is the error that the copy
//!   gave, which fails the same query batch with the same text. The check keeps the number of
//!   sequences, the longest, the sum of the lengths and the first length (`XChecked`), not a
//!   length per sequence (a `u32` each was 529 MB on core_nt, held through the search; NCBI maps
//!   the index): a search reads the length of a sequence where it reads the sequence
//!   (`XDbStream::length`), a report from the index and the sequence's last byte
//!   (`XDbStream::read_length`), and TBLASTN's chunks of subjects are cut from the index's
//!   approximate lengths (`XDbStream::x_approx_lengths`).
//! - The blocks of 256 KiB of the sequence file are fixed once, from the index alone
//!   (`x_block_grid`), since `SeqDb::next_block` starts a block at the sequence asked for and
//!   threads that ask for different sequences would cut different blocks. The search of a short
//!   subject makes what it needs from the bytes of its block, without a store
//!   (`XDbStream::with_block`): each thread keeps the last block it used, and a block that another
//!   thread holds or reads is not read again. The threads that hand out the sequences one at a
//!   time read each block about once, and the error of a failed read (the text of a damaged block
//!   names the block's bytes) does not depend on the threads or their order. A thread that
//!   searches a chunk of subjects (`x_plan_chunk`) lets its block go with the chunk; the threads
//!   of the neighbouring chunks share the block that the chunk's end cuts across.
//! - A long sequence (a block of its own) is read once while anything holds it: the stages and
//!   the query batches that run at the same time share it, and a thread that asks for it while
//!   another reads it waits for that read (`XDbStream::read_shared`). With more than one thread
//!   it is read into a buffer that glibc maps (`utils::x_mapped`), taken from the pool
//!   `utils::x_mapped::X_BYTES` when it holds one, and its store is that buffer cut to the 2-bit
//!   data, not a copy; dropped, its pages go back instead of staying in the worker's arena. At one
//!   thread it is read into the heap, where a dropped store's pages serve the next read.
//! - While a stream is open (`x_frames_pooled`), a search may give its idle large buffers to
//!   `X_BYTES`. When the last stream of the process closes, the idle buffers are freed.
//! - TBLASTX makes the frames of a short subject from the bytes of its block, with the ambiguity
//!   runs that its store would keep (`XSequence::with_view`), and keeps the block for the frames
//!   of its letters when it has ambiguous bases; a long one is read alone (`read_shared`) for its
//!   frames, and again for its letters' frames when it has ambiguous bases and hits reach the
//!   re-evaluation. Each query batch reads the sequences it reports once each
//!   (`tblastx::database::read_reported`), when a report reads their letters.
//!
//! Every read packs the stored bytes as the copy did (`common::database::x_read_nucleotide_store`:
//! the packing of `PackedNucleotide::from_stored`), and a view of a block's bytes gives the letters
//! and bases of that store (`XSequence::with_view`), so every stage computes the same values on the
//! same subjects in the same order. A read that fails during the search (a database changed or
//! damaged after the check, an error of the file system) is an `XReadError` holding the error that
//! the check gives for the same read; the program decides how the batch that read it fails.
//! The check does not read the 2-bit data of a block of one sequence but its last byte, so an
//! error of the disk in the rest comes in the search's read of the sequence and not in the check:
//! TBLASTX fails the batch that reads it, if one does, as the copy failed the first batch with a
//! valid context, except that a batch without one that comes first and is reported with the batch
//! (`tblastx::report_stream`) is reported before LOSAT's refusal.
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
#[cfg(any(feature = "parallel", test))]
use std::sync::OnceLock;
use std::sync::{Arc, Condvar, Mutex, Weak};

use anyhow::Result;

use crate::algorithm::common::database as shared;
use crate::seqdb::{AmbiguityRun, SeqDb, SequenceBlock};
use crate::sequence::packed_subject::{x_stored_runs, PackedNucleotide, PackedRef};

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
    /// The number of sequences and their lengths summed up (`x_check_nucleotide_stores`).
    checked: shared::XChecked,
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
    /// The cuts made of the sequences, by their key (`x_kept_cut`).
    #[cfg(any(feature = "parallel", test))]
    cuts: Mutex<Vec<(u64, Arc<[u32]>)>>,
    /// The blocks read (`with_block`).
    #[cfg(test)]
    block_reads: AtomicUsize,
    /// The sequences whose reads fail, and whether with LOSAT's refusal (`x_failing`).
    #[cfg(test)]
    x_failing: (Vec<usize>, bool),
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as the copy of the database was read (`common::database::x_check_nucleotide_stores`)
    /// and fixes the blocks (`x_block_grid`). The errors name `program`.
    pub(crate) fn new(db: &'a SeqDb, threads: usize, program: &'static str) -> Result<Self> {
        static STREAMS: AtomicU64 = AtomicU64::new(0);
        let checked = shared::x_check_nucleotide_stores(db, program)?;
        let grid = x_block_grid(db, PLAN_BLOCK_BYTES);
        X_OPEN_STREAMS.fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            db,
            program,
            id: STREAMS.fetch_add(1, Ordering::Relaxed),
            checked,
            grid,
            mapped: threads > 1,
            stores: XSharedReads::new(),
            blocks: XSharedReads::new(),
            #[cfg(any(feature = "parallel", test))]
            cuts: Mutex::new(Vec::new()),
            #[cfg(test)]
            block_reads: AtomicUsize::new(0),
            #[cfg(test)]
            x_failing: (Vec::new(), false),
        })
    }

    /// The number of sequences.
    pub(crate) fn count(&self) -> usize {
        self.checked.count
    }

    /// The number of bases of the longest sequence and of all the sequences.
    pub(crate) fn longest_and_total(&self) -> (usize, u64) {
        (self.checked.longest as usize, self.checked.total)
    }

    /// The number of bases of the first sequence (none without one), from the check.
    pub(crate) fn first_length(&self) -> usize {
        self.checked.first as usize
    }

    /// The number of bases of sequence `ordinal`, read where the search reads the sequence: a short
    /// one's from its block, which this thread keeps for the `with_block` that follows; a long
    /// one's (a block of its own) from the index and its last byte (`read_length`), without its
    /// sequence. A read that fails is an `XReadError` with the error of the sequence's read
    /// (`with_block`, `read`). Only the length is read of the block: no ambiguity data.
    pub(crate) fn length(&self, ordinal: usize) -> Result<usize> {
        self.with_grid_block(ordinal, |held| match held {
            None => self.read_length(ordinal),
            Some((block, i)) => block
                .nucleotide_len(i)
                .map_err(|error| XReadError(shared::search_read_error(error, self.program)).into()),
        })
    }

    /// The number of bases of sequence `ordinal` read from the index and the last byte of its 2-bit
    /// data (`common::database::x_nucleotide_length`, NCBI's `GetSeqLength`), without a block or
    /// the sequence: the length of the sequence's store, for the reports that read no letters.
    /// A read that fails is an `XReadError`.
    pub(crate) fn read_length(&self, ordinal: usize) -> Result<usize> {
        #[cfg(test)]
        self.x_failed(ordinal..ordinal + 1)?;
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        shared::x_nucleotide_length(self.db, ordinal, self.program)
            .map_err(|error| XReadError(error).into())
    }

    /// Calls `visit` with the length of every sequence in ordinal order as the index gives it
    /// without a read of the sequence file (`SeqDb::x_approx_lengths`: NCBI's approximate
    /// length, at most 3 bases off the stored one, the low bits of the OID for the last byte's).
    /// When the index cannot be read (it changed after the check), the sequences from there on
    /// are given by the blocks of the grid: the first of a block with the bases of a whole block
    /// and the others with none. For a cut of the sequences that is scheduling only
    /// (`x_kept_cut`).
    #[cfg(any(feature = "parallel", test))]
    pub(crate) fn x_approx_lengths(&self, mut visit: impl FnMut(u32)) {
        let count = u32::try_from(self.count()).expect("a database ordinal id fits Uint4");
        let mut visited = 0u32;
        let _ = self.db.x_approx_lengths(|length| {
            if visited < count {
                visit(length);
                visited += 1;
            }
        });
        self.x_grid_lengths(visited, visit);
    }

    /// `x_approx_lengths` of the sequences from ordinal id `from` on that the index did not give:
    /// a block of the grid as its first sequence with the bases of a whole block and the others
    /// with none.
    #[cfg(any(feature = "parallel", test))]
    fn x_grid_lengths(&self, from: u32, mut visit: impl FnMut(u32)) {
        let count = u32::try_from(self.count()).expect("a database ordinal id fits Uint4");
        let whole = u32::try_from(PLAN_BLOCK_BYTES * 4).unwrap_or(u32::MAX);
        let mut next = self.grid.partition_point(|&first| first < from);
        for ordinal in from..count {
            let starts = self.grid.get(next) == Some(&ordinal);
            next += usize::from(starts);
            visit(if starts { whole } else { 0 });
        }
    }

    /// The cut of the sequences that `make` gives for `key` (the ends of TBLASTN's chunks of
    /// subjects for their bases, from `x_approx_lengths`), made by the first caller with `key`
    /// and kept for the run: every query batch of the search and every pool (each `key`) cuts
    /// the same sequences, which on core_nt read 1 GB of the index.
    #[cfg(any(feature = "parallel", test))]
    pub(crate) fn x_kept_cut(&self, key: u64, make: impl FnOnce() -> Vec<u32>) -> Arc<[u32]> {
        let mut cuts = self
            .cuts
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some((_, cut)) = cuts.iter().find(|(kept, _)| *kept == key) {
            return Arc::clone(cut);
        }
        let cut: Arc<[u32]> = make().into();
        cuts.push((key, Arc::clone(&cut)));
        cut
    }

    /// Sequence `ordinal`, packed as the copy of the database packed it, with one read (`XReadError` if it fails).
    pub(crate) fn read(&self, ordinal: usize) -> Result<PackedNucleotide> {
        #[cfg(test)]
        self.x_failed(ordinal..ordinal + 1)?;
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        shared::x_read_nucleotide_store(self.db, ordinal, self.program, self.mapped)
            .map_err(|error| XReadError(error).into())
    }

    /// Sequence `ordinal` (`read`), shared with whatever holds it already: a stage or a query batch
    /// that runs at the same time reads it once, and a thread that asks for it while another
    /// reads it waits for that read (the threads of a stage that first align the same subject at
    /// once made one copy each). It is read again once nothing holds it.
    pub(crate) fn read_shared(&self, ordinal: usize) -> Result<Arc<PackedNucleotide>> {
        self.stores.read_once(ordinal, || self.read(ordinal))
    }

    /// `make` of sequence `ordinal` (`XSequence`), from the bytes of its block (`x_block_grid`), with
    /// the errors of its store but without the store. The thread keeps the last block it used; a
    /// block that another thread holds or reads is shared (`XSharedReads::read_once`). A block of
    /// one sequence (a long one) is read alone and shared (`read_shared`). The thread's block is
    /// out of its cell while `make` runs: a nested call (a pool thread that takes another subject
    /// while it waits) does not borrow it. A thread that searches a chunk of subjects
    /// (`x_plan_chunk`) lets its block go with the chunk.
    pub(crate) fn with_block<R>(
        &self,
        ordinal: usize,
        make: impl FnOnce(XSequence<'_>) -> Result<R>,
    ) -> Result<R> {
        self.with_grid_block(ordinal, |held| match held {
            None => make(XSequence::Store(&self.read_shared(ordinal)?)),
            Some((block, i)) => {
                let (length, runs) = self.parts_in(block, i, ordinal)?;
                make(XSequence::Block {
                    block,
                    i,
                    length,
                    runs,
                })
            }
        })
    }

    /// `f(None)` when sequence `ordinal` is a block of its own (a long one), else `f` of its block
    /// and its place in the block, with the block kept by this thread as `with_block` keeps it:
    /// out of the thread's cell while `f` runs, and back in it after.
    fn with_grid_block<R>(
        &self,
        ordinal: usize,
        f: impl FnOnce(Option<(&Arc<SequenceBlock>, usize)>) -> Result<R>,
    ) -> Result<R> {
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        let mut kept = X_PLAN_BLOCK.with(std::cell::Cell::take);
        let holds = |kept: &Option<XPlanBlock>| {
            kept.as_ref().is_some_and(|kept| {
                kept.stream == self.id
                    && (kept.first..kept.first + kept.block.len() as u32).contains(&ordinal)
            })
        };
        if !holds(&kept) {
            // Release the thread's previous block before the next one (or a long sequence) is
            // read, which keeps the peak memory down.
            drop(kept.take());
            let number = self.grid.partition_point(|&first| first <= ordinal) - 1;
            let (first, end) = (self.grid[number], self.grid[number + 1]);
            if end - first == 1 {
                return f(None);
            }
            // LOSAT_X_TBNPLANCUT: the block goes with this thread's chunk.
            crate::utils::xreach::x_reach!(
                "LOSAT_X_TBNPLANCUT",
                X_PLAN_CHUNKS.with(std::cell::Cell::get) > 0
            );
            kept = Some(XPlanBlock {
                stream: self.id,
                first,
                block: self
                    .blocks
                    .read_once(number, || self.read_grid_block(first, end))?,
            });
        }
        let block = kept.as_ref().expect("a block holds the sequence");
        let made = f(Some((&block.block, (ordinal - block.first) as usize)));
        X_PLAN_BLOCK.with(|cell| cell.set(kept));
        made
    }

    /// The stream with the blocks of `max_bytes` (`x_block_grid`), for the tests of a program:
    /// the fixtures are smaller than a block of the search.
    #[cfg(test)]
    pub(crate) fn x_regrid(mut self, max_bytes: u64) -> Self {
        self.grid = x_block_grid(self.db, max_bytes);
        self
    }

    /// The stream whose reads of the sequences `ordinals` fail (a block's read when it holds one of
    /// them) with NCBI's read error, or with LOSAT's refusal: a read that fails after the check,
    /// which the files of a test cannot give (the check keeps their handles and lengths).
    #[cfg(test)]
    pub(crate) fn x_failing(mut self, ordinals: &[usize], refusal: bool) -> Self {
        self.x_failing = (ordinals.to_vec(), refusal);
        self
    }

    /// The error of a read of the sequences `ordinals` that `x_failing` names one of.
    #[cfg(test)]
    fn x_failed(&self, ordinals: std::ops::Range<usize>) -> Result<()> {
        let (failing, refusal) = &self.x_failing;
        let error = match failing.iter().find(|ordinal| ordinals.contains(ordinal)) {
            None => return Ok(()),
            Some(&ordinal) if *refusal => shared::refusal(
                &format!(
                    "the read of database sequence {} failed",
                    self.db.oid(ordinal as u32)
                ),
                self.program,
            ),
            Some(_) => shared::PreliminaryReadError("Database memory map file error").into(),
        };
        Err(XReadError(error).into())
    }

    /// The block of the sequences `first..end` of the grid (`XReadError` if it cannot be read).
    fn read_grid_block(&self, first: u32, end: u32) -> Result<SequenceBlock> {
        #[cfg(test)]
        self.block_reads.fetch_add(1, Ordering::Relaxed);
        #[cfg(test)]
        self.x_failed(first as usize..end as usize)?;
        let range = self
            .db
            .x_block_range(first, end)
            .expect("a block of the grid is a range of the database");
        self.db
            .read_block(range)
            .map_err(|error| XReadError(shared::search_read_error(error, self.program)).into())
    }

    /// The length and the ambiguity runs of sequence `ordinal`, sequence `i` of `block`, with the
    /// errors of its store (`XReadError`).
    pub(crate) fn parts_in(
        &self,
        block: &SequenceBlock,
        i: usize,
        ordinal: usize,
    ) -> Result<(usize, Vec<AmbiguityRun>)> {
        shared::x_nucleotide_parts_in(block, i, || self.db.oid(ordinal as u32), self.program)
            .map_err(|error| XReadError(error).into())
    }
}

/// A sequence that `XDbStream::with_block` gives.
pub(crate) enum XSequence<'b> {
    /// Sequence `i` of a block of several, with its length and its ambiguity runs as stored.
    Block {
        block: &'b Arc<SequenceBlock>,
        i: usize,
        length: usize,
        runs: Vec<AmbiguityRun>,
    },
    /// A long sequence (a block of its own), read alone and shared (`XDbStream::read_shared`).
    Store(&'b Arc<PackedNucleotide>),
}

impl XSequence<'_> {
    /// The packed bytes, as the store's (`PackedNucleotide::packed`).
    pub(crate) fn packed(&self) -> &[u8] {
        match self {
            Self::Block { block, i, .. } => block.nucleotide_packed(*i),
            Self::Store(store) => store.packed(),
        }
    }

    /// The number of bases.
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Block { length, .. } => *length,
            Self::Store(store) => store.len(),
        }
    }

    /// `f` of the sequence as its store reads it (`PackedRef`): the store's view, or the block's
    /// bytes with the runs that the store keeps (`x_stored_runs`).
    pub(crate) fn with_view<R>(self, f: impl FnOnce(PackedRef<'_>) -> R) -> R {
        match self {
            Self::Block {
                block,
                i,
                length,
                mut runs,
            } => {
                let packed = block.nucleotide_packed(i);
                x_stored_runs(packed, length, &mut runs);
                f(PackedRef::x_stored(packed, length, &runs))
            }
            Self::Store(store) => f(store.view()),
        }
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
    /// The chunks of subjects this thread searches now (`x_plan_chunk`; they nest).
    static X_PLAN_CHUNKS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

// No NCBI counterpart (memory only; every sequence is made from the same bytes): reads the
// LOSAT_X_TBNPLANCUT variable once.
/// LOSAT_X_TBNPLANCUT, in the strict set (on unless `=0`): whether a thread's block goes with the
/// chunk of subjects the thread searches (`x_plan_chunk`).
#[cfg(any(feature = "parallel", test))]
fn x_plan_cut() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNPLANCUT"))
}

/// A thread's search of a chunk of subjects (`x_plan_chunk`), until it is dropped.
#[cfg(any(feature = "parallel", test))]
pub(crate) struct XPlanChunk(());

#[cfg(any(feature = "parallel", test))]
impl Drop for XPlanChunk {
    fn drop(&mut self) {
        X_PLAN_CHUNKS.with(|chunks| chunks.set(chunks.get() - 1));
        let _ = X_PLAN_BLOCK.try_with(std::cell::Cell::take);
    }
}

/// The calling thread searches a chunk of subjects until the value returned is dropped, and the
/// block it holds then goes with the chunk. A chunk holds about a sixteenth of the sequence file
/// over the threads and a thread's next chunk is about a sixteenth of the file further on, so past
/// 4 MiB of file the next chunk is outside the block; without this a thread kept its last block
/// until its next read, through the later stages of the search. The block is the grid's
/// (`x_block_grid`), not cut at the chunk's end: the threads of the neighbouring chunks, searched
/// at the same time, take the rest of it from the shared table (`XDbStream::with_block`), where a
/// block cut at a chunk's end would be the thread's own read. Chunks nest when a thread that waits
/// runs another flow's loop (query batches side by side): the inner one's end lets the outer block
/// go, and the outer chunk takes the block of its next sequence again. `None` with
/// `LOSAT_X_TBNPLANCUT=0`.
#[cfg(any(feature = "parallel", test))]
pub(crate) fn x_plan_chunk() -> Option<XPlanChunk> {
    x_plan_cut().then(|| {
        X_PLAN_CHUNKS.with(|chunks| chunks.set(chunks.get() + 1));
        XPlanChunk(())
    })
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

    /// A thread that waits for a read that fails makes the read itself, whenever it comes.
    #[test]
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
        let dir =
            std::env::temp_dir().join(format!("losat-x-db-stream-{tag}-{}", std::process::id()));
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

    /// A read that fails after the check (`x_failing`) is the error of every thread that asks for
    /// it, the one that reads and those that wait for its read, for a block of several sequences
    /// and for a sequence alone, with NCBI's read error and with LOSAT's refusal; the other blocks
    /// are read.
    #[test]
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
}

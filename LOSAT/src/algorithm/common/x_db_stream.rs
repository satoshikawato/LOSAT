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
//! - Each thread keeps the last block of 256 KiB of the sequence file it read
//!   (`XDbStream::plan_with`), from whose bytes the search of the next short subjects makes what
//!   it needs, without a store.
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
use crate::seqdb::SeqDb;
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

/// The sequences of `XDbStream::read_shared`.
#[derive(Default)]
struct XShared {
    /// The stores that something holds, by ordinal id.
    held: HashMap<usize, Weak<PackedNucleotide>>,
    /// The sequences that a thread reads now.
    reading: HashSet<usize>,
}

/// The read of a sequence by `XDbStream::read_shared`; dropped
/// (the read made, failed or unwound), the sequence is no longer read and the waiters wake.
struct XReading<'s, 'a> {
    stream: &'s XDbStream<'a>,
    oid: usize,
}

impl Drop for XReading<'_, '_> {
    fn drop(&mut self) {
        let mut shared = self.stream.shared.lock().unwrap_or_else(|e| e.into_inner());
        shared.reading.remove(&self.oid);
        drop(shared);
        self.stream.read_done.notify_all();
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
    lengths: Vec<usize>,
    /// Whether a sequence is read into a buffer that glibc maps (`x_read_nucleotide_store`): with
    /// more than one thread. The freed stores of the workers stayed resident in their arenas; at
    /// one thread the heap's freed store serves the next read without faults (C. elegans genome:
    /// 19 reads of its 7 sequences, against one read of the copy).
    mapped: bool,
    /// The sequences that something holds or a thread reads (`read_shared`).
    shared: Mutex<XShared>,
    /// Woken when a read of `read_shared` ends.
    read_done: Condvar,
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as the copy of the database was read (`common::database::x_check_nucleotide_stores`)
    /// and keeps the lengths. The errors name `program`.
    pub(crate) fn new(db: &'a SeqDb, threads: usize, program: &'static str) -> Result<Self> {
        static STREAMS: AtomicU64 = AtomicU64::new(0);
        let lengths = shared::x_check_nucleotide_stores(db, program)?;
        X_OPEN_STREAMS.fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            db,
            program,
            id: STREAMS.fetch_add(1, Ordering::Relaxed),
            lengths,
            mapped: threads > 1,
            shared: Mutex::new(XShared::default()),
            read_done: Condvar::new(),
        })
    }

    /// The number of sequences.
    pub(crate) fn count(&self) -> usize {
        self.lengths.len()
    }

    /// The number of bases of sequence `oid`.
    pub(crate) fn length(&self, oid: usize) -> usize {
        self.lengths[oid]
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
        let mut shared = self.shared.lock().unwrap_or_else(|e| e.into_inner());
        loop {
            if let Some(store) = shared.held.get(&oid).and_then(Weak::upgrade) {
                return Ok(store);
            }
            if !shared.reading.contains(&oid) {
                break;
            }
            shared = self
                .read_done
                .wait(shared)
                .unwrap_or_else(|e| e.into_inner());
        }
        shared.reading.insert(oid);
        drop(shared);
        let reading = XReading { stream: self, oid };
        let store = Arc::new(self.read(oid)?);
        let mut shared = self.shared.lock().unwrap_or_else(|e| e.into_inner());
        if shared.held.len() >= 1024 && shared.held.len().is_power_of_two() {
            shared.held.retain(|_, store| store.strong_count() > 0);
        }
        shared.held.insert(oid, Arc::downgrade(&store));
        drop(shared);
        drop(reading);
        Ok(store)
    }

    /// `make` of the packed bytes and the length of sequence `oid`. The thread keeps the last
    /// block of sequences it read (`SeqDb::next_block` of `PLAN_BLOCK_BYTES` from the sequence
    /// asked for) and gives the next sequences of the block from its bytes, with the errors of
    /// their stores but without the stores: a search that asks for the subjects of a range in
    /// ascending order on each thread gets most short subjects without a read or a copy. A block
    /// of one sequence (a long one) is read alone and shared (`read_shared`).
    pub(crate) fn plan_with<R>(
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
            match self.db.next_block(ordinal, PLAN_BLOCK_BYTES) {
                Some(range) if range.len() > 1 => {
                    let block = self.db.read_block(range).map_err(|error| {
                        XReadError(shared::search_read_error(error, self.program))
                    })?;
                    kept = Some(XPlanBlock {
                        stream: self.id,
                        first: range.first_ordinal,
                        block,
                    });
                }
                _ => {
                    let store = self.read_shared(oid)?;
                    return make(store.packed(), store.len());
                }
            }
        }
        let block = kept.as_ref().expect("a block holds the sequence");
        let i = (ordinal - block.first) as usize;
        let made = self
            .length_in(&block.block, i, oid)
            .and_then(|length| make(block.block.nucleotide_packed(i), length));
        X_PLAN_BLOCK.with(|cell| cell.set(kept));
        made
    }

    /// The length of sequence `oid`, sequence `i` of `block`, with the errors of its store
    /// (`XReadError`).
    fn length_in(
        &self,
        block: &crate::seqdb::SequenceBlock,
        i: usize,
        oid: usize,
    ) -> Result<usize> {
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

/// No NCBI counterpart: the bytes of the sequence file in a block that a thread keeps
/// (`XDbStream::plan_with`; the blocks of BLASTX's `-db` stream are 256 KiB too).
const PLAN_BLOCK_BYTES: u64 = 256 << 10;

/// A block of sequences that a thread keeps (`XDbStream::plan_with`), and its stream.
struct XPlanBlock {
    /// `XDbStream::id` of the stream that read it.
    stream: u64,
    first: u32,
    block: crate::seqdb::SequenceBlock,
}

thread_local! {
    /// This thread's block.
    static X_PLAN_BLOCK: std::cell::Cell<Option<XPlanBlock>> = const { std::cell::Cell::new(None) };
}

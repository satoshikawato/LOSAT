//! EXPERIMENT (LOSAT_X_BPDBKEEP): the preliminary stage of a `-db` search with one search scratch
//! and one block buffer per worker for the whole stage.
//!
//! Without the switch, every task of a round of `search_blocks_of` makes a search scratch
//! (`new_scratch`: the offset pairs, 256 KiB, the diagonal table, up to 512 KiB with a query block
//! of 50,000 letters, and smaller vectors) and reads its block into a new buffer (up to 1 MiB), and
//! drops both when it ends. glibc maps an allocation of 128 KiB or more, and raises its mmap
//! threshold to the size of each mapped block that is freed (up to 32 MiB on 64-bit): after the
//! first such free the scratches and blocks of the worker threads come from their arenas, and the
//! space they leave there stays resident (the check of `IndexFile::open`, which read 1 MiB at a
//! time, raised it 0.8 s into the search). BLASTP of 100 UniRef50 queries against UniRef50 at 32
//! threads peaked at 208.7 MiB, against 136 MiB with the threshold held at 128 KiB (a diagnostic
//! `GLIBC_TUNABLES` run; v4-MF t3, v4-D12 s15).
//!
//! With the switch, each worker of the pool keeps one scratch and one block buffer for the whole
//! stage, the one-thread loop reads every block into one buffer, and the check of
//! `IndexFile::open` reads 64 KiB at a time (`seqdb/index.rs`). The workers take the blocks one at
//! a time in ordinal order (the block and the one-thread order's diagonal offset of its first
//! sequence, under one lock) and put each block's results into a window of `blocks_per_thread` x
//! workers slots; the worker that ends the oldest block not yet merged merges it and every later
//! block already ended, in ordinal order. A worker takes block j only once the blocks before j - W
//! are merged (W, the window), and otherwise waits, so the results held for the merge stay within
//! W blocks, as within a round of `search_blocks_of`, without the end of a round where every worker
//! waits for the round's last block and the merge. NCBI's threads allocate their search structures
//! once each (`BlastCoreAuxStruct`, `aux_struct->offset_pairs`, quoted at `BlastpSubjectScratch`)
//! and read the sequences in the mapped files.
//!
//! Why the reports do not change: each block is searched as a task of `search_blocks_of` searches
//! it, from `prepare_independent_subject` at the one-thread order's diagonal offset of its first
//! sequence, with a scratch that may have searched an earlier block (as a task's scratch does for
//! its next block); the blocks, their reads and checks are those of `search_blocks_of`, and the
//! results of the blocks go to the hit lists in ordinal order, as the merge of its rounds puts
//! them. A search that fails returns the first error in ordinal order, the error of the rounds of
//! `search_blocks_of` (each round's merge stops at its first error, and the rounds before had
//! none): the merge stops at the first block that failed, after every block before it, and no
//! block is taken after it (blocks already taken end unmerged). No floating-point operation
//! changes order.
//!
//! Why the waits end: a worker waits only while block j - W, the oldest not merged, is taken and
//! not ended, so another worker is searching it; that search waits for nothing a waiting worker
//! holds and runs no rayon work (the kernel `process_subject` and the block's read are sequential),
//! so no waiting task can sit below it on its thread, and its worker merges it. A worker that panics stops the search and
//! wakes the others (`XStopOnPanic`); the scope then raises the panic.
//!
//! `LOSAT_X_BPDBKEEP` is a switch of the strict set: on unless it is `0`. It acts for every task of
//! BLASTP with `-db`; the chains of BLASTP's lookup builds (`tblastx/lookup/backbone.rs`) with
//! `-subject` too; the index check and the header pieces of the title checks for every program that
//! opens a database.

use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};

use anyhow::Result;

use crate::seqdb::{BlockRange, SeqDb, SequenceBlock};
use crate::utils::threading::SearchPool;

use super::{
    block_subject, diag_offset_after, prepare_search_block, search_read_error, BlastpDatabase,
    BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, PrelimSink,
};

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The largest block of the preliminary stage with the switch (`x_db_cut.rs`
/// `x_db_cut_block_bytes`): a worker's kept block buffer stays below glibc's first mmap threshold
/// (128 KiB) and so in the worker's arena, where the space it takes is the space a task's block
/// took without the switch, instead of a mapping of its own beside the arena (1 MiB per worker:
/// BLASTP of 100 UniRef50 queries against UniRef50 at 32 threads peaked 12 MiB above the search
/// without the switch with blocks of 1 MiB, v4-D12 s17).
pub(super) const X_KEEP_BLOCK_BYTES: u64 = 96 << 10;

// No NCBI counterpart: reads the LOSAT_X_BPDBKEEP switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): a strict-set switch, on unless it is `0`.
pub(super) fn x_bp_db_keep() -> bool {
    crate::seqdb::x_bp_db_keep()
}

/// `read_search_block` with the bytes in `buf`.
fn x_read_search_block_into(db: &SeqDb, range: BlockRange, buf: Vec<u8>) -> Result<SequenceBlock> {
    let mut block = db
        .x_read_block_into(range, buf)
        .map_err(search_read_error)?;
    prepare_search_block(&mut block)?;
    Ok(block)
}

// No NCBI counterpart: where LOSAT keeps its scratch and its block bytes (NCBI's
// `BlastCoreAuxStruct` per thread, quoted at `BlastpSubjectScratch`).
/// A worker's search scratch and block buffer.
struct XKit {
    scratch: BlastpSubjectScratch,
    bytes: Vec<u8>,
}

/// The task of `search_blocks_of` for block `range` with the worker's `kit`: the block read
/// into the kit's buffer (given back to it), its sequences searched from the diagonal offset
/// `first_diag_offset`, and the results of the sequences with HSP lists.
fn x_search_block<K>(
    db: &SeqDb,
    range: BlockRange,
    first_diag_offset: i32,
    kit: &mut XKit,
    process_subject: &K,
) -> Result<Vec<BlastpSubjectResult>>
where
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let block = x_read_search_block_into(db, range, std::mem::take(&mut kit.bytes))?;
    kit.scratch.prepare_independent_subject(first_diag_offset);
    let mut with_hits = Vec::new();
    let mut subject_result = BlastpSubjectResult::new(0);
    for i in 0..block.len() {
        let s_idx = range.first_ordinal as usize + i;
        process_subject(
            &mut kit.scratch,
            s_idx,
            block_subject(&block, i),
            &mut subject_result,
        );
        if !subject_result.x_compact_hsp_lists.is_empty() {
            with_hits.push(std::mem::replace(
                &mut subject_result,
                BlastpSubjectResult::new(0),
            ));
        }
    }
    kit.bytes = block.x_into_data();
    Ok(with_hits)
}

// No NCBI counterpart: the state of the window's scheduling and merge; it changes no value NCBI
// computes.
/// The window of `x_search_window`, under its lock.
struct XWindow<'a, P> {
    /// The ordinal id after the last block taken, and the one-thread order's diagonal offset of
    /// the first sequence from there.
    next_ordinal: u32,
    diag_offset: i32,
    /// The blocks taken (the index of the next block to take).
    taken: usize,
    /// The blocks merged (the index of the next block to merge). The slot of this block is empty:
    /// the worker that ends it merges it.
    merged: usize,
    /// The results of the blocks taken and not merged, block j in slot j mod the window.
    slots: Vec<Option<Result<Vec<BlastpSubjectResult>>>>,
    sink: &'a mut P,
    /// The first error in ordinal order, once the merge has reached it.
    error: Option<anyhow::Error>,
    /// The search stopped: an error merged, or a worker panicked.
    stopped: bool,
}

/// The window and its condition: blocks merged, or the search stopped.
struct XShared<'a, P> {
    window: Mutex<XWindow<'a, P>>,
    merged: Condvar,
}

impl<'a, P> XShared<'a, P> {
    fn lock(&self) -> MutexGuard<'_, XWindow<'a, P>> {
        self.window.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

/// Stops the search and wakes the waiting workers when the worker that holds it panics, so that
/// the scope can end and raise the panic.
struct XStopOnPanic<'s, 'a, P>(&'s XShared<'a, P>);

impl<P> Drop for XStopOnPanic<'_, '_, P> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.lock().stopped = true;
            self.0.merged.notify_all();
        }
    }
}

/// A worker of `x_search_window`: takes the next block once the window has room for it, searches
/// it with the worker's kit (made for its first block), and puts its results into its slot,
/// merging them and the later blocks already ended when it is the oldest block not merged. Ends
/// when no block is left or the search stopped.
fn x_window_worker<S, K, P>(
    db: &SeqDb,
    window: i32,
    block_bytes: u64,
    shared: &XShared<'_, P>,
    new_scratch: &S,
    process_subject: &K,
) where
    S: Fn() -> BlastpSubjectScratch,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink,
{
    let _stop_on_panic = XStopOnPanic(shared);
    let mut kit: Option<XKit> = None;
    let mut state = shared.lock();
    loop {
        if state.stopped {
            return;
        }
        let width = state.slots.len();
        if state.taken >= state.merged + width {
            state = shared
                .merged
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
            continue;
        }
        let Some(range) = db.next_block(state.next_ordinal, block_bytes) else {
            return;
        };
        let index = state.taken;
        let first_diag_offset = state.diag_offset;
        state.taken += 1;
        state.diag_offset = diag_offset_after(db, range, state.diag_offset, window);
        state.next_ordinal = range.end_ordinal();
        drop(state);
        let kit = kit.get_or_insert_with(|| XKit {
            scratch: new_scratch(),
            bytes: Vec::new(),
        });
        let found = x_search_block(db, range, first_diag_offset, kit, process_subject);
        state = shared.lock();
        if state.stopped {
            return;
        }
        state.slots[index % width] = Some(found);
        if index != state.merged {
            continue;
        }
        // The subjects without HSP lists change no hit list, so the order of the updates is
        // the one-thread order.
        let merge = &mut *state;
        while let Some(found) = merge.slots[merge.merged % width].take() {
            merge.merged += 1;
            match found {
                Ok(block_results) => {
                    for mut subject_result in block_results {
                        merge.sink.add(&mut subject_result);
                    }
                }
                Err(error) => {
                    merge.error = Some(error);
                    merge.stopped = true;
                    break;
                }
            }
        }
        shared.merged.notify_all();
    }
}

// No NCBI counterpart: scheduling of the preliminary loop (NCBI's loop is quoted above
// `search_blocks`); it changes no value NCBI computes.
/// The parallel stage of `x_search_blocks_kept`: the blocks of at most `block_bytes` searched by
/// one task per worker, in a window of `blocks_per_thread` blocks per worker (module header).
#[allow(clippy::too_many_arguments)]
fn x_search_window<S, K, P>(
    db: &SeqDb,
    pool: &SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
{
    let workers = pool.threads().max(1);
    let width = workers.saturating_mul(blocks_per_thread).max(1);
    let shared = XShared {
        window: Mutex::new(XWindow {
            next_ordinal: 0,
            diag_offset: window,
            taken: 0,
            merged: 0,
            slots: std::iter::repeat_with(|| None).take(width).collect(),
            sink,
            error: None,
            stopped: false,
        }),
        merged: Condvar::new(),
    };
    {
        let shared = &shared;
        pool.install(|| {
            rayon::scope(|scope| {
                for _ in 0..workers {
                    scope.spawn(move |_| {
                        x_window_worker(
                            db,
                            window,
                            block_bytes,
                            shared,
                            new_scratch,
                            process_subject,
                        );
                    });
                }
            });
        });
    }
    let state = shared
        .window
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner);
    if let Some(error) = state.error {
        return Err(error);
    }
    debug_assert_eq!(state.merged, state.taken, "every block taken is merged");
    Ok(())
}

// No NCBI counterpart: scheduling of the preliminary loop and where its buffers are kept (NCBI's
// loop is quoted above `search_blocks`); it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): `search_blocks_of` (the same blocks of at most `block_bytes`)
/// with one scratch and one block buffer per worker for the whole stage, the blocks taken in a
/// window of `blocks_per_thread` blocks per worker (module header).
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_blocks_kept<S, K, P>(
    database: &BlastpDatabase,
    pool: &SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
{
    crate::utils::xreach::x_reach!("LOSAT_X_BPDBKEEP", true);
    let db = &database.db;
    if super::super::blastp_should_parallelize_stage(
        "prelim",
        pool.threads(),
        pool.threads(),
        db.num_seqs() as usize,
        db.num_seqs() as usize,
        0,
    ) && pool.enabled()
    {
        return x_search_window(
            db,
            pool,
            window,
            block_bytes,
            blocks_per_thread,
            new_scratch,
            process_subject,
            sink,
        );
    }

    let mut scratch = new_scratch();
    let mut subject_result = BlastpSubjectResult::new(0);
    let mut bytes = Vec::new();
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block = x_read_search_block_into(db, range, std::mem::take(&mut bytes))?;
        for i in 0..block.len() {
            process_subject(
                &mut scratch,
                range.first_ordinal as usize + i,
                block_subject(&block, i),
                &mut subject_result,
            );
            sink.add(&mut subject_result);
        }
        next_ordinal = range.end_ordinal();
        bytes = block.x_into_data();
    }
    Ok(())
}

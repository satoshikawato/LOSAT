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
//! sequence, under one lock) and put each block's results into a window of slots (`XWindowSize`);
//! the worker that ends the oldest block not yet merged merges it and every later block already
//! ended, in ordinal order. A worker takes block j only once the blocks before j - W are merged (W,
//! the window), and otherwise waits, so the results held for the merge stay within W blocks, without
//! the end of a round of `search_blocks_of` where every worker waits for the round's last block and
//! the merge. NCBI's threads allocate their search structures
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
//! Why the waits end: a worker waits only while the oldest block not merged is taken and not ended
//! (W is at least 1, and the bytes narrow it to at least 1), so another worker is searching it; that search waits for nothing a waiting worker
//! holds and runs no rayon work (the kernel `process_subject` and the block's read are sequential),
//! so no waiting task can sit below it on its thread, and its worker merges it. A worker that panics stops the search and
//! wakes the others (`XStopOnPanic`); the scope then raises the panic.
//!
//! EXPERIMENT (LOSAT_X_BPDBWIDE, a switch of the strict set: on unless it is `0`): the window is
//! `X_BPDBWIDE_BLOCKS_PER_THREAD` x workers blocks, and a worker takes a block past
//! `blocks_per_thread` x workers only while at most `X_BPDBWIDE_WAITING_BYTES_PER_THREAD` x workers
//! bytes of results wait for the merge (`XWindow::full`). A window of `blocks_per_thread` per worker
//! was shorter than the blocks at the start of a pass, where a database sorted longest first keeps
//! its slow blocks (UniRef50, blastp-fast, 100 queries, 8 threads: 64 blocks are about 38 ms of the
//! wall clock, a block of the first tenth takes 39 ms at p99 and up to 1.2 s): the workers ran out
//! of room behind them. The width and the bytes decide only when a block is taken, not what it
//! finds or the order of the merge, so they change none of the reports. The oldest block not merged
//! is always within `blocks_per_thread` x workers, so the waits end as below.
//!
//! `LOSAT_X_BPDBKEEP` is a switch of the strict set: on unless it is `0`. It acts for every task of
//! BLASTP with `-db`; the chains of BLASTP's lookup builds (`tblastx/lookup/backbone.rs`) with
//! `-subject` too; the index check and the header pieces of the title checks for every program that
//! opens a database.

use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};

use anyhow::Result;

use crate::seqdb::{BlockRange, SeqDb, SequenceBlock};
use crate::utils::threading::SearchPool;

use super::super::super::hsp::BlastpHspList;
use super::super::super::x_prelim_compact::XPrelimHsp;
use super::x_db_cut::x_sequence_bytes;
use super::{
    block_subject, diag_offset_after, prepare_search_block, BlastpDatabase, BlastpSubjectRef,
    BlastpSubjectResult, BlastpSubjectScratch, PrelimSink,
};

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The largest block of the preliminary stage with the switch (`x_db_cut.rs`
/// `x_db_cut_block_bytes`): a worker's kept block buffer stays below glibc's first mmap threshold
/// (128 KiB) and so in the worker's arena, where the space it takes is the space a task's block
/// took without the switch, instead of a mapping of its own beside the arena (1 MiB per worker:
/// BLASTP of 100 UniRef50 queries against UniRef50 at 32 threads peaked 12 MiB above the search
/// without the switch with blocks of 1 MiB, v4-D12 s17).
pub(super) const X_KEEP_BLOCK_BYTES: u64 = 96 << 10;

// No NCBI counterpart (scheduling only): the window; it changes no value NCBI computes.
/// Blocks per worker of the window with LOSAT_X_BPDBWIDE. Measured on UniRef50 with blastp-fast at
/// 16 threads, 100 queries: 8 per worker waited 130 of 840 thread-seconds, 32 waited 44, 64 waited
/// 15, 128 waited 3 and took the same time as 64 with more memory.
const X_BPDBWIDE_BLOCKS_PER_THREAD: usize = 64;

// No NCBI counterpart (scheduling only): as the window.
/// Bytes of waiting results per worker above which LOSAT_X_BPDBWIDE's window narrows to
/// `blocks_per_thread` per worker.
const X_BPDBWIDE_WAITING_BYTES_PER_THREAD: usize = 1 << 20;

// No NCBI counterpart (scheduling only): as the window.
/// The window of `x_search_window`: past `width` blocks after the oldest not merged no block is
/// taken, and past `floor` none while more than `max_waiting_bytes` of results wait.
#[derive(Clone, Copy, Debug)]
pub(super) struct XWindowSize {
    pub(super) floor: usize,
    pub(super) width: usize,
    pub(super) max_waiting_bytes: usize,
}

impl XWindowSize {
    /// The window of `workers` workers: `blocks_per_thread` per worker, or with LOSAT_X_BPDBWIDE
    /// (`wide`) the wider one that narrows to it.
    fn of(workers: usize, blocks_per_thread: usize, wide: bool) -> Self {
        let floor = workers.saturating_mul(blocks_per_thread).max(1);
        Self {
            floor,
            width: if wide {
                workers
                    .saturating_mul(X_BPDBWIDE_BLOCKS_PER_THREAD)
                    .max(floor)
            } else {
                floor
            },
            max_waiting_bytes: if wide {
                workers.saturating_mul(X_BPDBWIDE_WAITING_BYTES_PER_THREAD)
            } else {
                usize::MAX
            },
        }
    }
}

// No NCBI counterpart (scheduling only): reads the LOSAT_X_BPDBWIDE switch once.
/// Whether LOSAT_X_BPDBWIDE is on (a switch of the strict set).
fn x_bp_db_wide() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDBWIDE"))
}

// No NCBI counterpart: reads the LOSAT_X_BPDBKEEP switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): a strict-set switch, on unless it is `0`.
pub(super) fn x_bp_db_keep() -> bool {
    crate::seqdb::x_bp_db_keep()
}

/// `read_search_block` with the bytes in `buf`; the errors name `program`.
fn x_read_search_block_into(
    db: &SeqDb,
    range: BlockRange,
    buf: Vec<u8>,
    program: &str,
) -> Result<SequenceBlock> {
    let mut block = db
        .x_read_block_into(range, buf)
        .map_err(|error| crate::algorithm::common::database::search_read_error(error, program))?;
    prepare_search_block(&mut block, program)?;
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
    program: &str,
    range: BlockRange,
    first_diag_offset: i32,
    kit: &mut XKit,
    process_subject: &K,
) -> Result<Vec<BlastpSubjectResult>>
where
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let block = x_read_search_block_into(db, range, std::mem::take(&mut kit.bytes), program)?;
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

/// No NCBI counterpart (scheduling only): the heap bytes of a block's results while they wait.
fn x_kept_bytes(found: &Result<Vec<BlastpSubjectResult>>) -> usize {
    use std::mem::size_of;
    let Ok(results) = found else {
        return 0;
    };
    results
        .iter()
        .map(|result| {
            let lists = &result.x_compact_hsp_lists;
            lists
                .iter()
                .map(|(_, list)| list.hsps.capacity() * size_of::<XPrelimHsp>())
                .sum::<usize>()
                + lists.capacity() * size_of::<(usize, BlastpHspList<XPrelimHsp>)>()
        })
        .sum::<usize>()
        + results.capacity() * size_of::<BlastpSubjectResult>()
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
    /// The results of the blocks taken and not merged and their heap bytes, block j in slot j mod
    /// the window.
    slots: Vec<Option<(Result<Vec<BlastpSubjectResult>>, usize)>>,
    /// The blocks within this many of the oldest not merged are taken whatever the bytes.
    floor: usize,
    /// The heap bytes of the results in `slots`, and the most of them above which only the blocks
    /// within `floor` are taken.
    waiting_bytes: usize,
    max_waiting_bytes: usize,
    sink: &'a mut P,
    /// The first error in ordinal order, once the merge has reached it.
    error: Option<anyhow::Error>,
    /// The search stopped: an error merged, or a worker panicked.
    stopped: bool,
}

impl<P> XWindow<'_, P> {
    /// Whether the next block must wait: past the window, or past `floor` with too many bytes
    /// waiting (the module header).
    fn full(&self) -> bool {
        let past = self.taken - self.merged;
        past >= self.slots.len()
            || (past >= self.floor && self.waiting_bytes > self.max_waiting_bytes)
    }
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
    program: &str,
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
        if state.full() {
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
        let found = x_search_block(db, program, range, first_diag_offset, kit, process_subject);
        let bytes = x_kept_bytes(&found);
        state = shared.lock();
        if state.stopped {
            return;
        }
        state.waiting_bytes += bytes;
        state.slots[index % width] = Some((found, bytes));
        if index != state.merged {
            continue;
        }
        // The subjects without HSP lists change no hit list, so the order of the updates is
        // the one-thread order.
        let merge = &mut *state;
        while let Some((found, bytes)) = merge.slots[merge.merged % width].take() {
            merge.merged += 1;
            merge.waiting_bytes -= bytes;
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
/// one task per worker, in the window `size` (module header).
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_window<S, K, P>(
    database: &BlastpDatabase,
    pool: &SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    size: XWindowSize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
{
    let db = &database.db;
    let program = database.program;
    let workers = pool.threads().max(1);
    let shared = XShared {
        window: Mutex::new(XWindow {
            next_ordinal: 0,
            diag_offset: window,
            taken: 0,
            merged: 0,
            slots: std::iter::repeat_with(|| None).take(size.width).collect(),
            floor: size.floor,
            waiting_bytes: 0,
            max_waiting_bytes: size.max_waiting_bytes,
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
                            program,
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
    debug_assert_eq!(
        state.waiting_bytes, 0,
        "every block merged gave back its bytes"
    );
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
        let wide = x_bp_db_wide();
        let size = XWindowSize::of(pool.threads().max(1), blocks_per_thread, wide);
        if wide {
            // A database of `floor` blocks or fewer behaves as without the switch.
            crate::utils::xreach::x_reach!(
                "LOSAT_X_BPDBWIDE",
                x_sequence_bytes(database) / block_bytes.max(1) > size.floor as u64
            );
        }
        return x_search_window(
            database,
            pool,
            window,
            block_bytes,
            size,
            new_scratch,
            process_subject,
            sink,
        );
    }

    if x_bp_db_wide() {
        // The blocks are searched one after another, no window.
        crate::utils::xreach::x_reach!("LOSAT_X_BPDBWIDE", false);
    }
    let mut scratch = new_scratch();
    let mut subject_result = BlastpSubjectResult::new(0);
    let mut bytes = Vec::new();
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block =
            x_read_search_block_into(db, range, std::mem::take(&mut bytes), database.program)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    struct NoSink;
    impl PrelimSink for NoSink {
        fn add(&mut self, _result: &mut BlastpSubjectResult) {}
    }

    fn window(sink: &mut NoSink, width: usize, floor: usize, max: usize) -> XWindow<'_, NoSink> {
        XWindow {
            next_ordinal: 0,
            diag_offset: 0,
            taken: 0,
            merged: 0,
            slots: std::iter::repeat_with(|| None).take(width).collect(),
            floor,
            waiting_bytes: 0,
            max_waiting_bytes: max,
            sink,
            error: None,
            stopped: false,
        }
    }

    // Past the width nothing is taken; past the floor only while few bytes wait; within the floor
    // always, so the oldest block not merged always runs.
    #[test]
    fn the_window_narrows_to_the_floor_while_many_bytes_wait() {
        let mut sink = NoSink;
        let mut w = window(&mut sink, 64, 8, 1000);
        w.merged = 5;
        for (past, bytes, full) in [
            (0, 5000, false),
            (7, 5000, false),
            (8, 5000, true),
            (8, 1000, false),
            (63, 1000, false),
            (63, 1001, true),
            (64, 0, true),
        ] {
            w.taken = w.merged + past;
            w.waiting_bytes = bytes;
            assert_eq!(w.full(), full, "past {past} bytes {bytes}");
        }
    }

    #[test]
    fn the_size_of_the_window_without_and_with_the_switch() {
        let narrow = XWindowSize::of(4, 8, false);
        assert_eq!((narrow.floor, narrow.width), (32, 32));
        assert_eq!(narrow.max_waiting_bytes, usize::MAX);
        let wide = XWindowSize::of(4, 8, true);
        assert_eq!((wide.floor, wide.width), (32, 256));
        assert_eq!(wide.max_waiting_bytes, 4 << 20);
    }

    #[test]
    fn the_bytes_of_a_block_count_its_vectors() {
        use std::mem::size_of;
        let mut result = BlastpSubjectResult::new(0);
        let list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: Vec::with_capacity(10),
            best_evalue: 0.0,
        };
        let hsps = list.hsps.capacity() * size_of::<XPrelimHsp>();
        result.x_compact_hsp_lists = Vec::with_capacity(2);
        let lists =
            result.x_compact_hsp_lists.capacity() * size_of::<(usize, BlastpHspList<XPrelimHsp>)>();
        result.x_compact_hsp_lists.push((0, list));
        let block = vec![result];
        let blocks = block.capacity() * size_of::<BlastpSubjectResult>();
        assert_eq!(x_kept_bytes(&Ok(block)), hsps + lists + blocks);
        assert_eq!(x_kept_bytes(&Err(anyhow::anyhow!("failed"))), 0);
    }
}

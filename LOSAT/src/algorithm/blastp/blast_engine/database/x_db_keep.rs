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
//! stage and takes the blocks of each round one at a time from a shared counter (as the `-subject`
//! runs do, `x_subject_runs.rs`), the one-thread loop reads every block into one buffer, and the
//! check of `IndexFile::open` reads 64 KiB at a time (`seqdb/index.rs`). NCBI's threads allocate
//! their search structures once each (`BlastCoreAuxStruct`, `aux_struct->offset_pairs`, quoted at
//! `BlastpSubjectScratch`) and read the sequences in the mapped files.
//!
//! Why the reports do not change: each block is searched as a task of `search_blocks_of` searches
//! it, from `prepare_independent_subject` at the one-thread order's diagonal offset of its first
//! sequence, with a scratch that may have searched an earlier block (as a task's scratch does for
//! its next block); the blocks, their reads and checks, and the merge of each round's results in
//! ordinal order (with the round's first error in that order) are those of `search_blocks_of`. No
//! floating-point operation changes order.
//!
//! `LOSAT_X_BPDBKEEP` is a switch of the strict set: on unless it is `0`. It acts for every task of
//! BLASTP with `-db`; the index check acts for every program that opens a database.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};

use anyhow::Result;

use crate::seqdb::{BlockRange, SeqDb, SequenceBlock};
use crate::utils::threading::SearchPool;

use super::{
    block_subject, diag_offset_after, prepare_search_block, search_read_error, BlastpDatabase,
    BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, PrelimSink,
};

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
        // (EXPERIMENT LOSAT_X_BPPCOMPACT: the HSP lists may all be compact ones.)
        if !subject_result.query_hsp_lists.is_empty()
            || !subject_result.x_compact_hsp_lists.is_empty()
        {
            with_hits.push(std::mem::replace(
                &mut subject_result,
                BlastpSubjectResult::new(0),
            ));
        }
    }
    kit.bytes = block.x_into_data();
    Ok(with_hits)
}

// No NCBI counterpart: scheduling of the preliminary loop and where its buffers are kept (NCBI's
// loop is quoted above `search_blocks`); it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): `search_blocks_of` (the same blocks of at most `block_bytes`,
/// rounds of `blocks_per_thread` blocks per thread) with one scratch and one block buffer per
/// worker for the whole stage (module header).
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
    P: PrelimSink,
{
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
        let workers = pool.threads().max(1);
        let round_blocks = workers.saturating_mul(blocks_per_thread).max(1);
        let mut kits: Vec<Option<XKit>> = std::iter::repeat_with(|| None).take(workers).collect();
        let mut next_ordinal = 0u32;
        let mut diag_offset = window;
        loop {
            // The blocks of the round and the diagonal offset of each one's first sequence in
            // the one-thread order (as `search_blocks_of`).
            let mut round: Vec<(BlockRange, i32)> = Vec::with_capacity(round_blocks);
            while round.len() < round_blocks {
                let Some(range) = db.next_block(next_ordinal, block_bytes) else {
                    break;
                };
                round.push((range, diag_offset));
                diag_offset = diag_offset_after(db, range, diag_offset, window);
                next_ordinal = range.end_ordinal();
            }
            if round.is_empty() {
                break;
            }
            let found: Vec<Mutex<Option<Result<Vec<BlastpSubjectResult>>>>> =
                round.iter().map(|_| Mutex::new(None)).collect();
            let next = AtomicUsize::new(0);
            {
                let (round, found, next) = (&round, &found, &next);
                pool.install(|| {
                    rayon::scope(|scope| {
                        for kit in kits.iter_mut().take(round.len()) {
                            scope.spawn(move |_| {
                                let kit = kit.get_or_insert_with(|| XKit {
                                    scratch: new_scratch(),
                                    bytes: Vec::new(),
                                });
                                loop {
                                    let index = next.fetch_add(1, Ordering::Relaxed);
                                    let Some(&(range, first_diag_offset)) = round.get(index) else {
                                        break;
                                    };
                                    let result = x_search_block(
                                        db,
                                        range,
                                        first_diag_offset,
                                        kit,
                                        process_subject,
                                    );
                                    *found[index].lock().unwrap_or_else(PoisonError::into_inner) =
                                        Some(result);
                                }
                            });
                        }
                    });
                });
            }
            // The subjects without HSP lists change no hit list, so the order of the
            // updates is the one-thread order.
            for slot in found {
                let block_results = slot
                    .into_inner()
                    .unwrap_or_else(PoisonError::into_inner)
                    .expect("every block of the round is searched");
                for mut subject_result in block_results? {
                    sink.add(&mut subject_result);
                }
            }
        }
        return Ok(());
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

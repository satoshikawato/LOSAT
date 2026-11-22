//! EXPERIMENT (LOSAT_X_BNDBCUT): the blocks of a `-db` search's
//! preliminary stage cut by the size of the database and the thread count (BLASTN; BLASTP's
//! `LOSAT_X_BPDBCUT`, `blastp/blast_engine/database/x_db_cut.rs`, cuts its blocks by the same
//! rule).
//!
//! Without the switch, `search_blocks` reads the database in blocks of at most `BLOCK_BYTES`
//! (1 MiB) of its sequence file, and with more than one thread searches rounds of up to
//! `BLOCKS_PER_THREAD` blocks per thread, one block per task. A database of a few MiB is a handful
//! of blocks: the 10.3 MB of 28,036 16S rRNA genes are 10 blocks, so 10 of 16 threads search, and
//! their first 2,000 (0.72 MB) are one block, so one thread searches (megablast, 100 16S queries,
//! 16 threads on the workstation: `-db` took x1.34 and x3.0 the wall time of `-subject`, line G
//! `PROBE_1.md` #3).
//!
//! With the switch on, a search whose preliminary stage runs its rounds on the pool (`parallel`:
//! more than one thread) cuts the blocks at
//! clamp(sequence bytes / (`X_BLOCKS_PER_THREAD` x threads), `X_MIN_BLOCK_BYTES`, `BLOCK_BYTES`)
//! bytes, where the sequence bytes are the 2-bit bytes of the letters and the last byte of each
//! sequence (`x_sequence_bytes`, the sequence file without the ambiguity data): about 4 blocks per
//! thread, of at least 4 KiB. A database of at least 4 x threads MiB keeps blocks of 1 MiB; at one
//! thread nothing changes. The rounds, the tasks and the merge are those of `search_blocks_of`;
//! only its `block_bytes` differs.
//!
//! Why the reports do not change (`search_blocks`'s comment, which holds for any block size):
//!
//! - the cut only moves where blocks end: `SeqDb::next_block` gives the consecutive sequences of a
//!   volume that fit the size, at least one (a longer sequence is a block of its own), and the
//!   blocks follow each other in ordinal order, so every sequence is searched once, in the order of
//!   its ordinal id, as in NCBI's loop with one thread (blast_engine.c:1407-1426, cited above
//!   `search_blocks`);
//! - the result of a sequence does not depend on what the worker's scratch searched before it: the
//!   rounds of 1 MiB blocks at two threads or more already rely on it (`map_init` hands a scratch a
//!   run of consecutive blocks of the round, a run that changes from one execution to the next),
//!   as the parallel `-subject` path does. The cut gives a scratch the same kind of history, a new
//!   scratch and then consecutive sequences in ordinal order, from other first sequences;
//! - the round's results are merged in ordinal order (blast_engine.c:1553-1554; the sequences that
//!   keep nothing change no hit list), so the hit lists receive the one-thread order's updates. No
//!   floating-point operation changes order: each sequence's kernel and each update are the
//!   one-thread ones;
//! - a search that fails returns the first error of the one-thread order (`x_one_thread_error`):
//!   with another cut a damaged database can meet its first error elsewhere (the byte range of a
//!   read past the end of the file names the block, and a block reads the index entries of all its
//!   sequences before their data, so a later sequence's index error and an earlier sequence's
//!   read error come in the other order), so the blocks of `BLOCK_BYTES` are read again in order,
//!   with each sequence's length, until the first error.
//!
//! Memory: a worker holds one block, at most the cut (or one longer sequence), and a round covers at
//! most `BLOCKS_PER_THREAD` x threads blocks, so the sequence memory does not grow; a database below
//! 4 x threads MiB is one round with and without the switch, so the round holds the same results.
//! More threads search at once, each with its scratch, as in the parallel `-subject` path.
//!
//! No NCBI counterpart for the cut itself: NCBI's search threads take chunks of
//! at least 1 (or `sequences / 100`) ordinal ids from one iterator (`BlastSeqSrcIteratorNewEx`,
//! blast_engine.c:1409; seqdbimpl.cpp:328-350, cited above `SeqDb::next_block`), whatever the
//! thread count; LOSAT's blocks are cut by bytes, because one block is one read of one file.
//!
//! With LOSAT_X_BNDBSTREAM (`x_db_stream.rs`, in the strict set) the parallel stage has no rounds:
//! the threads claim the cut's blocks in ordinal order and the results are handed over in that
//! order.
//!
//! EXPERIMENT (LOSAT_X_BNDBSMALL, a switch of the strict set: on unless it is `0`): the cut's largest
//! block is `X_BNDBSMALL_BLOCK_BYTES` (256 KiB) instead of `BLOCK_BYTES`. A thread holds the block it
//! reads and the results of the blocks that wait for the owner, and those go on its arena when the
//! search starts, before the owner has handed anything over. Where most sequences hit most queries,
//! 1 MiB of 16S genes gave about 1 MiB of preliminary HSPs (`v4-D12/STATE.md`, s29). The reports do not
//! change, for the reasons above, which hold for any cut. The window keeps its bytes
//! (`x_db_stream.rs` `x_wide_blocks_per_thread`): 16 blocks per thread of 256 KiB are shorter than a
//! slow block, and core_nt took x1.16 at 128 KiB without that. The blocks are listed with the index
//! read in runs (`x_db_stream.rs` `blocks_of`, `SeqDb::x_next_blocks`): listed one at a time on one
//! thread while the others waited, four times as many blocks made ref_prok_rep_genomes x1.003-1.011
//! slower (s30).
//!
//! `LOSAT_X_BNDBCUT` is a switch of the strict set: on unless it is `0`.
//! It acts for every BLASTN task with `-db` (megablast, blastn, dc-megablast,
//! blastn-short).

use std::sync::OnceLock;

use anyhow::Result;

use crate::seqdb::SeqDb;
use crate::utils::threading::SearchPool;

use super::{
    block_subject, search_blocks_of, search_read_error, BlastnDatabase, BlastnSubjectRef,
    BLOCK_BYTES,
};

// No NCBI counterpart: reads the LOSAT_X_BNDBCUT switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BNDBCUT): a strict-set switch, on unless it is `0`.
pub(super) fn x_bn_db_cut() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNDBCUT"))
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// Blocks per thread that the cut aims at: a block costs a task, and four per thread let the
/// threads end their round together (BLASTP's measured choice, `LOSAT_X_BPDBCUT`).
pub(super) const X_BLOCKS_PER_THREAD: u64 = 4;

// No NCBI counterpart: as `X_BLOCKS_PER_THREAD`.
/// The smallest block of the cut, so that a small database is not cut into blocks of a few
/// sequences.
pub(super) const X_MIN_BLOCK_BYTES: u64 = 4 << 10;

// No NCBI counterpart: reads the LOSAT_X_BNDBSMALL switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BNDBSMALL): a strict-set switch, on unless it is `0`.
pub(super) fn x_bn_db_small() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNDBSMALL"))
}

// No NCBI counterpart: as `X_BLOCKS_PER_THREAD`.
/// The largest block of the cut with LOSAT_X_BNDBSMALL. Measured with 1,000 amplicons against
/// SILVA at 16 threads, peak anonymous memory: 94 MiB with 1 MiB, 32 with 256 KiB, 27 with 128 KiB;
/// 128 KiB took x1.03 on ref_prok_rep_genomes at 32 threads and x1.009 on core_nt with the blocks
/// listed one at a time (module documentation).
const X_BNDBSMALL_BLOCK_BYTES: u64 = 256 << 10;

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The block size of a search with `threads` threads of a database of `sequence_bytes`:
/// `BLOCK_BYTES` at one thread, else `sequence_bytes / (X_BLOCKS_PER_THREAD x threads)` between
/// `X_MIN_BLOCK_BYTES` and `largest`.
pub(super) fn x_cut_block_bytes(sequence_bytes: u64, threads: usize, largest: u64) -> u64 {
    if threads <= 1 {
        return BLOCK_BYTES;
    }
    let parts = X_BLOCKS_PER_THREAD.saturating_mul(threads as u64);
    (sequence_bytes / parts).clamp(X_MIN_BLOCK_BYTES, largest)
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The bytes of the database's sequences in its sequence files without their ambiguity data: four
/// letters per byte and the last byte of each sequence, which holds the count of its letters
/// (seqdbvol.cpp:1805-1838, cited above `block_subject`). A sequence of `n` letters takes
/// `n / 4 + 1` bytes, so this is the sum up to 3 / 4 byte per sequence. The sequences are the
/// included ones (`SeqDb::included_length`); an index that cannot be read counts as the largest
/// database, whose blocks are the cut's largest (`BLOCK_BYTES`, or `X_BNDBSMALL_BLOCK_BYTES` with
/// LOSAT_X_BNDBSMALL on). The search then fails, and `x_one_thread_error` reads the blocks of
/// `BLOCK_BYTES` again to return the error of the one-thread order.
pub(super) fn x_sequence_bytes(db: &SeqDb) -> u64 {
    let length = db.included_length().unwrap_or(u64::MAX);
    (length / 4).saturating_add(u64::from(db.num_seqs()))
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The block size of the search of `db` on `pool`: the cut when the preliminary stage runs its
/// rounds on the pool (`parallel`), else `BLOCK_BYTES`.
pub(super) fn x_db_cut_block_bytes(db: &SeqDb, pool: &SearchPool<'_>, parallel: bool) -> u64 {
    let threads = if parallel { pool.threads() } else { 1 };
    if threads <= 1 {
        return BLOCK_BYTES;
    }
    let sequence_bytes = x_sequence_bytes(db);
    let cut = x_cut_block_bytes(sequence_bytes, threads, BLOCK_BYTES);
    let small = x_bn_db_small().then_some(X_BNDBSMALL_BLOCK_BYTES);
    let block_bytes = x_cut_block_bytes(sequence_bytes, threads, small.unwrap_or(BLOCK_BYTES));
    crate::utils::xreach::x_reach!("LOSAT_X_BNDBSMALL", block_bytes < cut);
    block_bytes
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
// The one-thread loop of `search_blocks_of` reads its blocks of `BLOCK_BYTES` in ordinal order and
// stops at the first block that cannot be read or the first sequence whose length is rejected
// (`block_subject`); the kernel itself returns no error. The parallel rounds of the same blocks
// return the same error (the round's first in ordinal order).
/// The error that the search of `db` with blocks of `BLOCK_BYTES` returns, found by reading those
/// blocks in order and each sequence's length; `error` (the cut's) when every block and length
/// reads (an error that does not come back on a second read).
pub(super) fn x_one_thread_error(db: &SeqDb, error: anyhow::Error) -> anyhow::Error {
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, BLOCK_BYTES) {
        let block = match db.read_block(range) {
            Ok(block) => block,
            Err(first) => return search_read_error(first),
        };
        for i in 0..block.len() {
            if let Err(first) = block_subject(&block, &range, i) {
                return first;
            }
        }
        next_ordinal = range.end_ordinal();
    }
    error
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1413
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
/// EXPERIMENT (LOSAT_X_BNDBCUT): `search_blocks` with the blocks of the cut
/// (`x_db_cut_block_bytes`) and, on failure, the error of the one-thread order.
///
/// The output is that of the blocks of `BLOCK_BYTES` (the one-thread order's) for any cut: the cut
/// only moves where blocks end, a sequence's result does not depend on what the worker's scratch
/// searched before it, and the round's results are merged in ordinal order (module header). No
/// floating-point operation changes order.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_blocks_cut<T, R, I, K, M>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    blocks_per_thread: usize,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send + super::XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    let db = &database.db;
    let block_bytes = x_db_cut_block_bytes(db, pool, parallel);
    crate::utils::xreach::x_reach!("LOSAT_X_BNDBCUT", block_bytes != BLOCK_BYTES);
    let one_thread_error = |error| {
        if block_bytes == BLOCK_BYTES {
            error
        } else {
            x_one_thread_error(db, error)
        }
    };
    search_blocks_of(
        database,
        pool,
        parallel,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        search,
        merge,
    )
    .map_err(one_thread_error)
}

#[cfg(test)]
#[path = "x_db_cut_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "x_db_cut_fixture_tests.rs"]
mod fixture_tests;

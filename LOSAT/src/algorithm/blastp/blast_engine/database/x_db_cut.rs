//! EXPERIMENT (LOSAT_X_BPDBCUT): the blocks of a `-db` search's preliminary stage cut by the size
//! of the database and the thread count.
//!
//! Without the switch, `search_blocks` reads the database in blocks of at most `BLOCK_BYTES`
//! (1 MiB) of its sequence file, and with more than one thread searches rounds of up to
//! `BLOCKS_PER_THREAD` blocks per thread, one block per task. A database of a few MiB is a
//! handful of blocks: the 2.92 MB of `BlastDBSpecies6` (an OrthoFinder species database) are 3
//! blocks, so 3 of 16 threads search and the stage took 12.85 s of wall for 33.7 s of CPU on the
//! workstation (perf-mb-20261014 `logs/split2_report.md`).
//!
//! With the switch on, a search with more than one thread on the pool cuts the blocks at
//! clamp(sequence bytes / (`X_BLOCKS_PER_THREAD` x threads), `X_MIN_BLOCK_BYTES`, `BLOCK_BYTES`)
//! bytes, where the sequence bytes are the letters, one sentinel byte per sequence and one before
//! the first (the size of a one-volume sequence file; `x_sequence_bytes`): about 4 blocks per
//! thread, of at least 4 KiB. A database of at least 4 x threads MiB keeps blocks of 1 MiB; at one
//! thread nothing changes. The rounds, the tasks and the merge are those of `search_blocks_of`; only
//! its `block_bytes` differs.
//!
//! Why the reports do not change (the integer argument of `search_blocks`'s comment, which already
//! holds for any block size):
//!
//! - the cut only moves where blocks end: `SeqDb::next_block` gives the consecutive sequences of a
//!   volume that fit the size, at least one (a longer sequence is a block of its own), and the
//!   blocks follow each other in ordinal order, so every sequence is searched once, in the order of
//!   its ordinal id, as in NCBI's loop with one thread (blast_engine.c:1407-1413, quoted above
//!   `search_blocks`);
//! - every block starts from a cleared diagonal table at the offset that the one-thread order gives
//!   its first sequence (`prepare_independent_subject`, `diag_offset_after`, the recurrence of
//!   `Blast_ExtendWordExit`, blast_extend.c:166-173 quoted above `diag_offset_after`) and advances it
//!   inside the block as the one-thread loop does, so each sequence is searched with the decisions
//!   of the one-thread loop (an entry left by an earlier sequence reads as a cleared entry;
//!   `x_subject_runs.rs` gives the argument in full, and its one exception, a sequence of more than
//!   2^29 letters in one-hit mode, is a block of its own with either size, so the sequence after it
//!   starts a block with and without the switch);
//! - the round's results are put into the hit lists in ordinal order (blast_engine.c:1553-1554;
//!   the sequences without HSP lists change no list), so every hit list receives the one-thread
//!   order's updates. No floating-point operation changes order: each sequence's kernel and each
//!   update are the one-thread ones;
//! - a search that fails returns the first error of the one-thread order (`x_one_thread_error`):
//!   with another cut a damaged database can meet its first error elsewhere (the byte range of a
//!   read past the end of the file names the block, and a read error and a rejected sequence of one
//!   1 MiB block come in the other order), so the blocks of `BLOCK_BYTES` are read and prepared again
//!   in order, until the first error.
//!
//! Memory: a worker holds one block, at most the cut (or one longer sequence), and a round covers at
//! most `BLOCKS_PER_THREAD` x threads blocks, so the sequence memory does not grow. More threads search
//! at once, each with its search scratch (diagonal table, offset pairs, alignment state), as in the
//! `-subject` path: exact blastp of 3,000 human proteins against `BlastDBSpecies6` at 16 threads
//! peaked at 42.9-44.1 MiB against 37.2 MiB without the switch (workstation, perf-mb-20261014
//! `logs/dbcut_report.md`).
//!
//! No NCBI counterpart for the cut itself: NCBI's search threads take chunks of
//! MAX(sequences / 100, 1) ordinal ids from one iterator (`BlastSeqSrcIteratorNewEx`,
//! blast_engine.c:1409; seqdbimpl.cpp:328-350, quoted above `SeqDb::next_block`), whatever the
//! thread count; LOSAT's blocks are cut by bytes, because one block is one read of one file.
//!
//! `LOSAT_X_BPDBCUT` is a switch of the strict set (on by default since 2026-10-17): on unless it is `0`.
//! It acts for every task of BLASTP with `-db` (blastp, blastp-fast, blastp-short), the only program whose search
//! reaches this driver.

use std::sync::OnceLock;

use anyhow::Result;

use crate::seqdb::SeqDb;
use crate::utils::threading::SearchPool;

use super::{
    read_search_block, search_blocks_of, BlastpDatabase, BlastpHitList, BlastpSubjectRef,
    BlastpSubjectResult, BlastpSubjectScratch, PrelimSink, XPrelimHsp, BLOCK_BYTES,
};

// No NCBI counterpart: reads the LOSAT_X_BPDBCUT switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDBCUT): a strict-set switch, on unless it is `0`.
pub(super) fn x_bp_db_cut() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDBCUT"))
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// Blocks per thread that the cut aims at: a block costs a task and a scratch, and four per thread
/// let the threads end their round together. Measured on `BlastDBSpecies6` at 16 threads: four per
/// thread 7.91 s, sixteen 8.70 s (perf-mb-20261014 `logs/split2_report.md`).
pub(super) const X_BLOCKS_PER_THREAD: u64 = 4;

// No NCBI counterpart: as `X_BLOCKS_PER_THREAD`.
/// The smallest block of the cut, so that a small database is not cut into blocks of a few
/// sequences.
pub(super) const X_MIN_BLOCK_BYTES: u64 = 4 << 10;

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The block size of a search with `threads` threads of a database of `sequence_bytes`:
/// `BLOCK_BYTES` at one thread, else `sequence_bytes / (X_BLOCKS_PER_THREAD x threads)` between
/// `X_MIN_BLOCK_BYTES` and `BLOCK_BYTES`.
pub(super) fn x_cut_block_bytes(sequence_bytes: u64, threads: usize) -> u64 {
    if threads <= 1 {
        return BLOCK_BYTES;
    }
    let parts = X_BLOCKS_PER_THREAD.saturating_mul(threads as u64);
    (sequence_bytes / parts).clamp(X_MIN_BLOCK_BYTES, BLOCK_BYTES)
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The bytes of the database's sequences as a one-volume sequence file holds them: the letters, a
/// sentinel byte after each sequence and one before the first (each further volume has one more).
pub(super) fn x_sequence_bytes(database: &BlastpDatabase) -> u64 {
    database
        .total_letters()
        .saturating_add(u64::from(database.num_sequences()))
        .saturating_add(1)
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The block size of the search of `database` on `pool`: the cut when the pool searches the rounds
/// with more than one thread, else `BLOCK_BYTES`.
pub(super) fn x_db_cut_block_bytes(database: &BlastpDatabase, pool: &SearchPool<'_>) -> u64 {
    let threads = if pool.enabled() { pool.threads() } else { 1 };
    let cut = x_cut_block_bytes(x_sequence_bytes(database), threads);
    // EXPERIMENT (LOSAT_X_BPDBKEEP): blocks that a worker's kept buffer holds below glibc's mmap
    // threshold (`x_db_keep.rs`).
    if crate::seqdb::x_bp_db_small_blocks() {
        cut.min(super::x_db_keep::X_KEEP_BLOCK_BYTES)
    } else {
        cut
    }
}

// No NCBI counterpart: counts blocks for the unit tests; it changes no value NCBI computes.
/// The number of blocks of at most `block_bytes` that `search_blocks_of` reads.
#[cfg(test)]
pub(super) fn x_count_blocks(db: &SeqDb, block_bytes: u64) -> u64 {
    let mut blocks = 0;
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        blocks += 1;
        next_ordinal = range.end_ordinal();
    }
    blocks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
// ```c
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// The one-thread loop of `search_blocks_of` reads its blocks of `BLOCK_BYTES` in ordinal order and
// stops at the first that cannot be read or prepared (`read_search_block`); the kernel itself
// returns no error. The parallel rounds of the same blocks return the same error (the round's first
// in ordinal order).
/// The error that the search of `db` with blocks of `BLOCK_BYTES` returns, found by reading and
/// preparing those blocks in order; `error` (the cut's) when every block reads (an error that does
/// not come back on a second read).
pub(super) fn x_one_thread_error(db: &SeqDb, error: anyhow::Error) -> anyhow::Error {
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, BLOCK_BYTES) {
        if let Err(first) = read_search_block(db, range) {
            return first;
        }
        next_ordinal = range.end_ordinal();
    }
    error
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1413
// ```c
//     db_length = BlastSeqSrcGetTotLen(seq_src);
//
//     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
// ```c
//           /* Save the results. */
//           status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-173
// ```c
//     if (ewp->diag_table) {
//         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//             ewp->diag_table->offset = ewp->diag_table->window;
//             s_BlastDiagClear(ewp->diag_table);
//         } else {
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//         }
// ```
/// EXPERIMENT (LOSAT_X_BPDBCUT): `search_blocks` with the blocks of the cut
/// (`x_db_cut_block_bytes`) and, on failure, the error of the one-thread order.
///
/// The output is that of the blocks of `BLOCK_BYTES` (the one-thread order's, apart from the
/// exception in the module header) for any cut: the cut only moves where blocks end, every
/// block starts at the one-thread order's diagonal offset of its first sequence from a cleared
/// table, and the round's results are merged in ordinal order (module header). No floating-point
/// operation changes order.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_blocks_cut<S, K, P>(
    database: &BlastpDatabase,
    pool: &SearchPool<'_>,
    window: i32,
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
    let db = &database.db;
    let block_bytes = x_db_cut_block_bytes(database, pool);
    crate::utils::xreach::x_reach!("LOSAT_X_BPDBCUT", block_bytes != BLOCK_BYTES);
    // EXPERIMENT (LOSAT_X_BPDBKEEP=3, a diagnostic): rounds of as many bytes as the cut's.
    let blocks_per_thread = if crate::seqdb::x_bp_db_byte_rounds() {
        let threads = if pool.enabled() { pool.threads() } else { 1 };
        let cut = x_cut_block_bytes(x_sequence_bytes(database), threads);
        blocks_per_thread * cut.div_ceil(block_bytes) as usize
    } else {
        blocks_per_thread
    };
    search_blocks_of(
        database,
        pool,
        window,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        process_subject,
        sink,
    )
    .map_err(|error| {
        if block_bytes == BLOCK_BYTES {
            error
        } else {
            x_one_thread_error(db, error)
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // The cut: one thread keeps `BLOCK_BYTES`; more threads a quarter of the bytes per thread,
    // between 4 KiB and 1 MiB.
    #[test]
    fn the_cut_is_a_quarter_of_the_bytes_per_thread_between_4_kib_and_1_mib() {
        for bytes in [0u64, 1, 812, 2_920_000, 64 << 20, u64::MAX] {
            assert_eq!(x_cut_block_bytes(bytes, 0), BLOCK_BYTES, "{bytes}");
            assert_eq!(x_cut_block_bytes(bytes, 1), BLOCK_BYTES, "{bytes}");
        }
        // BlastDBSpecies6 (2.92 MB) at 16 threads: 64 blocks of 45,625 bytes.
        assert_eq!(x_cut_block_bytes(2_920_000, 16), 45_625);
        assert_eq!(x_cut_block_bytes(2_920_000, 2), 365_000);
        // Small databases: the floor.
        assert_eq!(x_cut_block_bytes(0, 2), X_MIN_BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(812, 7), X_MIN_BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(8 * 4096 - 1, 2), X_MIN_BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(8 * 4096 + 8, 2), 4097);
        // From 4 x threads MiB on: `BLOCK_BYTES`.
        assert_eq!(x_cut_block_bytes((64 << 20) - 64, 16), BLOCK_BYTES - 1);
        assert_eq!(x_cut_block_bytes(64 << 20, 16), BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(u64::MAX, 7), BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(u64::MAX, usize::MAX), X_MIN_BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(1 << 40, usize::MAX), X_MIN_BLOCK_BYTES);
        // Never above one thread's blocks, never below the floor; smaller with more threads.
        for bytes in [0u64, 5_000, 300_000, 2_920_000, 40 << 20, 200 << 20] {
            let mut previous = BLOCK_BYTES;
            for threads in 1..=130 {
                let cut = x_cut_block_bytes(bytes, threads);
                assert!((X_MIN_BLOCK_BYTES..=BLOCK_BYTES).contains(&cut));
                assert!(cut <= previous, "{bytes} bytes, {threads} threads");
                previous = cut;
            }
        }
    }
}

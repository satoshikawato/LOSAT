//! EXPERIMENT (LOSAT_X_BPDBCUT / LOSAT_X_BPDBCUTSHADOW): the blocks of a `-db` search's
//! preliminary stage cut by the size of the database and the thread count.
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
//!   its ordinal id, as in NCBI's loop with one thread (blast_engine.c:1407-1426, quoted above
//!   `search_blocks`);
//! - every block starts from a cleared diagonal table at the offset that the one-thread order gives
//!   its first sequence (`prepare_independent_subject`, `diag_offset_after`, the recurrence of
//!   `Blast_ExtendWordExit`, blast_extend.c:160-172 quoted above `diag_offset_after`) and advances it
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
//! Memory: a worker holds one block, at most the cut, and a round covers at most
//! `BLOCKS_PER_THREAD` x threads blocks, so neither grows.
//!
//! No NCBI counterpart for the cut itself: NCBI's search threads take chunks of
//! MAX(sequences / 100, 1) ordinal ids from one iterator (`BlastSeqSrcIteratorNewEx`,
//! blast_engine.c:1407-1426; seqdbimpl.cpp:328-350, quoted above `SeqDb::next_block`), whatever the
//! thread count; LOSAT's blocks are cut by bytes, because one block is one read of one file.
//!
//! Off unless `LOSAT_X_BPDBCUT` is set (to any value); not in the strict set. It acts for every
//! task of BLASTP with `-db` (blastp, blastp-fast), the only program whose search reaches this
//! driver.
//!
//! `LOSAT_X_BPDBCUTSHADOW=1` (shadow mode, implies the switch): after each search with the cut, the
//! search runs again as without the switch (`search_blocks_of` with `BLOCK_BYTES`, the same pool,
//! rounds and kernel) into separate hit lists, and the two are compared: every preliminary and
//! compact hit list of every query, every field (their `Debug` text), or the texts of the two
//! errors. The first difference stops the run (panic) after printing the line that is otherwise
//! printed at exit: `[X_SHADOW] BPDBCUT searches= blocks= blocks_ref= lists_compared= identical`
//! (`blocks`: the blocks of the cut; `blocks_ref`: the blocks of `BLOCK_BYTES`; `lists_compared`:
//! hit lists compared, preliminary and compact, that are not empty on either side).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use anyhow::Result;

use crate::seqdb::SeqDb;
use crate::utils::threading::SearchPool;

use super::{
    read_search_block, search_blocks_of, BlastpDatabase, BlastpHitList, BlastpSubjectRef,
    BlastpSubjectResult, BlastpSubjectScratch, XPrelimHsp, BLOCK_BYTES,
};

/// EXPERIMENT (LOSAT_X_BPDBCUT): the mode the variables select.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XDbCut {
    Off,
    On,
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BPDBCUT / LOSAT_X_BPDBCUTSHADOW switches once; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPDBCUT): off unless set to any value; `LOSAT_X_BPDBCUTSHADOW` (any value)
/// selects the shadow mode.
pub(super) fn x_bp_db_cut() -> XDbCut {
    static MODE: OnceLock<XDbCut> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPDBCUTSHADOW").is_some() {
            XDbCut::Shadow
        } else if std::env::var_os("LOSAT_X_BPDBCUT").is_some() {
            XDbCut::On
        } else {
            XDbCut::Off
        }
    })
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
    x_cut_block_bytes(x_sequence_bytes(database), threads)
}

// No NCBI counterpart: counts blocks for the shadow's line; it changes no value NCBI computes.
/// The number of blocks of at most `block_bytes` that `search_blocks_of` reads.
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

static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_BLOCKS: AtomicU64 = AtomicU64::new(0);
static X_BLOCKS_REF: AtomicU64 = AtomicU64::new(0);
static X_LISTS_COMPARED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BPDBCUT searches={} blocks={} blocks_ref={} lists_compared={}{}",
        X_SEARCHES.load(Ordering::Relaxed),
        X_BLOCKS.load(Ordering::Relaxed),
        X_BLOCKS_REF.load(Ordering::Relaxed),
        X_LISTS_COMPARED.load(Ordering::Relaxed),
        if identical {
            " identical"
        } else {
            " DIFFERENT"
        },
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPDBCUTSHADOW): the line `main` prints at exit.
pub(super) fn x_bp_db_cut_print_summary() {
    if x_bp_db_cut() == XDbCut::Shadow {
        x_shadow_print_line(true);
    }
}

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    // The panic ends the run before `main` prints the line at exit: it is printed here.
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BPDBCUT {message}");
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// Compares the hit lists of the cut (`lists`) with those of `BLOCK_BYTES` (`expected`) by their
/// `Debug` text, one query at a time; returns the lists compared (those not empty on either side).
fn x_compare_lists<H: std::fmt::Debug>(
    kind: &str,
    lists: &[Option<BlastpHitList<H>>],
    expected: &[Option<BlastpHitList<H>>],
) -> u64 {
    if lists.len() != expected.len() {
        x_shadow_difference(format!(
            "{kind} hit lists: {} queries, {} without the cut",
            lists.len(),
            expected.len()
        ));
    }
    let mut compared = 0;
    for (query, (list, expected)) in lists.iter().zip(expected).enumerate() {
        if list.is_none() && expected.is_none() {
            continue;
        }
        let (text, expected_text) = (format!("{list:?}"), format!("{expected:?}"));
        if text != expected_text {
            x_shadow_difference(format!(
                "{kind} hit list of query {query} differs from the one without the cut: {text} / {expected_text}"
            ));
        }
        compared += 1;
    }
    compared
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:160-172
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
/// (`x_db_cut_block_bytes`) and, on failure, the error of the one-thread order; in `mode`
/// `Shadow`, the search again with `BLOCK_BYTES` and the comparison of the two (module header).
///
/// The output is the one-thread order's for any cut: the cut only moves where blocks end, every
/// block starts at the one-thread order's diagonal offset of its first sequence from a cleared
/// table, and the round's results are merged in ordinal order (module header). No floating-point
/// operation changes order.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_blocks_cut<S, K>(
    mode: XDbCut,
    database: &BlastpDatabase,
    pool: &SearchPool<'_>,
    window: i32,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    prelim_hitlist_size: usize,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let db = &database.db;
    let block_bytes = x_db_cut_block_bytes(database, pool);
    let shadow = mode == XDbCut::Shadow;
    if shadow
        && !(preliminary_hit_lists.iter().all(Option::is_none)
            && x_compact_hit_lists.iter().all(Option::is_none))
    {
        x_shadow_difference("the hit lists are not empty before the search".to_string());
    }
    let result = search_blocks_of(
        database,
        pool,
        window,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        process_subject,
        preliminary_hit_lists,
        x_compact_hit_lists,
        prelim_hitlist_size,
    )
    .map_err(|error| {
        if block_bytes == BLOCK_BYTES {
            error
        } else {
            x_one_thread_error(db, error)
        }
    });
    if !shadow {
        return result;
    }
    let mut expected_lists: Vec<Option<BlastpHitList>> = std::iter::repeat_with(|| None)
        .take(preliminary_hit_lists.len())
        .collect();
    let mut expected_compact: Vec<Option<BlastpHitList<XPrelimHsp>>> =
        std::iter::repeat_with(|| None)
            .take(x_compact_hit_lists.len())
            .collect();
    let expected = search_blocks_of(
        database,
        pool,
        window,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        process_subject,
        &mut expected_lists,
        &mut expected_compact,
        prelim_hitlist_size,
    );
    let compared = match (&result, &expected) {
        (Ok(()), Ok(())) => {
            x_compare_lists("preliminary", preliminary_hit_lists, &expected_lists)
                + x_compare_lists("compact", x_compact_hit_lists, &expected_compact)
        }
        (Err(error), Err(expected)) => {
            let (text, expected_text) = (format!("{error:#}"), format!("{expected:#}"));
            if text != expected_text {
                x_shadow_difference(format!(
                    "the search's error differs from the one without the cut: {text} / {expected_text}"
                ));
            }
            0
        }
        (Ok(()), Err(expected)) => x_shadow_difference(format!(
            "the search succeeded; without the cut it failed: {expected:#}"
        )),
        (Err(error), Ok(())) => x_shadow_difference(format!(
            "the search failed, without the cut it succeeded: {error:#}"
        )),
    };
    X_SEARCHES.fetch_add(1, Ordering::Relaxed);
    X_BLOCKS.fetch_add(x_count_blocks(db, block_bytes), Ordering::Relaxed);
    X_BLOCKS_REF.fetch_add(x_count_blocks(db, BLOCK_BYTES), Ordering::Relaxed);
    X_LISTS_COMPARED.fetch_add(compared, Ordering::Relaxed);
    result
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

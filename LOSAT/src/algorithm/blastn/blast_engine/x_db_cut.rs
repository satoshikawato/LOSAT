//! EXPERIMENT (LOSAT_X_BNDBCUT / LOSAT_X_BNDBCUTSHADOW): the blocks of a `-db` search's
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
//!   its ordinal id, as in NCBI's loop with one thread (blast_engine.c:1407-1426, quoted above
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
//! MAX(sequences / 100, 1) ordinal ids from one iterator (`BlastSeqSrcIteratorNewEx`,
//! blast_engine.c:1409; seqdbimpl.cpp:328-350, quoted above `SeqDb::next_block`), whatever the
//! thread count; LOSAT's blocks are cut by bytes, because one block is one read of one file.
//!
//! With LOSAT_X_BNDBSTREAM (`x_db_stream.rs`, in the strict set) the parallel stage has no rounds:
//! the threads claim the cut's blocks in ordinal order and the results are handed over in that
//! order; the shadow's second search with blocks of 1 MiB runs the same way.
//!
//! `LOSAT_X_BNDBCUT` is a switch of the strict set: on unless it is `0`, or `LOSAT_X_DEFAULT=0`
//! and it is not set. It acts for every BLASTN task with `-db` (megablast, blastn, dc-megablast,
//! blastn-short).
//!
//! `LOSAT_X_BNDBCUTSHADOW=1` (shadow mode, implies the switch): each `-db` search (at one thread too,
//! where both have blocks of `BLOCK_BYTES`) keeps what it hands to the merge, holding the kept
//! results of the whole database, the search runs again as without the switch (`search_blocks_of`
//! with `BLOCK_BYTES`, the same pool, rounds and kernel), and the two are compared: the ordinal ids in
//! order and each kept result (its `Debug` text) and the good initial extensions each counted, or
//! the texts of the two errors. The count of good initial extensions, which sizes the next query
//! batch (NCBI's batch size mixer), is put back to what the search with the cut left, and the kept
//! results go to the merge after the comparison: the run continues as without the shadow (only
//! `LOSAT_TIMING` and `LOSAT_X_STATS` count the preliminary stage twice). The first difference
//! stops the run (panic) after printing the line that is otherwise printed at exit:
//! `[X_SHADOW] BNDBCUT searches= blocks= blocks_ref= merged= identical` (`blocks`: the blocks of
//! the cut; `blocks_ref`: the blocks of `BLOCK_BYTES`; `merged`: kept results compared).

use std::fmt::Debug;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use anyhow::Result;

use crate::seqdb::SeqDb;
use crate::utils::threading::SearchPool;

use super::{
    block_subject, search_blocks_of, search_read_error, BlastnDatabase, BlastnSubjectRef,
    BLOCK_BYTES,
};

/// EXPERIMENT (LOSAT_X_BNDBCUT): the mode the variables select.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XDbCut {
    Off,
    On,
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BNDBCUT / LOSAT_X_BNDBCUTSHADOW switches once; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BNDBCUT): a strict-set switch, on unless it is `0`, or `LOSAT_X_DEFAULT=0`
/// and it is not set; `LOSAT_X_BNDBCUTSHADOW` (any value) selects the shadow mode.
pub(super) fn x_bn_db_cut() -> XDbCut {
    static MODE: OnceLock<XDbCut> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BNDBCUTSHADOW").is_some() {
            XDbCut::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BNDBCUT") {
            XDbCut::On
        } else {
            XDbCut::Off
        }
    })
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// Blocks per thread that the cut aims at: a block costs a task, and four per thread let the
/// threads end their round together (BLASTP's measured choice, `LOSAT_X_BPDBCUT`).
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
/// The bytes of the database's sequences in its sequence files without their ambiguity data: four
/// letters per byte and the last byte of each sequence, which holds the count of its letters
/// (seqdbvol.cpp:1805-1838, quoted above `block_subject`). A sequence of `n` letters takes
/// `n / 4 + 1` bytes, so this is the sum up to 3 / 4 byte per sequence.
pub(super) fn x_sequence_bytes(db: &SeqDb) -> u64 {
    (db.total_length() / 4).saturating_add(u64::from(db.num_seqs()))
}

// No NCBI counterpart: the size of the blocks (scheduling); it changes no value NCBI computes.
/// The block size of the search of `db` on `pool`: the cut when the preliminary stage runs its
/// rounds on the pool (`parallel`), else `BLOCK_BYTES`.
pub(super) fn x_db_cut_block_bytes(db: &SeqDb, pool: &SearchPool<'_>, parallel: bool) -> u64 {
    let threads = if parallel { pool.threads() } else { 1 };
    x_cut_block_bytes(x_sequence_bytes(db), threads)
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
// ```c
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
//        Int4 stat_length;
//        if (seq_arg.oid == BLAST_SEQSRC_ERROR) {
//            status = BLASTERR_SEQSRC;
//            break;
//        }
// ...
//        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//            continue;
// ```
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
            if let Err(first) = block_subject(&block, i) {
                return first;
            }
        }
        next_ordinal = range.end_ordinal();
    }
    error
}

static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_BLOCKS: AtomicU64 = AtomicU64::new(0);
static X_BLOCKS_REF: AtomicU64 = AtomicU64::new(0);
static X_MERGED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BNDBCUT searches={} blocks={} blocks_ref={} merged={}{}",
        X_SEARCHES.load(Ordering::Relaxed),
        X_BLOCKS.load(Ordering::Relaxed),
        X_BLOCKS_REF.load(Ordering::Relaxed),
        X_MERGED.load(Ordering::Relaxed),
        if identical {
            " identical"
        } else {
            " DIFFERENT"
        },
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BNDBCUTSHADOW): the line `main` prints at exit.
pub(super) fn x_bn_db_cut_print_summary() {
    if x_bn_db_cut() == XDbCut::Shadow {
        x_shadow_print_line(true);
    }
}

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    // The panic ends the run before `main` prints the line at exit: it is printed here.
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BNDBCUT {message}");
}

// No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
/// Compares what the search with the cut handed to the merge (`found`) with what the search of
/// `BLOCK_BYTES` handed to it (`expected`, the `Debug` texts); returns the results compared.
fn x_compare_merged<R: Debug>(found: &[(u32, R)], expected: &[(u32, String)]) -> u64 {
    if found.len() != expected.len() {
        x_shadow_difference(format!(
            "{} results merged, {} without the cut",
            found.len(),
            expected.len()
        ));
    }
    for (index, ((ordinal, result), (expected_ordinal, expected_text))) in
        found.iter().zip(expected).enumerate()
    {
        if ordinal != expected_ordinal {
            x_shadow_difference(format!(
                "result {index} merged for sequence {ordinal}, without the cut for sequence {expected_ordinal}"
            ));
        }
        let text = format!("{result:?}");
        if &text != expected_text {
            x_shadow_difference(format!(
                "the result of sequence {ordinal} differs from the one without the cut: {text} / {expected_text}"
            ));
        }
    }
    found.len() as u64
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
/// EXPERIMENT (LOSAT_X_BNDBCUT): `search_blocks` with the blocks of the cut
/// (`x_db_cut_block_bytes`) and, on failure, the error of the one-thread order; in `mode`
/// `Shadow`, the search again with `BLOCK_BYTES` and the comparison of the two (module header).
///
/// The output is that of the blocks of `BLOCK_BYTES` (the one-thread order's) for any cut: the cut
/// only moves where blocks end, a sequence's result does not depend on what the worker's scratch
/// searched before it, and the round's results are merged in ordinal order (module header). No
/// floating-point operation changes order. `extends`: the count of good initial extensions that
/// `search` adds to; the shadow mode leaves it as the search with the cut leaves it.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_blocks_cut<T, R, I, K, M>(
    mode: XDbCut,
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    blocks_per_thread: usize,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
    extends: &AtomicU64,
) -> Result<()>
where
    R: Send + Debug,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    let db = &database.db;
    let block_bytes = x_db_cut_block_bytes(db, pool, parallel);
    let one_thread_error = |error| {
        if block_bytes == BLOCK_BYTES {
            error
        } else {
            x_one_thread_error(db, error)
        }
    };
    if mode != XDbCut::Shadow {
        return search_blocks_of(
            database,
            pool,
            parallel,
            block_bytes,
            blocks_per_thread,
            new_scratch,
            search,
            merge,
            extends,
        )
        .map_err(one_thread_error);
    }
    let mut found: Vec<(u32, R)> = Vec::new();
    let before = extends.load(Ordering::Relaxed);
    let result = search_blocks_of(
        database,
        pool,
        parallel,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        search,
        &mut |ordinal, result| found.push((ordinal, result)),
        extends,
    )
    .map_err(one_thread_error);
    let counted = extends.load(Ordering::Relaxed);
    let mut expected: Vec<(u32, String)> = Vec::new();
    let expected_result = search_blocks_of(
        database,
        pool,
        parallel,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        search,
        &mut |ordinal, result| expected.push((ordinal, format!("{result:?}"))),
        extends,
    );
    // The second search's extensions do not count: the next query batch is sized as without the
    // shadow.
    let counted_ref = extends.load(Ordering::Relaxed).wrapping_sub(counted);
    extends.store(counted, Ordering::Relaxed);
    let compared = match (&result, &expected_result) {
        (Ok(()), Ok(())) => {
            let counted_cut = counted.wrapping_sub(before);
            if counted_cut != counted_ref {
                x_shadow_difference(format!(
                    "{counted_cut} good initial extensions counted, {counted_ref} without the cut"
                ));
            }
            x_compare_merged(&found, &expected)
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
    X_MERGED.fetch_add(compared, Ordering::Relaxed);
    result?;
    for (ordinal, kept) in found {
        merge(ordinal, kept);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open(name: &str) -> BlastnDatabase {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/db")
            .join(name);
        BlastnDatabase::open(
            path.to_str().expect("fixture path is UTF-8"),
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    // The shadow mode runs the kernel twice: the count of good initial extensions that the kernel
    // adds to (it sizes the next query batch) is left as one search leaves it, and the kept results
    // are merged once, in ordinal order, as with the switch on.
    #[test]
    fn the_shadow_leaves_the_extension_count_of_one_search_and_merges_once() {
        let database = open("nucl_mv");
        let sequences = u64::from(database.db.num_seqs());
        assert_eq!(sequences, 12);
        for mode in [XDbCut::On, XDbCut::Shadow] {
            let extends = AtomicU64::new(5);
            let mut merged = Vec::new();
            crate::utils::threading::with_search_pool(2, "blastn", |pool| {
                x_search_blocks_cut(
                    mode,
                    &database,
                    pool,
                    true,
                    2,
                    &|| 0u32,
                    &|calls: &mut u32, ordinal: u32, _subject: BlastnSubjectRef<'_>| {
                        *calls += 1;
                        extends.fetch_add(2, Ordering::Relaxed);
                        (ordinal % 3 == 0).then_some(ordinal)
                    },
                    &mut |ordinal, kept| {
                        assert_eq!(ordinal, kept);
                        merged.push(ordinal);
                    },
                    &extends,
                )
            })
            .unwrap();
            assert_eq!(
                extends.load(Ordering::Relaxed),
                5 + 2 * sequences,
                "{mode:?}"
            );
            assert_eq!(merged, [0, 3, 6, 9], "{mode:?}");
        }
    }

    // The cut: one thread keeps `BLOCK_BYTES`; more threads a quarter of the bytes per thread,
    // between 4 KiB and 1 MiB.
    #[test]
    fn the_cut_is_a_quarter_of_the_bytes_per_thread_between_4_kib_and_1_mib() {
        for bytes in [0u64, 1, 812, 10_300_000, 64 << 20, u64::MAX] {
            assert_eq!(x_cut_block_bytes(bytes, 0), BLOCK_BYTES, "{bytes}");
            assert_eq!(x_cut_block_bytes(bytes, 1), BLOCK_BYTES, "{bytes}");
        }
        // 28,036 16S rRNA genes (10.3 MB) at 16 threads: 64 blocks of 160,937 bytes; their first
        // 2,000 (0.72 MB) at 16 threads: blocks of 11,250 bytes.
        assert_eq!(x_cut_block_bytes(10_300_000, 16), 160_937);
        assert_eq!(x_cut_block_bytes(720_000, 16), 11_250);
        assert_eq!(x_cut_block_bytes(720_000, 2), 90_000);
        // Small databases: the floor.
        assert_eq!(x_cut_block_bytes(0, 2), X_MIN_BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(8 * 4096 - 1, 2), X_MIN_BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(8 * 4096 + 8, 2), 4097);
        // From 4 x threads MiB on: `BLOCK_BYTES`.
        assert_eq!(x_cut_block_bytes((64 << 20) - 64, 16), BLOCK_BYTES - 1);
        assert_eq!(x_cut_block_bytes(64 << 20, 16), BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(u64::MAX, 7), BLOCK_BYTES);
        assert_eq!(x_cut_block_bytes(u64::MAX, usize::MAX), X_MIN_BLOCK_BYTES);
        // Never above one thread's blocks, never below the floor; smaller with more threads.
        for bytes in [0u64, 5_000, 300_000, 10_300_000, 40 << 20, 200 << 20] {
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

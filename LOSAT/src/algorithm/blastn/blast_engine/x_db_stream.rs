//! EXPERIMENT (LOSAT_X_BNDBSTREAM / LOSAT_X_BNDBSTREAMSHADOW): the parallel preliminary stage of
//! a `-db` search without rounds. The pool threads claim the blocks once each in ordinal order,
//! and one of them (the owner) hands each block's results to the hit lists in that order while
//! the others go on with the next blocks (`x_bnstream::search_in_order`, the driver of the
//! `-subject` path's LOSAT_X_BNSTREAM). A thread claims a block only while it is fewer than
//! `X_BNDBSTREAM_BLOCKS_PER_THREAD` x threads blocks past the last one handed over, so that at
//! most that many blocks' results wait for the owner.
//!
//! Without the switch, `search_blocks_of` searches rounds of `BLOCKS_PER_THREAD` x threads blocks
//! (`search_rounds`) and merges a round's results when the whole round is searched, so the
//! preliminary HSP lists of every sequence of the round are in memory together. A database of
//! fewer than `BLOCKS_PER_THREAD` x threads blocks is one round; with LOSAT_X_BNDBCUT (about 4
//! blocks per thread), every database below 8 x 4 x threads MiB is. Megablast of 100 16S queries
//! against the 28,036 16S rRNA genes on the workstation, one round from 2 threads on: at 16
//! threads the round held 69 MiB of preliminary HSP lists (the lookup table is 16 MiB), and the
//! arenas kept the freed space to the end of the run: max RSS 34 MiB at 1 thread, 108 MiB at 2,
//! 126 MiB at 16 (NCBI `-db` at 16 threads: 109 MiB with outfmt 8). v4-B session 10, probe of the
//! allocations by thread and region (`v4-B/results/b14probe`).
//!
//! Why the reports do not change (`search_blocks`'s comment holds for the order of the merges):
//!
//! - every block is claimed once, in ordinal order, and searched once, with the blocks cut as
//!   before (`block_bytes`), so every sequence is searched once;
//! - the owner hands the blocks' results over in ordinal order, and a block's results in the
//!   order of its sequences (blast_engine.c:1553-1554), so the hit lists receive the one-thread
//!   order's updates. No floating-point operation changes order;
//! - the result of a sequence does not depend on what the thread's scratch searched before it
//!   (the rounds' `map_init` already hands a scratch runs of blocks that change from one execution
//!   to the next); here a thread makes one scratch (`new_scratch`) when it claims its first block
//!   and keeps it for the blocks it claims after, in increasing order;
//! - a block that fails stops the hand-over: the first error in ordinal order is returned after
//!   the results of every block before it are merged (as a round returns the first error of the
//!   round after merging the blocks before it), and blocks that are not claimed yet are not
//!   searched. The caller maps the error as before (`x_one_thread_error` when the blocks are cut);
//! - the count of good initial extensions (the next query batch's size) adds up over every block
//!   searched once, as in the rounds; a failed search ends the run.
//!
//! Memory: at most the window's blocks of results wait for the owner, and each thread holds one
//! block and one scratch (the rounds make a scratch per `map_init` run, a few hundred per thread
//! in the case above). The blocks are listed before the search: one `BlockRange` (16 bytes) per
//! block, 1.6 MB for a database of 100 GB of sequence bytes in blocks of 1 MiB.
//!
//! No NCBI counterpart for the window itself: NCBI's search threads take chunks of ordinal ids
//! from one iterator and write each sequence's HSP list to the HSP stream when it is done
//! (blast_engine.c:1409, 1553-1554, quoted above `search_blocks`); the hand-over in ordinal order
//! keeps the order of its one-thread loop.
//!
//! `LOSAT_X_BNDBSTREAM` is a switch of the strict set: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. It acts for every BLASTN task with `-db` at more than
//! one thread; at one thread the blocks are searched one after another, as without it.
//!
//! `LOSAT_X_BNDBSTREAMSHADOW=1` (shadow mode, implies the switch): each parallel preliminary stage
//! runs twice, with the window and with the rounds, and compares every result handed to the hit
//! lists (ordinal and `Debug` text, in order), the errors' texts and the counts of good initial
//! extensions; the count is left as one search leaves it and the results are merged once. A
//! difference prints `[X_SHADOW] BNDBSTREAM ... DIFFERENT` and stops the run (panic); otherwise
//! `[X_SHADOW] BNDBSTREAM searches=<parallel stages> blocks=<blocks> merged=<results compared>
//! identical` at exit.

use std::fmt::Debug;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::OnceLock;

use anyhow::Result;

use crate::seqdb::{BlockRange, SeqDb};
use crate::utils::threading::SearchPool;

use super::{block_subject, search_read_error, search_rounds, BlastnSubjectRef};

/// No NCBI counterpart: the modes of LOSAT_X_BNDBSTREAM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum XDbStream {
    Off,
    On,
    Shadow,
}

/// No NCBI counterpart (scheduling only): the mode of LOSAT_X_BNDBSTREAM.
pub(super) fn x_bn_db_stream() -> XDbStream {
    static MODE: OnceLock<XDbStream> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BNDBSTREAMSHADOW").is_some() {
            XDbStream::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BNDBSTREAM") {
            XDbStream::On
        } else {
            XDbStream::Off
        }
    })
}

// No NCBI counterpart (scheduling only): the window; it changes no value NCBI computes.
/// Blocks per thread of the window: a thread finds a block to claim while the owner waits for
/// the oldest one.
const X_BNDBSTREAM_BLOCKS_PER_THREAD: usize = 2;

/// The blocks of at most `block_bytes` of `db`, in ordinal order.
fn blocks_of(db: &SeqDb, block_bytes: u64) -> Vec<BlockRange> {
    let mut blocks = Vec::new();
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        blocks.push(range);
        next_ordinal = range.end_ordinal();
    }
    blocks
}

/// The preliminary stage of `search_blocks_of` at more than one thread with LOSAT_X_BNDBSTREAM
/// (`mode` On or Shadow): the blocks of at most `block_bytes` searched by the threads of `pool`
/// and their results passed to `merge` in ordinal order (module documentation). The shadow
/// mode compares the results with those of the rounds of `blocks_per_thread` blocks per thread.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_search_blocks_stream<T, R, I, K, M>(
    mode: XDbStream,
    db: &SeqDb,
    pool: &SearchPool<'_>,
    block_bytes: u64,
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
    let blocks = blocks_of(db, block_bytes);
    if mode != XDbStream::Shadow {
        return search_window(db, pool, &blocks, new_scratch, search, merge);
    }
    let mut found: Vec<(u32, R)> = Vec::new();
    let before = extends.load(Ordering::Relaxed);
    let result = search_window(db, pool, &blocks, new_scratch, search, &mut |ordinal, kept| {
        found.push((ordinal, kept))
    });
    let counted = extends.load(Ordering::Relaxed);
    let mut expected: Vec<(u32, String)> = Vec::new();
    let expected_result = search_rounds(
        db,
        pool,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        search,
        &mut |ordinal, kept| expected.push((ordinal, format!("{kept:?}"))),
    );
    // The rounds searched the database once more: leave the count as one search leaves it.
    let counted_ref = extends.load(Ordering::Relaxed).wrapping_sub(counted);
    extends.store(counted, Ordering::Relaxed);
    let compared = match (&result, &expected_result) {
        (Ok(()), Ok(())) => {
            let counted_window = counted.wrapping_sub(before);
            if counted_window != counted_ref {
                x_shadow_difference(format!(
                    "{counted_window} good initial extensions counted, {counted_ref} in rounds"
                ));
            }
            x_compare_merged(&found, &expected)
        }
        (Err(error), Err(expected)) => {
            let (text, expected_text) = (format!("{error:#}"), format!("{expected:#}"));
            if text != expected_text {
                x_shadow_difference(format!(
                    "the search's error differs from the one in rounds: {text} / {expected_text}"
                ));
            }
            0
        }
        (Ok(()), Err(expected)) => x_shadow_difference(format!(
            "the search succeeded; in rounds it failed: {expected:#}"
        )),
        (Err(error), Ok(())) => x_shadow_difference(format!(
            "the search failed, in rounds it succeeded: {error:#}"
        )),
    };
    X_SEARCHES.fetch_add(1, Ordering::Relaxed);
    X_BLOCKS.fetch_add(blocks.len() as u64, Ordering::Relaxed);
    X_MERGED.fetch_add(compared, Ordering::Relaxed);
    result?;
    for (ordinal, kept) in found {
        merge(ordinal, kept);
    }
    Ok(())
}

/// The window: `blocks` claimed by the threads of `pool` in order, and each block's results
/// handed to `merge` in ordinal order by the owner; the first error in ordinal order is returned
/// after the blocks before it are merged.
fn search_window<T, R, I, K, M>(
    db: &SeqDb,
    pool: &SearchPool<'_>,
    blocks: &[BlockRange],
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    let stop = AtomicBool::new(false);
    let search_block = |scratch: &mut T, i: usize| -> Result<Vec<(u32, R)>> {
        if stop.load(Ordering::Relaxed) {
            return Ok(Vec::new());
        }
        let range = blocks[i];
        let block = db.read_block(range).map_err(search_read_error)?;
        let mut kept = Vec::new();
        for j in 0..block.len() {
            let ordinal = range.first_ordinal + j as u32;
            if let Some(result) = search(scratch, ordinal, block_subject(&block, j)?) {
                kept.push((ordinal, result));
            }
        }
        Ok(kept)
    };
    let window = pool
        .threads()
        .max(1)
        .saturating_mul(X_BNDBSTREAM_BLOCKS_PER_THREAD);
    let mut first_error: Option<anyhow::Error> = None;
    let mut hand_over = |_: usize, found: Result<Vec<(u32, R)>>| {
        if first_error.is_some() {
            return;
        }
        match found {
            Ok(kept) => {
                for (ordinal, result) in kept {
                    merge(ordinal, result);
                }
            }
            Err(error) => {
                // The blocks not claimed yet are not searched; those claimed are dropped.
                stop.store(true, Ordering::Relaxed);
                first_error = Some(error);
            }
        }
    };
    pool.install(|| {
        super::super::x_bnstream::search_in_order(
            blocks.len(),
            window,
            new_scratch,
            &search_block,
            &mut hand_over,
        )
    });
    match first_error {
        Some(error) => Err(error),
        None => Ok(()),
    }
}

static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_BLOCKS: AtomicU64 = AtomicU64::new(0);
static X_MERGED: AtomicU64 = AtomicU64::new(0);

fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BNDBSTREAM searches={} blocks={} merged={}{}",
        X_SEARCHES.load(Ordering::Relaxed),
        X_BLOCKS.load(Ordering::Relaxed),
        X_MERGED.load(Ordering::Relaxed),
        if identical {
            " identical"
        } else {
            " DIFFERENT"
        },
    );
}

/// The shadow mode's line at exit.
pub(super) fn x_bn_db_stream_print_summary() {
    if x_bn_db_stream() == XDbStream::Shadow {
        x_shadow_print_line(true);
    }
}

fn x_shadow_difference(message: String) -> ! {
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BNDBSTREAM {message}");
}

fn x_compare_merged<R: Debug>(found: &[(u32, R)], expected: &[(u32, String)]) -> u64 {
    if found.len() != expected.len() {
        x_shadow_difference(format!(
            "{} results merged, {} in rounds",
            found.len(),
            expected.len()
        ));
    }
    for (index, ((ordinal, result), (expected_ordinal, expected_text))) in
        found.iter().zip(expected).enumerate()
    {
        if ordinal != expected_ordinal {
            x_shadow_difference(format!(
                "result {index} merged for sequence {ordinal}, in rounds for sequence {expected_ordinal}"
            ));
        }
        let text = format!("{result:?}");
        if &text != expected_text {
            x_shadow_difference(format!(
                "the result of sequence {ordinal} differs from the one in rounds: {text} / {expected_text}"
            ));
        }
    }
    found.len() as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::BlastnDatabase;

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

    // Uneven work per sequence, so that the blocks end out of order.
    fn work(ordinal: u32) -> u64 {
        let rounds = if ordinal % 5 == 1 { 200_000 } else { 100 };
        (0..rounds).fold(u64::from(ordinal), |acc, x| {
            acc.wrapping_mul(31).wrapping_add(x)
        })
    }

    // Every sequence is searched once, the kept results reach `merge` once each in ordinal
    // order, the count of good initial extensions is one search's, and a thread makes one scratch
    // at most: at 1 to 8 threads, blocks of one sequence (1 byte) to the whole database, with the
    // window and in the shadow mode.
    #[test]
    fn the_window_merges_in_ordinal_order_once() {
        let database = open("nucl_mv");
        let db = &database.db;
        let sequences = db.num_seqs();
        assert_eq!(sequences, 12);
        let want: Vec<(u32, u64)> = (0..sequences)
            .filter(|ordinal| ordinal % 3 != 2)
            .map(|ordinal| (ordinal, work(ordinal)))
            .collect();
        for mode in [XDbStream::On, XDbStream::Shadow] {
            for threads in 1..=8 {
                for block_bytes in [1u64, 40, 1 << 20] {
                    let extends = AtomicU64::new(7);
                    let scratches = AtomicU64::new(0);
                    let mut merged: Vec<(u32, u64)> = Vec::new();
                    crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                        x_search_blocks_stream(
                            mode,
                            db,
                            pool,
                            block_bytes,
                            2,
                            &|| {
                                scratches.fetch_add(1, Ordering::Relaxed);
                                Vec::<u32>::new()
                            },
                            &|seen: &mut Vec<u32>, ordinal: u32, subject: BlastnSubjectRef<'_>| {
                                assert!(subject.length > 0);
                                // a scratch sees its sequences in increasing order
                                assert!(seen.last().is_none_or(|&last| last < ordinal));
                                seen.push(ordinal);
                                extends.fetch_add(3, Ordering::Relaxed);
                                (ordinal % 3 != 2).then(|| work(ordinal))
                            },
                            &mut |ordinal, kept| merged.push((ordinal, kept)),
                            &extends,
                        )
                    })
                    .unwrap();
                    let case = format!("{mode:?} threads {threads} block_bytes {block_bytes}");
                    assert_eq!(merged, want, "{case}");
                    assert_eq!(
                        extends.load(Ordering::Relaxed),
                        7 + 3 * u64::from(sequences),
                        "{case}"
                    );
                    if mode == XDbStream::On {
                        assert!(
                            (1..=threads as u64).contains(&scratches.load(Ordering::Relaxed)),
                            "{case}: {} scratches",
                            scratches.load(Ordering::Relaxed)
                        );
                    }
                }
            }
        }
    }

    // The list of blocks is the order `next_block` gives: consecutive, every sequence once.
    #[test]
    fn the_blocks_cover_the_database_in_order() {
        let database = open("nucl_mv");
        let db = &database.db;
        for block_bytes in [1u64, 40, 100, 1 << 20] {
            let blocks = blocks_of(db, block_bytes);
            let mut next = 0u32;
            for range in &blocks {
                assert_eq!(range.first_ordinal, next, "{block_bytes}");
                assert!(range.end_ordinal() > next, "{block_bytes}");
                next = range.end_ordinal();
            }
            assert_eq!(next, db.num_seqs(), "{block_bytes}");
        }
        assert_eq!(blocks_of(db, 1).len(), 12);
    }
}

//! EXPERIMENT (LOSAT_X_BNDBSTREAM): the parallel preliminary stage of
//! a `-db` search without rounds. The pool threads claim the blocks once each in ordinal order,
//! and one of them (the owner) hands each block's results to the hit lists in that order while
//! the others go on with the next blocks (`x_bnstream::search_in_order`, the driver of the
//! `-subject` path's LOSAT_X_BNSTREAM). A thread claims a block only while it is fewer than
//! `X_BNDBSTREAM_BLOCKS_PER_THREAD` x threads blocks past the last one handed over, so that at
//! most that many blocks' results wait for the owner. A thread that waits (for room in the window,
//! or the owner for the next block) runs the pool's pending jobs first (`search_in_order`'s
//! `assist`): the chunks of a long sequence that a search hands out (a nested parallel stage) are
//! searched by the waiting threads too, as the idle workers of the rounds' `par_iter` search them.
//!
//! EXPERIMENT (LOSAT_X_BNDBWIDE, a switch of the strict set: on unless it is `0`): the window is
//! `X_BNDBWIDE_BLOCKS_PER_THREAD` x threads blocks, and it narrows to
//! `X_BNDBSTREAM_BLOCKS_PER_THREAD` x threads while more than
//! `X_BNDBWIDE_WAITING_BYTES_PER_THREAD` x threads bytes of results wait for the owner
//! (`WaitingBound`). A window of 2 blocks per thread was shorter than the blocks themselves: on
//! core_nt at 16 threads its 32 blocks were about 24 ms of the wall clock, against 54 ms for a block
//! at p99, or for the owner's own search of one. The threads waited for room 3,344 s of 10,860 s,
//! evenly over the run (`dropin-survey-20261029/BIGDB/STATE.md`, 2026-10-31). The bytes keep the
//! bound on memory: the results that wait are at most the larger of the old window's and those
//! bytes, plus the blocks in flight. They act behind a stall on a database of many hits: 1,600
//! mutated Sakai copies in 1 kb contigs with slow stretches, draft200 at 16 threads, 25.6 s and 114
//! MiB of max RSS with them, 24.7 s and 133 MiB without (`dropin-survey-20261029/BIGDB/s3/dense/`).
//!
//! Without the switch, `search_blocks_of` searches rounds of `BLOCKS_PER_THREAD` x threads blocks
//! (`search_rounds`) and merges a round's results when the whole round is searched, so the
//! preliminary HSP lists of every sequence of the round are in memory together. A database of at
//! most `BLOCKS_PER_THREAD` x threads blocks is one round: blocks hold at most 1 MiB, with or
//! without LOSAT_X_BNDBCUT, so every database below 8 x threads MiB of sequence bytes is (128 MiB
//! at 16 threads). Megablast of 100 16S queries against the 28,036 16S rRNA genes on the
//! workstation (one round from 2 threads on): at 16 threads the round held 69 MiB of preliminary
//! HSP lists at the peak of the live allocations (the lookup table is 16 MiB), and the arenas kept
//! the freed space to the end of the run (v4-B session 10, probe of the allocations by thread and
//! region, `v4-B/results/b14probe`); max RSS with outfmt 8, 33 MiB at one thread, 117-120 MiB at
//! 16 threads in rounds and 81-82 MiB with the window, against NCBI `-db`'s 108 MiB
//! (`v4-B/results/b14ab/mem.tsv`).
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
//! The window's width and the bytes decide only when a block is claimed, so they change none of
//! the above.
//!
//! Memory: at most the window's blocks of results wait for the owner, and each thread holds one
//! block and one scratch (the rounds make a scratch per `map_init` run, a few hundred per thread
//! in the case above). The blocks are listed `X_SEGMENT_BLOCKS` at a time, one `BlockRange` (24
//! bytes on 64-bit targets) each, at most 1.5 MiB however many blocks the database has (a damaged
//! index makes blocks of one sequence); the threads end a segment together before the next is
//! listed.
//!
//! No NCBI counterpart for the window itself: NCBI's search threads take chunks of ordinal ids
//! from one iterator and write each sequence's HSP list to the HSP stream when it is done
//! (blast_engine.c:1409, 1553-1554, quoted above `search_blocks`); the hand-over in ordinal order
//! keeps the order of its one-thread loop.
//!
//! `LOSAT_X_BNDBSTREAM` is a switch of the strict set: on unless it is `0`.
//! It acts for every BLASTN task with `-db` at more than
//! one thread; at one thread the blocks are searched one after another, as without it.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::OnceLock;

use anyhow::Result;

use crate::seqdb::{BlockRange, SeqDb};
use crate::utils::threading::SearchPool;

use super::super::x_bnstream::{search_in_order_bounded, WaitingBound};
use super::{block_subject, search_read_error, BlastnSubjectRef, XWaitingBytes};

/// No NCBI counterpart (scheduling only): whether LOSAT_X_BNDBSTREAM is on (a strict-set switch:
/// on unless it is `0`).
pub(super) fn x_bn_db_stream() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNDBSTREAM"))
}

/// No NCBI counterpart (scheduling only): whether LOSAT_X_BNDBWIDE is on (a strict-set switch:
/// on unless it is `0`).
fn x_bn_db_wide() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNDBWIDE"))
}

// No NCBI counterpart (scheduling only): the window; it changes no value NCBI computes.
/// Blocks per thread of the window that bytes of waiting results cannot narrow: a thread finds a
/// block to claim while the owner waits for the oldest one.
const X_BNDBSTREAM_BLOCKS_PER_THREAD: usize = 2;

// No NCBI counterpart (scheduling only): as the window.
/// Blocks per thread of the window with LOSAT_X_BNDBWIDE. Measured on core_nt (575,156 blocks,
/// 10.7 ms at the median, 54 ms at p99): 32 blocks are about 24 ms of the wall clock at 16 threads,
/// less than one slow block; 256 are about 190 ms (module documentation).
const X_BNDBWIDE_BLOCKS_PER_THREAD: usize = 16;

// No NCBI counterpart (scheduling only): as the window.
/// Bytes of waiting results per thread above which LOSAT_X_BNDBWIDE's window narrows to
/// `X_BNDBSTREAM_BLOCKS_PER_THREAD`.
const X_BNDBWIDE_WAITING_BYTES_PER_THREAD: usize = 1 << 20;

/// No NCBI counterpart (scheduling only): the heap bytes of a block's kept results.
fn kept_bytes<R: XWaitingBytes>(kept: &Vec<(u32, R)>) -> usize {
    kept.iter()
        .map(|(_, result)| result.x_waiting_bytes())
        .sum::<usize>()
        .saturating_add(kept.x_waiting_bytes())
}

// No NCBI counterpart (scheduling only): as the window.
/// Blocks listed at a time.
const X_SEGMENT_BLOCKS: usize = 1 << 16;

// No NCBI counterpart (scheduling only): `SeqDb::next_block`'s blocks, which follow each other in
// ordinal order (the rounds list them alike).
/// Up to `max` blocks of at most `block_bytes` of `db` from ordinal id `first_ordinal` on, in
/// ordinal order.
fn blocks_of(db: &SeqDb, first_ordinal: u32, block_bytes: u64, max: usize) -> Vec<BlockRange> {
    let mut blocks = Vec::new();
    let mut next_ordinal = first_ordinal;
    while blocks.len() < max {
        let Some(range) = db.next_block(next_ordinal, block_bytes) else {
            break;
        };
        blocks.push(range);
        next_ordinal = range.end_ordinal();
    }
    blocks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413,1553-1554
// ```c
//     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ...
//           /* Save the results. */
//           status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
/// The window over the whole database: its blocks of at most `block_bytes`, `segment` blocks
/// after `segment` blocks (`X_SEGMENT_BLOCKS`), through `search_window`; the first error in
/// ordinal order is returned after the blocks before it are merged.
fn stream_blocks<T, R, I, K, M>(
    db: &SeqDb,
    pool: &SearchPool<'_>,
    block_bytes: u64,
    segment: usize,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send + XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    let mut next_ordinal = 0u32;
    loop {
        let blocks = blocks_of(db, next_ordinal, block_bytes, segment.max(1));
        let Some(last) = blocks.last() else {
            return Ok(());
        };
        next_ordinal = last.end_ordinal();
        search_window(db, pool, &blocks, new_scratch, search, merge)?;
    }
}

/// The preliminary stage of `search_blocks_of` at more than one thread with LOSAT_X_BNDBSTREAM:
/// the blocks of at most `block_bytes` searched by the threads of `pool` and their results passed
/// to `merge` in ordinal order (module documentation).
pub(super) fn x_search_blocks_stream<T, R, I, K, M>(
    db: &SeqDb,
    pool: &SearchPool<'_>,
    block_bytes: u64,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send + XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    stream_blocks(
        db,
        pool,
        block_bytes,
        X_SEGMENT_BLOCKS,
        new_scratch,
        search,
        merge,
    )
}

// No NCBI counterpart (scheduling only): the window (module documentation).
/// The window over `blocks`: claimed by the threads of `pool` in order, and each block's results
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
    R: Send + XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    let stop = AtomicBool::new(false);
    let waiting_bytes = AtomicUsize::new(0);
    let search_block = |scratch: &mut T, i: usize| -> Result<Vec<(u32, R)>> {
        if stop.load(Ordering::Relaxed) {
            return Ok(Vec::new());
        }
        let range = blocks[i];
        let block = super::x_read_block(db, range).map_err(search_read_error)?;
        let mut kept = Vec::new();
        for j in 0..block.len() {
            let ordinal = range.first_ordinal + j as u32;
            if let Some(result) = search(scratch, ordinal, block_subject(&block, j)?) {
                kept.push((ordinal, result));
            }
        }
        super::x_done_block(block);
        waiting_bytes.fetch_add(kept_bytes(&kept), Ordering::Relaxed);
        Ok(kept)
    };
    let threads = pool.threads().max(1);
    let floor = threads.saturating_mul(X_BNDBSTREAM_BLOCKS_PER_THREAD);
    let wide = x_bn_db_wide();
    let window = if wide {
        threads.saturating_mul(X_BNDBWIDE_BLOCKS_PER_THREAD)
    } else {
        floor
    };
    let bound = wide.then_some(WaitingBound {
        floor,
        waiting_bytes: &waiting_bytes,
        max_bytes: threads.saturating_mul(X_BNDBWIDE_WAITING_BYTES_PER_THREAD),
    });
    if wide {
        // A list of `floor` blocks or fewer behaves as without the switch.
        crate::utils::xreach::x_reach!("LOSAT_X_BNDBWIDE", blocks.len() > floor);
    }
    let mut first_error: Option<anyhow::Error> = None;
    let mut hand_over = |_: usize, found: Result<Vec<(u32, R)>>| {
        if first_error.is_some() {
            return;
        }
        match found {
            Ok(kept) => {
                let bytes = kept_bytes(&kept);
                for (ordinal, result) in kept {
                    merge(ordinal, result);
                }
                waiting_bytes.fetch_sub(bytes, Ordering::Relaxed);
            }
            Err(error) => {
                // The blocks not claimed yet are not searched; those claimed are dropped.
                stop.store(true, Ordering::Relaxed);
                first_error = Some(error);
            }
        }
    };
    pool.install(|| {
        search_in_order_bounded(
            blocks.len(),
            window,
            true,
            bound,
            new_scratch,
            &search_block,
            &mut hand_over,
        )
    });
    match first_error {
        Some(error) => Err(error),
        None => {
            // Every block was handed over: a count that does not return to 0 would pin the window
            // at its floor, or lift the bound, without changing a report.
            debug_assert_eq!(waiting_bytes.load(Ordering::Relaxed), 0);
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::{search_rounds, BlastnDatabase};
    use super::*;
    use std::sync::atomic::AtomicU64;

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
    // The bytes that wait are the results' lists and the list of the block's kept results.
    #[test]
    fn the_kept_bytes_are_the_capacities_of_the_lists() {
        let mut kept: Vec<(u32, Vec<u64>)> = Vec::with_capacity(3);
        kept.push((0, Vec::with_capacity(100)));
        kept.push((1, Vec::with_capacity(1000)));
        let outer = kept.capacity() * std::mem::size_of::<(u32, Vec<u64>)>();
        let inner = (kept[0].1.capacity() + kept[1].1.capacity()) * std::mem::size_of::<u64>();
        assert_eq!(kept_bytes(&kept), outer + inner);
    }

    fn work(ordinal: u32) -> u64 {
        let rounds = if ordinal % 5 == 1 { 200_000 } else { 100 };
        (0..rounds).fold(u64::from(ordinal), |acc, x| {
            acc.wrapping_mul(31).wrapping_add(x)
        })
    }

    // Every sequence is searched once, the kept results reach `merge` once each in ordinal
    // order, and a thread makes one scratch at most: at 1 to 8 threads, blocks of one sequence
    // (1 byte) to the whole database.
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
        for threads in 1..=8 {
            for block_bytes in [1u64, 40, 1 << 20] {
                let searched = AtomicU64::new(0);
                let scratches = AtomicU64::new(0);
                let mut merged: Vec<(u32, u64)> = Vec::new();
                crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                    x_search_blocks_stream(
                        db,
                        pool,
                        block_bytes,
                        &|| {
                            scratches.fetch_add(1, Ordering::Relaxed);
                            Vec::<u32>::new()
                        },
                        &|seen: &mut Vec<u32>, ordinal: u32, subject: BlastnSubjectRef<'_>| {
                            assert!(subject.length > 0);
                            // a scratch sees its sequences in increasing order
                            assert!(seen.last().is_none_or(|&last| last < ordinal));
                            seen.push(ordinal);
                            searched.fetch_add(1, Ordering::Relaxed);
                            (ordinal % 3 != 2).then(|| work(ordinal))
                        },
                        &mut |ordinal, kept| merged.push((ordinal, kept)),
                    )
                })
                .unwrap();
                let case = format!("threads {threads} block_bytes {block_bytes}");
                assert_eq!(merged, want, "{case}");
                assert_eq!(
                    searched.load(Ordering::Relaxed),
                    u64::from(sequences),
                    "{case}"
                );
                assert!(
                    (1..=threads as u64).contains(&scratches.load(Ordering::Relaxed)),
                    "{case}: {} scratches",
                    scratches.load(Ordering::Relaxed)
                );
            }
        }
    }

    // The list of blocks is the order `next_block` gives: consecutive, every sequence once.
    #[test]
    fn the_blocks_cover_the_database_in_order() {
        let database = open("nucl_mv");
        let db = &database.db;
        for block_bytes in [1u64, 40, 100, 1 << 20] {
            let blocks = blocks_of(db, 0, block_bytes, usize::MAX);
            let mut next = 0u32;
            for range in &blocks {
                assert_eq!(range.first_ordinal, next, "{block_bytes}");
                assert!(range.end_ordinal() > next, "{block_bytes}");
                next = range.end_ordinal();
            }
            assert_eq!(next, db.num_seqs(), "{block_bytes}");
        }
        assert_eq!(blocks_of(db, 0, 1, usize::MAX).len(), 12);
        // Segments: up to `max` blocks from the given ordinal id on.
        let first = blocks_of(db, 0, 1, 5);
        assert_eq!(first.len(), 5);
        let rest = blocks_of(db, first[4].end_ordinal(), 1, usize::MAX);
        assert_eq!(rest.len(), 7);
        assert_eq!(rest[0].first_ordinal, 5);
    }

    // Segments of 1, 5 and 64 blocks: every kept result once, in ordinal order.
    #[test]
    fn segments_merge_in_ordinal_order_once() {
        let database = open("nucl_mv");
        let db = &database.db;
        let want: Vec<u32> = (0..db.num_seqs()).filter(|o| o % 2 == 0).collect();
        for threads in [2usize, 5] {
            for segment in [1usize, 5, 64] {
                let mut merged = Vec::new();
                let result = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                    Ok(stream_blocks(
                        db,
                        pool,
                        1,
                        segment,
                        &|| (),
                        &|_: &mut (), ordinal: u32, _: BlastnSubjectRef<'_>| {
                            (ordinal % 2 == 0).then_some(ordinal)
                        },
                        &mut |ordinal, kept| {
                            assert_eq!(ordinal, kept);
                            merged.push(ordinal);
                        },
                    ))
                })
                .unwrap();
                result.unwrap();
                assert_eq!(merged, want, "threads {threads} segment {segment}");
            }
        }
    }

    // A volume in the middle of the database cut short: the window returns the rounds' error (the
    // first failing block in ordinal order) after merging the same blocks before it, and nothing
    // after it, at 2 to 8 threads with blocks of one sequence, of 40 bytes and of 1 MiB.
    #[test]
    fn a_failing_block_in_the_middle_stops_the_window_as_it_stops_the_rounds() {
        let dir = std::env::temp_dir().join(format!(
            "losat-bndbstream-middle-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db");
        for entry in std::fs::read_dir(&fixtures).unwrap() {
            let entry = entry.unwrap();
            let name = entry.file_name().into_string().unwrap();
            if name.starts_with("nucl_mv.") {
                std::fs::copy(entry.path(), dir.join(&name)).unwrap();
            }
        }
        let nsq = dir.join("nucl_mv.02.nsq");
        let bytes = std::fs::read(&nsq).unwrap();
        std::fs::write(&nsq, &bytes[..bytes.len() / 2]).unwrap();
        let database = BlastnDatabase::open(
            dir.join("nucl_mv").to_str().unwrap(),
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .unwrap();
        let db = &database.db;
        let kernel = |_: &mut (), ordinal: u32, _: BlastnSubjectRef<'_>| Some(ordinal);
        for threads in 2..=8 {
            for block_bytes in [1u64, 40, 1 << 20] {
                let (mut rounds, mut window) = (Vec::new(), Vec::new());
                let (rounds_result, window_result) =
                    crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                        let rounds_result = search_rounds(
                            db,
                            pool,
                            block_bytes,
                            2,
                            &|| (),
                            &kernel,
                            &mut |o, k| {
                                assert_eq!(o, k);
                                rounds.push(o);
                            },
                        );
                        let window_result = stream_blocks(
                            db,
                            pool,
                            block_bytes,
                            3,
                            &|| (),
                            &kernel,
                            &mut |o, k| {
                                assert_eq!(o, k);
                                window.push(o);
                            },
                        );
                        Ok((rounds_result, window_result))
                    })
                    .unwrap();
                let case = format!("threads {threads} block_bytes {block_bytes}");
                let rounds_error = format!("{:#}", rounds_result.expect_err(&case));
                let window_error = format!("{:#}", window_result.expect_err(&case));
                assert_eq!(window_error, rounds_error, "{case}");
                assert_eq!(window, rounds, "{case}");
                assert!(window.len() < db.num_seqs() as usize, "{case}");
                assert_eq!(
                    window,
                    (0..window.len() as u32).collect::<Vec<_>>(),
                    "{case}"
                );
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}

//! EXPERIMENT (LOSAT_X_QPREPPAR): the BLASTNA encoding of a query
//! batch and the ungapped Karlin blocks of its contexts made on the search pool, as ordered
//! parallel maps. The switch also makes the batch's DUST masks on the pool
//! (`algorithm/blastn/x_query_dust_par.rs`, formerly LOSAT_X_QDUSTPAR, folded in on 2026-11-01).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2754-2784
//! NCBI encodes the queries of a batch one after the other into one buffer and computes the
//! ungapped Karlin block of each context one after the other. `search_query_batch` does the
//! same, serially, before the lookup table is built (about 0.03 s for a batch of 5 Mb of
//! queries, while no other thread works).
//!
//! With the switch on, and a LOSAT search pool of two or more threads current, both steps are
//! ordered parallel maps run with `SearchPool::install`, so on the pool of the search that calls
//! them. Rayon's global pool is never used (`SearchPool::install` runs the work on the pool the
//! search owns). A search without a pool (one thread, a build without the `parallel` feature,
//! WebAssembly without `wasm-threads`) runs today's serial code, and so does
//! a batch of fewer than two queries (encoding) or two contexts (Karlin blocks).
//!
//! - Encoding: the buffer holds the contexts of each query between its first and last context.
//!   When these spans are in query order and disjoint (always, in the layout of
//!   `search_query_batch`; if not, the serial function runs), the buffer is cut into the spans with
//!   `split_at_mut` and each query writes its own span. A query's letters are written by exactly
//!   the two calls the serial function makes (plus context, then minus context), into the same
//!   bytes; the bytes between and outside the spans stay the sentinels, as in the serial buffer.
//! - Karlin blocks: `par_iter().map().collect()` over the contexts, which returns the results in
//!   the order of the input; rayon splits the contexts into ranges, so the blocks are per context
//!   but the jobs are not.
//!
//! Why the results are the same: the encoding of a query is a function of its letters, its two
//! context offsets and the length of the block only. The functions it calls read constant tables
//! (`IUPACNA_TO_BLASTNA_FULL`) and the table of `blastna_reverse_complement_table`, a
//! `OnceLock` that every caller initialises to the same 256 bytes. The Karlin block of a context
//! is a function of the context's bytes, `reward` and `penalty` only: `query_ungapped_karlin_blastna`
//! counts the letters into a local array, builds the BLASTNA matrix and the score frequency
//! profile as locals and calls `compute_karlin_params_ungapped`, which has no global
//! state of its own. Neither uses a cache, a memo table, a counter or the random generator, and
//! no environment-variable path. Every block is
//! computed whole by one thread from its own context, in the operation order of the serial code, and the
//! blocks are collected in context order, so no floating-point operation changes order (the sums
//! are inside one context, never across contexts). The verbose lines, the totals and the empty
//! cases stay in `search_query_batch`, in the same order.
//!
//! Programs that reach it: only BLASTN (megablast, dc-megablast, blastn, blastn-short), through
//! `search_query_batch` (every query batch of the loop; the `-db` search of BLASTN uses the same
//! function). BLASTP, BLASTX, TBLASTN and TBLASTX do not call it.
//!
//! With `LOSAT_X_STATS`, `main` prints
//! `[X_STATS] QPREPPAR encode_s= karlin_s= fallbacks=`, the summed wall times of the two steps over
//! every batch and the number of batches with two or more queries on a pool whose spans were not
//! in query order and disjoint, which the serial function encoded (a fallback; the layout of
//! `search_query_batch` never causes one), with or without the switch.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::algorithm::blastn::coordination::QueryLetters;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastn::lookup::encode_reverse_complement_blastna_into;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastn::pairwise::query_ungapped_karlin_blastna;
use crate::algorithm::blastn::scoring::context_ungapped_blocks;
use crate::config::NuclScoringSpec;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::core::blast_encoding::encode_iupac_to_blastna_into;
use crate::stats::KarlinParams;
use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switch once (LOSAT_X_QPREPPAR is a strict-set switch since
// 2026-10-17: on unless it is `0`). It chooses between two LOSAT code paths with the same result
// and changes no value NCBI computes.
/// Whether LOSAT_X_QPREPPAR is on.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_QPREPPAR"))
}

// No NCBI counterpart: the fallbacks and the summed wall times of the two steps (LOSAT_X_STATS);
// they change no value NCBI computes.
static FALLBACKS: AtomicU64 = AtomicU64::new(0);
static ENCODE_NANOS: AtomicU64 = AtomicU64::new(0);
static KARLIN_NANOS: AtomicU64 = AtomicU64::new(0);
static CALLS: AtomicU64 = AtomicU64::new(0);

/// Whether the two steps are timed (the `[X_STATS]` line costs four clock reads per batch).
fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_STATS").is_some())
}

/// Whether the work of `count` items (queries or contexts) runs on `pool`: two or more items and
/// a LOSAT search pool of two or more threads.
fn on_pool(pool: &SearchPool<'_>, count: usize) -> bool {
    count > 1 && pool.enabled() && pool.threads() > 1
}

/// The serial buffer: today's `fill_query_blastna_buffer`, into `buffer` (the batch loop's).
fn serial_fill<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    buffer: &mut Vec<u8>,
) {
    super::fill_query_blastna_buffer(query_letters, context_offsets, query_concat_length, buffer)
}

/// `serial_fill` into a buffer of its own (the tests).
#[cfg(all(test, feature = "parallel"))]
fn serial_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
) -> Vec<u8> {
    let mut buffer = Vec::new();
    serial_fill(
        query_letters,
        context_offsets,
        query_concat_length,
        &mut buffer,
    );
    buffer
}

/// The serial blocks: today's `context_ungapped_blocks` over `contexts`.
fn serial_blocks(contexts: &[&[u8]], spec: &NuclScoringSpec) -> Vec<Option<KarlinParams>> {
    context_ungapped_blocks(contexts.iter().copied(), spec)
}

/// The range of a context of `length` letters at `offset` in the buffer with its first sentinel
/// (`None` for a context that is empty or does not end inside the block, which stays sentinels).
// Copy of the `place` closure in `fill_query_blastna_buffer` (run.rs): the serial function stays
// untouched with the switch off, so a change of its rule must be made here too; the random
// pool-vs-serial tests compare the two.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn place(offset: i32, length: usize, query_concat_length: usize) -> Option<std::ops::Range<usize>> {
    let start = offset.max(0) as usize;
    let end = start.saturating_add(length);
    (start < end && end <= query_concat_length).then(|| start + 1..end + 1)
}

/// The parallel buffer, into `buffer`; the same bytes as `serial_fill`. Only called when
/// `on_pool`. The flag is `true` when the queries were encoded on the pool and `false` when the
/// spans were not in query order and disjoint and the serial function made the buffer (a
/// fallback). `probe` is called once per query, on the thread that encodes it; production passes
/// `|| {}`, which the compiler removes (a unit test passes a closure that looks at the thread).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel_fill<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    probe: &P,
    buffer: &mut Vec<u8>,
) -> bool {
    debug_assert_eq!(2 * query_letters.count(), context_offsets.len());
    // (first byte, end, query) of the span of each query with a placed context, in query order.
    let mut spans: Vec<(usize, usize, usize)> = Vec::with_capacity(query_letters.count());
    for (q, offsets) in context_offsets.chunks_exact(2).enumerate() {
        let length = query_letters.length(q);
        let r0 = place(offsets[0], length, query_concat_length);
        let r1 = place(offsets[1], length, query_concat_length);
        let (a, b) = match (&r0, &r1) {
            (Some(x), Some(y)) => (x.start.min(y.start), x.end.max(y.end)),
            (Some(x), None) | (None, Some(x)) => (x.start, x.end),
            (None, None) => continue,
        };
        spans.push((a, b, q));
    }
    // The spans of different queries must not overlap and must follow the query order (the serial
    // function writes in query order, so overlapping contexts would make its bytes depend on it).
    if !spans.windows(2).all(|w| w[0].1 <= w[1].0) {
        serial_fill(query_letters, context_offsets, query_concat_length, buffer);
        return false;
    }
    super::x_query_sentinels(buffer, query_concat_length + 2);
    let mut parts: Vec<(&mut [u8], usize, usize)> = Vec::with_capacity(spans.len());
    let mut rest: &mut [u8] = buffer;
    let mut base = 0usize;
    for &(a, b, q) in &spans {
        let (_, tail) = std::mem::take(&mut rest).split_at_mut(a - base);
        let (part, tail) = tail.split_at_mut(b - a);
        parts.push((part, a, q));
        rest = tail;
        base = b;
    }
    pool.install(|| {
        parts.par_iter_mut().for_each(|(part, a, q)| {
            probe();
            let offsets = &context_offsets[2 * *q..2 * *q + 2];
            // The same two calls as `fill_query_blastna_buffer` (run.rs), in its order, into the
            // same bytes (`range` is relative to the buffer, `a` is the buffer index of the part's
            // first byte). A copy: a change of the serial calls must be made here too.
            query_letters.read(*q, |letters| {
                if let Some(range) = place(offsets[0], letters.len(), query_concat_length) {
                    encode_iupac_to_blastna_into(
                        letters,
                        &mut part[range.start - *a..range.end - *a],
                    );
                }
                if let Some(range) = place(offsets[1], letters.len(), query_concat_length) {
                    encode_reverse_complement_blastna_into(
                        letters,
                        &mut part[range.start - *a..range.end - *a],
                    );
                }
            });
        });
    });
    true
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel_fill<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    _pool: &SearchPool<'_>,
    _probe: &P,
    buffer: &mut Vec<u8>,
) -> bool {
    serial_fill(query_letters, context_offsets, query_concat_length, buffer);
    false
}

/// `parallel_fill` into a buffer of its own (the tests).
#[cfg(all(test, feature = "parallel"))]
fn parallel_buffer<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    probe: &P,
) -> (Vec<u8>, bool) {
    let mut buffer = Vec::new();
    let on_the_pool = parallel_fill(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        probe,
        &mut buffer,
    );
    (buffer, on_the_pool)
}

/// The parallel blocks; the same blocks as `serial_blocks`. Only called when `on_pool`. `probe` is
/// called once per context, on the thread that computes it (see `parallel_fill`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel_blocks<P: Fn() + Sync>(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
    probe: &P,
) -> Vec<Option<KarlinParams>> {
    let (reward, penalty) = (spec.reward, spec.penalty);
    pool.install(|| {
        contexts
            .par_iter()
            .map(|context| {
                probe();
                // A copy of the body of `context_ungapped_blocks` (scoring.rs): a change of the
                // serial function must be made here too.
                if context.is_empty() {
                    None
                } else {
                    query_ungapped_karlin_blastna(context, reward, penalty)
                }
            })
            .collect()
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel_blocks<P: Fn() + Sync>(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    _pool: &SearchPool<'_>,
    _probe: &P,
) -> Vec<Option<KarlinParams>> {
    serial_blocks(contexts, spec)
}

/// Whether two blocks are the same bit for bit (`KarlinParams` has no `PartialEq`, and `==` on
/// floats would equate `0.0` and `-0.0` and differ on a NaN).
#[cfg(test)]
fn same_bits(a: &Option<KarlinParams>, b: &Option<KarlinParams>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.lambda.to_bits() == b.lambda.to_bits()
                && a.k.to_bits() == b.k.to_bits()
                && a.h.to_bits() == b.h.to_bits()
                && a.alpha.to_bits() == b.alpha.to_bits()
                && a.beta.to_bits() == b.beta.to_bits()
        }
        _ => false,
    }
}

/// The query buffer as the switch `on` makes it, into `buffer`. Separate from the environment so
/// that tests choose the switch.
pub(super) fn fill_with_mode<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    on: bool,
    buffer: &mut Vec<u8>,
) {
    fill_with_probe(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        on,
        &|| {},
        buffer,
    )
}

/// `fill_with_mode` with a `probe` that the parallel map calls once per query.
fn fill_with_probe<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    on: bool,
    probe: &P,
    buffer: &mut Vec<u8>,
) {
    if on {
        crate::utils::xreach::x_reach!("LOSAT_X_QPREPPAR", on_pool(pool, query_letters.count()));
    }
    if !on || !on_pool(pool, query_letters.count()) {
        serial_fill(query_letters, context_offsets, query_concat_length, buffer);
        return;
    }
    let on_the_pool = parallel_fill(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        probe,
        buffer,
    );
    if !on_the_pool {
        // The serial function made it (spans not in query order and disjoint): not an encoding on
        // a pool, nothing to compare.
        FALLBACKS.fetch_add(1, Ordering::Relaxed);
    }
}

/// `fill_with_mode` into a buffer of its own (the tests).
#[cfg(all(test, feature = "parallel"))]
fn buffer_with_mode<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    on: bool,
) -> Vec<u8> {
    let mut buffer = Vec::new();
    fill_with_mode(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        on,
        &mut buffer,
    );
    buffer
}

/// `fill_with_probe` into a buffer of its own (the tests).
#[cfg(all(test, feature = "parallel"))]
fn buffer_with_probe<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    on: bool,
    probe: &P,
) -> Vec<u8> {
    let mut buffer = Vec::new();
    fill_with_probe(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        on,
        probe,
        &mut buffer,
    );
    buffer
}

/// The ungapped blocks of `contexts` as the switch `on` makes them (see `fill_with_mode`).
pub(super) fn blocks_with_mode(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
    on: bool,
) -> Vec<Option<KarlinParams>> {
    blocks_with_probe(contexts, spec, pool, on, &|| {})
}

/// `blocks_with_mode` with a `probe` that the parallel map calls once per context.
fn blocks_with_probe<P: Fn() + Sync>(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
    on: bool,
    probe: &P,
) -> Vec<Option<KarlinParams>> {
    if on {
        crate::utils::xreach::x_reach!("LOSAT_X_QPREPPAR", on_pool(pool, contexts.len()));
    }
    if !on || !on_pool(pool, contexts.len()) {
        return serial_blocks(contexts, spec);
    }
    parallel_blocks(contexts, spec, pool, probe)
}

/// The BLASTNA buffer of a query batch with its sentinels, as `fill_query_blastna_buffer` makes it,
/// into `buffer` (the batch loop's): the environment's switch decides.
pub(super) fn query_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    buffer: &mut Vec<u8>,
) {
    let start = timing_on().then(std::time::Instant::now);
    fill_with_mode(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        on(),
        buffer,
    );
    if let Some(start) = start {
        ENCODE_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        CALLS.fetch_add(1, Ordering::Relaxed);
    }
}

/// The ungapped Karlin block of each context, as `context_ungapped_blocks` makes them: the
/// environment's switch decides.
pub(super) fn ungapped_blocks(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
) -> Vec<Option<KarlinParams>> {
    let start = timing_on().then(std::time::Instant::now);
    let blocks = blocks_with_mode(contexts, spec, pool, on());
    if let Some(start) = start {
        KARLIN_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
    blocks
}

/// With LOSAT_X_STATS the summed wall times of the two steps (printed by `main` at exit), only when
/// at least one batch was prepared.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if timing_on() && CALLS.load(Ordering::Relaxed) > 0 {
        eprintln!(
            "[X_STATS] QPREPPAR encode_s={:.3} karlin_s={:.3} fallbacks={}",
            ENCODE_NANOS.load(Ordering::Relaxed) as f64 / 1e9,
            KARLIN_NANOS.load(Ordering::Relaxed) as f64 / 1e9,
            FALLBACKS.load(Ordering::Relaxed)
        );
    }
}

#[cfg(test)]
#[path = "x_query_prep_par_tests.rs"]
mod tests;

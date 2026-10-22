//! EXPERIMENT (LOSAT_X_QPREPPAR, shadow LOSAT_X_QPREPPARSHADOW): the BLASTNA encoding of a query
//! batch and the ungapped Karlin blocks of its contexts made on the search pool, as ordered
//! parallel maps.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
//! ```c++
//! int buflen = QueryInfo_GetSeqBufLen(qinfo);
//! TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
//! ...
//! sequence = queries.GetBlastSequence(index, encoding, strand, eSentinels);
//! ...
//! memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2754-2784
//! ```c
//!     for (context = query_info->first_context;
//!          context <= query_info->last_context; context++) {
//! ...
//!         Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
//!         sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
//!         BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
//!         sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
//!         loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
//! ```
//! NCBI encodes the queries of a batch one after the other into one buffer and computes the
//! ungapped Karlin block of each context one after the other. `search_query_batch` does the
//! same, serially, before the lookup table is built (about 0.03 s for a batch of 5 Mb of
//! queries, while no other thread works).
//!
//! With the switch on, and a LOSAT search pool of two or more threads current, both steps are
//! ordered parallel maps run with `SearchPool::install`, so on the pool of the search that calls
//! them: the whole search pool, or the private pool of a LOSAT_X_QBATCHPAR slot with two or more
//! threads. Rayon's global pool is never used (`SearchPool::install` runs the work on the pool the
//! search owns). A search without a pool (one thread, a slot of one thread, a build without the
//! `parallel` feature, WebAssembly without `wasm-threads`) runs today's serial code, and so does
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
//! profile as locals and calls `compute_karlin_params_ungapped`, whose only global is the
//! `LOSAT_X_KARLINFAST` flag (`OnceLock`, read once, never written again). Neither uses a cache,
//! a memo table, a counter or the random generator, and no `LOSAT_X_ENVCACHE` path. Every block is
//! computed whole by one thread from its own context, in the operation order of the serial code, and the
//! blocks are collected in context order, so no floating-point operation changes order (the sums
//! are inside one context, never across contexts). The verbose lines, the totals and the empty
//! cases stay in `search_query_batch`, in the same order.
//!
//! Programs that reach it: only BLASTN (megablast, dc-megablast, blastn, blastn-short), through
//! `search_query_batch` (every query batch, in the loop or in a LOSAT_X_QBATCHPAR slot; the `-db`
//! search of BLASTN uses the same function). BLASTP, BLASTX, TBLASTN and TBLASTX do not call it.
//!
//! The shadow also encodes and computes the blocks serially and compares them: the encodings
//! byte for byte, the blocks field by field with the bits of every float. On a difference it
//! prints `[X_SHADOW] QPREPPAR ... DIFFERENT` and panics with a message that names what differs
//! (the first differing byte of the buffer, or the context and the field of the block); otherwise
//! `main` prints
//! `[X_SHADOW] QPREPPAR batches= queries= contexts= fallbacks= diffs=0 identical` at exit
//! (`batches=` and `queries=` count the encodings made on a pool, `contexts=` the contexts whose
//! blocks were made on a pool; `fallbacks=` counts the batches with two or more queries on a pool
//! whose spans were not in query order and disjoint, which the serial function encoded: they are
//! not counted in `batches=`/`queries=` and not compared, the layout of `search_query_batch` never
//! causes one; `diffs=` counts the differences found, also in a batch that a LOSAT_X_QBATCHPAR slot
//! made and the loop discarded; the line ends in `identical` only when it is 0, else in
//! `DIFFERENT`). With `LOSAT_X_STATS`, `main` prints
//! `[X_STATS] QPREPPAR encode_s= karlin_s= fallbacks=`, the summed wall times of the two steps over
//! every batch (so with LOSAT_X_QBATCHPAR they include the slots' batches, which run concurrently,
//! and those the loop discards; with the shadow they include the serial reference pass) and the
//! fallbacks (counted whenever the switch is on), with or without the switch.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::algorithm::blastn::lookup::encode_reverse_complement_blastna_into;
use crate::algorithm::blastn::pairwise::query_ungapped_karlin_blastna;
use crate::algorithm::blastn::scoring::context_ungapped_blocks;
use crate::config::NuclScoringSpec;
use crate::core::blast_encoding::encode_iupac_to_blastna_into;
use crate::stats::KarlinParams;
use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switches once (LOSAT_X_QPREPPAR is a strict-set switch since
// 2026-10-17: on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set; the shadow, off
// unless set to any value, implies the switch). They choose between two LOSAT code paths with
// the same result and change no value NCBI computes.
/// 0 = off, 1 = LOSAT_X_QPREPPAR, 2 = LOSAT_X_QPREPPARSHADOW.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_QPREPPARSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_QPREPPAR") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: counters of the shadow's summary line and the summed wall times of the two
// steps (LOSAT_X_STATS); they change no value NCBI computes.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static QUERIES: AtomicU64 = AtomicU64::new(0);
static CONTEXTS: AtomicU64 = AtomicU64::new(0);
static FALLBACKS: AtomicU64 = AtomicU64::new(0);
static DIFFS: AtomicU64 = AtomicU64::new(0);
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

/// The serial buffer: today's `build_query_blastna_buffer`.
fn serial_buffer(
    query_letters: &[&[u8]],
    context_offsets: &[i32],
    query_concat_length: usize,
) -> Vec<u8> {
    super::build_query_blastna_buffer(query_letters, context_offsets, query_concat_length)
}

/// The serial blocks: today's `context_ungapped_blocks` over `contexts`.
fn serial_blocks(contexts: &[&[u8]], spec: &NuclScoringSpec) -> Vec<Option<KarlinParams>> {
    context_ungapped_blocks(contexts.iter().copied(), spec)
}

/// The range of a context of `length` letters at `offset` in the buffer with its first sentinel
/// (`None` for a context that is empty or does not end inside the block, which stays sentinels).
// Copy of the `place` closure in `build_query_blastna_buffer` (run.rs): the serial function stays
// untouched with the switch off, so a change of its rule must be made here too; the random
// pool-vs-serial tests and the shadow compare the two.
fn place(offset: i32, length: usize, query_concat_length: usize) -> Option<std::ops::Range<usize>> {
    let start = offset.max(0) as usize;
    let end = start.saturating_add(length);
    (start < end && end <= query_concat_length).then(|| start + 1..end + 1)
}

/// The parallel buffer; the same bytes as `serial_buffer`. Only called when `on_pool`. The flag is
/// `true` when the queries were encoded on the pool and `false` when the spans were not in query
/// order and disjoint and the serial function made the buffer (a fallback). `probe` is called once
/// per query, on the thread that encodes it; production passes `|| {}`, which the compiler removes
/// (a unit test passes a closure that looks at the thread).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel_buffer<P: Fn() + Sync>(
    query_letters: &[&[u8]],
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    probe: &P,
) -> (Vec<u8>, bool) {
    debug_assert_eq!(2 * query_letters.len(), context_offsets.len());
    // (first byte, end, query) of the span of each query with a placed context, in query order.
    let mut spans: Vec<(usize, usize, usize)> = Vec::with_capacity(query_letters.len());
    for (q, (letters, offsets)) in query_letters
        .iter()
        .zip(context_offsets.chunks_exact(2))
        .enumerate()
    {
        let r0 = place(offsets[0], letters.len(), query_concat_length);
        let r1 = place(offsets[1], letters.len(), query_concat_length);
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
        return (
            serial_buffer(query_letters, context_offsets, query_concat_length),
            false,
        );
    }
    let mut with_sentinels = vec![super::NUCL_SENTINEL; query_concat_length + 2];
    let mut parts: Vec<(&mut [u8], usize, usize)> = Vec::with_capacity(spans.len());
    let mut rest: &mut [u8] = &mut with_sentinels;
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
            let letters = query_letters[*q];
            let offsets = &context_offsets[2 * *q..2 * *q + 2];
            // The same two calls as `build_query_blastna_buffer` (run.rs), in its order, into the
            // same bytes (`range` is relative to the buffer, `a` is the buffer index of the part's
            // first byte). A copy: a change of the serial calls must be made here too.
            if let Some(range) = place(offsets[0], letters.len(), query_concat_length) {
                encode_iupac_to_blastna_into(letters, &mut part[range.start - *a..range.end - *a]);
            }
            if let Some(range) = place(offsets[1], letters.len(), query_concat_length) {
                encode_reverse_complement_blastna_into(
                    letters,
                    &mut part[range.start - *a..range.end - *a],
                );
            }
        });
    });
    (with_sentinels, true)
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel_buffer<P: Fn() + Sync>(
    query_letters: &[&[u8]],
    context_offsets: &[i32],
    query_concat_length: usize,
    _pool: &SearchPool<'_>,
    _probe: &P,
) -> (Vec<u8>, bool) {
    (
        serial_buffer(query_letters, context_offsets, query_concat_length),
        false,
    )
}

/// The parallel blocks; the same blocks as `serial_blocks`. Only called when `on_pool`. `probe` is
/// called once per context, on the thread that computes it (see `parallel_buffer`).
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

/// What differs between the buffer made on the pool and the serial one, or `None` when they are
/// equal (length and every byte).
fn first_buffer_diff(fast: &[u8], slow: &[u8]) -> Option<String> {
    if fast.len() != slow.len() {
        return Some(format!(
            "query buffer: lengths {} / {}",
            fast.len(),
            slow.len()
        ));
    }
    fast.iter()
        .zip(slow)
        .position(|(fast, slow)| fast != slow)
        .map(|i| {
            format!(
                "query buffer: first differing byte at index {i} ({} / {})",
                fast[i], slow[i]
            )
        })
}

/// The name of the first field in which two blocks differ bit for bit ("none" when `same_bits`).
fn block_field_diff(fast: &Option<KarlinParams>, slow: &Option<KarlinParams>) -> &'static str {
    match (fast, slow) {
        (Some(a), Some(b)) => [
            ("lambda", a.lambda, b.lambda),
            ("k", a.k, b.k),
            ("h", a.h, b.h),
            ("alpha", a.alpha, b.alpha),
            ("beta", a.beta, b.beta),
        ]
        .into_iter()
        .find(|(_, x, y)| x.to_bits() != y.to_bits())
        .map_or("none", |(name, _, _)| name),
        (Some(_), None) => "Some / None",
        (None, Some(_)) => "None / Some",
        (None, None) => "none",
    }
}

/// What differs between the blocks made on the pool and the serial ones, or `None` when they are
/// equal (count and every block, field by field with the bits of the floats).
fn first_block_diff(
    fast: &[Option<KarlinParams>],
    slow: &[Option<KarlinParams>],
) -> Option<String> {
    if fast.len() != slow.len() {
        return Some(format!(
            "Karlin blocks: counts {} / {}",
            fast.len(),
            slow.len()
        ));
    }
    fast.iter()
        .zip(slow)
        .position(|(fast, slow)| !same_bits(fast, slow))
        .map(|i| {
            format!(
                "Karlin blocks: context {i}, field {} ({:?} / {:?})",
                block_field_diff(&fast[i], &slow[i]),
                fast[i],
                slow[i]
            )
        })
}

// Test-only hooks for the shadow's tests: the next buffer / block set made on the pool on this
// thread is perturbed once, after it is made and before the shadow compares it. Not compiled
// outside tests.
#[cfg(test)]
thread_local! {
    /// The index of the byte of the next pool buffer to change.
    static TAMPER_BUFFER: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
    /// The field (0 = lambda .. 4 = beta) of the first block of the next pool block set to change.
    static TAMPER_FIELD: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
}

#[cfg(test)]
fn tamper_buffer(buffer: &mut [u8]) {
    if let Some(index) = TAMPER_BUFFER.with(|cell| cell.take()) {
        buffer[index] ^= 1;
    }
}

#[cfg(not(test))]
#[inline(always)]
fn tamper_buffer(_buffer: &mut [u8]) {}

#[cfg(test)]
fn tamper_blocks(blocks: &mut [Option<KarlinParams>]) {
    if let Some(field) = TAMPER_FIELD.with(|cell| cell.take()) {
        let block = blocks
            .iter_mut()
            .flatten()
            .next()
            .expect("a block to change");
        let slot = match field {
            0 => &mut block.lambda,
            1 => &mut block.k,
            2 => &mut block.h,
            3 => &mut block.alpha,
            _ => &mut block.beta,
        };
        *slot = f64::from_bits(slot.to_bits() ^ 1);
    }
}

#[cfg(not(test))]
#[inline(always)]
fn tamper_blocks(_blocks: &mut [Option<KarlinParams>]) {}

/// The query buffer as the switch `mode` makes it (0 = off, 1 = on, 2 = on and compared with the
/// serial buffer). Separate from the environment so that tests choose the mode.
pub(super) fn buffer_with_mode(
    query_letters: &[&[u8]],
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    mode: u8,
) -> Vec<u8> {
    buffer_with_probe(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        mode,
        &|| {},
    )
}

/// `buffer_with_mode` with a `probe` that the parallel map calls once per query.
fn buffer_with_probe<P: Fn() + Sync>(
    query_letters: &[&[u8]],
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
    mode: u8,
    probe: &P,
) -> Vec<u8> {
    if mode != 0 {
        crate::utils::xreach::x_reach!("LOSAT_X_QPREPPAR", on_pool(pool, query_letters.len()));
    }
    if mode == 0 || !on_pool(pool, query_letters.len()) {
        return serial_buffer(query_letters, context_offsets, query_concat_length);
    }
    let (mut buffer, on_the_pool) = parallel_buffer(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        probe,
    );
    if !on_the_pool {
        // The serial function made it (spans not in query order and disjoint): not an encoding on
        // a pool, nothing to compare.
        FALLBACKS.fetch_add(1, Ordering::Relaxed);
        return buffer;
    }
    tamper_buffer(&mut buffer);
    if mode == 2 {
        let reference = serial_buffer(query_letters, context_offsets, query_concat_length);
        BATCHES.fetch_add(1, Ordering::Relaxed);
        QUERIES.fetch_add(query_letters.len() as u64, Ordering::Relaxed);
        if let Some(what) = first_buffer_diff(&buffer, &reference) {
            DIFFS.fetch_add(1, Ordering::Relaxed);
            eprintln!(
                "[X_SHADOW] QPREPPAR batches={} queries={} contexts={} fallbacks={} diffs={} DIFFERENT ({what})",
                BATCHES.load(Ordering::Relaxed),
                QUERIES.load(Ordering::Relaxed),
                CONTEXTS.load(Ordering::Relaxed),
                FALLBACKS.load(Ordering::Relaxed),
                DIFFS.load(Ordering::Relaxed),
            );
            panic!("LOSAT_X_QPREPPARSHADOW: the query buffer made on the pool differs from the serial one ({what})");
        }
    }
    buffer
}

/// The ungapped blocks of `contexts` as the switch `mode` makes them (see `buffer_with_mode`).
pub(super) fn blocks_with_mode(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
    mode: u8,
) -> Vec<Option<KarlinParams>> {
    blocks_with_probe(contexts, spec, pool, mode, &|| {})
}

/// `blocks_with_mode` with a `probe` that the parallel map calls once per context.
fn blocks_with_probe<P: Fn() + Sync>(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
    mode: u8,
    probe: &P,
) -> Vec<Option<KarlinParams>> {
    if mode != 0 {
        crate::utils::xreach::x_reach!("LOSAT_X_QPREPPAR", on_pool(pool, contexts.len()));
    }
    if mode == 0 || !on_pool(pool, contexts.len()) {
        return serial_blocks(contexts, spec);
    }
    let mut blocks = parallel_blocks(contexts, spec, pool, probe);
    tamper_blocks(&mut blocks);
    if mode == 2 {
        let reference = serial_blocks(contexts, spec);
        CONTEXTS.fetch_add(contexts.len() as u64, Ordering::Relaxed);
        if let Some(what) = first_block_diff(&blocks, &reference) {
            DIFFS.fetch_add(1, Ordering::Relaxed);
            eprintln!(
                "[X_SHADOW] QPREPPAR batches={} queries={} contexts={} fallbacks={} diffs={} DIFFERENT ({what})",
                BATCHES.load(Ordering::Relaxed),
                QUERIES.load(Ordering::Relaxed),
                CONTEXTS.load(Ordering::Relaxed),
                FALLBACKS.load(Ordering::Relaxed),
                DIFFS.load(Ordering::Relaxed),
            );
            panic!("LOSAT_X_QPREPPARSHADOW: the Karlin blocks made on the pool differ from the serial ones ({what})");
        }
    }
    blocks
}

/// The BLASTNA buffer of a query batch with its sentinels, as `build_query_blastna_buffer` makes it:
/// the environment's switch decides.
pub(super) fn query_buffer(
    query_letters: &[&[u8]],
    context_offsets: &[i32],
    query_concat_length: usize,
    pool: &SearchPool<'_>,
) -> Vec<u8> {
    let start = timing_on().then(std::time::Instant::now);
    let buffer = buffer_with_mode(
        query_letters,
        context_offsets,
        query_concat_length,
        pool,
        mode(),
    );
    if let Some(start) = start {
        ENCODE_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        CALLS.fetch_add(1, Ordering::Relaxed);
    }
    buffer
}

/// The ungapped Karlin block of each context, as `context_ungapped_blocks` makes them: the
/// environment's switch decides.
pub(super) fn ungapped_blocks(
    contexts: &[&[u8]],
    spec: &NuclScoringSpec,
    pool: &SearchPool<'_>,
) -> Vec<Option<KarlinParams>> {
    let start = timing_on().then(std::time::Instant::now);
    let blocks = blocks_with_mode(contexts, spec, pool, mode());
    if let Some(start) = start {
        KARLIN_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
    blocks
}

/// The shadow's summary line (printed by `main` at exit), only with LOSAT_X_QPREPPARSHADOW, and
/// with LOSAT_X_STATS the summed wall times of the two steps, only when at least one batch was
/// prepared.
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
    if mode() == 2 {
        let diffs = DIFFS.load(Ordering::Relaxed);
        eprintln!(
            "[X_SHADOW] QPREPPAR batches={} queries={} contexts={} fallbacks={} diffs={diffs} {}",
            BATCHES.load(Ordering::Relaxed),
            QUERIES.load(Ordering::Relaxed),
            CONTEXTS.load(Ordering::Relaxed),
            FALLBACKS.load(Ordering::Relaxed),
            if diffs == 0 { "identical" } else { "DIFFERENT" }
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::threading::with_search_pool;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// A random query with N, IUPAC letters, lower case, low-complexity runs, and now and then
    /// no letters at all.
    fn sequence(rng: &mut Rng) -> Vec<u8> {
        let len = match rng.below(6) {
            0 => rng.below(3) as usize,
            1 => 3 + rng.below(40) as usize,
            _ => 100 + rng.below(1500) as usize,
        };
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            match rng.below(8) {
                0 => {
                    let unit = [b"A".as_slice(), b"AT", b"CAG", b"GGGA"][rng.below(4) as usize];
                    for _ in 0..(3 + rng.below(30)) {
                        seq.extend_from_slice(unit);
                    }
                }
                1 => seq.extend(std::iter::repeat_n(b'N', 1 + rng.below(12) as usize)),
                2 => seq.push(b"RYKMSWBDHVnacgtrykmswbdhv-"[rng.below(26) as usize]),
                _ => seq.push(b"ACGT"[rng.below(4) as usize]),
            }
        }
        seq.truncate(len);
        seq
    }

    /// The context offsets and the length of the block that `search_query_batch` gives to
    /// `queries`: `context_room` for each context (its letters and one sentinel, nothing for an
    /// empty query), plus context then minus context of each query, and the
    /// `query_concat_length` formula of run.rs (the last context's offset and length plus 2, or
    /// plus 1 when it is empty, less the two end sentinels).
    fn production_layout(queries: &[Vec<u8>]) -> (Vec<i32>, usize) {
        let room = |length: usize| if length == 0 { 0 } else { length + 1 };
        let mut offsets = Vec::new();
        let mut at = 0usize;
        let mut last = (0usize, 0usize);
        for query in queries {
            let plus = at;
            let minus = plus + room(query.len());
            offsets.push(plus as i32);
            offsets.push(minus as i32);
            at = minus + room(query.len());
            last = (minus, query.len());
        }
        let total = if queries.is_empty() {
            0
        } else {
            (last.0 + last.1 + if last.1 == 0 { 1 } else { 2 }).saturating_sub(2)
        };
        (offsets, total)
    }

    /// A batch of `count` random queries in the production layout (`production_layout`): the
    /// letters, the context offsets and the length of the block.
    fn batch(rng: &mut Rng, count: usize) -> (Vec<Vec<u8>>, Vec<i32>, usize) {
        let queries: Vec<Vec<u8>> = (0..count).map(|_| sequence(rng)).collect();
        let (offsets, total) = production_layout(&queries);
        (queries, offsets, total)
    }

    /// A batch in a COMPACT layout, which is NOT the production one: the minus context follows
    /// the plus context with no sentinel between them, and the next query follows at once (`gap`
    /// is 0) or after `gap` sentinels. It tests spans that touch each other, which the production
    /// layout (a sentinel after every context) never has.
    fn compact_batch(rng: &mut Rng, count: usize, gap: usize) -> (Vec<Vec<u8>>, Vec<i32>, usize) {
        let queries: Vec<Vec<u8>> = (0..count).map(|_| sequence(rng)).collect();
        let mut offsets = Vec::new();
        let mut at = 0usize;
        for query in &queries {
            offsets.push(at as i32);
            offsets.push((at + query.len()) as i32);
            at += 2 * query.len() + gap;
        }
        (queries, offsets, at)
    }

    /// The contexts (slices of the block without its end sentinels) of a batch, in context order
    /// (an empty query gives two empty slices).
    fn contexts_of<'a>(block: &'a [u8], queries: &[Vec<u8>], offsets: &[i32]) -> Vec<&'a [u8]> {
        queries
            .iter()
            .zip(offsets.chunks_exact(2))
            .flat_map(|(query, pair)| {
                pair.iter().map(move |&offset| {
                    if query.is_empty() {
                        &[][..]
                    } else {
                        &block[offset as usize..offset as usize + query.len()]
                    }
                })
            })
            .collect()
    }

    /// The text of the panic that `work` raises (the test fails when it does not panic).
    fn panic_text<R>(work: impl FnOnce() -> R) -> String {
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)) {
            Ok(_) => panic!("expected a panic"),
            Err(payload) => payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|text| text.to_string()))
                .unwrap_or_default(),
        }
    }

    fn spec() -> NuclScoringSpec {
        NuclScoringSpec {
            reward: 1,
            penalty: -2,
            gap_open: 0,
            gap_extend: 0,
        }
    }

    /// `same_bits` compares every field by its bits: a block that differs in one field, in the sign
    /// of a zero, or between `None` and `Some` is different, and a NaN equals itself.
    #[test]
    fn same_bits_sees_every_field() {
        let base = KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 1.5,
            beta: -2.0,
        };
        assert!(same_bits(&None, &None));
        assert!(!same_bits(&None, &Some(base)));
        assert!(same_bits(&Some(base), &Some(base)));
        let variants = [
            KarlinParams {
                lambda: 0.626,
                ..base
            },
            KarlinParams { k: 0.42, ..base },
            KarlinParams { h: 0.79, ..base },
            KarlinParams { alpha: 1.6, ..base },
            KarlinParams { beta: -2.1, ..base },
        ];
        for variant in variants {
            assert!(!same_bits(&Some(base), &Some(variant)), "{variant:?}");
        }
        let zero = KarlinParams { beta: 0.0, ..base };
        let negative_zero = KarlinParams { beta: -0.0, ..base };
        assert!(!same_bits(&Some(zero), &Some(negative_zero)));
        let nan = KarlinParams {
            k: f64::NAN,
            ..base
        };
        assert!(same_bits(&Some(nan), &Some(nan)));
    }

    /// Pool == serial for the buffer at 1, 2, 3 and 7 threads on random batches (N, IUPAC, lower
    /// case, empty queries, one query, none) in the production layout, which never takes the
    /// fallback, and in the compact layouts (spans that touch, or with sentinel gaps).
    #[cfg(feature = "parallel")]
    #[test]
    fn pool_buffer_equals_serial_buffer() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        for threads in [1usize, 2, 3, 7] {
            with_search_pool(threads, "blastn", |pool| {
                for (count, gap) in [
                    (0, None),
                    (1, None),
                    (2, None),
                    (40, None),
                    (90, None),
                    (0, Some(0)),
                    (1, Some(0)),
                    (2, Some(0)),
                    (40, Some(0)),
                    (40, Some(3)),
                    (90, Some(0)),
                ] {
                    let (queries, offsets, total) = match gap {
                        None => batch(&mut rng, count),
                        Some(gap) => compact_batch(&mut rng, count, gap),
                    };
                    let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                    let want = serial_buffer(&letters, &offsets, total);
                    for mode in [1u8, 2] {
                        let got = buffer_with_mode(&letters, &offsets, total, pool, mode);
                        assert_eq!(got, want, "threads {threads} count {count} gap {gap:?}");
                    }
                    if gap.is_none() && count > 1 && pool.threads() > 1 {
                        let (got, on_the_pool) =
                            parallel_buffer(&letters, &offsets, total, pool, &|| {});
                        assert!(on_the_pool, "the production layout took the fallback");
                        assert_eq!(got, want);
                    }
                }
                Ok(())
            })
            .unwrap();
        }
    }

    /// Contexts that do not end inside the block, empty contexts and overlapping spans give the
    /// serial buffer (the serial function stays the reference for them).
    #[cfg(feature = "parallel")]
    #[test]
    fn irregular_layouts_equal_serial_buffer() {
        let mut rng = Rng(0x1234_5678_9abc_def1);
        with_search_pool(3, "blastn", |pool| {
            let (queries, offsets, total) = batch(&mut rng, 30);
            let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
            for mutate in 0..4 {
                let mut offsets = offsets.clone();
                match mutate {
                    // A context past the block (stays sentinels).
                    0 => offsets[7] = total as i32 + 5,
                    // A negative offset (clamped to 0, as the serial function does).
                    1 => offsets[10] = -4,
                    // Queries out of order and overlapping.
                    2 => offsets.swap(6, 20),
                    // A minus context of one query overlapping the plus context of the next.
                    _ => offsets[9] = offsets[10] - 1,
                }
                let want = serial_buffer(&letters, &offsets, total);
                let got = buffer_with_mode(&letters, &offsets, total, pool, 2);
                assert_eq!(got, want, "mutation {mutate}");
                // Out of order or overlapping spans are the fallback (the serial function
                // made the buffer, nothing was encoded on the pool).
                let (got, on_the_pool) = parallel_buffer(&letters, &offsets, total, pool, &|| {});
                assert_eq!(got, want, "mutation {mutate}");
                assert_eq!(on_the_pool, mutate < 2, "mutation {mutate}");
            }
            Ok(())
        })
        .unwrap();
    }

    /// Pool == serial for the Karlin blocks, bit for bit, at 1, 2, 3 and 7 threads on the
    /// contexts of random batches (BLASTNA letters from the buffer, strands in both orders).
    #[cfg(feature = "parallel")]
    #[test]
    fn pool_blocks_equal_serial_blocks() {
        let mut rng = Rng(0xdead_beef_cafe_f00d);
        let (queries, offsets, total) = batch(&mut rng, 70);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let buffer = serial_buffer(&letters, &offsets, total);
        let contexts = contexts_of(&buffer[1..total + 1], &queries, &offsets);
        let mut compared = 0;
        for (reward, penalty) in [(1, -2), (2, -3), (1, -3), (4, -5)] {
            let spec = NuclScoringSpec {
                reward,
                penalty,
                ..spec()
            };
            let want = serial_blocks(&contexts, &spec);
            compared += want.iter().filter(|block| block.is_some()).count();
            for threads in [1usize, 2, 3, 7] {
                with_search_pool(threads, "blastn", |pool| {
                    for take in [contexts.len(), 3, 2, 1, 0] {
                        let part = &contexts[..take];
                        for mode in [1u8, 2] {
                            let got = blocks_with_mode(part, &spec, pool, mode);
                            assert_eq!(got.len(), take);
                            for (i, (got, want)) in got.iter().zip(&want).enumerate() {
                                assert!(
                                    same_bits(got, want),
                                    "threads {threads} take {take} mode {mode} context {i}: {got:?} / {want:?}"
                                );
                            }
                        }
                    }
                    Ok(())
                })
                .unwrap();
            }
        }
        assert!(compared > 200, "blocks compared: {compared}");
    }

    /// The parallel path acts only where a search pool of two or more threads is current.
    #[cfg(feature = "parallel")]
    #[test]
    fn only_a_search_pool_of_two_threads_acts() {
        with_search_pool(1, "blastn", |pool| {
            assert!(!on_pool(pool, 10));
            Ok(())
        })
        .unwrap();
        with_search_pool(3, "blastn", |pool| {
            assert!(on_pool(pool, 10));
            assert!(!on_pool(pool, 1));
            assert!(!on_pool(pool, 0));
            Ok(())
        })
        .unwrap();
    }

    /// With the switch on, both steps run on the threads of the search pool, and on more than one
    /// of them: the probe that each parallel map calls once per item sees a rayon worker of a pool
    /// of three threads, and the items are spread over at least two workers. (Each probe call waits
    /// until a second worker has called it, the first call that waits 10 s in vain ends all
    /// waiting, so that the spread does not depend on the scheduler; a serial map has no second
    /// worker and fails after those 10 s.) With the switch off, on a pool of one thread, or for one
    /// item, the probe is never called.
    #[cfg(feature = "parallel")]
    #[test]
    fn steps_run_on_the_pool_threads() {
        use std::collections::BTreeSet;
        use std::sync::atomic::{AtomicBool, AtomicUsize};
        use std::sync::Mutex;
        use std::time::{Duration, Instant};

        let mut rng = Rng(0x0bad_cafe_1234_5678);
        let (queries, offsets, total) = batch(&mut rng, 30);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let want_buffer = serial_buffer(&letters, &offsets, total);
        let contexts = contexts_of(&want_buffer[1..total + 1], &queries, &offsets);
        let spec = spec();
        let want_blocks = serial_blocks(&contexts, &spec);

        let calls = AtomicUsize::new(0);
        let workers = Mutex::new(BTreeSet::new());
        let in_pool_of_three = AtomicUsize::new(0);
        let gave_up = AtomicBool::new(false);
        let probe = || {
            calls.fetch_add(1, Ordering::SeqCst);
            if let Some(index) = rayon::current_thread_index() {
                workers.lock().unwrap().insert(index);
                if rayon::current_num_threads() == 3 {
                    in_pool_of_three.fetch_add(1, Ordering::SeqCst);
                }
            }
            let start = Instant::now();
            while !gave_up.load(Ordering::SeqCst) && workers.lock().unwrap().len() < 2 {
                if start.elapsed() > Duration::from_secs(10) {
                    gave_up.store(true, Ordering::SeqCst);
                }
                std::thread::yield_now();
            }
        };
        let reset = || {
            calls.store(0, Ordering::SeqCst);
            in_pool_of_three.store(0, Ordering::SeqCst);
            workers.lock().unwrap().clear();
            gave_up.store(false, Ordering::SeqCst);
        };
        // The queries that have a placed context: the probe is called once for each of them.
        let placed = letters
            .iter()
            .zip(offsets.chunks_exact(2))
            .filter(|(letters, pair)| {
                pair.iter()
                    .any(|&offset| place(offset, letters.len(), total).is_some())
            })
            .count();
        with_search_pool(3, "blastn", |pool| {
            for mode in [1u8, 2] {
                reset();
                let got = buffer_with_probe(&letters, &offsets, total, pool, mode, &probe);
                assert_eq!(got, want_buffer, "mode {mode}");
                assert_eq!(calls.load(Ordering::SeqCst), placed, "mode {mode}");
                assert_eq!(
                    in_pool_of_three.load(Ordering::SeqCst),
                    placed,
                    "mode {mode}"
                );
                {
                    let workers = workers.lock().unwrap();
                    assert!(
                        workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                        "mode {mode}: workers {workers:?}"
                    );
                }
                reset();
                let got = blocks_with_probe(&contexts, &spec, pool, mode, &probe);
                assert!(
                    got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)),
                    "mode {mode}"
                );
                assert_eq!(calls.load(Ordering::SeqCst), contexts.len(), "mode {mode}");
                assert_eq!(
                    in_pool_of_three.load(Ordering::SeqCst),
                    contexts.len(),
                    "mode {mode}"
                );
                let workers = workers.lock().unwrap();
                assert!(
                    workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                    "mode {mode}: workers {workers:?}"
                );
            }
            // Not acted on: switch off, one query, one context.
            reset();
            assert_eq!(
                buffer_with_probe(&letters, &offsets, total, pool, 0, &probe),
                want_buffer
            );
            let got = blocks_with_probe(&contexts, &spec, pool, 0, &probe);
            assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
            let one = buffer_with_probe(&letters[..1], &offsets[..2], total, pool, 1, &probe);
            assert_eq!(one, serial_buffer(&letters[..1], &offsets[..2], total));
            let got = blocks_with_probe(&contexts[..1], &spec, pool, 1, &probe);
            assert!(same_bits(&got[0], &want_blocks[0]));
            assert_eq!(calls.load(Ordering::SeqCst), 0);
            Ok(())
        })
        .unwrap();
        with_search_pool(1, "blastn", |pool| {
            reset();
            assert_eq!(
                buffer_with_probe(&letters, &offsets, total, pool, 1, &probe),
                want_buffer
            );
            let got = blocks_with_probe(&contexts, &spec, pool, 1, &probe);
            assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
            assert_eq!(calls.load(Ordering::SeqCst), 0);
            Ok(())
        })
        .unwrap();
    }

    /// A random query of 0 to `max` letters: often empty, one letter or very short, with bytes of
    /// every value, IUPAC letters and lower case (the reviewer's generator, for the production
    /// layout tests).
    fn varied_sequence(rng: &mut Rng, max: usize) -> Vec<u8> {
        let len = match rng.below(7) {
            0 => 0,
            1 => 1,
            2 => rng.below(4) as usize,
            _ => rng.below(max as u64 + 1) as usize,
        };
        (0..len)
            .map(|_| match rng.below(10) {
                0 => rng.below(256) as u8,
                1 => b"NnRYKMSWBDHVrykmswbdhvUu-*"[rng.below(26) as usize],
                2 => b"acgt"[rng.below(4) as usize],
                _ => b"ACGT"[rng.below(4) as usize],
            })
            .collect()
    }

    /// The production layout (offsets of `search_query_batch`, empty queries first or last,
    /// bytes 0..=255) on pools of 2, 3, 5 and 16 threads, in modes 1 and 2, with two score pairs:
    /// the buffer byte for byte, the blocks bit for bit, and never the fallback.
    #[cfg(feature = "parallel")]
    #[test]
    fn production_layout_buffers_and_blocks() {
        let mut rng = Rng(0x5151_7e57_0dd5_eed5);
        let mut cases = 0;
        for threads in [2usize, 3, 5, 16] {
            with_search_pool(threads, "blastn", |pool| {
                for (count, max) in [
                    (2usize, 50usize),
                    (3, 5),
                    (7, 3000),
                    (64, 400),
                    (300, 60),
                    (2_000, 30),
                ] {
                    for _ in 0..2 {
                        let mut queries: Vec<Vec<u8>> =
                            (0..count).map(|_| varied_sequence(&mut rng, max)).collect();
                        if rng.below(3) == 0 {
                            queries.push(Vec::new());
                        }
                        if rng.below(3) == 0 {
                            queries.insert(0, Vec::new());
                        }
                        if queries.iter().all(|query| query.is_empty()) {
                            queries.push(b"ACGT".to_vec());
                        }
                        let (offsets, total) = production_layout(&queries);
                        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                        let want = serial_buffer(&letters, &offsets, total);
                        let (got, on_the_pool) =
                            parallel_buffer(&letters, &offsets, total, pool, &|| {});
                        assert!(on_the_pool, "the production layout took the fallback");
                        assert_eq!(got, want, "threads {threads} count {count}");
                        for mode in [1u8, 2] {
                            let got = buffer_with_mode(&letters, &offsets, total, pool, mode);
                            assert_eq!(got, want, "threads {threads} count {count} mode {mode}");
                        }
                        let contexts = contexts_of(&want[1..total + 1], &queries, &offsets);
                        for (reward, penalty) in [(1, -2), (2, -3)] {
                            let spec = NuclScoringSpec {
                                reward,
                                penalty,
                                gap_open: 0,
                                gap_extend: 0,
                            };
                            let want = serial_blocks(&contexts, &spec);
                            for mode in [1u8, 2] {
                                let got = blocks_with_mode(&contexts, &spec, pool, mode);
                                assert_eq!(got.len(), want.len());
                                assert!(
                                    got.iter().zip(&want).all(|(a, b)| same_bits(a, b)),
                                    "threads {threads} count {count} r {reward} p {penalty}"
                                );
                            }
                        }
                        cases += 1;
                    }
                }
                Ok(())
            })
            .unwrap();
        }
        assert!(cases >= 40, "{cases}");
    }

    /// Called from a plain thread that is not a worker of the pool it is given: the work must still
    /// run on that pool (this is what `SearchPool::install` adds), never on rayon's global pool.
    /// Without `install` in the two parallel steps the probe sees a thread outside the pool, and
    /// this test fails.
    #[cfg(feature = "parallel")]
    #[test]
    fn caller_outside_the_pool_still_runs_on_it() {
        use std::sync::atomic::AtomicUsize;
        let mut rng = Rng(0x7777_1234_abcd_0001);
        let queries: Vec<Vec<u8>> = (0..64).map(|_| varied_sequence(&mut rng, 500)).collect();
        let (offsets, total) = production_layout(&queries);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let want = serial_buffer(&letters, &offsets, total);
        let contexts = contexts_of(&want[1..total + 1], &queries, &offsets);
        let spec = spec();
        with_search_pool(3, "blastn", |pool| {
            std::thread::scope(|scope| {
                scope
                    .spawn(|| {
                        assert!(rayon::current_thread_index().is_none());
                        let wrong = AtomicUsize::new(0);
                        let probe = || {
                            if rayon::current_num_threads() != 3
                                || rayon::current_thread_index().is_none()
                            {
                                wrong.fetch_add(1, Ordering::SeqCst);
                            }
                        };
                        let got = buffer_with_probe(&letters, &offsets, total, pool, 1, &probe);
                        assert_eq!(got, want);
                        let got = blocks_with_probe(&contexts, &spec, pool, 1, &probe);
                        assert_eq!(got.len(), contexts.len());
                        assert_eq!(
                            wrong.load(Ordering::SeqCst),
                            0,
                            "work ran off the search pool"
                        );
                    })
                    .join()
                    .unwrap();
            });
            Ok(())
        })
        .unwrap();
    }

    /// A LOSAT_X_QBATCHPAR-like slot: a plain thread with a private pool of 2 while the outer pool
    /// of 4 exists, and the steps called from inside a parallel loop of the outer pool (nested).
    #[cfg(feature = "parallel")]
    #[test]
    fn slot_pool_and_nested_calls() {
        use rayon::prelude::*;
        use std::sync::atomic::AtomicUsize;
        let mut rng = Rng(0x0123_4567_89ab_cdef);
        let queries: Vec<Vec<u8>> = (0..200).map(|_| varied_sequence(&mut rng, 800)).collect();
        let (offsets, total) = production_layout(&queries);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let want = serial_buffer(&letters, &offsets, total);
        let contexts = contexts_of(&want[1..total + 1], &queries, &offsets);
        let spec = spec();
        let want_blocks = serial_blocks(&contexts, &spec);
        with_search_pool(4, "blastn", |outer| {
            std::thread::scope(|scope| {
                let slot_thread = scope.spawn(|| {
                    with_search_pool(2, "blastn", |slot| {
                        let wrong = AtomicUsize::new(0);
                        let probe = || {
                            if rayon::current_num_threads() != 2 {
                                wrong.fetch_add(1, Ordering::SeqCst);
                            }
                        };
                        let got = buffer_with_probe(&letters, &offsets, total, slot, 2, &probe);
                        assert_eq!(got, want);
                        let got = blocks_with_probe(&contexts, &spec, slot, 2, &probe);
                        assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
                        assert_eq!(wrong.load(Ordering::SeqCst), 0);
                        Ok(())
                    })
                    .unwrap()
                });
                outer.install(|| {
                    (0..6).into_par_iter().for_each(|_| {
                        let got = buffer_with_mode(&letters, &offsets, total, outer, 2);
                        assert_eq!(got, want);
                        let got = blocks_with_mode(&contexts, &spec, outer, 2);
                        assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
                    })
                });
                slot_thread.join().unwrap();
            });
            Ok(())
        })
        .unwrap();
    }

    /// The comparison functions the shadow uses see a length mismatch, one changed byte, and in
    /// the blocks a count mismatch, `None` against `Some` and a change in each of the five fields.
    #[test]
    fn diff_functions_name_what_differs() {
        assert_eq!(first_buffer_diff(&[1, 2, 3], &[1, 2, 3]), None);
        assert!(first_buffer_diff(&[1, 2, 3], &[1, 2])
            .unwrap()
            .contains("lengths 3 / 2"));
        let text = first_buffer_diff(&[1, 2, 3], &[1, 7, 3]).unwrap();
        assert!(text.contains("index 1 (2 / 7)"), "{text}");

        let base = KarlinParams {
            lambda: 0.625,
            k: 0.41,
            h: 0.78,
            alpha: 1.5,
            beta: -2.0,
        };
        assert_eq!(
            first_block_diff(&[Some(base), None], &[Some(base), None]),
            None
        );
        assert!(first_block_diff(&[None], &[None, None])
            .unwrap()
            .contains("counts 1 / 2"));
        let text = first_block_diff(&[None, Some(base)], &[None, None]).unwrap();
        assert!(
            text.contains("context 1") && text.contains("Some / None"),
            "{text}"
        );
        let text = first_block_diff(&[None], &[Some(base)]).unwrap();
        assert!(text.contains("None / Some"), "{text}");
        let variants = [
            (
                "lambda",
                KarlinParams {
                    lambda: 0.626,
                    ..base
                },
            ),
            ("k", KarlinParams { k: 0.42, ..base }),
            ("h", KarlinParams { h: 0.79, ..base }),
            ("alpha", KarlinParams { alpha: 1.6, ..base }),
            ("beta", KarlinParams { beta: -2.1, ..base }),
        ];
        for (name, variant) in variants {
            let text = first_block_diff(&[None, Some(variant)], &[None, Some(base)]).unwrap();
            assert!(
                text.contains("context 1") && text.contains(&format!("field {name} ")),
                "{text}"
            );
        }
    }

    /// The shadow panics, naming the first differing byte, when the buffer made on the pool
    /// differs from the serial one (the test hook changes one encoded byte of the pool's buffer);
    /// without the hook the same call is silent, and with the switch only (mode 1) nothing is
    /// compared.
    #[cfg(feature = "parallel")]
    #[test]
    fn shadow_panics_on_a_changed_buffer_byte() {
        let mut rng = Rng(0x0bad_0001_0002_0003);
        let (queries, offsets, total) = batch(&mut rng, 40);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let want = serial_buffer(&letters, &offsets, total);
        with_search_pool(3, "blastn", |pool| {
            assert_eq!(buffer_with_mode(&letters, &offsets, total, pool, 2), want);
            for index in [0, total / 2, total + 1] {
                TAMPER_BUFFER.with(|cell| cell.set(Some(index)));
                let text = panic_text(|| buffer_with_mode(&letters, &offsets, total, pool, 2));
                assert!(
                    text.contains("LOSAT_X_QPREPPARSHADOW")
                        && text.contains("query buffer")
                        && text.contains(&format!("first differing byte at index {index} (")),
                    "{text}"
                );
            }
            // Mode 1 compares nothing: the changed byte goes through.
            TAMPER_BUFFER.with(|cell| cell.set(Some(total / 2)));
            let got = buffer_with_mode(&letters, &offsets, total, pool, 1);
            assert_ne!(got, want);
            assert_eq!(buffer_with_mode(&letters, &offsets, total, pool, 2), want);
            Ok(())
        })
        .unwrap();
    }

    /// The shadow panics, naming the context and the field, when one field of a block made on the
    /// pool differs by one bit from the serial one (the test hook changes it), for each of the five
    /// fields.
    #[cfg(feature = "parallel")]
    #[test]
    fn shadow_panics_on_a_changed_karlin_field() {
        let mut rng = Rng(0x0bad_0004_0005_0006);
        let (queries, offsets, total) = batch(&mut rng, 40);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let buffer = serial_buffer(&letters, &offsets, total);
        let contexts = contexts_of(&buffer[1..total + 1], &queries, &offsets);
        let spec = spec();
        let want = serial_blocks(&contexts, &spec);
        let first = want.iter().position(Option::is_some).expect("a block");
        with_search_pool(3, "blastn", |pool| {
            let got = blocks_with_mode(&contexts, &spec, pool, 2);
            assert!(got.iter().zip(&want).all(|(a, b)| same_bits(a, b)));
            for (field, name) in ["lambda", "k", "h", "alpha", "beta"]
                .into_iter()
                .enumerate()
            {
                TAMPER_FIELD.with(|cell| cell.set(Some(field)));
                let text = panic_text(|| blocks_with_mode(&contexts, &spec, pool, 2));
                assert!(
                    text.contains("LOSAT_X_QPREPPARSHADOW")
                        && text.contains(&format!("context {first}, field {name} (")),
                    "{text}"
                );
            }
            Ok(())
        })
        .unwrap();
    }
}

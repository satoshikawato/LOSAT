//! EXPERIMENT (LOSAT_X_QPREPPAR): the DUST masks of a query batch made on the search pool, the
//! per-query map split over the pool's threads and collected in query order. Until 2026-11-01 this
//! map had a switch of its own (LOSAT_X_QDUSTPAR) with the same rule as LOSAT_X_QPREPPAR's maps
//! of the same batch (a pool of two or more threads, two or more queries); it was folded into
//! LOSAT_X_QPREPPAR then (line SF S3, family bn-query-read), which makes a batch's setup on the
//! pool: these masks, the encoding and the Karlin blocks (`blast_engine/x_query_prep_par.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:92-102,171-180
//! NCBI makes a fresh `CSymDustMasker` for each query of the batch and masks the queries one
//! after the other. `apply_dust_masking` does the same, serially, in query order. For a batch of
//! thousands of queries that is most of the time before the lookup table is built, and no other
//! thread works meanwhile.
//!
//! With the switch on, and a LOSAT search pool of two or more threads current, the per-query map
//! `DustMasker::new(level, window, linker).mask_sequence(seq)` is an ordered parallel map
//! (`par_iter().map().collect()`, which returns the results in the order of the input; rayon
//! splits the queries into ranges, so the masks are per query but the jobs are not) run with
//! `SearchPool::install`, so on the pool of the search that calls it. Rayon's global pool is never
//! used. A search without a pool (one thread, a build without the `parallel` feature, WebAssembly
//! without `wasm-threads`), and a batch of fewer than two queries, run the serial map. With `-dust no` there is no map at all (`apply_dust_masking`
//! makes no masks).
//!
//! Why the masks are the same: the mask of a query is a function of that query's letters and of
//! the three DUST parameters only. `DustMasker::new` builds `level`, `window`, `linker`, `low_k`
//! and a threshold table from the three parameters and keeps no other state; every call of
//! `mask_sequence` makes its own `NcbiRandom` (the random numbers that resolve `N`),
//! windows and perfect-interval lists as locals, which is the reset of NCBI's fresh
//! `CSymDustMasker` for each query. `utils/dust.rs` has no statics that a call writes. So no query's mask depends on another query or on the order in which they are
//! masked, and the results are collected in query order. The verbose lines, the totals and the
//! empty case stay in `apply_dust_masking`, in the same order. No floating-point operation
//! changes order.
//!
//! Programs that reach it: only BLASTN (megablast, dc-megablast, blastn, blastn-short with
//! `-dust yes`), through `coordination::query_masks` (every query batch of `search_query_batch`)
//! and `chunk_query_masks` (the parts of a split query). The `-db` search of BLASTN uses the same
//! functions, so it is reached too. BLASTP, BLASTX, TBLASTN and TBLASTX do not call
//! `apply_dust_masking`; they mask with SEG.
//!
//! With `LOSAT_TIMING`, `main` prints `[TIMING] query_dust: <seconds>s`, and with `LOSAT_X_STATS`
//! `[X_STATS] QDUSTPAR dust_s=<seconds>`, only when at least one batch was masked (so not for
//! BLASTP, BLASTX, TBLASTN, TBLASTX or `-dust no`), with or without the switch. The seconds are
//! the sum of the wall times of every DUST call of the process (every batch's).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use super::coordination::QueryLetters;
use crate::utils::dust::{DustMasker, MaskedInterval};
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

/// The DUST masks of `seqs`, one list per sequence, in order: the serial map of NCBI's loop.
fn serial<L: QueryLetters + ?Sized>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
) -> Vec<Vec<MaskedInterval>> {
    (0..seqs.count())
        .map(|index| {
            seqs.read(index, |seq| {
                DustMasker::new(level, window, linker).mask_sequence(seq)
            })
        })
        .collect()
}

/// Whether the masks of `seq_count` queries are made on `pool`: two or more queries and a LOSAT
/// search pool of two or more threads.
fn on_pool(pool: &SearchPool<'_>, seq_count: usize) -> bool {
    seq_count > 1 && pool.enabled() && pool.threads() > 1
}

/// The parallel map on `pool`; the same lists as `serial`. Only called when `on_pool`. `probe` is
/// called once per query, on the thread that masks it; production passes `|| {}`, which the
/// compiler removes (a unit test passes a closure that looks at the thread).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    probe: &P,
) -> Vec<Vec<MaskedInterval>> {
    pool.install(|| {
        (0..seqs.count())
            .into_par_iter()
            .map(|index| {
                probe();
                seqs.read(index, |seq| {
                    DustMasker::new(level, window, linker).mask_sequence(seq)
                })
            })
            .collect()
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    _pool: &SearchPool<'_>,
    _probe: &P,
) -> Vec<Vec<MaskedInterval>> {
    serial(seqs, level, window, linker)
}

/// The masks of `seqs` as the switch `on` makes them. Separate from the environment so that tests
/// choose the switch.
pub(crate) fn masks_with_mode<L: QueryLetters + ?Sized>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    on: bool,
) -> Vec<Vec<MaskedInterval>> {
    masks_with_probe(seqs, level, window, linker, pool, on, &|| {})
}

/// `masks_with_mode` with a `probe` that the parallel map calls once per query (see `parallel`).
fn masks_with_probe<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    on: bool,
    probe: &P,
) -> Vec<Vec<MaskedInterval>> {
    if on {
        crate::utils::xreach::x_reach!("LOSAT_X_QPREPPAR", on_pool(pool, seqs.count()));
    }
    if !on || !on_pool(pool, seqs.count()) {
        return serial(seqs, level, window, linker);
    }
    DUST_POOL_CALLS.fetch_add(1, Ordering::Relaxed);
    parallel(seqs, level, window, linker, pool, probe)
}

// No NCBI counterpart: the summed wall time of the DUST calls and their number, for LOSAT_TIMING
// (any switch state), so that a run with the switch and one without can be compared; they change
// no value NCBI computes.
static DUST_NANOS: AtomicU64 = AtomicU64::new(0);
static DUST_CALLS: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: the DUST calls whose masks were made on the pool (the `pool=` of the
// `[X_STATS]` line); it changes no value NCBI computes.
static DUST_POOL_CALLS: AtomicU64 = AtomicU64::new(0);

/// Whether the DUST step is timed: with LOSAT_TIMING (the `[TIMING]` line) or LOSAT_X_STATS (the
/// `[X_STATS]` line, which costs no more than two clock reads per batch; LOSAT_TIMING also times
/// every ungapped extension on shared counters).
fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var_os("LOSAT_TIMING").is_some() || std::env::var_os("LOSAT_X_STATS").is_some()
    })
}

/// The masks of `seqs` as `apply_dust_masking` makes them: the environment's switch decides.
pub(crate) fn masks<L: QueryLetters + ?Sized>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
) -> Vec<Vec<MaskedInterval>> {
    let start = timing_on().then(std::time::Instant::now);
    let masks = masks_with_mode(seqs, level, window, linker, pool, on());
    if let Some(start) = start {
        DUST_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        DUST_CALLS.fetch_add(1, Ordering::Relaxed);
    }
    masks
}

/// With LOSAT_TIMING / LOSAT_X_STATS the summed wall time of every DUST call (not the wall time of
/// the step: see the module doc), printed by `main` at exit, only when at least one batch was
/// masked.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if timing_on() && DUST_CALLS.load(Ordering::Relaxed) > 0 {
        let seconds = DUST_NANOS.load(Ordering::Relaxed) as f64 / 1e9;
        if std::env::var_os("LOSAT_TIMING").is_some() {
            eprintln!("[TIMING] query_dust: {seconds:.3}s");
        }
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!(
                "[X_STATS] QDUSTPAR dust_s={seconds:.3} pool={}",
                DUST_POOL_CALLS.load(Ordering::Relaxed)
            );
        }
    }
}

#[cfg(test)]
#[path = "x_query_dust_par_tests.rs"]
mod tests;

//! EXPERIMENT (LOSAT_X_SEEDBUCKET / LOSAT_X_SEEDBUCKETSHADOW): the two-hit
//! stage of the TBLASTX word finder with the hits grouped by diagonal range.
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-614
//! (`BlastAaWordFinder_TwoHit`): the subject is scanned left to right, and
//! every hit `(q_off, s_off)` is tested against — and updates — exactly one
//! cell of the diagonal table, `diag_array[(q_off - s_off) & diag_mask]`.
//! Nothing else is read or written by the test: the extension that a passing
//! hit triggers reads only the two sequences and the matrix, and the result is
//! appended to the hit list.  So the sequence of states a diagonal cell goes
//! through, and the extensions it triggers, depend only on the hits of that
//! diagonal, taken in scan order.
//!
//! This module keeps that order per diagonal but changes the order *across*
//! diagonals: the hits of a scan window are first appended to buckets, one
//! per range of `2^shift` consecutive diagonals, and then each bucket is
//! processed in turn.  Each bucket touches a contiguous slice of the diagonal
//! table that fits in the first-level cache, instead of the whole table (4 MiB
//! for a 290 kb query, 8 MiB for 660 kb: a last-level-cache access per hit).
//! The idea is DIAMOND's: make the seed-join cache-local by reordering it
//! (double indexing), not by changing what is compared.
//!
//! The hit list the reference produces is in scan order (a hit's extension is
//! appended when the hit is processed).  Every hit carries its position in the
//! scan stream, and the saved HSPs are put back into that order before the
//! stable score sort (`Blast_InitHitListSortByScore`), so that the list the
//! rest of the pipeline sees is the reference's list element for element.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-538,588-606
//! ```c
//! Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//! Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ...
//! diag_coord = (query_offset - subject_offset) & diag_mask;
//! ...
//! if (diag_array[diag_coord].flag) {
//! ...
//!     last_hit = diag_array[diag_coord].last_hit - diag_offset;
//!     diff = subject_offset - last_hit;
//! ...
//!     if (score >= cutoffs->cutoff_score)
//!         BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//!                          query_offset, subject_offset, hsp_len,
//!                          score);
//! ```
//! This is the per-hit body of BlastAaWordFinder_TwoHit that the bucketed order applies. The Rust
//! code runs the same statements on a hit in the same way; only the order across diagonals
//! (which cell is visited next) differs. Within one diagonal the hits arrive in scan order, so
//! each cell goes through the same states and starts the same extensions. This is by
//! construction (integer operations only); LOSAT_X_SEEDBUCKETSHADOW checks it on real runs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-309
//! ```c
//! void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
//! {
//!     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//!           sizeof(BlastInitHSP), score_compare_match);
//! ```
//! NCBI sorts the hit list by score after the scan (called at aa_ungapped.c:282). Equal scores
//! occur, and the order of the result among them follows the order of the input, so the saved HSPs
//! must be in scan order before the sort. `restore_order` puts them back.
//!
//! `LOSAT_X_SEEDBUCKETSHADOW=1` runs both orders on every subject chunk (the
//! bucketed one on a copy of the diagonal table) and aborts when the HSP
//! lists or the diagonal tables differ.
//!
//! The reordering only pays when the table is large (see `min_cells`): with
//! `LOSAT_X_SEEDBUCKET=1` it is used for tables of at least
//! `LOSAT_X_SEEDBUCKET_MIN_CELLS` cells (default 2^21, 8 MiB), and in BLASTX
//! only with `LOSAT_X_BXSEEDBUCKET=1` as well (see `blastx_mode`).

use super::blast_gapalign::InitHSP;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// Subject chunks (TBLASTX) / subjects (BLASTX) compared in shadow mode, and
/// the HSPs they produced.
pub(crate) static SHADOW_CHUNKS: AtomicU64 = AtomicU64::new(0);
pub(crate) static SHADOW_HITS: AtomicU64 = AtomicU64::new(0);

/// Shadow-mode summary (printed by `main` at exit).
// No NCBI counterpart: prints the shadow-mode counters; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] SEEDBUCKET chunks_compared={} hsps_compared={} (HSP lists and diagonal tables identical)",
            SHADOW_CHUNKS.load(Ordering::Relaxed),
            SHADOW_HITS.load(Ordering::Relaxed)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// ```c
// while (scan_range[1] <= scan_range[2]) {
// ...
//     for (i = 0; i < hits; ++i) {
//         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//         diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// Dispatch point: mode 0 runs the scan-order loop (the port of the C loop above), mode 1 the
// bucketed order, mode 2 both, compared.
/// 0 = scan order (reference), 1 = bucketed, 2 = both, compared.
pub(crate) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_SEEDBUCKETSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_SEEDBUCKET").is_some() {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: a size threshold for choosing the order. It changes speed only; it does
// not change any value NCBI computes.
/// Smallest diagonal table, in cells, for which the bucketed order is used
/// when `LOSAT_X_SEEDBUCKET=1` (LOSAT_X_SEEDBUCKET_MIN_CELLS, default 2^21:
/// an 8 MiB table of 4-byte cells).  Below that the table is served well
/// enough by the cache hierarchy that the extra pass over the hits costs more
/// than it saves (measured on this VM: a 4 MiB table, tblastx p03, is 50 %
/// slower bucketed; an 8 MiB table, d04, twice as fast).  The threshold is a
/// property of the machine (STLB reach, L3 size), not of the result.
pub(crate) fn min_cells() -> usize {
    static MIN: OnceLock<usize> = OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_SEEDBUCKET_MIN_CELLS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1 << 21)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// ```c
// while (scan_range[1] <= scan_range[2]) {
// ...
//     for (i = 0; i < hits; ++i) {
//         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//         diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// Dispatch point: chooses between the scan-order loop (the port of the C loop above) and the
// bucketed order, by table size.
/// The order to use for a diagonal table of `cells` cells: the shadow mode
/// compares whatever the size, the plain bucketed mode applies only from
/// `min_cells()` up.
#[inline]
pub(crate) fn mode_for(cells: usize) -> u8 {
    match mode() {
        1 if cells < min_cells() => 0,
        m => m,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// ```c
// while (scan_range[1] <= scan_range[2]) {
// ...
//     for (i = 0; i < hits; ++i) {
//         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//         diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// Dispatch point for BLASTX: the same choice between the scan-order loop and the bucketed order.
/// BLASTX: the bucketed order needs `LOSAT_X_BXSEEDBUCKET=1` as well.  The
/// subjects are single proteins, so a subject's hits are far fewer than the
/// rows of the (query-sized) table and the grouping buys nothing; measured
/// slightly slower (blastx_ap 8.4 → 8.7 s).  The shadow mode still compares.
pub(crate) fn blastx_mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| match mode() {
        1 if std::env::var_os("LOSAT_X_BXSEEDBUCKET").is_none() => 0,
        m => m,
    })
}

// No NCBI counterpart: how many hits are buffered before they are processed. Buffering changes
// only when a hit is processed; the hits of one diagonal stay in scan order.
/// Hits buffered before a flush (LOSAT_X_SEEDBUCKET_BUDGET, default 2^21:
/// 24 MiB of buffered hits).
pub(crate) fn budget() -> usize {
    static BUDGET: OnceLock<usize> = OnceLock::new();
    *BUDGET.get_or_init(|| {
        std::env::var("LOSAT_X_SEEDBUCKET_BUDGET")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1 << 21)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516
// ```c
// diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// A bucket is a range of consecutive values of `diag_coord`, the index of the diagonal table.
/// Diagonal cells per bucket (LOSAT_X_SEEDBUCKET_CELLS, default 32768: a
/// 128 KiB slice of the diagonal table).
fn cells_per_bucket() -> u32 {
    static CELLS: OnceLock<u32> = OnceLock::new();
    *CELLS.get_or_init(|| {
        std::env::var("LOSAT_X_SEEDBUCKET_CELLS")
            .ok()
            .and_then(|v| v.parse().ok())
            .filter(|&v: &u32| v.is_power_of_two() && v >= 256)
            .unwrap_or(32768)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
// ```c
// Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
// Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
// diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// A buffered hit: the (q_off, s_off) of one offset pair, and its position `seq` in the scan
// stream (used to restore the scan order of the saved HSPs).
#[derive(Clone, Copy, Default)]
struct Hit {
    q_off: u32,
    s_off: u32,
    seq: u32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
// ```c
// while (scan_range[1] <= scan_range[2]) {
// ...
//     for (i = 0; i < hits; ++i) {
//         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//         diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// The buckets hold the hits of one scan window (one call of scansub in the C loop), grouped by
// range of diag_coord. Within a bucket the hits are in scan order.
pub(crate) struct SeedBuckets {
    diag_mask: u32,
    shift: u32,
    buckets: Vec<Vec<Hit>>,
    pub(crate) total: usize,
}

impl SeedBuckets {
    /// Buckets of `cells_per_bucket()` consecutive diagonals of a table of
    /// `diag_array_size` cells (a power of two).
    pub(crate) fn new(diag_array_size: u32, diag_mask: u32) -> Self {
        debug_assert!(diag_array_size.is_power_of_two());
        let cells = cells_per_bucket().min(diag_array_size);
        let shift = cells.trailing_zeros();
        let count = (diag_array_size / cells).max(1) as usize;
        let reserve = (budget() / count * 5 / 4).max(1024);
        SeedBuckets {
            diag_mask,
            shift,
            buckets: (0..count).map(|_| Vec::with_capacity(reserve)).collect(),
            total: 0,
        }
    }

    /// Whether a flush is due (budget reached).
    #[inline(always)]
    pub(crate) fn flush_due(&self, budget: usize) -> bool {
        self.total >= budget
    }

    #[inline(always)]
    pub(crate) fn push(&mut self, q_off: u32, s_off: u32, seq: u32) {
        let diag = q_off.wrapping_sub(s_off) & self.diag_mask;
        let b = (diag >> self.shift) as usize;
        // SAFETY: diag <= diag_mask < diag_array_size, so b < count.
        unsafe { self.buckets.get_unchecked_mut(b) }.push(Hit { q_off, s_off, seq });
        self.total += 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
    // ```c
    // for (i = 0; i < hits; ++i) {
    //     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
    //     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
    // ...
    //     diag_coord = (query_offset - subject_offset) & diag_mask;
    // ```
    // The C loop processes the offset pairs in scan order. This processes them bucket by bucket; the
    // hits of one diagonal share a bucket and keep their order, so every diagonal cell sees the
    // same sequence of hits.
    /// Processes every buffered hit, bucket by bucket, each bucket in scan
    /// order, and empties the buckets.
    #[inline(always)]
    pub(crate) fn flush<F: FnMut(u32, u32, u32)>(&mut self, mut process: F) {
        for bucket in self.buckets.iter_mut() {
            for hit in bucket.iter() {
                process(hit.q_off, hit.s_off, hit.seq);
            }
            bucket.clear();
        }
        self.total = 0;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-591
// ```c
// if (score >= cutoffs->cutoff_score)
//     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                      query_offset, subject_offset, hsp_len,
//                      score);
// ```
// NCBI saves an HSP when the hit that produced it is processed, so the saved list is in scan
// order (the sort at blast_extend.c:306-309 follows). This puts the list back into that order.
/// Puts `items` back into increasing `seq_keys` order (the scan order).
pub(crate) fn restore_order<T: Copy>(items: &mut Vec<T>, seq_keys: &[u32]) {
    debug_assert_eq!(items.len(), seq_keys.len());
    let n = items.len();
    if n < 2 {
        return;
    }
    let mut order: Vec<u32> = (0..n as u32).collect();
    // The keys are distinct (one per hit), so a plain sort is a total order.
    order.sort_unstable_by_key(|&i| seq_keys[i as usize]);
    let reordered: Vec<T> = order.iter().map(|&i| items[i as usize]).collect();
    *items = reordered;
}

/// Puts `init_hsps` back into increasing `seq_keys` order (the scan order).
pub(crate) fn restore_scan_order(init_hsps: &mut Vec<InitHSP>, seq_keys: &[u32]) {
    restore_order(init_hsps, seq_keys);
}

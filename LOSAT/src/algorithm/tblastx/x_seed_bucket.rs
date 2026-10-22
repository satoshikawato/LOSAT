//! EXPERIMENT (SEEDBUCKET, graduated on 2026-10-22: no variable): the two-hit
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
//! appended when the hit is processed).  Every buffered hit carries its
//! position among the hits buffered since the last flush, and when a flush
//! ends the HSPs it saved are put back into that order (`restore_flush_order`).
//! The flushes follow each other in scan order (a flush processes every
//! buffered hit before the scan goes on), so the list is in scan order when
//! the scan ends, before the stable score sort
//! (`Blast_InitHitListSortByScore`): the list the rest of the pipeline sees
//! is the reference's list element for element.  The position restarts at
//! every flush, so it is bounded by the budget plus the hits of one scan call
//! (see `budget`) whatever the number of hits of the whole scan.  (Until
//! 2026-10 the position ran over the whole scan in a `u32`, which a query and
//! a subject of a few megabases pass: 2^32 hits per subject-frame scan.)
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
//! construction (integer operations only); `tests/tblastx_seed_bucket.rs` compares real runs with
//! NCBI's output.
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
//! must be in scan order before the sort. `restore_flush_order` puts them back, flush by flush.
//!
//! SEEDBUCKET is a rule without a variable (graduated on 2026-10-22; `utils/xswitch.rs`
//! `GRADUATED`).  The reordering
//! only pays when the table is large (see `min_cells`): it is used for tables
//! of at least `LOSAT_X_SEEDBUCKET_MIN_CELLS` cells (default 2^21, 8 MiB), and
//! never for a traced or debugged run, which keeps the scan order (`scan_chunk_plain`).
//! TBLASTX reaches the scan loop this module acts in only
//! where the unit path (PAIRPAR) does not take the batch (its gates
//! and exactness guards, `x_pair_par.rs` `decide`).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

// No NCBI counterpart: counters for LOSAT_X_STATS that tell which order ran; they do not change
// any value NCBI computes.
/// Scans run in the bucketed order (a subject chunk of one frame), counted where their buckets
/// are made, and the flushes of those scans.
static BUCKETED_SCANS: AtomicU64 = AtomicU64::new(0);
static BUCKET_FLUSHES: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit): with LOSAT_X_STATS the number of scans that ran
/// in the bucketed order and of their flushes.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] SEEDBUCKET bucketed_scans={} flushes={}",
            BUCKETED_SCANS.load(Ordering::Relaxed),
            BUCKET_FLUSHES.load(Ordering::Relaxed)
        );
    }
}

// No NCBI counterpart: a size threshold for choosing the order. It changes speed only; it does
// not change any value NCBI computes.
/// Smallest diagonal table, in cells, for which the bucketed order is used
/// (LOSAT_X_SEEDBUCKET_MIN_CELLS, default 2^21:
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
/// Whether a diagonal table of `cells` cells takes the bucketed order: from `min_cells()` up.
#[inline]
pub(crate) fn bucketed(cells: usize) -> bool {
    let bucketed = cells >= min_cells();
    crate::utils::xreach::x_reach!("LOSAT_X_SEEDBUCKET", bucketed);
    bucketed
}

// No NCBI counterpart: how many hits are buffered before they are processed. Buffering changes
// only when a hit is processed; the hits of one diagonal stay in scan order.
/// Hits buffered before a flush (LOSAT_X_SEEDBUCKET_BUDGET, default 2^21:
/// 24 MiB of buffered hits), at most 2^31.
///
/// The cap bounds a hit's position within its flush (`Hit::seq`, a `u32`): a
/// flush is due as soon as `budget` hits are buffered, and the callers test
/// that after each scan call, whose hit count is an `Int4`.  So at most
/// `budget - 1 + (2^31 - 1) < 2^32` hits are buffered at a time.  The cap is
/// about that arithmetic only: memory (12 bytes per buffered hit, reserved
/// up front) limits a real setting long before it.
pub(crate) fn budget() -> usize {
    static BUDGET: OnceLock<usize> = OnceLock::new();
    *BUDGET.get_or_init(|| {
        capped_budget(
            std::env::var("LOSAT_X_SEEDBUCKET_BUDGET")
                .ok()
                .and_then(|v| v.parse().ok()),
        )
    })
}

// No NCBI counterpart: the largest number of buffered hits the buckets accept as a budget. It
// bounds the bookkeeping number `Hit::seq`; it does not change any value NCBI computes.
/// Largest budget (see `budget`).
const MAX_BUDGET: usize = 1 << 31;

// No NCBI counterpart: the budget for a given setting (none: 2^21), capped at `MAX_BUDGET`.
fn capped_budget(setting: Option<usize>) -> usize {
    setting.unwrap_or(1 << 21).min(MAX_BUDGET)
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
// A buffered hit: the (q_off, s_off) of one offset pair, and `seq`, its position among the hits
// buffered since the last flush (used to restore the scan order of the HSPs that flush saves).
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
// The buckets hold the hits buffered since the last flush (one or more calls of scansub in the C
// loop), grouped by range of diag_coord. Within a bucket the hits are in scan order.
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
        let reserve = ((budget() / count).saturating_mul(5) / 4).max(1024);
        BUCKETED_SCANS.fetch_add(1, Ordering::Relaxed);
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

    /// Appends one hit.  The hits must be pushed in scan order: a hit's `seq`
    /// is the number of hits buffered before it since the last flush.
    #[inline(always)]
    pub(crate) fn push(&mut self, q_off: u32, s_off: u32) {
        let diag = q_off.wrapping_sub(s_off) & self.diag_mask;
        let b = (diag >> self.shift) as usize;
        // No NCBI counterpart: the hit's position among the buffered hits, a bookkeeping number
        // for `restore_flush_order`; it does not change any value NCBI computes. `total` stays
        // below 2^32 between two flushes (see `budget`).
        debug_assert!(self.total < u32::MAX as usize);
        let seq = self.total as u32;
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
    /// order, and empties the buckets.  `process` gets `(q_off, s_off, seq)`;
    /// the caller records `seq` for each HSP it saves and calls
    /// `restore_flush_order` when the flush returns.
    #[inline(always)]
    pub(crate) fn flush<F: FnMut(u32, u32, u32)>(&mut self, mut process: F) {
        BUCKET_FLUSHES.fetch_add(1, Ordering::Relaxed);
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
// order (the sort at blast_extend.c:306-309 follows). A flush saves its HSPs bucket by bucket;
// this puts them back into scan order. Every hit of a flush comes after every hit of the flush
// before it in the scan, so restoring each flush's HSPs in place restores the whole list.
/// Puts `items[start..]`, the HSPs saved by the flush that just ended, back
/// into scan order: `seq_keys[k]` is the position within that flush
/// (`Hit::seq`) of the hit that saved `items[start + k]`.  Empties `seq_keys`
/// for the next flush.
pub(crate) fn restore_flush_order<T: Copy>(items: &mut [T], start: usize, seq_keys: &mut Vec<u32>) {
    let saved = &mut items[start..];
    // One key per saved HSP: the callers push the two together.
    assert_eq!(saved.len(), seq_keys.len());
    let n = saved.len();
    if n >= 2 {
        // A flush holds fewer than 2^32 hits (see `budget`), so `n` fits a `u32`.
        let mut order: Vec<u32> = (0..n as u32).collect();
        // The keys are distinct (a hit saves at most one HSP), so a plain sort is a total order.
        order.sort_unstable_by_key(|&i| seq_keys[i as usize]);
        let reordered: Vec<T> = order.iter().map(|&i| saved[i as usize]).collect();
        saved.copy_from_slice(&reordered);
    }
    seq_keys.clear();
}

#[cfg(test)]
mod tests {
    use super::*;

    // No NCBI counterpart: checks the order bookkeeping of the bucketed two-hit stage; it does
    // not compute any value NCBI computes.
    #[test]
    fn restore_flush_order_sorts_only_the_last_flush() {
        // Three items of earlier flushes, already in scan order, then four saved in bucket order.
        let mut items: Vec<u32> = vec![100, 101, 102, 13, 10, 12, 11];
        let mut keys: Vec<u32> = vec![7, 0, 5, 2];
        restore_flush_order(&mut items, 3, &mut keys);
        assert_eq!(items, vec![100, 101, 102, 10, 11, 12, 13]);
        assert!(keys.is_empty());

        // Keys up to the bound of a flush (2^32 - 3, see `budget`) order as unsigned numbers.
        let mut items: Vec<u32> = vec![1, 3, 2, 0];
        let mut keys: Vec<u32> = vec![3_000_000_000, u32::MAX - 2, 4_000_000_000, 5];
        restore_flush_order(&mut items, 0, &mut keys);
        assert_eq!(items, vec![0, 1, 2, 3]);

        // Nothing or one item saved: unchanged, keys emptied.
        let mut items: Vec<u32> = vec![5, 4];
        let mut keys: Vec<u32> = Vec::new();
        restore_flush_order(&mut items, 2, &mut keys);
        assert_eq!(items, vec![5, 4]);
        let mut keys: Vec<u32> = vec![9];
        restore_flush_order(&mut items, 1, &mut keys);
        assert_eq!(items, vec![5, 4]);
        assert!(keys.is_empty());
    }

    // No NCBI counterpart: drives the buckets the way the two callers do (push the hits of a
    // scan call, flush when the budget is reached, restore after each flush) and checks that the
    // saved items come out in scan order although the position restarts at every flush.
    #[test]
    fn saved_items_are_in_scan_order_across_flushes() {
        // 2^20 cells: several buckets (32 with the default bucket size), so a flush hands the
        // hits out in an order that is not the scan order.
        let diag_array_size: u32 = 1 << 20;
        let mut buckets = SeedBuckets::new(diag_array_size, diag_array_size - 1);
        assert!(
            buckets.buckets.len() > 1,
            "the test needs more than one bucket"
        );
        let budget = 1000usize;
        // Flushes whose saved items were out of scan order before the restore.
        let mut flushes_out_of_order = 0u32;
        // (scan position, per-flush position) of every saved hit.
        let mut saved: Vec<(u64, u32)> = Vec::new();
        let mut seq_keys: Vec<u32> = Vec::new();
        let mut flushes = 0u32;
        let mut largest_seq = 0u32;
        let mut scan_position = 0u64;
        // A deterministic stream: 300 scan calls of 0 to 96 hits on scattered diagonals.
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // The hits buffered since the last flush, in scan order, with their scan position.
        let mut pending: Vec<(u32, u32, u64)> = Vec::new();
        for call in 0..300 {
            let hits = (next() % 97) as usize;
            for _ in 0..hits {
                let q_off = (next() % (1 << 20)) as u32;
                let s_off = (next() % (1 << 20)) as u32;
                buckets.push(q_off, s_off);
                pending.push((q_off, s_off, scan_position));
                scan_position += 1;
            }
            let last_call = call == 299;
            if buckets.flush_due(budget) || last_call {
                let start = saved.len();
                let flushed = buckets.total;
                let mut seen = 0usize;
                buckets.flush(|q, s, seq| {
                    // The per-flush number is the index into this flush's hits.
                    let (pq, ps, position) = pending[seq as usize];
                    assert_eq!((pq, ps), (q, s));
                    seen += 1;
                    largest_seq = largest_seq.max(seq);
                    // Save about a third of the hits.
                    if position % 3 == 0 {
                        saved.push((position, seq));
                        seq_keys.push(seq);
                    }
                });
                assert_eq!(seen, flushed);
                if saved[start..].windows(2).any(|w| w[0].0 > w[1].0) {
                    flushes_out_of_order += 1;
                }
                restore_flush_order(&mut saved, start, &mut seq_keys);
                pending.clear();
                flushes += 1;
            }
        }
        assert!(flushes > 5, "the stream must span several flushes");
        assert!(
            flushes_out_of_order > 5,
            "the buckets must have changed the order within several flushes"
        );
        // The position restarted: it never reached the number of hits of the whole stream.
        assert!(u64::from(largest_seq) < scan_position / 4);
        assert!((largest_seq as usize) < budget + 96);
        // Every third hit was saved, in scan order.
        let positions: Vec<u64> = saved.iter().map(|&(p, _)| p).collect();
        let expected: Vec<u64> = (0..scan_position).filter(|p| p % 3 == 0).collect();
        assert_eq!(positions, expected);
    }

    // No NCBI counterpart: the cap that keeps the per-flush position within a u32.
    #[test]
    fn budget_cap_keeps_the_flush_position_in_u32() {
        assert_eq!(capped_budget(None), 1 << 21);
        assert_eq!(capped_budget(Some(4096)), 4096);
        assert_eq!(capped_budget(Some(MAX_BUDGET)), MAX_BUDGET);
        assert_eq!(capped_budget(Some(usize::MAX)), MAX_BUDGET);
        // budget - 1 buffered hits, then one scan call of at most i32::MAX hits.
        let most_buffered = (capped_budget(Some(usize::MAX)) as u64 - 1) + i32::MAX as u64;
        assert!(most_buffered < u64::from(u32::MAX));
    }
}

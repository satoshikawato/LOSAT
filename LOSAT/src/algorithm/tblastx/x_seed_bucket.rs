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
pub fn print_shadow_stats() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] SEEDBUCKET chunks_compared={} hsps_compared={} (HSP lists and diagonal tables identical)",
            SHADOW_CHUNKS.load(Ordering::Relaxed),
            SHADOW_HITS.load(Ordering::Relaxed)
        );
    }
}

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

#[derive(Clone, Copy, Default)]
struct Hit {
    q_off: u32,
    s_off: u32,
    seq: u32,
}

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

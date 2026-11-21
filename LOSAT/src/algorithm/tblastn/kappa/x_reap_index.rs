//! No NCBI counterpart: `LOSAT_X_TBNREAPIDX`. The reap of the HSPs that one composition redo
//! gives a match (`s_HitlistReapContained`, `blast_kappa.c` 224-273; `reap_contained_postredo_hsps`
//! in `kappa.rs`) answered from an index of the HSPs kept so far by subject start.
//!
//! NCBI tests every HSP against every earlier HSP still kept, so a match with tens of thousands of
//! HSPs (a protein against a human chromosome) costs the square of that: the B11h probe (35 human
//! proteins x the human genome `-db`, 8 threads) counted 21.11 G iterations in 73.6 thread-s, 10 %
//! of the search's user time, and the index of this file visited 6,326 candidates in 0.30
//! thread-s with the same answers. The answer for an HSP is whether some earlier HSP still kept,
//! of the same context and frame, contains both of its end points (`CONTAINED_IN_HSP`,
//! `blast_hits_priv.h` 65-70) and scores at least as much; it does not depend on the order the
//! earlier HSPs are met in, so the index may meet them in any order. Lists shorter than
//! `X_REAP_INDEX_FROM` HSPs keep NCBI's loop.
use std::collections::BTreeSet;
use std::sync::OnceLock;

use super::GappedHsp;

/// The list length from which the reap takes the index.
pub(super) const X_REAP_INDEX_FROM: usize = 32;

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// Whether the reap takes the index (`LOSAT_X_TBNREAPIDX`, read once; a unit test sets it per test
/// thread).
fn x_reap_idx_on() -> bool {
    #[cfg(test)]
    if let Some(on) = tests::X_ON.get() {
        return on;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNREAPIDX"))
}

/// Whether the reap of a list of `len` HSPs takes the index (the switch on and `len` at least
/// `X_REAP_INDEX_FROM`); records the side with `LOSAT_X_REACH`.
pub(super) fn x_indexed(len: usize) -> bool {
    if !x_reap_idx_on() {
        return false;
    }
    let indexed = len >= X_REAP_INDEX_FROM;
    crate::utils::xreach::x_reach!("LOSAT_X_TBNREAPIDX", indexed);
    indexed
}

/// `CONTAINED_IN_HSP` of both end points of `inner` in `outer`, and `inner.score <= outer.score`
/// (the test of `reap_keep` in `kappa.rs`, without the context and frame).
#[inline]
fn contains(outer: &GappedHsp, inner: &GappedHsp) -> bool {
    outer.q_start <= inner.q_start
        && outer.q_end >= inner.q_start
        && outer.s_start <= inner.s_start
        && outer.s_end >= inner.s_start
        && outer.q_start <= inner.q_end
        && outer.q_end >= inner.q_end
        && outer.s_start <= inner.s_end
        && outer.s_end >= inner.s_end
        && inner.score <= outer.score
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:224-273
// The same test, over the earlier HSPs still kept of the same context and frame whose subject start
// lies in [max(s_start, s_end) - dmax, min(s_start, s_end)] only (dmax: the largest `s_end -
// s_start` of an HSP kept so far, at least 0): an HSP `outer` that contains both end points has
// s_start <= both of them and s_end >= both, so 0 <= s_end - s_start <= dmax and s_start >=
// max(s_start, s_end) of the HSP - dmax. Every HSP in that range is tested with the whole
// condition, the score included.
/// Whether each HSP of `hsps` is kept by the reap (`reap_keep` in `kappa.rs`), from the index.
#[inline(never)]
pub(super) fn x_reap_keep(hsps: &[(usize, GappedHsp)]) -> Vec<bool> {
    let mut keep = vec![true; hsps.len()];
    // (context, frame, subject start, position in `hsps`) of each HSP kept so far.
    let mut by_start: BTreeSet<(usize, i8, i32, usize)> = BTreeSet::new();
    let mut dmax = 0i64;
    for (read, &(context, inner)) in hsps.iter().enumerate() {
        let lo = (i64::from(inner.s_start.max(inner.s_end)) - dmax)
            .clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
        let hi = inner.s_start.min(inner.s_end);
        // Lowest start first, the order the B11h probe counted; the answer does not depend on the
        // order.
        let contained = lo <= hi
            && by_start
                .range((context, inner.frame, lo, 0)..=(context, inner.frame, hi, usize::MAX))
                .any(|&(_, _, _, earlier)| contains(&hsps[earlier].1, &inner));
        if contained {
            keep[read] = false;
        } else {
            by_start.insert((context, inner.frame, inner.s_start, read));
            dmax = dmax.max(i64::from(inner.s_end) - i64::from(inner.s_start));
        }
    }
    keep
}

#[cfg(test)]
#[path = "x_reap_index_tests.rs"]
mod tests;

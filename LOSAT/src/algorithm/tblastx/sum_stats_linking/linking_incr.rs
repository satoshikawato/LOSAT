//! Incremental replacement of the per-group link_hsps kernel
//! (`LOSAT_LINK_FAST=2`; the default path is unchanged). It ports
//! s_BlastEvenGapLinkHSPs for one group as `linking_index.rs` does - the same
//! rounds, "current max" selection, add-back, E-values and removal - but a
//! recompute pass after the first computes again only the HSPs whose values
//! can have changed since the previous pass.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,584-589
//! ```c
//! s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
//! ...
//!       first_pass=1;    /* do full search */
//!       path_changed=1;
//! ...
//!       while (number_of_hsps > 0)
//! ```
//! The rounds of s_BlastEvenGapLinkHSPs for one group (one query frame and one subject frame
//! sign), with the same control flow as link_hsp_group_ncbi.
//!
//! # What a pass computes
//!
//! A recompute pass gives every remaining HSP H, for each ordering method,
//! the sum `score(H) - cutoff + max(0, largest candidate sum)`, its number of
//! HSPs and xsum, and the link to the candidate with the largest (sum, list
//! index) (`linking_index.rs`, "What a search returns"; NCBI's index-1 reuse,
//! start value and `next_larger` select the same HSP). The sums are added in
//! i64 ("Sums" below), so these are the exact chain maxima over the remaining
//! HSPs.
//!
//! # Which HSPs a later pass visits
//!
//! Between two passes HSPs are only removed. The candidates of H can only
//! disappear and, by induction in list order, no value increases
//! (`linking_index.rs`, step 2). Let P be the HSP that H selected in the
//! previous pass.
//!
//! 1. P remains and its sum did not change: every other candidate has a sum
//!    at most its previous one, which was at most P's, with a smaller list
//!    index on a tie, so H selects P again. H's values change only if P's
//!    number of HSPs, xsum or link changed.
//! 2. H selected nothing: no candidate had a positive sum, and none has now.
//!
//! So a pass computes again only the HSPs whose selected HSP was removed since
//! the previous pass and, recursively, those whose selected HSP changed in
//! this pass, each method on its own (the two methods do not read each
//! other's values). A selected HSP precedes H in the list, so visiting these
//! HSPs in list order reads final values; in case 1 H keeps P without a
//! search. Every other remaining HSP keeps its values and its link.
//!
//! # NCBI's other state
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:685,764-766,964-970
//! ```c
//!       H->linked_to = 0;
//! ...
//!             if(H_hsp_link)
//!                ((LinkHSPStruct*)H_hsp_link)->linked_to++;
//! ...
//! if (best[ordering_method]->linked_to>0) path_changed=1;
//! ...
//!    if (H->linked_to>1) path_changed=1;
//!    H->linked_to=-1000;
//!    H->hsp_link.changed=1;
//! ```
//! A pass sets `linked_to` of every remaining HSP to the number of links of
//! the remaining HSPs into it; removals do not change it until the next pass.
//! `indeg` keeps that number: a pass takes out the links of the HSPs removed
//! since the previous pass and moves the links that change. NCBI's `changed`
//! flags only decide its index-1 reuse, which selects the same HSPs (above).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:901-937
//! ```c
//! best[0]->hsp_link.sum[0] +=
//!    (best[0]->hsp_link.num[0])*cutoff[0];
//! ```
//! The add-back changes a stored sum until the next pass computes it again,
//! and the stored sums are read only by the "current max" scan. Here the
//! stored sums are the leaves of the key tree; `sum` holds the pass values,
//! which the searches read, and a pass writes the pass value back into the
//! leaf of every HSP that received an add-back.
//!
//! # Index 1
//!
//! Every pass sweeps the remaining HSPs with score > cutoff[1] by decreasing
//! `q_end_trim` and fills the prefix-maximum Fenwick tree of
//! `linking_index.rs`, but computes again only the HSPs that need it. A
//! selected HSP has a larger `q_end_trim` than the HSPs that select it, so
//! it is final when they are visited. (The removal of a chain changes the
//! values of a large share of the HSPs downstream of it at index 1, so a
//! search per HSP in a dynamic 2-d structure costs more than the sweep.)
//!
//! # Sums
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:753,868
//! ```c
//! Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
//! ```
//! NCBI adds the chain sums in Int4 without a range check; beyond Int4 the C
//! addition is undefined behaviour (in practice it wraps, and a large chain
//! can lose to a small one). This kernel adds them in i64 and does not
//! reproduce the wrap (Owner decision 2026-10-10). Within Int4 the two agree;
//! the largest pass value of the 16 tblastx comparison rows is 0.04% of
//! INT4_MAX.
//!
//! The trees keep a value and a list index in one key (`LinkKey`): a u64 with
//! the value in 40 bits and the list index + 1 in 24 when the group has fewer
//! than 2^24 - 1 HSPs and every pass value is below 2^39 - 1 (`pass_bound`),
//! else a u128 with 64 bits each. Both run the same kernel, for any group
//! below NCBI's Int4 count of HSPs. Pass values and the values of HSPs without
//! a link (`score - cutoff`, above -2^32) fit the key. A value after the
//! add-back is clamped to the range of the key, which does not change what the
//! kernel selects:
//!
//! - Only best[index] receives an add-back, and between two passes no other
//!   value grows. A value reaches the top only while its HSP is the maximum;
//!   it then stays above every other value until the HSP is removed or the
//!   next pass writes its pass value back, so no two values are at the top.
//! - A value at the bottom is below `-cutoff[index]`, the start value of the
//!   "current max" scan, so its HSP is not selected again before the next
//!   pass.
//! - The predecessor searches read the pass values (`sum`), not these values.
//!
//! # Checking
//!
//! `LOSAT_LINK_FAST_VERIFY=1` compares, after every pass, the values, links
//! and `indeg` of every remaining HSP with a pass computed from scratch by
//! NCBI's scans. `LOSAT_LINK_FAST_SHADOW=1` (in `linking.rs`) also runs the
//! literal NCBI port on every group and compares the complete result.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::sync::OnceLock;

use rustc_hash::FxHashMap;

use crate::algorithm::tblastx::chaining::UngappedHit;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::stats::sum_statistics::{gap_decay_divisor, ncbi_large_gap_sum_e, small_gap_sum_e};

use super::linking_index::{stats_enabled, MaxTree, NONE, TRIM_SIZE, WINDOW_SIZE};
use super::params::LinkHspCutoffs;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:610-623
// ```c
// Int4 sum0=H->hsp_link.sum[0];
// Int4 sum1=H->hsp_link.sum[1];
// if(sum0>=max0)
// {
//    max0=sum0;
//    best[0]=H;
// }
// ```
// NCBI keeps `best[index]` as the last HSP in list order with the largest sum (`>=`). A key
// (value, index + 1) has the same order: larger value first, then larger list index.
/// A tree key: a value and a list index packed so that the integer order is
/// larger value first, then larger index. `Default` (0) is "no entry". A key
/// holds values in `MIN_VALUE..=MAX_VALUE` (module documentation, "Sums").
pub(super) trait LinkKey: Copy + Ord + Default {
    const MIN_VALUE: i64;
    const MAX_VALUE: i64;
    fn new(value: i64, idx: usize) -> Self;
    fn idx(self) -> usize;
    fn value(self) -> i64;
}

/// A 40-bit value (offset by 2^39) and a 24-bit list index + 1.
impl LinkKey for u64 {
    const MIN_VALUE: i64 = -(1 << 39);
    const MAX_VALUE: i64 = (1 << 39) - 1;

    #[inline(always)]
    fn new(value: i64, idx: usize) -> u64 {
        (((value - Self::MIN_VALUE) as u64) << 24) | ((idx as u64) + 1)
    }

    #[inline(always)]
    fn idx(self) -> usize {
        ((self & 0xff_ffff) - 1) as usize
    }

    #[inline(always)]
    fn value(self) -> i64 {
        (self >> 24) as i64 + Self::MIN_VALUE
    }
}

/// A 64-bit value (sign bit flipped) and a 64-bit list index + 1.
impl LinkKey for u128 {
    const MIN_VALUE: i64 = i64::MIN;
    const MAX_VALUE: i64 = i64::MAX;

    #[inline(always)]
    fn new(value: i64, idx: usize) -> u128 {
        ((((value as u64) ^ (1 << 63)) as u128) << 64) | ((idx as u128) + 1)
    }

    #[inline(always)]
    fn idx(self) -> usize {
        ((self as u64) - 1) as usize
    }

    #[inline(always)]
    fn value(self) -> i64 {
        (((self >> 64) as u64) ^ (1 << 63)) as i64
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:753,868
// ```c
// Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
// ```
// No NCBI counterpart for the bound: it only chooses the key of this kernel (module documentation,
// "Sums").
/// A bound of every pass value: a chain sum adds `score - cutoff > 0` over
/// distinct HSPs, so it is at most the sum of the positive `score - cutoff`
/// (for the method with the larger sum).
pub(super) fn pass_bound(group_hits: &[UngappedHit], cutoffs: &LinkHspCutoffs) -> i64 {
    [cutoffs.cutoff_small_gap, cutoffs.cutoff_big_gap]
        .iter()
        .map(|&c| {
            group_hits
                .iter()
                .map(|hit| (i64::from(hit.raw_score) - i64::from(c)).max(0))
                .sum::<i64>()
        })
        .max()
        .unwrap_or(0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:607-623,692,772
// ```c
// max0 = -cutoff[0];
// max1 = -cutoff[1];
// for (H=hp_start->next; H!=NULL; H=H->next) {
// ...
//    if(sum0>=max0)
//    {
//       max0=sum0;
//       best[0]=H;
//    }
// ...
// maxscore = -cutoff[index];
// ```
// NCBI starts both maxima at -cutoff[index], so best[index] stays NULL when every sum is below it.
// The root of the key tree is the largest (value, index); it is best[index] when its value reaches
// the start value.
/// `best[0]` and `best[1]` from the root of the key tree (`None` when the
/// largest value is below `-cutoff[index]`).
#[inline]
fn best_of_root<K: LinkKey>(
    root: [K; 2],
    c: [i32; 2],
    ignore_small_gaps: bool,
) -> [Option<usize>; 2] {
    let mut best = [None, None];
    for m in usize::from(ignore_small_gaps)..2 {
        if root[m] != K::default() && root[m].value() >= -i64::from(c[m]) {
            best[m] = Some(root[m].idx());
        }
    }
    best
}

// No NCBI counterpart: options of the check (verify) and of the key (wide_keys).
/// How the kernel runs. `from_env` is what a search uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct IncrOptions {
    /// Compare every pass with a pass computed from scratch.
    pub verify: bool,
    /// Use the u128 key whatever the group (`LOSAT_LINK_WIDE_KEYS=1`, to
    /// measure its cost; the result is the same).
    pub wide_keys: bool,
}

impl IncrOptions {
    pub(super) fn from_env() -> Self {
        static OPTIONS: OnceLock<IncrOptions> = OnceLock::new();
        *OPTIONS.get_or_init(|| IncrOptions {
            verify: std::env::var_os("LOSAT_LINK_FAST_VERIFY").is_some(),
            wide_keys: std::env::var_os("LOSAT_LINK_WIDE_KEYS").is_some(),
        })
    }
}

// No NCBI counterpart: counters of what the kernel did; they do not change any value NCBI computes.
/// What one call did. Tests assert on it; `LOSAT_LINK_STATS=1` prints it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct IncrStats {
    pub rounds: u64,
    pub passes: u64,
    /// HSPs removed before a later pass (their links are taken out then).
    pub removed: u64,
    /// Index 0 / index 1: HSPs computed again, of which kept their selected
    /// HSP without a search, and searched.
    pub visited0: u64,
    pub kept0: u64,
    pub searched0: u64,
    pub visited1: u64,
    pub kept1: u64,
    pub searched1: u64,
    /// Searches that selected the previous choice again.
    pub same0: u64,
    pub same1: u64,
    /// HSPs read by the index-1 sweeps of the passes after the first.
    pub swept: u64,
    /// Sweep searches answered by NCBI's scan (the tree returned a later HSP).
    pub fallbacks: u64,
    pub verified_passes: u64,
    /// Groups linked with the u128 key, and add-backs clamped to the range of
    /// the key (module documentation, "Sums").
    pub wide_groups: u64,
    pub clamped: u64,
    /// Largest number of HSPs of a pass value (exact, before NCBI's Int2).
    pub max_num: i64,
    pub max_pass_sum: i64,
}

impl IncrStats {
    #[cfg(test)]
    pub(super) fn add(&mut self, o: &IncrStats) {
        self.rounds += o.rounds;
        self.passes += o.passes;
        self.removed += o.removed;
        self.visited0 += o.visited0;
        self.kept0 += o.kept0;
        self.searched0 += o.searched0;
        self.visited1 += o.visited1;
        self.kept1 += o.kept1;
        self.searched1 += o.searched1;
        self.same0 += o.same0;
        self.same1 += o.same1;
        self.swept += o.swept;
        self.fallbacks += o.fallbacks;
        self.verified_passes += o.verified_passes;
        self.wide_groups += o.wide_groups;
        self.clamped += o.clamped;
        self.max_num = self.max_num.max(o.max_num);
        self.max_pass_sum = self.max_pass_sum.max(o.max_pass_sum);
    }
}

// No NCBI counterpart: reverse links. NCBI keeps only the forward link of each HSP; these lists
// find the HSPs that selected a given HSP.
/// For one ordering method, the HSPs whose link is a given HSP, as
/// intrusive doubly linked lists.
struct Children {
    first: Vec<u32>,
    next: Vec<u32>,
    prev: Vec<u32>,
}

impl Children {
    fn new(n: usize) -> Self {
        Self {
            first: vec![NONE; n],
            next: vec![NONE; n],
            prev: vec![NONE; n],
        }
    }

    #[inline]
    fn insert(&mut self, parent: usize, child: usize) {
        let head = self.first[parent];
        self.next[child] = head;
        self.prev[child] = NONE;
        if head != NONE {
            self.prev[head as usize] = child as u32;
        }
        self.first[parent] = child as u32;
    }

    #[inline]
    fn remove(&mut self, parent: usize, child: usize) {
        let p = self.prev[child];
        let x = self.next[child];
        if p != NONE {
            self.next[p as usize] = x;
        } else {
            self.first[parent] = x;
        }
        if x != NONE {
            self.prev[x as usize] = p;
        }
        self.next[child] = NONE;
        self.prev[child] = NONE;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,589
// ```c
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
// ...
//       while (number_of_hsps > 0)
// ```
// link_hsp_group_incr only wraps link_hsp_group_incr_with: it adds the LOSAT_LINK_STATS line and
// changes no value of the result.
/// The kernel `linking.rs` calls under `LOSAT_LINK_FAST=2`.
#[allow(clippy::too_many_arguments)]
pub(super) fn link_hsp_group_incr(
    group_hits: Vec<UngappedHit>,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
) -> Vec<UngappedHit> {
    let n = group_hits.len();
    let stats_on = stats_enabled() && n > 0;
    let bound = if stats_on {
        pass_bound(&group_hits, cutoffs)
    } else {
        0
    };
    let started = stats_on.then(std::time::Instant::now);
    let (group_hits, stats) = link_hsp_group_incr_with(
        group_hits,
        cutoffs,
        gap_decay_rate,
        subject_len_nucl,
        query_contexts,
        length_adj_per_context,
        eff_searchsp_per_context,
        log_k_by_ctx,
        IncrOptions::from_env(),
    );
    if let Some(started) = started {
        eprintln!(
            "[LINK_INCR_STATS] n={} us={} bound={} max_pass_sum={} max_num={} rounds={} passes={} removed={} visited0={} kept0={} searched0={} visited1={} kept1={} searched1={} same0={} same1={} swept={} fallbacks={} verified_passes={} wide={} clamped={}",
            n,
            started.elapsed().as_micros(),
            bound,
            stats.max_pass_sum,
            stats.max_num,
            stats.rounds,
            stats.passes,
            stats.removed,
            stats.visited0,
            stats.kept0,
            stats.searched0,
            stats.visited1,
            stats.kept1,
            stats.searched1,
            stats.same0,
            stats.same1,
            stats.swept,
            stats.fallbacks,
            stats.verified_passes,
            stats.wide_groups,
            stats.clamped
        );
    }
    group_hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:702-745
// ```c
// if (H->hsp->score > cutoff[index]) {
// ...
//    for (H2_index=H_index-1; H2_index>1; H2_index=H2_index-1)
// ...
//       b1 = q_off_t <= H_query_etrim;
//       b2 = s_off_t <= H_sub_etrim;
// ...
//       if (b1|b2|b5|b4) continue;
// ...
//       if (sum>H_hsp_sum)
// ```
// Check only: one pass computed from scratch with NCBI's scans (nearest first, strictly larger
// sum, start value 0), for the remaining HSPs in list order, and the number of links into each.
/// `LOSAT_LINK_FAST_VERIFY`: (sum, num, xsum bits, link) of both methods and
/// the number of links into every HSP, for the remaining HSPs, by NCBI's
/// scans.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
fn verify_pass(
    alive: &[bool],
    score: &[i32],
    qo: &[i32],
    so: &[i32],
    qe: &[i32],
    se: &[i32],
    sl: &[f64],
    lk: &[f64],
    c: [i32; 2],
    ignore_small_gaps: bool,
) -> ([Vec<(i64, i16, u64, u32)>; 2], Vec<i32>) {
    let n = alive.len();
    let mut out: [Vec<(i64, i16, u64, u32)>; 2] = [
        (0..n)
            .map(|i| {
                (
                    i64::from(score[i]) - i64::from(c[0]),
                    1,
                    (sl[i] - lk[i]).to_bits(),
                    NONE,
                )
            })
            .collect(),
        (0..n)
            .map(|i| {
                (
                    i64::from(score[i]) - i64::from(c[1]),
                    1,
                    (sl[i] - lk[i]).to_bits(),
                    NONE,
                )
            })
            .collect(),
    ];
    let mut indeg = vec![0i32; n];
    for m in usize::from(ignore_small_gaps)..2 {
        for i in 0..n {
            if !alive[i] || score[i] <= c[m] {
                continue;
            }
            let (mut h_sum, mut h_num, mut h_xsum, mut h_link) = (0i64, 0i16, 0.0f64, NONE);
            for j in (0..i).rev() {
                if !alive[j] {
                    continue;
                }
                let ok = if m == 0 {
                    qo[j] > qe[i]
                        && so[j] > se[i]
                        && qo[j] <= qe[i] + WINDOW_SIZE
                        && so[j] <= se[i] + WINDOW_SIZE
                } else {
                    qo[j] > qe[i] && so[j] > se[i]
                };
                if ok && out[m][j].0 > h_sum {
                    h_sum = out[m][j].0;
                    h_num = out[m][j].1;
                    h_xsum = f64::from_bits(out[m][j].2);
                    h_link = j as u32;
                }
            }
            let new_sum = h_sum + (i64::from(score[i]) - i64::from(c[m]));
            out[m][i] = (
                new_sum,
                h_num.wrapping_add(1),
                (h_xsum + sl[i] - lk[i]).to_bits(),
                h_link,
            );
            if h_link != NONE {
                indeg[h_link as usize] += 1;
            }
        }
    }
    (out, indeg)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,584-589
// ```c
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
// ...
//       first_pass=1;    /* do full search */
//       path_changed=1;
// ...
//       while (number_of_hsps > 0)
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:158
// ```c
//    Int4 hspcnt; /**< Number of HSPs saved */
// ```
// NCBI counts the HSPs of a list in Int4. Below that count every pass value fits i64 (at most
// n * 2^32 < 2^63). The u64 key when the group fits it, else the u128 key (module documentation,
// "Sums"); both run the same kernel.
/// Links one group with the u64 key when it fits, else with the u128 key.
#[allow(clippy::too_many_arguments)]
pub(super) fn link_hsp_group_incr_with(
    group_hits: Vec<UngappedHit>,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
    options: IncrOptions,
) -> (Vec<UngappedHit>, IncrStats) {
    let n = group_hits.len();
    assert!(
        n < 1 << 31,
        "link_hsps: a group of {n} HSPs (NCBI counts HSPs in Int4)"
    );
    let wide = options.wide_keys
        || n >= (1 << 24) - 1
        || pass_bound(&group_hits, cutoffs) >= <u64 as LinkKey>::MAX_VALUE;
    if wide {
        let (group_hits, mut stats) = link_group::<u128>(
            group_hits,
            cutoffs,
            gap_decay_rate,
            subject_len_nucl,
            query_contexts,
            length_adj_per_context,
            eff_searchsp_per_context,
            log_k_by_ctx,
            options,
        );
        stats.wide_groups = 1;
        (group_hits, stats)
    } else {
        link_group::<u64>(
            group_hits,
            cutoffs,
            gap_decay_rate,
            subject_len_nucl,
            query_contexts,
            length_adj_per_context,
            eff_searchsp_per_context,
            log_k_by_ctx,
            options,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,584-589
// ```c
// s_BlastEvenGapLinkHSPs(EBlastProgramType program_number, BlastHSPList* hsp_list,
// ...
//       first_pass=1;    /* do full search */
//       path_changed=1;
// ...
//       while (number_of_hsps > 0)
// ```
// The replacement for the per-group body of s_BlastEvenGapLinkHSPs: same rounds, same selection
// rule, same removal; a pass after the first computes again only the HSPs described in the
// module documentation.
#[allow(clippy::too_many_arguments)]
fn link_group<K: LinkKey>(
    mut group_hits: Vec<UngappedHit>,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
    options: IncrOptions,
) -> (Vec<UngappedHit>, IncrStats) {
    let mut stats = IncrStats::default();
    let n = group_hits.len();
    if n == 0 {
        return (group_hits, stats);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:559-571
    // ```c
    // query_length = query_info->contexts[query_context].query_length;
    // query_length = MAX(query_length - length_adjustment, 1);
    // ...
    // if (Blast_SubjectIsTranslated(program_number))
    // {
    //    length_adjustment /= CODON_LENGTH;
    //    subject_length /= CODON_LENGTH;
    // }
    // subject_length = MAX(subject_length - length_adjustment, 1);
    // ```
    // Query and subject lengths and the length adjustment, as in link_hsp_group_ncbi.
    let query_context = group_hits[0].ctx_idx;
    let query_len_aa = query_contexts[query_context].aa_len as i64;
    let subject_len_aa = (subject_len_nucl / 3).max(1);
    let length_adjustment = length_adj_per_context[query_context];
    let eff_search_space = eff_searchsp_per_context[query_context];
    let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;
    let length_adj_for_subject = length_adjustment / 3;
    let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;

    let c = [cutoffs.cutoff_small_gap, cutoffs.cutoff_big_gap];
    let gap_prob = cutoffs.gap_prob;
    let ignore_small_gaps = cutoffs.ignore_small_gaps;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:545-550
    // ```c
    // q_length = (hsp->query.end - hsp->query.offset) / 4;
    // s_length = (hsp->subject.end - hsp->subject.offset) / 4;
    // H->q_offset_trim = hsp->query.offset + MIN(q_length, trim_size);
    // H->q_end_trim = hsp->query.end - MIN(q_length, trim_size);
    // H->s_offset_trim = hsp->subject.offset + MIN(s_length, trim_size);
    // H->s_end_trim = hsp->subject.end - MIN(s_length, trim_size);
    // ```
    // The trimmed offsets, with the same expressions.
    let mut score = Vec::with_capacity(n);
    let mut qo = Vec::with_capacity(n);
    let mut so = Vec::with_capacity(n);
    let mut qe = Vec::with_capacity(n);
    let mut se = Vec::with_capacity(n);
    let mut sl = Vec::with_capacity(n);
    let mut lk = Vec::with_capacity(n);
    for hit in group_hits.iter() {
        let q_off = hit.q_aa_start as i32;
        let q_end = hit.q_aa_end as i32;
        let s_off = hit.s_aa_start as i32;
        let s_end = hit.s_aa_end as i32;
        let qt = TRIM_SIZE.min((q_end - q_off) / 4);
        let st = TRIM_SIZE.min((s_end - s_off) / 4);
        score.push(hit.raw_score);
        qo.push(q_off + qt);
        so.push(s_off + st);
        qe.push(q_end - qt);
        se.push(s_end - st);
        sl.push((hit.raw_score as f64) * query_contexts[hit.ctx_idx].karlin_params.lambda);
        lk.push(log_k_by_ctx[hit.ctx_idx]);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:750-757,764-766
    // ```c
    // double new_xsum =
    //   H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
    //   kbp[H->hsp->context]->logK;
    // Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
    // ...
    // H->hsp_link.sum[index] = new_sum;
    // H->hsp_link.num[index] = H_hsp_num+1;
    // H->hsp_link.link[index] = H_hsp_link;
    // ...
    // H->hsp_link.xsum[index] = new_xsum;
    // if(H_hsp_link)
    //    ((LinkHSPStruct*)H_hsp_link)->linked_to++;
    // ```
    // The pass values of LinkHSPStruct (hsp_link.sum, num, xsum, link) per method, initialised to
    // the values of an HSP without a predecessor; `indeg` is linked_to.
    let mut sum: [Vec<i64>; 2] = [
        (0..n)
            .map(|i| i64::from(score[i]) - i64::from(c[0]))
            .collect(),
        (0..n)
            .map(|i| i64::from(score[i]) - i64::from(c[1]))
            .collect(),
    ];
    let base_xsum: Vec<f64> = (0..n).map(|i| sl[i] - lk[i]).collect();
    let mut xsum: [Vec<f64>; 2] = [base_xsum.clone(), base_xsum];
    let mut num: [Vec<i16>; 2] = [vec![1; n], vec![1; n]];
    let mut link: [Vec<u32>; 2] = [vec![NONE; n], vec![NONE; n]];
    let mut indeg: Vec<i32> = vec![0; n];
    let mut alive: Vec<bool> = vec![true; n];
    let mut next_active: Vec<u32> = (0..n)
        .map(|i| if i + 1 < n { (i + 1) as u32 } else { NONE })
        .collect();
    let mut prev_active: Vec<u32> = (0..n)
        .map(|i| if i > 0 { (i - 1) as u32 } else { NONE })
        .collect();
    let mut active_head: u32 = 0;
    let mut children = [Children::new(n), Children::new(n)];
    // Per pass: queued for a method, sum changed in this pass, leaf to write.
    let mut queued: [Vec<bool>; 2] = [vec![false; n], vec![false; n]];
    let mut sum_changed: [Vec<bool>; 2] = [vec![false; n], vec![false; n]];
    let mut sum_changed_list: Vec<(u8, u32)> = Vec::new();
    let mut leaf_dirty: Vec<bool> = vec![false; n];
    let mut leaf_list: Vec<u32> = Vec::new();
    let mut seeds: [Vec<u32>; 2] = [Vec::new(), Vec::new()];
    let mut addback_hsps: Vec<u32> = Vec::new();
    let mut removed: Vec<u32> = Vec::new();
    let mut heap: BinaryHeap<Reverse<u32>> = BinaryHeap::new();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:702-745
    // ```c
    // if (H->hsp->score > cutoff[index]) {
    // ...
    //    for (H2_index=H_index-1; H2_index>1; H2_index=H2_index-1)
    // ...
    //       b1 = q_off_t <= H_query_etrim;
    //       b2 = s_off_t <= H_sub_etrim;
    // ...
    //       if (b1|b2|b5|b4) continue;
    // ...
    //       if (sum>H_hsp_sum)
    // ```
    // The grid of linking_index.rs for the index-0 scan: W x W cells over (q_off_t, s_off_t) of
    // the HSPs with score > cutoff[0], each listing its HSPs by increasing list index.
    let w = WINDOW_SIZE;
    let cell_of = |q: i32, s: i32| -> u64 {
        (((q.div_euclid(w)) as u32 as u64) << 32) | ((s.div_euclid(w)) as u32 as u64)
    };
    let mut grid_index: FxHashMap<u64, (u32, u32)> = FxHashMap::default();
    let mut grid_items: Vec<u32> = Vec::new();
    if !ignore_small_gaps {
        let mut cells: Vec<(u64, u32)> = (0..n)
            .filter(|&i| score[i] > c[0])
            .map(|i| (cell_of(qo[i], so[i]), i as u32))
            .collect();
        cells.sort_unstable();
        grid_items.reserve(cells.len());
        let mut start = 0usize;
        while start < cells.len() {
            let cell = cells[start].0;
            let mut end = start;
            while end < cells.len() && cells[end].0 == cell {
                grid_items.push(cells[end].1);
                end += 1;
            }
            grid_index.insert(cell, (start as u32, end as u32));
            start = end;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-861
    // ```c
    // b0 = sum <= H_hsp_sum;
    // ...
    // b1 = q_off_t <= H_query_etrim;
    // b2 = s_off_t <= H_sub_etrim;
    // ...
    // if (!(b0|b1|b2) )
    // ...
    //    H2 = H2_helper->ptr;
    // ```
    // The structures of the index-1 scan: the orders, ranks and Fenwick tree of the sweep of
    // linking_index.rs, over the HSPs with score > cutoff[1].
    let elig1: Vec<u32> = (0..n)
        .filter(|&i| score[i] > c[1])
        .map(|i| i as u32)
        .collect();
    let mut qorder: Vec<u32> = elig1.clone();
    qorder.sort_unstable_by(|&a, &b| qe[b as usize].cmp(&qe[a as usize]).then_with(|| a.cmp(&b)));
    let mut iorder: Vec<u32> = elig1.clone();
    iorder.sort_unstable_by(|&a, &b| qo[b as usize].cmp(&qo[a as usize]).then_with(|| a.cmp(&b)));
    let mut so_vals: Vec<i32> = elig1.iter().map(|&i| so[i as usize]).collect();
    so_vals.sort_unstable();
    so_vals.dedup();
    let msz = so_vals.len();
    let mut rr: Vec<u32> = vec![0; n];
    let mut kq: Vec<u32> = vec![0; n];
    for &i in elig1.iter() {
        let i = i as usize;
        let rank = so_vals.partition_point(|&v| v < so[i]);
        rr[i] = (msz - rank) as u32;
        let le = so_vals.partition_point(|&v| v <= se[i]);
        kq[i] = (msz - le) as u32;
    }
    let mut fen: Vec<K> = vec![K::default(); msz + 1];

    let mut tree = MaxTree::<K>::new(n);

    let mut remaining = n;
    let mut first_pass = true;
    let mut path_changed = true;
    let int4_max = i32::MAX as f64;

    while remaining > 0 {
        stats.rounds += 1;
        let mut best: [Option<usize>; 2] = [None, None];
        let mut use_current_max = false;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:603-652
        // ```c
        // use_current_max=0;
        // if (!first_pass){
        // ...
        //    if(path_changed==0){
        //       /* No path was changed, use these max sums. */
        //       use_current_max=1;
        // ...
        //       use_current_max=1;
        //       if(!ignore_small_gaps){
        //          for (H=best[0]; H!=NULL; H=H->hsp_link.link[0])
        //             if (H->linked_to==-1000) {use_current_max=0; break;}
        //       }
        //       if(use_current_max)
        //          for (H=best[1]; H!=NULL; H=H->hsp_link.link[1])
        //             if (H->linked_to==-1000) {use_current_max=0; break;}
        // ```
        // Same rule; a removed HSP is one that is not alive.
        if !first_pass {
            best = best_of_root(tree.root(), c, ignore_small_gaps);
            use_current_max = true;
            if path_changed {
                for (m, b) in best.iter().enumerate() {
                    if (m == 0 && ignore_small_gaps) || !use_current_max {
                        continue;
                    }
                    if let Some(bi) = *b {
                        let mut cur = bi as u32;
                        while cur != NONE {
                            if !alive[cur as usize] {
                                use_current_max = false;
                                break;
                            }
                            cur = link[m][cur as usize];
                        }
                    }
                }
            }
        }

        if !use_current_max {
            stats.passes += 1;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:685
            // ```c
            // H->linked_to = 0;
            // ```
            // NCBI counts linked_to again from the remaining HSPs: the links of the HSPs removed since
            // the previous pass are taken out, and the HSPs that selected them are computed again.
            if first_pass {
                for m in usize::from(ignore_small_gaps)..2 {
                    seeds[m].clear();
                    for i in 0..n {
                        if score[i] > c[m] {
                            queued[m][i] = true;
                            seeds[m].push(i as u32);
                        }
                    }
                }
            } else {
                stats.removed += removed.len() as u64;
                for &y in removed.iter() {
                    let y = y as usize;
                    for m in 0..2 {
                        let p = link[m][y];
                        if p != NONE && alive[p as usize] {
                            children[m].remove(p as usize, y);
                            indeg[p as usize] -= 1;
                        }
                        let mut ch = children[m].first[y];
                        while ch != NONE {
                            let cu = ch as usize;
                            let next = children[m].next[cu];
                            children[m].next[cu] = NONE;
                            children[m].prev[cu] = NONE;
                            if alive[cu] && !queued[m][cu] {
                                queued[m][cu] = true;
                                seeds[m].push(ch);
                            }
                            ch = next;
                        }
                        children[m].first[y] = NONE;
                    }
                }
                removed.clear();
                for &x in addback_hsps.iter() {
                    let x = x as usize;
                    if alive[x] && !leaf_dirty[x] {
                        leaf_dirty[x] = true;
                        leaf_list.push(x as u32);
                    }
                }
            }
            addback_hsps.clear();

            // Index 0 in list order (link_hsps.c:691-768), index 1 by the sweep in decreasing
            // q_end_trim order (link_hsps.c:771-896).
            for m in usize::from(ignore_small_gaps)..2 {
                let sweep = m == 1;
                if sweep {
                    seeds[1].clear();
                    for f in fen.iter_mut() {
                        *f = K::default();
                    }
                }
                let ilen = iorder.len();
                let mut p = 0usize;
                let mut qw = 0usize;
                let mut qi = 0usize;
                if !sweep {
                    heap.clear();
                    heap.extend(seeds[m].drain(..).map(Reverse));
                }
                loop {
                    // The next HSP to visit: from the heap in list order, or the next remaining HSP of
                    // the sweep (which also fills the Fenwick tree with the HSPs before it).
                    let x = if sweep {
                        let mut found = NONE;
                        while qi < qorder.len() {
                            let i = qorder[qi] as usize;
                            qi += 1;
                            if !alive[i] {
                                continue;
                            }
                            qorder[qw] = i as u32;
                            qw += 1;
                            let h_qe = qe[i];
                            while p < ilen {
                                let j = iorder[p] as usize;
                                if qo[j] <= h_qe {
                                    break;
                                }
                                p += 1;
                                if alive[j] {
                                    let k = K::new(sum[1][j], j);
                                    let mut f = rr[j] as usize;
                                    while f <= msz {
                                        if fen[f] >= k {
                                            break;
                                        }
                                        fen[f] = k;
                                        f += f & f.wrapping_neg();
                                    }
                                }
                            }
                            if !first_pass {
                                stats.swept += 1;
                            }
                            if queued[1][i] {
                                found = i as u32;
                                break;
                            }
                        }
                        if found == NONE {
                            break;
                        }
                        found as usize
                    } else {
                        match heap.pop() {
                            Some(Reverse(x)) => x as usize,
                            None => break,
                        }
                    };
                    queued[m][x] = false;
                    if !alive[x] || score[x] <= c[m] {
                        continue;
                    }
                    if m == 0 {
                        stats.visited0 += 1;
                    } else {
                        stats.visited1 += 1;
                    }
                    let old = link[m][x];
                    let keep = !first_pass
                        && old != NONE
                        && alive[old as usize]
                        && !sum_changed[m][old as usize];
                    let (h_sum, h_num, h_xsum, h_link) = if keep {
                        if m == 0 {
                            stats.kept0 += 1;
                        } else {
                            stats.kept1 += 1;
                        }
                        let o = old as usize;
                        (sum[m][o], num[m][o], xsum[m][o], old)
                    } else {
                        let mut bestk = K::default();
                        if m == 0 {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:702-745
                            // ```c
                            // if (H->hsp->score > cutoff[index]) {
                            // ...
                            //       b1 = q_off_t <= H_query_etrim;
                            //       b2 = s_off_t <= H_sub_etrim;
                            // ...
                            //       if (b1|b2|b5|b4) continue;
                            // ...
                            //       if (sum>H_hsp_sum)
                            // ```
                            // The grid walk of linking_index.rs: the remaining HSP before H in the
                            // window with the largest positive (sum, list index).
                            stats.searched0 += 1;
                            let h_qe = qe[x];
                            let h_se = se[x];
                            let h_qe_gap = h_qe + w;
                            let h_se_gap = h_se + w;
                            for cq in (h_qe + 1).div_euclid(w)..=h_qe_gap.div_euclid(w) {
                                for cs in (h_se + 1).div_euclid(w)..=h_se_gap.div_euclid(w) {
                                    let cell = ((cq as u32 as u64) << 32) | (cs as u32 as u64);
                                    if let Some(&(a, b)) = grid_index.get(&cell) {
                                        for &j in &grid_items[a as usize..b as usize] {
                                            let j = j as usize;
                                            if j >= x {
                                                break;
                                            }
                                            if !alive[j] {
                                                continue;
                                            }
                                            let (qj, sj) = (qo[j], so[j]);
                                            if qj <= h_qe
                                                || sj <= h_se
                                                || qj > h_qe_gap
                                                || sj > h_se_gap
                                            {
                                                continue;
                                            }
                                            let s = sum[0][j];
                                            if s > 0 {
                                                bestk = bestk.max(K::new(s, j));
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-861
                            // ```c
                            // for (H2_index=H_index-1; H2_index>1;)
                            // ...
                            //    b1 = q_off_t <= H_query_etrim;
                            //    b2 = s_off_t <= H_sub_etrim;
                            // ...
                            //    if (!(b0|b1|b2) )
                            // ...
                            //       H2 = H2_helper->ptr;
                            // ```
                            // The prefix maximum of the sweep; NCBI's scan when it returns an HSP
                            // after H (linking_index.rs, "Index 1").
                            stats.searched1 += 1;
                            let mut f = kq[x] as usize;
                            while f > 0 {
                                bestk = bestk.max(fen[f]);
                                f &= f - 1;
                            }
                            if bestk != K::default() && bestk.idx() > x {
                                stats.fallbacks += 1;
                                bestk = K::default();
                                let mut j = prev_active[x];
                                while j != NONE {
                                    let jj = j as usize;
                                    if score[jj] > c[1] && qo[jj] > qe[x] && so[jj] > se[x] {
                                        bestk = bestk.max(K::new(sum[1][jj], jj));
                                    }
                                    j = prev_active[jj];
                                }
                            }
                        }
                        if bestk != K::default() && bestk.idx() as u32 == old {
                            if m == 0 {
                                stats.same0 += 1;
                            } else {
                                stats.same1 += 1;
                            }
                        }
                        if bestk != K::default() {
                            let j = bestk.idx();
                            (sum[m][j], num[m][j], xsum[m][j], j as u32)
                        } else {
                            (0, 0, 0.0, NONE)
                        }
                    };

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:750-757,764-766
                    // ```c
                    // double new_xsum =
                    //   H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
                    //   kbp[H->hsp->context]->logK;
                    // Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
                    // ...
                    // H->hsp_link.sum[index] = new_sum;
                    // H->hsp_link.num[index] = H_hsp_num+1;
                    // H->hsp_link.link[index] = H_hsp_link;
                    // ...
                    // H->hsp_link.xsum[index] = new_xsum;
                    // if(H_hsp_link)
                    //    ((LinkHSPStruct*)H_hsp_link)->linked_to++;
                    // ```
                    // Same statements, the sum in i64 (module documentation, "Sums"); the counter and the
                    // reverse links follow a changed link, and the HSPs that selected H are computed
                    // again when H's values change.
                    let new_sum = h_sum + (i64::from(score[x]) - i64::from(c[m]));
                    stats.max_pass_sum = stats.max_pass_sum.max(new_sum);
                    stats.max_num = stats.max_num.max(i64::from(h_num) + 1);
                    let new_num = h_num.wrapping_add(1);
                    let new_xsum = h_xsum + sl[x] - lk[x];
                    if h_link != old {
                        if old != NONE && alive[old as usize] {
                            children[m].remove(old as usize, x);
                            indeg[old as usize] -= 1;
                        }
                        if h_link != NONE {
                            children[m].insert(h_link as usize, x);
                            indeg[h_link as usize] += 1;
                        }
                    }
                    let sum_moved = new_sum != sum[m][x];
                    let changed = sum_moved
                        || new_num != num[m][x]
                        || new_xsum.to_bits() != xsum[m][x].to_bits()
                        || h_link != old;
                    if changed {
                        sum[m][x] = new_sum;
                        num[m][x] = new_num;
                        xsum[m][x] = new_xsum;
                        link[m][x] = h_link;
                        if sum_moved {
                            sum_changed[m][x] = true;
                            sum_changed_list.push((m as u8, x as u32));
                            if !leaf_dirty[x] {
                                leaf_dirty[x] = true;
                                leaf_list.push(x as u32);
                            }
                        }
                        let mut ch = children[m].first[x];
                        while ch != NONE {
                            let cu = ch as usize;
                            if alive[cu] && !queued[m][cu] {
                                queued[m][cu] = true;
                                if !sweep {
                                    heap.push(Reverse(ch));
                                }
                            }
                            ch = children[m].next[cu];
                        }
                    }
                }
                if sweep {
                    qorder.truncate(qw);
                    iorder.retain(|&j| alive[j as usize]);
                }
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:759-763,887-891
            // ```c
            // if (new_sum >= maxscore)
            // {
            //    maxscore=new_sum;
            //    best[index]=H;
            // }
            // ```
            // best[] after the pass from the key tree; the leaves of the HSPs whose pass sum changed or
            // that received an add-back are written.
            if first_pass {
                for i in 0..n {
                    tree.set(i, [K::new(sum[0][i], i), K::new(sum[1][i], i)]);
                }
            } else {
                for &x in leaf_list.iter() {
                    let x = x as usize;
                    if alive[x] {
                        tree.set(x, [K::new(sum[0][x], x), K::new(sum[1][x], x)]);
                    }
                }
            }
            for &x in leaf_list.iter() {
                leaf_dirty[x as usize] = false;
            }
            leaf_list.clear();
            for &(m, x) in sum_changed_list.iter() {
                sum_changed[m as usize][x as usize] = false;
            }
            sum_changed_list.clear();

            if options.verify {
                let (expected, expected_indeg) = verify_pass(
                    &alive,
                    &score,
                    &qo,
                    &so,
                    &qe,
                    &se,
                    &sl,
                    &lk,
                    c,
                    ignore_small_gaps,
                );
                let mut cur = active_head;
                while cur != NONE {
                    let i = cur as usize;
                    for m in usize::from(ignore_small_gaps)..2 {
                        let actual = (sum[m][i], num[m][i], xsum[m][i].to_bits(), link[m][i]);
                        assert_eq!(
                            expected[m][i], actual,
                            "LOSAT_LINK_FAST_VERIFY: pass {}, method {m}, HSP {i} of {n}",
                            stats.passes
                        );
                    }
                    assert_eq!(
                        expected_indeg[i], indeg[i],
                        "LOSAT_LINK_FAST_VERIFY: pass {}, linked_to of HSP {i} of {n}",
                        stats.passes
                    );
                    cur = next_active[i];
                }
                stats.verified_passes += 1;
            }
            best = best_of_root(tree.root(), c, ignore_small_gaps);
            path_changed = false;
            first_pass = false;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:901-937
        // ```c
        // best[0]->hsp_link.sum[0] +=
        //    (best[0]->hsp_link.num[0])*cutoff[0];
        // ...
        // prob[0] = BLAST_SmallGapSumE(window_size,
        // ...
        // if( best[0]->hsp_link.num[0] > 1 ) {
        //   if( gap_prob == 0 || (prob[0] /= gap_prob) > INT4_MAX ) {
        // ...
        // prob[1] = BLAST_LargeGapSumE(best[1]->hsp_link.num[1],
        // ...
        // if( best[1]->hsp_link.num[1] > 1 ) {
        //   if( 1 - gap_prob == 0 || (prob[1] /= 1 - gap_prob) > INT4_MAX ) {
        // ...
        // ordering_method =
        //    prob[0]<=prob[1] ? eLinkSmallGaps : eLinkLargeGaps;
        // ```
        // Same statements. The add-back changes the stored sum, which is the leaf of the key tree.
        let mut prob = [f64::MAX, f64::MAX];
        // The value is clamped to the range of the key (module documentation, "Sums").
        let addback = |tree: &mut MaxTree<K>, m: usize, bi: usize, stats: &mut IncrStats| {
            let mut v = tree.leaf(bi);
            let added = v[m]
                .value()
                .saturating_add(i64::from(num[m][bi]) * i64::from(c[m]));
            let clamped = added.clamp(K::MIN_VALUE, K::MAX_VALUE);
            stats.clamped += u64::from(clamped != added);
            v[m] = K::new(clamped, bi);
            tree.set(bi, v);
        };
        if !ignore_small_gaps {
            if let Some(bi) = best[0] {
                addback(&mut tree, 0, bi, &mut stats);
                addback_hsps.push(bi as u32);
                let k = num[0][bi] as usize;
                let divisor = gap_decay_divisor(gap_decay_rate, k);
                prob[0] = small_gap_sum_e(
                    WINDOW_SIZE,
                    k as i16,
                    xsum[0][bi],
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space,
                    divisor,
                );
                if k > 1 {
                    if gap_prob == 0.0 || (prob[0] / gap_prob) > int4_max {
                        prob[0] = int4_max;
                    } else {
                        prob[0] /= gap_prob;
                    }
                }
            }
            if let Some(bi) = best[1] {
                let k = num[1][bi] as usize;
                let divisor = gap_decay_divisor(gap_decay_rate, k);
                prob[1] = ncbi_large_gap_sum_e(
                    k as i16,
                    xsum[1][bi],
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space,
                    divisor,
                );
                if k > 1 {
                    let denom = 1.0 - gap_prob;
                    if denom == 0.0 || (prob[1] / denom) > int4_max {
                        prob[1] = int4_max;
                    } else {
                        prob[1] /= denom;
                    }
                }
            }
        } else if let Some(bi) = best[1] {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:939-953
            // ```c
            // best[1]->hsp_link.sum[1] +=
            //    (best[1]->hsp_link.num[1])*cutoff[1];
            // ...
            // prob[1] = BLAST_LargeGapSumE(
            //              best[1]->hsp_link.num[1],
            // ...
            // ordering_method = eLinkLargeGaps;
            // ```
            // Only large gaps are considered when cutoff[0] == 0.
            addback(&mut tree, 1, bi, &mut stats);
            addback_hsps.push(bi as u32);
            let k = num[1][bi] as usize;
            let divisor = gap_decay_divisor(gap_decay_rate, k);
            prob[1] = ncbi_large_gap_sum_e(
                k as i16,
                xsum[1][bi],
                eff_query_len as i32,
                eff_subject_len as i32,
                eff_search_space,
                divisor,
            );
            if k > 1 {
                let denom = 1.0 - gap_prob;
                if denom == 0.0 || (prob[1] / denom) > int4_max {
                    prob[1] = int4_max;
                } else {
                    prob[1] /= denom;
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:936-937
        // ```c
        // ordering_method =
        //    prob[0]<=prob[1] ? eLinkSmallGaps : eLinkLargeGaps;
        // ```
        // Same choice of ordering method.
        let ordering = if !ignore_small_gaps && prob[0] <= prob[1] {
            0
        } else {
            1
        };
        let best_i = if let Some(i) = best[ordering] {
            i
        } else if let Some(i) = best[1 - ordering] {
            i
        } else {
            break;
        };
        let evalue = prob[ordering];

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:955-980
        // ```c
        // best[ordering_method]->start_of_chain = TRUE;
        // ...
        // if (best[ordering_method]->linked_to>0) path_changed=1;
        // for (H=best[ordering_method]; H!=NULL;
        //      H=H->hsp_link.link[ordering_method])
        // ...
        //    if (H->linked_to>1) path_changed=1;
        //    H->linked_to=-1000;
        //    H->hsp_link.changed=1;
        // ...
        //    H->linked_set = linked_set;
        //    H->ordering_method = ordering_method;
        //    H->hsp->evalue = prob[ordering_method];
        //    if (H->next)
        //       (H->next)->prev=H->prev;
        //    if (H->prev)
        //       (H->prev)->next=H->next;
        //    number_of_hsps--;
        // ```
        // Same removal of the chosen chain and the same fields set on it, written to the group as
        // the chain is removed; `indeg` is linked_to as of the last pass.
        if indeg[best_i] > 0 {
            path_changed = true;
        }
        let linked_set = link[ordering][best_i] != NONE;
        group_hits[best_i].hsp_link_num = num[ordering][best_i];
        let mut cur = best_i;
        loop {
            if indeg[cur] > 1 {
                path_changed = true;
            }
            alive[cur] = false;
            removed.push(cur as u32);
            let prev_idx = prev_active[cur];
            let next_idx = next_active[cur];
            if prev_idx != NONE {
                next_active[prev_idx as usize] = next_idx;
            } else {
                active_head = next_idx;
            }
            if next_idx != NONE {
                prev_active[next_idx as usize] = prev_idx;
            }
            tree.set(cur, [K::default(); 2]);
            let next = link[ordering][cur];
            let chain_next_link_id = (next != NONE).then(|| group_hits[next as usize].link_id);
            let hit = &mut group_hits[cur];
            hit.ordering_method = ordering as u8;
            hit.e_value = evalue;
            hit.linked_set = linked_set;
            hit.chain_next_link_id = chain_next_link_id;
            hit.start_of_chain = cur == best_i;
            remaining -= 1;
            if next == NONE {
                break;
            }
            cur = next as usize;
        }
    }

    (group_hits, stats)
}

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
//! start value and `next_larger` select the same HSP). While no pass value
//! leaves Int4 (`linking_index.rs`, "Int4 range") these are the exact chain
//! maxima over the remaining HSPs.
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
//! The candidates of H are the remaining HSPs before H in the list with
//! score > cutoff[1], `q_off_t > q_end_trim(H)` and `s_off_t > s_end_trim(H)`
//! (link_hsps.c:827-861). `DomTree` (a 2-d tree over (q_off_t, s_off_t) with
//! the largest key and the range of list indices of every subtree) answers
//! "largest key among them" for one H. When many HSPs need index 1 in a pass,
//! the pass instead sweeps the remaining HSPs by decreasing `q_end_trim` with
//! the prefix-maximum Fenwick tree of `linking_index.rs`, computing again only
//! the HSPs that need it.
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

use super::linking_index::{
    best_of_root, int4_sum, key, key_idx, key_sum, stats_enabled, sum_bound, MaxTree, NONE,
    TRIM_SIZE, WINDOW_SIZE,
};
use super::params::LinkHspCutoffs;

// No NCBI counterpart: how index 1 is searched in a pass after the first. Every choice selects
// the same HSP; it does not change any value NCBI computes.
/// How a pass after the first searches index 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Index1Search {
    /// The tree when few HSPs need index 1, the sweep otherwise.
    Auto,
    Tree,
    Sweep,
}

// No NCBI counterpart: options of the check (verify), the Int4 check, and the index-1 search.
/// How the kernel runs. `from_env` is what a search uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct IncrOptions {
    /// Compare every pass with a pass computed from scratch.
    pub verify: bool,
    /// Stop at the first pass value outside Int4 (`linking_index.rs`, "Int4
    /// range"). Only tests turn it off.
    pub check_int4: bool,
    pub index1: Index1Search,
    /// `Auto` sweeps when (HSPs needing index 1 at the start of the pass) *
    /// this >= remaining HSPs with score > cutoff[1].
    pub sweep_factor: u32,
}

impl IncrOptions {
    pub(super) fn from_env() -> Self {
        static OPTIONS: OnceLock<IncrOptions> = OnceLock::new();
        *OPTIONS.get_or_init(|| IncrOptions {
            verify: std::env::var_os("LOSAT_LINK_FAST_VERIFY").is_some(),
            check_int4: true,
            index1: match std::env::var("LOSAT_LINK_INCR_INDEX1").as_deref() {
                Ok("tree") => Index1Search::Tree,
                Ok("sweep") => Index1Search::Sweep,
                _ => Index1Search::Auto,
            },
            sweep_factor: std::env::var("LOSAT_LINK_INCR_SWEEP")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(16),
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
    /// Nodes of `DomTree` read by the index-1 searches.
    pub tree_visits: u64,
    /// Passes after the first that swept index 1, and HSPs those sweeps read.
    pub sweeps: u64,
    pub swept: u64,
    /// Sweep searches answered by NCBI's scan (the tree returned a later HSP).
    pub fallbacks: u64,
    pub verified_passes: u64,
    /// Largest number of HSPs of a pass value (exact, before NCBI's Int2).
    pub max_num: i64,
    pub max_pass_sum: i64,
    pub max_addback_sum: i64,
    /// The first pass with a pass value outside Int4 (1 = the first pass).
    pub int4_overflow_pass: u64,
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
        self.tree_visits += o.tree_visits;
        self.sweeps += o.sweeps;
        self.swept += o.swept;
        self.fallbacks += o.fallbacks;
        self.verified_passes += o.verified_passes;
        self.max_num = self.max_num.max(o.max_num);
        self.max_pass_sum = self.max_pass_sum.max(o.max_pass_sum);
        self.max_addback_sum = self.max_addback_sum.max(o.max_addback_sum);
        self.int4_overflow_pass = self.int4_overflow_pass.max(o.int4_overflow_pass);
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
// The tree answers the index-1 scan for one H: the remaining HSP before H with q_off_t >
// H_query_etrim and s_off_t > H_sub_etrim and the largest (sum, list index), the HSP NCBI's
// downward scan ends with.
#[derive(Clone, Copy)]
struct DomNode {
    max: u64,
    own: u64,
    qo_lo: i32,
    qo_hi: i32,
    so_lo: i32,
    so_hi: i32,
    idx_lo: u32,
    idx_hi: u32,
    left: u32,
    right: u32,
    pt: u32,
    pt_qo: i32,
    pt_so: i32,
}

/// A 2-d tree over the HSPs with score > cutoff[1] at (q_off_t, s_off_t):
/// every node holds one HSP, the bounding box and list-index range of its
/// subtree, and the largest key in it (0 for removed HSPs).
struct DomTree {
    nodes: Vec<DomNode>,
    parent: Vec<u32>,
    node_of: Vec<u32>,
    root: u32,
}

impl DomTree {
    fn build(points: &[u32], qo: &[i32], so: &[i32], n: usize) -> Self {
        let mut t = DomTree {
            nodes: Vec::with_capacity(points.len()),
            parent: Vec::with_capacity(points.len()),
            node_of: vec![NONE; n],
            root: NONE,
        };
        let mut pts = points.to_vec();
        if !pts.is_empty() {
            t.root = t.build_rec(&mut pts, qo, so, 0, NONE);
        }
        t
    }

    fn build_rec(&mut self, pts: &mut [u32], qo: &[i32], so: &[i32], depth: u32, par: u32) -> u32 {
        let mid = pts.len() / 2;
        if depth % 2 == 0 {
            pts.select_nth_unstable_by_key(mid, |&p| (qo[p as usize], p));
        } else {
            pts.select_nth_unstable_by_key(mid, |&p| (so[p as usize], p));
        }
        let pt = pts[mid];
        let mut node = DomNode {
            max: 0,
            own: 0,
            qo_lo: i32::MAX,
            qo_hi: i32::MIN,
            so_lo: i32::MAX,
            so_hi: i32::MIN,
            idx_lo: u32::MAX,
            idx_hi: 0,
            left: NONE,
            right: NONE,
            pt,
            pt_qo: qo[pt as usize],
            pt_so: so[pt as usize],
        };
        for &p in pts.iter() {
            let (q, s) = (qo[p as usize], so[p as usize]);
            node.qo_lo = node.qo_lo.min(q);
            node.qo_hi = node.qo_hi.max(q);
            node.so_lo = node.so_lo.min(s);
            node.so_hi = node.so_hi.max(s);
            node.idx_lo = node.idx_lo.min(p);
            node.idx_hi = node.idx_hi.max(p);
        }
        let id = self.nodes.len() as u32;
        self.nodes.push(node);
        self.parent.push(par);
        self.node_of[pt as usize] = id;
        let (lo, rest) = pts.split_at_mut(mid);
        let hi = &mut rest[1..];
        if !lo.is_empty() {
            let l = self.build_rec(lo, qo, so, depth + 1, id);
            self.nodes[id as usize].left = l;
        }
        if !hi.is_empty() {
            let r = self.build_rec(hi, qo, so, depth + 1, id);
            self.nodes[id as usize].right = r;
        }
        id
    }

    #[inline]
    fn subtree_max(&self, v: usize) -> u64 {
        let nd = &self.nodes[v];
        let mut m = nd.own;
        if nd.left != NONE {
            m = m.max(self.nodes[nd.left as usize].max);
        }
        if nd.right != NONE {
            m = m.max(self.nodes[nd.right as usize].max);
        }
        m
    }

    /// Sets the key of every HSP in the tree (children follow their parent in
    /// `nodes`, so one backward walk fills the maxima).
    fn fill(&mut self, key_of: impl Fn(usize) -> u64) {
        for v in 0..self.nodes.len() {
            self.nodes[v].own = key_of(self.nodes[v].pt as usize);
        }
        for v in (0..self.nodes.len()).rev() {
            self.nodes[v].max = self.subtree_max(v);
        }
    }

    fn set(&mut self, hsp: usize, k: u64) {
        let mut v = self.node_of[hsp];
        if v == NONE {
            return;
        }
        self.nodes[v as usize].own = k;
        while v != NONE {
            let m = self.subtree_max(v as usize);
            if m == self.nodes[v as usize].max {
                break;
            }
            self.nodes[v as usize].max = m;
            v = self.parent[v as usize];
        }
    }

    /// The largest key among the HSPs with q_off_t > `a`, s_off_t > `b` and a
    /// list index below `before` (0: none), and the nodes read.
    fn query(&self, a: i32, b: i32, before: u32, stack: &mut Vec<u32>) -> (u64, u64) {
        let mut best = 0u64;
        let mut visits = 0u64;
        if self.root == NONE {
            return (0, 0);
        }
        stack.clear();
        stack.push(self.root);
        while let Some(v) = stack.pop() {
            visits += 1;
            let nd = &self.nodes[v as usize];
            if nd.max <= best || nd.qo_hi <= a || nd.so_hi <= b || nd.idx_lo >= before {
                continue;
            }
            if nd.qo_lo > a && nd.so_lo > b && nd.idx_hi < before {
                best = nd.max;
                continue;
            }
            if nd.own > best && nd.pt_qo > a && nd.pt_so > b && nd.pt < before {
                best = nd.own;
            }
            let (l, r) = (nd.left, nd.right);
            if l != NONE && r != NONE {
                if self.nodes[l as usize].max >= self.nodes[r as usize].max {
                    stack.push(r);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(r);
                }
            } else if l != NONE {
                stack.push(l);
            } else if r != NONE {
                stack.push(r);
            }
        }
        (best, visits)
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
/// The kernel `linking.rs` calls under `LOSAT_LINK_FAST=2`. `Err` returns the
/// group unchanged when a pass value left the Int4 range; `linking.rs` then
/// links it with the literal port.
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
) -> Result<Vec<UngappedHit>, Vec<UngappedHit>> {
    let n = group_hits.len();
    let stats_on = stats_enabled() && n > 0;
    let bound = if stats_on {
        sum_bound(&group_hits, cutoffs)
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
            "[LINK_INCR_STATS] n={} us={} bound={} max_pass_sum={} max_addback_sum={} max_num={} int4_overflow_pass={} rounds={} passes={} removed={} visited0={} kept0={} searched0={} visited1={} kept1={} searched1={} tree_visits={} sweeps={} swept={} fallbacks={} verified_passes={}",
            n,
            started.elapsed().as_micros(),
            bound,
            stats.max_pass_sum,
            stats.max_addback_sum,
            stats.max_num,
            stats.int4_overflow_pass,
            stats.rounds,
            stats.passes,
            stats.removed,
            stats.visited0,
            stats.kept0,
            stats.searched0,
            stats.visited1,
            stats.kept1,
            stats.searched1,
            stats.tree_visits,
            stats.sweeps,
            stats.swept,
            stats.fallbacks,
            stats.verified_passes
        );
    }
    if stats.int4_overflow_pass != 0 {
        Err(group_hits)
    } else {
        Ok(group_hits)
    }
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
) -> ([Vec<(i32, i16, u64, u32)>; 2], Vec<i32>) {
    let n = alive.len();
    let mut out: [Vec<(i32, i16, u64, u32)>; 2] = [
        (0..n)
            .map(|i| {
                (
                    score[i].wrapping_sub(c[0]),
                    1,
                    (sl[i] - lk[i]).to_bits(),
                    NONE,
                )
            })
            .collect(),
        (0..n)
            .map(|i| {
                (
                    score[i].wrapping_sub(c[1]),
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
            let (mut h_sum, mut h_num, mut h_xsum, mut h_link) = (0i32, 0i16, 0.0f64, NONE);
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
            let (new_sum, _, _) = int4_sum(h_sum, score[i], c[m]);
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
// The replacement for the per-group body of s_BlastEvenGapLinkHSPs: same rounds, same selection
// rule, same removal; a pass after the first computes again only the HSPs described in the
// module documentation.
#[allow(clippy::too_many_arguments)]
pub(super) fn link_hsp_group_incr_with(
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
    let mut sum: [Vec<i32>; 2] = [
        (0..n).map(|i| score[i].wrapping_sub(c[0])).collect(),
        (0..n).map(|i| score[i].wrapping_sub(c[1])).collect(),
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:955,970-974
    // ```c
    // best[ordering_method]->start_of_chain = TRUE;
    // ...
    //    H->linked_set = linked_set;
    //    H->ordering_method = ordering_method;
    //    H->hsp->evalue = prob[ordering_method];
    // ```
    // The fields NCBI sets on a removed chain, collected and written to the group at the end.
    let mut out_ordering: Vec<u8> = vec![0; n];
    let mut out_evalue: Vec<f64> = vec![0.0; n];
    let mut out_linked_set: Vec<bool> = vec![false; n];
    let mut out_head_num: Vec<Option<i16>> = vec![None; n];
    let mut out_next: Vec<u32> = vec![NONE; n];

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
    // linking_index.rs, and the 2-d tree, over the HSPs with score > cutoff[1].
    let elig1: Vec<u32> = (0..n)
        .filter(|&i| score[i] > c[1])
        .map(|i| i as u32)
        .collect();
    let mut alive1 = elig1.len();
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
    let mut fen: Vec<u64> = vec![0; msz + 1];
    let mut dom = DomTree::build(&elig1, &qo, &so, n);
    let mut stack: Vec<u32> = Vec::new();

    let mut tree = MaxTree::new(n);

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
            let mut pass_int4_overflow = false;

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

            // Index 0 in list order (link_hsps.c:691-768), and index 1 by the tree in list order or
            // by the sweep in decreasing q_end_trim order (link_hsps.c:771-896).
            for m in usize::from(ignore_small_gaps)..2 {
                let sweep = m == 1
                    && (first_pass
                        || match options.index1 {
                            Index1Search::Tree => false,
                            Index1Search::Sweep => true,
                            Index1Search::Auto => {
                                seeds[1].len() as u64 * u64::from(options.sweep_factor)
                                    >= alive1 as u64
                            }
                        });
                if sweep {
                    if !first_pass {
                        stats.sweeps += 1;
                    }
                    seeds[1].clear();
                    for f in fen.iter_mut() {
                        *f = 0;
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
                                    let k = key(sum[1][j], j);
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
                        let mut bestk = 0u64;
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
                                                bestk = bestk.max(key(s, j));
                                            }
                                        }
                                    }
                                }
                            }
                        } else if sweep {
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
                            if bestk != 0 && key_idx(bestk) > x {
                                stats.fallbacks += 1;
                                bestk = 0;
                                let mut j = prev_active[x];
                                while j != NONE {
                                    let jj = j as usize;
                                    if score[jj] > c[1] && qo[jj] > qe[x] && so[jj] > se[x] {
                                        bestk = bestk.max(key(sum[1][jj], jj));
                                    }
                                    j = prev_active[jj];
                                }
                            }
                        } else {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-861
                            // ```c
                            // b0 = sum <= H_hsp_sum;
                            // ...
                            // b1 = q_off_t <= H_query_etrim;
                            // b2 = s_off_t <= H_sub_etrim;
                            // ...
                            // if (!(b0|b1|b2) )
                            // ```
                            // The 2-d tree: the remaining HSPs before H with q_off_t > q_end_trim and
                            // s_off_t > s_end_trim, largest (sum, list index).
                            stats.searched1 += 1;
                            let (k, visits) = dom.query(qe[x], se[x], x as u32, &mut stack);
                            stats.tree_visits += visits;
                            bestk = k;
                        }
                        if bestk != 0 {
                            let j = key_idx(bestk);
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
                    // Same statements; the counter and the reverse links follow a changed link, and the
                    // HSPs that selected H are computed again when H's values change.
                    let (new_sum, exact, outside) = int4_sum(h_sum, score[x], c[m]);
                    stats.max_pass_sum = stats.max_pass_sum.max(exact);
                    stats.max_num = stats.max_num.max(i64::from(h_num) + 1);
                    pass_int4_overflow |= outside;
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
                            if m == 1 && !first_pass {
                                dom.set(x, key(new_sum, x));
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
                    tree.set(i, [key(sum[0][i], i), key(sum[1][i], i)]);
                }
                dom.fill(|i| key(sum[1][i], i));
            } else {
                for &x in leaf_list.iter() {
                    let x = x as usize;
                    if alive[x] {
                        tree.set(x, [key(sum[0][x], x), key(sum[1][x], x)]);
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

            if pass_int4_overflow && stats.int4_overflow_pass == 0 {
                stats.int4_overflow_pass = stats.passes;
            }
            if options.check_int4 && pass_int4_overflow {
                return (group_hits, stats);
            }
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
        let addback = |tree: &mut MaxTree, m: usize, bi: usize, stats: &mut IncrStats| {
            let leaf = tree.leaf(bi);
            let added = i64::from(key_sum(leaf[m])) + i64::from(num[m][bi]) * i64::from(c[m]);
            stats.max_addback_sum = stats.max_addback_sum.max(added);
            let mut v = leaf;
            v[m] = key(added as i32, bi);
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
        //    if (H->next)
        //       (H->next)->prev=H->prev;
        //    if (H->prev)
        //       (H->prev)->next=H->next;
        //    number_of_hsps--;
        // ```
        // Same removal of the chosen chain; `indeg` is linked_to as of the last pass.
        if indeg[best_i] > 0 {
            path_changed = true;
        }
        let linked_set = link[ordering][best_i] != NONE;
        out_head_num[best_i] = Some(num[ordering][best_i]);
        let mut cur = best_i;
        loop {
            if indeg[cur] > 1 {
                path_changed = true;
            }
            alive[cur] = false;
            removed.push(cur as u32);
            if score[cur] > c[1] {
                alive1 -= 1;
                dom.set(cur, 0);
            }
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
            tree.set(cur, [0, 0]);
            out_ordering[cur] = ordering as u8;
            out_evalue[cur] = evalue;
            out_linked_set[cur] = linked_set;
            let next = link[ordering][cur];
            out_next[cur] = next;
            remaining -= 1;
            if next == NONE {
                break;
            }
            cur = next as usize;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:955,970-974
    // ```c
    // best[ordering_method]->start_of_chain = TRUE;
    // ...
    //    H->linked_set = linked_set;
    //    H->ordering_method = ordering_method;
    //    H->hsp->evalue = prob[ordering_method];
    // ```
    // The collected fields of every removed HSP, as linking_index.rs writes them.
    for i in 0..n {
        if alive[i] {
            continue;
        }
        let next = out_next[i];
        let chain_next_link_id = (next != NONE).then(|| group_hits[next as usize].link_id);
        let hit = &mut group_hits[i];
        hit.ordering_method = out_ordering[i];
        hit.e_value = out_evalue[i];
        hit.linked_set = out_linked_set[i];
        hit.chain_next_link_id = chain_next_link_id;
        hit.start_of_chain = out_head_num[i].is_some();
        if let Some(k) = out_head_num[i] {
            hit.hsp_link_num = k;
        }
    }

    (group_hits, stats)
}

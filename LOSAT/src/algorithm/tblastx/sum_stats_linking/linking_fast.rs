//! Incremental, result-identical implementation of NCBI's even-gap HSP linking
//! (`s_BlastEvenGapLinkHSPs`, ncbi-blast/c++/src/algo/blast/core/link_hsps.c:414-1091).
//!
//! NCBI's procedure, per frame group, is: repeat { compute the best-predecessor DP for
//! every remaining HSP (two ordering methods), pick the best chain per method, choose the
//! method by its sum E-value, mark and remove the chain } until no HSP remains. NCBI
//! recomputes the DP with an O(n^2) scan (`for (H2_index=H_index-1; H2_index>1;)`,
//! link_hsps.c:710 and :827) and repeats it after chain removals, which is cubic in the
//! worst case.
//!
//! This module computes exactly the same sequence of (best chain, ordering method,
//! E-value) decisions with a different evaluation strategy:
//!
//! * The DP value of an HSP is `argmax (sum, list index)` over the valid earlier HSPs
//!   (link_hsps.c:736-745 and :838-860, where ties go to the first entry met while
//!   scanning from `H_index-1` downwards, i.e. the largest list index). That query is
//!   answered with a 2-D range-maximum tree over (list index, `s_off_trim`) plus a
//!   bounded scan of the few HSPs whose `q_off_trim` can differ from `q_off` by at most
//!   `trim_size` (link_hsps.c:728-734).
//! * After a chain is removed only the HSPs whose best chain passed through a removed HSP
//!   can change (DP values are monotone non-increasing under removals); they are
//!   recomputed in list order, exactly as a full NCBI pass would recompute them
//!   (link_hsps.c:781-795 reuses an unchanged predecessor for the same reason).
//! * NCBI's bookkeeping that decides *when* it recomputes (`path_changed`,
//!   `use_current_max`, the `linked_to` counters, link_hsps.c:603-652 and :958-981) and
//!   the persistent `best[0]->hsp_link.sum[0] += num*cutoff[0]` add-back
//!   (link_hsps.c:907-908) are reproduced on a snapshot ("stored") copy of the values,
//!   because they influence which HSP NCBI selects as `best[]` in the rounds between
//!   two recomputations.
//!
//! The chain selection, the ordering-method choice, the E-values (same floating-point
//! expressions in the same order), and the output flags are therefore identical to
//! `link_hsp_group_ncbi`, which remains available (`LOSAT_LINKING_LEGACY=1`) as the
//! reference implementation.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::algorithm::tblastx::chaining::UngappedHit;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::stats::sum_statistics::{
    defaults::{GAP_SIZE, OVERLAP_SIZE},
    gap_decay_divisor, ncbi_large_gap_sum_e, small_gap_sum_e,
};

use super::params::LinkHspCutoffs;

const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
const TRIM_SIZE: i32 = (OVERLAP_SIZE + 1) / 2;
const NONE: u32 = u32::MAX;

/// Lexicographic (sum, list index) candidate value; `NEG` never wins.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Cand {
    sum: i32,
    idx: u32,
}

impl Cand {
    const NEG: Cand = Cand {
        sum: i32::MIN,
        idx: 0,
    };
}

// ---------------------------------------------------------------------------
// Static balanced k-d tree over (list index, s_off_trim) with subtree maxima
// per ordering method. O(n) memory; point values change (activate / decrease /
// remove) and rectangle maximum queries prune on the subtree maxima.
// ---------------------------------------------------------------------------

struct KdNode {
    point: u32,
    left: u32,
    right: u32,
    // bounding box of the subtree (inclusive)
    x_lo: u32,
    x_hi: u32,
    y_lo: i32,
    y_hi: i32,
    max: [Cand; 2],
}

struct KdTree {
    nodes: Vec<KdNode>,
    parent: Vec<u32>,
    node_of_point: Vec<u32>,
    root: u32,
}

impl KdTree {
    fn build(ys: &[i32]) -> Self {
        let n = ys.len();
        let mut nodes: Vec<KdNode> = Vec::with_capacity(n);
        let mut parent = vec![NONE; n];
        let mut node_of_point = vec![NONE; n];
        let mut pts: Vec<u32> = (0..n as u32).collect();
        let root = if n == 0 {
            NONE
        } else {
            Self::build_rec(
                &mut nodes,
                &mut parent,
                &mut node_of_point,
                &mut pts,
                ys,
                0,
                NONE,
            )
        };
        KdTree {
            nodes,
            parent,
            node_of_point,
            root,
        }
    }

    fn build_rec(
        nodes: &mut Vec<KdNode>,
        parent: &mut [u32],
        node_of_point: &mut [u32],
        pts: &mut [u32],
        ys: &[i32],
        depth: u32,
        par: u32,
    ) -> u32 {
        let len = pts.len();
        let mid = len / 2;
        if depth % 2 == 0 {
            pts.select_nth_unstable(mid);
        } else {
            pts.select_nth_unstable_by_key(mid, |&p| (ys[p as usize], p));
        }
        let point = pts[mid];
        let id = nodes.len() as u32;
        let mut x_lo = point;
        let mut x_hi = point;
        let mut y_lo = ys[point as usize];
        let mut y_hi = ys[point as usize];
        for &p in pts.iter() {
            x_lo = x_lo.min(p);
            x_hi = x_hi.max(p);
            y_lo = y_lo.min(ys[p as usize]);
            y_hi = y_hi.max(ys[p as usize]);
        }
        nodes.push(KdNode {
            point,
            left: NONE,
            right: NONE,
            x_lo,
            x_hi,
            y_lo,
            y_hi,
            max: [Cand::NEG; 2],
        });
        parent[point as usize] = par;
        node_of_point[point as usize] = id;
        let (lo, rest) = pts.split_at_mut(mid);
        let hi = &mut rest[1..];
        if !lo.is_empty() {
            let l = Self::build_rec(nodes, parent, node_of_point, lo, ys, depth + 1, id);
            nodes[id as usize].left = l;
        }
        if !hi.is_empty() {
            let r = Self::build_rec(nodes, parent, node_of_point, hi, ys, depth + 1, id);
            nodes[id as usize].right = r;
        }
        id
    }

    /// Set the value of `point` for method `m` and refresh the maxima on its root path.
    #[inline]
    fn update(&mut self, point: u32, m: usize, value: Cand, own: &mut [Cand]) {
        own[point as usize] = value;
        let mut v = self.node_of_point[point as usize];
        while v != NONE {
            let node = &self.nodes[v as usize];
            let mut mx = own[node.point as usize];
            if node.left != NONE {
                mx = mx.max(self.nodes[node.left as usize].max[m]);
            }
            if node.right != NONE {
                mx = mx.max(self.nodes[node.right as usize].max[m]);
            }
            let node = &mut self.nodes[v as usize];
            if node.max[m] == mx {
                break;
            }
            node.max[m] = mx;
            v = self.parent[node.point as usize];
        }
    }

    /// Maximum value for method `m` over points with list index in `[x_lo, x_hi]`
    /// and y in `[y_lo, y_hi]` (all inclusive), starting from `best`.
    fn query(
        &self,
        m: usize,
        x_lo: u32,
        x_hi: u32,
        y_lo: i32,
        y_hi: i32,
        own: &[Cand],
        ys: &[i32],
        best: &mut Cand,
    ) -> u64 {
        let mut visits = 0u64;
        if self.root == NONE || x_lo > x_hi || y_lo > y_hi {
            return visits;
        }
        let mut stack = [0u32; 128];
        let mut sp = 0usize;
        stack[sp] = self.root;
        sp += 1;
        while sp > 0 {
            sp -= 1;
            let v = stack[sp];
            visits += 1;
            let node = &self.nodes[v as usize];
            if node.max[m] <= *best {
                continue;
            }
            if node.x_hi < x_lo || node.x_lo > x_hi || node.y_hi < y_lo || node.y_lo > y_hi {
                continue;
            }
            if node.x_lo >= x_lo && node.x_hi <= x_hi && node.y_lo >= y_lo && node.y_hi <= y_hi {
                *best = node.max[m];
                continue;
            }
            let p = node.point;
            let val = own[p as usize];
            if val > *best && p >= x_lo && p <= x_hi {
                let y = ys[p as usize];
                if y >= y_lo && y <= y_hi {
                    *best = val;
                }
            }
            // Best-first: explore the child with the larger subtree maximum first so
            // that the other child is more likely to be pruned.
            let (l, r) = (node.left, node.right);
            if l != NONE && r != NONE {
                if self.nodes[l as usize].max[m] > self.nodes[r as usize].max[m] {
                    stack[sp] = r;
                    stack[sp + 1] = l;
                } else {
                    stack[sp] = l;
                    stack[sp + 1] = r;
                }
                sp += 2;
            } else if l != NONE {
                stack[sp] = l;
                sp += 1;
            } else if r != NONE {
                stack[sp] = r;
                sp += 1;
            }
        }
        visits
    }
}

// ---------------------------------------------------------------------------
// 1-D maximum tree over list index for NCBI's `best[]` selection
// (link_hsps.c:610-623 and :759-763/:887-891: `>=`, so the last list entry wins ties).
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Dual {
    sum: [i32; 2],
    idx: [u32; 2],
}

impl Dual {
    const INACTIVE: Dual = Dual {
        sum: [0; 2],
        idx: [NONE; 2],
    };

    #[inline]
    fn reduce(l: Dual, r: Dual) -> Dual {
        let mut out = l;
        for c in 0..2 {
            if r.idx[c] != NONE
                && (l.idx[c] == NONE
                    || r.sum[c] > l.sum[c]
                    || (r.sum[c] == l.sum[c] && r.idx[c] > l.idx[c]))
            {
                out.sum[c] = r.sum[c];
                out.idx[c] = r.idx[c];
            }
        }
        out
    }
}

struct DualTree {
    cap: usize,
    nodes: Vec<Dual>,
}

impl DualTree {
    fn new(n: usize) -> Self {
        let cap = n.max(1).next_power_of_two();
        DualTree {
            cap,
            nodes: vec![Dual::INACTIVE; cap * 2],
        }
    }

    #[inline]
    fn update(&mut self, i: usize, leaf: Dual) {
        let mut p = self.cap + i;
        self.nodes[p] = leaf;
        while p > 1 {
            p /= 2;
            self.nodes[p] = Dual::reduce(self.nodes[2 * p], self.nodes[2 * p + 1]);
        }
    }

    #[inline]
    fn maxima(&self, cutoffs: [i32; 2], ignore_small_gaps: bool) -> ([Option<usize>; 2], [i32; 2]) {
        let root = self.nodes[1];
        let mut best = [None; 2];
        let mut sums = [-cutoffs[0], -cutoffs[1]];
        for c in usize::from(ignore_small_gaps)..2 {
            if root.idx[c] != NONE && root.sum[c] >= sums[c] {
                best[c] = Some(root.idx[c] as usize);
                sums[c] = root.sum[c];
            }
        }
        (best, sums)
    }
}

// ---------------------------------------------------------------------------
// Per-group linking state
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
struct Tuple {
    sum: i32,
    num: i16,
    xsum: f64,
    link: u32,
}

#[derive(Default, Debug)]
struct Stats {
    rounds: u64,
    recompute_events: u64,
    propagated: u64,
    queries: u64,
    kd_visits: u64,
    zone_scanned: u64,
    chain_len_total: u64,
}

struct Linker {
    stats: Stats,
    // static per-HSP data in list order (link_hsps.c:540-551)
    q_off: Vec<i32>,
    q_off_t: Vec<i32>,
    s_off_t: Vec<i32>,
    q_end_t: Vec<i32>,
    s_end_t: Vec<i32>,
    score: Vec<i32>,
    lambda: Vec<f64>,
    log_k: Vec<f64>,
    cutoff: [i32; 2],
    ignore_small_gaps: bool,
    // current ("true") DP values, equal to a from-scratch NCBI pass on the active set
    cur: [Vec<Tuple>; 2],
    active: Vec<bool>,
    indeg: Vec<i32>,
    // NCBI's stored values: a snapshot taken at the last recomputation pass
    version: u32,
    snap_version: Vec<u32>,
    snap: [Vec<Tuple>; 2],
    snap_indeg: Vec<i32>,
    snapped: Vec<u32>,
    // reverse links (children) per method, intrusive doubly linked lists
    first_child: [Vec<u32>; 2],
    next_sib: [Vec<u32>; 2],
    prev_sib: [Vec<u32>; 2],
    // range-maximum structures
    kd: KdTree,
    kd_val: [Vec<Cand>; 2],
    stored_tree: DualTree,
    // propagation work list
    heap: BinaryHeap<Reverse<u32>>,
    queued: Vec<bool>,
}

impl Linker {
    #[inline]
    fn stored(&self, i: usize, m: usize) -> Tuple {
        if self.snap_version[i] == self.version {
            self.snap[m][i]
        } else {
            self.cur[m][i]
        }
    }

    #[inline]
    fn stored_indeg(&self, i: usize) -> i32 {
        if self.snap_version[i] == self.version {
            self.snap_indeg[i]
        } else {
            self.indeg[i]
        }
    }

    /// Copy-on-write snapshot before the first modification after a recomputation.
    #[inline]
    fn ensure_snapshot(&mut self, i: usize) {
        if self.snap_version[i] != self.version {
            self.snap_version[i] = self.version;
            self.snap[0][i] = self.cur[0][i];
            self.snap[1][i] = self.cur[1][i];
            self.snap_indeg[i] = self.indeg[i];
            self.snapped.push(i as u32);
        }
    }

    /// NCBI recomputation pass (link_hsps.c:659-899): stored values become the current ones.
    fn recompute_event(&mut self) {
        self.stats.recompute_events += 1;
        self.version += 1;
        let snapped = std::mem::take(&mut self.snapped);
        for &i in &snapped {
            let i = i as usize;
            if self.active[i] {
                let leaf = Dual {
                    sum: [self.cur[0][i].sum, self.cur[1][i].sum],
                    idx: [i as u32, i as u32],
                };
                self.stored_tree.update(i, leaf);
            }
        }
        self.snapped = snapped;
        self.snapped.clear();
    }

    #[inline]
    fn kd_value(&self, i: usize, m: usize) -> Cand {
        if self.active[i] && self.cur[m][i].sum > 0 {
            Cand {
                sum: self.cur[m][i].sum,
                idx: i as u32,
            }
        } else {
            Cand::NEG
        }
    }

    #[inline]
    fn set_kd(&mut self, i: usize, m: usize) {
        let v = self.kd_value(i, m);
        let mut own = std::mem::take(&mut self.kd_val[m]);
        self.kd.update(i as u32, m, v, &mut own);
        self.kd_val[m] = own;
    }

    #[inline]
    fn child_insert(&mut self, m: usize, parent: usize, child: usize) {
        let head = self.first_child[m][parent];
        self.next_sib[m][child] = head;
        self.prev_sib[m][child] = NONE;
        if head != NONE {
            self.prev_sib[m][head as usize] = child as u32;
        }
        self.first_child[m][parent] = child as u32;
    }

    #[inline]
    fn child_remove(&mut self, m: usize, parent: usize, child: usize) {
        let prev = self.prev_sib[m][child];
        let next = self.next_sib[m][child];
        if prev != NONE {
            self.next_sib[m][prev as usize] = next;
        } else {
            self.first_child[m][parent] = next;
        }
        if next != NONE {
            self.prev_sib[m][next as usize] = prev;
        }
        self.next_sib[m][child] = NONE;
        self.prev_sib[m][child] = NONE;
    }

    /// First list index whose `q_off` is `<= x` (the list is non-increasing in `q_off`).
    #[inline]
    fn first_index_with_q_off_le(&self, x: i32) -> usize {
        self.q_off.partition_point(|&q| q > x)
    }

    /// Best predecessor of `i` for method `m`: `argmax (sum, index)` over active `j < i`
    /// with `sum > 0` and the method's coordinate conditions
    /// (link_hsps.c:719-745 for small gaps, :838-860 for large gaps).
    fn best_predecessor(&self, i: usize, m: usize) -> (Cand, u64, u64) {
        let a = self.q_end_t[i];
        let b = self.s_end_t[i];
        let mut best = Cand::NEG;
        let mut visits = 0u64;
        let mut zone = 0u64;
        // Core: HSPs whose q_off alone guarantees the q_off_trim condition.
        // q_off_trim in [q_off, q_off + TRIM_SIZE].
        let p = self.first_index_with_q_off_le(a); // q_off > a  <=>  index < p
        if m == 1 {
            if p > 0 {
                visits += self.kd.query(
                    1,
                    0,
                    (p - 1) as u32,
                    b + 1,
                    i32::MAX,
                    &self.kd_val[1],
                    &self.s_off_t,
                    &mut best,
                );
            }
        } else {
            let p_lo = self.first_index_with_q_off_le(a + WINDOW_SIZE - TRIM_SIZE);
            if p_lo < p {
                visits += self.kd.query(
                    0,
                    p_lo as u32,
                    (p - 1) as u32,
                    b + 1,
                    b + WINDOW_SIZE,
                    &self.kd_val[0],
                    &self.s_off_t,
                    &mut best,
                );
            }
        }
        // Boundary zones where q_off_trim may or may not satisfy the condition.
        let mut scan_zone = |lo_excl: i32, hi_incl: i32, best: &mut Cand, zone: &mut u64| {
            // indices with q_off in (lo_excl, hi_incl]
            let start = self.first_index_with_q_off_le(hi_incl);
            let end = self.first_index_with_q_off_le(lo_excl);
            for j in start..end.min(i) {
                *zone += 1;
                if !self.active[j] {
                    continue;
                }
                let sum = self.cur[m][j].sum;
                if sum <= 0 {
                    continue;
                }
                let c = Cand { sum, idx: j as u32 };
                if c <= *best {
                    continue;
                }
                let qo = self.q_off_t[j];
                let so = self.s_off_t[j];
                let ok = if m == 1 {
                    qo > a && so > b
                } else {
                    qo > a && so > b && qo <= a + WINDOW_SIZE && so <= b + WINDOW_SIZE
                };
                if ok {
                    *best = c;
                }
            }
        };
        scan_zone(a - TRIM_SIZE, a, &mut best, &mut zone);
        if m == 0 {
            scan_zone(
                a + WINDOW_SIZE - TRIM_SIZE,
                a + WINDOW_SIZE,
                &mut best,
                &mut zone,
            );
        }
        (best, visits, zone)
    }

    /// DP value of `i` for method `m` from the current values of the earlier HSPs
    /// (link_hsps.c:748-757 and :863-872).
    fn compute(&mut self, i: usize, m: usize) -> Tuple {
        let score = self.score[i];
        let mut h_sum = 0i32;
        let mut h_num = 0i16;
        let mut h_xsum = 0.0f64;
        let mut h_link = NONE;
        if score > self.cutoff[m] {
            let (c, visits, zone) = self.best_predecessor(i, m);
            self.stats.queries += 1;
            self.stats.kd_visits += visits;
            self.stats.zone_scanned += zone;
            if c != Cand::NEG {
                let j = c.idx as usize;
                let t = self.cur[m][j];
                h_sum = t.sum;
                h_num = t.num;
                h_xsum = t.xsum;
                h_link = c.idx;
            }
        }
        // NCBI link_hsps.c:750-753 / :865-868, same left-to-right evaluation order.
        let new_xsum = h_xsum + (score as f64) * self.lambda[i] - self.log_k[i];
        let new_sum = h_sum + (score - self.cutoff[m]);
        Tuple {
            sum: new_sum,
            num: h_num + 1,
            xsum: new_xsum,
            link: h_link,
        }
    }

    /// Install a new DP value for `i`/`m`, maintaining reverse links, in-degrees,
    /// range structures and snapshots. Returns true when the value changed.
    fn install(&mut self, i: usize, m: usize, t: Tuple, initial: bool) -> bool {
        let old = self.cur[m][i];
        if !initial && old == t {
            return false;
        }
        if !initial {
            self.ensure_snapshot(i);
        }
        if old.link != t.link {
            if !initial && old.link != NONE {
                let p = old.link as usize;
                self.child_remove(m, p, i);
                if !initial {
                    self.ensure_snapshot(p);
                }
                self.indeg[p] -= 1;
            }
            if t.link != NONE {
                let p = t.link as usize;
                self.child_insert(m, p, i);
                if !initial {
                    self.ensure_snapshot(p);
                }
                self.indeg[p] += 1;
            }
        }
        self.cur[m][i] = t;
        if initial || old.sum != t.sum {
            self.set_kd(i, m);
        }
        true
    }

    /// Recompute every HSP whose chain passed through a removed HSP, in list order.
    fn propagate(&mut self) {
        while let Some(Reverse(x)) = self.heap.pop() {
            let x = x as usize;
            self.queued[x] = false;
            if !self.active[x] {
                continue;
            }
            self.stats.propagated += 1;
            for m in 0..2 {
                let t = self.compute(x, m);
                if self.install(x, m, t, false) {
                    let mut c = self.first_child[m][x];
                    while c != NONE {
                        let cu = c as usize;
                        if !self.queued[cu] && self.active[cu] {
                            self.queued[cu] = true;
                            self.heap.push(Reverse(c));
                        }
                        c = self.next_sib[m][cu];
                    }
                }
            }
        }
    }

    fn remove(&mut self, i: usize) {
        self.active[i] = false;
        self.stored_tree.update(i, Dual::INACTIVE);
        for m in 0..2 {
            self.set_kd(i, m);
            let mut c = self.first_child[m][i];
            while c != NONE {
                let cu = c as usize;
                if !self.queued[cu] && self.active[cu] {
                    self.queued[cu] = true;
                    self.heap.push(Reverse(c));
                }
                c = self.next_sib[m][cu];
            }
            let p = self.cur[m][i].link;
            if p != NONE {
                let p = p as usize;
                self.child_remove(m, p, i);
                self.ensure_snapshot(p);
                self.indeg[p] -= 1;
            }
        }
    }

    /// NCBI link_hsps.c:644-650: walk the stored chain; false if a removed HSP is met.
    fn stored_path_intact(&self, start: usize, m: usize) -> bool {
        let mut x = start;
        loop {
            if !self.active[x] {
                return false;
            }
            let next = self.stored(x, m).link;
            if next == NONE {
                return true;
            }
            x = next as usize;
        }
    }
}

/// Drop-in replacement for `link_hsp_group_ncbi` (same inputs, same output fields).
#[allow(clippy::too_many_arguments)]
pub(super) fn link_hsp_group_fast(
    mut group_hits: Vec<UngappedHit>,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
) -> Vec<UngappedHit> {
    if group_hits.is_empty() {
        return group_hits;
    }
    let n = group_hits.len();
    // LOSAT_TIMING only adds the time of the group to the line printed below.
    let started = std::env::var_os("LOSAT_TIMING")
        .is_some()
        .then(std::time::Instant::now);

    // Effective lengths and search space, as in link_hsps.c:559-571.
    let query_context = group_hits[0].ctx_idx;
    let query_len_aa = query_contexts[query_context].aa_len as i64;
    let subject_len_aa = (subject_len_nucl / 3).max(1);
    let length_adjustment = length_adj_per_context[query_context];
    let eff_search_space = eff_searchsp_per_context[query_context];
    let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;
    let length_adj_for_subject = length_adjustment / 3;
    let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;

    let cutoff = [cutoffs.cutoff_small_gap, cutoffs.cutoff_big_gap];
    let gap_prob = cutoffs.gap_prob;
    let ignore_small_gaps = cutoffs.ignore_small_gaps;

    // Static per-HSP data (link_hsps.c:540-551).
    let mut q_off = Vec::with_capacity(n);
    let mut q_off_t = Vec::with_capacity(n);
    let mut s_off_t = Vec::with_capacity(n);
    let mut q_end_t = Vec::with_capacity(n);
    let mut s_end_t = Vec::with_capacity(n);
    let mut score = Vec::with_capacity(n);
    let mut lambda = Vec::with_capacity(n);
    let mut log_k = Vec::with_capacity(n);
    for hit in &group_hits {
        let (qo, qe, so, se) = (
            hit.q_aa_start as i32,
            hit.q_aa_end as i32,
            hit.s_aa_start as i32,
            hit.s_aa_end as i32,
        );
        let qt = TRIM_SIZE.min((qe - qo) / 4);
        let st = TRIM_SIZE.min((se - so) / 4);
        q_off.push(qo);
        q_off_t.push(qo + qt);
        s_off_t.push(so + st);
        q_end_t.push(qe - qt);
        s_end_t.push(se - st);
        score.push(hit.raw_score);
        lambda.push(query_contexts[hit.ctx_idx].karlin_params.lambda);
        log_k.push(log_k_by_ctx[hit.ctx_idx]);
    }
    debug_assert!(q_off.windows(2).all(|w| w[0] >= w[1]));

    let kd = KdTree::build(&s_off_t);
    let empty = Tuple {
        sum: 0,
        num: 0,
        xsum: 0.0,
        link: NONE,
    };
    let mut lk = Linker {
        stats: Stats::default(),
        q_off,
        q_off_t,
        s_off_t,
        q_end_t,
        s_end_t,
        score,
        lambda,
        log_k,
        cutoff,
        ignore_small_gaps,
        cur: [vec![empty; n], vec![empty; n]],
        active: vec![true; n],
        indeg: vec![0; n],
        version: 0,
        snap_version: vec![u32::MAX; n],
        snap: [vec![empty; n], vec![empty; n]],
        snap_indeg: vec![0; n],
        snapped: Vec::new(),
        first_child: [vec![NONE; n], vec![NONE; n]],
        next_sib: [vec![NONE; n], vec![NONE; n]],
        prev_sib: [vec![NONE; n], vec![NONE; n]],
        kd,
        kd_val: [vec![Cand::NEG; n], vec![Cand::NEG; n]],
        stored_tree: DualTree::new(n),
        heap: BinaryHeap::new(),
        queued: vec![false; n],
    };

    // First NCBI pass (link_hsps.c:691-896): every HSP in list order.
    for i in 0..n {
        for m in 0..2 {
            let t = lk.compute(i, m);
            lk.install(i, m, t, true);
        }
        let leaf = Dual {
            sum: [lk.cur[0][i].sum, lk.cur[1][i].sum],
            idx: [i as u32, i as u32],
        };
        lk.stored_tree.update(i, leaf);
    }
    lk.version = 1; // stored == current after the first pass

    let int4_max = i32::MAX as f64;
    let mut remaining = n;
    let mut first_pass = true;
    let mut path_changed = true;

    while remaining > 0 {
        let mut best: [Option<usize>; 2] = [None, None];
        let mut use_current_max = false;
        if !first_pass {
            // link_hsps.c:603-652 on NCBI's stored sums
            let (b, _) = lk.stored_tree.maxima(cutoff, ignore_small_gaps);
            best = b;
            if !path_changed {
                use_current_max = true;
            } else {
                use_current_max = true;
                if !ignore_small_gaps {
                    if let Some(b0) = best[0] {
                        if !lk.stored_path_intact(b0, 0) {
                            use_current_max = false;
                        }
                    }
                }
                if use_current_max {
                    if let Some(b1) = best[1] {
                        if !lk.stored_path_intact(b1, 1) {
                            use_current_max = false;
                        }
                    }
                }
            }
        }
        if !use_current_max {
            // link_hsps.c:659-899: NCBI recomputes every stored value. The current
            // values of the HSPs whose chains were cut since the last pass are brought
            // up to date here (they are not needed earlier: between two passes NCBI
            // only reads the stored values).
            if !first_pass {
                lk.propagate();
                lk.recompute_event();
            }
            let (b, _) = lk.stored_tree.maxima(cutoff, ignore_small_gaps);
            best = b;
            path_changed = false;
            first_pass = false;
        }

        // link_hsps.c:901-953
        let mut prob = [f64::MAX, f64::MAX];
        if !ignore_small_gaps {
            if let Some(bi) = best[0] {
                // link_hsps.c:907-908, a persistent modification of the stored sum
                lk.ensure_snapshot(bi);
                let num_i = lk.snap[0][bi].num as i32;
                lk.snap[0][bi].sum += num_i * cutoff[0];
                let leaf = Dual {
                    sum: [lk.snap[0][bi].sum, lk.snap[1][bi].sum],
                    idx: [bi as u32, bi as u32],
                };
                lk.stored_tree.update(bi, leaf);
                let t = lk.stored(bi, 0);
                let num = t.num as usize;
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[0] = small_gap_sum_e(
                    WINDOW_SIZE,
                    num as i16,
                    t.xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space,
                    divisor,
                );
                if num > 1 {
                    if gap_prob == 0.0 || (prob[0] / gap_prob) > int4_max {
                        prob[0] = int4_max;
                    } else {
                        prob[0] /= gap_prob;
                    }
                }
            }
            if let Some(bi) = best[1] {
                let t = lk.stored(bi, 1);
                let num = t.num as usize;
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[1] = ncbi_large_gap_sum_e(
                    num as i16,
                    t.xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space,
                    divisor,
                );
                if num > 1 {
                    let denom = 1.0 - gap_prob;
                    if denom == 0.0 || (prob[1] / denom) > int4_max {
                        prob[1] = int4_max;
                    } else {
                        prob[1] /= denom;
                    }
                }
            }
        } else if let Some(bi) = best[1] {
            lk.ensure_snapshot(bi);
            let num_i = lk.snap[1][bi].num as i32;
            lk.snap[1][bi].sum += num_i * cutoff[1];
            let leaf = Dual {
                sum: [lk.snap[0][bi].sum, lk.snap[1][bi].sum],
                idx: [bi as u32, bi as u32],
            };
            lk.stored_tree.update(bi, leaf);
            let t = lk.stored(bi, 1);
            let num = t.num as usize;
            let divisor = gap_decay_divisor(gap_decay_rate, num);
            prob[1] = ncbi_large_gap_sum_e(
                num as i16,
                t.xsum,
                eff_query_len as i32,
                eff_subject_len as i32,
                eff_search_space,
                divisor,
            );
            if num > 1 {
                let denom = 1.0 - gap_prob;
                if denom == 0.0 || (prob[1] / denom) > int4_max {
                    prob[1] = int4_max;
                } else {
                    prob[1] /= denom;
                }
            }
        }

        let ordering = if !ignore_small_gaps && prob[0] <= prob[1] {
            0
        } else {
            1
        };
        let best_i = match best[ordering].or(best[1 - ordering]) {
            Some(i) => i,
            None => break,
        };
        let evalue = prob[ordering];

        // link_hsps.c:955-981
        let head = lk.stored(best_i, ordering);
        let linked_set = head.link != NONE;
        if lk.stored_indeg(best_i) > 0 {
            path_changed = true;
        }
        group_hits[best_i].start_of_chain = true;
        group_hits[best_i].hsp_link_num = head.num;

        let mut chain: Vec<usize> = Vec::new();
        let mut cur = best_i;
        loop {
            chain.push(cur);
            let next = lk.stored(cur, ordering).link;
            if next == NONE {
                break;
            }
            cur = next as usize;
        }
        for (k, &cur) in chain.iter().enumerate() {
            if lk.stored_indeg(cur) > 1 {
                path_changed = true;
            }
            group_hits[cur].ordering_method = ordering as u8;
            group_hits[cur].e_value = evalue;
            group_hits[cur].linked_set = linked_set;
            group_hits[cur].chain_next_link_id =
                chain.get(k + 1).map(|&next| group_hits[next].link_id);
            if k > 0 {
                group_hits[cur].start_of_chain = false;
            }
        }
        for &cur in &chain {
            lk.remove(cur);
        }
        remaining -= chain.len();
        lk.stats.rounds += 1;
        lk.stats.chain_len_total += chain.len() as u64;
    }
    if let Some(started) = started {
        eprintln!(
            "[TIMING] linking_fast group n={} us={} {:?}",
            n,
            started.elapsed().as_micros(),
            lk.stats
        );
    }

    group_hits
}

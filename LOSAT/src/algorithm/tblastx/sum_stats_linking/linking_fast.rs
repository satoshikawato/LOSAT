//! Index-backed replacement of the per-group link_hsps kernel
//! (`LOSAT_LINK_FAST=1`; the default path is unchanged).
//!
//! It keeps NCBI's round structure (first_pass / path_changed /
//! use_current_max / linked_to / changed) and the same selection, statistics
//! and removal statements as `link_hsp_group_ncbi`, and answers the two
//! predecessor searches with index structures instead of linear scans.
//!
//! # What a search returns
//!
//! NCBI scans the remaining HSPs before H in list order, nearest first, and
//! keeps a candidate when `sum > H_hsp_sum` (link_hsps.c:727, 857). The result
//! is the candidate with the largest sum and, among equal sums, the one
//! nearest to H: the largest post-sort index. `key(sum, index)` orders
//! candidates exactly that way, so every search here is "largest key".
//!
//! # Index 0 (small gaps, link_hsps.c:691-768)
//!
//! The candidates of H are the remaining HSPs before H whose trimmed start
//! lies in (qe, qe+W] x (se, se+W] and whose sum is positive. They are found
//! in a uniform grid of W x W cells; each cell lists its HSPs by increasing
//! index and the walk stops at H. NCBI's early `break` (line 717) only skips
//! HSPs that fail `q_off_t > H_q_et_gap` and has no counterpart.
//!
//! NCBI searches again in every recomputation pass. This kernel keeps the
//! choice of the previous pass when that choice still stands unchanged
//! (on unless `LOSAT_LINK_FAST_REUSE0=0`), by the argument NCBI uses for
//! index 1 (lines 781-795). Write S_r(j) for sum[0] of j as pass r computes
//! it (the add-back of line 907 is overwritten for every remaining HSP by the
//! next pass, lines 755-758).
//!
//! 1. The window test is fixed and HSPs are only removed, so the candidates
//!    of H in pass r are a subset of those in pass r-1.
//! 2. S_r(j) <= S_{r-1}(j): S(j) = score(j) - cutoff + max(0, max of S(k)
//!    over the candidates k of j); by induction in list order this is a
//!    maximum over fewer and not larger values. (No sum leaves the Int4
//!    range: `sums_fit_int4` is checked before this kernel is used.)
//! 3. `changed0[p] == false` means p kept its own choice in pass r, whose sum
//!    is unchanged by induction, so S_r(p) = S_{r-1}(p).
//!
//! If p, the choice of pass r-1, remains with `changed0[p] == false`, every
//! other candidate k has S_r(k) <= S_{r-1}(k) <= S_{r-1}(p) = S_r(p), and
//! where S_{r-1}(k) = S_{r-1}(p) the scan of pass r-1 met p first; the scan
//! of pass r therefore selects p again and reads the same num/sum/xsum from
//! it. If pass r-1 selected nothing, no candidate had a positive sum then,
//! and none has now. In every other case the search runs.
//!
//! # Index 1 (large gaps, link_hsps.c:771-896)
//!
//! The candidates of H are the remaining HSPs before H with qo > qe(H) and
//! so > se(H). NCBI's own reuse rule (lines 781-795) is kept as is. A search
//! visits the HSPs by decreasing qe: every candidate of H has qe >= qo >
//! qe(H), so it has been given its sum for this pass before H is visited, as
//! in NCBI's list-order loop. "qo > qe(H)" is then a prefix of the HSPs by
//! decreasing qo, and a prefix-maximum Fenwick tree over the ranks of so
//! answers "so > se(H)". NCBI's starting value `sum - 1` (lines 812-823)
//! only lowers the threshold below a candidate that is still present and
//! does not change the selected HSP.
//!
//! The tree also holds HSPs after H in list order. One of those can pass the
//! coordinate test only if H spans at most 4 residues: qo(j) is at most
//! q_off(H) - 1 + 5, and qe(H) is at least q_off(H) + 4 from 5 residues on.
//! Such an H appears when `Blast_HSPReevaluateWithAmbiguitiesUngapped` trims
//! an HSP. When the tree answers with an HSP after H, the search falls back
//! to NCBI's scan.
//!
//! # Checking
//!
//! `LOSAT_LINK_FAST_VERIFY=1` repeats every choice - searched or kept, both
//! indexes - with a plain NCBI-order scan and compares the selected HSP, its
//! sum, its number of HSPs and the bits of its xsum.
//! `LOSAT_LINK_FAST_SHADOW=1` (in `linking.rs`) also runs the NCBI kernel on
//! every group and compares the complete result.

use std::sync::OnceLock;

use rustc_hash::FxHashMap;

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

pub(super) fn link_fast_enabled() -> bool {
    static FLAG: OnceLock<bool> = OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("LOSAT_LINK_FAST").is_some())
}

/// `LOSAT_LINK_FAST_SHADOW=1`: `linking.rs` also runs the NCBI kernel and
/// compares the complete result of every group.
pub(super) fn link_fast_shadow_enabled() -> bool {
    static FLAG: OnceLock<bool> = OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("LOSAT_LINK_FAST_SHADOW").is_some())
}

fn stats_enabled() -> bool {
    static FLAG: OnceLock<bool> = OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("LOSAT_LINK_STATS").is_some())
}

/// How the kernel runs. `from_env` is what a search uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LinkFastOptions {
    /// Repeat every choice with a plain NCBI-order scan and compare.
    pub verify: bool,
    /// Keep an unchanged index-0 choice without searching again. `false`
    /// searches in every pass, as NCBI does.
    pub reuse_index0: bool,
}

impl LinkFastOptions {
    pub(super) fn from_env() -> Self {
        static OPTIONS: OnceLock<LinkFastOptions> = OnceLock::new();
        *OPTIONS.get_or_init(|| LinkFastOptions {
            verify: std::env::var_os("LOSAT_LINK_FAST_VERIFY").is_some(),
            reuse_index0: std::env::var("LOSAT_LINK_FAST_REUSE0").map_or(true, |v| v != "0"),
        })
    }
}

/// What one call did. Tests assert on it; `LOSAT_LINK_STATS=1` prints it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct LinkFastStats {
    pub rounds: u64,
    pub recompute_rounds: u64,
    pub recompute_alive_sum: u64,
    /// Index-0 choices made by a grid search.
    pub idx0_searched: u64,
    /// Grid entries read by those searches.
    pub idx0_scanned: u64,
    /// Index-0 choices kept from the previous pass: an HSP / nothing.
    pub idx0_kept_link: u64,
    pub idx0_kept_none: u64,
    /// Index-1 choices made by a tree search.
    pub idx1_searched: u64,
    /// Index-1 choices kept from the previous pass: an HSP / nothing.
    pub idx1_kept_link: u64,
    pub idx1_kept_none: u64,
    /// Index-1 searches answered by NCBI's scan because the tree returned an
    /// HSP after H in list order.
    pub fallbacks: u64,
    /// Choices compared with the plain scan (`verify`).
    pub verified: u64,
}

impl LinkFastStats {
    #[cfg(test)]
    pub(super) fn add(&mut self, other: &LinkFastStats) {
        self.rounds += other.rounds;
        self.recompute_rounds += other.recompute_rounds;
        self.recompute_alive_sum += other.recompute_alive_sum;
        self.idx0_searched += other.idx0_searched;
        self.idx0_scanned += other.idx0_scanned;
        self.idx0_kept_link += other.idx0_kept_link;
        self.idx0_kept_none += other.idx0_kept_none;
        self.idx1_searched += other.idx1_searched;
        self.idx1_kept_link += other.idx1_kept_link;
        self.idx1_kept_none += other.idx1_kept_none;
        self.fallbacks += other.fallbacks;
        self.verified += other.verified;
    }
}

/// True when no chain sum of this group can leave the Int4 range, which the
/// index-0 reuse argument (module documentation, step 2) assumes. A chain sum
/// adds `score - cutoff` over distinct HSPs, so the sum of all positive parts
/// bounds it. `linking.rs` uses the NCBI kernel for a group that fails this.
pub(super) fn sums_fit_int4(group_hits: &[UngappedHit], cutoffs: &LinkHspCutoffs) -> bool {
    let slack = i64::from(cutoffs.cutoff_small_gap.min(0).unsigned_abs())
        + i64::from(cutoffs.cutoff_big_gap.min(0).unsigned_abs());
    let mut bound = 0i64;
    for hit in group_hits {
        bound += i64::from(hit.raw_score.max(0)) + slack;
        if bound > i64::from(i32::MAX) {
            return false;
        }
    }
    true
}

/// (sum, index) packed so that the integer order is: larger sum first, then
/// larger index. 0 is "no entry".
#[inline(always)]
fn key(sum: i32, idx: usize) -> u64 {
    ((((sum as i64) + (1i64 << 31)) as u64) << 32) | ((idx as u64) + 1)
}

#[inline(always)]
fn key_idx(k: u64) -> usize {
    ((k & 0xffff_ffff) - 1) as usize
}

/// Flat maximum tree over the stable post-sort indices, one channel per
/// ordering method (the existing `DualMaximumTree`, with packed keys and an
/// update that stops at the first unchanged parent).
struct MaxTree {
    cap: usize,
    nodes: Vec<[u64; 2]>,
}

impl MaxTree {
    fn new(n: usize) -> Self {
        let cap = n.max(1).next_power_of_two();
        Self {
            cap,
            nodes: vec![[0u64; 2]; cap * 2],
        }
    }

    #[inline]
    fn leaf(&self, i: usize) -> [u64; 2] {
        self.nodes[self.cap + i]
    }

    #[inline]
    fn set(&mut self, i: usize, v: [u64; 2]) {
        let mut pos = self.cap + i;
        self.nodes[pos] = v;
        while pos > 1 {
            pos >>= 1;
            let l = self.nodes[pos * 2];
            let r = self.nodes[pos * 2 + 1];
            let m = [l[0].max(r[0]), l[1].max(r[1])];
            if self.nodes[pos] == m {
                break;
            }
            self.nodes[pos] = m;
        }
    }

    fn build(&mut self) {
        for pos in (1..self.cap).rev() {
            let l = self.nodes[pos * 2];
            let r = self.nodes[pos * 2 + 1];
            self.nodes[pos] = [l[0].max(r[0]), l[1].max(r[1])];
        }
    }

    #[inline]
    fn root(&self) -> [u64; 2] {
        self.nodes[1]
    }
}

/// What a predecessor search hands to link_hsps.c:750-767 / 866-894: the
/// selected HSP (or `NONE`) and the three values read from it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Choice {
    link: u32,
    sum: i32,
    num: i16,
    xsum_bits: u64,
}

/// `LOSAT_LINK_FAST_VERIFY`: the index-0 choice of link_hsps.c:706-742, by
/// NCBI's own scan over the remaining HSPs before `i`, nearest first.
#[allow(clippy::too_many_arguments)]
fn verify_scan_index0(
    i: usize,
    prev_active: &[u32],
    qo: &[i32],
    so: &[i32],
    sum0: &[i32],
    num0: &[i16],
    xsum0: &[f64],
    h_qe: i32,
    h_se: i32,
) -> Choice {
    let h_qe_gap = h_qe + WINDOW_SIZE;
    let h_se_gap = h_se + WINDOW_SIZE;
    let mut choice = Choice {
        link: NONE,
        sum: 0,
        num: 0,
        xsum_bits: 0.0f64.to_bits(),
    };
    let mut j = prev_active[i];
    while j != NONE {
        let jj = j as usize;
        j = prev_active[jj];
        let q_off_t = qo[jj];
        let s_off_t = so[jj];
        // link_hsps.c:717
        if q_off_t > h_qe_gap + TRIM_SIZE {
            break;
        }
        // link_hsps.c:719-724
        if q_off_t <= h_qe || s_off_t <= h_se || q_off_t > h_qe_gap || s_off_t > h_se_gap {
            continue;
        }
        // link_hsps.c:727
        if sum0[jj] > choice.sum {
            choice = Choice {
                link: jj as u32,
                sum: sum0[jj],
                num: num0[jj],
                xsum_bits: xsum0[jj].to_bits(),
            };
        }
    }
    choice
}

/// `LOSAT_LINK_FAST_VERIFY`: the index-1 choice of link_hsps.c:827-863, by a
/// plain scan over the remaining HSPs before `i`, nearest first, starting
/// from `start_sum` (NCBI's `H_hsp_sum` before the loop). `next_larger` only
/// skips HSPs that fail `sum > H_hsp_sum` and has no counterpart.
#[allow(clippy::too_many_arguments)]
fn verify_scan_index1(
    i: usize,
    prev_active: &[u32],
    qo: &[i32],
    so: &[i32],
    sum1: &[i32],
    num1: &[i16],
    xsum1: &[f64],
    h_qe: i32,
    h_se: i32,
    start_sum: i32,
) -> Choice {
    let mut choice = Choice {
        link: NONE,
        sum: start_sum,
        num: 0,
        xsum_bits: 0.0f64.to_bits(),
    };
    let mut j = prev_active[i];
    while j != NONE {
        let jj = j as usize;
        j = prev_active[jj];
        // link_hsps.c:840-857: b0, b1, b2
        if sum1[jj] <= choice.sum || qo[jj] <= h_qe || so[jj] <= h_se {
            continue;
        }
        choice = Choice {
            link: jj as u32,
            sum: sum1[jj],
            num: num1[jj],
            xsum_bits: xsum1[jj].to_bits(),
        };
    }
    choice
}

/// The kernel `linking.rs` calls under `LOSAT_LINK_FAST=1`.
#[allow(clippy::too_many_arguments)]
pub(super) fn link_hsp_group_fast(
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
    let (group_hits, stats) = link_hsp_group_fast_with(
        group_hits,
        cutoffs,
        gap_decay_rate,
        subject_len_nucl,
        query_contexts,
        length_adj_per_context,
        eff_searchsp_per_context,
        log_k_by_ctx,
        LinkFastOptions::from_env(),
    );
    if stats_enabled() && n > 0 {
        eprintln!(
            "[LINK_FAST_STATS] n={} rounds={} recompute_rounds={} recompute_alive_sum={} idx0_searched={} idx0_scanned={} idx0_kept_link={} idx0_kept_none={} idx1_searched={} idx1_kept_link={} idx1_kept_none={} fallbacks={} verified={}",
            n,
            stats.rounds,
            stats.recompute_rounds,
            stats.recompute_alive_sum,
            stats.idx0_searched,
            stats.idx0_scanned,
            stats.idx0_kept_link,
            stats.idx0_kept_none,
            stats.idx1_searched,
            stats.idx1_kept_link,
            stats.idx1_kept_none,
            stats.fallbacks,
            stats.verified
        );
    }
    group_hits
}

#[allow(clippy::too_many_arguments)]
pub(super) fn link_hsp_group_fast_with(
    mut group_hits: Vec<UngappedHit>,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
    options: LinkFastOptions,
) -> (Vec<UngappedHit>, LinkFastStats) {
    let mut stats = LinkFastStats::default();
    let n = group_hits.len();
    if n == 0 {
        return (group_hits, stats);
    }
    let verify = options.verify;
    let reuse_index0 = options.reuse_index0;

    // link_hsps.c:559-571 (same expressions as link_hsp_group_ncbi)
    let query_context = group_hits[0].ctx_idx;
    let query_len_aa = query_contexts[query_context].aa_len as i64;
    let subject_len_aa = (subject_len_nucl / 3).max(1);
    let length_adjustment = length_adj_per_context[query_context];
    let eff_search_space = eff_searchsp_per_context[query_context];
    let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;
    let length_adj_for_subject = length_adjustment / 3;
    let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;

    let c0 = cutoffs.cutoff_small_gap;
    let c1 = cutoffs.cutoff_big_gap;
    let gap_prob = cutoffs.gap_prob;
    let ignore_small_gaps = cutoffs.ignore_small_gaps;

    // ---- static per-HSP data (link_hsps.c:545-549) ----
    let mut score = Vec::with_capacity(n);
    let mut qo = Vec::with_capacity(n);
    let mut so = Vec::with_capacity(n);
    let mut qe = Vec::with_capacity(n);
    let mut se = Vec::with_capacity(n);
    let mut sl = Vec::with_capacity(n); // score * Lambda
    let mut lk = Vec::with_capacity(n); // logK
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
        let lambda = query_contexts[hit.ctx_idx].karlin_params.lambda;
        sl.push((hit.raw_score as f64) * lambda);
        lk.push(log_k_by_ctx[hit.ctx_idx]);
    }

    // ---- dynamic per-HSP data ----
    let mut sum0: Vec<i32> = (0..n).map(|i| score[i] - c0).collect();
    let mut sum1: Vec<i32> = (0..n).map(|i| score[i] - c1).collect();
    let mut xsum0: Vec<f64> = (0..n).map(|i| sl[i] - lk[i]).collect();
    let mut xsum1: Vec<f64> = xsum0.clone();
    let mut num0: Vec<i16> = vec![1; n];
    let mut num1: Vec<i16> = vec![1; n];
    let mut link0: Vec<u32> = vec![NONE; n];
    let mut link1: Vec<u32> = vec![NONE; n];
    let mut changed0: Vec<bool> = vec![true; n];
    let mut changed1: Vec<bool> = vec![true; n];
    let mut linked_to: Vec<i32> = vec![0; n];
    let mut alive: Vec<bool> = vec![true; n];
    let mut next_active: Vec<u32> = (0..n)
        .map(|i| if i + 1 < n { (i + 1) as u32 } else { NONE })
        .collect();
    let mut prev_active: Vec<u32> = (0..n)
        .map(|i| if i > 0 { (i - 1) as u32 } else { NONE })
        .collect();
    let mut active_head: u32 = 0;

    // ---- index 0: grid of W x W cells over (qo, so), HSPs with score > c0 ----
    // cell entries are ascending post-sort indices.
    let w = WINDOW_SIZE;
    let cell_of = |q: i32, s: i32| -> u64 {
        (((q.div_euclid(w)) as u32 as u64) << 32) | ((s.div_euclid(w)) as u32 as u64)
    };
    let mut grid_index: FxHashMap<u64, (u32, u32)> = FxHashMap::default();
    let mut grid_items: Vec<u32> = Vec::new();
    if !ignore_small_gaps {
        let mut cells: Vec<(u64, u32)> = (0..n)
            .filter(|&i| score[i] > c0)
            .map(|i| (cell_of(qo[i], so[i]), i as u32))
            .collect();
        cells.sort_unstable();
        grid_items.reserve(cells.len());
        let mut start = 0usize;
        while start < cells.len() {
            let c = cells[start].0;
            let mut end = start;
            while end < cells.len() && cells[end].0 == c {
                grid_items.push(cells[end].1);
                end += 1;
            }
            grid_index.insert(c, (start as u32, end as u32));
            start = end;
        }
    }

    // ---- index 1: orders and ranks over the HSPs with score > c1 ----
    let elig1: Vec<u32> = (0..n)
        .filter(|&i| score[i] > c1)
        .map(|i| i as u32)
        .collect();
    let mut qorder: Vec<u32> = elig1.clone();
    qorder.sort_unstable_by(|&a, &b| qe[b as usize].cmp(&qe[a as usize]).then_with(|| a.cmp(&b)));
    let mut iorder: Vec<u32> = elig1.clone();
    iorder.sort_unstable_by(|&a, &b| qo[b as usize].cmp(&qo[a as usize]).then_with(|| a.cmp(&b)));
    let mut so_vals: Vec<i32> = elig1.iter().map(|&i| so[i as usize]).collect();
    so_vals.sort_unstable();
    so_vals.dedup();
    let m = so_vals.len();
    // rr = position (1-based) in decreasing-so order; "so > B" is a prefix.
    let mut rr: Vec<u32> = vec![0; n];
    let mut kq: Vec<u32> = vec![0; n];
    for &i in elig1.iter() {
        let i = i as usize;
        let rank = so_vals.partition_point(|&v| v < so[i]); // 0-based ascending
        rr[i] = (m - rank) as u32; // 1..=m
        let le = so_vals.partition_point(|&v| v <= se[i]);
        kq[i] = (m - le) as u32; // number of distinct so values > se[i]
    }
    let mut fen: Vec<u64> = vec![0; m + 1];

    let mut tree = MaxTree::new(n);

    let mut remaining = n;
    let mut first_pass = true;
    let mut path_changed = true;
    let int4_max = i32::MAX as f64;

    while remaining > 0 {
        stats.rounds += 1;
        let mut best: [Option<usize>; 2] = [None, None];
        let mut use_current_max = false;

        // link_hsps.c:603-652
        if !first_pass {
            let root = tree.root();
            for ch in usize::from(ignore_small_gaps)..2 {
                if root[ch] != 0 {
                    best[ch] = Some(key_idx(root[ch]));
                }
            }
            if !path_changed {
                use_current_max = true;
            } else {
                use_current_max = true;
                if !ignore_small_gaps {
                    if let Some(bi) = best[0] {
                        let mut cur = bi;
                        loop {
                            if linked_to[cur] == -1000 {
                                use_current_max = false;
                                break;
                            }
                            let p = link0[cur];
                            if p == NONE {
                                break;
                            }
                            cur = p as usize;
                        }
                    }
                }
                if use_current_max {
                    if let Some(bi) = best[1] {
                        let mut cur = bi;
                        loop {
                            if linked_to[cur] == -1000 {
                                use_current_max = false;
                                break;
                            }
                            let p = link1[cur];
                            if p == NONE {
                                break;
                            }
                            cur = p as usize;
                        }
                    }
                }
            }
        }

        if !use_current_max {
            stats.recompute_rounds += 1;

            // link_hsps.c:685: H->linked_to = 0;
            let mut cur = active_head;
            let mut alive_count = 0u64;
            while cur != NONE {
                linked_to[cur as usize] = 0;
                cur = next_active[cur as usize];
                alive_count += 1;
            }
            stats.recompute_alive_sum += alive_count;

            // ---------------- index 0 (link_hsps.c:691-768) ----------------
            if !ignore_small_gaps {
                let mut cur = active_head;
                while cur != NONE {
                    let i = cur as usize;
                    cur = next_active[i];
                    let i_score = score[i];
                    let mut h_sum = 0i32;
                    let mut h_xsum = 0.0f64;
                    let mut h_num = 0i16;
                    let mut h_link = NONE;
                    if i_score > c0 {
                        let prev = link0[i];
                        let keep = reuse_index0
                            && !first_pass
                            && (prev == NONE || !changed0[prev as usize]);
                        if keep {
                            // The previous choice was neither removed nor
                            // re-chosen: it is still the best choice (module
                            // documentation, "Index 0").
                            if prev != NONE {
                                let p = prev as usize;
                                h_num = num0[p];
                                h_sum = sum0[p];
                                h_xsum = xsum0[p];
                                stats.idx0_kept_link += 1;
                            } else {
                                stats.idx0_kept_none += 1;
                            }
                            h_link = prev;
                            changed0[i] = false;
                        } else {
                            changed0[i] = true;
                            stats.idx0_searched += 1;
                            let h_qe = qe[i];
                            let h_se = se[i];
                            let h_qe_gap = h_qe + w;
                            let h_se_gap = h_se + w;
                            let mut bestk = 0u64;
                            let cq_lo = (h_qe + 1).div_euclid(w);
                            let cq_hi = h_qe_gap.div_euclid(w);
                            let cs_lo = (h_se + 1).div_euclid(w);
                            let cs_hi = h_se_gap.div_euclid(w);
                            for cq in cq_lo..=cq_hi {
                                for cs in cs_lo..=cs_hi {
                                    let c = ((cq as u32 as u64) << 32) | (cs as u32 as u64);
                                    if let Some(&(a, b)) = grid_index.get(&c) {
                                        for &j in &grid_items[a as usize..b as usize] {
                                            let j = j as usize;
                                            if j >= i {
                                                break;
                                            }
                                            stats.idx0_scanned += 1;
                                            if !alive[j] {
                                                continue;
                                            }
                                            let qj = qo[j];
                                            let sj = so[j];
                                            // link_hsps.c:719-724
                                            if qj <= h_qe
                                                || sj <= h_se
                                                || qj > h_qe_gap
                                                || sj > h_se_gap
                                            {
                                                continue;
                                            }
                                            let s = sum0[j];
                                            if s > 0 {
                                                let k = key(s, j);
                                                if k > bestk {
                                                    bestk = k;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if bestk != 0 {
                                let j = key_idx(bestk);
                                h_num = num0[j];
                                h_sum = sum0[j];
                                h_xsum = xsum0[j];
                                h_link = j as u32;
                            }
                        }
                        if verify {
                            let expected = verify_scan_index0(
                                i,
                                &prev_active,
                                &qo,
                                &so,
                                &sum0,
                                &num0,
                                &xsum0,
                                qe[i],
                                se[i],
                            );
                            let actual = Choice {
                                link: h_link,
                                sum: h_sum,
                                num: h_num,
                                xsum_bits: h_xsum.to_bits(),
                            };
                            assert_eq!(
                                expected, actual,
                                "LOSAT_LINK_FAST_VERIFY: index 0, HSP {i} of {n}, kept={keep}"
                            );
                            stats.verified += 1;
                        }
                    }
                    // link_hsps.c:750-767
                    let new_sum = h_sum + (i_score - c0);
                    let new_xsum = h_xsum + sl[i] - lk[i];
                    sum0[i] = new_sum;
                    num0[i] = h_num + 1;
                    link0[i] = h_link;
                    xsum0[i] = new_xsum;
                    if h_link != NONE {
                        linked_to[h_link as usize] += 1;
                    }
                }
            }

            // ---------------- index 1 (link_hsps.c:771-896) ----------------
            for f in fen.iter_mut() {
                *f = 0;
            }
            let ilen = iorder.len();
            let mut p = 0usize;
            let mut qw = 0usize;
            for qi in 0..qorder.len() {
                let i = qorder[qi] as usize;
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
                        let k = key(sum1[j], j);
                        let mut x = rr[j] as usize;
                        while x <= m {
                            if fen[x] >= k {
                                break;
                            }
                            fen[x] = k;
                            x += x & x.wrapping_neg();
                        }
                    }
                }

                let mut h_sum = 0i32;
                let mut h_xsum = 0.0f64;
                let mut h_num = 0i16;
                let mut h_link = NONE;
                let prev = link1[i];
                // link_hsps.c:781-795
                let keep = !first_pass && (prev == NONE || !changed1[prev as usize]);
                if keep {
                    if prev != NONE {
                        let pj = prev as usize;
                        h_num = num1[pj];
                        h_sum = sum1[pj];
                        h_xsum = xsum1[pj];
                        stats.idx1_kept_link += 1;
                    } else {
                        stats.idx1_kept_none += 1;
                    }
                    h_link = prev;
                    changed1[i] = false;
                } else {
                    changed1[i] = true;
                    stats.idx1_searched += 1;
                    let mut bestk = 0u64;
                    let mut x = kq[i] as usize;
                    while x > 0 {
                        if fen[x] > bestk {
                            bestk = fen[x];
                        }
                        x &= x - 1;
                    }
                    if bestk != 0 && key_idx(bestk) > i {
                        // The tree answered with an HSP after H in list
                        // order (H spans at most 4 residues): NCBI's scan
                        // over the HSPs before H decides.
                        stats.fallbacks += 1;
                        let h_se = se[i];
                        bestk = 0;
                        let mut j = prev_active[i];
                        while j != NONE {
                            let jj = j as usize;
                            if score[jj] > c1 && qo[jj] > h_qe && so[jj] > h_se {
                                let k = key(sum1[jj], jj);
                                if k > bestk {
                                    bestk = k;
                                }
                            }
                            j = prev_active[jj];
                        }
                    }
                    if bestk != 0 {
                        let j = key_idx(bestk);
                        h_num = num1[j];
                        h_sum = sum1[j];
                        h_xsum = xsum1[j];
                        h_link = j as u32;
                    }
                }
                if verify {
                    // link_hsps.c:812-823: a search starts one below the sum
                    // of the previous choice when that HSP remains. A kept
                    // choice is compared with a search from 0.
                    let start_sum =
                        if !keep && !first_pass && prev != NONE && linked_to[prev as usize] >= 0 {
                            sum1[prev as usize] - 1
                        } else {
                            0
                        };
                    let expected = verify_scan_index1(
                        i,
                        &prev_active,
                        &qo,
                        &so,
                        &sum1,
                        &num1,
                        &xsum1,
                        h_qe,
                        se[i],
                        start_sum,
                    );
                    let actual = Choice {
                        link: h_link,
                        sum: h_sum,
                        num: h_num,
                        xsum_bits: h_xsum.to_bits(),
                    };
                    assert_eq!(
                        expected, actual,
                        "LOSAT_LINK_FAST_VERIFY: index 1, HSP {i} of {n}, kept={keep}"
                    );
                    stats.verified += 1;
                }
                // link_hsps.c:866-894
                let new_sum = h_sum + (score[i] - c1);
                let new_xsum = h_xsum + sl[i] - lk[i];
                sum1[i] = new_sum;
                num1[i] = h_num + 1;
                link1[i] = h_link;
                xsum1[i] = new_xsum;
                if h_link != NONE {
                    linked_to[h_link as usize] += 1;
                }
            }
            qorder.truncate(qw);
            iorder.retain(|&j| alive[j as usize]);

            // ---- maxima (link_hsps.c:759-763, 887-891) and the tree ----
            if first_pass {
                for i in 0..n {
                    tree.nodes[tree.cap + i] = [key(sum0[i], i), key(sum1[i], i)];
                }
                tree.build();
            } else {
                let mut cur = active_head;
                while cur != NONE {
                    let i = cur as usize;
                    cur = next_active[i];
                    let v = [key(sum0[i], i), key(sum1[i], i)];
                    if tree.leaf(i) != v {
                        tree.set(i, v);
                    }
                }
            }
            let root = tree.root();
            best = [None, None];
            for ch in usize::from(ignore_small_gaps)..2 {
                if root[ch] != 0 {
                    best[ch] = Some(key_idx(root[ch]));
                }
            }

            path_changed = false;
            first_pass = false;
        }

        // ---- ordering method (link_hsps.c:901-952) ----
        let mut prob = [f64::MAX, f64::MAX];
        if !ignore_small_gaps {
            if let Some(bi) = best[0] {
                sum0[bi] += (num0[bi] as i32) * c0;
                tree.set(bi, [key(sum0[bi], bi), key(sum1[bi], bi)]);
                let num = num0[bi] as usize;
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[0] = small_gap_sum_e(
                    WINDOW_SIZE,
                    num as i16,
                    xsum0[bi],
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
                let num = num1[bi] as usize;
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[1] = ncbi_large_gap_sum_e(
                    num as i16,
                    xsum1[bi],
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
            sum1[bi] += (num1[bi] as i32) * c1;
            tree.set(bi, [key(sum0[bi], bi), key(sum1[bi], bi)]);
            let num = num1[bi] as usize;
            let divisor = gap_decay_divisor(gap_decay_rate, num);
            prob[1] = ncbi_large_gap_sum_e(
                num as i16,
                xsum1[bi],
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
        let best_i = if let Some(i) = best[ordering] {
            i
        } else if let Some(i) = best[1 - ordering] {
            i
        } else {
            break;
        };
        let evalue = prob[ordering];

        // ---- remove the chain (link_hsps.c:955-980) ----
        if linked_to[best_i] > 0 {
            path_changed = true;
        }
        let head_link = if ordering == 0 {
            link0[best_i]
        } else {
            link1[best_i]
        };
        let linked_set = head_link != NONE;
        group_hits[best_i].start_of_chain = true;
        group_hits[best_i].hsp_link_num = if ordering == 0 {
            num0[best_i]
        } else {
            num1[best_i]
        };

        let mut is_first = true;
        let mut cur = best_i;
        loop {
            if linked_to[cur] > 1 {
                path_changed = true;
            }
            linked_to[cur] = -1000;
            changed0[cur] = true;
            changed1[cur] = true;
            alive[cur] = false;

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
            // Keep prev_active[cur] / next_active[cur]: they are not read again.
            tree.set(cur, [0, 0]);

            group_hits[cur].ordering_method = ordering as u8;
            group_hits[cur].e_value = evalue;
            group_hits[cur].linked_set = linked_set;
            let next = if ordering == 0 {
                link0[cur]
            } else {
                link1[cur]
            };
            group_hits[cur].chain_next_link_id = if next == NONE {
                None
            } else {
                Some(group_hits[next as usize].link_id)
            };
            if is_first {
                is_first = false;
            } else {
                group_hits[cur].start_of_chain = false;
            }
            remaining -= 1;
            if next == NONE {
                break;
            }
            cur = next as usize;
        }
    }

    (group_hits, stats)
}

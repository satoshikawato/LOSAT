//! Core HSP linking algorithm
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c

use crate::stats::sum_statistics::{
    defaults::{GAP_SIZE, OVERLAP_SIZE},
    gap_decay_divisor, ncbi_large_gap_sum_e, normalize_score, small_gap_sum_e,
};
use crate::stats::KarlinParams;
use crate::utils::glibc_fma_libm::GlibcLibm;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::cell::RefCell;
use std::cmp::Ordering;
use std::sync::OnceLock;
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1297
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::algorithm::tblastx::chaining::UngappedHit;
use crate::algorithm::tblastx::diagnostics::diagnostics_enabled;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::tracing::int4_nucleotide_coords;
use crate::algorithm::tblastx::ungapped_hit_sort::sort_ungapped_hits_by_index_replay;

use super::params::LinkingParams;

// The HSP list sorted, linked as sub-slices and replayed without a second whole copy alive; a
// child module of this file so that it uses the comparators, per-group functions and replay below
// unchanged (see tblastx/x_hsp_inplace.rs).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-480
// NCBI works on one array of pointers to the HSPs.
#[path = "../x_hsp_inplace.rs"]
pub(crate) mod x_hsp_inplace;
pub use x_hsp_inplace::print_summary as x_hsp_inplace_print_summary;

// The index sorts of
// `x_hsp_inplace` as one final sort instead of two and, for long lists, as stable parallel sorts;
// a child module of this file so that it names the comparators below unchanged (see
// tblastx/x_link_sort.rs).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// NCBI's qsort is not stable; the module returns the order of LOSAT's stable sorts.
#[path = "../x_link_sort.rs"]
mod x_link_sort;
pub use x_link_sort::print_summary as x_link_sort_print_summary;

// The groups of the linking call of `x_hsp_inplace` linked as disjoint sub-slices of the one
// sorted list instead of Vecs cut off it and put together again; a child module of this file so
// that it calls the per-group functions below unchanged (see tblastx/x_link_slices.rs).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// NCBI's groups are sublists of the one sorted array, linked in place.
#[path = "../x_link_slices.rs"]
mod x_link_slices;
pub use x_link_slices::print_summary as x_link_slices_print_summary;
pub(crate) use x_link_slices::{x_helped as x_link_helped, XHelped};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
const TRIM_SIZE: i32 = (OVERLAP_SIZE + 1) / 2;

/// Sentinel index indicating end of active list or "not in list"
const SENTINEL_IDX: usize = usize::MAX;

// ---------------------------------------------------------------------------
// HSP tracing (opt-in via env var)
// ---------------------------------------------------------------------------

/// Trace a single HSP through linking by matching final outfmt coordinates
/// (qstart,qend,sstart,send).
///
/// Enable by setting:
/// - `LOSAT_TRACE_HSP="qstart,qend,sstart,send"`
/// - `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` for chain-selection
///   summaries without the per-candidate DP trace.
#[derive(Clone, Copy, Debug)]
struct TraceHspTarget {
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
}

static TRACE_HSP_TARGET: OnceLock<Option<TraceHspTarget>> = OnceLock::new();
static TRACE_CHAIN_TARGET: OnceLock<Option<TraceHspTarget>> = OnceLock::new();
static TRACE_LINK_SELECTIONS: OnceLock<bool> = OnceLock::new();

#[inline(always)]
fn trace_hsp_target() -> Option<TraceHspTarget> {
    if !crate::utils::xdiag::ON {
        return None;
    }
    *TRACE_HSP_TARGET.get_or_init(|| parse_trace_hsp_target_env("LOSAT_TRACE_HSP"))
}

#[inline(always)]
fn trace_chain_target() -> Option<TraceHspTarget> {
    if !crate::utils::xdiag::ON {
        return None;
    }
    *TRACE_CHAIN_TARGET.get_or_init(|| parse_trace_hsp_target_env("LOSAT_TRACE_CHAIN_HSP"))
}

#[inline(always)]
fn trace_link_selection_summaries() -> bool {
    if !crate::utils::xdiag::ON {
        return false;
    }
    *TRACE_LINK_SELECTIONS
        .get_or_init(|| crate::utils::xdiag::var_os("LOSAT_TRACE_LINK_SELECTIONS").is_some())
}

/// Use the literal O(n^2)-per-pass port (`link_hsp_group_ncbi`) instead of the
/// result-identical incremental linker: requested explicitly, or when any linking
/// trace/debug output is enabled (only the literal port prints those traces).
fn legacy_linking_requested() -> bool {
    static LEGACY: OnceLock<bool> = OnceLock::new();
    *LEGACY.get_or_init(|| {
        std::env::var_os("LOSAT_LINKING_LEGACY").is_some()
            || crate::utils::xdiag::var_os("LOSAT_TRACE_HSP").is_some()
            || crate::utils::xdiag::var_os("LOSAT_TRACE_CHAIN_HSP").is_some()
            || crate::utils::xdiag::var_os("LOSAT_TRACE_LINK_SELECTIONS").is_some()
            || crate::utils::xdiag::var_os("LOSAT_DEBUG_CHAINING").is_some()
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419
// No NCBI counterpart for the switch: with it set, the literal port of s_BlastEvenGapLinkHSPs
// also links every group and the complete results are compared. It does not change any value NCBI
// computes.
/// `LOSAT_LINK_FAST_SHADOW=1`: also link every group with the literal port and
/// compare the complete results (`link_hsp_group_incremental`).
fn link_shadow_enabled() -> bool {
    static FLAG: OnceLock<bool> = OnceLock::new();
    *FLAG.get_or_init(|| std::env::var_os("LOSAT_LINK_FAST_SHADOW").is_some())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:901-982
// The parser is diagnostic-only; it selects the NCBI chain-removal boundary to
// print and does not alter the linking algorithm or runtime output.
#[inline(always)]
fn parse_trace_hsp_target_env(var_name: &str) -> Option<TraceHspTarget> {
    let raw = crate::utils::xdiag::var(var_name).ok()?;
    let parts: Vec<&str> = raw
        .split(|c: char| c == ',' || c == ':' || c == ';' || c.is_whitespace())
        .filter(|p| !p.is_empty())
        .collect();
    if parts.len() != 4 {
        eprintln!(
            "[TRACE_HSP] invalid {} (expected 4 integers): {:?}",
            var_name, raw
        );
        return None;
    }
    let q_start: i32 = parts[0].parse().ok()?;
    let q_end: i32 = parts[1].parse().ok()?;
    let s_start: i32 = parts[2].parse().ok()?;
    let s_end: i32 = parts[3].parse().ok()?;
    Some(TraceHspTarget {
        q_start,
        q_end,
        s_start,
        s_end,
    })
}

// No NCBI counterpart (diagnostics): whether an HSP's nucleotide coordinates, its ends read as
// NCBI's `Int4` (`tracing::int4_nucleotide_coords`), are the trace target's.
#[inline(always)]
fn trace_hit_matches_target(target: TraceHspTarget, h: &UngappedHit) -> bool {
    int4_nucleotide_coords(h)
        == (
            i64::from(target.q_start),
            i64::from(target.q_end),
            i64::from(target.s_start),
            i64::from(target.s_end),
        )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:343-349
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:76-77
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1373-1376
// Every query has six contexts in NCBI, also a frame without a residue, so
// `context/(NUM_FRAMES/2)` is `2 * query + (frame < 0)`. LOSAT's `ctx_idx` counts only the
// frames that have a residue (none for 1-2 nt, two for 3 nt, four for 4 nt), so it would
// shift the groups of the queries after such a query in the batch.
#[inline]
fn translated_context_group(hit: &UngappedHit) -> usize {
    hit.q_idx as usize * 2 + usize::from(hit.q_frame < 0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:330-376
fn rev_compare_hsps_tbx(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    translated_context_group(a)
        .cmp(&translated_context_group(b))
        .then_with(|| a.s_frame.signum().cmp(&b.s_frame.signum()))
        .then_with(|| b.q_aa_start.cmp(&a.q_aa_start))
        .then_with(|| b.q_aa_end.cmp(&a.q_aa_end))
        .then_with(|| b.s_aa_start.cmp(&a.s_aa_start))
        .then_with(|| b.s_aa_end.cmp(&a.s_aa_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:359-374
// `s_RevCompareHSPsTbx` with the ends read as NCBI's `Int4` (`UngappedHit::q_end_int4`), for a
// list of the one-hit word finder that holds an HSP of length 0 or less
// (`link_one_hit_with_int4_ends`). On ends of 0 or more it is `rev_compare_hsps_tbx`.
fn rev_compare_hsps_tbx_int4(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    translated_context_group(a)
        .cmp(&translated_context_group(b))
        .then_with(|| a.s_frame.signum().cmp(&b.s_frame.signum()))
        .then_with(|| b.q_aa_start.cmp(&a.q_aa_start))
        .then_with(|| b.q_end_int4().cmp(&a.q_end_int4()))
        .then_with(|| b.s_aa_start.cmp(&a.s_aa_start))
        .then_with(|| b.s_end_int4().cmp(&a.s_end_int4()))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:233-264
fn rev_compare_hsps_transl(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    translated_context_group(a)
        .cmp(&translated_context_group(b))
        .then_with(|| b.q_aa_start.cmp(&a.q_aa_start))
        .then_with(|| b.s_aa_start.cmp(&a.s_aa_start))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:157-186
fn fwd_compare_hsps_transl(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    translated_context_group(a)
        .cmp(&translated_context_group(b))
        .then_with(|| a.q_aa_start.cmp(&b.q_aa_start))
        .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-486
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
#[inline]
fn sort_hsps_by_ncbi_link_order<C>(hits: &mut [UngappedHit], compare: C)
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
{
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    // NCBI's translated-query comparators are partial and NCBI delegates
    // comparator-equal rows to the platform qsort implementation. Do not add
    // frame or LOSAT insertion-order tie-breakers here.
    sort_ungapped_hits_by_index_replay(hits, compare);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// Keep the completed groups' payloads stationary during both final sorts.
#[inline]
fn sort_linked_hsp_references(
    hits: &mut [&mut UngappedHit],
    compare: fn(&UngappedHit, &UngappedHit) -> Ordering,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    // NCBI's translated-query comparators are partial and NCBI delegates
    // comparator-equal rows to the platform qsort implementation. Do not add
    // frame or LOSAT insertion-order tie-breakers here.
    // Keep every payload in the caller's fixed Vec. Stable ordering retains
    // the same comparator-equal order as the previous index-replay sort.
    hits.sort_by(|lhs, rhs| compare(lhs, rhs));
}

// ---------------------------------------------------------------------------
// Core structures
// ---------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:104-105
// NCBI keeps the sums in Int4; LOSAT keeps them in i64 and does not reproduce NCBI's wrap beyond
// Int4 (Owner decision 2026-10-10; linking_incr.rs, "Sums").
/// lh_helper structure (NCBI link_hsps.c:100-109)
/// This is the "hot" data used in the inner loop
#[derive(Clone)]
struct LhHelper {
    hsp_idx: usize, // Original HSP index (NCBI: ptr field)
    q_off_trim: i32,
    s_off_trim: i32,
    sum: [i64; 2],      // sum for both indices (small gap, large gap)
    next_larger: usize, // index of next HSP with larger sum[1]
    maxsum1: i64,       // NCBI threshold for stopping link attempts (disabled in NCBI source)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:608-633,654-658:
// Update each channel's maximum (`max0`, `max1`) and winning HSP (`best`) whenever a
// candidate's sum (`sum0`, `sum1`) is at least the current maximum so later ties win.
// Each channel retains its own (score, stable post-sort group index) winner.
// SENTINEL_IDX means inactive, never a real contender at cutoff equality.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct DualMaximum {
    sum: [i64; 2],
    index: [usize; 2],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:610-623:
// Each channel updates its running maximum and winning HSP on greater-than-or-equal sums.
impl DualMaximum {
    const INACTIVE: Self = Self {
        sum: [0; 2],
        index: [SENTINEL_IDX; 2],
    };

    #[inline]
    fn active(index: usize, sum: [i64; 2]) -> Self {
        Self {
            sum,
            index: [index; 2],
        }
    }

    #[inline]
    fn reduce(left: Self, right: Self) -> Self {
        let mut result = left;
        for channel in 0..2 {
            if right.index[channel] != SENTINEL_IDX
                && (left.index[channel] == SENTINEL_IDX
                    || right.sum[channel] > left.sum[channel]
                    || (right.sum[channel] == left.sum[channel]
                        && right.index[channel] > left.index[channel]))
            {
                result.sum[channel] = right.sum[channel];
                result.index[channel] = right.index[channel];
            }
        }
        result
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:603-659,897-908,942-978:
// After the first pass, locate the current maximum stored sums; when the cached maximum cannot
// be reused, recompute stored sums, and before scoring add back `num` times the small-gap cutoff
// to the small-gap winner's stored sum.
// A group-local flat tree mirrors stored sums and physical membership only.
// It does not validate paths or recompute survivors after predecessor removal.
struct DualMaximumTree {
    capacity: usize,
    nodes: Vec<DualMaximum>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:608-633,755-758,869-873,975-978:
// Maintain the running maximum and winning HSP for both gap channels using non-strict ties.
impl DualMaximumTree {
    fn new(group_len: usize) -> Self {
        let capacity = group_len.max(1).next_power_of_two();
        Self {
            capacity,
            nodes: vec![DualMaximum::INACTIVE; capacity * 2],
        }
    }

    fn rebuild(&mut self, links: &[HspLink], mut active: usize) {
        self.nodes[self.capacity..].fill(DualMaximum::INACTIVE);
        while active != SENTINEL_IDX {
            self.nodes[self.capacity + active] = DualMaximum::active(active, links[active].sum);
            active = links[active].next_active;
        }
        for parent in (1..self.capacity).rev() {
            self.nodes[parent] =
                DualMaximum::reduce(self.nodes[parent * 2], self.nodes[parent * 2 + 1]);
        }
    }

    #[inline]
    fn update(&mut self, index: usize, leaf: DualMaximum) {
        let mut position = self.capacity + index;
        self.nodes[position] = leaf;
        while position > 1 {
            position /= 2;
            self.nodes[position] =
                DualMaximum::reduce(self.nodes[position * 2], self.nodes[position * 2 + 1]);
        }
    }

    #[inline]
    fn maxima(&self, cutoffs: [i32; 2], ignore_small_gaps: bool) -> ([Option<usize>; 2], [i64; 2]) {
        let root = self.nodes[1];
        let mut best = [None; 2];
        let mut sums = [-i64::from(cutoffs[0]), -i64::from(cutoffs[1])];
        for channel in usize::from(ignore_small_gaps)..2 {
            if root.index[channel] != SENTINEL_IDX && root.sum[channel] >= sums[channel] {
                best[channel] = Some(root.index[channel]);
                sums[channel] = root.sum[channel];
            }
        }
        (best, sums)
    }
}

/// HSP link information (NCBI BlastHSPLink structure)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:63-71
#[derive(Clone)]
struct HspLink {
    score: i32,
    ctx_idx: usize, // Query context index (`H->hsp->context`, used for `kbp` lookup)
    q_off_trim: i32,
    s_off_trim: i32,
    q_end_trim: i32,
    s_end_trim: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:67
    // i64, not Int4: NCBI's wrap beyond Int4 is not reproduced (Owner decision 2026-10-10;
    // linking_incr.rs, "Sums").
    sum: [i64; 2],  // sum for both indices
    xsum: [f64; 2], // normalized sum for both indices
    num: [i16; 2],  // number of HSPs in chain for both indices
    // Best previous HSP to link with for both indices (NCBI: LinkHSPStruct* or NULL).
    // We store an index into `hsp_links`, using SENTINEL_IDX as the NULL marker.
    link: [usize; 2],
    changed: bool,        // Has link been changed? (NCBI: hsp_link.changed)
    linked_to: i32,       // How many HSPs link TO this one (NCBI: linked_to)
    start_of_chain: bool, // Is this HSP the start of a chain? (for output filtering)
    linked_set: bool,     // Is this HSP part of a multi-HSP chain? (NCBI: linked_set)
    // Intrusive linked list for active HSPs (NCBI uses hp_start->next linked list)
    // SENTINEL_IDX indicates end of list or not in list
    next_active: usize, // Next active HSP in list
    prev_active: usize, // Previous active HSP in list
}

/// Thread-local buffer pools for memory reuse across frame groups
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:452-454
/// Allocates `lh_helper` once with size `max(1024, hspcnt + 5)`
/// and reuses it across all frame groups (lines 553-983)
struct BufferPools {
    lh_helpers: Vec<LhHelper>,
    hsp_links: Vec<HspLink>,
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
pub fn apply_sum_stats_even_gap_linking_with_parallel(
    mut hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
) -> Vec<UngappedHit> {
    if hits.is_empty() {
        return hits;
    }

    // Dispatch: the steps of this function on one list (sort in place, groups linked as
    // sub-slices of the list, output replay in place) with the same per-group calls; see
    // x_hsp_inplace.rs for why it is the same list. With LOSAT_DEBUG_CUTOFFS or
    // LOSAT_DIAGNOSTICS, the body below runs on the untouched input.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486,990-994
    // NCBI sorts and regroups one array of pointers and never copies the HSPs.
    match x_hsp_inplace::linking_dispatch(
        hits,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
    ) {
        Ok(linked) => return linked,
        Err(untouched) => hits = untouched,
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-480
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1773-1776
    for (link_id, hit) in hits.iter_mut().enumerate() {
        hit.link_id = u32::try_from(link_id)
            .expect("link_id: a position in the linking list fits u32 (NCBI's hspcnt is Int4)");
        hit.chain_next_link_id = None;
        hit.num = 1;
    }

    let diag_enabled = diagnostics_enabled();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:468-493
    let cutoffs = linking_params.cutoffs;

    // Debug output for cutoffs (controlled by LOSAT_DEBUG_CUTOFFS env var)
    static DEBUG_CUTOFFS_PRINTED: std::sync::atomic::AtomicBool =
        std::sync::atomic::AtomicBool::new(false);
    if crate::utils::xenv::debug_cutoffs_is_ok()
        && !DEBUG_CUTOFFS_PRINTED.swap(true, std::sync::atomic::Ordering::SeqCst)
    {
        eprintln!("=== LOSAT Linking Cutoffs Debug ===");
        eprintln!("  subject_len_nucl: {}", linking_params.subject_len_nucl);
        eprintln!("  lambda: {:.6}", params.lambda);
        eprintln!("  K: {:.6}", params.k);
        eprintln!("  H: {:.6}", params.h);
        eprintln!("  cutoff_small_gap: {}", cutoffs.cutoff_small_gap);
        eprintln!("  cutoff_big_gap: {}", cutoffs.cutoff_big_gap);
        eprintln!("  gap_prob: {:.6}", cutoffs.gap_prob);
        eprintln!("  ignore_small_gaps: {}", cutoffs.ignore_small_gaps());
        eprintln!("  WINDOW_SIZE: {}", WINDOW_SIZE);
        eprintln!("  TRIM_SIZE: {}", TRIM_SIZE);
        eprintln!("===================================");
    }

    // =======================================================================
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c EXACT REPLICATION:
    // 1. Sort ALL HSPs by s_RevCompareHSPsTbx (lines 484-490)
    // 2. Scan sorted list and detect frame boundaries (lines 514-534)
    // 3. Process each frame group SEQUENTIALLY (lines 553-982)
    // =======================================================================

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2722-2725
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:750-752,865-867
    // Keep per-context logK available for the exact NCBI xsum expression below.
    let log_k_by_ctx: Vec<f64> = query_contexts
        .iter()
        .map(|ctx| ctx.karlin_params.k.glibc_ln())
        .collect();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-486
    sort_hsps_by_ncbi_link_order(&mut hits, rev_compare_hsps_tbx);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:510-531
    let total_hits = hits.len(); // Save for pool sizing (NCBI `hsp_list->hspcnt`)
    let mut frame_groups: Vec<Vec<UngappedHit>> = Vec::new();
    let mut current_group: Vec<UngappedHit> = Vec::new();

    for hit in hits {
        let ctx_group = translated_context_group(&hit);
        let s_sign = hit.s_frame.signum();

        if let Some(prev) = current_group.last() {
            let prev_ctx_group = translated_context_group(prev);
            let prev_s_sign = prev.s_frame.signum();

            if prev_ctx_group != ctx_group || prev_s_sign != s_sign {
                // Frame switch - start new group
                if !current_group.is_empty() {
                    frame_groups.push(std::mem::take(&mut current_group));
                }
            }
        }
        current_group.push(hit);
    }
    if !current_group.is_empty() {
        frame_groups.push(current_group);
    }

    // Step 3: Process each frame group (parallel when enabled, sequential otherwise).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
    // NCBI processes groups in frame_index order (0, 1, 2, ...), which matches
    // the order of frame_groups in LOSAT (created from sorted HSP list).
    //
    // NCBI memory pooling: c++/src/algo/blast/core/link_hsps.c:452-454 allocates `lh_helper` once with
    // size `max(1024, hspcnt + 5)` and reuses it across all frame groups.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
    crate::utils::threading::report_stage(
        "tblastx",
        "linking",
        frame_groups.len(),
        allow_parallel
            && frame_groups.len() > 1
            && cfg!(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
    );
    let initial_lh_size = (total_hits + 5).max(1024); // NCBI: the larger of 1024 and `hspcnt + 5`

    // The incremental kernel of linking_incr.rs links every group with the same result as the
    // literal port `link_hsp_group_ncbi`, which stays for LOSAT_LINKING_LEGACY, the linking traces
    // and as the reference of the tests.
    let use_legacy = legacy_linking_requested();

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let mut results: Vec<UngappedHit> = if allow_parallel && frame_groups.len() > 1 {
        // Use thread_local! macro for thread-local storage (rayon doesn't have ThreadLocal)
        thread_local! {
            static POOLS: RefCell<Option<BufferPools>> = const { RefCell::new(None) };
        }

        let mut indexed_results: Vec<(usize, Vec<UngappedHit>)> = frame_groups
            .into_par_iter()
            .enumerate()
            .map(|(group_idx, group_hits)| {
                POOLS.with(|pools| {
                    let mut pools_opt = pools.borrow_mut();
                    let pool = pools_opt.get_or_insert_with(|| BufferPools {
                        lh_helpers: Vec::with_capacity(initial_lh_size),
                        hsp_links: Vec::with_capacity(total_hits),
                    });
                    let processed = if use_legacy {
                        link_hsp_group_ncbi(
                            group_hits,
                            params,
                            &cutoffs,
                            linking_params.gap_decay_rate,
                            diag_enabled,
                            linking_params.subject_len_nucl,
                            query_contexts,
                            subject_frame_bases,
                            length_adj_per_context,
                            eff_searchsp_per_context,
                            &log_k_by_ctx,
                            &mut pool.lh_helpers,
                            &mut pool.hsp_links,
                        )
                    } else {
                        link_hsp_group_incremental(
                            group_hits,
                            params,
                            &cutoffs,
                            linking_params.gap_decay_rate,
                            linking_params.subject_len_nucl,
                            query_contexts,
                            subject_frame_bases,
                            length_adj_per_context,
                            eff_searchsp_per_context,
                            &log_k_by_ctx,
                            &mut pool.lh_helpers,
                            &mut pool.hsp_links,
                        )
                    };
                    (group_idx, processed)
                })
            })
            .collect();

        // Sort by group index to preserve NCBI processing order
        indexed_results.sort_by_key(|(idx, _)| *idx);

        indexed_results
            .into_iter()
            .flat_map(|(_, hits)| hits)
            .collect()
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
        let mut pools = BufferPools {
            lh_helpers: Vec::with_capacity(initial_lh_size),
            hsp_links: Vec::with_capacity(total_hits),
        };
        let mut ordered: Vec<UngappedHit> = Vec::new();
        for group_hits in frame_groups.into_iter() {
            let processed = if use_legacy {
                link_hsp_group_ncbi(
                    group_hits,
                    params,
                    &cutoffs,
                    linking_params.gap_decay_rate,
                    diag_enabled,
                    linking_params.subject_len_nucl,
                    query_contexts,
                    subject_frame_bases,
                    length_adj_per_context,
                    eff_searchsp_per_context,
                    &log_k_by_ctx,
                    &mut pools.lh_helpers,
                    &mut pools.hsp_links,
                )
            } else {
                link_hsp_group_incremental(
                    group_hits,
                    params,
                    &cutoffs,
                    linking_params.gap_decay_rate,
                    linking_params.subject_len_nucl,
                    query_contexts,
                    subject_frame_bases,
                    length_adj_per_context,
                    eff_searchsp_per_context,
                    &log_k_by_ctx,
                    &mut pools.lh_helpers,
                    &mut pools.hsp_links,
                )
            };
            ordered.extend(processed);
        }
        ordered
    };

    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let mut results: Vec<UngappedHit> = {
        let mut pools = BufferPools {
            lh_helpers: Vec::with_capacity(initial_lh_size),
            hsp_links: Vec::with_capacity(total_hits),
        };
        let mut ordered: Vec<UngappedHit> = Vec::new();
        for group_hits in frame_groups.into_iter() {
            let processed = if use_legacy {
                link_hsp_group_ncbi(
                    group_hits,
                    params,
                    &cutoffs,
                    linking_params.gap_decay_rate,
                    diag_enabled,
                    linking_params.subject_len_nucl,
                    query_contexts,
                    subject_frame_bases,
                    length_adj_per_context,
                    eff_searchsp_per_context,
                    &log_k_by_ctx,
                    &mut pools.lh_helpers,
                    &mut pools.hsp_links,
                )
            } else {
                link_hsp_group_incremental(
                    group_hits,
                    params,
                    &cutoffs,
                    linking_params.gap_decay_rate,
                    linking_params.subject_len_nucl,
                    query_contexts,
                    subject_frame_bases,
                    length_adj_per_context,
                    eff_searchsp_per_context,
                    &log_k_by_ctx,
                    &mut pools.lh_helpers,
                    &mut pools.hsp_links,
                )
            };
            ordered.extend(processed);
        }
        ordered
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994,1080-1085
    // Preserve the existing empty/singleton replay shortcut. The owned group
    // kernel and its contiguous payload layout are unchanged; only these
    // adjacent final sorts and replay borrow stationary result records.
    if results.len() <= 1 {
        return results;
    }
    let mut references: Vec<&mut UngappedHit> = results.iter_mut().collect();
    sort_linked_hsp_references(&mut references, rev_compare_hsps_transl);
    sort_linked_hsp_references(&mut references, fwd_compare_hsps_transl);

    replay_ncbi_output_list(references)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873 (the first
// `BLAST_LinkHsps` of `s_BlastSearchEngineCore`) and c++/src/algo/blast/core/blast_hits.c:2703-2731
// (`Blast_HSPListReevaluateUngapped`, called at blast_engine.c:1492-1497)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:687,699 (in
// `Blast_HSPReevaluateWithAmbiguitiesUngapped`) and :452-453 (in `s_UpdateReevaluatedHSP`)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1480-1497 and
// c++/src/algo/blast/core/blast_hits.c:2651-2665 (the re-evaluation's fetch of the subject)
// The first linking call of a subject whose list from the one-hit word finder (`-window_size 0`)
// holds an HSP of length 0 or less (`UngappedHit::is_int4_inverted`), with its ends as NCBI's
// `Int4` (`chaining.rs`, the impl of `UngappedHit`). The list is first sorted by score with the
// `Int4` comparator: NCBI's frame lists are combined by that sort (blast_hits.c:2758-2766), which
// LOSAT's `append_frame_hits` does with the `u32` one. The two comparators order an end below 0
// differently but have the same ties (two `i32` are equal exactly when their `u32` bits are), so
// the stable `u32` sorts of the appends keep every tie class of the `Int4` comparator in the
// order of the concatenated frame lists, and this stable sort then gives NCBI's list (the stable
// sort of that concatenation; `x_link_sort.rs`, (S)). The link is then `link_in_place`'s with
// `rev_compare_hsps_tbx_int4` and the literal port `link_hsp_group_ncbi` for every group
// (`x_hsp_inplace::link_in_place_int4`): the incremental kernel assumes HSPs whose end is after
// their start (`linking_incr.rs`), while the literal port computes with `Int4` coordinates as NCBI
// does. NCBI then sorts by score, reaps by the preliminary E-value and re-evaluates, whose fetch
// of the subject comes first for every list that still holds an HSP, and which scores an HSP of length 0 or less 0 (its loop does not run),
// below every context's cutoff (at least 1), and deletes it. So those HSPs stay in the list with
// the E-values this call gives them, through the sort and the reap, up to that fetch; they leave
// after it (`subject_link::drop_unread_int4_inverted`). This path runs only for such a subject;
// the same subject without such an HSP links with `apply_sum_stats_even_gap_linking_with_parallel`.
/// The first linking call of a one-hit list holding an HSP of length 0 or less, which keeps those
/// HSPs.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn link_one_hit_with_int4_ends(
    mut hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    allow_parallel: bool,
) -> Vec<UngappedHit> {
    if hits.is_empty() {
        return hits;
    }
    crate::algorithm::tblastx::blast_gapalign::sort_ungapped_hits_by_score_int4(&mut hits);
    x_hsp_inplace::link_in_place_int4(
        hits,
        params,
        linking_params,
        query_contexts,
        subject_frame_bases,
        length_adj_per_context,
        eff_searchsp_per_context,
        allow_parallel,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1012-1028
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1059
fn replay_ncbi_output_list(mut sorted_hits: Vec<&mut UngappedHit>) -> Vec<UngappedHit> {
    if sorted_hits.len() <= 1 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1085
        // Retain the existing singleton path, including a non-head singleton.
        // Materialize once at this Rust ownership boundary.
        return sorted_hits.into_iter().map(|hit| hit.clone()).collect();
    }

    let max_link_id = sorted_hits
        .iter()
        .map(|h| h.link_id as usize)
        .max()
        .unwrap_or(0);
    let mut id_to_index = vec![SENTINEL_IDX; max_link_id + 1];
    for (index, hit) in sorted_hits.iter().enumerate() {
        if (hit.link_id as usize) < id_to_index.len() {
            id_to_index[hit.link_id as usize] = index;
        }
    }

    let mut replayed_indices = Vec::with_capacity(sorted_hits.len());
    for index in 0..sorted_hits.len() {
        let hit = &sorted_hits[index];
        if hit.linked_set && !hit.start_of_chain {
            continue;
        }

        let mut current = index;
        let mut hops = 0usize;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1061
        // Every HSP of a chain gets the head's count; an HSP without a link keeps 1.
        let num_links = sorted_hits[index]
            .chain_next_link_id
            .map(|_| i32::from(sorted_hits[index].hsp_link_num));
        loop {
            replayed_indices.push(current);
            if let Some(num_links) = num_links {
                sorted_hits[current].num = num_links;
            }
            let Some(next_id) = sorted_hits[current].chain_next_link_id else {
                break;
            };
            let next_id = next_id as usize;
            if next_id >= id_to_index.len() {
                break;
            }
            let next_index = id_to_index[next_id];
            if next_index == SENTINEL_IDX {
                break;
            }
            current = next_index;
            hops += 1;
            if hops > sorted_hits.len() {
                break;
            }
        }
    }

    replayed_indices
        .into_iter()
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1080-1085
        // Preserve the complete replay, including all chain members. Only
        // this final boundary copies payloads back into an owned result Vec.
        .map(|index| (*sorted_hits[index]).clone())
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-863
// When the score guard `b0` triggers, jump the predecessor index to `next_larger`; otherwise
// adopt the predecessor's `sum`, `xsum`, `num`, and link when all geometry predicates pass.
// Both cursor drivers expand this single selection/trace body. All caller
// state is explicit; continue targets the surrounding predecessor loop.
macro_rules! visit_large_gap_predecessor {
    (@advance indexed, $cursor:ident, $next_cursor:expr) => { $cursor -= 1; };
    (@advance prefix, $cursor:ident, $next_cursor:expr) => { $cursor = $next_cursor; };
    (@jump indexed, $cursor:ident, $next:expr) => { $cursor = $next; };
    (@jump prefix, $cursor:ident, $next:expr) => {
        $cursor = &$cursor[..$next.saturating_sub(1)];
    };
    ($mode:ident, $current_idx:ident, $helper:ident, $cursor:ident, $next_cursor:expr;
     $h_sum:ident, $h_num:ident, $h_xsum:ident, $h_link:ident,
     $is_target_hsp:ident, $h_qe:ident, $h_se:ident,
     $pool_hsp_links:ident, $group_hits:ident) => {{
        let sum = $helper.sum[1];
        let next_larger = $helper.next_larger;
        let b0 = sum <= $h_sum;

        // Step backward to the previous predecessor index (link_hsps.c:841).
        visit_large_gap_predecessor!(@advance $mode, $cursor, $next_cursor);
        if b0 {
            visit_large_gap_predecessor!(@jump $mode, $cursor, next_larger);
            if $is_target_hsp {
                eprintln!(
                    "  [SKIP] current_idx={}, sum={} <= h_sum={}, jump to next_larger={}",
                    $current_idx, sum, $h_sum, next_larger
                );
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:838-859
            // When `b0` holds, jump the predecessor index to `next_larger`; a predecessor can only
            // be selected when `b0`, `b1`, and `b2` are all false.
            // A true `b0` cannot select this predecessor. Keep the
            // same jump, but defer coordinate loads to contenders.
            // Tracing still evaluates all three reported predicates.
            if !$is_target_hsp {
                continue;
            }
        }

        let qo = $helper.q_off_trim;
        let so = $helper.s_off_trim;
        let b1 = qo <= $h_qe;
        let b2 = so <= $h_se;

        // Select this predecessor only when all three rejection flags are false (link_hsps.c:852).
        if !(b0 || b1 || b2) {
            // Use saved $current_idx to retrieve the predecessor HSP index from the helper.
            let j = $helper.hsp_idx;
            if $is_target_hsp {
                eprintln!(
                    "  [LINK] current_idx={}, j={}, sum={} > h_sum={}, linking to HSP {}",
                    $current_idx, j, sum, $h_sum, j
                );
                eprintln!(
                    "    Linked HSP: q_aa={}-{}, s_aa={}-{}, score={}, num={}",
                    $group_hits[j].q_aa_start,
                    $group_hits[j].q_end_int4(),
                    $group_hits[j].s_aa_start,
                    $group_hits[j].s_end_int4(),
                    $pool_hsp_links[j].score,
                    $pool_hsp_links[j].num[1]
                );
            }
            $h_num = $pool_hsp_links[j].num[1];
            $h_sum = $pool_hsp_links[j].sum[1];
            $h_xsum = $pool_hsp_links[j].xsum[1];
            $h_link = j;
        } else if $is_target_hsp {
            eprintln!("  [SKIP] current_idx={}, conditions: b0={} (sum<={}), b1={} (qo<={}), b2={} (so<={})",
                $current_idx, b0, $h_sum, b1, $h_qe, b2, $h_se);
        }
    }};
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-863
// Walk large-gap predecessors backward, jumping to `next_larger` when `b0` holds and recording
// the best predecessor's count, sum, scaled sum, and link index.
// A call per HSP provides a repeated real function-entry boundary within a long
// threaded-WASI group invocation. Native tests can also check the indexed driver.
// The loop, trace timing, initial state and returned values are unchanged;
// all slices are borrowed, with no allocation.
#[cfg(any(test, all(losat_wasi_threads, feature = "wasm-threads")))]
#[cfg_attr(all(losat_wasi_threads, feature = "wasm-threads"), inline(never))]
#[cfg_attr(not(all(losat_wasi_threads, feature = "wasm-threads")), inline(always))]
#[allow(clippy::too_many_arguments)]
fn scan_large_gap_predecessors(
    h_lh_idx: usize,
    pool_lh_helpers: &[LhHelper],
    pool_hsp_links: &[HspLink],
    group_hits: &[UngappedHit],
    h_qe: i32,
    h_se: i32,
    is_target_hsp: bool,
    mut h_sum: i64,
    mut h_num: i16,
    mut h_xsum: f64,
    mut h_link: usize,
) -> (i64, i16, f64, usize) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:675-684,876-885
    // Follow `next_larger` backward while the current sum is at least the previous entry's sum
    // and the previous index remains positive.
    // Every jump is backward. A borrowed prefix excludes sentinels 0/1;
    // its last element is the same predecessor and truncation is metadata-only.
    // One macro owns the unchanged selection and tracing statements for both
    // cursor representations. Native and plain WASI retain indexed traversal.

    #[cfg(all(losat_wasi_threads, feature = "wasm-threads"))]
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827
        // Preserve the no-visit case before borrowing any helper slice.
        if h_lh_idx - 1 > 1 {
            let mut prefix = &pool_lh_helpers[2..h_lh_idx];
            while let Some((helper, tail)) = prefix.split_last() {
                let current_idx = prefix.len() + 1;
                visit_large_gap_predecessor!(prefix, current_idx, helper, prefix, tail;
                h_sum, h_num, h_xsum, h_link, is_target_hsp, h_qe, h_se, pool_hsp_links, group_hits);
            }
        }
    }
    #[cfg(not(all(losat_wasi_threads, feature = "wasm-threads")))]
    {
        let mut j_lh_idx = h_lh_idx - 1;
        while j_lh_idx > 1 {
            let current_idx = j_lh_idx;
            let helper = &pool_lh_helpers[current_idx];
            visit_large_gap_predecessor!(indexed, current_idx, helper, j_lh_idx, j_lh_idx - 1;
                h_sum, h_num, h_xsum, h_link, is_target_hsp, h_qe, h_se, pool_hsp_links, group_hits);
        }
    }
    (h_sum, h_num, h_xsum, h_link)
}

// ---------------------------------------------------------------------------
// Per-group linking (internal)
// ---------------------------------------------------------------------------

use super::params::LinkHspCutoffs;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:974,1049,1061,1080-1081
// The shadow check compares what s_BlastEvenGapLinkHSPs leaves in the list: the order of the
// HSPs, evalue, num and xsum of each. The E-value is compared bit for bit.
/// `LOSAT_LINK_FAST_SHADOW`: the two kernels must return the same HSPs in the
/// same order with every field equal, E-values bit for bit.
fn assert_same_linking(expected: &[UngappedHit], actual: &[UngappedHit]) {
    assert_eq!(
        expected.len(),
        actual.len(),
        "LOSAT_LINK_FAST_SHADOW: group size"
    );
    for (index, (a, b)) in expected.iter().zip(actual).enumerate() {
        assert!(
            a.e_value.to_bits() == b.e_value.to_bits() && format!("{a:?}") == format!("{b:?}"),
            "LOSAT_LINK_FAST_SHADOW: HSP {index} of {}\n  literal port:       {a:?}\n  incremental kernel: {b:?}",
            expected.len()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,589
// The default per-group port of s_BlastEvenGapLinkHSPs: the incremental kernel of linking_incr.rs,
// for every group and also in a run with diagnostics on. LOSAT_LINK_FAST_SHADOW only adds a run of
// the literal port and a comparison; it changes no value of the result.
/// Links one group with the incremental kernel of `linking_incr.rs`, which
/// returns the same result as the literal port `link_hsp_group_ncbi`. Under
/// `LOSAT_LINK_FAST_SHADOW` the literal port also links the group (its
/// diagnostics output off) and the two results must be equal.
#[allow(clippy::too_many_arguments)]
fn link_hsp_group_incremental(
    group_hits: Vec<UngappedHit>,
    params: &KarlinParams,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
    pool_lh_helpers: &mut Vec<LhHelper>,
    pool_hsp_links: &mut Vec<HspLink>,
) -> Vec<UngappedHit> {
    // No NCBI counterpart: the incremental kernel's large-gap sweep assumes every end is after its
    // start (`linking_incr.rs`). A list that holds an HSP of the one-hit word finder of length 0 or
    // less links with the literal port (`x_hsp_inplace::link_in_place_int4`), so none arrives here
    // (debug builds check it).
    debug_assert!(
        group_hits.iter().all(|hit| !hit.is_int4_inverted()),
        "an HSP of length 0 or less reached the incremental linking kernel"
    );
    let expected = link_shadow_enabled().then(|| {
        link_hsp_group_ncbi(
            group_hits.clone(),
            params,
            cutoffs,
            gap_decay_rate,
            false,
            subject_len_nucl,
            query_contexts,
            subject_frame_bases,
            length_adj_per_context,
            eff_searchsp_per_context,
            log_k_by_ctx,
            pool_lh_helpers,
            pool_hsp_links,
        )
    });
    let linked = super::linking_incr::link_hsp_group_incr(
        group_hits,
        cutoffs,
        gap_decay_rate,
        subject_len_nucl,
        query_contexts,
        length_adj_per_context,
        eff_searchsp_per_context,
        log_k_by_ctx,
    );
    if let Some(expected) = expected {
        assert_same_linking(&expected, &linked);
    }
    linked
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:414-419,589
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
// `link_hsp_group_incremental` for `x_link_slices.rs` on the
// group's records where they lie, a sub-slice of the linking call's list, as NCBI links a sublist
// of its one array. The same kernel with the same arguments; under LOSAT_LINK_FAST_SHADOW the
// literal port links a copy of the group (`to_vec`, as `clone` above) and the two results must be
// equal. It changes no value of the result.
/// `link_hsp_group_incremental` on the group's records in place.
#[allow(clippy::too_many_arguments)]
fn link_hsp_group_incremental_slice(
    group_hits: &mut [UngappedHit],
    params: &KarlinParams,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],
    pool_lh_helpers: &mut Vec<LhHelper>,
    pool_hsp_links: &mut Vec<HspLink>,
) {
    // No NCBI counterpart: the incremental kernel's large-gap sweep assumes every end is after its
    // start (`linking_incr.rs`). A list that holds an HSP of the one-hit word finder of length 0 or
    // less links with the literal port (`x_hsp_inplace::link_in_place_int4`), so none arrives here
    // (debug builds check it).
    debug_assert!(
        group_hits.iter().all(|hit| !hit.is_int4_inverted()),
        "an HSP of length 0 or less reached the incremental linking kernel"
    );
    let expected = link_shadow_enabled().then(|| {
        link_hsp_group_ncbi(
            group_hits.to_vec(),
            params,
            cutoffs,
            gap_decay_rate,
            false,
            subject_len_nucl,
            query_contexts,
            subject_frame_bases,
            length_adj_per_context,
            eff_searchsp_per_context,
            log_k_by_ctx,
            pool_lh_helpers,
            pool_hsp_links,
        )
    });
    super::linking_incr::link_hsp_group_incr_slice(
        group_hits,
        cutoffs,
        gap_decay_rate,
        subject_len_nucl,
        query_contexts,
        length_adj_per_context,
        eff_searchsp_per_context,
        log_k_by_ctx,
    );
    if let Some(expected) = expected {
        assert_same_linking(&expected, group_hits);
    }
}

/// NCBI-style HSP linking for moderate-sized groups
///
/// # Arguments
/// * `group_hits` - HSPs in this group (same q_idx, s_idx, q_strand, s_strand)
/// * `params` - Karlin-Altschul parameters
/// * `cutoffs` - Pre-computed NCBI-style cutoffs for this subject
/// * `gap_decay_rate` - Gap decay rate for E-value calculation
/// * `diag_enabled` - Whether diagnostics are enabled
/// * `subject_len_nucl` - Subject length in nucleotides
/// * `query_contexts` - Query context information
/// * `subject_frame_bases` - Subject frame base offsets (unused, kept for API)
/// * `length_adj_per_context` - Pre-computed length adjustment per context
///   (NCBI: `query_info->contexts[ctx].length_adjustment`)
/// * `eff_searchsp_per_context` - Pre-computed effective search space per context
///   (NCBI: `query_info->contexts[ctx].eff_searchsp`)
/// * `log_k_by_ctx` - Pre-computed log(K) per context (NCBI: `kbp[ctx]->logK`)
/// * `pool_lh_helpers` - Thread-local pool for lh_helper buffer reuse
///   (NCBI: link_hsps.c:452-454 allocates once with size `max(1024, hspcnt + 5)`)
/// * `pool_hsp_links` - Thread-local pool for hsp_link buffer reuse
fn link_hsp_group_ncbi(
    mut group_hits: Vec<UngappedHit>,
    params: &KarlinParams,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    diag_enabled: bool,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    _subject_frame_bases: &[i32], // Kept for API compatibility; no longer used after frame-relative coord fix
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
    log_k_by_ctx: &[f64],                // Pre-computed log(K) per context
    pool_lh_helpers: &mut Vec<LhHelper>, // Thread-local pool for lh_helper reuse
    pool_hsp_links: &mut Vec<HspLink>,   // Thread-local pool for hsp_link reuse
) -> Vec<UngappedHit> {
    if group_hits.is_empty() {
        return group_hits;
    }

    let n = group_hits.len();

    // ========================================================================
    // NCBI parity: Use frame-relative coordinates for sorting, NOT concatenated
    // absolute offsets.
    //
    // NCBI's hsp->query.offset and hsp->subject.offset are frame-relative
    // coordinates (0-indexed) after s_AdjustInitialHSPOffsets subtracts the
    // context offset. LOSAT pre-groups by strand, so within a group, we compare
    // frame-relative coordinates directly.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:330-376 (`s_RevCompareHSPsTbx`, the
    // comparator of tblastx's reverse sort)
    // ========================================================================
    #[inline]
    fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
        // NCBI parity: hsp->query.offset is frame-relative (0-indexed after
        // s_AdjustInitialHSPOffsets subtracts context offset).
        // LOSAT's hit.q_aa_start is already 0-indexed frame-relative. The ends are NCBI's `Int4`
        // (`UngappedHit::q_end_int4`): a one-hit list may hold an HSP of length 0 or less, whose
        // end lies at or before its start (`x_hsp_inplace::link_in_place_int4`).
        (
            hit.q_aa_start as i32,
            hit.q_end_int4(),
            hit.s_aa_start as i32,
            hit.s_end_int4(),
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:483-559
    // NCBI does not re-sort within each frame group; it preserves the initial
    // qsort order. Therefore we must not reorder `group_hits` here.

    // NCBI uses the first HSP in this frame/strand group to select the query context
    // used for effective length/search-space in sum-statistics.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:559-562
    let query_context = group_hits[0].ctx_idx as usize;
    let query_len_aa = query_contexts[query_context].aa_len as i64;
    let subject_len_aa = (subject_len_nucl / 3).max(1);

    // ========================================================================
    // NCBI Parity: Use pre-computed length_adjustment and eff_searchsp
    // ========================================================================
    // NCBI stores these in query_info->contexts[ctx] via BLAST_CalcEffLengths
    // (blast_setup.c:700-850) and references them in link_hsps.c:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:560-562
    //
    // These are the SAME pre-computed values as the rest of the search uses (utils.rs), for NCBI
    // parity.
    // ========================================================================
    let length_adjustment = length_adj_per_context[query_context];
    let eff_search_space = eff_searchsp_per_context[query_context];

    // Clamp the effective query length to at least 1 after subtracting the length adjustment.
    let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;

    // LOCAL effective subject length for BLAST_SmallGapSumE/LargeGapSumE arguments
    // (uses 1/3 adjustment for subject per link_hsps.c:566-571)
    let length_adj_for_subject = length_adjustment / 3; // integer division
    let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;

    // Use pre-computed cutoffs from NCBI algorithm
    let cutoff_small = cutoffs.cutoff_small_gap;
    let cutoff_big = cutoffs.cutoff_big_gap;
    let gap_prob = cutoffs.gap_prob;
    let ignore_small_gaps = cutoffs.ignore_small_gaps();

    // Statistics for filtering analysis (long sequences only)
    let stats_index0_filtered = 0usize;
    let stats_index0_passed = 0usize;
    let mut stats_index1_filtered = 0usize;
    let mut stats_index1_passed = 0usize;
    let is_long_sequence = subject_len_nucl > 600_000;

    // Debug output for first group only (to avoid spam)
    static DEBUG_ONCE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if diag_enabled && !DEBUG_ONCE.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[DEBUG] link_hsp_group_ncbi first group:");
        eprintln!(
            "  query_context={}, q_aa_len={}, subject_len_nucl={}, subject_aa_len={}",
            query_context, query_len_aa, subject_len_nucl, subject_len_aa
        );
        eprintln!(
            "  length_adjustment={} (pre-computed), length_adj_for_subject={}",
            length_adjustment, length_adj_for_subject
        );
        eprintln!(
            "  eff_query_len={:.2}, eff_subject_len={:.2} (local)",
            eff_query_len, eff_subject_len
        );
        eprintln!(
            "  eff_search_space={} (pre-computed, Int8)",
            eff_search_space
        );
        eprintln!(
            "  lambda={}, k={}, logK={:.4}",
            params.lambda,
            params.k,
            params.k.glibc_ln()
        );
    }

    // Debug: print cutoffs for first group
    static DEBUG_CUTOFFS: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if diag_enabled && !DEBUG_CUTOFFS.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!(
            "[DEBUG] Linking cutoffs: small_gap={}, big_gap={}, ignore_small_gaps={}, gap_prob={}",
            cutoff_small, cutoff_big, ignore_small_gaps, gap_prob
        );
    }

    // NOTE: We already sorted group_hits above using NCBI's reverse comparator semantics.

    // NCBI `c++/src/algo/blast/core/link_hsps.c` does NOT pre-initialize per-HSP e-values here.
    // Every HSP receives its (possibly linked-set) e-value during the
    // extraction phase (NCBI lines 955-980).

    // Debug: trace a specific HSP (via LOSAT_TRACE_HSP) through linking.
    // Also supports legacy `LOSAT_DEBUG_CHAINING` to enable debug prints.
    let trace_target = trace_hsp_target();
    let debug_chaining =
        trace_target.is_some() || crate::utils::xdiag::var("LOSAT_DEBUG_CHAINING").is_ok();
    let trace_chain_target = trace_chain_target();
    let trace_link_selections = trace_link_selection_summaries();

    // Initialize HspLink array (indices 0..n correspond to HSPs)
    // NCBI uses 1-based indexing with lh_helper[0] as sentinel
    // We use 0-based but handle the logic appropriately
    // NCBI memory pooling: reuse pool_hsp_links buffer instead of allocating fresh
    // Resize to exact group size (n), clearing and reusing existing capacity
    pool_hsp_links.clear();
    pool_hsp_links.reserve(n);
    for (i, hit) in group_hits.iter().enumerate() {
        // NCBI parity: Use frame-relative coordinates for trim calculations
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:545-549
        let (q_off, q_end, s_off, s_end) = frame_relative_coords(hit);
        let q_len_quarter = (q_end - q_off) / 4;
        let s_len_quarter = (s_end - s_off) / 4;
        let qt = TRIM_SIZE.min(q_len_quarter);
        let st = TRIM_SIZE.min(s_len_quarter);
        let score = hit.raw_score;
        let ctx_idx = hit.ctx_idx as usize;
        // Per-context `Lambda` and `logK` from `kbp` for this HSP's context (link_hsps.c:750-752).
        let ctx_lambda = query_contexts[ctx_idx].karlin_params.lambda;
        let ctx_log_k = log_k_by_ctx[ctx_idx];
        let xscore = normalize_score(score, ctx_lambda, ctx_log_k);

        pool_hsp_links.push(HspLink {
            score,
            ctx_idx,
            q_off_trim: q_off + qt,
            s_off_trim: s_off + st,
            q_end_trim: q_end - qt,
            s_end_trim: s_end - st,
            sum: [
                i64::from(score) - i64::from(cutoff_small),
                i64::from(score) - i64::from(cutoff_big),
            ],
            xsum: [xscore, xscore],
            num: [1, 1],
            link: [SENTINEL_IDX, SENTINEL_IDX],
            changed: true,
            linked_to: 0,
            start_of_chain: false, // Will be set to true for chain heads
            linked_set: false,     // Will be set to true for multi-HSP chains
            // Initialize intrusive linked list: 0 -> 1 -> 2 -> ... -> n-1 -> SENTINEL
            next_active: if i + 1 < n { i + 1 } else { SENTINEL_IDX },
            prev_active: if i > 0 { i - 1 } else { SENTINEL_IDX },
        });
    }
    // Now pool_hsp_links has exactly n elements
    // Use pool_hsp_links directly as hsp_links throughout the function
    // Since pool_hsp_links is &mut Vec<HspLink>, we use it directly
    // All references to hsp_links below should use pool_hsp_links instead

    // Resolve the traced target HSP (if any) within this group (post-sort order).
    let target_hsp_idx: Option<usize> = trace_target.and_then(|t| {
        group_hits
            .iter()
            .position(|h| trace_hit_matches_target(t, h))
    });
    let chain_target_hsp_idx: Option<usize> = trace_chain_target.and_then(|t| {
        group_hits
            .iter()
            .position(|h| trace_hit_matches_target(t, h))
    });

    if debug_chaining {
        if let Some(i) = target_hsp_idx {
            let hit = &group_hits[i];
            let (q_start, q_end, s_start, s_end) = int4_nucleotide_coords(hit);
            let link = &pool_hsp_links[i];
            eprintln!(
                "[TRACE_LINKING] target_in_group idx={} raw_score={} q_frame={} s_frame={} q={}-{} s={}-{} q_aa={}-{} s_aa={}-{} q_trim={}..{} s_trim={}..{} cutoff_small={} cutoff_big={}",
                i,
                hit.raw_score,
                hit.q_frame,
                hit.s_frame,
                q_start,
                q_end,
                s_start,
                s_end,
                hit.q_aa_start,
                hit.q_end_int4(),
                hit.s_aa_start,
                hit.s_end_int4(),
                link.q_off_trim,
                link.q_end_trim,
                link.s_off_trim,
                link.s_end_trim,
                cutoff_small,
                cutoff_big,
            );
        }
    }

    // Head of active HSP list (first active HSP index, or SENTINEL_IDX if empty)
    let mut active_head: usize = if n > 0 { 0 } else { SENTINEL_IDX };

    // NCBI lh_helper array structure (lines 579-583):
    // slot 0 is an unused extra end marker
    // lh_helper[1] = hp_start = empty entry
    // slots from 2 on hold the HSPs, shifted by two
    // So: pool_hsp_links[j] <-> pool_lh_helpers[j+2]
    //
    // We allocate n+2 entries and use direct index calculation:
    // - pool_lh_helpers[0], pool_lh_helpers[1]: sentinels
    // - pool_lh_helpers[i+2]: corresponds to pool_hsp_links[i]
    // NCBI `c++/src/algo/blast/core/link_hsps.c` allocates `lh_helper` with `calloc`, so the sentinel
    // entries are zero-initialized. We mirror that exactly (sums=0, next_larger=0).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:576 sets `maxsum1` of the
    // first sentinel to -10000.
    // NCBI memory pooling: reuse pool_lh_helpers buffer instead of allocating fresh
    // Resize to n+2 (group size + sentinels), clearing and reusing existing capacity
    // NCBI uses calloc to allocate n+2 entries, all zero-initialized.
    // We use resize to create n+2 default entries, then set the sentinels.
    pool_lh_helpers.clear();
    pool_lh_helpers.resize(
        n + 2,
        LhHelper {
            hsp_idx: SENTINEL_IDX,
            q_off_trim: 0,
            s_off_trim: 0,
            sum: [0, 0], // NCBI: calloc zero-initializes
            next_larger: 0,
            maxsum1: 0, // Default; sentinels will override
        },
    );
    // Initialize sentinels (indices 0 and 1)
    pool_lh_helpers[0] = LhHelper {
        hsp_idx: SENTINEL_IDX,
        q_off_trim: 0,
        s_off_trim: 0,
        sum: [0, 0], // NCBI: calloc zero-initializes
        next_larger: 0,
        maxsum1: -10000, // NCBI line 576
    };
    pool_lh_helpers[1] = LhHelper {
        hsp_idx: SENTINEL_IDX,
        q_off_trim: 0,
        s_off_trim: 0,
        sum: [0, 0], // NCBI: calloc zero-initializes
        next_larger: 0,
        maxsum1: -10000, // Will be reset to -10000 in each pass (line 1191)
    };
    // Remaining entries (2..n+2) will be populated during linking passes
    // Use pool_lh_helpers directly throughout the function (all lh_helpers references replaced)

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:654-658:
    // Allocate once per sorted group; stable leaves survive all cached rounds.
    let mut maximum_tree = DualMaximumTree::new(n);

    let mut remaining = n;
    let mut selection_round = 0usize;
    let mut first_pass = true;
    // NCBI: path_changed tracks whether any removed HSP had linked_to > 0
    // When `path_changed` is zero, we can skip chain validation entirely.
    let mut path_changed = true;
    // NCBI `c++/src/algo/blast/core/link_hsps.c` takes `gap_prob` from the link parameters (:468)
    // regardless of `ignore_small_gaps` (small-gap cutoff 0, :495). With small gaps the E-values of
    // sets of more than one HSP are divided by `gap_prob` and `1-gap_prob` (:918-935); with
    // large gaps only they are not divided (:939-952). gap_prob is 0 for small search spaces
    // (blast_parameters.c:1076).

    while remaining > 0 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:603-652
        // Diagnostic logs preserve the incoming `path_changed` value so we can
        // distinguish NCBI's cached-path selection from a full DP recomputation.
        let path_changed_at_loop_start = path_changed;
        // NCBI lines 607-625: Initialize max sums each pass
        // CRITICAL: must reset to -cutoff each pass to allow all HSPs to be candidates
        // (`best_sum` is reset where the sums are recomputed, below; with `use_current_max`
        // only `best` is read, as NCBI's max0/max1 are local to the block for `first_pass` false)
        let mut best: [Option<usize>; 2] = [None, None];
        let mut best_sum: [i64; 2];
        let mut use_current_max = false;

        // NCBI lines 603-652: Try to reuse previous max
        if !first_pass {
            // NCBI lines 607-625: Find the current max sums
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:610-623
            // Update `max0`/`best` and `max1`/`best` whenever a candidate's stored sum is at least
            // the current maximum.
            // NCBI does NOT check `linked_to` or `changed` flags here.
            // Processed HSPs are physically removed from the linked list (lines 975-978),
            // so they cannot appear during active-list traversal.
            // In LOSAT, we also unlink processed HSPs (lines 1623-1642), so they cannot
            // appear in `active_head` traversal either. No defensive checks needed.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:608-633:
            // Read the exact stored-score maxima for both channels before unchanged path validation.
            (best, _) = maximum_tree.maxima([cutoff_small, cutoff_big], ignore_small_gaps);

            // Reuse the current maxima directly when no path changed (link_hsps.c:635).
            if !path_changed {
                // No path was changed, use these max sums directly
                use_current_max = true;
            } else {
                // NCBI lines 639-650: If max path hasn't changed, we can use it
                // Walk down `best`, giving up if a removed item (`linked_to` equal to -1000) is encountered.
                use_current_max = true;
                if !ignore_small_gaps {
                    if let Some(bi) = best[0] {
                        // NCBI line 645: check if best[0] itself or its chain has processed HSP
                        if pool_hsp_links[bi].linked_to == -1000 {
                            use_current_max = false;
                        } else {
                            let mut cur = bi;
                            loop {
                                let p = pool_hsp_links[cur].link[0];
                                if p == SENTINEL_IDX {
                                    break;
                                }
                                if pool_hsp_links[p].linked_to == -1000 {
                                    use_current_max = false;
                                    break;
                                }
                                cur = p;
                            }
                        }
                    }
                }
                if use_current_max {
                    if let Some(bi) = best[1] {
                        // NCBI line 648: check if best[1] itself or its chain has processed HSP
                        if pool_hsp_links[bi].linked_to == -1000 {
                            use_current_max = false;
                        } else {
                            let mut cur = bi;
                            loop {
                                let p = pool_hsp_links[cur].link[1];
                                if p == SENTINEL_IDX {
                                    break;
                                }
                                if pool_hsp_links[p].linked_to == -1000 {
                                    use_current_max = false;
                                    break;
                                }
                                cur = p;
                            }
                        }
                    }
                }
            }
        }

        if !use_current_max {
            // NCBI `c++/src/algo/blast/core/link_hsps.c` (s_BlastEvenGapLinkHSPs) lines 659-686:
            //
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:659-685
            //
            // Performance note: NCBI reuses a pre-allocated `lh_helper` buffer.
            // We do the same (no truncate/resize) and only overwrite the active prefix.

            // Build lh_helpers for the current active list in ONE pass (also gives active_count)
            // pool_lh_helpers[0] and [1] are sentinels; active HSPs occupy [2..lh_len).
            let mut lh_len = 2usize;
            let mut cur = active_head;
            // NCBI lines 669-671: track max sum[1] for maxsum1 calculation
            // Since LOSAT groups by frame, all HSPs have the same frame sign, so we track one max
            let mut running_max = -10000i64;
            while cur != SENTINEL_IDX {
                let hsp_idx = cur;

                let next_cur = pool_hsp_links[hsp_idx].next_active;

                // Reset the incoming-link counter for each active HSP (link_hsps.c:685).
                pool_hsp_links[hsp_idx].linked_to = 0;

                // NCBI lines 664-668: copy data + store ptr
                pool_lh_helpers[lh_len].hsp_idx = hsp_idx;
                pool_lh_helpers[lh_len].q_off_trim = pool_hsp_links[hsp_idx].q_off_trim;
                pool_lh_helpers[lh_len].s_off_trim = pool_hsp_links[hsp_idx].s_off_trim;
                pool_lh_helpers[lh_len].sum = [
                    pool_hsp_links[hsp_idx].sum[0],
                    pool_hsp_links[hsp_idx].sum[1],
                ];

                // NCBI lines 669-671: update maxsum1
                running_max = running_max.max(pool_hsp_links[hsp_idx].sum[1]);
                pool_lh_helpers[lh_len].maxsum1 = running_max;

                // Compute `next_larger` by following prior links while `cur_sum` is at least
                // `prev_sum` and `prev` stays positive (link_hsps.c:675-684).
                let cur_sum = pool_lh_helpers[lh_len].sum[1];
                let mut prev = lh_len - 1;
                let mut prev_sum = pool_lh_helpers[prev].sum[1];
                while cur_sum >= prev_sum && prev > 0 {
                    prev = pool_lh_helpers[prev].next_larger;
                    prev_sum = pool_lh_helpers[prev].sum[1];
                }
                pool_lh_helpers[lh_len].next_larger = prev;

                lh_len += 1;
                cur = next_cur;
            }
            // Initialize sentinel slot 1's `maxsum1` to -10000 (link_hsps.c:688).
            pool_lh_helpers[1].maxsum1 = -10000;
            let _active_count = lh_len - 2;

            best = [None, None];
            best_sum = [-i64::from(cutoff_small), -i64::from(cutoff_big)];

            // INDEX 0 LOOP (small gaps) - NCBI lines 691-768
            if !ignore_small_gaps {
                // Walk active HSPs starting at helper index 2 (`pool_lh_helpers[2..]`).
                for h_lh_idx in 2..lh_len {
                    let i = pool_lh_helpers[h_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)

                    let mut h_sum = 0i64;
                    let mut h_xsum = 0.0f64;
                    let mut h_num = 0i16;
                    let mut h_link: usize = SENTINEL_IDX;

                    // Only attempt linking when the HSP's score strictly exceeds the channel cutoff (link_hsps.c:702).
                    let is_target_hsp = if debug_chaining {
                        target_hsp_idx == Some(i)
                    } else {
                        false
                    };
                    let i_score = pool_hsp_links[i].score;
                    if i_score > cutoff_small {
                        let h_qe = pool_hsp_links[i].q_end_trim;
                        let h_se = pool_hsp_links[i].s_end_trim;
                        let h_qe_gap = h_qe + WINDOW_SIZE;
                        let h_se_gap = h_se + WINDOW_SIZE;

                        if is_target_hsp {
                            eprintln!(
                                "[DEBUG_CHAINING] INDEX 0 (small gap) linking for target HSP:"
                            );
                            eprintln!(
                                "  h_qe: {}, h_se: {}, h_qe_gap: {}, h_se_gap: {}",
                                h_qe, h_se, h_qe_gap, h_se_gap
                            );
                        }

                        // Inner loop: NCBI lines 706-742
                        for j_lh_idx in (2..h_lh_idx).rev() {
                            let helper = &pool_lh_helpers[j_lh_idx];
                            let qo = helper.q_off_trim;
                            let so = helper.s_off_trim;
                            let sum = helper.sum[0];

                            // NCBI line 717
                            if qo > h_qe_gap + TRIM_SIZE {
                                if is_target_hsp {
                                    eprintln!(
                                        "  [BREAK] j_lh_idx={}, qo={} > h_qe_gap+TRIM_SIZE={}",
                                        j_lh_idx,
                                        qo,
                                        h_qe_gap + TRIM_SIZE
                                    );
                                }
                                break;
                            }
                            // NCBI lines 719-724
                            let skip = qo <= h_qe || so <= h_se || qo > h_qe_gap || so > h_se_gap;
                            if skip {
                                if is_target_hsp {
                                    eprintln!("  [SKIP] j_lh_idx={}, qo={}, so={}, conditions: qo<={}? {} so<={}? {} qo>{}? {} so>{}? {}",
                                        j_lh_idx, qo, so, h_qe, qo <= h_qe, h_se, so <= h_se,
                                        h_qe_gap, qo > h_qe_gap, h_se_gap, so > h_se_gap);
                                }
                                continue;
                            }

                            // NCBI line 727: update when `sum` strictly exceeds `H_hsp_sum`.
                            if sum > h_sum {
                                let j = helper.hsp_idx; // Use stored hsp_idx (NCBI: ptr)
                                if is_target_hsp {
                                    eprintln!("  [LINK] j_lh_idx={}, j={}, sum={} > h_sum={}, linking to HSP {}",
                                        j_lh_idx, j, sum, h_sum, j);
                                    eprintln!(
                                        "    Linked HSP: q_aa={}-{}, s_aa={}-{}, score={}, num={}",
                                        group_hits[j].q_aa_start,
                                        group_hits[j].q_end_int4(),
                                        group_hits[j].s_aa_start,
                                        group_hits[j].s_end_int4(),
                                        pool_hsp_links[j].score,
                                        pool_hsp_links[j].num[0]
                                    );
                                }
                                h_num = pool_hsp_links[j].num[0];
                                h_sum = pool_hsp_links[j].sum[0];
                                h_xsum = pool_hsp_links[j].xsum[0];
                                h_link = j;
                            }
                        }
                    } else if is_target_hsp {
                        eprintln!("[DEBUG_CHAINING] INDEX 0 (small gap): score {} <= cutoff_small {}, skipping",
                            i_score, cutoff_small);
                    }

                    // NCBI lines 750-767: Update this HSP's link info (the sum in i64, see LhHelper)
                    let new_sum = h_sum + (i64::from(i_score) - i64::from(cutoff_small));
                    let ctx_idx = pool_hsp_links[i].ctx_idx;
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:750-752
                    // Preserve C's left-associative evaluation order; adding `score` times `Lambda`
                    // less `logK` as one group changes 1-ULP chain E-value ties.
                    let new_xsum = h_xsum
                        + (i_score as f64) * query_contexts[ctx_idx].karlin_params.lambda
                        - log_k_by_ctx[ctx_idx];

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:756
                    // `H_hsp_num` is an Int4 and `num` an Int2: the store keeps the low 16 bits, so
                    // the count of a set of more than 32,767 HSPs wraps (as in `core::link_hsps`).
                    let new_num = h_num.wrapping_add(1);

                    // is_target_hsp is already defined at line 1085, reuse it here
                    if is_target_hsp {
                        eprintln!(
                            "[DEBUG_CHAINING] INDEX 0 update: new_sum={}, num={}, link={}",
                            new_sum,
                            new_num,
                            if h_link != SENTINEL_IDX {
                                h_link
                            } else {
                                999999
                            }
                        );
                    }

                    // NCBI lines 755-758: Update sum, num, link, lh_helper
                    pool_hsp_links[i].sum[0] = new_sum;
                    pool_hsp_links[i].num[0] = new_num;
                    pool_hsp_links[i].link[0] = h_link;
                    pool_lh_helpers[h_lh_idx].sum[0] = new_sum;

                    // NCBI line 759-763: update best
                    if new_sum >= best_sum[0] {
                        best_sum[0] = new_sum;
                        best[0] = Some(i);
                        if is_target_hsp {
                            eprintln!(
                                "[DEBUG_CHAINING] INDEX 0: new best_sum={}, best={}",
                                new_sum, i
                            );
                        }
                    }
                    // Store `new_xsum` after updating the channel's best candidate (link_hsps.c:764).
                    pool_hsp_links[i].xsum[0] = new_xsum;
                    // Increment the chosen predecessor's `linked_to` count when linked (link_hsps.c:765-766).
                    if h_link != SENTINEL_IDX {
                        pool_hsp_links[h_link].linked_to += 1;
                    }
                }
            }

            // INDEX 1 LOOP (large gaps) with next_larger optimization - NCBI lines 771-896
            for h_lh_idx in 2..lh_len {
                let i = pool_lh_helpers[h_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)

                let mut h_sum = 0i64;
                let mut h_xsum = 0.0f64;
                let mut h_num = 0i16;
                let mut h_link: usize = SENTINEL_IDX;

                let i_score = pool_hsp_links[i].score;
                // Mark the current HSP as changed before checking the fast-path skip (link_hsps.c:781).
                pool_hsp_links[i].changed = true;
                let prev_link = pool_hsp_links[i].link[1];

                // Define is_target_hsp at loop start for use in all blocks
                let is_target_hsp = if debug_chaining {
                    target_hsp_idx == Some(i)
                } else {
                    false
                };

                // NCBI `c++/src/algo/blast/core/link_hsps.c:781-782` lines 781-795:
                //
                // This is the original fast-path: if the previous best choice exists and
                // was not changed in the last pass, it is still the best.
                // NCBI only checks that the previous link's `changed` flag is zero (no `linked_to`
                // check needed). When an HSP is processed, both `linked_to` (-1000) and `changed`
                // (1) are set, so a cleared `changed` flag suffices to know it is still active.
                let can_skip_ncbi = !first_pass
                    && (prev_link == SENTINEL_IDX || !pool_hsp_links[prev_link].changed);

                if can_skip_ncbi {
                    if prev_link != SENTINEL_IDX {
                        h_num = pool_hsp_links[prev_link].num[1];
                        h_sum = pool_hsp_links[prev_link].sum[1];
                        h_xsum = pool_hsp_links[prev_link].xsum[1];
                    }
                    h_link = prev_link;
                    pool_hsp_links[i].changed = false;
                } else {
                    if is_long_sequence {
                        if i_score > cutoff_big {
                            stats_index1_passed += 1;
                        } else {
                            stats_index1_filtered += 1;
                        }
                    }
                    if i_score > cutoff_big {
                        let h_qe = pool_hsp_links[i].q_end_trim;
                        let h_se = pool_hsp_links[i].s_end_trim;
                        if is_target_hsp {
                            eprintln!(
                                "[DEBUG_CHAINING] INDEX 1 (large gap) linking for target HSP:"
                            );
                            eprintln!(
                                "  h_qe: {}, h_se: {}, prev_link: {}",
                                h_qe,
                                h_se,
                                if prev_link != SENTINEL_IDX {
                                    prev_link
                                } else {
                                    999999
                                }
                            );
                        }

                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:812-823 (the
                        // condition at 812 is the constant 1, so the branch always runs).
                        //
                        // This sets the initial best score to slightly less than the previous
                        // best's real value so ties preserve original ordering.
                        if !first_pass
                            && prev_link != SENTINEL_IDX
                            && pool_hsp_links[prev_link].linked_to >= 0
                        {
                            h_sum = pool_hsp_links[prev_link].sum[1] - 1;
                            if is_target_hsp {
                                eprintln!("  Initial h_sum from prev_link: {}", h_sum);
                            }
                        }

                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:812-863
                        // Enter once per original HSP, after its initial-best setup.
                        #[cfg(all(losat_wasi_threads, feature = "wasm-threads"))]
                        {
                            (h_sum, h_num, h_xsum, h_link) = scan_large_gap_predecessors(
                                h_lh_idx,
                                pool_lh_helpers,
                                pool_hsp_links,
                                &group_hits,
                                h_qe,
                                h_se,
                                is_target_hsp,
                                h_sum,
                                h_num,
                                h_xsum,
                                h_link,
                            );
                        }
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:827-863
                        // Scan predecessor helper entries backward from `h_lh_idx - 1` while the index stays above sentinel 1.
                        // Keep the indexed scan in the caller without a function/tuple boundary.
                        #[cfg(not(all(losat_wasi_threads, feature = "wasm-threads")))]
                        {
                            let mut j_lh_idx = h_lh_idx - 1;
                            while j_lh_idx > 1 {
                                let current_idx = j_lh_idx;
                                let helper = &pool_lh_helpers[current_idx];
                                visit_large_gap_predecessor!(indexed, current_idx, helper, j_lh_idx, j_lh_idx - 1;
                                    h_sum, h_num, h_xsum, h_link, is_target_hsp, h_qe, h_se, pool_hsp_links, group_hits);
                            }
                        }
                    } else if is_target_hsp {
                        eprintln!("[DEBUG_CHAINING] INDEX 1 (large gap): score {} <= cutoff_big {}, skipping",
                            i_score, cutoff_big);
                    }
                }

                // NCBI lines 863-895: Update this HSP's link info (the sum in i64, see LhHelper)
                let new_sum = h_sum + (i64::from(i_score) - i64::from(cutoff_big));
                let ctx_idx = pool_hsp_links[i].ctx_idx;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:865-867
                // Preserve C's left-associative evaluation order; adding `score` times `Lambda`
                // less `logK` as one group changes 1-ULP chain E-value ties.
                let new_xsum = h_xsum
                    + (i_score as f64) * query_contexts[ctx_idx].karlin_params.lambda
                    - log_k_by_ctx[ctx_idx];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:871
                // The Int4 sum stored into the Int2 count wraps, as at index 0.
                let new_num = h_num.wrapping_add(1);

                // is_target_hsp is already defined at line 1218, reuse it here
                if is_target_hsp {
                    eprintln!(
                        "[DEBUG_CHAINING] INDEX 1 update: new_sum={}, num={}, link={}",
                        new_sum,
                        new_num,
                        if h_link != SENTINEL_IDX {
                            h_link
                        } else {
                            999999
                        }
                    );
                }

                // NCBI lines 870-873: Update sum, num, link, lh_helper
                pool_hsp_links[i].sum[1] = new_sum;
                pool_hsp_links[i].num[1] = new_num;
                pool_hsp_links[i].link[1] = h_link;
                pool_lh_helpers[h_lh_idx].sum[1] = new_sum;

                // Update `maxsum1` to the greater of the previous helper's `maxsum1` and `new_sum` (link_hsps.c:874).
                // Note: `maxsum1` is disabled in NCBI's inner loop (link_hsps.c:850), but we compute it for strict parity.
                pool_lh_helpers[h_lh_idx].maxsum1 =
                    pool_lh_helpers[h_lh_idx - 1].maxsum1.max(new_sum);

                // Update next_larger for this entry (NCBI lines 877-884)
                //
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:877-884 (s_BlastEvenGapLinkHSPs)
                let cur_sum = new_sum;
                let mut prev = h_lh_idx - 1;
                let mut prev_sum = pool_lh_helpers[prev].sum[1];
                while cur_sum >= prev_sum && prev > 0 {
                    prev = pool_lh_helpers[prev].next_larger;
                    prev_sum = pool_lh_helpers[prev].sum[1];
                }
                pool_lh_helpers[h_lh_idx].next_larger = prev;

                // NCBI line 887-891: update best
                if new_sum >= best_sum[1] {
                    best_sum[1] = new_sum;
                    best[1] = Some(i);
                    if is_target_hsp {
                        eprintln!(
                            "[DEBUG_CHAINING] INDEX 1: new best_sum={}, best={}",
                            new_sum, i
                        );
                    }
                }
                // Store `new_xsum` after updating the channel's best candidate (link_hsps.c:892).
                pool_hsp_links[i].xsum[1] = new_xsum;
                // Increment the chosen predecessor's `linked_to` count when linked (link_hsps.c:893-894).
                if h_link != SENTINEL_IDX {
                    pool_hsp_links[h_link].linked_to += 1;
                }
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:755-758,869-898:
            // Rebuild only after both existing DP stages finish. Preserve their
            // selected best/sums, including any incoming cached maximum boundary.
            maximum_tree.rebuild(pool_hsp_links, active_head);

            // Clear `path_changed` after the full DP pass finishes (link_hsps.c:897).
            path_changed = false;
            first_pass = false;
        }

        // Select best ordering method (NCBI lines 901-952)
        let mut prob = [f64::MAX, f64::MAX];
        let int4_max = i32::MAX as f64;

        if !ignore_small_gaps {
            if let Some(bi) = best[0] {
                // Add back `num` times the small-gap cutoff that was subtracted during DP (link_hsps.c:907-908).
                pool_hsp_links[bi].sum[0] +=
                    i64::from(pool_hsp_links[bi].num[0]) * i64::from(cutoff_small);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:907-908:
                // Mirror the persistent small-gap cutoff addback in the tree even if another chain is selected.
                maximum_tree.update(bi, DualMaximum::active(bi, pool_hsp_links[bi].sum));

                // The count is NCBI's Int2: a wrapped (negative) count is not "> 1"
                // (link_hsps.c:918). BLAST_GapDecayDivisor takes it as an unsigned and raises the
                // rate to the Int4 `nsegs` less 1; `as usize` sign-extends and gap_decay_divisor
                // truncates `num_segments - 1` to the same i32.
                let num = pool_hsp_links[bi].num[0];
                let xsum = pool_hsp_links[bi].xsum[0];
                let divisor = gap_decay_divisor(gap_decay_rate, num as usize);
                prob[0] = small_gap_sum_e(
                    WINDOW_SIZE,
                    num,
                    xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space,
                    divisor,
                );
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:918-922: with more
                // than one HSP, `prob` becomes `INT4_MAX` when `gap_prob` is 0; otherwise it is
                // divided by `gap_prob` and capped at `INT4_MAX`.
                if num > 1 {
                    if gap_prob == 0.0 || (prob[0] / gap_prob) > int4_max {
                        prob[0] = int4_max;
                    } else {
                        prob[0] /= gap_prob;
                    }
                }
            }

            // NCBI lines 924-935: Also compute prob[1] for large gaps when small gaps enabled
            if let Some(bi) = best[1] {
                let num = pool_hsp_links[bi].num[1];
                let xsum = pool_hsp_links[bi].xsum[1];
                let divisor = gap_decay_divisor(gap_decay_rate, num as usize);
                prob[1] = ncbi_large_gap_sum_e(
                    num,
                    xsum,
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
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:939-952
            // Only large gaps: the E-value is not divided by `1 - gap_prob`.
            if let Some(bi) = best[1] {
                pool_hsp_links[bi].sum[1] +=
                    i64::from(pool_hsp_links[bi].num[1]) * i64::from(cutoff_big);
                // Mirror the persistent addback even if another chain is selected.
                maximum_tree.update(bi, DualMaximum::active(bi, pool_hsp_links[bi].sum));

                let num = pool_hsp_links[bi].num[1];
                let xsum = pool_hsp_links[bi].xsum[1];
                let divisor = gap_decay_divisor(gap_decay_rate, num as usize);
                prob[1] = ncbi_large_gap_sum_e(
                    num,
                    xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space,
                    divisor,
                );
            }
        }

        // NCBI line 940: Select best ordering method
        let ordering = if !ignore_small_gaps && prob[0] <= prob[1] {
            0
        } else {
            1
        };
        // NCBI never has a null `best` here as long as `number_of_hsps > 0`:
        // the DP loops always assign `best[index]` while scanning the list.
        let best_i = if let Some(i) = best[ordering] {
            i
        } else if let Some(i) = best[1 - ordering] {
            i
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:589-982
            // In NCBI, the outer loop over remaining HSPs (`number_of_hsps > 0`) processes ALL HSPs.
            // If active_head becomes SENTINEL_IDX but remaining > 0, this indicates
            // a bug in the linking logic. NCBI would never reach this state.
            // However, we must handle this gracefully to avoid panic.
            if is_long_sequence && remaining > 0 {
                eprintln!("[DEBUG LINKING] WARNING: active_head=SENTINEL but remaining={}, this should not happen in NCBI", remaining);
            }
            // Break out of the loop to avoid panic
            break;
        };

        let evalue = prob[ordering];
        let remaining_before_chain = remaining;

        if diag_enabled {
            // Debug: count chain statistics
            static CHAIN_STATS: std::sync::atomic::AtomicBool =
                std::sync::atomic::AtomicBool::new(false);
            if !CHAIN_STATS.swap(true, std::sync::atomic::Ordering::Relaxed) {
                // Count chain size for first group only
                let chain_len = pool_hsp_links[best_i].num[ordering];
                eprintln!(
                    "[DEBUG] First chain: len={}, e-value={:.2e}, ordering={}",
                    chain_len, evalue, ordering
                );
            }
        }

        // NCBI lines 964-968: If best has linked_to > 0, set path_changed
        // This means some other HSP was linking to this chain
        if pool_hsp_links[best_i].linked_to > 0 {
            path_changed = true;
        }

        // NCBI lines 959-963: Determine if this is a linked set (multi-HSP chain)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:959-962
        // If best has a link, it's a linked set (multiple HSPs)
        let linked_set = pool_hsp_links[best_i].link[ordering] != SENTINEL_IDX;

        // Mark the selected chain head with `start_of_chain` true (link_hsps.c:955).
        // CRITICAL: NCBI sets `start_of_chain` true for ALL selected HSPs (both single and chain heads)
        // before the loop that processes chain members.
        pool_hsp_links[best_i].start_of_chain = true;
        group_hits[best_i].start_of_chain = true;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1045
        // The selected HSPs leave the list (link_hsps.c:965-980), so the head's count is
        // not changed again before the HSPs are hooked up (`replay_ncbi_output_list`).
        group_hits[best_i].hsp_link_num = pool_hsp_links[best_i].num[ordering];

        // Mark chain head and remove chain from consideration (NCBI lines 961-1000)
        // IMPORTANT: Set the chain's E-value to ALL members of the chain,
        // not just the head. This matches NCBI's behavior where all HSPs
        // in a chain share the same E-value.
        let mut is_first = true;
        let mut cur = best_i;
        let collect_chain_members =
            debug_chaining || chain_target_hsp_idx.is_some() || trace_link_selections;
        let mut chain_members = Vec::new();
        loop {
            if debug_chaining && target_hsp_idx == Some(cur) {
                eprintln!("[DEBUG_CHAINING] Target HSP is in chain! ordering={}, evalue={:.2e}, is_first={}",
                    ordering, evalue, is_first);
                eprintln!(
                    "  Chain head: best_i={}, chain length: {}",
                    best_i, pool_hsp_links[best_i].num[ordering]
                );
            }
            // Avoid heap allocation in normal runs: chain_members is only used for debug output.
            if collect_chain_members {
                chain_members.push(cur);
            }
            // Mark `path_changed` when more than one HSP linked to this chain member (link_hsps.c:968).
            if pool_hsp_links[cur].linked_to > 1 {
                path_changed = true;
            }
            pool_hsp_links[cur].linked_to = -1000;
            pool_hsp_links[cur].changed = true;

            // === O(1) UNLINK from active list (NCBI: physically remove from linked list) ===
            {
                let prev_idx = pool_hsp_links[cur].prev_active;
                let next_idx = pool_hsp_links[cur].next_active;

                if prev_idx != SENTINEL_IDX {
                    pool_hsp_links[prev_idx].next_active = next_idx;
                } else {
                    // cur was the head, update active_head
                    active_head = next_idx;
                }

                if next_idx != SENTINEL_IDX {
                    pool_hsp_links[next_idx].prev_active = prev_idx;
                }

                // Mark as not in list
                pool_hsp_links[cur].next_active = SENTINEL_IDX;
                pool_hsp_links[cur].prev_active = SENTINEL_IDX;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:965-978:
            // Remove every physical chain member without renumbering survivors.
            maximum_tree.update(cur, DualMaximum::INACTIVE);

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:972
            // Record `linked_set` for ALL HSPs in the chain (both head and members).
            pool_hsp_links[cur].linked_set = linked_set;

            // Record which ordering method (0 for small gap, 1 for large gap) selected this chain
            // for output (link_hsps.c:973).
            group_hits[cur].ordering_method = ordering as u8;

            // Set E-value for ALL HSPs in the chain (not just head)
            // This is key: chain members inherit the chain's E-value
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:974
            group_hits[cur].e_value = evalue;

            // Transfer `linked_set` to `UngappedHit` for output filtering (link_hsps.c:972).
            group_hits[cur].linked_set = linked_set;

            let next = pool_hsp_links[cur].link[ordering];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1027-1028
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1046-1056
            group_hits[cur].chain_next_link_id = if next == SENTINEL_IDX {
                None
            } else {
                Some(group_hits[next].link_id)
            };

            if is_first {
                // Chain head: `start_of_chain` was already set to true before the loop (link_hsps.c:955).
                is_first = false;
            } else {
                // Non-head HSP in chain: mark `start_of_chain` false; only the head keeps it true
                // while subsequent chain members remain false from initialization.
                pool_hsp_links[cur].start_of_chain = false;
                group_hits[cur].start_of_chain = false;
            }
            remaining -= 1;

            if next == SENTINEL_IDX {
                break;
            }
            cur = next;
        }

        let chain_contains_trace_target = chain_target_hsp_idx
            .map(|target| chain_members.contains(&target))
            .unwrap_or(false);

        if trace_link_selections || chain_contains_trace_target {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:901-982
            let describe_best = |candidate: Option<usize>| -> String {
                let Some(idx) = candidate else {
                    return "NA".to_string();
                };
                let hit = &group_hits[idx];
                let link = &pool_hsp_links[idx];
                format!(
                    "idx={} raw={} num0={} sum0={} xsum0={:.6e} num1={} sum1={} xsum1={:.6e} q_aa={}..{} s_aa={}..{} q_frame={} s_frame={}",
                    idx,
                    hit.raw_score,
                    link.num[0],
                    link.sum[0],
                    link.xsum[0],
                    link.num[1],
                    link.sum[1],
                    link.xsum[1],
                    hit.q_aa_start,
                    hit.q_end_int4(),
                    hit.s_aa_start,
                    hit.s_end_int4(),
                    hit.q_frame,
                    hit.s_frame
                )
            };
            let head = &group_hits[best_i];
            let (head_q_start, head_q_end, head_s_start, head_s_end) = int4_nucleotide_coords(head);
            eprintln!(
                "[TRACE_LINK_SELECTION] context={} subject_sign={} round={} remaining_before={} used_current_max={} path_changed_in={} selected={} evalue={:.6e} chain_len={} linked_set={} target_in_chain={} prob_small={:.6e} prob_large={:.6e} head_idx={} head_raw={} head_q={}-{} head_s={}-{} head_q_frame={} head_s_frame={} best_small={} best_large={}",
                query_context,
                group_hits[0].s_frame.signum(),
                selection_round,
                remaining_before_chain,
                use_current_max,
                path_changed_at_loop_start,
                ordering,
                evalue,
                chain_members.len(),
                linked_set,
                chain_contains_trace_target,
                prob[0],
                prob[1],
                best_i,
                head.raw_score,
                head_q_start,
                head_q_end,
                head_s_start,
                head_s_end,
                head.q_frame,
                head.s_frame,
                describe_best(best[0]),
                describe_best(best[1])
            );
        }

        if chain_contains_trace_target {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1046-1056
            // Print the replay chain that NCBI follows from a selected chain
            // head, without enabling the per-candidate DP trace.
            eprintln!(
                "[TRACE_LINK_SELECTION] target_chain_members head_idx={} ordering={} evalue={:.6e}",
                best_i, ordering, evalue
            );
            for &idx in &chain_members {
                let hit = &group_hits[idx];
                let (q_start, q_end, s_start, s_end) = int4_nucleotide_coords(hit);
                let mark = if chain_target_hsp_idx == Some(idx) {
                    "*"
                } else {
                    " "
                };
                eprintln!(
                    "  {}idx={} raw_score={} q={}-{} s={}-{} q_frame={} s_frame={} q_aa={}-{} s_aa={}-{} start_of_chain={}",
                    mark,
                    idx,
                    hit.raw_score,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                    hit.q_frame,
                    hit.s_frame,
                    hit.q_aa_start,
                    hit.q_end_int4(),
                    hit.s_aa_start,
                    hit.s_end_int4(),
                    hit.start_of_chain,
                );
            }
        }

        if debug_chaining && chain_members.iter().any(|&idx| target_hsp_idx == Some(idx)) {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:901-952
            let describe_best = |candidate: Option<usize>| -> String {
                let Some(idx) = candidate else {
                    return "NA".to_string();
                };
                let hit = &group_hits[idx];
                let link = &pool_hsp_links[idx];
                format!(
                    "idx={} raw={} num0={} sum0={} xsum0={:.6e} num1={} sum1={} xsum1={:.6e} q_aa={}..{} s_aa={}..{}",
                    idx,
                    hit.raw_score,
                    link.num[0],
                    link.sum[0],
                    link.xsum[0],
                    link.num[1],
                    link.sum[1],
                    link.xsum[1],
                    hit.q_aa_start,
                    hit.q_end_int4(),
                    hit.s_aa_start,
                    hit.s_end_int4()
                )
            };
            eprintln!(
                "[TRACE_LINKING] ordering_choice prob_small={:.6e} prob_large={:.6e} selected={} best_small={} best_large={}",
                prob[0],
                prob[1],
                ordering,
                describe_best(best[0]),
                describe_best(best[1])
            );
            eprintln!(
                "[TRACE_LINKING] target_chain head_idx={} chain_len={} linked_set={} ordering={} evalue={:.3e}",
                best_i,
                pool_hsp_links[best_i].num[ordering],
                linked_set,
                ordering,
                evalue
            );
            for &idx in &chain_members {
                let hit = &group_hits[idx];
                let (q_start, q_end, s_start, s_end) = int4_nucleotide_coords(hit);
                let next = pool_hsp_links[idx].link[ordering];
                let next_str = if next == SENTINEL_IDX {
                    "END".to_string()
                } else {
                    next.to_string()
                };
                let mark = if target_hsp_idx == Some(idx) {
                    "*"
                } else {
                    " "
                };
                eprintln!(
                    "  {}idx={} raw_score={} q={}-{} s={}-{} q_aa={}-{} s_aa={}-{} start_of_chain={} next={}",
                    mark,
                    idx,
                    hit.raw_score,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                    hit.q_aa_start,
                    hit.q_end_int4(),
                    hit.s_aa_start,
                    hit.s_end_int4(),
                    hit.start_of_chain,
                    next_str,
                );
            }
        }

        selection_round += 1;
    }

    // NOTE: NCBI link_hsps.c:1016-1020 filters chain members (linked_set=true,
    // start_of_chain=false) during OUTPUT phase, not here. The flags are set
    // so that downstream output code can perform the filtering.
    //
    // Output filtering should happen in the format/output code, not here.
    // This function returns ALL HSPs with their linked_set/start_of_chain flags
    // properly set for downstream use.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:589-982
    // In NCBI, the outer loop over remaining HSPs (`number_of_hsps > 0`) processes ALL HSPs.
    // HSPs that are not selected as best in one iteration are still processed
    // in subsequent iterations. Therefore, there are NO unprocessed HSPs in NCBI.
    // If active_head becomes SENTINEL_IDX but remaining > 0, this indicates
    // a bug in the linking logic, not a valid case for E-value calculation.

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:492-496
    // NCBI applies the cutoffs without unconditional stderr diagnostics.
    if (debug_chaining || trace_link_selections || std::env::var_os("LOSAT_DIAGNOSTICS").is_some())
        && is_long_sequence
        && (stats_index0_filtered > 0
            || stats_index0_passed > 0
            || stats_index1_filtered > 0
            || stats_index1_passed > 0)
    {
        let total_index0 = stats_index0_filtered + stats_index0_passed;
        let total_index1 = stats_index1_filtered + stats_index1_passed;
        eprintln!(
            "[DEBUG LINKING_FILTER] subject_len_nucl={}, group_size={}",
            subject_len_nucl, n
        );
        eprintln!("[DEBUG LINKING_FILTER] INDEX 0 (small gap): cutoff={}, filtered={}, passed={}, total={}, filter_rate={:.2}%",
            cutoff_small, stats_index0_filtered, stats_index0_passed, total_index0,
            if total_index0 > 0 { (stats_index0_filtered as f64 / total_index0 as f64) * 100.0 } else { 0.0 });
        eprintln!("[DEBUG LINKING_FILTER] INDEX 1 (large gap): cutoff={}, filtered={}, passed={}, total={}, filter_rate={:.2}%",
            cutoff_big, stats_index1_filtered, stats_index1_passed, total_index1,
            if total_index1 > 0 { (stats_index1_filtered as f64 / total_index1 as f64) * 100.0 } else { 0.0 });
    }

    group_hits
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "linking_tests.rs"]
mod tests;

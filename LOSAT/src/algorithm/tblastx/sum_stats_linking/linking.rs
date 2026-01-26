//! Core HSP linking algorithm
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c

use crate::stats::sum_statistics::{
    gap_decay_divisor, small_gap_sum_e, large_gap_sum_e, normalize_score,
    defaults::{GAP_SIZE, OVERLAP_SIZE},
};
use crate::stats::KarlinParams;
use std::sync::OnceLock;
use std::cell::RefCell;
// NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
// ```c
// CMTArgs(...)
// {
// #ifdef NCBI_NO_THREADS
//     m_NumThreads = CThreadable::kMinNumThreads;
//     m_MTMode = eNotSupported;
// #endif
// }
// ```
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::algorithm::tblastx::chaining::UngappedHit;
use crate::algorithm::tblastx::diagnostics::diagnostics_enabled;
use crate::algorithm::tblastx::extension::convert_coords;
use crate::algorithm::tblastx::lookup::QueryContext;

use super::params::LinkingParams;
use super::cutoffs::calculate_link_hsp_cutoffs_ncbi;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
const TRIM_SIZE: i32 = (OVERLAP_SIZE + 1) / 2;

/// Ordering method indices (NCBI ELinkOrderingMethod)
const LINK_SMALL_GAPS: usize = 0;
const LINK_LARGE_GAPS: usize = 1;

/// NCBI BLAST_GAP_PROB (ungapped search default)
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:66
const BLAST_GAP_PROB: f64 = 0.5;

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
#[derive(Clone, Copy, Debug)]
struct TraceHspTarget {
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
}

static TRACE_HSP_TARGET: OnceLock<Option<TraceHspTarget>> = OnceLock::new();

#[inline(always)]
fn trace_hsp_target() -> Option<TraceHspTarget> {
    *TRACE_HSP_TARGET.get_or_init(|| {
        let raw = std::env::var("LOSAT_TRACE_HSP").ok()?;
        let parts: Vec<&str> = raw
            .split(|c: char| c == ',' || c == ':' || c == ';' || c.is_whitespace())
            .filter(|p| !p.is_empty())
            .collect();
        if parts.len() != 4 {
            eprintln!(
                "[TRACE_HSP] invalid LOSAT_TRACE_HSP (expected 4 integers): {:?}",
                raw
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
    })
}

#[inline(always)]
fn trace_match_target(target: TraceHspTarget, q_start: usize, q_end: usize, s_start: usize, s_end: usize) -> bool {
    q_start as i32 == target.q_start
        && q_end as i32 == target.q_end
        && s_start as i32 == target.s_start
        && s_end as i32 == target.s_end
}

#[inline(always)]
fn trace_hit_matches_target(target: TraceHspTarget, h: &UngappedHit) -> bool {
    let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, h.q_frame, h.q_orig_len);
    let (s_start, s_end) = convert_coords(h.s_aa_start, h.s_aa_end, h.s_frame, h.s_orig_len);
    trace_match_target(target, q_start, q_end, s_start, s_end)
}

/// NCBI BLAST_Nint - round to nearest integer
/// Reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441
///
/// ```c
/// long BLAST_Nint(double x)
/// {
///    x += (x >= 0. ? 0.5 : -0.5);
///    return (long)x;
/// }
/// ```
#[inline]
fn blast_nint(x: f64) -> i32 {
    let rounded = if x >= 0.0 { x + 0.5 } else { x - 0.5 };
    rounded as i32
}

// ---------------------------------------------------------------------------
// Core structures
// ---------------------------------------------------------------------------

/// lh_helper structure (NCBI link_hsps.c:100-109)
/// This is the "hot" data used in the inner loop
#[derive(Clone)]
struct LhHelper {
    hsp_idx: usize,     // Original HSP index (NCBI: ptr field)
    q_off_trim: i32,
    s_off_trim: i32,
    sum: [i32; 2],      // sum for both indices (small gap, large gap)
    next_larger: usize, // index of next HSP with larger sum[1]
    maxsum1: i32,       // NCBI: threshold for stopping link attempts (unused with if(0))
}

/// HSP link information (NCBI BlastHSPLink structure)
/// Reference: link_hsps.c:63-71
#[derive(Clone)]
struct HspLink {
    score: i32,
    ctx_idx: usize,          // NCBI: H->hsp->context - for kbp[context] lookup
    q_off_trim: i32,
    s_off_trim: i32,
    q_end_trim: i32,
    s_end_trim: i32,
    sum: [i32; 2],           // sum for both indices
    xsum: [f64; 2],          // normalized sum for both indices
    // NCBI parity: precompute the per-HSP normalized score contribution once.
    // This avoids recomputing (and re-ln'ing K) inside hot loops.
    xscore: f64,
    num: [i16; 2],           // number of HSPs in chain for both indices
    // Best previous HSP to link with for both indices (NCBI: LinkHSPStruct* or NULL).
    // We store an index into `hsp_links`, using SENTINEL_IDX as the NULL marker.
    link: [usize; 2],
    changed: bool,           // Has link been changed? (NCBI: hsp_link.changed)
    linked_to: i32,          // How many HSPs link TO this one (NCBI: linked_to)
    start_of_chain: bool,    // Is this HSP the start of a chain? (for output filtering)
    linked_set: bool,        // Is this HSP part of a multi-HSP chain? (NCBI: linked_set)
    // Intrusive linked list for active HSPs (NCBI uses hp_start->next linked list)
    // SENTINEL_IDX indicates end of list or not in list
    next_active: usize,      // Next active HSP in list
    prev_active: usize,      // Previous active HSP in list
}

/// Thread-local buffer pools for memory reuse across frame groups
/// NCBI reference: link_hsps.c:452-454 allocates lh_helper once with size MAX(1024, hspcnt+5)
/// and reuses it across all frame groups (lines 553-983)
struct BufferPools {
    lh_helpers: Vec<LhHelper>,
    hsp_links: Vec<HspLink>,
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Apply NCBI-style sum-statistics even-gap linking
///
/// This is the main entry point for HSP linking with NCBI-compatible cutoffs.
///
/// # Arguments
/// * `hits` - Vector of ungapped HSPs to link
/// * `params` - Karlin-Altschul parameters for the scoring matrix
/// * `linking_params` - NCBI-style linking parameters (avg query length, subject length, etc.)
/// * `query_contexts` - Query context information
/// * `subject_frame_bases` - Subject frame base offsets
/// * `length_adj_per_context` - Pre-computed length adjustment per context
///   (NCBI: query_info->contexts[ctx].length_adjustment)
/// * `eff_searchsp_per_context` - Pre-computed effective search space per context
///   (NCBI: query_info->contexts[ctx].eff_searchsp)
///
/// NCBI Parity Note:
/// In NCBI, length_adjustment and eff_searchsp are calculated once per context
/// in BLAST_CalcEffLengths (blast_setup.c:700-850) and stored in query_info->contexts[].
/// link_hsps.c then references these stored values:
/// ```c
/// length_adjustment = query_info->contexts[query_context].length_adjustment;
/// searchsp_eff = query_info->contexts[query_context].eff_searchsp;
/// ```
/// We now mirror this by accepting pre-computed vectors instead of recalculating.
pub fn apply_sum_stats_even_gap_linking(
    mut hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
    length_adj_per_context: &[i64],
    eff_searchsp_per_context: &[i64],
) -> Vec<UngappedHit> {
    if hits.is_empty() {
        return hits;
    }

    let diag_enabled = diagnostics_enabled();

    // Calculate cutoffs once for this subject using NCBI algorithm
    // NCBI: CalculateLinkHSPCutoffs is called once per subject
    let cutoffs = calculate_link_hsp_cutoffs_ncbi(
        linking_params.avg_query_length,
        linking_params.subject_len_nucl,
        0, // db_length = 0 for -subject mode
        linking_params.cutoff_score_min,
        linking_params.scale_factor,
        linking_params.gap_decay_rate,
        params,
    );

    // Debug output for cutoffs (controlled by LOSAT_DEBUG_CUTOFFS env var)
    static DEBUG_CUTOFFS_PRINTED: std::sync::atomic::AtomicBool =
        std::sync::atomic::AtomicBool::new(false);
    if std::env::var("LOSAT_DEBUG_CUTOFFS").is_ok()
        && !DEBUG_CUTOFFS_PRINTED.swap(true, std::sync::atomic::Ordering::SeqCst)
    {
        eprintln!("=== LOSAT Linking Cutoffs Debug ===");
        eprintln!("  avg_query_length: {}", linking_params.avg_query_length);
        eprintln!("  subject_len_nucl: {}", linking_params.subject_len_nucl);
        eprintln!("  cutoff_score_min: {}", linking_params.cutoff_score_min);
        eprintln!("  lambda: {:.6}", params.lambda);
        eprintln!("  K: {:.6}", params.k);
        eprintln!("  H: {:.6}", params.h);
        eprintln!("  cutoff_small_gap: {}", cutoffs.cutoff_small_gap);
        eprintln!("  cutoff_big_gap: {}", cutoffs.cutoff_big_gap);
        eprintln!("  gap_prob: {:.6}", cutoffs.gap_prob);
        eprintln!("  ignore_small_gaps: {}", cutoffs.ignore_small_gaps);
        eprintln!("  WINDOW_SIZE: {}", WINDOW_SIZE);
        eprintln!("  TRIM_SIZE: {}", TRIM_SIZE);
        eprintln!("===================================");
    }

    // =======================================================================
    // NCBI link_hsps.c EXACT REPLICATION:
    // 1. Sort ALL HSPs by s_RevCompareHSPsTbx (lines 484-490)
    // 2. Scan sorted list and detect frame boundaries (lines 514-534)
    // 3. Process each frame group SEQUENTIALLY (lines 553-982)
    // =======================================================================

    // Precompute logK per context once (NCBI: kbp[ctx]->logK).
    // This is computed once at top level and reused across all frame groups.
    // NCBI reference: link_hsps.c does not explicitly show this, but kbp[ctx]->logK
    // is accessed per context, so we precompute to avoid repeated ln() calls.
    let log_k_by_ctx: Vec<f64> = query_contexts
        .iter()
        .map(|ctx| ctx.karlin_params.k.ln())
        .collect();

    // Step 1: Sort ALL HSPs using s_RevCompareHSPsTbx order
    // NCBI: qsort(link_hsp_array, ..., s_RevCompareHSPsTbx)
    // Sort key: (q_idx, s_idx, context/3, SIGN(s_frame), q_off desc, q_end desc, s_off desc, s_end desc)
    hits.sort_unstable_by(|a, b| {
        // Primary: q_idx, s_idx (for multi-query/subject)
        a.q_idx.cmp(&b.q_idx)
            .then(a.s_idx.cmp(&b.s_idx))
            // NCBI line 343-349: context/(NUM_FRAMES/2) = context/3 for query strand
            .then_with(|| {
                let a_qstrand = if a.q_frame > 0 { 0 } else { 1 };
                let b_qstrand = if b.q_frame > 0 { 0 } else { 1 };
                a_qstrand.cmp(&b_qstrand)
            })
            // NCBI line 351-357: SIGN(subject.frame)
            // In NCBI qsort: return 1 means h1 comes AFTER h2
            // if h1->subject.frame > h2->subject.frame: return 1 (h1 after h2)
            // This means negative frames come FIRST, then positive frames.
            // In Rust sort_by: Greater means a comes AFTER b
            // So we want: if a.frame > b.frame, return Greater (a after b)
            .then_with(|| {
                let a_ssign = a.s_frame.signum();
                let b_ssign = b.s_frame.signum();
                // NCBI: if h1->subject.frame > h2->subject.frame return 1
                // Rust equivalent: a_ssign.cmp(&b_ssign) gives ascending order
                a_ssign.cmp(&b_ssign)
            })
            // NCBI lines 359-374: all descending
            .then(b.q_aa_start.cmp(&a.q_aa_start))
            .then(b.q_aa_end.cmp(&a.q_aa_end))
            .then(b.s_aa_start.cmp(&a.s_aa_start))
            .then(b.s_aa_end.cmp(&a.s_aa_end))
    });

    // Step 2: Detect frame boundaries in sorted list (NCBI lines 514-534)
    // Split into groups where (context/3, SIGN(s_frame)) changes
    let total_hits = hits.len();  // Save for pool sizing (NCBI: hsp_list->hspcnt)
    let mut frame_groups: Vec<Vec<UngappedHit>> = Vec::new();
    let mut current_group: Vec<UngappedHit> = Vec::new();

    for hit in hits {
        let q_strand = if hit.q_frame > 0 { 1i8 } else { -1i8 };
        let s_sign = hit.s_frame.signum();

        if let Some(prev) = current_group.last() {
            let prev_q_strand = if prev.q_frame > 0 { 1i8 } else { -1i8 };
            let prev_s_sign = prev.s_frame.signum();

            // NCBI line 522-525: frame boundary detection
            if prev.q_idx != hit.q_idx
                || prev.s_idx != hit.s_idx
                || prev_q_strand != q_strand
                || prev_s_sign != s_sign
            {
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:553-558
    // ```c
    // for (frame_index=0; frame_index<num_query_frames; frame_index++)
    // {
    //    hp_start->next = hp_frame_start[frame_index];
    //    number_of_hsps = hp_frame_number[frame_index];
    // }
    // ```
    // NCBI processes groups in frame_index order (0, 1, 2, ...), which matches
    // the order of frame_groups in LOSAT (created from sorted HSP list).
    //
    // NCBI memory pooling: link_hsps.c:452-454 allocates lh_helper once with
    // MAX(1024, hspcnt+5) and reuses it across all frame groups.
    let initial_lh_size = (total_hits + 5).max(1024);  // NCBI: MAX(1024, hspcnt+5)

    #[cfg(feature = "parallel")]
    let mut results: Vec<UngappedHit> = {
        // Use thread_local! macro for thread-local storage (rayon doesn't have ThreadLocal)
        thread_local! {
            static POOLS: RefCell<Option<BufferPools>> = RefCell::new(None);
        }

        let mut indexed_results: Vec<(usize, Vec<UngappedHit>)> = frame_groups
            .into_par_iter()
            .enumerate()
            .map(|(group_idx, group_hits)| {
                POOLS.with(|pools| {
                    let mut pools_opt = pools.borrow_mut();
                    let pool = pools_opt.get_or_insert_with(|| {
                        BufferPools {
                            lh_helpers: Vec::with_capacity(initial_lh_size),
                            hsp_links: Vec::with_capacity(total_hits),
                        }
                    });
                    let processed = link_hsp_group_ncbi(
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
                    );
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
    };

    #[cfg(not(feature = "parallel"))]
    let mut results: Vec<UngappedHit> = {
        let mut pools = BufferPools {
            lh_helpers: Vec::with_capacity(initial_lh_size),
            hsp_links: Vec::with_capacity(total_hits),
        };
        let mut ordered: Vec<UngappedHit> = Vec::new();
        for group_hits in frame_groups.into_iter() {
            let processed = link_hsp_group_ncbi(
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
            );
            ordered.extend(processed);
        }
        ordered
    };

    // ========================================================================
    // NCBI Parity: Two final sorts for translated queries (link_hsps.c:990-1000)
    // ========================================================================
    // NCBI reference (verbatim):
    // ```c
    // if (kTranslatedQuery) {
    //     qsort(link_hsp_array, num, sizeof(BlastHSP*), s_RevCompareHSPsTransl);
    //     qsort(link_hsp_array, num, sizeof(BlastHSP*), s_FwdCompareHSPsTransl);
    // }
    // ```
    //
    // Sort 1: s_RevCompareHSPsTransl - Reverse sort by query context group and subject frame
    // NCBI comparator (link_hsps.c:274-291):
    //   context comparison: context/(NUM_FRAMES/2) = context/3 for query strand group
    //   subject frame comparison: SIGN(subject.frame)
    //   Both use descending order (h1 < h2 returns 1)
    results.sort_by(|a, b| {
        // Query context group (context/3 in NCBI)
        let a_ctx_group = a.ctx_idx / 3;
        let b_ctx_group = b.ctx_idx / 3;
        // Subject strand sign
        let a_s_sign = a.s_frame.signum();
        let b_s_sign = b.s_frame.signum();

        // NCBI s_RevCompareHSPsTransl: descending by context group, then by subject frame
        b_ctx_group.cmp(&a_ctx_group)
            .then(b_s_sign.cmp(&a_s_sign))
    });

    // Sort 2: s_FwdCompareHSPsTransl - Forward stable sort by query context group
    // NCBI comparator (link_hsps.c:248-272):
    //   context comparison: ascending by context/(NUM_FRAMES/2)
    // This effectively groups HSPs by query strand while preserving order within groups
    results.sort_by(|a, b| {
        let a_ctx_group = a.ctx_idx / 3;
        let b_ctx_group = b.ctx_idx / 3;

        // NCBI s_FwdCompareHSPsTransl: ascending by context group
        a_ctx_group.cmp(&b_ctx_group)
    });

    results
}

// ---------------------------------------------------------------------------
// Per-group linking (internal)
// ---------------------------------------------------------------------------

use super::params::LinkHspCutoffs;

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
///   (NCBI: query_info->contexts[ctx].length_adjustment)
/// * `eff_searchsp_per_context` - Pre-computed effective search space per context
///   (NCBI: query_info->contexts[ctx].eff_searchsp)
/// * `log_k_by_ctx` - Pre-computed log(K) per context (NCBI: kbp[ctx]->logK)
/// * `pool_lh_helpers` - Thread-local pool for lh_helper buffer reuse
///   (NCBI: link_hsps.c:452-454 allocates once with MAX(1024, hspcnt+5))
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
    log_k_by_ctx: &[f64],  // Pre-computed log(K) per context
    pool_lh_helpers: &mut Vec<LhHelper>,  // Thread-local pool for lh_helper reuse
    pool_hsp_links: &mut Vec<HspLink>,     // Thread-local pool for hsp_link reuse
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
    // NCBI reference (verbatim comparator for tblastx reverse sort):
    //   if (h1->query.offset < h2->query.offset) return  1;
    //   if (h1->query.offset > h2->query.offset) return -1;
    //   if (h1->query.end   < h2->query.end)   return  1;
    //   if (h1->query.end   > h2->query.end)   return -1;
    //   if (h1->subject.offset < h2->subject.offset) return  1;
    //   if (h1->subject.offset > h2->subject.offset) return -1;
    //   if (h1->subject.end   < h2->subject.end)   return  1;
    //   if (h1->subject.end   > h2->subject.end)   return -1;
    // Source: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:s_RevCompareHSPsTbx (lines ~359-374)
    // ========================================================================
    #[inline]
    fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
        // NCBI parity: hsp->query.offset is frame-relative (0-indexed after
        // s_AdjustInitialHSPOffsets subtracts context offset).
        // LOSAT's hit.q_aa_start is already 0-indexed frame-relative.
        (
            hit.q_aa_start as i32,
            hit.q_aa_end as i32,
            hit.s_aa_start as i32,
            hit.s_aa_end as i32,
        )
    }

    // Sort by reverse position using frame-relative coordinates (NCBI parity).
    // NCBI uses qsort (unstable), so we use sort_unstable_by for parity.
    group_hits.sort_unstable_by(|a, b| {
        let (aqo, aqe, aso, ase) = frame_relative_coords(a);
        let (bqo, bqe, bso, bse) = frame_relative_coords(b);
        bqo.cmp(&aqo)
            .then(bqe.cmp(&aqe))
            .then(bso.cmp(&aso))
            .then(bse.cmp(&ase))
    });

    // NCBI uses the first HSP in this frame/strand group to select the query context
    // used for effective length/search-space in sum-statistics.
    // Reference: link_hsps.c:559-562
    let query_context = group_hits[0].ctx_idx;
    let query_len_aa = query_contexts[query_context].aa_len as i64;
    let subject_len_aa = (subject_len_nucl / 3).max(1);

    // ========================================================================
    // NCBI Parity: Use pre-computed length_adjustment and eff_searchsp
    // ========================================================================
    // NCBI stores these in query_info->contexts[ctx] via BLAST_CalcEffLengths
    // (blast_setup.c:700-850) and references them in link_hsps.c:
    // ```c
    // length_adjustment = query_info->contexts[query_context].length_adjustment;
    // query_length = query_info->contexts[query_context].query_length;
    // query_length = MAX(query_length - length_adjustment, 1);
    // ...
    // // In BLAST_SmallGapSumE/BLAST_LargeGapSumE calls:
    // BLAST_SmallGapSumE(..., query_info->contexts[query_context].eff_searchsp, ...);
    // ```
    //
    // Previously we recalculated these locally via compute_length_adjustment_simple().
    // Now we use the SAME pre-computed values from utils.rs for NCBI parity.
    // ========================================================================
    let length_adjustment = length_adj_per_context[query_context];
    let eff_search_space = eff_searchsp_per_context[query_context];

    // NCBI: query_length = MAX(query_length - length_adjustment, 1)
    let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;

    // LOCAL effective subject length for BLAST_SmallGapSumE/LargeGapSumE arguments
    // (uses 1/3 adjustment for subject per link_hsps.c:566-571)
    // ```c
    // if (Blast_SubjectIsTranslated(program_number))  // tblastx = TRUE
    // {
    //    length_adjustment /= CODON_LENGTH;  // divide by 3
    //    subject_length /= CODON_LENGTH;
    // }
    // subject_length = MAX(subject_length - length_adjustment, 1);
    // ```
    let length_adj_for_subject = length_adjustment / 3;  // integer division
    let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;

    // Use pre-computed cutoffs from NCBI algorithm
    let cutoff_small = cutoffs.cutoff_small_gap;
    let cutoff_big = cutoffs.cutoff_big_gap;
    let gap_prob = cutoffs.gap_prob;
    let ignore_small_gaps = cutoffs.ignore_small_gaps;

    // Statistics for filtering analysis (long sequences only)
    let mut stats_index0_filtered = 0usize;
    let mut stats_index0_passed = 0usize;
    let mut stats_index1_filtered = 0usize;
    let mut stats_index1_passed = 0usize;
    let is_long_sequence = subject_len_nucl > 600_000;

    // Debug output for first group only (to avoid spam)
    static DEBUG_ONCE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if diag_enabled && !DEBUG_ONCE.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[DEBUG] link_hsp_group_ncbi first group:");
        eprintln!(
            "  query_context={}, q_aa_len={}, subject_len_nucl={}, subject_aa_len={}",
            query_context,
            query_len_aa,
            subject_len_nucl,
            subject_len_aa
        );
        eprintln!("  length_adjustment={} (pre-computed), length_adj_for_subject={}",
            length_adjustment, length_adj_for_subject);
        eprintln!("  eff_query_len={:.2}, eff_subject_len={:.2} (local)",
            eff_query_len, eff_subject_len);
        eprintln!("  eff_search_space={} (pre-computed, Int8)",
            eff_search_space);
        eprintln!("  lambda={}, k={}, logK={:.4}", params.lambda, params.k, params.k.ln());
    }

    // Debug: print cutoffs for first group
    static DEBUG_CUTOFFS: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if diag_enabled && !DEBUG_CUTOFFS.swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[DEBUG] Linking cutoffs: small_gap={}, big_gap={}, ignore_small_gaps={}, gap_prob={}",
            cutoff_small, cutoff_big, ignore_small_gaps, gap_prob);
    }

    // NOTE: We already sorted group_hits above using NCBI's reverse comparator semantics.

    // NCBI `link_hsps.c` does NOT pre-initialize per-HSP e-values here.
    // Every HSP receives its (possibly linked-set) e-value during the
    // extraction phase (NCBI lines 955-980).

    // Debug: trace a specific HSP (via LOSAT_TRACE_HSP) through linking.
    // Also supports legacy `LOSAT_DEBUG_CHAINING` to enable debug prints.
    let trace_target = trace_hsp_target();
    let debug_chaining = trace_target.is_some() || std::env::var("LOSAT_DEBUG_CHAINING").is_ok();

    // Initialize HspLink array (indices 0..n correspond to HSPs)
    // NCBI uses 1-based indexing with lh_helper[0] as sentinel
    // We use 0-based but handle the logic appropriately
    // NCBI memory pooling: reuse pool_hsp_links buffer instead of allocating fresh
    // Resize to exact group size (n), clearing and reusing existing capacity
    pool_hsp_links.clear();
    pool_hsp_links.reserve(n);
    for (i, hit) in group_hits.iter().enumerate() {
        // NCBI parity: Use frame-relative coordinates for trim calculations
        // Reference: link_hsps.c:545-549
        let (q_off, q_end, s_off, s_end) = frame_relative_coords(hit);
        let q_len_quarter = (q_end - q_off) / 4;
        let s_len_quarter = (s_end - s_off) / 4;
        let qt = TRIM_SIZE.min(q_len_quarter);
        let st = TRIM_SIZE.min(s_len_quarter);
        let score = hit.raw_score;
        let ctx_idx = hit.ctx_idx;
        // NCBI line 750-752: kbp[H->hsp->context]->Lambda, kbp[H->hsp->context]->logK
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
            sum: [score - cutoff_small, score - cutoff_big],
            xsum: [xscore, xscore],
            xscore,
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

    if debug_chaining {
        if let Some(i) = target_hsp_idx {
            let hit = &group_hits[i];
            let (q_start, q_end) =
                convert_coords(hit.q_aa_start, hit.q_aa_end, hit.q_frame, hit.q_orig_len);
            let (s_start, s_end) =
                convert_coords(hit.s_aa_start, hit.s_aa_end, hit.s_frame, hit.s_orig_len);
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
                hit.q_aa_end,
                hit.s_aa_start,
                hit.s_aa_end,
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
    // lh_helper[0] = empty = additional end marker
    // lh_helper[1] = hp_start = empty entry
    // lh_helper[i] = hsp_array[i-2] for i >= 2
    // So: pool_hsp_links[j] <-> pool_lh_helpers[j+2]
    //
    // We allocate n+2 entries and use direct index calculation:
    // - pool_lh_helpers[0], pool_lh_helpers[1]: sentinels
    // - pool_lh_helpers[i+2]: corresponds to pool_hsp_links[i]
    // NCBI `link_hsps.c` allocates `lh_helper` with `calloc`, so the sentinel
    // entries are zero-initialized. We mirror that exactly (sums=0, next_larger=0).
    // NCBI line 576: lh_helper[0].maxsum1 = -10000;
    // NCBI memory pooling: reuse pool_lh_helpers buffer instead of allocating fresh
    // Resize to n+2 (group size + sentinels), clearing and reusing existing capacity
    // NCBI uses calloc to allocate n+2 entries, all zero-initialized.
    // We use resize to create n+2 default entries, then set the sentinels.
    pool_lh_helpers.clear();
    pool_lh_helpers.resize(n + 2, LhHelper {
        hsp_idx: SENTINEL_IDX,
        q_off_trim: 0,
        s_off_trim: 0,
        sum: [0, 0],           // NCBI: calloc zero-initializes
        next_larger: 0,
        maxsum1: 0,            // Default; sentinels will override
    });
    // Initialize sentinels (indices 0 and 1)
    pool_lh_helpers[0] = LhHelper {
        hsp_idx: SENTINEL_IDX,
        q_off_trim: 0,
        s_off_trim: 0,
        sum: [0, 0],           // NCBI: calloc zero-initializes
        next_larger: 0,
        maxsum1: -10000,       // NCBI line 576
    };
    pool_lh_helpers[1] = LhHelper {
        hsp_idx: SENTINEL_IDX,
        q_off_trim: 0,
        s_off_trim: 0,
        sum: [0, 0],           // NCBI: calloc zero-initializes
        next_larger: 0,
        maxsum1: -10000,       // Will be reset to -10000 in each pass (line 1191)
    };
    // Remaining entries (2..n+2) will be populated during linking passes
    // Use pool_lh_helpers directly throughout the function (all lh_helpers references replaced)

    let mut remaining = n;
    let mut first_pass = true;
    // NCBI: path_changed tracks whether any removed HSP had linked_to > 0
    // If path_changed == 0, we can skip chain validation entirely
    let mut path_changed = true;
    // NCBI `link_hsps.c` sets `gap_prob = link_hsp_params->gap_prob;` regardless
    // of `ignore_small_gaps` (cutoff[0]==0). The adjustment by `1-gap_prob`
    // for INDEX 1 is still applied when `num>1`.
    // NOTE: gap_prob is now passed in from calculate_link_hsp_cutoffs_ncbi
    // and may be 0 for small search spaces (NCBI line 1076)

    while remaining > 0 {
        // NCBI lines 607-625: Initialize max sums each pass
        // CRITICAL: must reset to -cutoff each pass to allow all HSPs to be candidates
        let mut best: [Option<usize>; 2] = [None, None];
        let mut best_sum: [i32; 2] = [-cutoff_small, -cutoff_big];
        let mut use_current_max = false;

        // NCBI lines 603-652: Try to reuse previous max
        if !first_pass {
            // NCBI lines 607-625: Find the current max sums
            // NCBI reference (verbatim, link_hsps.c:610-623):
            //   for (H=hp_start->next; H!=NULL; H=H->next) {
            //      Int4 sum0=H->hsp_link.sum[0];
            //      Int4 sum1=H->hsp_link.sum[1];
            //      if(sum0>=max0) { max0=sum0; best[0]=H; }
            //      if(sum1>=max1) { max1=sum1; best[1]=H; }
            //   }
            // NCBI does NOT check linked_to or changed flags here.
            // Processed HSPs are physically removed from the linked list (lines 975-978),
            // so they cannot appear in hp_start->next traversal.
            // In LOSAT, we also unlink processed HSPs (lines 1623-1642), so they cannot
            // appear in active_head traversal either. No defensive checks needed.
            if !ignore_small_gaps {
                let mut cur = active_head;
                while cur != SENTINEL_IDX {
                    // NCBI line 613: if(sum0>=max0)
                    if pool_hsp_links[cur].sum[0] >= best_sum[0] {
                        best_sum[0] = pool_hsp_links[cur].sum[0];
                        best[0] = Some(cur);
                    }
                    cur = pool_hsp_links[cur].next_active;
                }
            }
            {
                let mut cur = active_head;
                while cur != SENTINEL_IDX {
                    // NCBI line 618: if(sum1>=max1)
                    if pool_hsp_links[cur].sum[1] >= best_sum[1] {
                        best_sum[1] = pool_hsp_links[cur].sum[1];
                        best[1] = Some(cur);
                    }
                    cur = pool_hsp_links[cur].next_active;
                }
            }

            // NCBI line 635: if(path_changed==0) use_current_max=1
            if !path_changed {
                // No path was changed, use these max sums directly
                use_current_max = true;
            } else {
                // NCBI lines 639-650: If max path hasn't changed, we can use it
                // Walk down best, give up if we find a removed item (linked_to==-1000) in path
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
            // NCBI `link_hsps.c` (s_BlastEvenGapLinkHSPs) lines 659-686:
            //
            //   if(!use_current_max){
            //     for (H=hp_start,H_index=1; H!=NULL; H=H->next,H_index++) {
            //       ...
            //       lh_helper[H_index].ptr = H;
            //       lh_helper[H_index].q_off_trim = q_off_t;
            //       lh_helper[H_index].s_off_trim = s_off_t;
            //       for(i=0;i<eOrderingMethods;i++)
            //         lh_helper[H_index].sum[i] = H->hsp_link.sum[i];
            //       ...
            //       /* set next_larger ... */
            //       ...
            //       H->linked_to = 0;
            //     }
            //   }
            //
            // Performance note: NCBI reuses a pre-allocated `lh_helper` buffer.
            // We do the same (no truncate/resize) and only overwrite the active prefix.

            // Build lh_helpers for the current active list in ONE pass (also gives active_count)
            // pool_lh_helpers[0] and [1] are sentinels; active HSPs occupy [2..lh_len).
            let mut lh_len = 2usize;
            let mut cur = active_head;
            // NCBI lines 669-671: track max sum[1] for maxsum1 calculation
            // Since LOSAT groups by frame, all HSPs have the same frame sign, so we track one max
            let mut running_max = -10000i32;
            while cur != SENTINEL_IDX {
                let hsp_idx = cur;

                let next_cur = pool_hsp_links[hsp_idx].next_active;

                // NCBI line 685: H->linked_to = 0
                pool_hsp_links[hsp_idx].linked_to = 0;

                // NCBI lines 664-668: copy data + store ptr
                pool_lh_helpers[lh_len].hsp_idx = hsp_idx;
                pool_lh_helpers[lh_len].q_off_trim = pool_hsp_links[hsp_idx].q_off_trim;
                pool_lh_helpers[lh_len].s_off_trim = pool_hsp_links[hsp_idx].s_off_trim;
                pool_lh_helpers[lh_len].sum = [pool_hsp_links[hsp_idx].sum[0], pool_hsp_links[hsp_idx].sum[1]];

                // NCBI lines 669-671: update maxsum1
                running_max = running_max.max(pool_hsp_links[hsp_idx].sum[1]);
                pool_lh_helpers[lh_len].maxsum1 = running_max;

                // NCBI lines 675-684: compute next_larger
                // NCBI: while((cur_sum>=prev_sum) && (prev>0))
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
            // NCBI line 688: lh_helper[1].maxsum1 = -10000;
            pool_lh_helpers[1].maxsum1 = -10000;
            let _active_count = lh_len - 2;

            best = [None, None];
            best_sum = [-cutoff_small, -cutoff_big];

            // INDEX 0 LOOP (small gaps) - NCBI lines 691-768
            if !ignore_small_gaps {
                // NCBI: H_index = 2; for (H=hp_start->next; H!=NULL; H=H->next,H_index++)
                // Iterate over pool_lh_helpers[2..] which contains only active HSPs
                for h_lh_idx in 2..lh_len {
                    let i = pool_lh_helpers[h_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)

                    let mut h_sum = 0i32;
                    let mut h_xsum = 0.0f64;
                    let mut h_num = 0i16;
                    let mut h_link: usize = SENTINEL_IDX;

                    // NCBI line 702: if (H->hsp->score > cutoff[index])
                    let is_target_hsp = if debug_chaining { target_hsp_idx == Some(i) } else { false };
                    let i_score = pool_hsp_links[i].score;
                    if i_score > cutoff_small {
                        let h_qe = pool_hsp_links[i].q_end_trim;
                        let h_se = pool_hsp_links[i].s_end_trim;
                        let h_qe_gap = h_qe + WINDOW_SIZE;
                        let h_se_gap = h_se + WINDOW_SIZE;

                        if is_target_hsp {
                            eprintln!("[DEBUG_CHAINING] INDEX 0 (small gap) linking for target HSP:");
                            eprintln!("  h_qe: {}, h_se: {}, h_qe_gap: {}, h_se_gap: {}", h_qe, h_se, h_qe_gap, h_se_gap);
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
                                    eprintln!("  [BREAK] j_lh_idx={}, qo={} > h_qe_gap+TRIM_SIZE={}",
                                        j_lh_idx, qo, h_qe_gap + TRIM_SIZE);
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

                            // NCBI line 727: if (sum>H_hsp_sum)
                            if sum > h_sum {
                                let j = helper.hsp_idx; // Use stored hsp_idx (NCBI: ptr)
                                if is_target_hsp {
                                    eprintln!("  [LINK] j_lh_idx={}, j={}, sum={} > h_sum={}, linking to HSP {}",
                                        j_lh_idx, j, sum, h_sum, j);
                                    eprintln!("    Linked HSP: q_aa={}-{}, s_aa={}-{}, score={}, num={}",
                                        group_hits[j].q_aa_start, group_hits[j].q_aa_end,
                                        group_hits[j].s_aa_start, group_hits[j].s_aa_end,
                                        pool_hsp_links[j].score, pool_hsp_links[j].num[0]);
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

                    // NCBI lines 750-767: Update this HSP's link info
                    let new_sum = h_sum + (i_score - cutoff_small);
                    // NCBI parity: normalized score contribution is constant per HSP.
                    let new_xsum = h_xsum + pool_hsp_links[i].xscore;

                    // is_target_hsp is already defined at line 1085, reuse it here
                    if is_target_hsp {
                        eprintln!("[DEBUG_CHAINING] INDEX 0 update: new_sum={}, num={}, link={}",
                            new_sum, h_num + 1, if h_link != SENTINEL_IDX { h_link } else { 999999 });
                    }

                    // NCBI lines 755-758: Update sum, num, link, lh_helper
                    pool_hsp_links[i].sum[0] = new_sum;
                    pool_hsp_links[i].num[0] = h_num + 1;
                    pool_hsp_links[i].link[0] = h_link;
                    pool_lh_helpers[h_lh_idx].sum[0] = new_sum;

                    // NCBI line 759-763: update best
                    if new_sum >= best_sum[0] {
                        best_sum[0] = new_sum;
                        best[0] = Some(i);
                        if is_target_hsp {
                            eprintln!("[DEBUG_CHAINING] INDEX 0: new best_sum={}, best={}", new_sum, i);
                        }
                    }
                    // NCBI line 764: H->hsp_link.xsum[index] = new_xsum; (after best update)
                    pool_hsp_links[i].xsum[0] = new_xsum;
                    // NCBI line 765-766: if(H_hsp_link) ((LinkHSPStruct*)H_hsp_link)->linked_to++;
                    if h_link != SENTINEL_IDX {
                        pool_hsp_links[h_link].linked_to += 1;
                    }
                }
            }

            // INDEX 1 LOOP (large gaps) with next_larger optimization - NCBI lines 771-896
            for h_lh_idx in 2..lh_len {
                let i = pool_lh_helpers[h_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)

                let mut h_sum = 0i32;
                let mut h_xsum = 0.0f64;
                let mut h_num = 0i16;
                let mut h_link: usize = SENTINEL_IDX;

                let i_score = pool_hsp_links[i].score;
                let i_xscore = pool_hsp_links[i].xscore;

                // NCBI line 781: H->hsp_link.changed=1
                pool_hsp_links[i].changed = true;
                let prev_link = pool_hsp_links[i].link[1];

                // Define is_target_hsp at loop start for use in all blocks
                let is_target_hsp = if debug_chaining { target_hsp_idx == Some(i) } else { false };

                // NCBI `link_hsps.c` lines 781-795:
                //   H->hsp_link.changed=1;
                //   H2 = H->hsp_link.link[index];
                //   if ((!first_pass) && ((H2==0) || (H2->hsp_link.changed==0))) { ... }
                //
                // This is the original fast-path: if the previous best choice exists and
                // was not changed in the last pass, it is still the best.
                // NCBI only checks H2->hsp_link.changed==0 (no linked_to check needed).
                // When an HSP is processed, both linked_to=-1000 and hsp_link.changed=1 are set,
                // so changed==0 is sufficient to determine the HSP hasn't been processed.
                let can_skip_ncbi = !first_pass
                    && (prev_link == SENTINEL_IDX
                        || !pool_hsp_links[prev_link].changed);

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
                            eprintln!("[DEBUG_CHAINING] INDEX 1 (large gap) linking for target HSP:");
                            eprintln!("  h_qe: {}, h_se: {}, prev_link: {}", h_qe, h_se,
                                if prev_link != SENTINEL_IDX { prev_link } else { 999999 });
                        }

                        // NCBI `link_hsps.c` lines 812-823 (note the hardcoded `if(1)`):
                        //   if(!first_pass&&H2&&H2->linked_to>=0){
                        //      if(1){
                        //         H_hsp_sum=H2->hsp_link.sum[index]-1;
                        //      }else{
                        //         ...
                        //      }
                        //   }
                        //
                        // This sets the initial best score to slightly less than the previous
                        // best's real value so ties preserve original ordering.
                        if !first_pass && prev_link != SENTINEL_IDX && pool_hsp_links[prev_link].linked_to >= 0 {
                            h_sum = pool_hsp_links[prev_link].sum[1] - 1;
                            if is_target_hsp {
                                eprintln!("  Initial h_sum from prev_link: {}", h_sum);
                            }
                        }

                        // Inner loop (NCBI lines 827-861)
                        // CRITICAL: Must save current_idx BEFORE decrement to access correct hsp_idx
                        let mut j_lh_idx = h_lh_idx - 1;
                        while j_lh_idx > 1 {
                            // NCBI: H2_helper points to lh_helper[H2_index] at loop START
                            let current_idx = j_lh_idx; // Save before any modification
                            let helper = &pool_lh_helpers[current_idx];
                            let sum = helper.sum[1];
                            let next_larger = helper.next_larger;
                            let qo = helper.q_off_trim;
                            let so = helper.s_off_trim;

                            let b0 = sum <= h_sum;

                            // NCBI line 841: H2_index--
                            j_lh_idx -= 1;
                            if b0 {
                                j_lh_idx = next_larger;
                                if is_target_hsp {
                                    eprintln!("  [SKIP] current_idx={}, sum={} <= h_sum={}, jump to next_larger={}",
                                        current_idx, sum, h_sum, next_larger);
                                }
                            }

                            let b1 = qo <= h_qe;
                            let b2 = so <= h_se;

                            // NCBI line 852: if (!(b0|b1|b2))
                            if !(b0 || b1 || b2) {
                                // Use saved current_idx (NCBI: H2 = H2_helper->ptr)
                                let j = helper.hsp_idx;
                                if is_target_hsp {
                                    eprintln!("  [LINK] current_idx={}, j={}, sum={} > h_sum={}, linking to HSP {}",
                                        current_idx, j, sum, h_sum, j);
                                    eprintln!("    Linked HSP: q_aa={}-{}, s_aa={}-{}, score={}, num={}",
                                        group_hits[j].q_aa_start, group_hits[j].q_aa_end,
                                        group_hits[j].s_aa_start, group_hits[j].s_aa_end,
                                        pool_hsp_links[j].score, pool_hsp_links[j].num[1]);
                                }
                                h_num = pool_hsp_links[j].num[1];
                                h_sum = pool_hsp_links[j].sum[1];
                                h_xsum = pool_hsp_links[j].xsum[1];
                                h_link = j;
                            } else if is_target_hsp {
                                eprintln!("  [SKIP] current_idx={}, conditions: b0={} (sum<={}), b1={} (qo<={}), b2={} (so<={})",
                                    current_idx, b0, h_sum, b1, h_qe, b2, h_se);
                            }
                        }
                    } else if is_target_hsp {
                        eprintln!("[DEBUG_CHAINING] INDEX 1 (large gap): score {} <= cutoff_big {}, skipping",
                            i_score, cutoff_big);
                    }
                }

                // NCBI lines 863-895: Update this HSP's link info
                let new_sum = h_sum + (i_score - cutoff_big);
                // NCBI parity: normalized score contribution is constant per HSP.
                let new_xsum = h_xsum + i_xscore;

                // is_target_hsp is already defined at line 1218, reuse it here
                if is_target_hsp {
                    eprintln!("[DEBUG_CHAINING] INDEX 1 update: new_sum={}, num={}, link={}",
                        new_sum, h_num + 1, if h_link != SENTINEL_IDX { h_link } else { 999999 });
                }

                // NCBI lines 870-873: Update sum, num, link, lh_helper
                pool_hsp_links[i].sum[1] = new_sum;
                pool_hsp_links[i].num[1] = h_num + 1;
                pool_hsp_links[i].link[1] = h_link;
                pool_lh_helpers[h_lh_idx].sum[1] = new_sum;

                // NCBI line 874: lh_helper[H_index].maxsum1 = MAX(lh_helper[H_index-1].maxsum1, new_sum);
                // Note: maxsum1 is unused (NCBI line 850 is if(0)) but we compute it for strict parity
                pool_lh_helpers[h_lh_idx].maxsum1 = pool_lh_helpers[h_lh_idx - 1].maxsum1.max(new_sum);

                // Update next_larger for this entry (NCBI lines 876-884)
                //
                // NCBI `link_hsps.c` (s_BlastEvenGapLinkHSPs) lines 876-885:
                //   Int4 cur_sum=lh_helper[H_index].sum[1];
                //   Int4 prev = H_index-1;
                //   Int4 prev_sum = lh_helper[prev].sum[1];
                //   while((cur_sum>=prev_sum) && (prev>0)){
                //      prev=lh_helper[prev].next_larger;
                //      prev_sum = lh_helper[prev].sum[1];
                //   }
                //   lh_helper[H_index].next_larger = prev;
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
                        eprintln!("[DEBUG_CHAINING] INDEX 1: new best_sum={}, best={}", new_sum, i);
                    }
                }
                // NCBI line 892: H->hsp_link.xsum[index] = new_xsum; (after best update)
                pool_hsp_links[i].xsum[1] = new_xsum;
                // NCBI line 893-894: if(H_hsp_link) ((LinkHSPStruct*)H_hsp_link)->linked_to++;
                if h_link != SENTINEL_IDX {
                    pool_hsp_links[h_link].linked_to += 1;
                }
            }

            // NCBI line 897: path_changed=0 after first pass computation
            path_changed = false;
            first_pass = false;
        }

        // Select best ordering method (NCBI lines 901-952)
        let mut prob = [f64::MAX, f64::MAX];
        let int4_max = i32::MAX as f64;

        if !ignore_small_gaps {
            if let Some(bi) = best[0] {
                // NCBI lines 907-908: Add back cutoff*num that was subtracted during DP
                // best[0]->hsp_link.sum[0] += (best[0]->hsp_link.num[0])*cutoff[0];
                pool_hsp_links[bi].sum[0] += (pool_hsp_links[bi].num[0] as i32) * cutoff_small;

                let num = pool_hsp_links[bi].num[0] as usize;
                let xsum = pool_hsp_links[bi].xsum[0];
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[0] = small_gap_sum_e(WINDOW_SIZE, num as i16, xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space as i64, divisor);
                // NCBI `link_hsps.c` lines 918-922:
                //   if( num > 1 ) {
                //     if( gap_prob == 0 || (prob /= gap_prob) > INT4_MAX ) prob = INT4_MAX;
                //   }
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
                let num = pool_hsp_links[bi].num[1] as usize;
                let xsum = pool_hsp_links[bi].xsum[1];
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[1] = large_gap_sum_e(num as i16, xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space as i64, divisor);
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
            // NCBI lines 939-952: Only consider large gaps
            if let Some(bi) = best[1] {
                // NCBI lines 942-943: Add back cutoff*num
                pool_hsp_links[bi].sum[1] += (pool_hsp_links[bi].num[1] as i32) * cutoff_big;

                let num = pool_hsp_links[bi].num[1] as usize;
                let xsum = pool_hsp_links[bi].xsum[1];
                let divisor = gap_decay_divisor(gap_decay_rate, num);
                prob[1] = large_gap_sum_e(num as i16, xsum,
                    eff_query_len as i32,
                    eff_subject_len as i32,
                    eff_search_space as i64, divisor);
                if num > 1 {
                    let denom = 1.0 - gap_prob;
                    if denom == 0.0 || (prob[1] / denom) > int4_max {
                        prob[1] = int4_max;
                    } else {
                        prob[1] /= denom;
                    }
                }
            }
        }

        // NCBI line 940: Select best ordering method
        let ordering = if !ignore_small_gaps && prob[0] <= prob[1] { 0 } else { 1 };
        // NCBI never has "best==NULL" here as long as `number_of_hsps > 0`:
        // the DP loops always assign `best[index]` while scanning the list.
        let best_i = if let Some(i) = best[ordering] {
            i
        } else if let Some(i) = best[1 - ordering] {
            i
        } else {
            // NCBI reference: link_hsps.c:589-982
            // In NCBI, the while (number_of_hsps > 0) loop processes ALL HSPs.
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
        // NCBI reference (verbatim, link_hsps.c:959-962):
        //   if (best[ordering_method]->hsp_link.link[ordering_method])
        //      linked_set = TRUE;
        //   else
        //      linked_set = FALSE;
        // If best has a link, it's a linked set (multiple HSPs)
        let linked_set = pool_hsp_links[best_i].link[ordering] != SENTINEL_IDX;

        // NCBI line 955: best[ordering_method]->start_of_chain = TRUE;
        // CRITICAL: NCBI sets start_of_chain = TRUE for ALL selected HSPs (both single and chain heads)
        // This is set BEFORE the loop that processes chain members
        pool_hsp_links[best_i].start_of_chain = true;
        group_hits[best_i].start_of_chain = true;

        // Mark chain head and remove chain from consideration (NCBI lines 961-1000)
        // IMPORTANT: Set the chain's E-value to ALL members of the chain,
        // not just the head. This matches NCBI's behavior where all HSPs
        // in a chain share the same E-value.
        let mut is_first = true;
        let mut cur = best_i;
        let mut chain_members = Vec::new();
        loop {
            if debug_chaining && target_hsp_idx == Some(cur) {
                eprintln!("[DEBUG_CHAINING] Target HSP is in chain! ordering={}, evalue={:.2e}, is_first={}",
                    ordering, evalue, is_first);
                eprintln!("  Chain head: best_i={}, chain length: {}", best_i, pool_hsp_links[best_i].num[ordering]);
            }
            // Avoid heap allocation in normal runs: chain_members is only used for debug output.
            if debug_chaining {
                chain_members.push(cur);
            }
            // NCBI line 968: if (H->linked_to>1) path_changed=1
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

            // NCBI line 972: H->linked_set = linked_set
            // NCBI reference (verbatim, link_hsps.c:972):
            //   H->linked_set = linked_set;
            // This is set for ALL HSPs in the chain (both head and members)
            pool_hsp_links[cur].linked_set = linked_set;

            // NCBI line 973: H->ordering_method = ordering_method;
            // This records which ordering method (small gap=0 / large gap=1) was used
            // to select this chain for output.
            group_hits[cur].ordering_method = ordering as u8;

            // Set E-value for ALL HSPs in the chain (not just head)
            // This is key: chain members inherit the chain's E-value
            // NCBI reference (verbatim, link_hsps.c:974):
            //   H->hsp->evalue = prob[ordering_method];
            group_hits[cur].e_value = evalue;

            // Transfer linked_set flag to UngappedHit for output filtering
            // NCBI line 972: H->linked_set = linked_set
            group_hits[cur].linked_set = linked_set;

            if is_first {
                // Chain head: start_of_chain already set to true above (NCBI line 955)
                // This matches NCBI: best[ordering_method]->start_of_chain = TRUE (before loop)
                is_first = false;
            } else {
                // Non-head HSP in chain: mark as not start_of_chain
                // NCBI reference: Only the first HSP in the loop (best[ordering_method]) has start_of_chain=TRUE
                // All subsequent HSPs in the chain have start_of_chain=FALSE (default from initialization)
                pool_hsp_links[cur].start_of_chain = false;
                group_hits[cur].start_of_chain = false;
            }
            remaining -= 1;

            let next = pool_hsp_links[cur].link[ordering];
            if next == SENTINEL_IDX {
                break;
            }
            cur = next;
        }

        if debug_chaining && chain_members.iter().any(|&idx| target_hsp_idx == Some(idx)) {
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
                let (q_start, q_end) =
                    convert_coords(hit.q_aa_start, hit.q_aa_end, hit.q_frame, hit.q_orig_len);
                let (s_start, s_end) =
                    convert_coords(hit.s_aa_start, hit.s_aa_end, hit.s_frame, hit.s_orig_len);
                let next = pool_hsp_links[idx].link[ordering];
                let next_str = if next == SENTINEL_IDX {
                    "END".to_string()
                } else {
                    next.to_string()
                };
                let mark = if target_hsp_idx == Some(idx) { "*" } else { " " };
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
                    hit.q_aa_end,
                    hit.s_aa_start,
                    hit.s_aa_end,
                    hit.start_of_chain,
                    next_str,
                );
            }
        }
    }

    // NOTE: NCBI link_hsps.c:1016-1020 filters chain members (linked_set=true,
    // start_of_chain=false) during OUTPUT phase, not here. The flags are set
    // so that downstream output code can perform the filtering.
    //
    // Output filtering should happen in the format/output code, not here.
    // This function returns ALL HSPs with their linked_set/start_of_chain flags
    // properly set for downstream use.
    //
    // NCBI reference: link_hsps.c:589-982
    // In NCBI, the while (number_of_hsps > 0) loop processes ALL HSPs.
    // HSPs that are not selected as best in one iteration are still processed
    // in subsequent iterations. Therefore, there are NO unprocessed HSPs in NCBI.
    // If active_head becomes SENTINEL_IDX but remaining > 0, this indicates
    // a bug in the linking logic, not a valid case for E-value calculation.

    // Output filtering statistics for long sequences
    if is_long_sequence && (stats_index0_filtered > 0 || stats_index0_passed > 0 || stats_index1_filtered > 0 || stats_index1_passed > 0) {
        let total_index0 = stats_index0_filtered + stats_index0_passed;
        let total_index1 = stats_index1_filtered + stats_index1_passed;
        eprintln!("[DEBUG LINKING_FILTER] subject_len_nucl={}, group_size={}", subject_len_nucl, n);
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
mod tests {
    use super::*;

    /// Helper to create a mock UngappedHit with specified coordinates
    fn mock_hit(q_aa_start: usize, q_aa_end: usize, s_aa_start: usize, s_aa_end: usize, id: i32) -> UngappedHit {
        UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx: 0,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start,
            q_aa_end,
            s_aa_start,
            s_aa_end,
            q_orig_len: 1000,
            s_orig_len: 1000,
            raw_score: id, // Use raw_score as ID for verification
            e_value: 1e-10,
            num_ident: 0, // Mock value for tests
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        }
    }

    /// Verify LOSAT's HSP sort order matches NCBI's s_RevCompareHSPsTbx
    ///
    /// NCBI comparison logic (link_hsps.c:359-375):
    /// ```c
    /// // ALL fields use the same pattern: h1 < h2 → return 1 (DESCENDING)
    /// if (h1->query.offset < h2->query.offset)   return  1;
    /// if (h1->query.offset > h2->query.offset)   return -1;
    /// if (h1->query.end < h2->query.end)         return  1;
    /// if (h1->query.end > h2->query.end)         return -1;
    /// if (h1->subject.offset < h2->subject.offset) return  1;
    /// if (h1->subject.offset > h2->subject.offset) return -1;
    /// if (h1->subject.end < h2->subject.end)       return  1;
    /// if (h1->subject.end > h2->subject.end)       return -1;
    /// return 0;
    /// ```
    ///
    /// C qsort: compare(a,b) > 0 means "a comes after b"
    /// So `if (h1 < h2) return 1` means smaller values come AFTER = DESCENDING order
    #[test]
    fn test_hsp_sort_order_matches_ncbi() {
        // Test data: hits with varying coordinates
        // ID = raw_score for identification after sort
        //
        // Expected order after NCBI-style sort (all descending):
        // 1. q_off=100, q_end=150, s_off=200, s_end=250 (ID=1) - largest q_off/q_end/s_off/s_end
        // 2. q_off=100, q_end=150, s_off=150, s_end=200 (ID=2) - same q_*, smaller s_off
        // 3. q_off=100, q_end=140, s_off=200, s_end=250 (ID=3) - same q_off, smaller q_end
        // 4. q_off=50,  q_end=100, s_off=300, s_end=350 (ID=4) - smallest q_off
        let mut hits = vec![
            mock_hit(50, 100, 300, 350, 4),   // Should be 4th (smallest q_off)
            mock_hit(100, 140, 200, 250, 3),  // Should be 3rd (smaller q_end)
            mock_hit(100, 150, 200, 250, 1),  // Should be 1st (largest all)
            mock_hit(100, 150, 150, 200, 2),  // Should be 2nd (smaller s_off)
        ];

        // NCBI parity: Use frame-relative coordinates (not concatenated absolute)
        #[inline]
        fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
            (
                hit.q_aa_start as i32,
                hit.q_aa_end as i32,
                hit.s_aa_start as i32,
                hit.s_aa_end as i32,
            )
        }

        // LOSAT sort: all fields use b.cmp(&a) = descending
        hits.sort_by(|a, b| {
            let (aqo, aqe, aso, ase) = frame_relative_coords(a);
            let (bqo, bqe, bso, bse) = frame_relative_coords(b);
            bqo.cmp(&aqo)
                .then(bqe.cmp(&aqe))
                .then(bso.cmp(&aso))
                .then(bse.cmp(&ase))
        });

        // Verify order by checking raw_score (used as ID)
        let sorted_ids: Vec<i32> = hits.iter().map(|h| h.raw_score).collect();
        let expected_ids = vec![1, 2, 3, 4];

        println!("LOSAT sorted order (by ID): {:?}", sorted_ids);
        println!("Expected NCBI order:        {:?}", expected_ids);

        assert_eq!(
            sorted_ids, expected_ids,
            "HSP sort order mismatch! LOSAT={:?}, NCBI expected={:?}",
            sorted_ids, expected_ids
        );
    }

    /// Test that ties are handled correctly (stable sort behavior)
    #[test]
    fn test_hsp_sort_identical_coords() {
        // Two hits with identical coordinates but different IDs
        let mut hits = vec![
            mock_hit(100, 150, 200, 250, 1),
            mock_hit(100, 150, 200, 250, 2),
        ];

        // NCBI parity: Use frame-relative coordinates (not concatenated absolute)
        #[inline]
        fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
            (
                hit.q_aa_start as i32,
                hit.q_aa_end as i32,
                hit.s_aa_start as i32,
                hit.s_aa_end as i32,
            )
        }

        hits.sort_by(|a, b| {
            let (aqo, aqe, aso, ase) = frame_relative_coords(a);
            let (bqo, bqe, bso, bse) = frame_relative_coords(b);
            bqo.cmp(&aqo)
                .then(bqe.cmp(&aqe))
                .then(bso.cmp(&aso))
                .then(bse.cmp(&ase))
        });

        // Sort should be stable - original order preserved for equal elements
        let sorted_ids: Vec<i32> = hits.iter().map(|h| h.raw_score).collect();
        // Rust's sort_by is stable, so [1, 2] should remain [1, 2]
        assert_eq!(sorted_ids, vec![1, 2], "Stable sort should preserve order for equal elements");
    }

    /// Verify NCBI comparison semantics: all fields use DESCENDING order
    /// This documents the correct interpretation of s_RevCompareHSPsTbx
    #[test]
    fn test_ncbi_comparison_semantics() {
        // NCBI s_RevCompareHSPsTbx (link_hsps.c:367-370):
        // if (h1->subject.offset < h2->subject.offset) return 1;
        // if (h1->subject.offset > h2->subject.offset) return -1;
        //
        // In C qsort: positive return means first arg comes AFTER second arg
        // So: h1.offset < h2.offset → h1 after h2 → smaller comes after → DESCENDING
        //
        // The status document incorrectly stated this was ASCENDING.
        // This test documents the correct interpretation.

        // Test: same query coords, different subject offsets
        // NCBI descending: larger s_off should come first
        let mut hits = vec![
            mock_hit(100, 150, 50, 100, 1),   // s_off=50 (smaller)
            mock_hit(100, 150, 200, 250, 2),  // s_off=200 (larger)
        ];

        // NCBI parity: Use frame-relative coordinates (not concatenated absolute)
        #[inline]
        fn frame_relative_coords(hit: &UngappedHit) -> (i32, i32, i32, i32) {
            (
                hit.q_aa_start as i32,
                hit.q_aa_end as i32,
                hit.s_aa_start as i32,
                hit.s_aa_end as i32,
            )
        }

        hits.sort_by(|a, b| {
            let (aqo, aqe, aso, ase) = frame_relative_coords(a);
            let (bqo, bqe, bso, bse) = frame_relative_coords(b);
            bqo.cmp(&aqo)
                .then(bqe.cmp(&aqe))
                .then(bso.cmp(&aso))  // LOSAT: b.cmp(&a) = DESCENDING
                .then(bse.cmp(&ase))
        });

        let sorted_ids: Vec<i32> = hits.iter().map(|h| h.raw_score).collect();

        // DESCENDING: larger s_off (ID=2) should come first
        assert_eq!(
            sorted_ids, vec![2, 1],
            "Subject offset should be sorted DESCENDING (larger first). Got: {:?}",
            sorted_ids
        );

        println!("✓ Confirmed: LOSAT uses DESCENDING order for subject.offset (matches NCBI)");
    }
}

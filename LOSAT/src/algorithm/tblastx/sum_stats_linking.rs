//! NCBI-style sum-statistics HSP linking with lh_helper optimization
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c
//! 
//! Key optimizations from NCBI:
//! 1. `changed` flag: skip recomputation when previous link is still valid
//! 2. `linked_to` counter: track how many HSPs link to this one
//! 3. `path_changed` flag: skip recomputation when no chains were affected
//! 4. `use_current_max`: reuse max chains when they weren't affected by removal
//! 5. `next_larger`: skip HSPs with too-small sum
//! 6. Dual-index approach: index=0 (small gaps), index=1 (large gaps)

use rustc_hash::FxHashMap;
use rayon::prelude::*;
use crate::stats::sum_statistics::{
    gap_decay_divisor, small_gap_sum_e, large_gap_sum_e, normalize_score,
    defaults::{GAP_SIZE, OVERLAP_SIZE},
};
use crate::stats::KarlinParams;
use crate::stats::length_adjustment::compute_length_adjustment_simple;

use super::chaining::UngappedHit;
use super::diagnostics::diagnostics_enabled;
use super::lookup::QueryContext;

// NCBI groups by `(context/3)` (query strand + query index) and `SIGN(subject.frame)` (subject strand).
// For LOSAT multi-sequence support we explicitly include `(q_idx, s_idx)` and store only the strand signs.
//
// Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:522-525
type ContextKey = (u32, u32, i8, i8); // (q_idx, s_idx, q_strand, s_strand)

const WINDOW_SIZE: i32 = GAP_SIZE + OVERLAP_SIZE + 1;
const TRIM_SIZE: i32 = (OVERLAP_SIZE + 1) / 2;

/// Ordering method indices (NCBI ELinkOrderingMethod)
const LINK_SMALL_GAPS: usize = 0;
const LINK_LARGE_GAPS: usize = 1;

/// NCBI BLAST_GAP_PROB (ungapped search default)
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:66
const BLAST_GAP_PROB: f64 = 0.5;

/// NCBI BLAST_GAP_DECAY_RATE (ungapped search default)  
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:68
const BLAST_GAP_DECAY_RATE: f64 = 0.5;

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

/// Parameters for HSP linking cutoffs
/// These are the outputs from calculate_link_hsp_cutoffs_ncbi()
#[derive(Debug, Clone, Copy)]
pub struct LinkHspCutoffs {
    /// Cutoff for small gap linking (index 0)
    pub cutoff_small_gap: i32,
    /// Cutoff for big gap linking (index 1)
    pub cutoff_big_gap: i32,
    /// Gap probability (may be set to 0 for small search spaces)
    pub gap_prob: f64,
    /// Whether to ignore small gaps (when gap_prob=0)
    pub ignore_small_gaps: bool,
}

/// Parameters required for NCBI-style linking
/// These are computed once per subject and passed to linking
#[derive(Debug, Clone)]
pub struct LinkingParams {
    /// Average query length in AA (NCBI CalculateLinkHSPCutoffs formula)
    pub avg_query_length: i32,
    /// Subject length in nucleotides
    pub subject_len_nucl: i64,
    /// Minimum cutoff score across all contexts for this subject
    pub cutoff_score_min: i32,
    /// Scale factor (typically 1.0 for standard BLOSUM62)
    pub scale_factor: f64,
    /// Gap decay rate
    pub gap_decay_rate: f64,
}

impl Default for LinkingParams {
    fn default() -> Self {
        Self {
            avg_query_length: 100,
            subject_len_nucl: 300,
            cutoff_score_min: 0,
            scale_factor: 1.0,
            gap_decay_rate: BLAST_GAP_DECAY_RATE,
        }
    }
}

/// Compute NCBI-style average query length from nucleotide lengths
///
/// This replicates NCBI's SetupQueryInfo_OMF + CalculateLinkHSPCutoffs logic:
/// - For each query, compute 6 context lengths using BLAST_GetTranslatedProteinLength
/// - Build context offsets using s_QueryInfo_SetContext pattern
/// - Apply NCBI's average formula: (last_offset + last_length - 1) / num_contexts
///
/// Reference: 
/// - blast_util.c:BLAST_GetTranslatedProteinLength
/// - blast_setup_cxx.cpp:s_QueryInfo_SetContext  
/// - blast_parameters.c:1023-1026
///
/// ```c
/// query_length =
///     (query_info->contexts[query_info->last_context].query_offset +
///     query_info->contexts[query_info->last_context].query_length - 1)
///     / (query_info->last_context + 1);
/// ```
pub fn compute_avg_query_length_ncbi(query_nucl_lengths: &[usize]) -> i32 {
    if query_nucl_lengths.is_empty() {
        return 1;
    }
    
    // Build context offsets and lengths like NCBI
    // Each query has 6 contexts (frames 0-5 in NCBI, which map to frames 1,2,3,-1,-2,-3)
    let mut contexts: Vec<(i32, i32)> = Vec::new(); // (query_offset, query_length)
    
    for &nucl_len in query_nucl_lengths {
        for context in 0..6u32 {
            // NCBI BLAST_GetTranslatedProteinLength:
            // return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
            let prot_len = if nucl_len == 0 || nucl_len <= (context % 3) as usize {
                0
            } else {
                ((nucl_len - (context % 3) as usize) / 3) as i32
            };
            
            // NCBI s_QueryInfo_SetContext offset calculation:
            // Uint4 shift = prev_len ? prev_len + 1 : 0;
            // qinfo->contexts[index].query_offset = prev_loc + shift;
            let (prev_offset, prev_len) = contexts.last().copied().unwrap_or((0, 0));
            let shift = if prev_len > 0 { prev_len + 1 } else { 0 };
            let new_offset = prev_offset + shift;
            
            contexts.push((new_offset, prot_len));
        }
    }
    
    if contexts.is_empty() {
        return 1;
    }
    
    // NCBI average formula: (last_offset + last_length - 1) / num_contexts
    let (last_offset, last_length) = contexts.last().copied().unwrap();
    let num_contexts = contexts.len() as i32;
    
    let avg = (last_offset + last_length - 1) / num_contexts;
    avg.max(1)
}

/// NCBI CalculateLinkHSPCutoffs - verbatim port
///
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
///
/// ```c
/// void
/// CalculateLinkHSPCutoffs(EBlastProgramType program, BlastQueryInfo* query_info, 
///    const BlastScoreBlk* sbp, BlastLinkHSPParameters* link_hsp_params, 
///    const BlastInitialWordParameters* word_params,
///    Int8 db_length, Int4 subject_length)
/// {
///     Blast_KarlinBlk* kbp = NULL;
///     double gap_prob, gap_decay_rate, x_variable, y_variable;
///     Int4 expected_length, window_size, query_length;
///     Int8 search_sp;
///     const double kEpsilon = 1.0e-9;
///
///     if (!link_hsp_params)
///         return;
///
///     /* Get KarlinBlk for context with smallest lambda (still greater than zero) */
///     s_BlastFindSmallestLambda(sbp->kbp, query_info, &kbp);
///     if (!kbp)
///         return;
///
///     window_size
///         = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
///     gap_prob = link_hsp_params->gap_prob = BLAST_GAP_PROB;
///     gap_decay_rate = link_hsp_params->gap_decay_rate;
///     /* Use average query length */
///     
///     query_length =
///         (query_info->contexts[query_info->last_context].query_offset +
///         query_info->contexts[query_info->last_context].query_length - 1)
///         / (query_info->last_context + 1);
///     
///     if (Blast_SubjectIsTranslated(program) || program == eBlastTypeRpsTblastn) {
///         /* Lengths in subsequent calculations should be on the protein scale */
///         subject_length /= CODON_LENGTH;
///         db_length /= CODON_LENGTH;
///     }
///
///     
///     /* Subtract off the expected score. */
///    expected_length = (Int4)BLAST_Nint(log(kbp->K*((double) query_length)*
///                                     ((double) subject_length))/(kbp->H));
///    query_length = query_length - expected_length;
///
///    subject_length = subject_length - expected_length;
///    query_length = MAX(query_length, 1);
///    subject_length = MAX(subject_length, 1);
///
///    /* If this is a database search, use database length, else the single 
///       subject sequence length */
///    if (db_length > subject_length) {
///       y_variable = log((double) (db_length)/(double) subject_length)*(kbp->K)/
///          (gap_decay_rate);
///    } else {
///       y_variable = log((double) (subject_length + expected_length)/
///                        (double) subject_length)*(kbp->K)/(gap_decay_rate);
///    }
///
///    search_sp = ((Int8) query_length)* ((Int8) subject_length);
///    x_variable = 0.25*y_variable*((double) search_sp);
///
///    /* To use "small" gaps the query and subject must be "large" compared to
///       the gap size. If small gaps may be used, then the cutoff values must be
///       adjusted for the "bayesian" possibility that both large and small gaps 
///       are being checked for. */
///
///    if (search_sp > 8*window_size*window_size) {
///       x_variable /= (1.0 - gap_prob + kEpsilon);
///       link_hsp_params->cutoff_big_gap = 
///          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
///       x_variable = y_variable*(window_size*window_size);
///       x_variable /= (gap_prob + kEpsilon);
///       link_hsp_params->cutoff_small_gap = 
///          MAX(word_params->cutoff_score_min, 
///              (Int4) floor((log(x_variable)/kbp->Lambda)) + 1);
///    } else {
///       link_hsp_params->cutoff_big_gap = 
///          (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
///       /* The following is equivalent to forcing small gap rule to be ignored
///          when linking HSPs. */
///       link_hsp_params->gap_prob = 0;
///       link_hsp_params->cutoff_small_gap = 0;
///    }	
///
///    link_hsp_params->cutoff_big_gap *= (Int4)sbp->scale_factor;
///    link_hsp_params->cutoff_small_gap *= (Int4)sbp->scale_factor;
/// }
/// ```
pub fn calculate_link_hsp_cutoffs_ncbi(
    avg_query_length: i32,
    subject_len_nucl: i64,
    db_length: i64, // 0 for -subject mode
    cutoff_score_min: i32,
    scale_factor: f64,
    gap_decay_rate: f64,
    params: &KarlinParams,
) -> LinkHspCutoffs {
    const K_EPSILON: f64 = 1.0e-9;
    
    // NCBI: gap_prob = link_hsp_params->gap_prob = BLAST_GAP_PROB;
    let mut gap_prob = BLAST_GAP_PROB;
    
    // NCBI: window_size = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
    let window_size = WINDOW_SIZE;
    
    // NCBI: query_length already provided as avg_query_length
    let mut query_length = avg_query_length;
    
    // NCBI: if (Blast_SubjectIsTranslated(program)) subject_length /= CODON_LENGTH;
    // For tblastx, subject is translated, so divide by 3
    let mut subject_length = (subject_len_nucl / 3) as i32;
    let db_length = db_length / 3; // Also scale db_length for tblastx
    
    // NCBI: expected_length = (Int4)BLAST_Nint(log(kbp->K*q*s)/(kbp->H));
    let expected_length = blast_nint(
        (params.k * (query_length as f64) * (subject_length as f64)).ln() / params.h
    );
    
    // NCBI: query_length = query_length - expected_length;
    query_length -= expected_length;
    // NCBI: subject_length = subject_length - expected_length;
    subject_length -= expected_length;
    // NCBI: query_length = MAX(query_length, 1);
    query_length = query_length.max(1);
    // NCBI: subject_length = MAX(subject_length, 1);
    subject_length = subject_length.max(1);
    
    // NCBI: y_variable calculation with db_length > subject_length branch
    // For -subject mode, db_length = 0, so we always use the else branch
    let y_variable = if db_length > subject_length as i64 {
        // Database search mode
        ((db_length as f64) / (subject_length as f64)).ln() * params.k / gap_decay_rate
    } else {
        // Subject mode (db_length == 0 or single sequence)
        (((subject_length + expected_length) as f64) / (subject_length as f64)).ln()
            * params.k / gap_decay_rate
    };
    
    // NCBI: search_sp = ((Int8) query_length)* ((Int8) subject_length);
    let search_sp = (query_length as i64) * (subject_length as i64);
    
    // NCBI: x_variable = 0.25*y_variable*((double) search_sp);
    let mut x_variable = 0.25 * y_variable * (search_sp as f64);
    
    let window_sq = (window_size * window_size) as i64;
    
    let (cutoff_big_gap, cutoff_small_gap, ignore_small_gaps) = if search_sp > 8 * window_sq {
        // NCBI: Large search space - use both small and big gap rules
        // x_variable /= (1.0 - gap_prob + kEpsilon);
        x_variable /= 1.0 - gap_prob + K_EPSILON;
        
        // cutoff_big_gap = (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
        let cutoff_big = (x_variable.ln() / params.lambda).floor() as i32 + 1;
        
        // x_variable = y_variable*(window_size*window_size);
        let x_small = y_variable * (window_sq as f64);
        // x_variable /= (gap_prob + kEpsilon);
        let x_small = x_small / (gap_prob + K_EPSILON);
        
        // cutoff_small_gap = MAX(word_params->cutoff_score_min, floor(log(x)/Lambda)+1);
        let cutoff_small_raw = (x_small.ln() / params.lambda).floor() as i32 + 1;
        let cutoff_small = cutoff_score_min.max(cutoff_small_raw);
        
        (cutoff_big, cutoff_small, false)
    } else {
        // NCBI: Small search space - disable small gap rule
        // cutoff_big_gap = (Int4) floor((log(x_variable)/kbp->Lambda)) + 1;
        let cutoff_big = (x_variable.ln() / params.lambda).floor() as i32 + 1;
        
        // link_hsp_params->gap_prob = 0;
        gap_prob = 0.0;
        // link_hsp_params->cutoff_small_gap = 0;
        let cutoff_small = 0;
        
        (cutoff_big, cutoff_small, true)
    };
    
    // NCBI: cutoff_big_gap *= (Int4)sbp->scale_factor;
    // NCBI: cutoff_small_gap *= (Int4)sbp->scale_factor;
    let scale = scale_factor as i32;
    
    LinkHspCutoffs {
        cutoff_small_gap: cutoff_small_gap * scale,
        cutoff_big_gap: cutoff_big_gap * scale,
        gap_prob,
        ignore_small_gaps,
    }
}

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
    q_off_trim: i32,
    s_off_trim: i32,
    q_end_trim: i32,
    s_end_trim: i32,
    sum: [i32; 2],           // sum for both indices
    xsum: [f64; 2],          // normalized sum for both indices
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

/// Sentinel index indicating end of active list or "not in list"
const SENTINEL_IDX: usize = usize::MAX;

/// Apply NCBI-style sum-statistics even-gap linking
///
/// This is the main entry point for HSP linking with NCBI-compatible cutoffs.
///
/// # Arguments
/// * `hits` - Vector of ungapped HSPs to link
/// * `params` - Karlin-Altschul parameters for the scoring matrix
/// * `linking_params` - NCBI-style linking parameters (avg query length, subject length, etc.)
pub fn apply_sum_stats_even_gap_linking(
    hits: Vec<UngappedHit>,
    params: &KarlinParams,
    linking_params: &LinkingParams,
    query_contexts: &[QueryContext],
    subject_frame_bases: &[i32],
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

    // Group by (q_idx, s_idx, q_strand, s_strand) to match NCBI link_hsps.c
    // NCBI groups by context/3 (query strand) and SIGN(subject.frame) (subject strand)
    // Reference: link_hsps.c lines 524-525
    let mut groups: FxHashMap<ContextKey, Vec<UngappedHit>> = FxHashMap::default();
    for hit in hits {
        let q_strand: i8 = if hit.q_frame > 0 { 1 } else { -1 };
        let s_strand: i8 = if hit.s_frame > 0 { 1 } else { -1 };
        let key = (hit.q_idx, hit.s_idx, q_strand, s_strand);
        groups.entry(key).or_default().push(hit);
    }

    // Process each group in parallel - all groups use NCBI-compliant linking
    let results: Vec<Vec<UngappedHit>> = groups
        .into_par_iter()
        .map(|(_, group_hits)| {
            link_hsp_group_ncbi(
                group_hits,
                params,
                &cutoffs,
                linking_params.gap_decay_rate,
                diag_enabled,
                linking_params.subject_len_nucl,
                query_contexts,
                subject_frame_bases,
            )
        })
        .collect();

    results.into_iter().flatten().collect()
}

/// NCBI-style HSP linking for moderate-sized groups
///
/// # Arguments
/// * `group_hits` - HSPs in this group (same q_idx, s_idx, q_strand, s_strand)
/// * `params` - Karlin-Altschul parameters
/// * `cutoffs` - Pre-computed NCBI-style cutoffs for this subject
/// * `gap_decay_rate` - Gap decay rate for E-value calculation
/// * `diag_enabled` - Whether diagnostics are enabled
fn link_hsp_group_ncbi(
    mut group_hits: Vec<UngappedHit>,
    params: &KarlinParams,
    cutoffs: &LinkHspCutoffs,
    gap_decay_rate: f64,
    diag_enabled: bool,
    subject_len_nucl: i64,
    query_contexts: &[QueryContext],
    _subject_frame_bases: &[i32], // Kept for API compatibility; no longer used after frame-relative coord fix
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
    group_hits.sort_by(|a, b| {
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
    // NCBI link_hsps.c:560-571 - Effective length calculation for tblastx
    // ========================================================================
    // ```c
    // length_adjustment = query_info->contexts[query_context].length_adjustment;
    // query_length = query_info->contexts[query_context].query_length;
    // query_length = MAX(query_length - length_adjustment, 1);
    // subject_length = subject_length_orig; /* in nucleotides even for tblast[nx] */
    // /* If subject is translated, length adjustment is given in nucleotide
    //    scale. */
    // if (Blast_SubjectIsTranslated(program_number))  // tblastx = TRUE
    // {
    //    length_adjustment /= CODON_LENGTH;  // divide by 3
    //    subject_length /= CODON_LENGTH;
    // }
    // subject_length = MAX(subject_length - length_adjustment, 1);
    // ```
    //
    // Key insight: For tblastx, NCBI applies length_adjustment differently:
    //   - query: subtract full length_adjustment
    //   - subject: subtract (length_adjustment / 3) after converting to AA
    // ========================================================================
    let length_adjustment = compute_length_adjustment_simple(
        query_len_aa,
        subject_len_aa,
        params,
    ).length_adjustment;
    
    // NCBI: query_length = MAX(query_length - length_adjustment, 1)
    let eff_query_len = (query_len_aa - length_adjustment).max(1) as f64;
    
    // NCBI: For tblastx (Blast_SubjectIsTranslated = TRUE):
    //   length_adjustment /= CODON_LENGTH (integer division by 3)
    //   subject_length = MAX(subject_length - length_adjustment, 1)
    let length_adj_for_subject = length_adjustment / 3;  // integer division
    let eff_subject_len = (subject_len_aa - length_adj_for_subject).max(1) as f64;
    
    let eff_search_space = eff_query_len * eff_subject_len;
    
    // Use pre-computed cutoffs from NCBI algorithm
    let cutoff_small = cutoffs.cutoff_small_gap;
    let cutoff_big = cutoffs.cutoff_big_gap;
    let gap_prob = cutoffs.gap_prob;
    let ignore_small_gaps = cutoffs.ignore_small_gaps;
    
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
        eprintln!("  length_adjustment={}, length_adj_for_subject={}", 
            length_adjustment, length_adj_for_subject);
        eprintln!("  eff_query_len={:.2}, eff_subject_len={:.2}, eff_search_space={:.2e}", 
            eff_query_len, eff_subject_len, eff_search_space);
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

    // Initialize HspLink array (indices 0..n correspond to HSPs)
    // NCBI uses 1-based indexing with lh_helper[0] as sentinel
    // We use 0-based but handle the logic appropriately
    let mut hsp_links: Vec<HspLink> = group_hits.iter().enumerate().map(|(i, hit)| {
        // NCBI parity: Use frame-relative coordinates for trim calculations
        // Reference: link_hsps.c:545-549
        let (q_off, q_end, s_off, s_end) = frame_relative_coords(hit);
        let q_len_quarter = (q_end - q_off) / 4;
        let s_len_quarter = (s_end - s_off) / 4;
        let qt = TRIM_SIZE.min(q_len_quarter);
        let st = TRIM_SIZE.min(s_len_quarter);
        let score = hit.raw_score;
        let xsum_val = normalize_score(score, params.lambda, params.k.ln());
        
        HspLink {
            score,
            q_off_trim: q_off + qt,
            s_off_trim: s_off + st,
            q_end_trim: q_end - qt,
            s_end_trim: s_end - st,
            sum: [score - cutoff_small, score - cutoff_big],
            xsum: [xsum_val, xsum_val],
            num: [1, 1],
            link: [SENTINEL_IDX, SENTINEL_IDX],
            changed: true,
            linked_to: 0,
            start_of_chain: false, // Will be set to true for chain heads
            linked_set: false,     // Will be set to true for multi-HSP chains
            // Initialize intrusive linked list: 0 -> 1 -> 2 -> ... -> n-1 -> SENTINEL
            next_active: if i + 1 < n { i + 1 } else { SENTINEL_IDX },
            prev_active: if i > 0 { i - 1 } else { SENTINEL_IDX },
        }
    }).collect();
    
    // Head of active HSP list (first active HSP index, or SENTINEL_IDX if empty)
    let mut active_head: usize = if n > 0 { 0 } else { SENTINEL_IDX };

    // NCBI lh_helper array structure (lines 579-583):
    // lh_helper[0] = empty = additional end marker
    // lh_helper[1] = hp_start = empty entry  
    // lh_helper[i] = hsp_array[i-2] for i >= 2
    // So: hsp_links[j] <-> lh_helpers[j+2]
    //
    // We allocate n+2 entries and use direct index calculation:
    // - lh_helpers[0], lh_helpers[1]: sentinels
    // - lh_helpers[i+2]: corresponds to hsp_links[i]
    // NCBI `link_hsps.c` allocates `lh_helper` with `calloc`, so the sentinel
    // entries are zero-initialized. We mirror that exactly (sums=0, next_larger=0).
    // NCBI line 576: lh_helper[0].maxsum1 = -10000;
    let mut lh_helpers: Vec<LhHelper> = vec![
        LhHelper {
            hsp_idx: SENTINEL_IDX,
            q_off_trim: 0,
            s_off_trim: 0,
            sum: [0, 0],           // NCBI: calloc zero-initializes
            next_larger: 0,
            maxsum1: -10000,       // NCBI line 576
        };
        n + 2
    ];

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
            // Using intrusive linked list traversal (O(active) instead of O(n))
            // This matches NCBI's hp_start->next traversal where processed HSPs are removed
            if !ignore_small_gaps {
                let mut cur = active_head;
                while cur != SENTINEL_IDX {
                    if hsp_links[cur].sum[0] >= best_sum[0] {
                        best_sum[0] = hsp_links[cur].sum[0];
                        best[0] = Some(cur);
                    }
                    cur = hsp_links[cur].next_active;
                }
            }
            {
                let mut cur = active_head;
                while cur != SENTINEL_IDX {
                    if hsp_links[cur].sum[1] >= best_sum[1] {
                        best_sum[1] = hsp_links[cur].sum[1];
                        best[1] = Some(cur);
                    }
                    cur = hsp_links[cur].next_active;
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
                        if hsp_links[bi].linked_to == -1000 { 
                            use_current_max = false; 
                        } else {
                            let mut cur = bi;
                            loop {
                                let p = hsp_links[cur].link[0];
                                if p == SENTINEL_IDX {
                                    break;
                                }
                                if hsp_links[p].linked_to == -1000 {
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
                        if hsp_links[bi].linked_to == -1000 {
                            use_current_max = false;
                        } else {
                            let mut cur = bi;
                            loop {
                                let p = hsp_links[cur].link[1];
                                if p == SENTINEL_IDX {
                                    break;
                                }
                                if hsp_links[p].linked_to == -1000 {
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
            // lh_helpers[0] and [1] are sentinels; active HSPs occupy [2..lh_len).
            let mut lh_len = 2usize;
            let mut cur = active_head;
            // NCBI lines 669-671: track max sum[1] for maxsum1 calculation
            // Since LOSAT groups by frame, all HSPs have the same frame sign, so we track one max
            let mut running_max = -10000i32;
            while cur != SENTINEL_IDX {
                let hsp_idx = cur;

                // NCBI line 685: H->linked_to = 0
                hsp_links[hsp_idx].linked_to = 0;

                // NCBI lines 664-668: copy data + store ptr
                lh_helpers[lh_len].hsp_idx = hsp_idx;
                lh_helpers[lh_len].q_off_trim = hsp_links[hsp_idx].q_off_trim;
                lh_helpers[lh_len].s_off_trim = hsp_links[hsp_idx].s_off_trim;
                lh_helpers[lh_len].sum = [hsp_links[hsp_idx].sum[0], hsp_links[hsp_idx].sum[1]];

                // NCBI lines 669-671: update maxsum1
                running_max = running_max.max(hsp_links[hsp_idx].sum[1]);
                lh_helpers[lh_len].maxsum1 = running_max;

                // NCBI lines 675-684: compute next_larger
                // NCBI: while((cur_sum>=prev_sum) && (prev>0))
                let cur_sum = lh_helpers[lh_len].sum[1];
                let mut prev = lh_len - 1;
                let mut prev_sum = lh_helpers[prev].sum[1];
                while cur_sum >= prev_sum && prev > 0 {
                    prev = lh_helpers[prev].next_larger;
                    prev_sum = lh_helpers[prev].sum[1];
                }
                lh_helpers[lh_len].next_larger = prev;

                lh_len += 1;
                cur = hsp_links[cur].next_active;
            }
            // NCBI line 688: lh_helper[1].maxsum1 = -10000;
            lh_helpers[1].maxsum1 = -10000;
            let active_count = lh_len - 2;

            best = [None, None];
            best_sum = [-cutoff_small, -cutoff_big];

            // INDEX 0 LOOP (small gaps) - NCBI lines 691-768
            if !ignore_small_gaps {
                // NCBI: H_index = 2; for (H=hp_start->next; H!=NULL; H=H->next,H_index++)
                // Iterate over lh_helpers[2..] which contains only active HSPs
                for h_lh_idx in 2..lh_len {
                    let i = lh_helpers[h_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)
                    
                    let mut h_sum = 0i32;
                    let mut h_xsum = 0.0f64;
                    let mut h_num = 0i16;
                    let mut h_link: usize = SENTINEL_IDX;
                    
                    // NCBI line 702: if (H->hsp->score > cutoff[index])
                    if hsp_links[i].score > cutoff_small {
                        let h_qe = hsp_links[i].q_end_trim;
                        let h_se = hsp_links[i].s_end_trim;
                        let h_qe_gap = h_qe + WINDOW_SIZE;
                        let h_se_gap = h_se + WINDOW_SIZE;
                        
                        // Inner loop: NCBI lines 706-742
                        for j_lh_idx in (2..h_lh_idx).rev() {
                            let qo = lh_helpers[j_lh_idx].q_off_trim;
                            let so = lh_helpers[j_lh_idx].s_off_trim;
                            let sum = lh_helpers[j_lh_idx].sum[0];
                            
                            // NCBI line 717
                            if qo > h_qe_gap + TRIM_SIZE { break; }
                            // NCBI lines 719-724
                            if qo <= h_qe || so <= h_se || qo > h_qe_gap || so > h_se_gap { continue; }
                            
                            // NCBI line 727: if (sum>H_hsp_sum)
                            if sum > h_sum {
                                let j = lh_helpers[j_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)
                                h_num = hsp_links[j].num[0];
                                h_sum = hsp_links[j].sum[0];
                                h_xsum = hsp_links[j].xsum[0];
                                h_link = j;
                            }
                        }
                    }
                    
                    // NCBI lines 750-767: Update this HSP's link info
                    let score = hsp_links[i].score;
                    let new_sum = h_sum + (score - cutoff_small);
                    let new_xsum = h_xsum + normalize_score(score, params.lambda, params.k.ln());
                    
                    hsp_links[i].sum[0] = new_sum;
                    hsp_links[i].num[0] = h_num + 1;
                    hsp_links[i].link[0] = h_link;
                    hsp_links[i].xsum[0] = new_xsum;
                    lh_helpers[h_lh_idx].sum[0] = new_sum;
                    
                    // NCBI line 759-763: update best
                    if new_sum >= best_sum[0] {
                        best_sum[0] = new_sum;
                        best[0] = Some(i);
                    }
                    // NCBI line 765-766
                    if h_link != SENTINEL_IDX {
                        hsp_links[h_link].linked_to += 1;
                    }
                }
            }

            // INDEX 1 LOOP (large gaps) with next_larger optimization - NCBI lines 771-896
            for h_lh_idx in 2..lh_len {
                let i = lh_helpers[h_lh_idx].hsp_idx; // Use stored hsp_idx (NCBI: ptr)
                
                let mut h_sum = 0i32;
                let mut h_xsum = 0.0f64;
                let mut h_num = 0i16;
                let mut h_link: usize = SENTINEL_IDX;
                
                // NCBI line 781: H->hsp_link.changed=1
                hsp_links[i].changed = true;
                let prev_link = hsp_links[i].link[1];
                
                // NCBI `link_hsps.c` lines 781-795:
                //   H->hsp_link.changed=1;
                //   H2 = H->hsp_link.link[index];
                //   if ((!first_pass) && ((H2==0) || (H2->hsp_link.changed==0))) { ... }
                //
                // This is the original fast-path: if the previous best choice exists and
                // was not changed in the last pass, it is still the best.
                // NCBI checks only changed==0, NOT linked_to>=0 here.
                let can_skip_ncbi = !first_pass
                    && (prev_link == SENTINEL_IDX || !hsp_links[prev_link].changed);
                
                if can_skip_ncbi {
                    if prev_link != SENTINEL_IDX {
                        h_num = hsp_links[prev_link].num[1];
                        h_sum = hsp_links[prev_link].sum[1];
                        h_xsum = hsp_links[prev_link].xsum[1];
                    }
                    h_link = prev_link;
                    hsp_links[i].changed = false;
                } else if hsp_links[i].score > cutoff_big {
                    let h_qe = hsp_links[i].q_end_trim;
                    let h_se = hsp_links[i].s_end_trim;
                    
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
                    if !first_pass && prev_link != SENTINEL_IDX && hsp_links[prev_link].linked_to >= 0 {
                        h_sum = hsp_links[prev_link].sum[1] - 1;
                    }
                    
                    // Inner loop (NCBI lines 827-861)
                    // CRITICAL: Must save current_idx BEFORE decrement to access correct hsp_idx
                    let mut j_lh_idx = h_lh_idx - 1;
                    while j_lh_idx > 1 {
                        // NCBI: H2_helper points to lh_helper[H2_index] at loop START
                        let current_idx = j_lh_idx; // Save before any modification
                        let sum = lh_helpers[current_idx].sum[1];
                        let next_larger = lh_helpers[current_idx].next_larger;
                        let qo = lh_helpers[current_idx].q_off_trim;
                        let so = lh_helpers[current_idx].s_off_trim;
                        
                        let b0 = sum <= h_sum;
                        
                        // NCBI line 841: H2_index--
                        j_lh_idx -= 1;
                        if b0 {
                            j_lh_idx = next_larger;
                        }
                        
                        let b1 = qo <= h_qe;
                        let b2 = so <= h_se;
                        
                        // NCBI line 852: if (!(b0|b1|b2))
                        if !(b0 || b1 || b2) {
                            // Use saved current_idx (NCBI: H2 = H2_helper->ptr)
                            let j = lh_helpers[current_idx].hsp_idx;
                            h_num = hsp_links[j].num[1];
                            h_sum = hsp_links[j].sum[1];
                            h_xsum = hsp_links[j].xsum[1];
                            h_link = j;
                        }
                    }
                }
                
                // NCBI lines 863-895: Update this HSP's link info
                let score = hsp_links[i].score;
                let new_sum = h_sum + (score - cutoff_big);
                let new_xsum = h_xsum + normalize_score(score, params.lambda, params.k.ln());
                
                hsp_links[i].sum[1] = new_sum;
                hsp_links[i].num[1] = h_num + 1;
                hsp_links[i].link[1] = h_link;
                hsp_links[i].xsum[1] = new_xsum;
                lh_helpers[h_lh_idx].sum[1] = new_sum;
                
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
                let mut prev_sum = lh_helpers[prev].sum[1];
                while cur_sum >= prev_sum && prev > 0 {
                    prev = lh_helpers[prev].next_larger;
                    prev_sum = lh_helpers[prev].sum[1];
                }
                lh_helpers[h_lh_idx].next_larger = prev;
                
                // NCBI line 874: lh_helper[H_index].maxsum1 = MAX(lh_helper[H_index-1].maxsum1, new_sum);
                // Note: maxsum1 is unused (NCBI line 850 is if(0)) but we compute it for strict parity
                lh_helpers[h_lh_idx].maxsum1 = lh_helpers[h_lh_idx - 1].maxsum1.max(new_sum);
                
                if new_sum >= best_sum[1] {
                    best_sum[1] = new_sum;
                    best[1] = Some(i);
                }
                if h_link != SENTINEL_IDX {
                    hsp_links[h_link].linked_to += 1;
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
                let num = hsp_links[bi].num[0] as usize;
                let xsum = hsp_links[bi].xsum[0];
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
        }
        
        if let Some(bi) = best[1] {
            let num = hsp_links[bi].num[1] as usize;
            let xsum = hsp_links[bi].xsum[1];
            let divisor = gap_decay_divisor(gap_decay_rate, num);
            prob[1] = large_gap_sum_e(num as i16, xsum,
                eff_query_len as i32,
                eff_subject_len as i32,
                eff_search_space as i64, divisor);
            // NCBI `link_hsps.c` lines 931-935:
            //   if( num > 1 ) {
            //     if( 1 - gap_prob == 0 || (prob /= 1 - gap_prob) > INT4_MAX ) prob = INT4_MAX;
            //   }
            if num > 1 {
                let denom = 1.0 - gap_prob;
                if denom == 0.0 || (prob[1] / denom) > int4_max {
                    prob[1] = int4_max;
                } else {
                    prob[1] /= denom;
                }
            }
            if diag_enabled {
                // Debug: log the first few E-value calculations
                static DEBUG_EVALUE_COUNT: std::sync::atomic::AtomicUsize =
                    std::sync::atomic::AtomicUsize::new(0);
                let count = DEBUG_EVALUE_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if count < 3 {
                    eprintln!(
                        "[DEBUG] E-value calc #{}: score={}, num={}, xsum={:.4}, divisor={:.4}, e={:.4e}",
                        count,
                        hsp_links[bi].score,
                        num,
                        xsum,
                        divisor,
                        prob[1]
                    );
                }
            }
        }

        // NCBI line 940: Select best ordering method
        let ordering = if !ignore_small_gaps && prob[0] <= prob[1] { 0 } else { 1 };
        // NCBI never has "best==NULL" here as long as `number_of_hsps > 0`:
        // the DP loops always assign `best[index]` while scanning the list.
        let best_i = if let Some(i) = best[ordering] {
            i
        } else {
            // Fallback for safety; should be unreachable with NCBI-equivalent logic.
            best[1 - ordering].expect("NCBI parity: best must exist when remaining>0")
        };
        
        let evalue = prob[ordering];
        
        if diag_enabled {
            // Debug: count chain statistics
            static CHAIN_STATS: std::sync::atomic::AtomicBool =
                std::sync::atomic::AtomicBool::new(false);
            if !CHAIN_STATS.swap(true, std::sync::atomic::Ordering::Relaxed) {
                // Count chain size for first group only
                let chain_len = hsp_links[best_i].num[ordering];
                eprintln!(
                    "[DEBUG] First chain: len={}, e-value={:.2e}, ordering={}",
                    chain_len, evalue, ordering
                );
            }
        }

        // NCBI lines 964-968: If best has linked_to > 0, set path_changed
        // This means some other HSP was linking to this chain
        if hsp_links[best_i].linked_to > 0 {
            path_changed = true;
        }

        // NCBI lines 959-963: Determine if this is a linked set (multi-HSP chain)
        // If best has a link, it's a linked set (multiple HSPs)
        let linked_set = hsp_links[best_i].link[ordering] != SENTINEL_IDX;
        
        // Mark chain head and remove chain from consideration (NCBI lines 961-1000)
        // IMPORTANT: Set the chain's E-value to ALL members of the chain,
        // not just the head. This matches NCBI's behavior where all HSPs
        // in a chain share the same E-value.
        let mut is_first = true;
        let mut cur = best_i;
        loop {
            // NCBI line 968: if (H->linked_to>1) path_changed=1
            if hsp_links[cur].linked_to > 1 {
                path_changed = true;
            }
            hsp_links[cur].linked_to = -1000;
            hsp_links[cur].changed = true;
            
            // === O(1) UNLINK from active list (NCBI: physically remove from linked list) ===
            {
                let prev_idx = hsp_links[cur].prev_active;
                let next_idx = hsp_links[cur].next_active;
                
                if prev_idx != SENTINEL_IDX {
                    hsp_links[prev_idx].next_active = next_idx;
                } else {
                    // cur was the head, update active_head
                    active_head = next_idx;
                }
                
                if next_idx != SENTINEL_IDX {
                    hsp_links[next_idx].prev_active = prev_idx;
                }
                
                // Mark as not in list
                hsp_links[cur].next_active = SENTINEL_IDX;
                hsp_links[cur].prev_active = SENTINEL_IDX;
            }
            
            // NCBI line 972: H->linked_set = linked_set
            hsp_links[cur].linked_set = linked_set;
            
            // NCBI line 973: H->ordering_method = ordering_method;
            // This records which ordering method (small gap=0 / large gap=1) was used
            // to select this chain for output.
            group_hits[cur].ordering_method = ordering as u8;
            
            // Set E-value for ALL HSPs in the chain (not just head)
            // This is key: chain members inherit the chain's E-value
            group_hits[cur].e_value = evalue;
            
            if is_first {
                // Chain head: mark as start_of_chain (NCBI line 955)
                hsp_links[cur].start_of_chain = true;
                is_first = false;
            } else {
                // Non-head HSP in chain: mark as not start_of_chain
                hsp_links[cur].start_of_chain = false;
            }
            remaining -= 1;

            let next = hsp_links[cur].link[ordering];
            if next == SENTINEL_IDX {
                break;
            }
            cur = next;
        }
    }

    // NCBI behavior: All HSPs are kept in the array, but chain information
    // (num, xsum, evalue) is used for filtering and reporting.
    group_hits
}

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
            ordering_method: 0,
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


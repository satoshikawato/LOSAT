//! HSP extraction and coordinate adjustment for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2384-2392, 4719-4775
//!
//! This module contains InitHSP (initial HSP with absolute coordinates) and
//! functions for converting to UngappedHit (context-relative coordinates).

use super::chaining::UngappedHit;
use super::extension::convert_coords;
use super::lookup::QueryContext;
use super::translation::QueryFrame;
use super::tracing::{trace_hsp_target, trace_match_target, trace_ungapped_hit_if_match};

/// NCBI BlastInitHSP equivalent - stores initial HSP with absolute coordinates
///
/// Reference: blast_hits.h:BlastInitHSP, blast_extend.c:360-375 BlastSaveInitHsp
///
/// This structure stores HSPs after extension but before coordinate conversion.
/// Coordinates are stored as absolute positions in the concatenated buffer.
#[derive(Clone, Copy)]
pub struct InitHSP {
    /// Query absolute coordinate (concatenated buffer, 0-based)
    /// Reference: blast_query_info.c:311-315, blast_util.c:112-116.
    pub q_start_absolute: i32,
    /// Query end absolute coordinate (concatenated buffer, 0-based)
    /// Reference: blast_query_info.c:311-315, blast_util.c:112-116.
    pub q_end_absolute: i32,
    /// Subject coordinate (frame-relative, 0-based)
    /// Reference: blast_aascan.c:110-113.
    pub s_start: i32,
    /// Subject end coordinate (frame-relative, 0-based)
    /// Reference: blast_aascan.c:110-113.
    pub s_end: i32,
    /// Raw score from extension
    pub score: i32,
    /// Query context index
    pub ctx_idx: usize,
    /// Subject frame index
    pub s_f_idx: usize,
    /// Query record index
    pub q_idx: u32,
    /// Subject record index
    pub s_idx: u32,
    /// Query frame
    pub q_frame: i8,
    /// Subject frame
    pub s_frame: i8,
    /// Query original length
    pub q_orig_len: usize,
    /// Subject original length
    pub s_orig_len: usize,
}

/// Trace an InitHSP if it matches the trace target.
///
/// This function converts absolute coordinates to nucleotide coordinates
/// for comparison with the trace target.
#[inline]
pub fn trace_init_hsp_if_match(stage: &str, init: &InitHSP, contexts: &[QueryContext]) {
    let Some(target) = trace_hsp_target() else {
        return;
    };

    let ctx = &contexts[init.ctx_idx];
    // Convert absolute query coords to context-relative (0-based).
    // Reference: blast_gapalign.c:4756-4768, blast_query_info.c:311-315.
    let q_start_rel = adjust_initial_hsp_offsets(init.q_start_absolute, ctx.frame_base);
    let q_end_rel = adjust_initial_hsp_offsets(init.q_end_absolute, ctx.frame_base);
    if q_start_rel < 0 || q_end_rel < 0 || init.s_start < 0 || init.s_end < 0 {
        return;
    }

    // Convert to logical AA coords (0-based, half-open).
    let q_aa_start = q_start_rel as usize;
    let q_aa_end = q_end_rel as usize;
    let s_aa_start = init.s_start as usize;
    let s_aa_end = init.s_end as usize;

    let (q_start, q_end) = convert_coords(q_aa_start, q_aa_end, init.q_frame, init.q_orig_len);
    let (s_start, s_end) = convert_coords(s_aa_start, s_aa_end, init.s_frame, init.s_orig_len);

    if trace_match_target(target, q_start, q_end, s_start, s_end) {
        eprintln!(
            "[TRACE_HSP] stage={} raw_score={} ctx_idx={} s_f_idx={} q_frame={} s_frame={} q={}-{} s={}-{}",
            stage,
            init.score,
            init.ctx_idx,
            init.s_f_idx,
            init.q_frame,
            init.s_frame,
            q_start,
            q_end,
            s_start,
            s_end
        );
    }
}

/// NCBI s_AdjustInitialHSPOffsets equivalent
///
/// Reference: blast_gapalign.c:2384-2392
///
/// NCBI code:
/// ```c
/// static NCBI_INLINE void
/// s_AdjustInitialHSPOffsets(BlastInitHSP* init_hsp, Int4 query_start)
/// {
///     init_hsp->offsets.qs_offsets.q_off -= query_start;
///     if (init_hsp->ungapped_data) {
///         init_hsp->ungapped_data->q_start -= query_start;
///     }
///     ASSERT(init_hsp->ungapped_data == NULL ||
///            init_hsp->ungapped_data->q_start >= 0);
/// }
/// ```
///
/// Converts absolute coordinates to context-relative coordinates.
/// This is called in BLAST_GetUngappedHSPList (blast_gapalign.c:4756-4758).
#[inline]
pub fn adjust_initial_hsp_offsets(
    hsp_q_absolute: i32, // Absolute coordinate in concatenated buffer
    frame_base: i32,     // Context start position (absolute)
) -> i32 {
    // NCBI: init_hsp->ungapped_data->q_start -= query_start;
    hsp_q_absolute - frame_base
}

/// NCBI BLAST_GetUngappedHSPList equivalent
///
/// Reference: blast_gapalign.c:4719-4775
///
/// NCBI code flow:
/// 1. Get context for each InitHSP (s_GetUngappedHSPContext)
/// 2. Adjust coordinates (s_AdjustInitialHSPOffsets)
/// 3. Initialize HSP (Blast_HSPInit)
/// 4. Add to hsp_list
///
/// This function converts InitHSPs (with absolute coordinates) to UngappedHits
/// (with context-relative coordinates). Reevaluation is performed separately
/// by Blast_HSPListReevaluateUngapped equivalent.
///
/// NCBI reference (verbatim, blast_gapalign.c:4756-4768):
/// ```c
/// context = s_GetUngappedHSPContext(query_info, init_hsp);
/// s_AdjustInitialHSPOffsets(init_hsp, query_info->contexts[context].query_offset);
/// ungapped_data = init_hsp->ungapped_data;
/// Blast_HSPInit(ungapped_data->q_start,
///               ungapped_data->length+ungapped_data->q_start,
///               ungapped_data->s_start,
///               ungapped_data->length+ungapped_data->s_start,
///               init_hsp->offsets.qs_offsets.q_off,
///               init_hsp->offsets.qs_offsets.s_off,
///               context, query_info->contexts[context].frame,
///               subject->frame, ungapped_data->score, NULL, &new_hsp);
/// Blast_HSPListSaveHSP(hsp_list, new_hsp);
/// ```
pub fn get_ungapped_hsp_list(
    init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
) -> Vec<UngappedHit> {
    let mut ungapped_hits = Vec::new();

    for init_hsp in init_hsps {
        let ctx = &contexts[init_hsp.ctx_idx];
        let _s_frame = &s_frames[init_hsp.s_f_idx];

        // NCBI: s_GetUngappedHSPContext equivalent
        // Context is already stored in init_hsp.ctx_idx

        // NCBI: s_AdjustInitialHSPOffsets (blast_gapalign.c:2384-2392)
        // Convert absolute coordinates to context-relative coordinates
        // NCBI: init_hsp->ungapped_data->q_start -= query_start;
        // where query_start = query_info->contexts[context].query_offset
        let q_start_relative = adjust_initial_hsp_offsets(init_hsp.q_start_absolute, ctx.frame_base);
        let q_end_relative = adjust_initial_hsp_offsets(init_hsp.q_end_absolute, ctx.frame_base);

        // NCBI: Blast_HSPInit (blast_hits.c:150-189)
        // query.offset = ungapped_data->q_start (context-relative, 0-based)
        // query.end = ungapped_data->length + ungapped_data->q_start
        //
        // In LOSAT, q_start_relative is context-relative coordinate (0-based).
        let qs = q_start_relative as usize;
        let qe = q_end_relative as usize;
        let ss = init_hsp.s_start as usize;
        let se = init_hsp.s_end as usize;
        let len_u = qe.saturating_sub(qs);

        if len_u == 0 {
            continue;
        }

        // Offsets are already 0-based in query/subject->sequence buffers.
        // Reference: blast_gapalign.c:4756-4768, blast_aascan.c:110-113.
        let (qs_l, qe_l) = (qs, qe);
        let (ss_l, se_l) = (ss, se);

        // Create UngappedHit with context-relative coordinates (before reevaluation)
        // NCBI: Blast_HSPInit creates HSP with original extension score
        let uh = UngappedHit {
            q_idx: init_hsp.q_idx,
            s_idx: init_hsp.s_idx,
            ctx_idx: init_hsp.ctx_idx,
            s_f_idx: init_hsp.s_f_idx,
            q_frame: init_hsp.q_frame,
            s_frame: init_hsp.s_frame,
            q_aa_start: qs_l,
            q_aa_end: qe_l,
            s_aa_start: ss_l,
            s_aa_end: se_l,
            q_orig_len: init_hsp.q_orig_len,
            s_orig_len: init_hsp.s_orig_len,
            raw_score: init_hsp.score, // Original extension score (before reevaluation)
            e_value: f64::INFINITY,
            num_ident: 0, // Will be computed during reevaluation
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        };
        trace_ungapped_hit_if_match("after_get_ungapped_hsp_list", &uh);
        ungapped_hits.push(uh);
    }

    // NCBI: Sort the HSP array by score
    // Reference: blast_gapalign.c:4772
    ungapped_hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));

    ungapped_hits
}

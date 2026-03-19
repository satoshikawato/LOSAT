//! HSP extraction and coordinate adjustment for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2384-2392, 4719-4775
//!
//! This module contains InitHSP (initial HSP with absolute coordinates) and
//! functions for converting to UngappedHit (context-relative coordinates).

use super::chaining::UngappedHit;
use super::extension::convert_coords;
use super::lookup::QueryContext;
use super::tracing::{trace_hsp_target, trace_match_target, trace_ungapped_hit_if_match};
use super::translation::QueryFrame;
use std::cmp::Ordering;
use std::ffi::c_void;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1366-1377
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       ScoreCompareHSPs);
// ```
unsafe extern "C" {
    fn qsort(
        base: *mut c_void,
        nmemb: usize,
        size: usize,
        compar: Option<unsafe extern "C" fn(*const c_void, *const c_void) -> i32>,
    );
}

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
    /// Seed query offset stored in `offsets.qs_offsets.q_off`.
    /// Reference: blast_extend.c:351-370, blast_gapalign.c:4756-4768.
    pub q_seed_absolute: i32,
    /// Seed subject offset stored in `offsets.qs_offsets.s_off`.
    /// Reference: blast_extend.c:351-370, blast_gapalign.c:4756-4768.
    pub s_seed: i32,
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
        let q_start_relative =
            adjust_initial_hsp_offsets(init_hsp.q_start_absolute, ctx.frame_base);
        let q_end_relative = adjust_initial_hsp_offsets(init_hsp.q_end_absolute, ctx.frame_base);
        let q_seed_relative = adjust_initial_hsp_offsets(init_hsp.q_seed_absolute, ctx.frame_base);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2384-2392
        // ```c
        // ASSERT(init_hsp->ungapped_data == NULL ||
        //        init_hsp->ungapped_data->q_start >= 0);
        // ```
        assert!(
            q_start_relative >= 0,
            "NCBI BLAST requires context-relative ungapped q_start >= 0"
        );
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4760-4767
        // ```c
        // Blast_HSPInit(ungapped_data->q_start,
        //               ungapped_data->length+ungapped_data->q_start,
        //               ungapped_data->s_start,
        //               ungapped_data->length+ungapped_data->s_start,
        //               init_hsp->offsets.qs_offsets.q_off,
        //               init_hsp->offsets.qs_offsets.s_off, ...);
        // ```
        assert!(
            q_end_relative > q_start_relative,
            "NCBI BLAST requires positive ungapped query span"
        );
        assert!(
            q_seed_relative >= 0,
            "NCBI BLAST requires non-negative ungapped query seed offset"
        );
        assert!(
            init_hsp.s_start >= 0,
            "NCBI BLAST requires non-negative ungapped subject start"
        );
        assert!(
            init_hsp.s_end > init_hsp.s_start,
            "NCBI BLAST requires positive ungapped subject span"
        );
        assert!(
            init_hsp.s_seed >= 0,
            "NCBI BLAST requires non-negative ungapped subject seed offset"
        );

        // NCBI: Blast_HSPInit (blast_hits.c:150-189)
        // query.offset = ungapped_data->q_start (context-relative, 0-based)
        // query.end = ungapped_data->length + ungapped_data->q_start
        //
        // In LOSAT, q_start_relative is context-relative coordinate (0-based).
        let qs = usize::try_from(q_start_relative)
            .expect("NCBI BLAST requires non-negative ungapped query start");
        let qe = usize::try_from(q_end_relative)
            .expect("NCBI BLAST requires non-negative ungapped query end");
        let ss = usize::try_from(init_hsp.s_start)
            .expect("NCBI BLAST requires non-negative ungapped subject start");
        let se = usize::try_from(init_hsp.s_end)
            .expect("NCBI BLAST requires non-negative ungapped subject end");
        let q_seed = usize::try_from(q_seed_relative)
            .expect("NCBI BLAST requires non-negative ungapped query seed");
        let s_seed = usize::try_from(init_hsp.s_seed)
            .expect("NCBI BLAST requires non-negative ungapped subject seed");

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
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4756-4768
            // ```c
            // Blast_HSPInit(...,
            //               init_hsp->offsets.qs_offsets.q_off,
            //               init_hsp->offsets.qs_offsets.s_off, ...);
            // ```
            q_seed_off: q_seed,
            s_seed_off: s_seed,
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4719-4775
    // ```c
    // Blast_HSPListSaveHSP(hsp_list, new_hsp);
    // ...
    // Blast_HSPListSortByScore(hsp_list);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
    // ```c
    // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
    //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
    //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
    //     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
    //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1349-1369
    // ```c
    // if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //           ScoreCompareHSPs);
    // }
    // ```
    if !ungapped_hits_is_sorted_by_score_ncbi(&ungapped_hits) {
        ncbi_qsort_ungapped_hits_by_score(&mut ungapped_hits);
    }
    ungapped_hits
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2)
// {
//     ...
// }
// ```
fn score_compare_ungapped_hits_ncbi(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
        .then_with(|| b.s_aa_end.cmp(&a.s_aa_end))
        .then_with(|| a.q_aa_start.cmp(&b.q_aa_start))
        .then_with(|| b.q_aa_end.cmp(&a.q_aa_end))
}

#[inline]
fn ncbi_qsort_ordering(ordering: Ordering) -> i32 {
    match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1377
// ```c
// if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//           ScoreCompareHSPs);
// }
// ```
unsafe extern "C" fn score_compare_ungapped_hits_qsort(
    left: *const c_void,
    right: *const c_void,
) -> i32 {
    let left = unsafe { &*(left as *const UngappedHit) };
    let right = unsafe { &*(right as *const UngappedHit) };
    ncbi_qsort_ordering(score_compare_ungapped_hits_ncbi(left, right))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1366-1377
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       ScoreCompareHSPs);
// ```
fn ncbi_qsort_ungapped_hits_by_score(hits: &mut [UngappedHit]) {
    if hits.len() <= 1 {
        return;
    }
    debug_assert!(!std::mem::needs_drop::<UngappedHit>());
    unsafe {
        qsort(
            hits.as_mut_ptr().cast::<c_void>(),
            hits.len(),
            std::mem::size_of::<UngappedHit>(),
            Some(score_compare_ungapped_hits_qsort),
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1355-1369
// ```c
// Boolean Blast_HSPListIsSortedByScore(const BlastHSPList* hsp_list)
// {
//     ...
//     if (ScoreCompareHSPs(&hsp_list->hsp_array[index],
//                          &hsp_list->hsp_array[index+1]) > 0) {
//         return FALSE;
//     }
// }
// ```
fn ungapped_hits_is_sorted_by_score_ncbi(hits: &[UngappedHit]) -> bool {
    if hits.len() <= 1 {
        return true;
    }

    for index in 0..hits.len() - 1 {
        if score_compare_ungapped_hits_ncbi(&hits[index], &hits[index + 1]) == Ordering::Greater {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::KarlinParams;

    fn make_context() -> QueryContext {
        QueryContext {
            q_idx: 0,
            f_idx: 0,
            frame: 1,
            aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
            aa_seq_nomask: None,
            aa_len: 5,
            orig_len: 15,
            frame_base: 0,
            karlin_params: KarlinParams::default(),
        }
    }

    fn make_subject_frame() -> QueryFrame {
        QueryFrame {
            frame: 1,
            aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
            aa_seq_nomask: None,
            aa_len: 5,
            orig_len: 15,
            seg_masks: Vec::new(),
        }
    }

    fn make_init_hsp(q_start: i32, q_end: i32, s_start: i32, s_end: i32, score: i32) -> InitHSP {
        InitHSP {
            q_start_absolute: q_start,
            q_end_absolute: q_end,
            s_start,
            s_end,
            q_seed_absolute: q_start,
            s_seed: s_start,
            score,
            ctx_idx: 0,
            s_f_idx: 0,
            q_idx: 0,
            s_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_orig_len: 15,
            s_orig_len: 15,
        }
    }

    #[test]
    fn test_get_ungapped_hsp_list_sorts_like_ncbi_score_compare_match() {
        let contexts = vec![make_context()];
        let subject_frames = vec![make_subject_frame()];
        let hits = get_ungapped_hsp_list(
            vec![
                make_init_hsp(4, 10, 30, 36, 100),
                make_init_hsp(2, 10, 20, 28, 100),
                make_init_hsp(1, 8, 20, 27, 100),
                make_init_hsp(3, 11, 20, 28, 120),
            ],
            &contexts,
            &subject_frames,
        );

        let ordered: Vec<(i32, usize, usize, usize)> = hits
            .iter()
            .map(|hit| {
                (
                    hit.raw_score,
                    hit.s_aa_start,
                    hit.q_aa_end.saturating_sub(hit.q_aa_start),
                    hit.q_aa_start,
                )
            })
            .collect();

        assert_eq!(
            ordered,
            vec![
                (120, 20, 8, 3),
                (100, 20, 8, 2),
                (100, 20, 7, 1),
                (100, 30, 6, 4)
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1355-1369
    // ```c
    // if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //     qsort(..., ScoreCompareHSPs);
    // }
    // ```
    #[test]
    fn test_get_ungapped_hsp_list_preserves_insertion_order_when_score_sorted() {
        let contexts = vec![
            QueryContext {
                q_idx: 0,
                f_idx: 0,
                frame: 1,
                aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
                aa_seq_nomask: None,
                aa_len: 5,
                orig_len: 15,
                frame_base: 0,
                karlin_params: KarlinParams::default(),
            },
            QueryContext {
                q_idx: 1,
                f_idx: 0,
                frame: 1,
                aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
                aa_seq_nomask: None,
                aa_len: 5,
                orig_len: 15,
                frame_base: 100,
                karlin_params: KarlinParams::default(),
            },
        ];
        let subject_frames = vec![make_subject_frame()];
        let hits = get_ungapped_hsp_list(
            vec![
                InitHSP {
                    q_start_absolute: 101,
                    q_end_absolute: 106,
                    s_start: 20,
                    s_end: 25,
                    q_seed_absolute: 101,
                    s_seed: 20,
                    score: 100,
                    ctx_idx: 1,
                    s_f_idx: 0,
                    q_idx: 1,
                    s_idx: 0,
                    q_frame: 1,
                    s_frame: 1,
                    q_orig_len: 15,
                    s_orig_len: 15,
                },
                InitHSP {
                    q_start_absolute: 1,
                    q_end_absolute: 6,
                    s_start: 20,
                    s_end: 25,
                    q_seed_absolute: 1,
                    s_seed: 20,
                    score: 100,
                    ctx_idx: 0,
                    s_f_idx: 0,
                    q_idx: 0,
                    s_idx: 0,
                    q_frame: 1,
                    s_frame: 1,
                    q_orig_len: 15,
                    s_orig_len: 15,
                },
            ],
            &contexts,
            &subject_frames,
        );

        let ordered_qidx: Vec<u32> = hits.iter().map(|hit| hit.q_idx).collect();
        assert_eq!(ordered_qidx, vec![1, 0]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2384-2392
    // ```c
    // ASSERT(init_hsp->ungapped_data == NULL ||
    //        init_hsp->ungapped_data->q_start >= 0);
    // ```
    #[test]
    #[should_panic(expected = "NCBI BLAST requires context-relative ungapped q_start >= 0")]
    fn test_get_ungapped_hsp_list_fails_fast_on_negative_context_relative_query_start() {
        let contexts = vec![QueryContext {
            q_idx: 0,
            f_idx: 0,
            frame: 1,
            aa_seq: vec![0, 1, 1, 1, 1, 1, 0],
            aa_seq_nomask: None,
            aa_len: 5,
            orig_len: 15,
            frame_base: 10,
            karlin_params: KarlinParams::default(),
        }];
        let subject_frames = vec![make_subject_frame()];
        let invalid = InitHSP {
            q_start_absolute: 5,
            q_end_absolute: 10,
            s_start: 20,
            s_end: 25,
            q_seed_absolute: 5,
            s_seed: 20,
            score: 100,
            ctx_idx: 0,
            s_f_idx: 0,
            q_idx: 0,
            s_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_orig_len: 15,
            s_orig_len: 15,
        };

        let _ = get_ungapped_hsp_list(vec![invalid], &contexts, &subject_frames);
    }
}

//! `blast_kappa.c`-style post-processing helpers for the default `blastp`
//! parity path.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c
//! Reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c

use anyhow::Result;

use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::core::blast_seg::{SegMasker, SegParams};
use crate::core::blast_stat::{e_to_p, p_to_e};
use crate::core::composition_adjustment::adjust_scores::{
    blast_overall_p_value, read_aa_composition, AdjustedProteinMatrix,
};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, blast_redo_one_match, build_query_word_hashes, test_near_identical,
    BlastCompoAlignment, BlastCompoAlignmentContext, BlastCompoMatchingSequence,
    BlastCompoQueryInfo, BlastCompoSequenceData, BlastCompoSequenceRange, BlastRedoAlignCallbacks,
    BlastRedoRangeResult, EMatrixAdjustRule,
};
use crate::stats::{blast_spouge_stoe, BlastGumbelBlk, KarlinParams};
use crate::utils::matrix::{ncbistdaa, protein_score};

use super::alignment::build_alignment_view_with_matrix;
use super::gapalign::{blast_gapped_alignment_with_traceback, BlastpPreliminaryHsp};
use super::hsp::{
    reap_hsplist_by_evalue, sort_hsplist_by_score, update_best_evalue, BlastpHsp, BlastpHspList,
};

pub(crate) use crate::core::composition_adjustment::redo_alignment::BlastRedoAlignParams;

#[derive(Debug)]
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3675-3716
// ```c
// s_HitlistEvaluateAndPurge(&best_score, &best_evalue, ...);
// if (best_evalue <= hitParams->options->expect_value) {
//     ...
//     BlastCompo_HeapInsert(..., best_evalue, best_score, ...);
// }
// ```
pub(crate) struct BlastpPostprocessResult {
    pub hits: Vec<Hit>,
    pub best_score: i32,
    pub best_evalue: f64,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1000-1004
// ```c
// /** NCBIstdaa encoding for 'X' character */
// #define BLASTP_MASK_RESIDUE 21
// /** Default instructions and mask residue for SEG filtering */
// #define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"
// ```
fn blastp_subject_seg_params() -> SegParams {
    SegParams::new(10, 1.8, 2.1)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:770-803
// ```c
// new_align =
//     BlastCompo_AlignmentNew((int) (hsp->score * localScalingFactor),
//                             eDontAdjustMatrix,
//                             hsp->query.offset, hsp->query.end, hsp->context,
//                             hsp->subject.offset, hsp->subject.end,
//                             hsp->subject.frame, hsp);
// ```
fn build_incoming_alignment_list(
    preliminary_hits: &[BlastpPreliminaryHsp],
    score_divisor: f64,
) -> Option<Box<BlastCompoAlignment>> {
    let mut head: Option<Box<BlastCompoAlignment>> = None;
    let mut tail = &mut head;
    for (prelim_index, preliminary_hsp) in preliminary_hits.iter().enumerate() {
        let node = blast_compo_alignment_new(
            ((preliminary_hsp.raw_score as f64) * score_divisor) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            preliminary_hsp.query_start,
            preliminary_hsp.query_end,
            preliminary_hsp.query_context,
            preliminary_hsp.subject_start,
            preliminary_hsp.subject_end,
            preliminary_hsp.subject_frame,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:770-803
            // ```c
            // new_align =
            //     BlastCompo_AlignmentNew(..., hsp->subject.frame, hsp);
            // ```
            //
            // Rust stores the stable preliminary-hit index instead of a raw
            // `BlastHSP*`, preserving the same lookup without borrowing the
            // caller-owned slice across the redo-alignment list.
            Some(BlastCompoAlignmentContext::PreliminaryHspIndex(
                prelim_index,
            )),
        );
        *tail = Some(node);
        tail = &mut tail.as_mut().expect("incoming align inserted").next;
    }
    head
}

#[derive(Debug, Clone)]
struct RedoneBlastpHit {
    hit: Hit,
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1748-1770
// ```c
// s_NewAlignmentFromGapAlign(..., BlastCompo_SequenceRange * query_range,
//                            BlastCompo_SequenceRange * subject_range, ...)
// {
//     queryStart = gap_align->query_start + query_range->begin;
//     queryEnd   = gap_align->query_stop + query_range->begin;
//     matchStart = gap_align->subject_start + subject_range->begin;
//     matchEnd   = gap_align->subject_stop  + subject_range->begin;
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1902-1945
// ```c
// q_start = hsp->query.gapped_start - query_range->begin;
// s_start = hsp->subject.gapped_start - subject_range->begin;
// status = BLAST_GappedAlignmentWithTraceback(..., q_start, s_start, ...);
// return s_NewAlignmentFromGapAlign(..., query_range, subject_range, ...);
// ```
fn redo_preliminary_blastp_hit(
    preliminary_hsp: &BlastpPreliminaryHsp,
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    matrix: ScoringMatrix,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    _params: &KarlinParams,
    _effective_space: f64,
) -> Option<RedoneBlastpHit> {
    let q_start = preliminary_hsp.gapped_query_start - query_range.begin;
    let s_start = preliminary_hsp.gapped_subject_start - subject_range.begin;
    if q_start < 0 || s_start < 0 {
        return None;
    }

    let aligned = blast_gapped_alignment_with_traceback(
        query_data.data(),
        subject_data.data(),
        q_start as usize,
        s_start as usize,
        matrix,
        adjusted_matrix,
        gap_open,
        gap_extend,
        x_drop,
    )?;

    let alignment_len = aligned
        .edit_script
        .iter()
        .map(|op| op.num() as usize)
        .sum::<usize>();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1757-1770
    // ```c
    // queryStart = gap_align->query_start + query_range->begin;
    // queryEnd   = gap_align->query_stop + query_range->begin;
    // matchStart = gap_align->subject_start + subject_range->begin;
    // matchEnd   = gap_align->subject_stop  + subject_range->begin;
    // obj = BlastCompo_AlignmentNew(..., queryStart, queryEnd, ...,
    //                               matchStart, matchEnd, ...);
    // ```
    if alignment_len == 0 {
        return None;
    }

    let query_origin = usize::try_from(query_range.begin).ok()?;
    let subject_origin = usize::try_from(subject_range.begin).ok()?;
    let query_start = query_origin.checked_add(aligned.query_start)?;
    let query_end = query_origin.checked_add(aligned.query_stop)?;
    let subject_start = subject_origin.checked_add(aligned.subject_start)?;
    let subject_end = subject_origin.checked_add(aligned.subject_stop)?;

    Some(RedoneBlastpHit {
        hit: Hit {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: query_start + 1,
            q_end: query_end,
            s_start: subject_start + 1,
            s_end: subject_end,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: preliminary_hsp.query_frame,
            query_length: preliminary_hsp.query_length,
            q_idx: preliminary_hsp.q_idx,
            s_idx: preliminary_hsp.s_idx,
            raw_score: aligned.score,
            gap_info: Some(aligned.edit_script),
            num_positives: 0,
        },
        query_start: query_start as i32,
        query_end: query_end as i32,
        subject_start: subject_start as i32,
        subject_end: subject_end as i32,
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:305-354
// ```c
// s_HSPListFromDistinctAlignments(BlastHSPList *hsp_list,
//                                 BlastCompo_Alignment ** alignments,
//                                 int oid,
//                                 const BlastQueryInfo* queryInfo,
//                                 int frame)
// {
//     ...
//     editScript = align->context;
//     ...
//     Blast_HSPInit(..., align->score, &editScript, &new_hsp);
// }
// ```
fn redone_hit_from_alignment(
    align: &BlastCompoAlignment,
    template_hit: &BlastpPreliminaryHsp,
) -> Option<RedoneBlastpHit> {
    let query_start = usize::try_from(align.query_start).ok()?;
    let query_end = usize::try_from(align.query_end).ok()?;
    let subject_start = usize::try_from(align.match_start).ok()?;
    let subject_end = usize::try_from(align.match_end).ok()?;
    let gap_info = match align.context.as_ref() {
        Some(BlastCompoAlignmentContext::EditScript(edit_script)) => Some(edit_script.clone()),
        None => None,
        Some(BlastCompoAlignmentContext::PreliminaryHspIndex(_)) => return None,
    };

    Some(RedoneBlastpHit {
        hit: Hit {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: query_start + 1,
            q_end: query_end,
            s_start: subject_start + 1,
            s_end: subject_end,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: template_hit.query_frame,
            query_length: template_hit.query_length,
            q_idx: template_hit.q_idx,
            s_idx: template_hit.s_idx,
            raw_score: align.score,
            gap_info,
            num_positives: 0,
        },
        query_start: align.query_start,
        query_end: align.query_end,
        subject_start: align.match_start,
        subject_end: align.match_end,
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:305-354
// ```c
// s_HSPListFromDistinctAlignments(BlastHSPList *hsp_list,
//                                 BlastCompo_Alignment ** alignments,
//                                 int oid,
//                                 const BlastQueryInfo* queryInfo,
//                                 int frame)
// {
//     ...
//     editScript = align->context;
//     ...
//     Blast_HSPInit(..., align->score, &editScript, &new_hsp);
// }
// ```
fn hits_from_distinct_alignments(
    alignments: &Option<Box<BlastCompoAlignment>>,
    template_hit: &BlastpPreliminaryHsp,
) -> Vec<RedoneBlastpHit> {
    let mut hits = Vec::new();
    let mut current = alignments.as_deref();
    while let Some(align) = current {
        if let Some(hit) = redone_hit_from_alignment(align, template_hit) {
            hits.push(hit);
        }
        current = align.next.as_deref();
    }
    hits
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:305-354
// ```c
// s_HSPListFromDistinctAlignments(BlastHSPList *hsp_list,
//                                 BlastCompo_Alignment ** alignments,
//                                 int oid,
//                                 const BlastQueryInfo* queryInfo,
//                                 int frame)
// {
//     ...
//     Blast_HSPListSortByScore(hsp_list);
// }
// ```
fn hsplist_from_redone_hits(
    redone_hits: &[RedoneBlastpHit],
    oid: u32,
    query_index: u32,
) -> BlastpHspList {
    let mut hsp_list = BlastpHspList {
        oid,
        query_index,
        hsps: redone_hits
            .iter()
            .map(|redone_hit| BlastpHsp::from_hit(redone_hit.hit.clone()))
            .collect(),
        best_evalue: i32::MAX as f64,
    };
    sort_hsplist_by_score(&mut hsp_list);
    hsp_list
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1510-1647
// ```c
// static int
// s_SequenceGetProteinRange(..., const Boolean shouldTestIdentical, ...)
// {
//     ...
//     if (compo_adjust_mode && ...) {
//         if ((!shouldTestIdentical) || ... ) {
//             status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
//                                          subject_maybe_biased);
//         }
//     }
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1682-1704
// ```c
// Uint1 *origData = query->data + q_range->begin;
// queryData->length = q_range->end - q_range->begin;
// ...
// queryData->data[idx] = (origData[idx] != 24) ? origData[idx] : 3;
// ```
fn build_query_range_data(
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
) -> BlastCompoSequenceData {
    let begin = query_range.begin.max(0) as usize;
    let end = query_range.end.max(query_range.begin) as usize;
    let mut query = orig_query.data()[begin..end].to_vec();
    for residue in &mut query {
        if *residue == ncbistdaa::U {
            *residue = ncbistdaa::C;
        }
    }
    BlastCompoSequenceData::from_ncbistdaa(&query)
}

fn build_subject_range_data(
    matching_seq: &BlastCompoMatchingSequence,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    query_words: Option<&[u64]>,
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    subject_maybe_biased_in: bool,
    params: &BlastRedoAlignParams,
) -> BlastRedoRangeResult {
    let subject_len = matching_seq.data.len();
    let begin = (subject_range.begin.max(0) as usize).min(subject_len);
    let end = (subject_range.end.max(subject_range.begin) as usize).min(subject_len);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1677-1709
    // ```c
    // origData = query->data + q_range->begin;
    // queryData->length = q_range->end - q_range->begin;
    // ...
    // queryData->data[idx] = (origData[idx] != 24) ? origData[idx] : 3;
    // return s_SequenceGetProteinRange(..., queryData, ...);
    // ```
    let query = build_query_range_data(orig_query, query_range);
    let mut subject = BlastCompoSequenceData::from_ncbistdaa(&matching_seq.data);
    let mut subject_maybe_biased = subject_maybe_biased_in;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1626-1636
    // ```c
    // if (compo_adjust_mode
    //     && (!subject_maybe_biased || *subject_maybe_biased)) {
    //     if ((!shouldTestIdentical)
    //          || (shouldTestIdentical
    //              && (!s_TestNearIdentical(seqData, 0, queryData,
    //                                       q_range->begin, query_words,
    //                                       align)))) {
    //         status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
    //                                      subject_maybe_biased);
    //     }
    // }
    // ```
    let should_apply_seg = params.uses_composition_based_stats()
        && subject_maybe_biased
        && (!should_test_identical
            || !test_near_identical(&subject, 0, &query, query_range.begin, query_words, align));

    if should_apply_seg {
        let seg_params = blastp_subject_seg_params();
        let masker = SegMasker::with_params(&seg_params);
        let intervals = masker.mask_sequence(subject.data());
        subject_maybe_biased = !intervals.is_empty();
        for interval in &intervals {
            for pos in interval.start..interval.end {
                subject.buffer[subject.data_offset + pos] = ncbistdaa::X;
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1638-1646
    // ```c
    // seqData->data    = &seqData->data[range->begin - 1];
    // *seqData->data++ = '\0';
    // seqData->length  = range->end - range->begin;
    // ```
    subject = BlastCompoSequenceData::from_ncbistdaa(&subject.data()[begin..end]);

    BlastRedoRangeResult {
        query,
        subject,
        subject_maybe_biased,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:224-260
// ```c
// static void
// s_HitlistReapContained(BlastHSP * hsp_array[], Int4 * hspcnt)
// {
//     ...
// }
// ```
fn reap_contained_hits(hsp_list: &mut BlastpHspList) {
    if hsp_list.hsps.len() < 2 {
        return;
    }

    let mut keep = vec![true; hsp_list.hsps.len()];
    for iread in 1..hsp_list.hsps.len() {
        if !keep[iread] {
            continue;
        }
        let hsp1 = &hsp_list.hsps[iread];
        let hsp1_query_start = hsp1.q_start.saturating_sub(1) as i32;
        let hsp1_query_end = hsp1.q_end as i32;
        let hsp1_subject_start = hsp1.s_start.saturating_sub(1) as i32;
        let hsp1_subject_end = hsp1.s_end as i32;
        for iread_back in 0..iread {
            if !keep[iread_back] {
                continue;
            }
            let hsp2 = &hsp_list.hsps[iread_back];
            let hsp2_query_start = hsp2.q_start.saturating_sub(1) as i32;
            let hsp2_query_end = hsp2.q_end as i32;
            let hsp2_subject_start = hsp2.s_start.saturating_sub(1) as i32;
            let hsp2_subject_end = hsp2.s_end as i32;
            if hsp_list.hsps[iread_back].query_frame != hsp_list.hsps[iread].query_frame {
                continue;
            }
            let start_contained = hsp2_query_start <= hsp1_query_start
                && hsp2_query_end >= hsp1_query_start
                && hsp2_subject_start <= hsp1_subject_start
                && hsp2_subject_end >= hsp1_subject_start;
            let end_contained = hsp2_query_start <= hsp1_query_end
                && hsp2_query_end >= hsp1_query_end
                && hsp2_subject_start <= hsp1_subject_end
                && hsp2_subject_end >= hsp1_subject_end;
            if start_contained
                && end_contained
                && hsp_list.hsps[iread].raw_score <= hsp_list.hsps[iread_back].raw_score
            {
                keep[iread] = false;
                break;
            }
        }
    }

    let mut hsps = Vec::with_capacity(hsp_list.hsps.len());
    for (read_index, keep_hit) in keep.into_iter().enumerate() {
        if keep_hit {
            hsps.push(hsp_list.hsps[read_index].clone());
        }
    }
    hsp_list.hsps = hsps;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-116
// ```c
// static void
// s_HSPListNormalizeScores(BlastHSPList * hsp_list,
//                          double lambda,
//                          double logK,
//                          double scoreDivisor)
// {
//     ...
// }
// ```
fn normalize_scores(hsp_list: &mut BlastpHspList, lambda: f64, log_k: f64, score_divisor: f64) {
    for hsp in &mut hsp_list.hsps {
        hsp.raw_score = blast_nint((hsp.raw_score as f64) / score_divisor);
        hsp.bit_score =
            (hsp.raw_score as f64 * lambda * score_divisor - log_k) / std::f64::consts::LN_2;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441
// ```c
// long BLAST_Nint(double x)
// {
//    x += (x >= 0. ? 0.5 : -0.5);
//    return (long)x;
// }
// ```
fn blast_nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2120-2133
// ```c
// for (i = 0;  i < num_queries;  i++) {
//     if (sbp->kbp_gap[i] != NULL) {
//         Blast_KarlinBlk * kbp = sbp->kbp_gap[i];
//         kbp->Lambda /= scale_factor;
//         kbp->logK = log(kbp->K);
//     }
// }
// ```
fn scaled_karlin_params(params: &KarlinParams, score_divisor: f64) -> KarlinParams {
    if score_divisor == 1.0 {
        return params.clone();
    }
    KarlinParams {
        lambda: params.lambda / score_divisor,
        k: params.k,
        h: params.h,
        alpha: params.alpha,
        beta: params.beta,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:147-158
// ```c
// int query_length =      query_context->query_length;
// int length_adjustment = query_context->length_adjustment;
// double query_eff   = MAX((query_length - length_adjustment), 1);
// double dblen_eff = (double) query_context->eff_searchsp / query_eff;
// ```
fn build_query_info(
    query: &[u8],
    query_length: i32,
    length_adjustment: i32,
    effective_space: f64,
) -> BlastCompoQueryInfo {
    BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(query),
        composition: read_aa_composition(query),
        query_length,
        length_adjustment,
        eff_search_space: effective_space,
        words: None,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
// ```c
// static Int2
// s_Blast_HSPGetNumIdentitiesAndPositives(const Uint1* query, const Uint1* subject,
//                                         const BlastHSP* hsp, Int4* num_ident_ptr,
//                                         Int4* align_length_ptr, const BlastScoreBlk* sbp,
//                                         Int4* num_pos_ptr)
// {
//     ...
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:518-529
// ```c
// query = query_blk->sequence + query_info->contexts[hsp->context].query_offset;
// query_nomask = query_blk->sequence_nomask +
//     query_info->contexts[hsp->context].query_offset;
// ...
// status = Blast_HSPGetNumIdentitiesAndPositives(query_nomask, subject, hsp,
//                                                scoring_options, 0, sbp);
// ```
fn recompute_num_identities(
    query_nomask: &[u8],
    subject: &[u8],
    hits: &mut [Hit],
    matrix: ScoringMatrix,
) {
    for hit in hits {
        if let Some(view) = build_alignment_view_with_matrix(query_nomask, subject, hit, matrix) {
            hit.length = view.alignment_len;
            hit.num_ident = view.num_ident;
            hit.num_positives = view.num_positives;
            hit.gapopen = view.gap_opens;
            hit.mismatch = view
                .alignment_len
                .saturating_sub(view.num_ident + view.gap_letters);
            hit.identity = if view.alignment_len > 0 {
                (view.num_ident as f64 * 100.0) / (view.alignment_len as f64)
            } else {
                0.0
            };
            continue;
        }

        if hit.gap_info.is_some() || hit.q_start == 0 || hit.s_start == 0 {
            continue;
        }

        let q_start = hit.q_start - 1;
        let q_end = hit.q_end;
        let s_start = hit.s_start.min(hit.s_end) - 1;
        let s_end = hit.s_start.max(hit.s_end);
        if q_end < q_start || s_end < s_start {
            continue;
        }
        let q_span = q_end.saturating_sub(q_start);
        let s_span = s_end.saturating_sub(s_start);
        if q_span != s_span || q_end > query_nomask.len() || s_end > subject.len() {
            continue;
        }

        let mut num_ident = 0usize;
        let mut num_positives = 0usize;
        for (&q, &s) in query_nomask[q_start..q_end]
            .iter()
            .zip(subject[s_start..s_end].iter())
        {
            if q == s {
                num_ident += 1;
                num_positives += 1;
            } else if protein_score(matrix, q, s) > 0 {
                num_positives += 1;
            }
        }
        let align_len = q_span;
        hit.length = align_len;
        hit.num_ident = num_ident;
        hit.num_positives = num_positives;
        hit.mismatch = align_len.saturating_sub(num_ident);
        hit.gapopen = 0;
        hit.identity = if align_len > 0 {
            (num_ident as f64 * 100.0) / (align_len as f64)
        } else {
            0.0
        };
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:135-193
// ```c
// static void
// s_AdjustEvaluesForComposition(BlastHSPList *hsp_list, double comp_p_value,
//                               ..., Int4 subject_length,
//                               const BlastContextInfo * query_context,
//                               ...)
// {
//     double query_eff   = MAX((query_length - length_adjustment), 1);
//     double subject_eff = MAX((subject_length - length_adjustment), 1.0);
//     double dblen_eff = (double) query_context->eff_searchsp / query_eff;
//     double db_to_sequence_scale = subject_eff / dblen_eff;
//     ...
// }
// ```
fn adjust_evalues_for_composition(
    hsp_list: &mut BlastpHspList,
    comp_p_value: f64,
    query_info: &BlastCompoQueryInfo,
    subject_length: i32,
    lambda_ratio: f64,
) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:184-191
    // ```c
    // /* suppress unused parameter warnings if diagnostics are not printed */
    // (void) seqSrc;
    // (void) query_length;
    // (void) LambdaRatio;
    // (void) subject_id;
    // ```
    let _ = lambda_ratio;
    let query_eff = (query_info.query_length - query_info.length_adjustment).max(1) as f64;
    let subject_eff = (subject_length - query_info.length_adjustment).max(1) as f64;
    let dblen_eff = query_info.eff_search_space / query_eff;
    if !(dblen_eff.is_finite() && dblen_eff > 0.0) {
        return;
    }
    let db_to_sequence_scale = subject_eff / dblen_eff;
    if !(db_to_sequence_scale.is_finite() && db_to_sequence_scale > 0.0) {
        return;
    }

    let mut best_evalue = f64::MAX;
    for hsp in &mut hsp_list.hsps {
        hsp.e_value *= db_to_sequence_scale;
        let align_p_value = e_to_p(hsp.e_value);
        let combined_p_value = blast_overall_p_value(comp_p_value, align_p_value);
        hsp.e_value = p_to_e(combined_p_value);
        hsp.e_value /= db_to_sequence_scale;
        if hsp.e_value < best_evalue {
            best_evalue = hsp.e_value;
        }
    }
    hsp_list.best_evalue = best_evalue;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:395-457
// ```c
// static int
// s_HitlistEvaluateAndPurge(...)
// {
//     ...
//     Blast_HSPListReapByEvalue(hsp_list, hitParams->options);
// }
// ```
fn evaluate_and_purge(
    hsp_list: &mut BlastpHspList,
    params: &KarlinParams,
    gumbel_params: &BlastGumbelBlk,
    query_info: &BlastCompoQueryInfo,
    subject_length: i32,
    expect_value: f64,
    pvalue_for_this_pair: Option<f64>,
    lambda_ratio: Option<f64>,
    score_divisor: f64,
) -> (i32, f64) {
    let mut best_evalue = f64::MAX;
    let mut best_score = 0;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:414-427
    // ```c
    // status =
    //     Blast_HSPListGetEvalues(program_number, queryInfo,
    //                             s_GetSubjectLength(subject_length, program_number),
    //                             hsp_list, TRUE, FALSE, sbp, 0.0, 1.0);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:429-435
    // ```c
    // if ((0 <= pvalueForThisPair) && (pvalueForThisPair <= 1)) {
    //     s_AdjustEvaluesForComposition(..., LambdaRatio, ...);
    // }
    // ```
    let eval_params = scaled_karlin_params(params, score_divisor);
    for hsp in &mut hsp_list.hsps {
        hsp.e_value = blast_spouge_stoe(
            hsp.raw_score,
            &eval_params,
            gumbel_params,
            query_info.query_length,
            subject_length,
        );
    }
    update_best_evalue(hsp_list);
    if let Some(comp_p_value) = pvalue_for_this_pair.filter(|p| (0.0..=1.0).contains(p)) {
        adjust_evalues_for_composition(
            hsp_list,
            comp_p_value,
            query_info,
            subject_length,
            lambda_ratio.unwrap_or(1.0),
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:439-446
    // ```c
    // Blast_HSPListReapByEvalue(hsp_list, hitParams->options);
    // if (hsp_list->hspcnt > 0) {
    //     *pbestEvalue = hsp_list->best_evalue;
    //     *pbestScore  = hsp_list->hsp_array[0]->score;
    // }
    // ```
    reap_hsplist_by_evalue(hsp_list, expect_value);
    if let Some(first_hsp) = hsp_list.hsps.first() {
        best_evalue = hsp_list.best_evalue;
        best_score = first_hsp.raw_score;
    }

    (best_score, best_evalue)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3658-3710
// ```c
// s_HSPListFromDistinctAlignments(...);
// if (hsp_list->hspcnt > 1) {
//     s_HitlistReapContained(...);
// }
// s_HitlistEvaluateAndPurge(...);
// if (best_evalue <= hitParams->options->expect_value) {
//     s_HSPListNormalizeScores(...);
//     s_ComputeNumIdentities(...);
// }
// ```
pub(crate) fn postprocess_preliminary_hits(
    preliminary_hits: Vec<BlastpPreliminaryHsp>,
    query: &[u8],
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    karlin_params: &KarlinParams,
    gumbel_params: &BlastGumbelBlk,
    query_length: i32,
    length_adjustment: i32,
    effective_space: f64,
    redo_align_params: &BlastRedoAlignParams,
    expect_value: f64,
) -> Result<BlastpPostprocessResult> {
    if preliminary_hits.is_empty() {
        return Ok(BlastpPostprocessResult {
            hits: Vec::new(),
            best_score: 0,
            best_evalue: f64::MAX,
        });
    }

    let incoming_aligns =
        build_incoming_alignment_list(&preliminary_hits, redo_align_params.score_divisor);
    let mut query_info = build_query_info(query, query_length, length_adjustment, effective_space);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2321-2331
    // ```c
    // query_info->words = NULL;
    // s_CreateWordArray(query_info->seq.data, query_info->seq.length,
    //                   &query_info->words);
    // ```
    let query_words = build_query_word_hashes(query);
    query_info.words = (!query_words.is_empty()).then_some(query_words);
    let matching_seq = BlastCompoMatchingSequence::new(preliminary_hits[0].s_idx as i32, subject);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2117-2133
    // ```c
    // if (sbp->kbp_gap[i] != NULL) {
    //     Blast_KarlinBlk * kbp = sbp->kbp_gap[i];
    //     kbp->Lambda /= scale_factor;
    // }
    // ```
    let redo_karlin_params = scaled_karlin_params(karlin_params, redo_align_params.score_divisor);

    let mut callbacks = BlastRedoAlignCallbacks {
        calc_lambda: None,
        get_range: Box::new(
            |matching_seq,
             subject_range,
             orig_query,
             query_range,
             query_words,
             align,
             should_test_identical,
             subject_maybe_biased,
             redo_params| {
                let range_result = build_subject_range_data(
                    matching_seq,
                    subject_range,
                    orig_query,
                    query_range,
                    query_words,
                    align,
                    should_test_identical,
                    subject_maybe_biased,
                    redo_params,
                );
                Ok(range_result)
            },
        ),
        redo_one_alignment: Box::new(
            |incoming_align,
             matrix_adjust_rule,
             adjusted_matrix,
             query_data,
             query_range,
             _ccat_query_length,
             subject_data,
             subject_range,
             _full_subject_length,
             redo_params| {
                let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(prelim_index)) =
                    incoming_align.context.as_ref()
                else {
                    return Ok(None);
                };
                let preliminary_hsp = &preliminary_hits[*prelim_index];
                let Some(redone_hit) = redo_preliminary_blastp_hit(
                    preliminary_hsp,
                    query_data,
                    query_range,
                    subject_data,
                    subject_range,
                    matrix,
                    adjusted_matrix,
                    redo_params.gapping_params.gap_open,
                    redo_params.gapping_params.gap_extend,
                    redo_params.gapping_params.x_dropoff,
                    &redo_karlin_params,
                    effective_space,
                ) else {
                    return Ok(None);
                };
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1902-1945
                // ```c
                // status = BLAST_GappedAlignmentWithTraceback(...);
                // return s_NewAlignmentFromGapAlign(gap_align, ...);
                // ```
                Ok(Some(blast_compo_alignment_new(
                    redone_hit.hit.raw_score,
                    matrix_adjust_rule,
                    redone_hit.query_start,
                    redone_hit.query_end,
                    query_range.context,
                    redone_hit.subject_start,
                    redone_hit.subject_end,
                    subject_range.context,
                    redone_hit
                        .hit
                        .gap_info
                        .clone()
                        .map(BlastCompoAlignmentContext::EditScript),
                )))
            },
        ),
        new_xdrop_align: None,
        free_align_traceback: None,
    };

    let distinct_alignments = blast_redo_one_match(
        &incoming_aligns,
        redo_align_params,
        &matching_seq,
        &query_info,
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2120-2128
        // ```c
        // if (sbp->kbp_gap[i] != NULL) {
        //     Blast_KarlinBlk * kbp = sbp->kbp_gap[i];
        //     kbp->Lambda /= scale_factor;
        // }
        // ```
        karlin_params.lambda / redo_align_params.score_divisor,
        &mut callbacks,
    )?;
    drop(callbacks);
    let query_index =
        usize::try_from(preliminary_hits[0].query_context).expect("blastp query context is non-negative");
    let kept_hits = distinct_alignments
        .alignments_by_query
        .get(query_index)
        .map(|alignments| hits_from_distinct_alignments(alignments, &preliminary_hits[0]))
        .unwrap_or_default();
    let mut hsp_list =
        hsplist_from_redone_hits(&kept_hits, matching_seq.index as u32, preliminary_hits[0].q_idx);

    if hsp_list.hsps.len() > 1 {
        reap_contained_hits(&mut hsp_list);
    }
    let pvalue_for_this_pair = distinct_alignments.pvalue_for_this_pair;
    let lambda_ratio = distinct_alignments.lambda_ratio;
    let (best_score, best_evalue) = evaluate_and_purge(
        &mut hsp_list,
        karlin_params,
        gumbel_params,
        &query_info,
        matching_seq.length,
        expect_value,
        pvalue_for_this_pair,
        lambda_ratio,
        redo_align_params.score_divisor,
    );
    if best_evalue <= expect_value && !hsp_list.hsps.is_empty() {
        normalize_scores(
            &mut hsp_list,
            redo_karlin_params.lambda,
            redo_align_params.log_k,
            redo_align_params.score_divisor,
        );
        let mut raw_hits: Vec<Hit> = hsp_list.hsps.drain(..).map(BlastpHsp::into_hit).collect();
        recompute_num_identities(query_nomask, subject, &mut raw_hits, matrix);
        hsp_list.hsps = raw_hits.into_iter().map(BlastpHsp::from_hit).collect();
    }
    Ok(BlastpPostprocessResult {
        hits: hsp_list.hsps.into_iter().map(BlastpHsp::into_hit).collect(),
        best_score,
        best_evalue,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_redo_align_params() -> BlastRedoAlignParams {
        BlastRedoAlignParams {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2216-2231
            // ```c
            // self->ungappedLambda = sbp->kbp_ideal->Lambda / scale_factor;
            // status = s_GetStartFreqRatios(self->startFreqRatios, matrixName);
            // ```
            matrix_info:
                crate::core::composition_adjustment::adjust_scores::build_matrix_info(
                    ScoringMatrix::Blosum62,
                    0.3176,
                )
                .unwrap(),
            gapping_params:
                crate::core::composition_adjustment::redo_alignment::BlastCompoGappingParams {
                    gap_open: 11,
                    gap_extend: 1,
                    decline_align: 0,
                    x_dropoff: 50,
                    context: None,
                },
            compo_adjust_mode:
                crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode::CompositionMatrixAdjust,
            alphsize: crate::utils::matrix::BLASTAA_SIZE as i32,
            composition_test_index: 0,
            unified_p: false,
            log_k: 0.0,
            score_divisor: 1.0,
            restricted_alignment: false,
            smith_waterman: false,
            is_same_adjustment: false,
            near_identical_cutoff: 0.0,
            position_based: false,
            re_matrix_adjustment_pseudocounts: 20,
            ccat_query_length: 80,
            query_is_translated: false,
            subject_is_translated: false,
            cutoff_score: 0,
            cutoff_evalue: 10.0,
            do_link_hsps: false,
        }
    }

    fn make_hit(score: i32, q_start: usize, q_end: usize, s_start: usize, s_end: usize) -> Hit {
        Hit {
            identity: 0.0,
            length: q_end.saturating_sub(q_start).saturating_add(1),
            mismatch: 0,
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: 0,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
            gap_info: None,
            num_positives: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:224-260
    // ```c
    // if (CONTAINED_IN_HSP
    //     (hsp2->query.offset, hsp2->query.end, hsp1->query.offset,
    //      hsp2->subject.offset, hsp2->subject.end, hsp1->subject.offset) &&
    //     CONTAINED_IN_HSP
    //     (hsp2->query.offset, hsp2->query.end, hsp1->query.end,
    //      hsp2->subject.offset, hsp2->subject.end, hsp1->subject.end) &&
    //     hsp1->score <= hsp2->score) {
    //     hsp1 = hsp_array[iread] = Blast_HSPFree(hsp_array[iread]);
    // }
    // ```
    #[test]
    fn test_reap_contained_hits_uses_internal_offsets_ncbi_style() {
        let outer = RedoneBlastpHit {
            hit: make_hit(200, 101, 150, 201, 250),
            query_start: 100,
            query_end: 150,
            subject_start: 200,
            subject_end: 250,
        };
        let inner = RedoneBlastpHit {
            hit: make_hit(180, 101, 149, 201, 249),
            query_start: 100,
            query_end: 149,
            subject_start: 200,
            subject_end: 249,
        };

        let kept_hits = vec![outer.clone(), inner];
        let mut hsp_list = hsplist_from_redone_hits(&kept_hits, 0, 0);
        reap_contained_hits(&mut hsp_list);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert_eq!(hsp_list.hsps[0].raw_score, outer.hit.raw_score);
        assert_eq!(hsp_list.hsps[0].q_start, outer.hit.q_start);
        assert_eq!(hsp_list.hsps[0].s_start, outer.hit.s_start);
    }

    fn make_redone_hit(
        score: i32,
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
    ) -> RedoneBlastpHit {
        RedoneBlastpHit {
            hit: make_hit(score, q_start, q_end, s_start, s_end),
            query_start: q_start.saturating_sub(1) as i32,
            query_end: q_end as i32,
            subject_start: s_start.saturating_sub(1) as i32,
            subject_end: s_end as i32,
        }
    }

    #[test]
    fn test_hits_from_distinct_alignments_uses_edit_script_context() {
        let template_hit = BlastpPreliminaryHsp {
            query_start: 0,
            query_end: 20,
            subject_start: 4,
            subject_end: 24,
            gapped_query_start: 0,
            gapped_subject_start: 4,
            raw_score: 80,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: 9,
        };
        let mut alignments = Some(blast_compo_alignment_new(
            80,
            EMatrixAdjustRule::DontAdjustMatrix,
            0,
            20,
            0,
            4,
            24,
            0,
            Some(BlastCompoAlignmentContext::EditScript(vec![
                crate::common::GapEditOp::Sub(10),
                crate::common::GapEditOp::Del(1),
                crate::common::GapEditOp::Sub(9),
            ])),
        ));
        alignments.as_mut().unwrap().next = Some(blast_compo_alignment_new(
            70,
            EMatrixAdjustRule::DontAdjustMatrix,
            1,
            18,
            0,
            7,
            24,
            0,
            Some(BlastCompoAlignmentContext::EditScript(vec![
                crate::common::GapEditOp::Sub(17),
            ])),
        ));

        let hits = hits_from_distinct_alignments(&alignments, &template_hit);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].hit.raw_score, 80);
        assert_eq!(hits[1].hit.raw_score, 70);
        assert_eq!(
            hits[0].hit.gap_info,
            Some(vec![
                crate::common::GapEditOp::Sub(10),
                crate::common::GapEditOp::Del(1),
                crate::common::GapEditOp::Sub(9),
            ])
        );
        assert_eq!(hits[0].hit.q_start, 1);
        assert_eq!(hits[0].hit.s_start, 5);
        assert_eq!(hits[0].hit.s_idx, 9);
    }

    #[test]
    fn test_build_incoming_alignment_list_uses_truncating_scaled_score_and_query_index() {
        let preliminary_hits = vec![BlastpPreliminaryHsp {
            query_start: 10,
            query_end: 40,
            subject_start: 20,
            subject_end: 50,
            gapped_query_start: 25,
            gapped_subject_start: 35,
            raw_score: 7,
            query_context: 7,
            query_frame: 0,
            subject_frame: -1,
            query_length: 100,
            q_idx: 3,
            s_idx: 9,
        }];

        let alignments =
            build_incoming_alignment_list(&preliminary_hits, 3.2).expect("incoming alignment list");
        assert_eq!(alignments.score, 22);
        assert_eq!(alignments.query_index, 7);
        assert_eq!(alignments.query_start, 10);
        assert_eq!(alignments.query_end, 40);
        assert_eq!(alignments.match_start, 20);
        assert_eq!(alignments.match_end, 50);
        assert_eq!(alignments.frame, -1);
        assert_eq!(
            alignments.context,
            Some(BlastCompoAlignmentContext::PreliminaryHspIndex(0))
        );
    }

    #[test]
    fn test_adjust_evalues_for_composition_matches_ncbi_sequence_scaling() {
        let query_info = BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(&vec![1u8; 100]),
            composition: read_aa_composition(&vec![1u8; 100]),
            query_length: 100,
            length_adjustment: 10,
            eff_search_space: 9000.0,
            words: None,
        };
        let mut hsp_list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(80, 1, 20, 5, 24))],
            best_evalue: f64::MAX,
        };
        hsp_list.hsps[0].e_value = 2.5;

        adjust_evalues_for_composition(&mut hsp_list, 0.25, &query_info, 40, 1.0);

        let query_eff = 90.0;
        let subject_eff = 30.0;
        let db_to_sequence_scale = subject_eff / (query_info.eff_search_space / query_eff);
        let expected = p_to_e(blast_overall_p_value(
            0.25,
            e_to_p(2.5 * db_to_sequence_scale),
        )) / db_to_sequence_scale;
        assert!((hsp_list.hsps[0].e_value - expected).abs() < 1e-12);
        assert!((hsp_list.best_evalue - expected).abs() < 1e-12);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-116
    // ```c
    // hsp->score = (Int4)BLAST_Nint(((double) hsp->score) / scoreDivisor);
    // hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;
    // ```
    #[test]
    fn test_normalize_scores_uses_blast_nint_and_scaled_bit_score() {
        let mut hsp_list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(9, 1, 5, 1, 5))],
            best_evalue: f64::MAX,
        };

        normalize_scores(&mut hsp_list, 0.25, 0.0, 2.0);

        assert_eq!(hsp_list.hsps[0].raw_score, 5);
        let expected = (5.0_f64 * 0.25 * 2.0) / std::f64::consts::LN_2;
        assert!((hsp_list.hsps[0].bit_score - expected).abs() < 1e-12);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1811-1893
    // ```c
    // if (sbp->gbp) {
    //     hsp->evalue = BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
    //         query_info->contexts[hsp->context].query_length, subject_length);
    // } else {
    //     hsp->evalue = BLAST_KarlinStoE_simple(...);
    // }
    // ```
    #[test]
    fn test_evaluate_and_purge_uses_spouge_for_gapped_blastp() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let gumbel = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.669720,
            g: 12.0,
            a: 1.9,
            alpha: 42.602800,
            sigma: 43.636200,
            a_un: 0.7916,
            alpha_un: 4.964660,
            b: 2.0 * 12.0 * (0.7916 - 1.9),
            beta: 2.0 * 12.0 * (4.964660 - 42.602800),
            tau: 2.0 * 12.0 * (4.964660 - 43.636200),
            db_length: 5000,
        };
        let mut hsp_list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(1062, 1, 10, 1, 10))],
            best_evalue: f64::MAX,
        };
        let query_info = build_query_info(&vec![1u8; 100], 100, 0, 100000.0);

        let (_best_score, best_evalue) = evaluate_and_purge(
            &mut hsp_list,
            &params,
            &gumbel,
            &query_info,
            200,
            10.0,
            None,
            None,
            1.0,
        );

        let expected = blast_spouge_stoe(1062, &params, &gumbel, 100, 200);
        assert_eq!(hsp_list.hsps.len(), 1);
        assert!((hsp_list.hsps[0].e_value - expected).abs() < 1e-12);
        assert!((best_evalue - expected).abs() < 1e-12);
    }

    #[test]
    fn test_evaluate_and_purge_uses_scaled_lambda_for_spouge() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let gumbel = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.669720,
            g: 12.0,
            a: 1.9,
            alpha: 42.602800,
            sigma: 43.636200,
            a_un: 0.7916,
            alpha_un: 4.964660,
            b: 2.0 * 12.0 * (0.7916 - 1.9),
            beta: 2.0 * 12.0 * (4.964660 - 42.602800),
            tau: 2.0 * 12.0 * (4.964660 - 43.636200),
            db_length: 5000,
        };
        let mut hsp_list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(1062, 1, 10, 1, 10))],
            best_evalue: f64::MAX,
        };
        let query_info = build_query_info(&vec![1u8; 100], 100, 0, 100000.0);

        let (_best_score, best_evalue) = evaluate_and_purge(
            &mut hsp_list,
            &params,
            &gumbel,
            &query_info,
            200,
            10.0,
            None,
            None,
            32.0,
        );

        assert_eq!(hsp_list.hsps.len(), 1);
        let expected = blast_spouge_stoe(
            1062,
            &scaled_karlin_params(&params, 32.0),
            &gumbel,
            100,
            200,
        );
        assert!((hsp_list.hsps[0].e_value - expected).abs() < 1e-12);
        assert!((best_evalue - expected).abs() < 1e-12);
    }

    #[test]
    fn test_redo_preliminary_blastp_hit_preserves_internal_offsets_with_ranges() {
        let full_query = BlastCompoSequenceData::from_ncbistdaa(&vec![1u8; 12]);
        let full_subject = BlastCompoSequenceData::from_ncbistdaa(&vec![1u8; 11]);
        let query_range = BlastCompoSequenceRange {
            begin: 4,
            end: 12,
            context: 0,
        };
        let subject_range = BlastCompoSequenceRange {
            begin: 3,
            end: 11,
            context: 0,
        };
        let query_data = full_query.slice_range(&query_range);
        let subject_data = full_subject.slice_range(&subject_range);
        let preliminary_hsp = BlastpPreliminaryHsp {
            query_start: 4,
            query_end: 12,
            subject_start: 3,
            subject_end: 11,
            gapped_query_start: 4,
            gapped_subject_start: 3,
            raw_score: 0,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 12,
            q_idx: 0,
            s_idx: 0,
        };
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        let redone = redo_preliminary_blastp_hit(
            &preliminary_hsp,
            &query_data,
            &query_range,
            &subject_data,
            &subject_range,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            50,
            &params,
            100.0,
        )
        .expect("traceback should produce a hit");

        assert_eq!(redone.query_start, 4);
        assert_eq!(redone.subject_start, 3);
        assert_eq!(redone.hit.q_start, 5);
        assert_eq!(redone.hit.s_start, 4);
    }

    #[test]
    fn test_subject_seg_masks_when_enabled_and_not_near_identical() {
        let matching_seq = BlastCompoMatchingSequence::new(0, &vec![1u8; 80]);
        let subject_range = BlastCompoSequenceRange {
            begin: 0,
            end: 80,
            context: 0,
        };
        let params = test_redo_align_params();

        let query = BlastCompoSequenceData::from_ncbistdaa(&vec![2u8; 80]);
        let query_range = BlastCompoSequenceRange {
            begin: 0,
            end: 80,
            context: 0,
        };
        let align = blast_compo_alignment_new(
            80,
            EMatrixAdjustRule::DontAdjustMatrix,
            0,
            80,
            0,
            0,
            80,
            0,
            None,
        );
        let range = build_subject_range_data(
            &matching_seq,
            &subject_range,
            &query,
            &query_range,
            Some(&build_query_word_hashes(query.data())),
            &align,
            false,
            true,
            &params,
        );
        assert!(range.subject_maybe_biased);
        assert!(range.subject.data().contains(&ncbistdaa::X));
    }

    #[test]
    fn test_subject_seg_skips_exact_near_identical_subject() {
        let residues = vec![1u8; 80];
        let matching_seq = BlastCompoMatchingSequence::new(0, &residues);
        let subject_range = BlastCompoSequenceRange {
            begin: 0,
            end: 80,
            context: 0,
        };
        let query = BlastCompoSequenceData::from_ncbistdaa(&residues);
        let query_range = BlastCompoSequenceRange {
            begin: 0,
            end: 80,
            context: 0,
        };
        let params = test_redo_align_params();
        let align = blast_compo_alignment_new(
            80,
            EMatrixAdjustRule::DontAdjustMatrix,
            0,
            80,
            0,
            0,
            80,
            0,
            None,
        );

        let range = build_subject_range_data(
            &matching_seq,
            &subject_range,
            &query,
            &query_range,
            Some(&build_query_word_hashes(query.data())),
            &align,
            true,
            true,
            &params,
        );
        assert!(range.subject_maybe_biased);
        assert!(!range.subject.data().contains(&ncbistdaa::X));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1626-1636
    // ```c
    // if (compo_adjust_mode
    //     && (!subject_maybe_biased || *subject_maybe_biased)) {
    //     ...
    //     status = s_DoSegSequenceData(..., subject_maybe_biased);
    // }
    // ```
    #[test]
    fn test_subject_seg_preserves_prior_not_biased_state() {
        let residues = vec![1u8; 80];
        let matching_seq = BlastCompoMatchingSequence::new(0, &residues);
        let subject_range = BlastCompoSequenceRange {
            begin: 0,
            end: 80,
            context: 0,
        };
        let query = BlastCompoSequenceData::from_ncbistdaa(&residues);
        let query_range = BlastCompoSequenceRange {
            begin: 0,
            end: 80,
            context: 0,
        };
        let params = test_redo_align_params();
        let align = blast_compo_alignment_new(
            80,
            EMatrixAdjustRule::DontAdjustMatrix,
            0,
            80,
            0,
            0,
            80,
            0,
            None,
        );

        let range = build_subject_range_data(
            &matching_seq,
            &subject_range,
            &query,
            &query_range,
            Some(&build_query_word_hashes(query.data())),
            &align,
            false,
            false,
            &params,
        );
        assert!(!range.subject_maybe_biased);
        assert!(!range.subject.data().contains(&ncbistdaa::X));
    }

    #[test]
    fn test_subject_range_slice_has_fresh_leading_sentinel() {
        let matching_seq = BlastCompoMatchingSequence::new(0, &[1u8, 2, 3, 4, 5, 6]);
        let subject_range = BlastCompoSequenceRange {
            begin: 2,
            end: 5,
            context: 0,
        };
        let query = BlastCompoSequenceData::from_ncbistdaa(&[2u8, 3, 4, 5, 6, 7]);
        let query_range = BlastCompoSequenceRange {
            begin: 0,
            end: 6,
            context: 0,
        };
        let params = test_redo_align_params();
        let align = blast_compo_alignment_new(
            42,
            EMatrixAdjustRule::DontAdjustMatrix,
            0,
            6,
            0,
            0,
            6,
            0,
            None,
        );

        let range = build_subject_range_data(
            &matching_seq,
            &subject_range,
            &query,
            &query_range,
            None,
            &align,
            false,
            true,
            &params,
        );

        assert_eq!(range.subject.data(), &[3u8, 4, 5]);
        assert_eq!(range.subject.buffer[range.subject.data_offset - 1], 0);
    }

    #[test]
    fn test_reap_contained_hits_discards_lower_scoring_inner_hit() {
        let hits = vec![
            RedoneBlastpHit {
                hit: make_hit(100, 11, 50, 101, 140),
                query_start: 10,
                query_end: 50,
                subject_start: 100,
                subject_end: 140,
            },
            RedoneBlastpHit {
                hit: make_hit(90, 21, 40, 111, 130),
                query_start: 20,
                query_end: 40,
                subject_start: 110,
                subject_end: 130,
            },
            RedoneBlastpHit {
                hit: make_hit(80, 61, 80, 151, 170),
                query_start: 60,
                query_end: 80,
                subject_start: 150,
                subject_end: 170,
            },
        ];
        let mut hsp_list = hsplist_from_redone_hits(&hits, 0, 0);
        reap_contained_hits(&mut hsp_list);
        assert_eq!(hsp_list.hsps.len(), 2);
        assert_eq!(hsp_list.hsps[0].raw_score, 100);
        assert_eq!(hsp_list.hsps[1].raw_score, 80);
    }

    #[test]
    fn test_recompute_num_identities_ungapped() {
        let mut hits = vec![make_hit(20, 1, 4, 1, 4)];
        recompute_num_identities(
            &[1u8, 2, 3, 4],
            &[1u8, 2, 5, 4],
            &mut hits,
            ScoringMatrix::Blosum62,
        );
        assert_eq!(hits[0].num_ident, 3);
        assert_eq!(hits[0].length, 4);
        assert_eq!(hits[0].mismatch, 1);
    }

    #[test]
    fn test_recompute_num_identities_uses_query_nomask() {
        let mut hits = vec![make_hit(20, 1, 4, 1, 4)];
        recompute_num_identities(
            &[1u8, 2, 3, 4],
            &[1u8, 2, 3, 4],
            &mut hits,
            ScoringMatrix::Blosum62,
        );
        assert_eq!(hits[0].num_ident, 4);
        assert_eq!(hits[0].num_positives, 4);
    }
}

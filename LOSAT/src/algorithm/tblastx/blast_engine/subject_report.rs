//! The report hits of a TBLASTX subject (`linked_hits_to_report`): its linked HSPs, reaped by
//! E-value, with their bit scores and identities in report coordinates, as `TblastxHsp`s, with the
//! stage dumps and the diagnostics of subject_debug.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1535-1541 (in
//! `BLAST_PreliminarySearchEngine`: `s_Blast_HSPListReapByPrelimEvalue`, then
//! `Blast_HSPListGetBitScores`, c++/src/algo/blast/core/blast_hits.c:1907-1932), and
//! c++/src/algo/blast/api/blast_seqalign.cpp:1397 (`BLASTUngappedHspListToSeqAlign`, which turns
//! the HSPs into the report's alignments).
//!
//! This module is a child of `run_impl` for the timer it keeps private (`x_search_timing`). Its
//! function is made of statements of `process_subject` moved here in their order; it takes the
//! linked list by value and reads a struct of references and `Copy` values.

use super::subject_debug;
use super::subject_link::SubjectStageArgs;
use super::x_search_timing::XLap;
use super::{
    calc_bit_score, convert_coords, stage_dump, trace_final_hit_if_match, AtomicOrdering, Hit,
    QueryFrame, TblastxHsp, TbxScore, UngappedHit, XStage,
};

/// The report hits of the subject's linked HSPs (`linked`): the HSPs whose E-value passes
/// `evalue_threshold`, with their identities (counted by the re-evaluation), bit scores and
/// coordinates on the query and on the subject (`s_len` nucleotides; `s_frames` are the frames
/// the HSPs were scored on, `s_frames_report` those the report reads).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1535 and :1541 (the reap
/// and `Blast_HSPListGetBitScores`, c++/src/algo/blast/core/blast_hits.c:1907-1932).
pub(super) fn linked_hits_to_report<S: TbxScore>(
    sc: S,
    linked: Vec<UngappedHit>,
    (s_frames, s_frames_report): (&[QueryFrame], &[QueryFrame]),
    s_len: usize,
    stage: &SubjectStageArgs<'_>,
    x_lap: &mut XLap<'_>,
) -> Vec<TblastxHsp> {
    let &SubjectStageArgs {
        contexts_ref,
        evalue_threshold,
        diag_enabled,
        diagnostics,
        debug_output_filter,
        ..
    } = stage;
    // NCBI's culling (`-culling_limit`) runs on the final HSPs of all subjects, after
    // the e-value reap (`report::culled_hit_order`).

    // LOSAT_DEBUG_OUTPUT_FILTER: the roles of the linked HSPs (`subject_debug::count_chain_roles`).
    let chain_roles = subject_debug::count_chain_roles(&linked, debug_output_filter);

    let mut final_hits: Vec<TblastxHsp> = Vec::new();
    let dump_output_stage = stage_dump::enabled();
    let mut output_snapshot_hits: Vec<UngappedHit> = Vec::new();
    let mut output_snapshot_pairs: Vec<(Hit, UngappedHit)> = Vec::new();
    let mut filtered_by_evalue = 0usize;
    for h in linked {
        // NCBI reference (verbatim, link_hsps.c:1018-1020):
        //   /* If this is not a single piece or the start of a chain, then Skip it. */
        //   if (H->linked_set == TRUE && H->start_of_chain == FALSE)
        //       continue;
        // NOTE: This "continue" in NCBI is NOT a filter - NCBI then walks the link
        // pointer from each chain head to include all chain members (lines 1047-1059).
        // LOSAT's linking already produces a flat list of ALL HSPs (including chain
        // members) so we output everything without this skip logic.

        // NCBI reference: E-value filtering is applied during output conversion
        // The exact timing may differ, but the threshold check is standard
        if h.e_value > evalue_threshold {
            filtered_by_evalue += 1;
            if diag_enabled {
                diagnostics
                    .ungapped_evalue_failed
                    .fetch_add(1, AtomicOrdering::Relaxed);
            }
            continue;
        }

        if diag_enabled {
            diagnostics
                .ungapped_evalue_passed
                .fetch_add(1, AtomicOrdering::Relaxed);
        }
        if dump_output_stage {
            output_snapshot_hits.push(h.clone());
        }

        let ctx = &contexts_ref[h.ctx_idx as usize];
        let s_score_frame = &s_frames[usize::from(h.s_f_idx)];
        let s_frame = &s_frames_report[usize::from(h.s_f_idx)];

        let len = (h.q_aa_end as usize).saturating_sub(h.q_aa_start as usize);
        let q0 = h.q_aa_start as usize;
        let s0 = h.s_aa_start as usize;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2708-2720
        // ```c
        // Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
        //     subject_start, hsp, score_params->options, &align_length, sbp);
        // delete_hsp = Blast_HSPTest(hsp, ...);
        // ```
        // `reevaluate_ungapped_hsp_list` already performs this NCBI step
        // using the unmasked query and reporting subject buffers, then stores
        // the identity count on the HSP. Reuse it here instead of rescanning
        // every final hit.
        let matches = h.num_ident as usize;
        let mismatch = len.saturating_sub(matches);
        let identity = if len > 0 {
            (matches as f64 / len as f64) * 100.0
        } else {
            0.0
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1918-1928
        // ```c
        // kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
        // hsp->bit_score =
        //    (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
        //    NCBIMATH_LN2;
        // ```
        let bit_params = &contexts_ref[h.ctx_idx as usize].karlin_params;
        let bit = calc_bit_score(h.raw_score, bit_params);
        let (q_start, q_end) = convert_coords(
            h.q_aa_start as usize,
            h.q_aa_end as usize,
            ctx.frame,
            ctx.orig_len,
        );
        let (s_start, s_end) = convert_coords(
            h.s_aa_start as usize,
            h.s_aa_end as usize,
            s_frame.frame,
            s_len,
        );

        // Diagnostic (LOSAT_TRACE_HSP_MASKS, diagnostic build): `subject_debug::trace_hsp_masks`.
        // The variable is read here, so that a build without diagnostics drops the call.
        if crate::utils::xdiag::var_os("LOSAT_TRACE_HSP_MASKS").is_some() {
            subject_debug::trace_hsp_masks(
                sc,
                &subject_debug::MaskTrace {
                    h: &h,
                    ctx,
                    s_frame,
                    s_score_frame,
                    len,
                    q0,
                    s0,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                },
            );
        }

        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        let out_hit = Hit {
            identity,
            length: len,
            mismatch: mismatch as i64,
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: h.e_value,
            bit_score: bit,
            num_ident: matches,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame: ctx.frame as i32,
            query_length: 0,
            q_idx: ctx.q_idx,
            s_idx: h.s_idx,
            raw_score: h.raw_score,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4754-4767
            // ```c
            // s_AdjustInitialHSPOffsets(init_hsp,
            //                           query_info->contexts[context].query_offset);
            // Blast_HSPInit(ungapped_data->q_start, ...,
            //               ungapped_data->s_start, ...,
            //               context, query_info->contexts[context].frame,
            //               subject->frame, ungapped_data->score, NULL, &new_hsp);
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:808-812
            // ```c
            // subject->sequence = translation_buffer + frame_offsets[context] + 1;
            // subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
            // ```
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1347-1353
            // ```c
            // BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)
            // BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)
            // BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset)
            // BLAST_CMP(hsp2->query.end, hsp1->query.end)
            // ```
            // TBLASTX ScoreCompareHSPs uses the context/frame-relative
            // BlastSeg offsets stored in Blast_HSPInit, not formatted
            // nucleotide coordinates or LOSAT's concatenated frame bases.
            sort_query_offset: h.q_aa_start as usize,
            sort_query_end: h.q_aa_end as usize,
            sort_subject_offset: h.s_aa_start as usize,
            sort_subject_end: h.s_aa_end as usize,
            has_sort_offsets: true,
            gap_info: None,
            num_positives: matches,
        };
        trace_final_hit_if_match("output_hit", &out_hit);
        if dump_output_stage {
            output_snapshot_pairs.push((out_hit.clone(), h.clone()));
        }
        // NCBI reference: c++/src/algo/blast/core/blast_gapalign.c:4760-4767
        // ```c
        //       Blast_HSPInit(ungapped_data->q_start,
        //                     ungapped_data->length+ungapped_data->q_start,
        //     ...
        //                     context, query_info->contexts[context].frame,
        //                     subject->frame, ungapped_data->score, NULL, &new_hsp);
        // ```
        // The HSP keeps its subject frame and its linked-set size (`num`) for the
        // pairwise report.
        final_hits.push(TblastxHsp {
            hit: out_hit,
            subject_frame: s_frame.frame,
            num: h.num,
        });
    }
    x_lap.lap(XStage::Convert);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1018-1059
    // ```c
    // if (H->linked_set == TRUE && H->start_of_chain == FALSE)
    //     continue;
    // ...
    //             link = H->hsp_link.link[ordering_method];
    //             while (link)
    // ```
    // The converted output hits are represented here by their original
    // internal HSP rows after e-value filtering and before common.rs
    // applies final HSP-list sorting for the selected output format.
    if dump_output_stage {
        stage_dump::dump_ungapped_hits(
            "after_output_conversion_before_final_sort",
            &output_snapshot_hits,
        );
        stage_dump::dump_final_output_order("after_final_output_sort", &output_snapshot_pairs);
    }

    if debug_output_filter {
        // `[DEBUG OUTPUT_FILTER]` (`subject_debug::print_output_filter`).
        subject_debug::print_output_filter(
            &chain_roles,
            evalue_threshold,
            filtered_by_evalue,
            &final_hits,
        );
    }

    final_hits
}

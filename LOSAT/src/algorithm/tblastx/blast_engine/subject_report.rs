//! The report hits of a TBLASTX subject (`linked_hits_to_report`): its linked HSPs, reaped by
//! E-value, with their bit scores and identities in report coordinates, as `TblastxHsp`s, with the
//! stage dumps and the diagnostics of subject_debug.rs; and the stages after the subject's
//! re-evaluation that end in them (`report_subject_hits`: the HSP-saving statistics, the second
//! linking call, the report hits and the subject's final HSPs).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1535-1541 (in
//! `BLAST_PreliminarySearchEngine`: `s_Blast_HSPListReapByPrelimEvalue`, then
//! `Blast_HSPListGetBitScores`, c++/src/algo/blast/core/blast_hits.c:1907-1932), and
//! c++/src/algo/blast/api/blast_seqalign.cpp:1397 (`BLASTUngappedHspListToSeqAlign`, which turns
//! the HSPs into the report's alignments); for `report_subject_hits`,
//! c++/src/algo/blast/core/blast_engine.c:1515-1520 and 1539-1541 (after
//! `Blast_HSPListReevaluateUngapped`; the write of the list at :1553-1554 is
//! `SubjectEmit::emit`'s).
//!
//! This module is a child of `run_impl` for the timer it keeps private (`x_search_timing`). Its
//! functions are made of statements of `process_subject` moved here in their order; they take the
//! HSP list by value and read structs of references and `Copy` values.

use super::subject_debug;
use super::subject_link::{self, SubjectLinkArgs, SubjectStageArgs};
use super::x_search_timing::XLap;
use super::{
    calc_bit_score, convert_coords, stage_dump, trace_final_hit_if_match, AtomicOrdering, Hit,
    QueryFrame, TblastxHsp, TbxScore, UngappedHit, XStage,
};
use crate::algorithm::tblastx::x_link_timing::XLinkTimer;

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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1020
        // NOTE: This "continue" in NCBI is NOT a filter - NCBI then walks the link
        // pointer from each chain head to include all chain members (lines 1047-1059).
        // LOSAT's linking already produces a flat list of ALL HSPs (including chain
        // members) so we output everything without this skip logic.

        // NCBI: E-value filtering is applied during output conversion
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2720
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1918-1928
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

        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
            query_frame: ctx.frame as i32,
            query_length: 0,
            q_idx: ctx.q_idx,
            s_idx: h.s_idx,
            raw_score: h.raw_score,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4754-4767
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:808-812
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4760-4767
        // The HSP keeps its subject frame and its linked-set size (`num`) for the
        // pairwise report.
        final_hits.push(TblastxHsp {
            hit: out_hit,
            subject_frame: s_frame.frame,
            num: h.num,
        });
    }
    x_lap.lap(XStage::Convert);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1018-1059
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

/// What the stages after a subject's re-evaluation read besides its HSP list
/// (`report_subject_hits`), made per subject by `process_subject`: the values of the linking
/// calls and of the stages, the timer of the second linking call, the subject's frames
/// (`s_frames` scored the HSPs, `s_frames_report` the report reads) and its length in
/// nucleotides.
pub(super) struct SubjectReportArgs<'a> {
    pub(super) link: &'a SubjectLinkArgs<'a>,
    pub(super) stage: &'a SubjectStageArgs<'a>,
    pub(super) x_link_timer: &'a XLinkTimer,
    pub(super) s_frames: &'a [QueryFrame],
    pub(super) s_frames_report: &'a [QueryFrame],
    pub(super) s_len: usize,
}

/// The stages of a subject after its re-evaluation, when HSPs are left (`ungapped_hits`): the
/// HSP-saving statistics of a long subject (`hsp_saving`, `subject_debug::print_hsp_saving_stats`),
/// the second linking call (`subject_link::relink_subject_hits`) and the report hits
/// (`linked_hits_to_report`); when some are left, their LOSAT_DIAGNOSTICS counts and the
/// subject's final HSPs in `hits` (the worker's, sent at the end of the subject); then the `Emit`
/// lap of LOSAT_TIMING.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1515-1541 and 1539-1541
/// (after `Blast_HSPListReevaluateUngapped` in `BLAST_PreliminarySearchEngine`; the write of the
/// list that follows at :1553-1554 is `SubjectEmit::emit`'s)
pub(super) fn report_subject_hits<S: TbxScore>(
    sc: S,
    ungapped_hits: Vec<UngappedHit>,
    hsp_saving: subject_debug::HspSavingStats<'_>,
    report: &SubjectReportArgs<'_>,
    x_lap: &mut XLap<'_>,
    hits: &mut Vec<TblastxHsp>,
) {
    let &SubjectReportArgs {
        link,
        stage,
        x_link_timer,
        s_frames,
        s_frames_report,
        s_len,
    } = report;
    let &SubjectStageArgs {
        diag_enabled,
        diagnostics,
        ..
    } = stage;
    // DEBUG: the HSP-saving statistics of a long subject and LOSAT_DIAGNOSTICS'
    // `hsps_before_chain` (`subject_debug::print_hsp_saving_stats`).
    subject_debug::print_hsp_saving_stats(hsp_saving, &ungapped_hits, diag_enabled, diagnostics);

    // The second linking call (`subject_link::relink_subject_hits`), then the subject's
    // report hits (`subject_report::linked_hits_to_report`).
    let linked = subject_link::relink_subject_hits(link, ungapped_hits, stage, x_link_timer, x_lap);
    let final_hits =
        linked_hits_to_report(sc, linked, (s_frames, s_frames_report), s_len, stage, x_lap);

    if !final_hits.is_empty() {
        if diag_enabled {
            diagnostics
                .base
                .hsps_after_overlap_filter
                .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
            diagnostics
                .output_from_ungapped
                .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
        }
        // The subject's final HSPs, sent at the end of the subject (`SubjectEmit::emit`).
        debug_assert!(hits.is_empty());
        *hits = final_hits;
    }
    x_lap.lap(XStage::Emit);
}

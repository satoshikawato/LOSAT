//! The Smith-Waterman redo of one position-based match (`-use_sw_tback` with composition-based
//! statistics 1 to 3; RPS-BLAST's, series PR, item PS11c-3): NCBI's
//! `Blast_RedoOneMatchSmithWaterman` for the one window of a protein query and subject, then the
//! steps of the X-drop redo after it (`RedoMatch::redo_match`). The kernels are
//! `core::composition_adjustment::smith_waterman`'s, on the rows of the adjusted PSSM.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
//! ```c
//!     Boolean smithWaterman =
//!         (Boolean) (extendParams->options->eTbackExt == eSmithWatermanTbck);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3621
//! ```c
//!                     if (smithWaterman) {
//!                         *pStatusCode =
//!                                 Blast_RedoOneMatchSmithWaterman(
//!                                         alignments,
//!                                         redo_align_params,
//!                                         incoming_aligns,
//!                                         numAligns[frame_index],
//!                                         kbp->Lambda,
//!                                         kbp->logK,
//!                                         &matchingSeq,
//!                                         query_info,
//!                                         numQueries,
//!                                         matrix,
//!                                         BLASTAA_SIZE,
//!                                         NRrecord,
//!                                         forbidden,
//!                                         redoneMatches,
//!                                         &pvalueForThisPair,
//!                                         compositionTestIndex,
//!                                         &LambdaRatio
//!                                 );
//! ```
//! The steps after it are the X-drop redo's (blast_kappa.c:3654-3707, quoted at
//! `RedoMatch::redo_match`); they are this module's own copy so that psiblast's X-drop redo keeps
//! its code (DECISIONS 11 and 16 of line PS).

use super::{
    blast_nint, contained_in_hsp, preliminary_test_near_identical, score_compare, Align, RedoMatch,
    RedoRun, Redone, SCALING_FACTOR,
};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{BlastCompoHeap, BlastpHspList};
use crate::algorithm::pssm_search::gapalign::pssm_sw_xdrop_traceback;
use crate::algorithm::pssm_search::traceback0::final_hsp;
use crate::core::composition_adjustment::adjust_scores::read_aa_composition;
use crate::core::composition_adjustment::compo_scale::pssm_composition_based_stats_sw;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::core::composition_adjustment::redo_alignment::{
    test_near_identical, BlastCompoAlignment, BlastCompoSequenceData,
};
use crate::core::composition_adjustment::smith_waterman::{
    smith_waterman_find_start, smith_waterman_score_only, ForbiddenRanges, SwScratch,
};

/// What the Smith-Waterman redo keeps from one match to the next: the forbidden ranges (NCBI's
/// thread-local `forbidden_tld`), the scans' memory and the subject's letters between NCBI's two
/// sentinel bytes (`s_SequenceGetProteinRange`'s buffer).
pub(crate) struct SwRedoScratch {
    forbidden: ForbiddenRanges,
    scan: SwScratch,
    subject: Vec<u8>,
}

impl SwRedoScratch {
    pub(crate) fn new() -> Self {
        Self {
            forbidden: ForbiddenRanges::new(),
            scan: SwScratch::new(),
            subject: Vec::new(),
        }
    }
}

impl Default for SwRedoScratch {
    fn default() -> Self {
        Self::new()
    }
}

/// NCBI's redo of a match whose Smith-Waterman loop does not end: a later alignment of score 0
/// passes the e-value threshold and repeats the one before it (`Blast_RedoOneMatchSmithWaterman`
/// pushes a range that forbids nothing new, then finds the same alignment), which NCBI keeps until
/// its memory runs out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwEndless;

/// `s_EvalueFromScore`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:943-947
/// ```c
/// static double
/// s_EvalueFromScore(int score, double Lambda, double logK, double searchsp)
/// {
///     return searchsp * exp(-(Lambda * score) + logK);
/// }
/// ```
fn evalue_from_score(score: i32, lambda: f64, log_k: f64, searchsp: f64) -> f64 {
    searchsp * (-(lambda * f64::from(score)) + log_k).exp()
}

impl RedoMatch<'_> {
    /// The Smith-Waterman redo of one match (the HSPs of the query and one subject):
    /// `Blast_RedoOneMatchSmithWaterman`, then the X-drop redo's steps (`redo_match`).
    /// `heap` is the query's heap of redone matches, which the first alignment's significance
    /// asks (`BlastCompo_HeapWouldInsert`); `Ok(None)` when no HSP is left, also where NCBI's
    /// redo of the match fails (an alignment without an edit script).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1290-1305
    /// ```c
    ///     function_level_cleanup:
    ///     if (status != 0) {
    ///         for (query_index = 0;  query_index < numQueries;  query_index++) {
    ///             BlastCompo_AlignmentsFree(&alignments[query_index],
    ///                     callbacks->free_align_traceback);
    ///         }
    ///     }
    /// ```
    // One function for the window loop and the steps after it, called by RPS-BLAST's traceback
    // per profile where `-use_sw_tback` asks for it (its X-drop redo keeps its own code).
    #[inline(never)]
    pub(crate) fn redo_match_smith_waterman(
        &self,
        run: &mut RedoRun<'_>,
        sw: &mut SwRedoScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Result<Option<(BlastpHspList, i32, f64)>, SwEndless> {
        let Some(alignments) =
            self.redo_one_match_smith_waterman(run, sw, heap, oid, hsps, subject)?
        else {
            return Ok(None);
        };
        Ok(self.evaluate(run, oid, hsps, subject, alignments))
    }

    /// `Blast_RedoOneMatchSmithWaterman` for the one window of a protein query and subject (the
    /// whole query, the whole subject, `s_WindowsFromProteinAligns`): the window's alignments, the
    /// list's head first (the last one made), or `None` for NCBI's failed redo.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1366-1430
    /// ```c
    ///         window = windows[window_index];
    ///         query_index = window->query_range.context;
    ///         query_composition = &query_info[query_index].composition;
    ///
    ///         nearIdenticalStatus = s_preliminaryTestNearIdentical(query_info,
    ///                                              window,
    ///                                              window->align,
    ///                                              params->near_identical_cutoff);
    ///
    ///         status =
    ///             callbacks->get_range(matchingSeq, &window->subject_range,
    ///                                  &subject,
    ///                                  &query_info[query_index].seq,
    ///                                  &window->query_range,
    ///                                  &query,
    ///                                  query_info[query_index].words,
    ///                                  window->align, nearIdenticalStatus,
    ///                                  compo_adjust_mode, TRUE,
    ///                                  NULL);
    ///         ...
    ///         searchsp = query_info[query_index].eff_search_space;
    ///
    ///         /* For Smith-Waterman alignments, adjust the search using the
    ///          * composition of the highest scoring alignment in window */
    ///         if (compo_adjust_mode != eNoCompositionBasedStats) {
    ///             Blast_AminoAcidComposition subject_composition;
    ///             s_GetComposition(&subject_composition, alphsize,
    ///                                     &subject, &window->subject_range,
    ///                                     window->align, FALSE, subject_is_translated);
    ///             adjust_search_failed =
    ///                 Blast_AdjustScores(matrix,
    ///                                    query_composition, query.length,
    ///                                    &subject_composition, subject.length,
    ///                                    scaledMatrixInfo, compo_adjust_mode,
    ///                                    RE_pseudocounts, NRrecord,
    ///                                    &matrix_adjust_rule, callbacks->calc_lambda,
    ///                                    pvalueForThisPair,
    ///                                    compositionTestIndex,
    ///                                    LambdaRatio);
    ///             ...
    ///         }
    ///         if ( !adjust_search_failed ) {
    /// ```
    /// The subject range's last argument is `NULL`: the subject is SEG-masked unless the window's
    /// first alignment is near identical (`subject_range_sw`). The composition adjustment is the
    /// X-drop redo's (`RedoMatch::redo_one_match`, `Blast_AdjustScores` of a position-based
    /// search), on the copy `pssm_composition_based_stats_sw`.
    fn redo_one_match_smith_waterman(
        &self,
        run: &mut RedoRun<'_>,
        sw: &mut SwRedoScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Result<Option<Vec<Align>>, SwEndless> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3250-3256
        // ```c
        //         if(smithWaterman) {
        // 	   forbidden_tld[i] = calloc(1,sizeof(Blast_ForbiddenRanges));
        //            status_code =
        //                Blast_ForbiddenRangesInitialize(forbidden_tld[i], queryInfo->max_length);
        // ```
        let query_length = run.query.len();
        sw.forbidden.initialize(query_length);
        let Some(first) = hsps.first() else {
            return Ok(Some(Vec::new()));
        };
        // The window's first alignment (`window->align`), scaled as the X-drop redo's incoming
        // alignments (blast_kappa.c:788-797, quoted at `RedoMatch::redo_one_match`).
        let window_align = Align {
            score: (first.raw_score as f64 * SCALING_FACTOR) as i32,
            rule: EMatrixAdjustRule::DontAdjustMatrix,
            query_start: first.query_start,
            query_end: first.query_end,
            match_start: first.subject_start,
            match_end: first.subject_end,
            frame: first.subject_frame,
            edit_script: Vec::new(),
        };
        let hspcnt = hsps.len();
        let near_identical = preliminary_test_near_identical(
            query_length as i32,
            &window_align,
            run.near_identical_cutoff,
        );
        let subject_data = subject_range_sw(run, subject, &window_align, near_identical);
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1436-1442
        // ```c
        //     if (query_composition->numTrueAminoAcids == 0 ||
        //         subject_composition->numTrueAminoAcids == 0) {
        //         /* Either the query or subject contains only ambiguity
        //            characters, most likely because the entire subject has been
        //            SEGed.  Compositional adjustment is meaningless. */
        //         return 1;
        //     }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1484-1529
        // ```c
        //     if (matrixInfo->positionBased ||
        //         composition_adjust_mode == eCompositionBasedStats) {
        //         /* Use old-style composition-based statistics unconditionally. */
        //         *matrix_adjust_rule =  eCompoScaleOldMatrix;
        //     ...
        //     return Blast_CompositionBasedStats(matrix, ratioToPassBack, matrixInfo,
        //                                        query_composition->prob,
        //                                        subject_composition->prob,
        //                                        calc_lambda,
        //                                        (compositionTestIndex > 0));
        // ```
        let subject_composition = read_aa_composition(&subject_data);
        if run.query_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0 {
            return Ok(Some(Vec::new()));
        }
        let matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
        pssm_composition_based_stats_sw(&mut run.matrix.rows, self.info, &subject_composition.prob);
        // RPS-BLAST's bound of the adjusted matrix (`RedoRun::with_bound_check`): above `INT4_MAX`
        // NCBI's scores overflow and the caller refuses the search, so the scans do not run.
        if run.bound_check && run.bound_exceeded.is_none() {
            let bound: i64 = run
                .matrix
                .rows
                .iter()
                .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
                .sum();
            if bound > i64::from(i32::MAX) {
                run.bound_exceeded = Some(bound);
            }
        }
        if run.bound_exceeded.is_some() {
            return Ok(Some(Vec::new()));
        }
        // The subject's letters after one sentinel and before another (NCBI's
        // `seqData->buffer`, `seqData->data = buffer + 1`, blast_kappa.c:1601-1606 and 1640-1643):
        // `ALIGN_EX` reads `subject->data - 1`.
        sw.subject.clear();
        sw.subject.push(0);
        sw.subject.extend_from_slice(&subject_data);
        sw.subject.push(0);
        let letters = &sw.subject[1..=subject_data.len()];
        let rows = &run.matrix.rows;
        let searchsp = self.stats.search_space as f64;
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1431-1525
        // ```c
        //             Blast_ForbiddenRangesClear(forbidden);
        //             do {
        //                 int matchEnd, queryEnd;    /* end points of the alignments
        //                                             * computed by the Smith-Waterman
        //                                             * algorithm. */
        //                 status =
        //                     Blast_SmithWatermanScoreOnly(&aSwScore, &matchEnd,
        //                                                  &queryEnd,
        //                                                  subject.data,
        //                                                  subject.length,
        //                                                  query.data,
        //                                                  query.length, matrix,
        //                                                  gap_open, gap_extend,
        //                                                  positionBased,
        //                                                  forbidden);
        //                 ...
        //                 if (do_link_hsps) {
        //                     alignment_is_significant = aSwScore >= params->cutoff_s;
        //                 } else {
        //                     double newSwEvalue;     /* evalue as computed by the
        //                                              * Smith-Waterman algorithm */
        //                     newSwEvalue =
        //                         s_EvalueFromScore(aSwScore, Lambda, logK, searchsp);
        //
        //                     alignment_is_significant = newSwEvalue < params->cutoff_e;
        //                     if (alignments[query_index] == NULL) {
        //                         /* this is the most significant alignment; if
        //                          * it will not be accepted, no alignments from
        //                          * this match will */
        //                         alignment_is_significant =
        //                             alignment_is_significant &&
        //                             BlastCompo_HeapWouldInsert(
        //                                 &significantMatches[query_index],
        //                                 newSwEvalue, aSwScore, matchingSeq->index);
        //                     }
        //                 }
        //                 if (alignment_is_significant) {
        //                     ...
        //                     status =
        //                         Blast_SmithWatermanFindStart(&updatedScore,
        //                                                      &matchStart,
        //                                                      &queryStart,
        //                                                      subject.data,
        //                                                      subject.length,
        //                                                      query.data,
        //                                                      matrix, gap_open,
        //                                                      gap_extend,
        //                                                      matchEnd,
        //                                                      queryEnd,
        //                                                      aSwScore,
        //                                                      positionBased,
        //                                                      forbidden);
        //                     ...
        //                     status =
        //                         callbacks->
        //                         new_xdrop_align(&newAlign, &queryEnd, &matchEnd,
        //                                         queryStart, matchStart, aSwScore,
        //                                         &query, &window->query_range,
        //                                         ccat_query_length,
        //                                         &subject, &window->subject_range,
        //                                         matchingSeq->length,
        //                                         gapping_params, matrix_adjust_rule);
        //                     if (status != 0) {
        //                         goto window_index_loop_cleanup;
        //                     }
        //                     newAlign->next = alignments[query_index];
        //                     alignments[query_index] = newAlign;
        //
        //                     if (window->hspcnt > 1) {
        //                         /* We may compute more alignments; make the range
        //                            of the current alignment forbidden */
        //                         status =
        //                             Blast_ForbiddenRangesPush(forbidden,
        //                                                       queryStart, queryEnd,
        //                                                       matchStart, matchEnd);
        //                     }
        //                     ...
        //                 }
        //                 /* end if the next local alignment is significant */
        //             } while (alignment_is_significant && window->hspcnt > 1);
        // ```
        // `searchsp` is the query's `eff_search_space`, `Lambda` and `logK` the redo's scaled
        // block (`kbp_gap` of NCBI's copy), `cutoff_e` the `-evalue` (blast_kappa.c:2458). RPS-BLAST
        // refuses `-sum_stats true` and psiblast does not run this redo, so `do_link_hsps` is
        // false here.
        debug_assert!(self.link.is_none());
        sw.forbidden.clear();
        let mut alignments: Vec<Align> = Vec::new();
        loop {
            let end = smith_waterman_score_only(
                letters,
                query_length,
                |pos| &rows[pos],
                run.gap_open,
                run.gap_extend,
                &sw.forbidden,
                &mut sw.scan,
            );
            let evalue = evalue_from_score(end.score, run.karlin.lambda, run.log_k, searchsp);
            let mut alignment_is_significant = evalue < self.expect_value;
            if alignments.is_empty() {
                // RPS-BLAST's `matchingSeq->index` is -1 (its subject is the user's query,
                // blast_kappa.c:3576-3579); the empty heap of each profile's redo does not read it.
                alignment_is_significant =
                    alignment_is_significant && heap.would_insert(evalue, end.score, oid);
            }
            if !alignment_is_significant {
                break;
            }
            let start = smith_waterman_find_start(
                letters,
                |pos| &rows[pos],
                run.gap_open,
                run.gap_extend,
                end.match_end,
                end.query_end,
                end.score,
                &sw.forbidden,
                &mut sw.scan,
            );
            let traced = pssm_sw_xdrop_traceback(
                &run.matrix,
                &sw.subject,
                (start.query_start, end.query_end),
                (start.match_start, end.match_end),
                run.gap_open,
                run.gap_extend,
                run.x_dropoff,
                end.score,
                &mut run.scratch,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1856-1873
            // ```c
            //     if (editScript != NULL) {
            //         ...
            //         obj = BlastCompo_AlignmentNew(newScore, matrix_adjust_rule,
            //                                       aqueryStart, aqueryEnd,
            //                                       query_range->context,
            //                                       amatchStart, amatchEnd,
            //                                       subject_range->context, editScript);
            //         ...
            //     }
            //     *pnewAlign = obj;
            //
            //     return obj != NULL ? 0 : -1;
            // ```
            // An X-drop alignment without operations has no edit script
            // (`GapEditScriptNew(0)` is NULL): the redo of the match fails and frees the window's
            // alignments. The ranges start at 0 (the protein window's `query_range` and
            // `subject_range`) and its frame is 0.
            if traced.edit_script.is_empty() {
                return Ok(None);
            }
            let query_start = start.query_start as i32;
            let match_start = start.match_start as i32;
            let query_end = query_start + traced.query_extent;
            let match_end = match_start + traced.match_extent;
            alignments.push(Align {
                score: traced.score,
                rule: matrix_adjust_rule,
                query_start,
                query_end,
                match_start,
                match_end,
                frame: 0,
                edit_script: traced.edit_script,
            });
            if hspcnt <= 1 {
                break;
            }
            // Each alignment of a positive score starts at a cell that the scans did not forbid
            // and that its range forbids, so the loop ends. An alignment of score 0 (the start
            // and end at (0, 0), whose positive cell an earlier range forbids) forbids nothing
            // new: NCBI makes it again and again.
            if end.score == 0 {
                return Err(SwEndless);
            }
            sw.forbidden
                .push(query_start, query_end, match_start, match_end);
        }
        alignments.reverse();
        Ok(Some(alignments))
    }

    /// The steps of the redo after the window's alignments: `s_HSPListFromDistinctAlignments`,
    /// `s_HitlistReapContained`, `s_HitlistEvaluateAndPurge` and, for a list whose best e-value
    /// passes `-evalue`, `s_HSPListNormalizeScores` and `s_ComputeNumIdentities`; `RedoMatch::
    /// redo_match`'s without sum statistics (its quotes give NCBI's lines).
    fn evaluate(
        &self,
        run: &RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
        alignments: Vec<Align>,
    ) -> Option<(BlastpHspList, i32, f64)> {
        let mut redone: Vec<Redone> = alignments
            .into_iter()
            .map(|align| Redone {
                align,
                e_value: 0.0,
                num: 1,
            })
            .collect();
        redone.sort_by(|a, b| score_compare(&a.align, &b.align));
        if redone.len() > 1 {
            reap_contained_sw(&mut redone);
        }
        let subject_length = subject.len() as i32;
        let (m, n) = if self.rps_roles {
            (subject_length, self.query_length)
        } else {
            (self.query_length, subject_length)
        };
        let mut best_evalue = i32::MAX as f64;
        for hsp in &mut redone {
            hsp.e_value = self.stats.e_value(hsp.align.score, &run.karlin, m, n);
            best_evalue = best_evalue.min(hsp.e_value);
        }
        redone.retain(|hsp| hsp.e_value <= self.expect_value);
        let first = redone.first()?;
        let best_score = first.align.score;
        let template = &hsps[0];
        let mut hsp_list = BlastpHspList {
            oid,
            query_index: template.q_idx,
            hsps: Vec::with_capacity(redone.len()),
            best_evalue,
        };
        for Redone {
            align,
            e_value,
            num,
        } in redone
        {
            let score = blast_nint(align.score as f64 / SCALING_FACTOR);
            let traced = BlastpPreliminaryHsp {
                query_start: align.query_start,
                query_end: align.query_end,
                subject_start: align.match_start,
                subject_end: align.match_end,
                gapped_query_start: align.query_start,
                gapped_subject_start: align.match_start,
                raw_score: score,
                query_context: template.q_idx as i32,
                query_frame: template.query_frame,
                subject_frame: 0,
                query_length: template.query_length,
                q_idx: template.q_idx,
                s_idx: template.s_idx,
            };
            let mut hsp = final_hsp(
                &traced,
                align.edit_script,
                run.query_nomask,
                subject,
                self.matrix,
            );
            hsp.e_value = e_value;
            hsp.num = num;
            hsp.bit_score = (score as f64 * run.karlin.lambda * SCALING_FACTOR - run.log_k)
                / std::f64::consts::LN_2;
            hsp.comp_adjustment_method = match align.rule {
                EMatrixAdjustRule::DontAdjustMatrix => 0,
                EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                _ => 2,
            };
            hsp_list.hsps.push(hsp);
        }
        Some((hsp_list, best_score, best_evalue))
    }
}

/// `s_SequenceGetProteinRange` for the whole subject with NCBI's `NULL` for `subject_maybe_biased`
/// (`Blast_RedoOneMatchSmithWaterman`): its letters, masked by SEG ("S 10 1.8 2.1") unless the
/// window's first alignment is near identical. The X-drop redo's `RedoMatch::subject_range` does
/// the same for its first range; this is its copy, so that the X-drop redo keeps it inlined as
/// its only caller (DECISIONS 11 of line PS).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1638
/// ```c
///     if ( !(KAPPA_BLASTP_NO_SEG_SEQUENCE) ) {
///         if (compo_adjust_mode
///             && (!subject_maybe_biased || *subject_maybe_biased)) {
///
///             if ( (!shouldTestIdentical)
///                  || (shouldTestIdentical
///                      && (!s_TestNearIdentical(seqData, 0, queryData,
///                                               q_range->begin, query_words,
///                                               align)))) {
///
///                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
///                                              subject_maybe_biased);
///             }
///         }
///     }
/// ```
fn subject_range_sw(
    run: &RedoRun<'_>,
    subject: &[u8],
    window_align: &Align,
    near_identical: bool,
) -> Vec<u8> {
    let mut data = subject.to_vec();
    let test_identical = near_identical && {
        let align = BlastCompoAlignment {
            score: window_align.score,
            matrix_adjust_rule: window_align.rule,
            query_index: 0,
            query_start: window_align.query_start,
            query_end: window_align.query_end,
            match_start: window_align.match_start,
            match_end: window_align.match_end,
            frame: window_align.frame,
            context: None,
            next: None,
        };
        test_near_identical(
            &BlastCompoSequenceData::from_ncbistdaa(subject),
            0,
            &run.query_data,
            0,
            Some(&run.query_words),
            &align,
        )
    };
    if !test_identical {
        let intervals = crate::utils::seg::SegMasker::new(10, 1.8, 2.1).mask_sequence(&data);
        for interval in &intervals {
            data[interval.start..interval.end].fill(crate::utils::matrix::ncbistdaa::X);
        }
    }
    data
}

/// `s_HitlistReapContained` (quoted at `reap_contained`, blast_kappa.c:233-262): an HSP whose
/// ends lie in an earlier HSP and whose score is at most that HSP's is deleted. A copy of
/// `reap_contained`, so that the X-drop redo keeps it inlined as its only caller (DECISIONS 11 of
/// line PS).
fn reap_contained_sw(redone: &mut Vec<Redone>) {
    let mut keep = vec![true; redone.len()];
    for iread in 1..redone.len() {
        let hsp1 = &redone[iread].align;
        for ireadback in 0..iread {
            if !keep[ireadback] {
                continue;
            }
            let hsp2 = &redone[ireadback].align;
            if contained_in_hsp(
                hsp2.query_start,
                hsp2.query_end,
                hsp1.query_start,
                hsp2.match_start,
                hsp2.match_end,
                hsp1.match_start,
            ) && contained_in_hsp(
                hsp2.query_start,
                hsp2.query_end,
                hsp1.query_end,
                hsp2.match_start,
                hsp2.match_end,
                hsp1.match_end,
            ) && hsp1.score <= hsp2.score
            {
                keep[iread] = false;
                break;
            }
        }
    }
    let mut index = 0;
    redone.retain(|_| {
        index += 1;
        keep[index - 1]
    });
}

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;

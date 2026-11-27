//! The Smith-Waterman redo of one position-based match (`-use_sw_tback`): NCBI's
//! `Blast_RedoOneMatchSmithWaterman` for the one window of a protein query and subject, then the
//! steps of the X-drop redo after it (`RedoMatch::redo_match`), for psiblast's PSSM (modes 0 to 3,
//! line PS-SW: `PssmSwStep`, the step of `PssmRedo1::redo_loop`) and RPS-BLAST's profiles (modes 1
//! to 3, series PR, item PS11c-3); for RPS-TBLASTN's translated subject (item PS12c) the same
//! function over the windows of the frames, then the steps of the X-drop redo of frames
//! (`RedoMatch::evaluate_windows`). The kernels are `core::composition_adjustment::smith_waterman`'s,
//! on the rows of the adjusted PSSM, or, without composition-based statistics, of the PSSM as the
//! search made it.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3621
//! The steps after it are the X-drop redo's (blast_kappa.c:3654-3707, cited at
//! `RedoMatch::redo_match`). The functions the two redos share run as this redo's instance
//! (`compo_scale::SwRedo`), so that the X-drop redo's instance keeps one caller and psiblast's
//! machine code (DECISIONS 11 and 63 of line PS); the evaluation is a copy (`evaluate`).

use super::{
    blast_nint, frame_range, frame_windows, link_redone, old_test_near_identical,
    preliminary_test_near_identical, reap_contained, score_compare, window_composition, Align,
    FrameSubject, MatchStep, PssmRedo1, PssmSwOptions, RedoMatch, RedoRun, Redone,
};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{BlastCompoHeap, BlastpHspList};
use crate::algorithm::pssm_search::gapalign::pssm_sw_xdrop_traceback;
use crate::algorithm::pssm_search::traceback0::final_hsp;
use crate::core::composition_adjustment::adjust_scores::read_aa_composition;
use crate::core::composition_adjustment::compo_scale::{
    pssm_composition_based_stats, SwFrameRedo, SwRedo,
};
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::core::composition_adjustment::smith_waterman::{
    smith_waterman_redo_window, SwRedoScratch, SwSignificance, SwWindow, SwWindowRedo,
};
use crate::core::composition_adjustment::smith_waterman_match::{
    SwBatchStatus, SwRedoEndless, SwRedoFailed,
};

/// The memory of the position-based Smith-Waterman redo from one match to the next (psiblast's
/// subjects, RPS-BLAST's profiles): the round loop's (`SwRedoScratch::rounds`) and the subject's
/// letters between NCBI's two sentinel bytes (`s_SequenceGetProteinRange`'s buffer;
/// `SwRedoScratch::window`).
pub(crate) type PssmSwScratch = SwRedoScratch<Vec<u8>>;

/// NCBI's redo of a match whose Smith-Waterman loop does not end: the score of 0 is significant,
/// the window has two or more HSPs and profile row 0 scores the window's first letter positive.
/// After every positive alignment NCBI makes the alignment of that one cell again and again; its
/// range changes no score that the scan sees, and the loop stops only when an allocation fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwEndless;

/// A position-based match whose Smith-Waterman redo gives no alignments
/// ([`RedoMatch::redo_match_smith_waterman`]): NCBI's failed redo (`Failed`: a significant
/// alignment of score 0 without an edit script, which frees the match's alignments and sets the
/// status) or its endless one (`Endless`, as `SwEndless`). Each caller maps it: RPS-BLAST drops a
/// failed profile and refuses an endless one (`for_rps`); psiblast keeps the batch's status
/// (`PssmSwStep`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SwMatchFailure {
    Failed,
    Endless,
}

impl SwMatchFailure {
    /// RPS-BLAST's reading of the failure: a failed profile has no list (NCBI drops its
    /// alignments; RPS-BLAST's redo keeps no status), an endless one is `SwEndless`.
    #[inline(always)]
    pub(crate) fn for_rps<T>(self) -> Result<Option<T>, SwEndless> {
        match self {
            SwMatchFailure::Failed => Ok(None),
            SwMatchFailure::Endless => Err(SwEndless),
        }
    }
}

/// psiblast's step of `PssmRedo1::redo_loop` with `-use_sw_tback`: the subject read, the
/// Smith-Waterman redo of the match (`RedoMatch::redo_match_smith_waterman`) and NCBI's status of
/// the batch (`SwBatchStatus`: set by each redone match, cleared at a subject that the fetch
/// excludes, `PssmSwOptions::on_query_done`; a failed status at the batch's end is its error,
/// `PssmSwOptions::on_batch_failed`). A match whose redo NCBI never ends is refused where it comes
/// (`PssmSwOptions::on_endless`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3646
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
pub(super) struct PssmSwStep<'a> {
    options: &'a PssmSwOptions<'a>,
    scratch: PssmSwScratch,
    status: SwBatchStatus,
    /// The position in the query's list of the last match whose redo set the status.
    status_position: Option<usize>,
}

impl<'a> PssmSwStep<'a> {
    pub(super) fn new(options: &'a PssmSwOptions<'a>) -> Self {
        Self {
            options,
            scratch: PssmSwScratch::new(),
            status: SwBatchStatus::default(),
            status_position: None,
        }
    }
}

impl MatchStep for PssmSwStep<'_> {
    fn redo(
        &mut self,
        redo: &PssmRedo1<'_>,
        run: &mut RedoRun<'_>,
        heaps: &mut [BlastCompoHeap],
        q_idx: usize,
        position: usize,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
    ) -> anyhow::Result<Option<(BlastpHspList, i32, f64)>> {
        let subject_read = redo.subjects.read(oid)?;
        let redone = redo
            .matching()
            .redo_match_smith_waterman(
                run,
                &mut self.scratch,
                &mut heaps[q_idx],
                oid,
                hsps,
                subject_read.residues(),
            )
            .map_err(|failure| match failure {
                SwMatchFailure::Failed => anyhow::Error::from(SwRedoFailed),
                SwMatchFailure::Endless => anyhow::Error::from(SwRedoEndless),
            });
        let on_endless = self.options.on_endless;
        let matched = self
            .status
            .redone(redone, || on_endless(&*heaps, position))?;
        self.status_position = Some(position);
        Ok(matched.flatten())
    }

    fn query_done(&mut self, heaps: &[BlastCompoHeap]) {
        (self.options.on_query_done)(self.status_position.take(), heaps, &mut self.status);
    }

    fn batch_done(&mut self, heaps: &[BlastCompoHeap]) -> anyhow::Result<()> {
        if self.status.failed() {
            return Err((self.options.on_batch_failed)(heaps));
        }
        Ok(())
    }
}

impl RedoMatch<'_> {
    /// The Smith-Waterman redo of one match (the HSPs of the query and one subject):
    /// `Blast_RedoOneMatchSmithWaterman`, then the X-drop redo's steps (`redo_match`), at the
    /// run's scale (`RedoRun::scale`). `heap` is the query's heap of redone matches, which the
    /// first alignment's significance asks (`BlastCompo_HeapWouldInsert`); `Ok(None)` when no HSP
    /// is left, `SwMatchFailure` where NCBI's redo of the match fails or does not end.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1542-1548
    // One function for the window loop and the steps after it, called by psiblast's redo loop and
    // by RPS-BLAST's traceback per profile where `-use_sw_tback` asks for it (their X-drop redos
    // keep their own code).
    #[inline(never)]
    pub(crate) fn redo_match_smith_waterman(
        &self,
        run: &mut RedoRun<'_>,
        sw: &mut PssmSwScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Result<Option<(BlastpHspList, i32, f64)>, SwMatchFailure> {
        let alignments = self.redo_one_match_smith_waterman(run, sw, heap, oid, hsps, subject)?;
        Ok(self.evaluate(run, oid, hsps, subject, alignments))
    }

    /// `Blast_RedoOneMatchSmithWaterman` for the one window of a protein query and subject (the
    /// whole query, the whole subject, `s_WindowsFromProteinAligns`): the window's alignments, the
    /// list's head first (the last one made). With composition-based statistics
    /// (`RedoScale::compo_adjust`) the subject range is SEG-masked unless the window's first
    /// alignment is near identical and the matrix is adjusted to the range's composition; without
    /// them (psiblast's `-comp_based_stats 0`) neither runs and the rows are the PSSM's as the run
    /// holds them (the search's, unscaled).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1366-1430
    /// The subject range's last argument is `NULL`: SEG runs (with composition-based statistics)
    /// unless the window's first alignment is near identical (`RedoMatch::subject_range`). The
    /// composition adjustment is the X-drop redo's (`RedoMatch::redo_one_match`,
    /// `Blast_AdjustScores` of a position-based search), with this redo's instance (`SwRedo`) of
    /// the shared steps.
    fn redo_one_match_smith_waterman(
        &self,
        run: &mut RedoRun<'_>,
        sw: &mut PssmSwScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Result<Vec<Align>, SwMatchFailure> {
        let query_length = run.query.len();
        let Some(first) = hsps.first() else {
            return Ok(Vec::new());
        };
        let hspcnt = hsps.len();
        let (matrix_adjust_rule, subject_masked);
        let subject_slice: &[u8] = if run.scale.compo_adjust {
            // The window's first alignment (`window->align`), scaled as the X-drop redo's incoming
            // alignments (blast_kappa.c:788-797, cited at `RedoMatch::redo_one_match`).
            let window_align = Align {
                score: (first.raw_score as f64 * run.scale.factor) as i32,
                rule: EMatrixAdjustRule::DontAdjustMatrix,
                query_start: first.query_start,
                query_end: first.query_end,
                match_start: first.subject_start,
                match_end: first.subject_end,
                frame: first.subject_frame,
                edit_script: Vec::new(),
            };
            let near_identical = preliminary_test_near_identical(
                query_length as i32,
                &window_align,
                run.near_identical_cutoff,
            );
            // NCBI's `NULL` for `subject_maybe_biased`: SEG runs whenever the pointer is null or
            // points to true, and its answer is not kept, as for the X-drop redo's first range.
            let mut subject_maybe_biased = true;
            subject_masked = self.subject_range::<SwRedo>(
                run,
                subject,
                &window_align,
                near_identical,
                &mut subject_maybe_biased,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1436-1444
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1484-1530
            let subject_composition = read_aa_composition(&subject_masked);
            if run.query_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0 {
                return Ok(Vec::new());
            }
            matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
            pssm_composition_based_stats::<SwRedo>(
                &mut run.matrix.rows,
                self.info,
                &subject_composition.prob,
            );
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
                return Ok(Vec::new());
            }
            &subject_masked
        } else {
            // Without composition-based statistics the range is the subject's letters (no SEG:
            // `compo_adjust_mode` is 0), no adjustment runs and the rule is `eDontAdjustMatrix`.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1625-1639
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1410-1411
            matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
            subject
        };
        // The subject's letters after one sentinel and before another (NCBI's
        // `seqData->buffer`, with `seqData->data` set one past it, blast_kappa.c:1601-1606 and
        // 1640-1643): `ALIGN_EX` reads from one byte before `subject->data`.
        let subject_buffer = &mut sw.window;
        subject_buffer.clear();
        subject_buffer.push(0);
        subject_buffer.extend_from_slice(subject_slice);
        subject_buffer.push(0);
        let rows = &run.matrix.rows;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3261-3264
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1462-1474
        let window = SwWindow {
            subject_from_sentinel: subject_buffer.as_slice(),
            subject_length: subject_slice.len(),
            query_length,
            rows: |pos| &rows[pos],
            gap_open: run.gap_open,
            gap_extend: run.gap_extend,
            hspcnt,
        };
        let significance = match self.link {
            Some(link) => SwSignificance::CutoffScore(link.cutoff_s),
            None => SwSignificance::Evalue {
                lambda: run.karlin.lambda,
                log_k: run.log_k,
                searchsp: self.stats.search_space as f64,
                cutoff_e: self.expect_value,
            },
        };
        let (matrix, gap, x_dropoff) = (&run.matrix, &mut run.scratch, run.x_dropoff);
        let (gap_open, gap_extend) = (run.gap_open, run.gap_extend);
        // One window per profile/query: the query has no alignment at its start.
        let redone = smith_waterman_redo_window(
            &window,
            significance,
            true,
            |evalue, score| heap.would_insert(evalue, score, oid),
            |query, subject_span, score| {
                pssm_sw_xdrop_traceback(
                    matrix,
                    window.subject_from_sentinel,
                    query,
                    subject_span,
                    gap_open,
                    gap_extend,
                    x_dropoff,
                    score,
                    gap,
                )
            },
            &mut sw.rounds,
        );
        // The protein window's ranges start at 0 (`query_range`, `subject_range`) and its frame
        // is 0; NCBI's list starts with the last alignment made.
        let mut alignments: Vec<Align> = match redone {
            SwWindowRedo::Alignments(alignments) => alignments
                .into_iter()
                .map(|align| Align {
                    score: align.score,
                    rule: matrix_adjust_rule,
                    query_start: align.query_start,
                    query_end: align.query_end,
                    match_start: align.match_start,
                    match_end: align.match_end,
                    frame: 0,
                    edit_script: align.edit_script,
                })
                .collect(),
            SwWindowRedo::Failed => return Err(SwMatchFailure::Failed),
            SwWindowRedo::Endless => return Err(SwMatchFailure::Endless),
        };
        alignments.reverse();
        Ok(alignments)
    }

    /// The Smith-Waterman redo of one match of a translated subject (rpstblastn's
    /// `-use_sw_tback`, series PR, item PS12c): `Blast_RedoOneMatchSmithWaterman` over the windows
    /// of the frames, then the steps of the X-drop redo of frames after it
    /// (`evaluate_windows`). `sw` is the round loop's memory (its `window` is not used: each window
    /// has its frame's buffer). `heap` is the query's heap of redone matches, which the match's first
    /// alignment asks (`BlastCompo_HeapWouldInsert`); `Ok(None)` when no HSP is left, also where
    /// NCBI's redo of the match fails (an alignment without an edit script, in any window).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1353-1548
    /// The windows are those of the X-drop redo of frames (`s_WindowsFromAligns` of a translated
    /// subject with a position-based matrix, `frame_windows`). Each window is redone alone: the
    /// near-identical test of its first alignment `window->align` (`s_preliminaryTestNearIdentical`
    /// with the window's HSP count; the old test where the cutoff is 0,
    /// `RedoRun::without_valid_context`), its frame's range (`frame_range`, SEG unless near
    /// identical, NCBI's `NULL` for `subject_maybe_biased`), the composition of the window around
    /// that alignment (`s_GetComposition`, `window_composition`) and the adjustment
    /// (`Blast_AdjustScores`; a window without a true amino acid makes no alignment), the forbidden
    /// ranges cleared and the round loop on the window's letters (`subject.data` after the zero
    /// byte that `frame_range` writes before the window, `subject.length` its `end - begin`).
    /// The alignments are shifted by the window's start and take its frame
    /// (`s_NewAlignmentUsingXdrop`). All windows prepend to the match's one
    /// list (`alignments[query_index]`), and the heap is asked only while that list is empty.
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn redo_windows_smith_waterman(
        &self,
        run: &mut RedoRun<'_>,
        sw: &mut PssmSwScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
        frames: FrameSubject<'_>,
    ) -> Result<Option<(BlastpHspList, i32, f64)>, SwEndless> {
        debug_assert!(self.link.is_none());
        let query_length = run.query.len();
        let matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
        // The match's alignments in the order the windows' rounds made them.
        let mut made: Vec<Align> = Vec::new();
        for (window, order) in
            frame_windows::<SwFrameRedo>(hsps, query_length as i32, subject.len() as i32)
        {
            let hspcnt = order.len();
            let window_align = Align::incoming::<SwFrameRedo>(&hsps[order[0]]);
            let near_identical = if run.near_identical_cutoff > 0.0 {
                preliminary_test_near_identical(
                    query_length as i32,
                    &window_align,
                    run.near_identical_cutoff,
                )
            } else {
                old_test_near_identical(query_length as i32, hspcnt, &window_align)
            };
            let mut subject_maybe_biased = true;
            let subject_data = frame_range::<SwFrameRedo>(
                run,
                subject,
                window,
                &window_align,
                near_identical,
                &mut subject_maybe_biased,
            );
            let subject_composition =
                window_composition::<SwFrameRedo>(&subject_data, window, &window_align);
            if run.query_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0 {
                continue;
            }
            pssm_composition_based_stats::<SwFrameRedo>(
                &mut run.matrix.rows,
                self.info,
                &subject_composition.prob,
            );
            // RPS-BLAST's bound of the adjusted matrix (`RedoRun::with_bound_check`), as in
            // `redo_one_match_smith_waterman`: above `INT4_MAX` the caller refuses the search.
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
                return Ok(None);
            }
            // A window that ends before it begins (`-strand plus|minus`, NCBI's six-frame
            // arithmetic on the three-frame buffer; `frame_range`) does not get here: its end is
            // the frame's translated length, it holds one HSP, and every letter its composition
            // reads lies past the frame's copy (zero bytes, no true amino acid), so the test above
            // skips it as NCBI's does. The check keeps `end - begin` from wrapping; were it
            // reached, NCBI's scan could not allocate a row of a negative length and would fail
            // the match (as `SwWindowRedo::Failed`).
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:681-682
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1430
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:74-77
            if window.end < window.begin {
                return Ok(None);
            }
            // `frame_range`'s buffer holds the frame's letters from byte 1 and a zero byte before
            // the window (`buffer[begin]`): NCBI's `seqData->data` after the fit, whose byte before
            // `ALIGN_EX` reads.
            let rows = &run.matrix.rows;
            let window_letters = SwWindow {
                subject_from_sentinel: &subject_data[window.begin as usize..],
                subject_length: (window.end - window.begin) as usize,
                query_length,
                rows: |pos| &rows[pos],
                gap_open: run.gap_open,
                gap_extend: run.gap_extend,
                hspcnt,
            };
            let significance = SwSignificance::Evalue {
                lambda: run.karlin.lambda,
                log_k: run.log_k,
                searchsp: self.stats.search_space as f64,
                cutoff_e: self.expect_value,
            };
            // NCBI asks the heap only for the match's first alignment, while
            // `alignments[query_index]` is still `NULL`: a window after one that made alignments
            // starts with the query's list not empty (`query_had_none` false), as TBLASTN's.
            let query_had_none = made.is_empty();
            let (matrix, gap, x_dropoff) = (&run.matrix, &mut run.scratch, run.x_dropoff);
            let (gap_open, gap_extend) = (run.gap_open, run.gap_extend);
            let redone = smith_waterman_redo_window(
                &window_letters,
                significance,
                query_had_none,
                |evalue, score| heap.would_insert(evalue, score, oid),
                |query, subject_span, score| {
                    pssm_sw_xdrop_traceback(
                        matrix,
                        window_letters.subject_from_sentinel,
                        query,
                        subject_span,
                        gap_open,
                        gap_extend,
                        x_dropoff,
                        score,
                        gap,
                    )
                },
                &mut sw.rounds,
            );
            match redone {
                SwWindowRedo::Alignments(alignments) => {
                    made.extend(alignments.into_iter().map(|align| Align {
                        score: align.score,
                        rule: matrix_adjust_rule,
                        query_start: align.query_start,
                        query_end: align.query_end,
                        match_start: align.match_start + window.begin,
                        match_end: align.match_end + window.begin,
                        frame: window.frame,
                        edit_script: align.edit_script,
                    }))
                }
                SwWindowRedo::Failed => return Ok(None),
                SwWindowRedo::Endless => return Err(SwEndless),
            }
        }
        // NCBI's list starts with the last alignment made.
        made.reverse();
        Ok(self.evaluate_windows::<SwFrameRedo, FrameSubject<'_>>(
            run, oid, hsps, subject, frames, made,
        ))
    }

    /// The steps of the redo after the window's alignments: `s_HSPListFromDistinctAlignments`,
    /// `s_HitlistReapContained`, `s_HitlistEvaluateAndPurge` (with psiblast's `-sum_stats true`,
    /// the link of the redone alignments, `link_redone`) and, for a list whose best e-value passes
    /// `-evalue`, `s_HSPListNormalizeScores` and `s_ComputeNumIdentities`, at the run's scale
    /// (`RedoRun::scale`); `RedoMatch::redo_match`'s (its quotes give NCBI's lines). These lines
    /// are a copy of `redo_match`'s, which keeps them in its body (moved into a function shared by
    /// both redos they changed psiblast's machine code; line PS `DECISIONS.md` 63); the test
    /// `the_evaluation_is_the_x_drop_redos` holds the two together, with and without the link.
    fn evaluate(
        &self,
        run: &RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
        alignments: Vec<Align>,
    ) -> Option<(BlastpHspList, i32, f64)> {
        let scaling_factor = run.scale.factor;
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
            reap_contained::<SwRedo>(&mut redone);
        }
        let subject_length = subject.len() as i32;
        let best_evalue = match &self.link {
            Some(link) => link_redone(&mut redone, link, run, subject_length)?,
            None => {
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
                best_evalue
            }
        };
        let best_score = self.reap_by_expect(&mut redone, best_evalue)?;
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
            let score = blast_nint(align.score as f64 / scaling_factor);
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
            hsp.bit_score = (score as f64 * run.karlin.lambda * scaling_factor - run.log_k)
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

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;

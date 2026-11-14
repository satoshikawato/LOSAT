//! The preparation of a TBLASTX subject before its frame loop: its preliminary frames
//! (`subject_frames`), and the unit path (PAIRPAR), which searches the subject's frames as units on
//! the search pool when the batch has one and then leaves the frame loop no frame
//! (`run_unit_path`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1429-1432 (the subject loop
//! of `BLAST_PreliminarySearchEngine` skips a subject shorter than `s_GetMinimumSubjSeqLen`,
//! :1318-1330) and :772-775 (`s_BlastSearchEngineCore` translates the subject with
//! `BLAST_GetAllTranslations`, c++/src/algo/blast/core/blast_util.c:1045). The unit path has no
//! NCBI counterpart: it returns the HSP list of NCBI's frame loop, blast_engine.c:805-855
//! (tblastx/x_pair_par.rs).
//!
//! This module is a child of `run_impl` for the modules it keeps private (`x_pair_par`,
//! `x_seed_release`, `x_batch_pipe`, `x_search_timing`). Its functions are statements of
//! `process_subject` moved here in their order. The subject's guards stay in `process_subject`,
//! which lends `x_lap` and `x_pipe_seed` and hands its units stage over by value; the returns for a
//! subject that is not searched stay there too.

use super::x_batch_pipe::{XSeedEnd, XSeedGuard};
use super::x_pair_par;
use super::x_search_timing::XLap;
use super::x_seed_release::XSeedStage;
use super::{
    GeneticCode, QueryContext, QueryFrame, Result, Subjects, TbxScore, UngappedHit, XStage,
};
use crate::algorithm::tblastx::subjects::XKept;
use crate::algorithm::tblastx::x_work_board::{XBatch, XBoard};
use crate::utils::threading::SearchPool;
use std::sync::Arc;

/// The frames of a subject long enough to be searched, which `process_subject` reads under the
/// same names: its number of bases, the frames of its letters (empty until they are made after
/// the seed stage, `x_letter_frames_later`), its preliminary frames when they are not those, and
/// what a `-db` search keeps for the frames of its letters.
pub(super) struct SubjectFrames {
    pub(super) s_len: usize,
    pub(super) s_frames: Vec<QueryFrame>,
    pub(super) s_frames_preliminary_owned: Option<Vec<QueryFrame>>,
    pub(super) x_kept: XKept,
    pub(super) x_letter_frames_later: bool,
}

/// The preliminary frames of subject `s_idx` (`Subjects::x_preliminary_frames`, the genetic code
/// `db_code`); `None` for a subject shorter than `wordsize * 3 + 2` bases, which is not searched.
/// A read of a `-db` subject that fails is the error.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
/// (`BLAST_GetAllTranslations`, c++/src/algo/blast/core/blast_util.c:1045), :1318-1330
/// (`s_GetMinimumSubjSeqLen`) and :1429-1432 (the subject loop skips a shorter subject).
pub(super) fn subject_frames(
    subjects: Subjects<'_>,
    s_idx: usize,
    db_code: &GeneticCode,
    wordsize: i32,
) -> Result<Option<SubjectFrames>> {
    // NCBI reference (translate all frames for subject sequences):
    // ncbi-blast/c++/src/algo/blast/core/blast_util.c:1296-1308
    // ```c
    // for (context = 0; context < num_frames; ++context) {
    //    int frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //    retval->translations[context] = (Uint1*) malloc(...);
    //    BLAST_GetTranslation(subject_blk->sequence_start, nucl_seq_rev,
    //       subject_blk->length, frame, retval->translations[context], gen_code_string);
    // }
    // ```
    // Only the frames the preliminary stage translates are made here, in windows
    // (`Subjects::x_preliminary_frames`). For a subject with ambiguous bases the frames of
    // its letters (`s_frames`) are made after the seed stage, once the preliminary frames
    // are freed, and only when HSPs are left for the re-evaluation (`x_letter_frames_later`,
    // below); until then `s_frames` is empty and the HSP lists index the preliminary frames,
    // which have the same number and lengths (`process_subject`'s `s_frames_listed`). For a
    // subject without ambiguous bases both are the same frames.
    let mut x_letter_frames_later = false;
    // `-db`: the block of a short sequence with ambiguous bases, kept for the frames of its
    // letters (`Subjects::x_letter_frames`).
    let x_kept;
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1318-1330,1429-1431
    // return word_length * 3 + 2;
    // if (subject->length < min_subj_seq_length) { ... continue; }
    // Statistics include all records, including ones too short to search: such a subject
    // gets no frame (`None`). A `-db` sequence's length is read with its block, in the
    // frames' one look-up (it holds no length per sequence).
    let s_len;
    let (s_frames, s_frames_preliminary_owned) = {
        let (length, frames, of_letters, kept) =
            match subjects.x_preliminary_frames(s_idx, db_code, (wordsize * 3 + 2) as usize)? {
                Some(made) => made,
                None => return Ok(None),
            };
        s_len = length;
        x_kept = kept;
        if of_letters {
            (frames, None)
        } else {
            x_letter_frames_later = true;
            (Vec::new(), Some(frames))
        }
    };
    Ok(Some(SubjectFrames {
        s_len,
        s_frames,
        s_frames_preliminary_owned,
        x_kept,
        x_letter_frames_later,
    }))
}

/// What the unit path of a subject reads, made per subject by `process_subject` under the names
/// it reads them by: the subject's frames (the preliminary ones, and the ones its HSP lists
/// index), the batch's contexts with their cutoffs and X-drops, the window and word size, the
/// subject's index and length, the search pool, the batch on the board and the batch pipe's
/// seed-stage end, and whether the batch's lookup table is PAIRLUT's stub.
pub(super) struct UnitPathArgs<'a> {
    pub(super) s_frames_preliminary: &'a [QueryFrame],
    pub(super) s_frames_listed: &'a [QueryFrame],
    pub(super) contexts_ref: &'a [QueryContext],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) x_dropoff_per_context: &'a [i32],
    pub(super) window: i32,
    pub(super) wordsize: i32,
    pub(super) s_idx: usize,
    pub(super) s_len: usize,
    pub(super) parallel_pool: &'a SearchPool<'a>,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): the search's board and the batch on it.
    pub(super) x_board_batch: &'a Option<(&'a XBoard, Arc<XBatch>)>,
    /// EXPERIMENT (LOSAT_X_TBXPIPE): the batch pipe's view of this search's seed stage.
    pub(super) x_seed_end: Option<&'a XSeedEnd<'a>>,
    /// PAIRLUT: the batch's lookup table is a stub without entries.
    pub(super) x_pair_lut_is_stub: bool,
}

/// Dispatch points 2 and 3 of PAIRPAR and of the seed release (per subject): when the batch has a
/// unit path (`x_units_stage`'s tables), the subject's units run on the search pool, the units
/// stage ends and their lists are merged into `combined_ungapped_hits`; the frames the frame loop
/// runs over are then none. Without a unit path the stage ends and the frame loop gets the
/// preliminary frames. Stops the run when the lookup table is PAIRLUT's stub and the unit path did
/// not take the subject. The laps `Prepare`, `Units`, `UnitsRelease` and `Merge` go to `x_lap`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-855 (the frame loop of
/// `s_BlastSearchEngineCore`, whose list the unit path returns) and :561-562
/// (`BLAST_GetUngappedHSPList`, the conversion that the merge makes). The split into units and
/// the stage end between them have no NCBI counterpart (tblastx/x_pair_par.rs,
/// tblastx/x_seed_release.rs).
pub(super) fn run_unit_path<'a, S: TbxScore>(
    args: &UnitPathArgs<'a>,
    x_units_stage: XSeedStage<'_, Option<x_pair_par::XPairPar<S>>>,
    x_lap: &mut XLap<'_>,
    x_pipe_seed: &mut XSeedGuard<'_>,
    combined_ungapped_hits: &mut Vec<UngappedHit>,
) -> &'a [QueryFrame] {
    let &UnitPathArgs {
        s_frames_preliminary,
        s_frames_listed,
        contexts_ref,
        cutoff_scores,
        x_dropoff_per_context,
        window,
        wordsize,
        s_idx,
        s_len,
        parallel_pool,
        x_board_batch,
        x_seed_end,
        x_pair_lut_is_stub,
    } = args;
    // The batch's unit path, borrowed from the subject's units stage (dispatch point 2 of the
    // seed release): it cannot be used after `x_units_stage.end()` (dispatch point 3).
    let x_pair_par_ref: Option<&x_pair_par::XPairPar<S>> = x_units_stage.tables().as_ref();
    // EXPERIMENT (PAIRPAR), dispatch point 2 (per subject): the unit path fills
    // `combined_ungapped_hits` and the frame slice this fn returns is empty, so the reference
    // frame loop (unchanged) runs zero times.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805,835-844
    // ```c
    // for (context=first_context; context<=last_context; context++) {
    // ...
    //     status = s_BlastSearchEngineOneContext(program_number, query, query_info,
    // ...
    //     if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
    // ```
    // The unit path returns the list this frame loop builds (same calls, same order of
    // the appended frame lists); see x_pair_par.rs for why it is the same list.
    //
    // Dispatch point 3 of the seed release (per subject): the unit path runs
    // in two calls, the units (`run_subject_units`, the only reader of the group tables)
    // and the merge of their lists (`merge_subject`, a function without `self`), and the
    // units stage of this subject ends between them. What is still needed of the unit path
    // is read before that: `x_pair_par_stub` for the assertion below. After
    // `x_units_stage.end()` neither `x` nor `x_pair_par_ref` can be used. The batch's last
    // units stage frees the unit path there, before the merge, where the seed stage's
    // memory peaks.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-562
    // ```c
    // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
    //         hit_params->options, &hsp_list);
    // ```
    // NCBI converts and merges its init HSPs with the lookup table still allocated (freed
    // with the search's internal data); that conversion reads no lookup table.
    let x_pair_par_stub = x_pair_par_ref.is_some_and(|x| x.stub());
    let s_frames_preliminary: &[QueryFrame] = match x_pair_par_ref {
        Some(x) => {
            let x_subject = x_pair_par::XSubject {
                s_frames_preliminary,
                s_frames: s_frames_listed,
                contexts: contexts_ref,
                cutoff_scores,
                x_dropoff_per_context,
                window,
                wordsize,
                s_idx,
                s_len,
            };
            x_lap.lap(XStage::Prepare);
            let x_units = x.run_subject_units(&x_subject, parallel_pool);
            x_lap.lap(XStage::Units);
            x_pipe_seed.units_done(x_units.records());
            // EXPERIMENT (LOSAT_X_TBXBOARD): every unit of the batch taken: the next batch
            // may build its tables (tblastx/x_batch_pipe.rs `XSeedEnd::units_taken`).
            if let (Some((_, batch)), Some(seed_end)) = (x_board_batch.as_ref(), x_seed_end) {
                if batch.drained() {
                    seed_end.units_taken();
                }
            }
            x_units_stage.end();
            x_lap.lap(XStage::UnitsRelease);
            let x_hits = x_pair_par::XPairPar::<S>::merge_subject(x_units, &x_subject);
            x_lap.lap(XStage::Merge);
            *combined_ungapped_hits = x_hits;
            &[]
        }
        None => {
            x_lap.lap(XStage::Prepare);
            x_units_stage.end();
            x_lap.lap(XStage::UnitsRelease);
            s_frames_preliminary
        }
    };
    // EXPERIMENT (PAIRLUT): with the direct group tables the batch's `lookup` is a
    // stub without entries, so a scan of it by the frame loop of `process_subject` would find no
    // hit. The unit path above must have taken the subject (empty frame slice); anything else
    // stops the run here.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
    // ```c
    //     /* scan the subject sequence for hits */
    //     hits = scansub(lookup_wrap, subject,
    //                               offset_pairs, array_size, scan_range);
    // ```
    if x_pair_lut_is_stub {
        assert!(
            s_frames_preliminary.is_empty() && x_pair_par_stub,
            "PAIRLUT: the reference frame loop would scan the stub lookup table (subject {s_idx})"
        );
    }
    s_frames_preliminary
}

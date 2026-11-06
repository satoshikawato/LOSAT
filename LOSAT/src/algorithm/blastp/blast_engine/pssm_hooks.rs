//! The stages of a psiblast `-in_pssm` search (series PS, item PS4a) that differ from blastp's,
//! called once per stage from `run_resolved_in_pool` when the search has a PSSM
//! (`ResolvedBlastpArgs::psi_search`): the lookup table and contexts, the preliminary search with
//! its cutoffs from the PSSM's Karlin blocks, the CBS-0 traceback, and the Karlin blocks of the
//! report. The kernels and stages are `algorithm::pssm_search`'s; this module hands them blastp's
//! per-run values and collects their hit lists where blastp's stages put theirs. Every function
//! is out of line and cold: blastp's default path only tests `psi_search`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:186-214
//! ```c++
//!     // Set up the score block (some of this code is duplicated in
//!     // blast_setup_cxx.cpp:s_InitScoreBlk, as well as here)
//!     ...
//!     if (pssm.NotEmpty()) {
//!         PsiBlastSetupScoreBlock(score_blk, pssm, messages, options);
//!     }
//!     ...
//!     // Create the lookup table
//!     ...
//!     lookup_table.Reset(LookupTableWrapInit(query_data, ... , score_blk, ...));
//! ```

use anyhow::Result;
use std::cell::RefCell;

use super::{
    blastp_gapped_x_dropoff_raw, database, BlastpRuntimeLookup, BlastpSubjectRef,
    BlastpSubjectResult, BlastpSubjectScratch, PrelimCollector,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::pssm_search::lookup::build_pssm_lookup;
use crate::algorithm::pssm_search::prelim::{PssmPrelim, PssmSubjectScratch};
use crate::algorithm::pssm_search::traceback0::PssmTraceback0;
use crate::algorithm::pssm_search::PsiSearch;
use crate::algorithm::tblastx::constants::{
    GAP_TRIGGER_BIT_SCORE, X_DROP_GAPPED_FINAL, X_DROP_GAPPED_PRELIM, X_DROP_UNGAPPED_BITS,
};
use crate::algorithm::tblastx::lookup::{
    prepare_blosum62_lookup_query_for_word_size, QueryContext,
};
use crate::algorithm::tblastx::ncbi_cutoffs::{gap_trigger_raw_score, x_drop_raw_score};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::stats::search_space::SearchSpace;
use crate::stats::{blast_spouge_etos, BlastGumbelBlk, KarlinParams};

/// The lookup table of the PSSM and the query's contexts (with `kbp_std`, the query's
/// composition with BLOSUM62, and its validity).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:86-95 (quoted at
/// `pssm_search::lookup`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2791
/// ```c
///       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
///       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
/// ```
#[cold]
#[inline(never)]
pub(super) fn pssm_lookup(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let psi = args
        .psi_search
        .as_deref()
        .expect("the lookup table of a PSSM search");
    let (_, _, contexts) =
        prepare_blosum62_lookup_query_for_word_size(query_frames, args.word_size, false);
    let table = build_pssm_lookup(
        &psi.pssm,
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
    );
    (BlastpRuntimeLookup::Aa(table), contexts)
}

/// The Karlin blocks of a PSSM search: `kbp_psi[0]` and `kbp_gap_psi[0]`.
#[derive(Clone, Copy)]
pub(super) struct PsiKarlin {
    pub(super) ungapped: KarlinParams,
    pub(super) gapped: KarlinParams,
}

/// `kbp_psi[0]` from the query's `kbp_std[0]` and `kbp_gap_psi[0]` from `kbp_gap_std[0]` (the
/// table's BLOSUM62 values, `gapped_params`), as `PsiBlastSetupScoreBlock` sets them.
pub(super) fn psi_karlin(
    psi: &PsiSearch,
    contexts: &[QueryContext],
    gapped_params: &KarlinParams,
) -> PsiKarlin {
    PsiKarlin {
        ungapped: psi.ungapped_karlin(&contexts[0].karlin_params),
        gapped: psi.gapped_karlin(gapped_params),
    }
}

thread_local! {
    /// The PSSM preliminary search's memory of this thread (`PssmSubjectScratch`), which the
    /// subject loops' scratch (`BlastpSubjectScratch`, blastp's) does not hold.
    static PSSM_SCRATCH: RefCell<PssmSubjectScratch> = RefCell::new(PssmSubjectScratch::new());
}

/// The per-run values of the preliminary search that `run_resolved_in_pool` computed for blastp.
pub(super) struct PssmPrelimRun<'a> {
    pub(super) contexts: &'a [QueryContext],
    pub(super) search_spaces: &'a [SearchSpace],
    pub(super) exact_search_spaces: &'a [Option<i64>],
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gumbel: Option<&'a BlastGumbelBlk>,
    pub(super) cutoff_subject_length: i32,
    pub(super) prelim_evalue: f64,
    pub(super) window: i32,
    pub(super) diag_mask: i32,
    pub(super) offset_array_size: i32,
    pub(super) diag_array_size: i32,
    pub(super) total_query_span: i32,
    pub(super) prelim_hitlist_size: usize,
}

/// The preliminary search of a PSSM query over the database `database` or the subject set
/// `subjects`, its HSP lists collected into `x_compact_hit_lists` as blastp's are.
///
/// The cutoffs are NCBI's with the PSSM's blocks: the ungapped X-drop from `kbp_psi`
/// (`sbp->kbp`), the gap trigger from `kbp_std`, the hit-saving cutoff and the gapped X-drops
/// from `kbp_gap_psi` (`sbp->kbp_gap`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:204-221
/// ```c
///       kbp= sbp->kbp[context];
///       ...
///         p->cutoffs[context].x_dropoff_init =
///             (Int4)(sbp->scale_factor *
///                    ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:338-345
/// ```c
///       if (sbp->kbp_std) {     /* this may not be set for gapped blastn */
///          kbp = sbp->kbp_std[context];
///          if (s_BlastKarlinBlkIsValid(kbp)) {
///             gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 +
///                                      kbp->logK) / kbp->Lambda);
///          }
///       }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:914-943
/// ```c
///          kbp = kbp_array[context];
///          ASSERT(s_BlastKarlinBlkIsValid(kbp));
///          searchsp = query_info->contexts[context].eff_searchsp;
///          ...
///          if (sbp->gbp && sbp->gbp->filled) {
///              /* If cbs greater than 1 (2 or 3), then increase expect value by 5 for preliminary search. */
///              int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
///              params->prelim_evalue = cbs_stretch*evalue;
///              new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
///                          query_info->contexts[context].query_length,
///                          avg_subject_length);
///          } else {
///              BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
///          }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-463
/// ```c
///    if (sbp->kbp_gap) {
///       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
///       params->gap_x_dropoff = (Int4)
///           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:152-169
/// ```c
///    if (ewp->diag_table->offset >= INT4_MAX / 4) {
///        ewp->diag_table->offset = ewp->diag_table->window;
///        s_BlastDiagClear(ewp->diag_table);
///    } else {
///        ewp->diag_table->offset += subject_length + ewp->diag_table->window;
///    }
/// ```
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn pssm_prelim(
    psi: &PsiSearch,
    args: &ResolvedBlastpArgs,
    lookup: &BlastpRuntimeLookup,
    run: &PssmPrelimRun<'_>,
    database: Option<&database::BlastpDatabase>,
    subjects: &[EncodedProtein],
    pool: &crate::utils::threading::SearchPool<'_>,
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> Result<()> {
    let BlastpRuntimeLookup::Aa(table) = lookup else {
        unreachable!("the lookup table of a PSSM search is `pssm_lookup`'s");
    };
    let ctx = &run.contexts[0];
    // A query whose Karlin-Altschul parameters cannot be computed is not searched (psiblast
    // stops at it before, `psiblast_stop_at_invalid_query`).
    if !ctx.is_valid {
        return Ok(());
    }
    let karlin = psi_karlin(psi, run.contexts, run.gapped_params);
    let x_dropoff = x_drop_raw_score(X_DROP_UNGAPPED_BITS, &karlin.ungapped, 1.0);
    let hit_cutoff = match run.gumbel {
        Some(gumbel) => blast_spouge_etos(
            run.prelim_evalue,
            &karlin.gapped,
            gumbel,
            ctx.aa_len as i32,
            run.cutoff_subject_length,
        ),
        None => crate::algorithm::tblastx::ncbi_cutoffs::cutoff_score_from_evalue(
            args.evalue,
            run.exact_search_spaces[0].unwrap_or(0),
            &karlin.gapped,
        )
        .max(1),
    };
    let word_cutoff =
        gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, &ctx.karlin_params).min(hit_cutoff);
    let window = run.window;
    let prelim = PssmPrelim {
        rows: &psi.pssm,
        lookup: table,
        word_size: args.word_size as i32,
        window,
        x_dropoff,
        word_cutoff,
        hit_cutoff,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        x_drop_gapped: blastp_gapped_x_dropoff_raw(
            X_DROP_GAPPED_PRELIM as f64,
            karlin.gapped.lambda,
        ),
        karlin: &karlin.gapped,
        gumbel: run.gumbel,
        search_space: run.search_spaces[0].effective_space,
        prelim_evalue: run.prelim_evalue,
    };
    let diag_mask = run.diag_mask;
    let process_subject = |scratch: &mut BlastpSubjectScratch,
                           s_idx: usize,
                           subject: BlastpSubjectRef<'_>,
                           subject_result: &mut BlastpSubjectResult| {
        subject_result.reset(s_idx);
        // A subject without letters is not searched, and the diagonal offset does not move
        // for it (blast_engine.c:1429-1432).
        if subject.aa_len == 0 {
            return;
        }
        let list = PSSM_SCRATCH.with(|own| {
            prelim.search_subject(
                &mut own.borrow_mut(),
                &mut scratch.offset_pairs,
                &mut scratch.diag_array,
                diag_mask,
                scratch.diag_offset,
                s_idx as u32,
                subject.aa_seq,
                subject.aa_len,
            )
        });
        if scratch.diag_offset >= i32::MAX / 4 {
            scratch.diag_offset = window;
            for diag in &mut scratch.diag_array {
                *diag = crate::algorithm::tblastx::blast_extend::DiagStruct::clear(window);
            }
        } else {
            scratch.diag_offset += subject.aa_len as i32 + window;
        }
        if let Some(list) = list {
            subject_result.x_compact_hsp_lists.push((0, list));
        }
    };
    let new_scratch = || {
        BlastpSubjectScratch::new(
            run.offset_array_size,
            run.diag_array_size,
            window,
            run.total_query_span,
        )
    };
    match database {
        // The preliminary stage of a database: blocks of sequences in the order of their ordinal
        // ids (`database::search_blocks`, the deterministic merge of every thread count).
        Some(database) => database::search_blocks(
            database,
            pool,
            window,
            &new_scratch,
            &process_subject,
            &mut PrelimCollector {
                x_compact_hit_lists,
                prelim_hitlist_size: run.prelim_hitlist_size,
            },
        ),
        // A subject set: each subject in turn (`blast_engine.c:1409-1498`), with one diagonal
        // table, at every thread count.
        None => {
            let mut scratch = new_scratch();
            let mut subject_result = BlastpSubjectResult::new(0);
            for (s_idx, subject) in subjects.iter().enumerate() {
                process_subject(
                    &mut scratch,
                    s_idx,
                    BlastpSubjectRef::of(subject),
                    &mut subject_result,
                );
                for (q_idx, list) in subject_result.x_compact_hsp_lists.drain(..) {
                    x_compact_hit_lists[q_idx]
                        .get_or_insert_with(|| BlastpHitList::new(run.prelim_hitlist_size))
                        .update(list);
                }
            }
            Ok(())
        }
    }
}

/// The CBS-0 traceback of a PSSM search (`pssm_search::traceback0`): the final hit lists of the
/// preliminary lists, which it takes. The final gapped X-drop and the e-values read
/// `kbp_gap_psi`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:464-467
/// ```c
///       params->gap_x_dropoff_final = (Int4)
///           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
/// ```
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn pssm_traceback(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    total_query_span: i32,
    subjects: XRedoSubjects<'_>,
    gapped_params: &KarlinParams,
    gumbel: Option<&BlastGumbelBlk>,
    search_spaces: &[SearchSpace],
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<Vec<Option<BlastpHitList>>> {
    let psi = args
        .psi_search
        .as_deref()
        .expect("the traceback of a PSSM search");
    let karlin = psi_karlin(psi, contexts, gapped_params);
    let x_drop_gapped =
        blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_PRELIM as f64, karlin.gapped.lambda);
    let x_drop_final =
        blastp_gapped_x_dropoff_raw(X_DROP_GAPPED_FINAL as f64, karlin.gapped.lambda)
            .max(x_drop_gapped);
    PssmTraceback0 {
        rows: &psi.pssm,
        context: &contexts[0],
        total_query_span,
        subjects,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        x_drop_final,
        karlin: &karlin.gapped,
        gumbel,
        search_space: search_spaces[0].effective_space,
        expect_value: args.evalue,
        hitlist_size: args.traceback_hitlist_size(),
    }
    .run(preliminary_hit_lists, compact_hit_lists, pool)
}

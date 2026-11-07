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
use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::pssm_search::lookup::build_pssm_lookup;
use crate::algorithm::pssm_search::prelim::{PssmPrelim, PssmSubjectScratch};
use crate::algorithm::pssm_search::redo1;
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
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::composition_adjustment::compo_scale::PosMatrixInfo;
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

/// The Gumbel block of a PSSM search (`sbp->gbp`), which the e-values and the hit-saving cutoff
/// read (Spouge's finite-size correction). blastp's search has none only with `-searchsp`, which
/// LOSAT's psiblast does not take (its parser refuses it, and `PsiblastArgs` sets
/// `BlastpArgs::searchsp` to `None`), so a PSSM search always has one.
fn psi_gumbel(gumbel: Option<&BlastGumbelBlk>) -> &BlastGumbelBlk {
    gumbel.expect("a psiblast search has a Gumbel block (LOSAT's psiblast takes no -searchsp)")
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
    pub(super) gapped_params: &'a KarlinParams,
    /// `sbp->gbp` (`psi_gumbel`).
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
    let gumbel = psi_gumbel(run.gumbel);
    let hit_cutoff = blast_spouge_etos(
        run.prelim_evalue,
        &karlin.gapped,
        gumbel,
        ctx.aa_len as i32,
        run.cutoff_subject_length,
    );
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
        gumbel,
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

/// The traceback of a PSSM search: mode 0's (`pssm_search::traceback0`), or the position-based
/// redo of modes 1 to 3 (`pssm_redo`). The final hit lists of the preliminary lists, which it
/// takes. The final gapped X-drop and the e-values read `kbp_gap_psi`.
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
    if args.comp_based_stats.mode != BlastpCompositionMode::NoCompositionBasedStats {
        return pssm_redo(
            args,
            psi,
            &contexts[0],
            subjects,
            gapped_params,
            &karlin,
            x_drop_final,
            psi_gumbel(gumbel),
            preliminary_hit_lists,
            compact_hit_lists,
        );
    }
    PssmTraceback0 {
        rows: &psi.pssm,
        context: &contexts[0],
        total_query_span,
        subjects,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        x_drop_final,
        karlin: &karlin.gapped,
        gumbel: psi_gumbel(gumbel),
        expect_value: args.evalue,
        hitlist_size: args.traceback_hitlist_size(),
    }
    .run(preliminary_hit_lists, compact_hit_lists, pool)
}

/// The position-based redo of a PSSM search (`pssm_search::redo1`, modes 1 to 3): the matrix
/// information of the PSSM (`PosMatrixInfo`, built when a preliminary list is there to redo),
/// then the redo of every list.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1497
/// ```c
///     if (ext_params->options->compositionBasedStats > 0 ||
///         ext_params->options->eTbackExt == eSmithWatermanTbck) {
///         status =
///             Blast_RedoAlignmentCore_MT(program_number, ...);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3254-3258
/// ```c
///         sbp_tld[i] = s_BlastScoreBlk_Copy(
///                 program_number,
///                 sbp,
///                 sbp->alphabet_code,
///                 sbp->number_of_contexts
///         );
/// ```
/// The copy's matrix information reads its frequency ratios (zero where the PSSM has none: then
/// every column takes BLOSUM62's), its ideal block and `kbp_gap_std[0]`'s K.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn pssm_redo(
    args: &ResolvedBlastpArgs,
    psi: &PsiSearch,
    context: &QueryContext,
    subjects: XRedoSubjects<'_>,
    gapped_params: &KarlinParams,
    karlin: &PsiKarlin,
    x_drop_final: i32,
    gumbel: &BlastGumbelBlk,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> Result<Vec<Option<BlastpHitList>>> {
    let has_lists = preliminary_hit_lists
        .iter()
        .zip(compact_hit_lists.iter())
        .any(|(full, compact)| {
            full.as_ref().is_some_and(|list| list.hsplist_count > 0)
                || compact.as_ref().is_some_and(|list| list.hsplist_count > 0)
        });
    if !has_lists {
        return Ok(preliminary_hit_lists.iter().map(|_| None).collect());
    }
    let sbp = crate::algorithm::psiblast::pssm::score_block(
        args.scoring.matrix,
        crate::core::psi::KarlinBlk {
            lambda: gapped_params.lambda,
            k: gapped_params.k,
            h: gapped_params.h,
        },
    )?;
    let query = &context.aa_seq[1..context.aa_seq.len() - 1];
    let zero_ratios;
    let freq_ratios: &[[f64; crate::core::blast_util::BLASTAA_SIZE]] = if psi.freq_ratios.is_empty()
    {
        zero_ratios = vec![[0.0; crate::core::blast_util::BLASTAA_SIZE]; query.len()];
        &zero_ratios
    } else {
        &psi.freq_ratios
    };
    let kbp_psi = crate::core::psi::KarlinBlk {
        lambda: karlin.ungapped.lambda,
        k: karlin.ungapped.k,
        h: karlin.ungapped.h,
    };
    let Some(info) = PosMatrixInfo::new(query, freq_ratios, &sbp, redo1::SCALING_FACTOR, kbp_psi)
    else {
        anyhow::bail!(
            "a PSSM whose position-based scaling for the composition-based statistics gives no Karlin-Altschul parameters (NCBI BLAST+ 2.17.0 goes on with a Lambda of -1) is not supported by LOSAT's PSIBLAST"
        );
    };
    redo1::PssmRedo1 {
        pssm: &psi.pssm,
        info: &info,
        context,
        subjects,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        karlin: &karlin.gapped,
        x_drop_final,
        gumbel,
        expect_value: args.evalue,
        hitlist_size: args.traceback_hitlist_size(),
        inclusion_ethresh: args.inclusion_ethresh,
    }
    .run(preliminary_hit_lists, compact_hit_lists)
}

/// The warning line of a PSSM query (psiblast `-in_pssm`): its label is the PSSM's Seq-id (the
/// record's local id, DECISIONS 12 of line PS), and its messages add those of
/// `PsiBlastSetupScoreBlock` (`PsiSearch::setup_messages`), which NCBI keeps with the query's
/// others (sorted, `report::query_warnings::sorted_query_messages`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:188-190
/// ```c++
///         CRef<CSearchMessage> sm(new CSearchMessage(eBlastSevWarning, 0,
///                                                    os.str()));
///         _ASSERT(messages.size() == 1); // PSI-BLAST only works with one query
///         messages.front().push_back(sm);
/// ```
#[cold]
#[inline(never)]
pub(super) fn pssm_query_warning(
    program: &str,
    psi: &PsiSearch,
    record: &FastaRecord,
    messages: &[String],
) -> Vec<u8> {
    let mut all = messages.to_vec();
    all.extend(psi.setup_messages.iter().cloned());
    crate::report::query_warnings::query_warning_with_seq_id(
        program,
        record.local_id.as_bytes(),
        record,
        &all,
    )
}

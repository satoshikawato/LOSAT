//! The stages of a psiblast `-in_pssm` search (series PS, item PS4a) that differ from blastp's,
//! called once per stage from `run_call` when the search has a PSSM
//! (`ResolvedBlastpArgs::psi_search`): the lookup table and contexts, the preliminary search with
//! its cutoffs from the PSSM's Karlin blocks, the CBS-0 traceback, and the Karlin blocks of the
//! report. The kernels and stages are `algorithm::pssm_search`'s; this module hands them blastp's
//! per-run values and collects their hit lists where blastp's stages put theirs. Every function
//! is out of line and cold: blastp's default path only tests `psi_search`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:186-214

use anyhow::Result;
use std::cell::RefCell;

use super::{
    database, BlastpRuntimeLookup, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch,
    PrelimCollector, PrelimFilter,
};
use crate::algorithm::blastp::args::{BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::lookup::build_pssm_lookup_in;
use crate::algorithm::pssm_search::prelim::{PssmPrelim, PssmSubjectScratch};
use crate::algorithm::pssm_search::redo1;
use crate::algorithm::pssm_search::traceback0::PssmTraceback0;
use crate::algorithm::pssm_search::PsiSearch;
use crate::algorithm::tblastx::lookup::{prepare_matrix_lookup_query_for_word_size, QueryContext};
use crate::algorithm::tblastx::ncbi_cutoffs::{gap_trigger_raw_score, x_drop_raw_score};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::composition_adjustment::compo_scale::PosMatrixInfo;
use crate::core::link_hsps::{EvenGapLinkParams, LinkContext};
use crate::stats::search_space::SearchSpace;
use crate::stats::{BlastGumbelBlk, KarlinParams};

/// The lookup table of the PSSM and the query's contexts ([`psi_contexts`]): the PSSM's words
/// over the query's lookup segments, the complement of its masks (`-seg`, `-lcase_masking`; with
/// `-soft_masking true` the masks keep words out of the table only), which blastp's query
/// preparation computes for its own table (`prepare_matrix_lookup_query_for_word_size`). The
/// query frames are blastp's (`blastp_query_frames`): a hard mask's letters are X there, which
/// `kbp_std` reads (`psi_contexts`), and the PSSM search reads them where NCBI reads the masked
/// query (the gap trigger, the position-based scaling of `pssm_redo`). The PSSM itself is made
/// of the query's letters unmasked (`psi_round::next_pssm`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:170-214
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-637 (the masked letters,
/// then the complement of the masks)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:86-95 (cited at
/// `pssm_search::lookup`)
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
    let (lookup_segments, mut contexts) = psi_query(args, query_frames);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215 (cited above:
    // no lookup table for a split query)
    // A later round of a query that NCBI splits: its chunks search without the PSSM
    // (`search_split_query_batches`), so the PSSM's table indexes no word, and the round reads
    // `kbp_std` of its set-up's query (`PsiSearch::split_setup`), whose validity the round driver
    // checked.
    let lookup_segments: Vec<(i64, i64)> = match psi.split_setup.as_ref() {
        Some(setup) => {
            let (_, setup_contexts) = psi_query(args, &[vec![setup.clone()]]);
            for (context, setup) in contexts.iter_mut().zip(setup_contexts) {
                context.karlin_params = setup.karlin_params;
                context.is_valid = setup.is_valid;
            }
            Vec::new()
        }
        None => lookup_segments
            .iter()
            .map(|&(left, right)| (i64::from(left), i64::from(right)))
            .collect(),
    };
    let table = build_pssm_lookup_in(
        &psi.pssm,
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
        &lookup_segments,
    );
    (BlastpRuntimeLookup::Aa(table), contexts)
}

/// The query's contexts of a PSSM search (and of psiblast's check of its queries before the
/// search, `psiblast_stop_at_invalid_query`): `kbp_std`, the query's composition with the matrix
/// of the search (`-matrix`: the score block's matrix, which the PSSM's scores replace only after
/// its Karlin blocks), and its validity.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:498-515
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2791
pub(super) fn psi_contexts(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> Vec<QueryContext> {
    psi_query(args, query_frames).1
}

/// The query's lookup segments (inclusive offsets, NCBI's `lookup_segments` of the query's one
/// context) and its contexts ([`psi_contexts`]), as blastp's query preparation computes them.
fn psi_query(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> (Vec<(i32, i32)>, Vec<QueryContext>) {
    let (_, lookup_segments, contexts) = prepare_matrix_lookup_query_for_word_size(
        query_frames,
        args.word_size,
        false,
        args.scoring.matrix,
        None,
    );
    (lookup_segments, contexts)
}

/// The Karlin blocks of a PSSM search: `kbp_psi[0]` and `kbp_gap_psi[0]`.
#[derive(Clone, Copy)]
pub(super) struct PsiKarlin {
    pub(super) ungapped: KarlinParams,
    pub(super) gapped: KarlinParams,
}

/// What the e-values and the hit-saving cutoff of a PSSM search read (`EvalueStats`): the Gumbel
/// block (`sbp->gbp`, Spouge's finite-size correction), and without it (`-searchsp`, `OLD_FSC`:
/// `ResolvedBlastpArgs::old_fsc`) the effective search space of the query's one context,
/// `exact_search_spaces[0]`, which `call_setup::call_search_spaces` computes for blastp.
///
/// Called for a valid context only (`pssm_prelim`, `pssm_traceback`): the space is then exact
/// (`exact_search_spaces` holds every context's), and its one `None`, a subject set without
/// letters, has been NCBI's engine error before the search ("The average subject length is too
/// short").
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1885-1890
fn psi_stats<'a>(
    gumbel: Option<&'a BlastGumbelBlk>,
    exact_search_spaces: &[Option<i64>],
) -> EvalueStats<'a> {
    match gumbel {
        // With the block nothing of a psiblast search reads the search space (`e_value` and
        // `hit_cutoff` do not; `-use_sw_tback`, whose Smith-Waterman redo would, is refused).
        // A reader that needs it with the block (series PR, item PS11c-3) fills it first.
        Some(gumbel) => EvalueStats::spouge(gumbel),
        None => EvalueStats::new(
            None,
            exact_search_spaces[0]
                .expect("a valid context's search space without the Gumbel block is exact"),
        ),
    }
}

/// `kbp_psi[0]` from the query's `kbp_std[0]` and `kbp_gap_psi[0]` from `kbp_gap_std[0]` (the
/// table's values of the matrix and gap costs, `gapped_params`), as `PsiBlastSetupScoreBlock`
/// sets them.
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

/// The per-run values of the preliminary search that `run_call` computed for blastp.
pub(super) struct PssmPrelimRun<'a> {
    pub(super) contexts: &'a [QueryContext],
    pub(super) gapped_params: &'a KarlinParams,
    /// `sbp->gbp` (`psi_stats`).
    pub(super) gumbel: Option<&'a BlastGumbelBlk>,
    /// Each context's effective search space where it is exact (`psi_stats`).
    pub(super) exact_search_spaces: &'a [Option<i64>],
    /// Each context's length adjustment and effective search space (`query_info->contexts`,
    /// computed with `kbp_gap_std`), which psiblast's `-sum_stats true` links with
    /// (`psi_link_context`).
    pub(super) search_spaces: &'a [SearchSpace],
    pub(super) cutoff_subject_length: i32,
    pub(super) prelim_evalue: f64,
    /// The statistics of each subject (`DbStatistics::per_subject`), `None` for the setup's.
    pub(super) per_subject: Option<PsiSubjectStatistics>,
    pub(super) window: i32,
    pub(super) diag_mask: i32,
    pub(super) offset_array_size: i32,
    pub(super) diag_array_size: i32,
    pub(super) total_query_span: i32,
    pub(super) prelim_hitlist_size: usize,
}

/// What the statistics of one subject read in a PSSM search (`BLAST_OneSubjectUpdateParameters`,
/// `stats::setup::subject_eff_lengths`): the query's preliminary option space, the setup's number
/// of sequences and `-dbsize`.
#[derive(Clone, Copy)]
pub(super) struct PsiSubjectStatistics {
    pub(super) option_space: i64,
    pub(super) setup_count: i32,
    pub(super) dbsize: i64,
}

/// The preliminary and final gapped X-drops of a PSSM search (raw scores): the search's X-drops in
/// bits at the smallest gapped Lambda of the query batch, the PSSM's `kbp_gap_psi` (one query).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:454-463 (cited at
/// `core::x_drop::gapped_x_dropoffs`)
fn psi_gapped_x_dropoffs(args: &ResolvedBlastpArgs, karlin: &PsiKarlin) -> (i32, i32) {
    crate::core::x_drop::gapped_x_dropoffs(
        args.x_drop_gapped_bits,
        args.x_drop_final_bits,
        crate::core::x_drop::extension_lambda(&karlin.gapped),
    )
}

/// The hit-saving cutoff of a PSSM search's context (`hit_params->cutoffs[0]`) at the subject
/// length `cutoff_subject_length` and the e-value `prelim_evalue`: `cutoff_score_max` with
/// `kbp_gap_psi`, which the initial words' cutoff reads, and `cutoff_score`, which psiblast's
/// `-sum_stats true` lowers to the sum statistics' cutoff at that length where that is smaller
/// (`blast_engine::sum_stats_hit_cutoffs` with `kbp_gap_psi`: NCBI's `kbp_array` is
/// `sbp->kbp_gap`); without sum statistics the two are the same. `psi_prelim_cutoffs` calls it
/// with the setup's subject length and preliminary e-value, or with a subject's length and
/// `-evalue` unstretched (the statistics of each subject). With the setup's (as `pssm_traceback`
/// calls it), its `cutoff_score` is the traceback's `cutoff_score_min` (one context): the
/// traceback's hit-saving parameters are made again from the same options
/// (`BLAST_GapAlignSetUp`), and the update for each subject does not reach them
/// (blast_traceback.c:1596-1628 updates the thread's copy).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:949-976 (cited at
/// `blast_engine::sum_stats_hit_cutoffs`)
fn psi_hit_cutoffs(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    karlin: &PsiKarlin,
    stats: &EvalueStats<'_>,
    cutoff_subject_length: i32,
    prelim_evalue: f64,
) -> (i32, i32) {
    let cutoff_score_max = stats.hit_cutoff(
        prelim_evalue,
        &karlin.gapped,
        contexts[0].aa_len as i32,
        cutoff_subject_length,
    );
    let cutoff_score = if args.sum_stats {
        // The one context's (a PSSM search has one).
        let sum_cutoff =
            super::sum_stats_hit_cutoffs(contexts, cutoff_subject_length, &karlin.gapped)[0];
        cutoff_score_max.min(sum_cutoff)
    } else {
        cutoff_score_max
    };
    (cutoff_score_max, cutoff_score)
}

/// The query's context as psiblast's `-sum_stats true` links it in a PSSM search: `kbp_gap_psi`
/// (`sbp->kbp_gap`) and the context's length adjustment and effective search space, which NCBI
/// computed with `kbp_gap_std` (blast_setup.c:766-768) as for the matrix's rounds.
fn psi_link_context(
    karlin: &PsiKarlin,
    context: &QueryContext,
    search_space: &SearchSpace,
) -> LinkContext {
    super::blastp_link_context(&karlin.gapped, context.aa_len, search_space)
}

/// The cutoffs of a PSSM search's preliminary stage at one subject length: the setup's
/// (`cutoff_subject_length`, the preliminary e-value, the setup's search space) or, with the
/// statistics of each subject, the subject's own (`BLAST_OneSubjectUpdateParameters`: its
/// length, `-evalue` unstretched, its search space). The initial words' cutoff is the gap trigger
/// (`kbp_std`) at most `cutoff_score_max` (`BlastInitialWordParametersUpdate`), the hit-saving
/// cutoff is `psi_hit_cutoffs`' `cutoff_score`, and with psiblast's `-sum_stats true` the
/// linking reads the words' cutoff as its small-gap cutoff (`BlastLinkHSPParametersUpdate`, one
/// context: `EvenGapLinkParams::preliminary`) and the context at `search_space`. The matrix's
/// counterpart is `blast_engine::blastp_prelim_cutoffs`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:625-646
#[allow(clippy::too_many_arguments)]
fn psi_prelim_cutoffs(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    karlin: &PsiKarlin,
    stats: &EvalueStats<'_>,
    subject_length: i32,
    evalue: f64,
    search_space: &SearchSpace,
) -> (i32, i32, Option<(EvenGapLinkParams, LinkContext)>) {
    let ctx = &contexts[0];
    let (cutoff_score_max, hit_cutoff) =
        psi_hit_cutoffs(args, contexts, karlin, stats, subject_length, evalue);
    let word_cutoff =
        gap_trigger_raw_score(args.gap_trigger_bits, &ctx.karlin_params).min(cutoff_score_max);
    let link = args.sum_stats.then(|| {
        (
            EvenGapLinkParams::preliminary(word_cutoff),
            psi_link_context(karlin, ctx, search_space),
        )
    });
    (word_cutoff, hit_cutoff, link)
}

/// The preliminary search of a PSSM query over the database `database` or the subject set
/// `subjects`, its HSP lists collected into `x_compact_hit_lists` as blastp's are.
///
/// The cutoffs are NCBI's with the PSSM's blocks: the ungapped X-drop from `kbp_psi`
/// (`sbp->kbp`), the gap trigger from `kbp_std`, the hit-saving cutoff and the gapped X-drops
/// from `kbp_gap_psi` (`sbp->kbp_gap`). The X-drops and the gap trigger in bits are the search's
/// (`-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`, `-gap_trigger` or the defaults,
/// `ResolvedBlastpArgs`), converted as blastp's are (the X-drops in `core::x_drop`, the gap trigger
/// in `tblastx::ncbi_cutoffs::gap_trigger_raw_score`); a gapped X-drop that NCBI's
/// kernels cannot size with the PSSM's Lambda is refused before the search, as blastp refuses it
/// with the matrix's (`core::x_drop::check_gapped_x_dropoffs`). With `-culling_limit` or
/// `-best_hit_*` the subjects' lists go to the option's writer `prelim_filter` in place of the
/// collector, as blastp's do.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:204-221
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:338-345
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:914-943
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-463
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:152-172
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
    prelim_filter: Option<&mut PrelimFilter>,
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
    let x_dropoff = x_drop_raw_score(args.x_drop_ungapped_bits, &karlin.ungapped, 1.0);
    let x_drops = psi_gapped_x_dropoffs(args, &karlin);
    crate::core::x_drop::check_gapped_x_dropoffs(
        args.app.name(),
        (args.x_drop_gapped_bits, args.x_drop_final_bits),
        x_drops,
        (args.scoring.gap_open, args.scoring.gap_extend),
        args.comp_based_stats.mode == BlastpCompositionMode::NoCompositionBasedStats,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-236
    // A later round of a query that NCBI splits: only its chunks are searched
    // (`search_split_query_batches`); the whole query scans no subject (its table is empty,
    // `pssm_lookup`). The X-drops above are the PSSM's, which NCBI's traceback converts.
    if psi.split_setup.is_some() {
        return Ok(());
    }
    let stats = psi_stats(run.gumbel, run.exact_search_spaces);
    // `run.prelim_evalue` is `-evalue` without the Gumbel block (`call_setup::prelim_cutoffs`,
    // blast_parameters.c:931-939).
    // psiblast's `-sum_stats true`: the preliminary stage links each subject's HSPs
    // (`psi_prelim_cutoffs`).
    let (word_cutoff, hit_cutoff, link) = psi_prelim_cutoffs(
        args,
        run.contexts,
        &karlin,
        &stats,
        run.cutoff_subject_length,
        run.prelim_evalue,
        &run.search_spaces[0],
    );
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
        x_drop_gapped: x_drops.0,
        karlin: &karlin.gapped,
        stats,
        prelim_evalue: run.prelim_evalue,
        link,
        subject_besthit: args.subject_besthit,
    };
    let diag_mask = run.diag_mask;
    let subject_space_overflow = crate::stats::setup::SubjectSpaceOverflow::default();
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443 and
        // c++/src/algo/blast/core/blast_setup.c:1015-1025 (cited at `call_setup`'s
        // `subject_parameters`): the statistics of each subject, with the subject's length in the
        // hit-saving cutoff (no stretch of the e-value; with psiblast's `-sum_stats true` lowered
        // at the subject's length) and in the effective lengths, which the linking reads with the
        // subject's small-gap cutoff (`psi_prelim_cutoffs`).
        let subject_prelim;
        let prelim = match run.per_subject {
            None => &prelim,
            Some(statistics) => {
                let lengths = match crate::stats::setup::subject_eff_lengths(
                    ctx.aa_len as i64,
                    subject.aa_len as i64,
                    statistics.setup_count,
                    statistics.dbsize,
                    statistics.option_space,
                    false,
                    run.gapped_params,
                    args.app.name(),
                ) {
                    Ok(lengths) => lengths,
                    Err(overflow) => {
                        subject_space_overflow.note(overflow);
                        return;
                    }
                };
                let stats = match run.gumbel {
                    Some(gumbel) => EvalueStats::spouge(gumbel),
                    None => EvalueStats::new(None, lengths.eff_searchsp),
                };
                let (word_cutoff, hit_cutoff, link) = psi_prelim_cutoffs(
                    args,
                    run.contexts,
                    &karlin,
                    &stats,
                    subject.aa_len as i32,
                    args.evalue,
                    &SearchSpace {
                        effective_space: lengths.eff_searchsp as f64,
                        length_adjustment: i64::from(lengths.length_adjustment),
                    },
                );
                subject_prelim = PssmPrelim {
                    word_cutoff,
                    hit_cutoff,
                    stats,
                    prelim_evalue: args.evalue,
                    link,
                    ..prelim
                };
                &subject_prelim
            }
        };
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
    let searched = match (database, prelim_filter) {
        // The preliminary stage of a database: blocks of sequences in the order of their ordinal
        // ids (`database::search_blocks`, the deterministic merge of every thread count).
        (Some(database), None) => database::search_blocks(
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
        (Some(database), Some(filter)) => database::search_blocks(
            database,
            pool,
            window,
            &new_scratch,
            &process_subject,
            filter,
        ),
        // A subject set: each subject in turn (`blast_engine.c:1409-1498`), with one diagonal
        // table, at every thread count.
        (None, mut prelim_filter) => {
            let mut scratch = new_scratch();
            let mut subject_result = BlastpSubjectResult::new(0);
            for (s_idx, subject) in subjects.iter().enumerate() {
                process_subject(
                    &mut scratch,
                    s_idx,
                    BlastpSubjectRef::of(subject),
                    &mut subject_result,
                );
                if let Some(filter) = prelim_filter.as_deref_mut() {
                    filter.add(&mut subject_result);
                    continue;
                }
                for (q_idx, list) in subject_result.x_compact_hsp_lists.drain(..) {
                    x_compact_hit_lists[q_idx]
                        .get_or_insert_with(|| BlastpHitList::new(run.prelim_hitlist_size))
                        .update(list);
                }
            }
            Ok(())
        }
    };
    searched?;
    subject_space_overflow.check()?;
    Ok(())
}

/// The traceback of a PSSM search: mode 0's (`pssm_search::traceback0`), or the position-based
/// redo of modes 1 to 3 (`pssm_redo`). The final hit lists of the preliminary lists, which it
/// takes. The final gapped X-drop and the e-values read `kbp_gap_psi`. With psiblast's
/// `-sum_stats true` both link each list's HSPs with large gaps only
/// (`EvenGapLinkParams::traceback`), and the redo keeps the alignments that score at least the
/// traceback's hit-saving cutoff (`psi_hit_cutoffs`, `redo1::RedoLink`), whose `-evalue`
/// and subject length (`cutoff_subject_length`) are the preliminary stage's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:462-467
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
    exact_search_spaces: &[Option<i64>],
    search_spaces: &[SearchSpace],
    (cutoff_subject_length, prelim_evalue): (i32, f64),
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    pool: &crate::utils::threading::SearchPool<'_>,
    redone: &dyn Fn(&[crate::algorithm::blastp::hsp::BlastCompoHeap]) -> anyhow::Result<()>,
    sorts_results: bool,
) -> Result<Vec<Option<BlastpHitList>>> {
    let psi = args
        .psi_search
        .as_deref()
        .expect("the traceback of a PSSM search");
    // A query whose Karlin-Altschul parameters cannot be computed has no preliminary list
    // (`pssm_prelim`).
    if !contexts[0].is_valid {
        return Ok(preliminary_hit_lists.iter().map(|_| None).collect());
    }
    let stats = psi_stats(gumbel, exact_search_spaces);
    let karlin = psi_karlin(psi, contexts, gapped_params);
    let (_, x_drop_final) = psi_gapped_x_dropoffs(args, &karlin);
    let link_context = args
        .sum_stats
        .then(|| psi_link_context(&karlin, &contexts[0], &search_spaces[0]));
    if args.comp_based_stats.mode != BlastpCompositionMode::NoCompositionBasedStats {
        let link = link_context.map(|context| {
            let (_, cutoff_score_min) = psi_hit_cutoffs(
                args,
                contexts,
                &karlin,
                &stats,
                cutoff_subject_length,
                prelim_evalue,
            );
            redo1::RedoLink::new(cutoff_score_min, context)
        });
        return pssm_redo(
            args,
            psi,
            &contexts[0],
            subjects,
            gapped_params,
            &karlin,
            x_drop_final,
            stats,
            link,
            preliminary_hit_lists,
            compact_hit_lists,
            redone,
        )
        .map(seqalign_scores_of_lists);
    }
    PssmTraceback0 {
        rows: &psi.pssm,
        context: &contexts[0],
        total_query_span,
        subjects,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        matrix: args.scoring.matrix,
        x_drop_final,
        karlin: &karlin.gapped,
        stats,
        expect_value: args.evalue,
        hitlist_size: args.traceback_hitlist_size(),
        sorts_results,
        link: link_context.map(|context| (EvenGapLinkParams::traceback(), context)),
    }
    .run(preliminary_hit_lists, compact_hit_lists, pool)
    .map(seqalign_scores_of_lists)
}

/// The final lists' E-values and bit scores as the Seq-align that NCBI makes of them carries them
/// (`seqalign_scores`): a NaN E-value stays as the one it lacks (`seqalign_evalue_absent`), and a
/// bit score below 0 (a PSSM whose gapped Karlin-Altschul block gives them) is the readers' -1.
/// The lists are in their final
/// order (the traceback's sorts read the E-values as the search made them).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1185-1197
fn seqalign_scores_of_lists(mut lists: Vec<Option<BlastpHitList>>) -> Vec<Option<BlastpHitList>> {
    for hsp in lists
        .iter_mut()
        .flatten()
        .flat_map(|list| list.hsplist_array.iter_mut())
        .flat_map(|hsp_list| hsp_list.hsps.iter_mut())
    {
        (hsp.e_value, hsp.bit_score) =
            crate::report::serial::seqalign::seqalign_scores(hsp.e_value, hsp.bit_score);
    }
    lists
}

/// The position-based redo of a PSSM search (`pssm_search::redo1`, modes 1 to 3): the matrix
/// information of the PSSM (`PosMatrixInfo`, built when a preliminary list is there to redo),
/// then the redo of every list.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1497
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3254-3258
/// The copy's matrix information reads its frequency ratios (zero where the PSSM has none: then
/// every column takes the matrix's), its ideal block and `kbp_gap_std[0]`'s K.
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
    stats: EvalueStats<'_>,
    link: Option<redo1::RedoLink>,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    redone: &dyn Fn(&[crate::algorithm::blastp::hsp::BlastCompoHeap]) -> anyhow::Result<()>,
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
    let sbp = crate::algorithm::pssm_search::in_pssm::score_block(
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
            "a PSSM whose position-based scaling for the composition-based statistics gives no Karlin-Altschul parameters (NCBI BLAST+ 2.17.0 goes on with a Lambda of -1) is not supported by LOSAT's {}",
            args.app.name()
        );
    };
    redo1::PssmRedo1 {
        pssm: &psi.pssm,
        info: &info,
        context,
        subjects,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        matrix: args.scoring.matrix,
        karlin: &karlin.gapped,
        x_drop_final,
        x_drop_final_bits: args.x_drop_final_bits,
        stats,
        expect_value: args.evalue,
        hitlist_size: args.traceback_hitlist_size(),
        inclusion_ethresh: args.inclusion_ethresh,
        link,
        redone,
    }
    .run(preliminary_hit_lists, compact_hit_lists)
}

/// The warning line of a PSSM query (psiblast `-in_pssm`, and the later rounds of a FASTA query):
/// its label is `GetSeqIdString()` of the query's Seq-id (the PSSM's with `-in_pssm`, the record's
/// local id, DECISIONS 12 of line PS; a FASTA query's later rounds keep its query factory, whose
/// Seq-loc id is the record's Seq-align id: the `warning_id` that `-parse_deflines` read, `P69905`
/// of `sp|P69905.2|`) and the title of `CBlastQuerySourceBioseqSet`
/// (`PsiSearch::warning_title`; a FASTA query's later rounds keep the record's title), and its
/// messages add those of
/// `PsiBlastSetupScoreBlock` (`PsiSearch::setup_messages`), which NCBI keeps with the query's
/// others (sorted, `report::query_warnings::sorted_query_messages`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:188-193
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
    let title = match psi.warning_title.as_deref() {
        Some(title) => std::borrow::Cow::Borrowed(title),
        None => crate::blastinput::fasta_reader::InputRecord::title_bytes(record),
    };
    crate::report::query_warnings::query_warning_with_seq_id(
        program,
        record.warning_id(),
        &title,
        &all,
    )
}

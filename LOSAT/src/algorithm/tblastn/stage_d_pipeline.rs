//! Internal local-subject TBLASTN Stage C to D composition pipeline.
//! The public CLI remains gated while Stage D and Stage E are incomplete.

use std::collections::HashMap;
#[cfg(feature = "parallel")]
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{ensure, Context, Result};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::kappa::{
    convert_distinct_alignments, normalize_postredo_scores, postredo_converted_stats,
    postredo_converted_stats_with_matrix, reap_contained_postredo_hsps, redo_preliminary_match,
    ConvertedKappaHsp,
};
use super::kappa_heap::{compo_early_termination, CompoHeap, CompoHeapRecord};
use super::search_gapped::{
    full_translation_traceback_with_matrix_and_events_with_mask_mode_owned,
    preliminary_protein_hsps_in_ncbi_order, GappedHsp, PreliminaryProfile, TargetTranslation,
    TracebackOwnedHsp,
};
use super::stage_d_kappa_params::{local_extension_final_xdrop, local_kappa_redo_params};
use super::stage_d_linking::{
    compare_preliminary_lists_for_kappa, link_preliminary_hsps, reap_by_evalue, LinkedHsp,
    LinkedHspList,
};
use super::stage_d_results::{KappaHspPayload, KappaResultHitList, KappaResultList};
use super::stage_d_stats::{
    local_parameters_for_call, LocalParameterCall, LocalParameterOptions, LocalSubjectParameters,
};
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::gapalign::GapAlignScratch;
use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::core::composition_adjustment::adjust_scores::{
    read_aa_composition, BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    build_query_word_hashes, BlastCompoAdjustMode, BlastCompoQueryInfo, BlastCompoSequenceData,
    EMatrixAdjustRule,
};
use crate::stats::karlin::evalue_from_raw_score;
use crate::stats::spouge::{blast_spouge_stoe, lookup_protein_gumbel_params, BlastGumbelBlk};
use crate::stats::tables::{
    lookup_protein_params, lookup_protein_params_ungapped, protein_scoring_supported, KarlinParams,
};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::aa_char_to_ncbistdaa;
use crate::utils::seg::SegParams;
use crate::utils::threading::{with_search_pool, SearchPool};

// NCBI c++/src/algo/blast/core/blast_options.c:673-675;
// c++/include/algo/blast/core/blast_options.h:128,144,163:
// gap_x_dropoff = BLAST_GAP_X_DROPOFF_PROT; /* 15 */
// gap_x_dropoff_final = BLAST_GAP_X_DROPOFF_FINAL_PROT; /* 25 */
// #define PSI_INCLUSION_ETHRESH 0.002
#[derive(Clone, Copy)]
pub(super) struct LocalStageDProfile<'a> {
    pub seg: Option<&'a SegParams>,
    pub soft_masking: bool,
    pub mask_lowercase: bool,
    pub genetic_code: u8,
    pub expect_value: f64,
    pub max_target_seqs: usize,
}

// NCBI c++/src/algo/blast/core/blast_options.c:903-936;
// core/blast_parameters.c:302-383,902-999:
// scoring matrix, gap costs, word size, threshold and two-hit window are
// independent options supplied before the preliminary and Stage D calls.
// NCBI c++/src/algo/blast/core/blast_parameters.c:457-463:
// gap_x_dropoff and gap_x_dropoff_final are converted from bits using lambda.
#[derive(Clone, Copy)]
pub(super) struct LocalStageDScoring {
    pub matrix: ScoringMatrix,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub word_size: usize,
    pub threshold: i32,
    pub window: i32,
    pub gap_xdrop_bits: f64,
    pub final_xdrop_bits: f64,
}

impl Default for LocalStageDScoring {
    fn default() -> Self {
        Self {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
            word_size: 3,
            threshold: 13,
            window: 40,
            gap_xdrop_bits: 15.0,
            final_xdrop_bits: 25.0,
        }
    }
}

// NCBI c++/src/algo/blast/core/blast_engine.c:870-905;
// core/blast_kappa.c:3525-3565; core/blast_traceback.c:717-719:
// comparison-only boundary snapshots expose the actual Stage C HSP list,
// Kappa incoming order and ordinary posttraceback stat_length to tests.
// They are never used to produce a search result.
#[derive(Default)]
struct StageDBoundaryTrace {
    parameters: Option<LocalSubjectParameters>,
    // NCBI c++/src/algo/blast/api/blast_results.cpp:72-115:
    // s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
    ungapped_karlin: Vec<KarlinParams>,
    query_validity: Vec<bool>,
    subject_nt_lengths: Vec<usize>,
    preliminary: Vec<(usize, Vec<(usize, GappedHsp)>)>,
    redo: Vec<(usize, usize, Vec<GappedHsp>)>,
    posttraceback_lengths: Vec<(usize, usize, i32)>,
}

// NCBI c++/src/algo/blast/core/blast_setup.c:1011-1024;
// c++/src/algo/blast/format/blast_format.cpp:445-477:
// BLAST_CalcEffLengths(...); BlastHitSavingParametersUpdate(...);
// x_PrintOneQueryFooter(*results.GetAncillaryData());
// Formatting consumes the same initial context lengths computed before search.
pub(super) fn run_local_for_report(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    composition_mode2: bool,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
) -> Result<(
    Vec<KappaResultHitList>,
    LocalSubjectParameters,
    Vec<KarlinParams>,
    Vec<bool>,
)> {
    run_local_for_report_threads(
        queries,
        subjects,
        profile,
        composition_mode2,
        do_sum_stats,
        scoring,
        1,
    )
}

// NCBI c++/src/algo/blast/api/prelim_stage.cpp:145-188:
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
// NCBI c++/src/algo/blast/core/blast_engine.c:1411-1414:
// /* iterate over all subject sequences */
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//    Int4 stat_length;
// The pool owns one search; ordered subject results enter the same Stage D stream.
pub(super) fn run_local_for_report_threads(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    composition_mode2: bool,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    threads: usize,
) -> Result<(
    Vec<KappaResultHitList>,
    LocalSubjectParameters,
    Vec<KarlinParams>,
    Vec<bool>,
)> {
    with_search_pool(threads, "tblastn", |pool| {
        let mut trace = StageDBoundaryTrace::default();
        let results = run_local_search_with_pool(
            queries,
            subjects,
            profile,
            composition_mode2,
            do_sum_stats,
            scoring,
            Some(&mut trace),
            Some(pool),
        )?;
        Ok((
            results,
            trace
                .parameters
                .context("missing TBLASTN initial parameters")?,
            trace.ungapped_karlin,
            trace.query_validity,
        ))
    })
}

// NCBI c++/src/algo/blast/core/blast_engine.c:870-905;
// c++/src/algo/blast/core/blast_kappa.c:3383-3427,3525-3736,2494-2515;
// c++/src/algo/blast/core/blast_traceback.c:1763-1776:
// BLAST_LinkHsps(...); s_Blast_HSPListReapByPrelimEvalue(...);
// BlastCompo_EarlyTermination(...); Blast_RedoOneMatch(...);
// s_HitlistReapContained(...); s_HitlistEvaluateAndPurge(...);
// BlastCompo_HeapInsert(...); Blast_HitListUpdate(...);
// Blast_HSPResultsReverseOrder(results);
// if (BlastSeqSrcGetTotLen(seq_src) > 0) Blast_HSPResultsSortByEvalue(results);
#[allow(dead_code)] // This internal path is compared before the public CLI is enabled.
pub(super) fn run_local_mode2_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        true,
        true,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:1481-1501,1644-1707:
// if (ext_params->options->compositionBasedStats > 0) Blast_RedoAlignmentCore_MT(...);
// else Blast_TracebackFromHSPList(...); Blast_HSPResultsInsertHSPList(...);
#[allow(dead_code)]
pub(super) fn run_local_mode0_no_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        false,
        false,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:1486-1501;
// c++/src/algo/blast/core/blast_kappa.c:411-427:
// compositionBasedStats selects redo independently of hitParams->do_sum_stats.
#[allow(dead_code)]
pub(super) fn run_local_mode0_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        false,
        true,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:1486-1501;
// c++/src/algo/blast/core/blast_kappa.c:411-427:
// compositionBasedStats selects redo independently of hitParams->do_sum_stats.
#[allow(dead_code)]
pub(super) fn run_local_mode2_no_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(
        queries,
        subjects,
        profile,
        true,
        false,
        LocalStageDScoring::default(),
        None,
    )
}

// NCBI c++/src/algo/blast/core/blast_parameters.c:457-463;
// core/blast_kappa.c:2372-2384: selected x-drop bits reach both Stage C
// gapped extension and composition redo.
#[allow(dead_code)]
pub(super) fn run_local_mode2_sum_stats_with_scoring(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    scoring: LocalStageDScoring,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(queries, subjects, profile, true, true, scoring, None)
}

// NCBI c++/src/algo/blast/core/blast_options.c:903-936;
// core/blast_traceback.c:198-250,1486-1501:
// the selected protein scoring options enter ordinary traceback and statistics.
#[allow(dead_code)]
pub(super) fn run_local_mode0_no_sum_stats_with_scoring(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    scoring: LocalStageDScoring,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search(queries, subjects, profile, false, false, scoring, None)
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:1481-1501:
// if (compositionBasedStats > 0) Blast_RedoAlignmentCore_MT(...);
// else Blast_TracebackFromHSPList(...);
fn run_local_search(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    composition_mode2: bool,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    observer: Option<&mut StageDBoundaryTrace>,
) -> Result<Vec<KappaResultHitList>> {
    run_local_search_with_pool(
        queries,
        subjects,
        profile,
        composition_mode2,
        do_sum_stats,
        scoring,
        observer,
        None,
    )
}

// NCBI c++/src/algo/blast/core/blast_engine.c:1411-1414,1469-1475:
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//      status =
//          s_BlastSearchEngineCore(program_number, query, query_info,
//                                  seq_arg.seq, lookup_wrap, gap_align,
// NCBI c++/src/algo/blast/core/blast_engine.c:870-905:
// BLAST_LinkHsps(...); s_Blast_HSPListReapByPrelimEvalue(...);
// Preliminary work uses one immutable subject; linking and stream writes follow OID order.
fn run_local_search_with_pool(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    composition_mode2: bool,
    do_sum_stats: bool,
    scoring: LocalStageDScoring,
    mut observer: Option<&mut StageDBoundaryTrace>,
    pool: Option<&SearchPool<'_>>,
) -> Result<Vec<KappaResultHitList>> {
    ensure!(
        !queries.is_empty(),
        "TBLASTN Stage D requires a protein query"
    );
    ensure!(
        !subjects.is_empty(),
        "TBLASTN Stage D requires a local subject"
    );
    ensure!(
        profile.max_target_seqs > 0,
        "TBLASTN hitlist size must be positive"
    );
    let code = GeneticCode::try_from_id(profile.genetic_code).map_err(anyhow::Error::msg)?;
    let total_nt: usize = subjects.iter().map(Vec::len).sum();
    let min_subject_length = subjects.iter().map(Vec::len).min().unwrap() / 3;
    let query_lengths: Vec<i32> = queries
        .iter()
        .map(|query| i32::try_from(query.len()))
        .collect::<std::result::Result<_, _>>()?;
    // NCBI c++/src/algo/blast/core/blast_options.c:908-936;
    // core/blast_stat.c:2547-2605: reject scoring schemes absent from
    // the pinned NCBI gapped table before any search.
    let spec = ProteinScoringSpec {
        matrix: scoring.matrix,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
    };
    ensure!(
        protein_scoring_supported(&spec),
        "unsupported TBLASTN protein scoring scheme"
    );
    ensure!(
        !composition_mode2 || scoring.matrix == ScoringMatrix::Blosum62,
        "TBLASTN composition redo matrix is unimplemented"
    );
    let gapped = lookup_protein_params(&spec);
    let ungapped = lookup_protein_params_ungapped(scoring.matrix);
    // NCBI c++/src/algo/blast/core/blast_stat.c:2778-2803;
    // c++/src/algo/blast/core/blast_setup.c:770-847:
    // if (Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]))
    //     contexts[context].is_valid = FALSE;
    // if (query_info->contexts[index].is_valid) BLAST_CalcEffLengths(...);
    // Reuse the Stage C score-frequency calculation on the query bytes after
    // BlastSetUp_MaskQuery, which precedes BlastSetup_ScoreBlkInit at
    // pinned blast_setup.c:614-654.
    let working_frames: Vec<_> = queries
        .iter()
        .map(|query| {
            let frame = if profile.soft_masking {
                encode_protein_query_frame_with_seg(query, None)
            } else {
                super::search_seed::encode_tblastn_lookup_query(
                    query,
                    profile.seg,
                    profile.mask_lowercase,
                )
            };
            vec![frame]
        })
        .collect();
    // NCBI c++/src/algo/blast/core/blast_stat.c:2778-2803;
    // c++/src/algo/blast/core/lookup_wrap.c:91-100:
    // Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // BlastAaLookupTableNew(..., lookup_options->word_size, ...);
    // The active matrix determines context validity; the active word size
    // determines the lookup state used by Stage C.
    let (_, contexts) = build_ncbi_lookup_for_profile(
        &working_frames,
        scoring.threshold,
        &ungapped,
        false,
        scoring.matrix,
        scoring.word_size,
    );
    let query_contexts: Vec<_> = queries
        .iter()
        .zip(&contexts)
        .map(|(query, context)| (query.len(), context.is_valid))
        .collect();
    let valid_contexts: Vec<_> = contexts.iter().map(|context| context.is_valid).collect();
    let gapped_by_context = vec![gapped; queries.len()];
    let ungapped_by_context = vec![ungapped; queries.len()];
    // NCBI c++/src/algo/blast/core/blast_stat.c:4091-4171;
    // core/blast_hits.c:1870-1895: gbp is optional; absent gbp selects
    // BLAST_KarlinStoE_simple with the context's effective search space.
    let gumbel = lookup_protein_gumbel_params(&spec, (total_nt / 3) as i64);

    // NCBI c++/src/algo/blast/core/blast_setup.c:964-985,1011-1024;
    // c++/src/algo/blast/core/blast_parameters.c:457-463:
    // BlastHitSavingParametersNew(..., min_subject_length,
    //     options->compositionBasedStats, ...);
    // params->gap_x_dropoff = (Int4)(options->gap_x_dropoff*LN2/min_lambda);
    let parameters = local_parameters_for_call(
        &query_contexts,
        total_nt,
        subjects.len(),
        &gapped_by_context,
        &ungapped_by_context,
        LocalParameterOptions {
            expect_value: profile.expect_value,
            do_sum_stats,
            max_intron_length: 0,
            gap_trigger_bits: 22.0,
            word_xdrop_bits: 7.0,
            scale_factor: 1.0,
            gumbel: gumbel.as_ref(),
        },
        LocalParameterCall::Initial {
            min_subject_length: i32::try_from(min_subject_length)?,
            composition_based_stats: if composition_mode2 { 2 } else { 0 },
        },
    );
    if let Some(ref mut trace) = observer {
        trace.parameters = Some(parameters.clone());
        // NCBI c++/src/algo/blast/api/blast_results.cpp:72-115:
        // s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
        trace.ungapped_karlin = contexts
            .iter()
            .map(|context| context.karlin_params)
            .collect();
        // NCBI blast_results.cpp:82-97: no valid context returns before
        // Karlin block initialization and leaves search space at zero.
        trace.query_validity = contexts.iter().map(|context| context.is_valid).collect();
        trace.subject_nt_lengths = subjects.iter().map(Vec::len).collect();
    }
    let link = parameters.link;
    let word_xdrop: Vec<_> = parameters.cutoffs.iter().map(|v| v.word_xdrop).collect();
    let word_cutoff: Vec<_> = parameters.cutoffs.iter().map(|v| v.word_cutoff).collect();
    let hit_cutoff: Vec<_> = parameters.cutoffs.iter().map(|v| v.hit_cutoff).collect();
    let search = PreliminaryProfile {
        seg: profile.seg,
        soft_masking: profile.soft_masking,
        threshold: scoring.threshold,
        window: scoring.window,
        word_xdrop: &word_xdrop,
        word_cutoff: &word_cutoff,
        mask_lowercase: profile.mask_lowercase,
        matrix: scoring.matrix,
        word_size: scoring.word_size,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        gap_xdrop: ((scoring.gap_xdrop_bits * std::f64::consts::LN_2) / gapped.lambda) as i32,
        gapped_cutoff: &hit_cutoff,
        hsp_num_max: i32::MAX as usize,
    };
    let query_refs: Vec<_> = queries.iter().map(Vec::as_slice).collect();
    // NCBI c++/src/algo/blast/core/blast_hits.c:43-70;
    // c++/src/algo/blast/core/hspfilter_collector.c:105-161:
    // if (compositionBasedStats) {
    //     if (hitlist_size <= 500) prelim_hitlist_size = 1050;
    //     else prelim_hitlist_size = hitlist_size*2 + 50;
    // }
    // Blast_HitListUpdate(results->hitlist_array[index], hsp_list);
    // NCBI c++/src/algo/blast/core/blast_traceback.c:1839-1846;
    // c++/src/algo/blast/core/blast_hspstream.c:816-864:
    // if (ADAPTIVE_CBS_ENV && compositionBasedStats && hitlist_size < 1000)
    //     BlastHSPCBSStreamClose(hsp_stream, hitlist_size);
    // That close can trim the preliminary stream, so this internal path
    // rejects the environment mode until the close operation is ported.
    ensure!(
        !composition_mode2 || std::env::var_os("ADAPTIVE_CBS").is_none(),
        "TBLASTN ADAPTIVE_CBS preliminary stream close is unimplemented"
    );
    // NCBI c++/src/algo/blast/core/blast_hits.c:43-70:
    // if (compositionBasedStats) ... 1050 or 2*hitlist+50;
    // else if (gapped_calculation)
    //   prelim_hitlist_size = MIN(MAX(2*hitlist_size,10), hitlist_size+50);
    let prelim_size = if composition_mode2 {
        if profile.max_target_seqs <= 500 {
            1050
        } else {
            profile
                .max_target_seqs
                .checked_mul(2)
                .and_then(|v| v.checked_add(50))
                .context("TBLASTN preliminary cap overflow")?
        }
    } else {
        profile
            .max_target_seqs
            .saturating_mul(2)
            .max(10)
            .min(profile.max_target_seqs.saturating_add(50))
    };
    let mut preliminary_hitlists: Vec<Option<KappaResultHitList>> =
        (0..queries.len()).map(|_| None).collect();
    // NCBI c++/src/algo/blast/core/blast_engine.c:1469-1475,872-905:
    //      status =
    //          s_BlastSearchEngineCore(program_number, query, query_info,
    //                                  seq_arg.seq, lookup_wrap, gap_align,
    // if (hit_params->link_hsp_params) {
    //     status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
    //               subject->length, gap_align->sbp, hit_params->link_hsp_params,
    //               score_options->gapped_calculation);
    // }
    // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
    // Finish search, link/direct E-value and preliminary reap inside one subject
    // job, before its OID-ordered list enters the shared collector.
    let subject_core = |subject: &[u8]| -> Result<(Vec<(usize, GappedHsp)>, LinkedHspList)> {
        let (preliminary, _) = preliminary_protein_hsps_in_ncbi_order(
            &query_refs,
            subject,
            profile.genetic_code,
            search,
        )?;
        // NCBI c++/src/algo/blast/core/blast_engine.c:870-905:
        // BLAST_LinkHsps(..., hsp_list_out, ...);
        // s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, ...);
        // Even an allocated empty subject list enters link/E-value.
        let mut linked = if let Some(link) = link {
            link_preliminary_hsps(
                &preliminary,
                &query_lengths,
                &parameters.lengths,
                i32::try_from(subject.len())?,
                &gapped_by_context,
                gumbel
                    .as_ref()
                    .context("TBLASTN sum-statistics Spouge state is missing")?,
                &link,
            )?
        } else {
            // NCBI c++/src/algo/blast/core/blast_engine.c:728-729,804-813,882-887:
            // if (subject->length > 0) stat_length = subject->length;
            // Blast_HSPListGetEvalues(..., stat_length, ...);
            let stat_length = initial_translated_stat_length(subject.len());
            direct_local_evalues(
                &preliminary,
                &query_lengths,
                i32::try_from(stat_length)?,
                &gapped_by_context,
                &parameters.lengths,
                gumbel.as_ref(),
            )
        };
        reap_by_evalue(&mut linked, parameters.prelim_evalue);
        Ok((preliminary, linked))
    };
    // NCBI c++/src/algo/blast/core/blast_engine.c:1411-1414,1469-1475:
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //      status =
    //          s_BlastSearchEngineCore(program_number, query, query_info,
    //                                  seq_arg.seq, lookup_wrap, gap_align,
    // NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
    // (*thread)->Run(); (*thread)->Join(&result);
    // Each subject keeps the NCBI search-to-link-to-reap order. Indexed
    // collection finishes before the OID-ordered shared collector below.
    let parallel_selected = pool.is_some_and(SearchPool::enabled) && subjects.len() > 1;
    crate::utils::threading::report_stage(
        "tblastn",
        "subject_core",
        subjects.len(),
        parallel_selected,
    );
    // NCBI c++/src/algo/blast/api/prelim_stage.cpp:145-145:
    // TBlastThreads the_threads(GetNumberOfThreads());
    // NCBI c++/src/algo/blast/core/blast_engine.c:1409-1414:
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    // Keep at most one pool-sized batch of subject results before OID reduction;
    // memory follows NCBI's bounded number of workers, not total subject count.
    #[cfg(feature = "parallel")]
    let mut parallel_batch = std::collections::VecDeque::new();
    #[cfg(feature = "parallel")]
    let diagnose = crate::utils::threading::diagnostics_enabled();
    #[cfg(feature = "parallel")]
    let active = AtomicUsize::new(0);
    #[cfg(feature = "parallel")]
    let peak = AtomicUsize::new(0);
    #[cfg(feature = "parallel")]
    let mut workers = std::collections::BTreeSet::new();
    for (oid, subject) in subjects.iter().enumerate() {
        // NCBI c++/src/algo/blast/core/blast_engine.c:870-905:
        // BLAST_LinkHsps(..., subject->length, ...);
        // s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        #[cfg(feature = "parallel")]
        let (preliminary, linked) = if parallel_selected {
            if parallel_batch.is_empty() {
                let pool = pool.expect("parallel pool selected");
                let end = oid.saturating_add(pool.threads()).min(subjects.len());
                let batch: Vec<_> = pool.install(|| {
                    subjects[oid..end]
                        .par_iter()
                        .map(|subject| {
                            if diagnose {
                                let now = active.fetch_add(1, Ordering::SeqCst) + 1;
                                peak.fetch_max(now, Ordering::SeqCst);
                            }
                            let result = subject_core(subject);
                            if diagnose {
                                active.fetch_sub(1, Ordering::SeqCst);
                            }
                            (result, rayon::current_thread_index())
                        })
                        .collect()
                });
                if diagnose {
                    workers.extend(batch.iter().filter_map(|(_, worker)| *worker));
                }
                parallel_batch.extend(batch);
            }
            parallel_batch
                .pop_front()
                .context("missing ordered TBLASTN subject result")?
                .0?
        } else {
            subject_core(subject)?
        };
        #[cfg(not(feature = "parallel"))]
        let (preliminary, linked) = subject_core(subject)?;
        // NCBI c++/src/algo/blast/core/blast_engine.c:870-905:
        // BLAST_LinkHsps(..., hsp_list_out, ...);
        // s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, ...);
        // The comparison-only trace retains pre-link HSP bytes by OID.
        if let Some(ref mut trace) = observer {
            trace.preliminary.push((oid, preliminary.clone()));
        }
        // NCBI c++/src/algo/blast/core/blast_engine.c:1553-1555:
        // status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
        // NCBI c++/src/algo/blast/core/blast_hspstream.c:289-319:
        // Query-indexed lists are consumed in
        // reverse-sorted order by BlastHSPStreamRead.
        for context in 0..queries.len() {
            let hsps: Vec<_> = linked
                .hsps
                .iter()
                .copied()
                .filter(|h| h.context == context)
                .collect();
            if !hsps.is_empty() {
                let best_evalue = hsps.iter().map(|h| h.evalue).reduce(f64::min).unwrap();
                // NCBI c++/src/algo/blast/core/hspfilter_collector.c:129-161:
                // Blast_HitListUpdate(results->hitlist_array[index], hsp_list);
                // Each context owns an independently capped preliminary list.
                let collector = preliminary_hitlists[context]
                    .get_or_insert(KappaResultHitList::new(prelim_size)?);
                collector.update(KappaResultList {
                    oid: i32::try_from(oid)?,
                    hsps: LinkedHspList { hsps, best_evalue },
                    payloads: Vec::new(),
                })?;
            }
        }
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:1659-1668:
    // /* Use a local diagnostics structure, because the one passed in an input
    //    argument can be shared between multiple threads */
    #[cfg(feature = "parallel")]
    if parallel_selected && diagnose {
        eprintln!(
            "[losat-thread-activity] program=tblastn stage=subject_core work_items={} worker_slots={:?} peak_active={}",
            subjects.len(), workers, peak.load(Ordering::SeqCst)
        );
    }
    // NCBI c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319:
    // Blast_HSPResultsReverseSort(results); stream reads from the end.
    let mut preliminary_lists: Vec<_> = preliminary_hitlists
        .iter()
        .enumerate()
        .flat_map(|(context, hitlist)| {
            hitlist.iter().flat_map(move |hitlist| {
                hitlist
                    .lists()
                    .iter()
                    .map(move |list| (context, list.oid as usize, list.hsps.clone()))
            })
        })
        .collect();
    // NCBI c++/src/algo/blast/core/blast_hspstream.c:91-96,158-204,
    // 569-610: the ordinary traceback stream sorts HSP lists by decreasing
    // OID, then each batch reads from the end; equal-OID query lists were
    // inserted in query order and are consumed in reverse query order.
    // Composition redo uses the score-sorted stream at lines 289-319.
    if composition_mode2 {
        preliminary_lists.sort_by(|(qa, oa, a), (qb, ob, b)| {
            qa.cmp(qb)
                .then_with(|| compare_preliminary_lists_for_kappa(*oa as i32, a, *ob as i32, b))
        });
    } else {
        preliminary_lists.sort_by(|(qa, oa, _), (qb, ob, _)| oa.cmp(ob).then_with(|| qb.cmp(qa)));
    }

    if !composition_mode2 {
        // NCBI c++/src/algo/blast/core/blast_traceback.c:1481-1501,
        // 1644-1707,1763-1776: the control path performs ordinary
        // traceback and posttraceback statistics, then result post-pipes.
        return finish_local_mode0_no_sum_stats(
            queries,
            subjects,
            profile,
            preliminary_lists,
            &query_lengths,
            &gapped_by_context,
            gumbel.as_ref(),
            &code,
            &parameters.lengths,
            link.as_ref(),
            scoring,
            observer.as_deref_mut(),
        );
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:2352-2384,2418-2479,
    // 3099-3110,3290-3302:
    // s_RecordInitialSearch(...); s_GetAlignParams(...);
    // BlastCompo_HeapInitialize(..., hitlist_size, inclusion_ethresh);
    // NCBI c++/src/algo/blast/core/blast_setup.c:614-624;
    // c++/src/algo/blast/core/blast_kappa.c:2309-2337,3244-3248,515-526:
    // if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
    // s_GetQueryInfo(queryBlk->sequence, queryInfo, ...);
    // For TBLASTN s_ComputeNumIdentities passes the working query sequence.
    // Stage C uses these same hard-masked bytes for extension and traceback.
    let encoded: Vec<Vec<u8>> = working_frames
        .iter()
        .map(|frames| frames[0].aa_seq[1..1 + frames[0].aa_len].to_vec())
        .collect();
    // NCBI c++/src/objtools/align_format/tabular.cpp:983-1021:
    // alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
    // if (m_QuerySeq[i] == m_SubjectSeq[i]) ++num_ident;
    // The formatter reads the original query, not the hard-masked Kappa buffer.
    let report_queries: Vec<Vec<u8>> = queries
        .iter()
        .map(|query| {
            query
                .iter()
                .copied()
                .map(|aa| aa_char_to_ncbistdaa(aa.to_ascii_uppercase()))
                .collect()
        })
        .collect();
    let query_infos: Vec<_> = encoded
        .iter()
        .enumerate()
        .map(|(index, seq)| BlastCompoQueryInfo {
            origin: 0,
            seq: BlastCompoSequenceData::from_ncbistdaa(seq),
            composition: read_aa_composition(seq),
            eff_search_space: parameters.lengths[index].eff_searchsp as f64,
            words: Some(build_query_word_hashes(seq)),
        })
        .collect();
    // NCBI c++/src/algo/blast/core/blast_kappa.c:2439-2479:
    // do_link_hsps and cutoff_s follow the actual sum-statistics option.
    let gumbel = gumbel
        .as_ref()
        .context("TBLASTN composition Spouge state is missing")?;
    let redo_params = local_kappa_redo_params(
        ScoringMatrix::Blosum62,
        11,
        1,
        &gapped_by_context,
        &valid_contexts,
        &parameters,
        *query_lengths.iter().max().unwrap(),
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        false,
        profile.expect_value,
        do_sum_stats,
        scoring.final_xdrop_bits,
        local_extension_final_xdrop(
            scoring.gap_xdrop_bits,
            scoring.final_xdrop_bits,
            gapped.lambda,
        )?,
    )?;
    let scaled = KarlinParams {
        lambda: gapped.lambda / 32.0,
        ..gapped
    };
    let scaled_by_context = vec![scaled; queries.len()];
    // NCBI c++/src/algo/blast/core/blast_kappa.c:411-427;
    // core/blast_parameters.c:1245-1253: postredo link, when enabled, uses
    // the updated zero small-gap cutoff; direct E-value uses raw subject length.
    let mut post_link = link;
    if let Some(link) = post_link.as_mut() {
        link.cutoff_small_gap = 0;
    }
    let mut heaps: Vec<_> = (0..queries.len())
        .map(|_| CompoHeap::new(profile.max_target_seqs, 0.002))
        .collect::<Result<_>>()?;
    let mut owned_lists = HashMap::new();
    let mut scratch = GapAlignScratch::new();
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    for (context, oid, list) in preliminary_lists {
        // NCBI c++/src/algo/blast/core/blast_kappa.c:3383-3427,3525-3736:
        // if (BlastCompo_EarlyTermination(localMatch->best_evalue,
        //     redoneMatches, numQueries)) { Blast_HSPListFree(localMatch); continue; }
        // query_index = localMatch->query_index;
        if compo_early_termination(list.best_evalue, &heaps) {
            continue;
        }
        let subject = &subjects[oid];
        if let Some(ref mut trace) = observer {
            trace
                .redo
                .push((context, oid, list.hsps.iter().map(|hsp| hsp.hsp).collect()));
        }
        let incoming: Vec<_> = list.hsps.iter().map(|hsp| hsp.hsp).collect();
        let mut redo = redo_preliminary_match(
            &incoming,
            i32::try_from(context)?,
            &query_infos,
            subject,
            profile.genetic_code,
            &redo_params,
            scaled.lambda,
            ScoringMatrix::Blosum62,
            &mut scratch,
            &mut workspace,
        )?;
        let converted = convert_distinct_alignments(&mut redo.alignments_by_query[context])?;
        if converted.is_empty() {
            continue;
        }
        let mut hsp_values: Vec<_> = converted.iter().map(|h| (h.context, h.hsp)).collect();
        // NCBI c++/src/algo/blast/core/blast_kappa.c:3661-3706:
        // s_HSPListFromDistinctAlignments(...); s_HitlistReapContained(...);
        // s_HitlistEvaluateAndPurge(...); s_HSPListNormalizeScores(...);
        // s_ComputeNumIdentities(...);
        let survivor_indices = reap_contained_postredo_hsps(&mut hsp_values);
        let mut slots: Vec<_> = converted.into_iter().map(Some).collect();
        let survivors: Vec<_> = survivor_indices
            .into_iter()
            .map(|index| slots[index].take().expect("unique containment survivor"))
            .collect();
        // NCBI c++/src/algo/blast/core/blast_kappa.c:411-427,364-368:
        // if (hitParams->do_sum_stats) BLAST_LinkHsps(...);
        // else Blast_HSPListGetEvalues(..., s_GetSubjectLength(...), ..., 1.0);
        let mut postredo = if let Some(link) = post_link.as_ref() {
            link_preliminary_hsps(
                &hsp_values,
                &query_lengths,
                &parameters.lengths,
                i32::try_from(subject.len())?,
                &scaled_by_context,
                gumbel,
                link,
            )?
        } else {
            direct_local_evalues(
                &hsp_values,
                &query_lengths,
                i32::try_from(subject.len())?,
                &scaled_by_context,
                &parameters.lengths,
                Some(gumbel),
            )
        };
        reap_by_evalue(&mut postredo, profile.expect_value);
        if postredo.hsps.is_empty() {
            continue;
        }
        let candidate = CompoHeapRecord {
            subject_index: i32::try_from(oid)?,
            best_evalue: postredo.best_evalue,
            best_score: postredo.hsps[0].hsp.score,
        };
        let bits = normalize_postredo_scores(&mut postredo, scaled.lambda, gapped.k.ln(), 32.0);
        let mut target = TargetTranslation::new(subject, &code);
        let mut owned: Vec<_> = survivors.into_iter().map(Some).collect();
        let payloads = postredo
            .hsps
            .iter()
            .enumerate()
            .map(|(index, hsp)| {
                let converted = owned[hsp.source_index]
                    .take()
                    .context("TBLASTN Kappa source owner reused")?;
                let (ident, positives, length, mismatches, gap_opens, gap_letters) =
                    postredo_converted_stats(&converted, &encoded[context], &mut target)?;
                // NCBI core/blast_kappa.c:515-526 uses the masked working query;
                // objtools/align_format/tabular.cpp:983-1021 recomputes the
                // report counts from the original aligned query sequence.
                let (report_ident, report_positives, _, report_mismatches, _, _) =
                    postredo_converted_stats(&converted, &report_queries[context], &mut target)?;
                Ok(KappaHspPayload {
                    bit_score: bits[index],
                    num_ident: i32::try_from(ident)?,
                    num_positives: positives,
                    report_num_ident: report_ident,
                    report_num_positives: report_positives,
                    report_mismatches,
                    align_length: length,
                    mismatches,
                    gap_opens,
                    gap_letters,
                    edit_script: converted.edit_script,
                    matrix_adjust_rule: converted.matrix_adjust_rule,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        if heaps[context].would_insert(candidate) {
            let discarded = heaps[context].insert(candidate);
            owned_lists.insert(
                (context, candidate.subject_index),
                KappaResultList {
                    oid: candidate.subject_index,
                    hsps: postredo,
                    payloads,
                },
            );
            if let Some(previous) = discarded {
                owned_lists.remove(&(context, previous.subject_index));
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:2494-2515;
    // c++/src/algo/blast/core/blast_hits.c:3243-3297,3383-3400,3420-3437:
    // while ((hsp_list = BlastCompo_HeapPop(heap)) != NULL)
    //     Blast_HitListUpdate(hitlist, hsp_list);
    // Blast_HSPResultsReverseOrder(results);
    // Blast_HSPResultsSortByEvalue(results); /* positive local TotLen */
    let mut results = Vec::with_capacity(queries.len());
    for (context, heap) in heaps.iter_mut().enumerate() {
        let mut hitlist = KappaResultHitList::new(profile.max_target_seqs)?;
        while let Some(popped) = heap.pop() {
            let list = owned_lists
                .remove(&(context, popped.subject_index))
                .context("TBLASTN composition heap lost its HSP payload")?;
            hitlist.update(list)?;
        }
        hitlist.reverse_order();
        if total_nt > 0 {
            hitlist.sort_by_evalue_for_positive_totlen();
        }
        results.push(hitlist);
    }
    Ok(results)
}

// NCBI c++/src/algo/blast/core/blast_engine.c:728-729,804-813,882-887:
// Int4 stat_length = subject->length;
// if (subject->length > 0) stat_length = subject->length;
// Blast_HSPListGetEvalues(..., stat_length, ...);
// NCBI c++/src/algo/blast/core/blast_util.c:508-527,1070-1101:
// if ((length-ABS(frame)+1) < CODON_LENGTH) return prot_length;
// frame_offsets[context+1] = offset after the translated frame and sentinel.
// The final nonempty frame in the six-frame loop supplies stat_length.
fn initial_translated_stat_length(nt_length: usize) -> usize {
    match nt_length {
        0..=2 => nt_length,
        n => ((n - 2) / 3).max(1),
    }
}

// NCBI c++/src/algo/blast/core/blast_hits.c:141-144,1811-1926:
// Blast_HSPNew() uses calloc, so num remains zero without BLAST_LinkHsps.
// hsp->evalue = BLAST_SpougeStoE(score, kbp, sbp->gbp,
//   query_info->contexts[hsp->context].query_length, subject_length);
// hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
fn direct_local_evalues(
    input: &[(usize, super::search_gapped::GappedHsp)],
    query_lengths: &[i32],
    subject_stat_length: i32,
    gapped_params: &[KarlinParams],
    lengths: &[super::stage_d_stats::LocalContextLength],
    gumbel: Option<&BlastGumbelBlk>,
) -> LinkedHspList {
    let hsps: Vec<_> = input
        .iter()
        .enumerate()
        .map(|(source_index, &(context, hsp))| LinkedHsp {
            context,
            hsp,
            num: 0,
            // NCBI c++/src/algo/blast/core/blast_hits.c:1870-1895:
            // if (sbp->gbp) BLAST_SpougeStoE(...);
            // else BLAST_KarlinStoE_simple(score, kbp, eff_searchsp);
            evalue: if let Some(gumbel) = gumbel {
                blast_spouge_stoe(
                    hsp.score,
                    &gapped_params[context],
                    gumbel,
                    query_lengths[context],
                    subject_stat_length,
                )
            } else {
                evalue_from_raw_score(
                    hsp.score,
                    &gapped_params[context],
                    lengths[context].eff_searchsp as f64,
                )
            },
            source_index,
        })
        .collect();
    let best_evalue = hsps
        .iter()
        .map(|h| h.evalue)
        .reduce(f64::min)
        .unwrap_or(0.0);
    LinkedHspList { hsps, best_evalue }
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:200-251,675-721,
// 1644-1707,1763-1776; core/blast_hits.c:1811-1935,3243-3297:
// Blast_TracebackFromHSPList(...); s_HSPListPostTracebackUpdate(...);
// Blast_HSPResultsInsertHSPList(...); Blast_HSPResultsSortByEvalue(...);
fn finish_local_mode0_no_sum_stats(
    queries: &[Vec<u8>],
    subjects: &[Vec<u8>],
    profile: LocalStageDProfile<'_>,
    preliminary_lists: Vec<(usize, usize, LinkedHspList)>,
    query_lengths: &[i32],
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    code: &GeneticCode,
    lengths: &[super::stage_d_stats::LocalContextLength],
    link: Option<&super::stage_d_stats::LocalLinkParameters>,
    scoring: LocalStageDScoring,
    mut observer: Option<&mut StageDBoundaryTrace>,
) -> Result<Vec<KappaResultHitList>> {
    let mut results: Vec<_> = (0..queries.len())
        .map(|_| KappaResultHitList::new(profile.max_target_seqs))
        .collect::<Result<_>>()?;
    // NCBI c++/src/algo/blast/core/blast_traceback.c:583-605;
    // c++/src/algo/blast/core/blast_encoding.c:110-121:
    // the traceback identity pass reads the unmasked amino-acid query;
    // lower-case FASTA residues represent the same protein letters.
    let query_nomask: Vec<_> = queries
        .iter()
        .map(|query| {
            query
                .iter()
                .copied()
                .map(|aa| aa_char_to_ncbistdaa(aa.to_ascii_uppercase()))
                .collect::<Vec<_>>()
        })
        .collect();
    for (context, oid, preliminary) in preliminary_lists {
        let subject = &subjects[oid];
        let input: Vec<_> = preliminary.hsps.iter().map(|linked| linked.hsp).collect();
        // NCBI c++/src/algo/blast/core/blast_traceback.c:503-536:
        // traceback uses the selected scoring matrix and gap costs.
        let (traced, stat_length) =
            full_translation_traceback_with_matrix_and_events_with_mask_mode_owned(
                &queries[context],
                subject,
                profile.genetic_code,
                &input,
                scoring.matrix,
                scoring.gap_open,
                scoring.gap_extend,
                ((scoring.final_xdrop_bits * std::f64::consts::LN_2)
                    / gapped_params[context].lambda) as i32,
                0.0,
                0,
                profile.seg,
                profile.soft_masking,
                profile.mask_lowercase,
                None,
                None,
                None,
                None,
            )?;
        if let Some(ref mut trace) = observer {
            trace
                .posttraceback_lengths
                .push((context, oid, stat_length));
        }
        if traced.is_empty() {
            continue;
        }
        let input: Vec<_> = traced.iter().map(|row| (context, row.hsp)).collect();
        // NCBI c++/src/algo/blast/core/blast_traceback.c:213-250,717-719:
        // if (hit_params->link_hsp_params) BLAST_LinkHsps(..., stat_length);
        // else Blast_HSPListGetEvalues(..., stat_length, ...);
        let mut post = if let Some(link) = link {
            let mut post_link = *link;
            post_link.cutoff_small_gap = 0;
            link_preliminary_hsps(
                &input,
                query_lengths,
                lengths,
                stat_length,
                gapped_params,
                gumbel.context("TBLASTN sum-statistics Spouge state is missing")?,
                &post_link,
            )?
        } else {
            direct_local_evalues(
                &input,
                query_lengths,
                stat_length,
                gapped_params,
                lengths,
                gumbel,
            )
        };
        reap_by_evalue(&mut post, profile.expect_value);
        if post.hsps.is_empty() {
            continue;
        }
        let mut target = TargetTranslation::new(subject, code);
        let mut owned: Vec<Option<TracebackOwnedHsp>> = traced.into_iter().map(Some).collect();
        let payloads = post
            .hsps
            .iter()
            .map(|hsp| {
                let original = owned[hsp.source_index]
                    .take()
                    .context("TBLASTN ordinary traceback source owner reused")?;
                let converted = ConvertedKappaHsp {
                    context,
                    hsp: original.hsp,
                    edit_script: original.edit_script,
                    matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
                };
                let (ident, positives, length, mismatches, gap_opens, gap_letters) =
                    postredo_converted_stats_with_matrix(
                        &converted,
                        &query_nomask[context],
                        &mut target,
                        scoring.matrix,
                    )?;
                Ok(KappaHspPayload {
                    bit_score: (hsp.hsp.score as f64 * gapped_params[context].lambda
                        - gapped_params[context].k.ln())
                        / std::f64::consts::LN_2,
                    num_ident: i32::try_from(ident)?,
                    num_positives: positives,
                    report_num_ident: ident,
                    report_num_positives: positives,
                    report_mismatches: mismatches,
                    align_length: length,
                    mismatches,
                    gap_opens,
                    gap_letters,
                    edit_script: converted.edit_script,
                    matrix_adjust_rule: converted.matrix_adjust_rule,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        results[context].update(KappaResultList {
            oid: i32::try_from(oid)?,
            hsps: post,
            payloads,
        })?;
    }
    // NCBI c++/src/algo/blast/core/blast_traceback.c:1763-1776:
    // if (BlastSeqSrcGetTotLen(seq_src) > 0)
    //     Blast_HSPResultsSortByEvalue(results);
    for result in &mut results {
        result.sort_by_evalue_for_positive_totlen();
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI c++/src/algo/blast/core/blast_engine.c:728-729,804-813;
    // c++/src/algo/blast/core/blast_util.c:508-527,1070-1101:
    // stat_length retains the raw nucleotide length until a nonempty frame
    // replaces it, and the final translated context is frame -3.
    #[test]
    fn initial_stat_length_uses_last_nonempty_translated_frame() {
        for (nt_length, expected) in [
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 1),
            (4, 1),
            (5, 1),
            (6, 1),
            (7, 1),
            (8, 2),
            (360, 119),
            (362, 120),
        ] {
            assert_eq!(initial_translated_stat_length(nt_length), expected);
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:870-905;
    // core/blast_kappa.c:3525-3565; core/blast_traceback.c:717-719:
    // compare the actual Rust C-to-D list and Kappa incoming order with the
    // saved function-entry snapshots, before comparing retained results.
    fn assert_local_stage_d_boundaries(
        trace: &StageDBoundaryTrace,
        calls: &str,
        kappa: Option<&str>,
        do_sum_stats: bool,
    ) {
        let initial_phase = if do_sum_stats {
            "link_before"
        } else {
            "evalue_before"
        };
        // NCBI core/blast_engine.c:870-905; blast_traceback.c:1644-1707:
        // composition redo begins only after all preliminary lists were
        // written. In the ordinary path, the first posttraceback E-value
        // input has recomputed positive num_ident; initial HSPs have zero.
        let boundary = if kappa.is_some() {
            calls
                .lines()
                .find(|line| line.starts_with("D_CALL\t") && line.contains("\tredo\t"))
                .and_then(|line| calls.find(line))
                .expect("NCBI redo call")
        } else {
            calls
                .lines()
                .find(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    f.len() > 16
                        && f[0] == "D_HSP"
                        && f[2] == initial_phase
                        && f[11].parse::<i32>().unwrap() > 0
                })
                .and_then(|line| {
                    let event = line.split('\t').nth(1).unwrap();
                    calls.find(&format!("D_CALL\t{event}\t"))
                })
                .unwrap_or(calls.len())
        };
        let preliminary_calls = &calls[..boundary];
        let actual: Vec<_> = trace
            .preliminary
            .iter()
            .flat_map(|(_, hsps)| hsps.iter())
            .map(|(context, hsp)| {
                (
                    *context as i32,
                    hsp.frame as i32,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.score,
                    hsp.s_gapped_start,
                )
            })
            .collect();
        let expected: Vec<_> = preliminary_calls
            .lines()
            .filter_map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                if f.len() > 16 && f[0] == "D_HSP" && f[2] == initial_phase {
                    Some((
                        f[4].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                        f[16].parse().unwrap(),
                    ))
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(actual, expected, "Stage C to D HSP input/order/count");
        // NCBI core/blast_engine.c:870-905; link_hsps.c:1765-1810:
        // each allocated subject list calls link with raw nucleotide length,
        // or direct E-values with translated length, before preliminary reap.
        let initial_call = if do_sum_stats { "link" } else { "evalue" };
        let call_rows: Vec<Vec<_>> = preliminary_calls
            .lines()
            .filter_map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (f.len() > 4 && f[0] == "D_CALL" && f[2] == initial_call).then_some(f)
            })
            .collect();
        assert_eq!(
            call_rows.len(),
            trace.preliminary.len(),
            "Stage C to D initial function call count"
        );
        for (row, (oid, _)) in call_rows.iter().zip(&trace.preliminary) {
            let nt_length = trace.subject_nt_lengths[*oid];
            assert_eq!(
                row[4].parse::<usize>().unwrap(),
                if do_sum_stats {
                    nt_length
                } else {
                    initial_translated_stat_length(nt_length)
                },
                "Stage C to D subject length at initial call"
            );
            if do_sum_stats {
                let link = trace.parameters.as_ref().unwrap().link.unwrap();
                assert_eq!(row[6].parse::<i32>().unwrap(), link.longest_intron);
                assert_eq!(
                    row[7].parse::<f64>().unwrap().to_bits(),
                    link.gap_decay_rate.to_bits()
                );
                assert_eq!(row[8].parse::<i32>().unwrap(), link.gap_size);
                assert_eq!(row[9].parse::<i32>().unwrap(), link.cutoff_small_gap);
            }
        }
        if let Some(params) = &trace.parameters {
            for line in calls.lines().filter(|line| line.starts_with("D_CONTEXT\t")) {
                let f: Vec<_> = line.split('\t').collect();
                let context: usize = f[2].parse().unwrap();
                assert_eq!(
                    params.lengths[context].eff_searchsp,
                    f[4].parse().unwrap(),
                    "D effective search space"
                );
                assert_eq!(
                    params.lengths[context].length_adjustment,
                    f[5].parse().unwrap(),
                    "D length adjustment"
                );
            }
        }
        if let Some(kappa) = kappa {
            let actual: Vec<_> = trace
                .redo
                .iter()
                .enumerate()
                .flat_map(|(event, (context, _, hsps))| {
                    hsps.iter().enumerate().map(move |(index, hsp)| {
                        (
                            event,
                            index,
                            hsp.score,
                            *context,
                            hsp.frame as i32,
                            hsp.q_start,
                            hsp.q_end,
                            hsp.q_gapped_start,
                            hsp.s_start,
                            hsp.s_end,
                            hsp.s_gapped_start,
                        )
                    })
                })
                .collect();
            let expected: Vec<_> = kappa
                .lines()
                .filter(|line| line.starts_with("K_TRACE_PRELIM\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                        f[3].parse().unwrap(),
                        f[4].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                        f[11].parse().unwrap(),
                    )
                })
                .collect();
            assert_eq!(actual, expected, "Kappa redo incoming HSPs/order");
        } else {
            // NCBI core/blast_traceback.c:717-719: stat_length is passed
            // to the posttraceback E-value call. D_HSP field 16 is instead
            // subject.gapped_start (ncbi_d_call_trace.c:113-118).
            let mut seen = std::collections::HashSet::new();
            let mut current_length = None;
            let mut expected = Vec::new();
            for line in calls.lines() {
                let f: Vec<_> = line.split('\t').collect();
                if f.len() > 4
                    && f[0] == "D_CALL"
                    && f[2] == if do_sum_stats { "link" } else { "evalue" }
                {
                    current_length = Some(f[4].parse::<i32>().unwrap());
                } else if f.len() > 16
                    && f[0] == "D_HSP"
                    && f[2] == "bits_after"
                    && seen.insert(f[1])
                {
                    expected.push((f[4].parse::<usize>().unwrap(), current_length.unwrap()));
                }
            }
            let actual: Vec<_> = trace
                .posttraceback_lengths
                .iter()
                .map(|(context, _, length)| (*context, *length))
                .collect();
            assert_eq!(
                actual, expected,
                "posttraceback translated stat_length/call order"
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:3525-3736,2494-2515;
    // c++/src/algo/blast/core/blast_traceback.c:1763-1776:
    // The 112-subject comparison-only local run supplies natural Stage C
    // candidates and two retained HSPs after Kappa and the result post-pipes.
    #[test]
    fn natural_112_subject_local_path_matches_retained_ncbi_fields() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<Vec<u8>> {
            let mut records: Vec<Vec<u8>> = Vec::new();
            for line in std::fs::read_to_string(format!("{root}{name}"))
                .unwrap()
                .lines()
            {
                if line.starts_with('>') {
                    records.push(Vec::new());
                } else {
                    records.last_mut().unwrap().extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let seg = SegParams::default();
        let mut boundary = StageDBoundaryTrace::default();
        let results = run_local_search(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 2,
            },
            true,
            true,
            LocalStageDScoring::default(),
            Some(&mut boundary),
        )
        .unwrap();
        let calls =
            std::fs::read_to_string(format!("{root}natural_c_d_early_20260925/calls.tsv")).unwrap();
        let kappa =
            std::fs::read_to_string(format!("{root}result_order_20260925/ncbi.trace")).unwrap();
        assert_local_stage_d_boundaries(&boundary, &calls, Some(&kappa), true);
        assert_eq!(results.len(), 1);
        let hitlists = results[0].lists();
        assert_eq!(hitlists.iter().map(|h| h.oid).collect::<Vec<_>>(), [10, 11]);
        let report =
            std::fs::read_to_string(format!("{root}report_payload_20260925/report_fields.tsv"))
                .unwrap();
        for (list, row) in hitlists.iter().zip(report.lines()) {
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(list.payloads.len(), 1);
            let hsp = &list.hsps.hsps[0];
            let payload = &list.payloads[0];
            assert_eq!(f[1], format!("weak_{}", list.oid + 7));
            assert_eq!(hsp.hsp.score, f[2].parse().unwrap());
            assert_eq!(payload.num_ident, f[3].parse().unwrap());
            assert_eq!(payload.num_positives, f[4].parse().unwrap());
            assert_eq!(payload.align_length, f[5].parse().unwrap());
            assert_eq!(payload.mismatches, f[6].parse().unwrap());
            assert_eq!(payload.gap_letters, f[7].parse().unwrap());
            assert_eq!(payload.gap_opens, f[8].parse().unwrap());
            assert_eq!(hsp.hsp.q_start + 1, f[9].parse().unwrap());
            assert_eq!(hsp.hsp.q_end, f[10].parse().unwrap());
            assert_eq!(
                3 * hsp.hsp.s_start + i32::from(hsp.hsp.frame),
                f[11].parse().unwrap()
            );
            assert_eq!(
                3 * hsp.hsp.s_end + i32::from(hsp.hsp.frame) - 1,
                f[12].parse().unwrap()
            );
            assert_eq!(hsp.hsp.frame, f[13].parse().unwrap());
            assert!(!payload.edit_script.is_empty());
        }
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:3577-3736;
    // c++/src/algo/blast/core/blast_traceback.c:1763-1776:
    // Three query contexts share one local nucleotide subject; the third
    // contains only X and NCBI leaves it invalid with no retained HSPs.
    #[test]
    fn natural_multi_query_local_path_matches_retained_ncbi_scores_and_coordinates() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let read_fasta = |name: &str| -> Vec<Vec<u8>> {
            let mut records: Vec<Vec<u8>> = Vec::new();
            for line in std::fs::read_to_string(format!("{root}{name}"))
                .unwrap()
                .lines()
            {
                if line.starts_with('>') {
                    records.push(Vec::new());
                } else {
                    records.last_mut().unwrap().extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let seg = SegParams::default();
        let mut boundary = StageDBoundaryTrace::default();
        let results = run_local_search(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 500,
            },
            true,
            true,
            LocalStageDScoring::default(),
            Some(&mut boundary),
        )
        .unwrap();
        let calls = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.trace"
        ))
        .unwrap();
        let kappa = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/multi_query_20260924_default.tsv"
        )).unwrap();
        assert_local_stage_d_boundaries(&boundary, &calls, Some(&kappa), true);
        let expected = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.out"
        ))
        .unwrap();
        let mut rows = Vec::new();
        for (query_index, hitlist) in results.iter().enumerate() {
            for list in hitlist.lists() {
                for hsp in &list.hsps.hsps {
                    let frame = i32::from(hsp.hsp.frame);
                    let subject_len = subjects[list.oid as usize].len() as i32;
                    let (sstart, send) = if frame > 0 {
                        (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                    } else {
                        (
                            subject_len - 3 * hsp.hsp.s_start + frame + 1,
                            subject_len - 3 * hsp.hsp.s_end + frame + 2,
                        )
                    };
                    rows.push((
                        query_index,
                        hsp.hsp.score,
                        hsp.hsp.q_start + 1,
                        hsp.hsp.q_end,
                        sstart,
                        send,
                        frame,
                    ));
                }
            }
        }
        let expected_rows: Vec<_> = expected
            .lines()
            .map(|row| {
                let f: Vec<_> = row.split('\t').collect();
                let query_index = match f[0] {
                    "full" => 0,
                    "internal70" => 1,
                    _ => 2,
                };
                (
                    query_index,
                    f[2].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                )
            })
            .collect();
        assert_eq!(rows, expected_rows);
        // NCBI c++/src/algo/blast/core/blast_kappa.c:3690-3713;
        // c++/src/algo/blast/core/blast_hits.c:3243-3297:
        // s_HSPListNormalizeScores(...); s_ComputeNumIdentities(...);
        // Blast_HitListUpdate(hitlist, hsp_list);
        // Compare the unformatted per-HSP values after the full local path.
        let trace = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/multi_query_20260924_default.tsv")).unwrap();
        let heap_rows: std::collections::HashMap<_, _> = trace
            .lines()
            .filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t"))
            .map(|row| {
                let f: Vec<_> = row.split('\t').collect();
                (
                    (
                        f[7].parse::<usize>().unwrap(),
                        f[8].parse::<i8>().unwrap(),
                        f[9].parse::<i32>().unwrap(),
                        f[10].parse::<i32>().unwrap(),
                        f[11].parse::<i32>().unwrap(),
                        f[12].parse::<i32>().unwrap(),
                    ),
                    f,
                )
            })
            .collect();
        for (context, hitlist) in results.iter().enumerate() {
            for list in hitlist.lists() {
                for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                    let key = (
                        context,
                        hsp.hsp.frame,
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                    );
                    let row = &heap_rows[&key];
                    assert_eq!(hsp.hsp.score, row[3].parse().unwrap());
                    assert_eq!(
                        payload.bit_score.to_bits(),
                        row[4].parse::<f64>().unwrap().to_bits()
                    );
                    assert_eq!(
                        hsp.evalue.to_bits(),
                        row[5].parse::<f64>().unwrap().to_bits()
                    );
                    assert_eq!(payload.num_ident, row[6].parse().unwrap());
                }
            }
        }
    }
    // NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:800-812;
    // core/blast_engine.c:1460-1466; core/blast_kappa.c:3687-3736:
    // The comparison-only API uses FindGeneticCode(32) with the CLI-calibrated
    // local seqSrc state and gives one retained HSP with full-precision values.
    #[test]
    fn code32_local_path_matches_cli_calibrated_api_oracle() {
        let input = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let oracle = concat!(env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/code32_local_api_20260925_cli_calibrated/code32.trace");
        let read = |name: &str| -> Vec<u8> {
            std::fs::read_to_string(format!("{input}{name}"))
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('>'))
                .flat_map(str::bytes)
                .collect()
        };
        let queries = vec![read("query.faa")];
        let subjects = vec![read("subject_code32.fna")];
        let seg = SegParams::default();
        let results = run_local_mode2_sum_stats(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 32,
                expect_value: 10.0,
                max_target_seqs: 500,
            },
        )
        .unwrap();
        let trace = std::fs::read_to_string(oracle).unwrap();
        let row: Vec<_> = trace
            .lines()
            .find(|row| row.starts_with("K_TRACE_HEAP_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].lists().len(), 1);
        let list = &results[0].lists()[0];
        assert_eq!(list.hsps.hsps.len(), 1);
        let hsp = &list.hsps.hsps[0];
        let payload = &list.payloads[0];
        assert_eq!(hsp.hsp.score, row[3].parse().unwrap());
        assert_eq!(
            payload.bit_score.to_bits(),
            row[4].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            hsp.evalue.to_bits(),
            row[5].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(payload.num_ident, row[6].parse().unwrap());
        assert_eq!(hsp.hsp.frame, row[8].parse().unwrap());
        assert_eq!(
            (
                hsp.hsp.q_start,
                hsp.hsp.q_end,
                hsp.hsp.s_start,
                hsp.hsp.s_end
            ),
            (
                row[9].parse().unwrap(),
                row[10].parse().unwrap(),
                row[11].parse().unwrap(),
                row[12].parse().unwrap()
            )
        );
        assert!(!payload.edit_script.is_empty());
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:200-251,
    // 1481-1501,1644-1707,1763-1776:
    // without composition statistics, ordinary traceback recomputes E-values,
    // reaps, assigns bit scores, and inserts HSP lists before result sorting.
    #[test]
    fn ordinary_traceback_controls_match_local_ncbi_full_precision_fields() {
        for case in [
            "seg_hard_query_20260924",
            "multi_query_20260924",
            "run_20260923",
        ] {
            let input_root = format!(
                "{}/../docs/evidence/tlosan_stage_c/{case}",
                env!("CARGO_MANIFEST_DIR")
            );
            let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
                let mut records: Vec<(String, Vec<u8>)> = Vec::new();
                for line in std::fs::read_to_string(format!("{input_root}/{name}"))
                    .unwrap()
                    .lines()
                {
                    if let Some(id) = line.strip_prefix('>') {
                        records
                            .push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                    } else {
                        records.last_mut().unwrap().1.extend(line.bytes());
                    }
                }
                records
            };
            let queries = read_fasta("query.faa");
            let subjects = read_fasta("subjects.fna");
            let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
            let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
            let seg = SegParams::default();
            let results = run_local_mode0_no_sum_stats(
                &query_sequences,
                &subject_sequences,
                LocalStageDProfile {
                    seg: Some(&seg),
                    soft_masking: false,
                    mask_lowercase: false,
                    genetic_code: 1,
                    expect_value: 10.0,
                    max_target_seqs: 500,
                },
            )
            .unwrap();
            let root = concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../docs/evidence/tlosan_stage_d/run_20260924/"
            );
            let trace = std::fs::read_to_string(format!("{root}{case}_control.trace")).unwrap();
            let expected: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("D_HSP\t") && line.contains("\tbits_after\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        (
                            f[4].parse::<usize>().unwrap(),
                            f[5].parse::<i8>().unwrap(),
                            f[6].parse::<i32>().unwrap(),
                            f[7].parse::<i32>().unwrap(),
                            f[8].parse::<i32>().unwrap(),
                            f[9].parse::<i32>().unwrap(),
                        ),
                        (
                            f[10].parse::<i32>().unwrap(),
                            f[11].parse::<i32>().unwrap(),
                            f[12].parse::<i32>().unwrap(),
                            f[13].parse::<f64>().unwrap().to_bits(),
                            f[14].parse::<f64>().unwrap().to_bits(),
                        ),
                    )
                })
                .collect();
            let mut actual: Vec<_> = results
                .iter()
                .enumerate()
                .flat_map(|(query_index, hitlist)| {
                    hitlist.lists().iter().flat_map(move |list| {
                        list.hsps
                            .hsps
                            .iter()
                            .zip(&list.payloads)
                            .map(move |(hsp, payload)| {
                                (
                                    (
                                        query_index,
                                        hsp.hsp.frame,
                                        hsp.hsp.q_start,
                                        hsp.hsp.q_end,
                                        hsp.hsp.s_start,
                                        hsp.hsp.s_end,
                                    ),
                                    (
                                        hsp.hsp.score,
                                        payload.num_ident,
                                        hsp.num,
                                        hsp.evalue.to_bits(),
                                        payload.bit_score.to_bits(),
                                    ),
                                )
                            })
                    })
                })
                .collect();
            let mut expected = expected;
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(actual, expected, "{case} full-precision ordinary traceback");
            // NCBI c++/src/algo/blast/core/blast_traceback.c:1763-1776;
            // c++/src/algo/blast/core/blast_hits.c:3383-3400:
            // positive local TotLen sorts result lists by best E-value.
            // Here the saved diagnostic outfmt-6 columns witness Stage D
            // order, coordinates and alignment length; Stage E byte output
            // remains gated.
            let expected_output =
                std::fs::read_to_string(format!("{root}{case}_control.out")).unwrap();
            let expected_order: Vec<_> = expected_output
                .lines()
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[0].to_owned(),
                        f[1].to_owned(),
                        f[2].parse::<i32>().unwrap(),
                        f[5].parse::<i32>().unwrap(),
                        f[6].parse::<i32>().unwrap(),
                        f[7].parse::<i32>().unwrap(),
                        f[8].parse::<i32>().unwrap(),
                        f[9].parse::<i8>().unwrap(),
                        f[10].parse::<usize>().unwrap(),
                    )
                })
                .collect();
            let mut actual_order = Vec::new();
            for (query_index, hitlist) in results.iter().enumerate() {
                for list in hitlist.lists() {
                    for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                        let frame = i32::from(hsp.hsp.frame);
                        let subject_len = subjects[list.oid as usize].1.len() as i32;
                        let (sstart, send) = if frame > 0 {
                            (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                        } else {
                            (
                                subject_len - 3 * hsp.hsp.s_start + frame + 1,
                                subject_len - 3 * hsp.hsp.s_end + frame + 2,
                            )
                        };
                        actual_order.push((
                            queries[query_index].0.clone(),
                            subjects[list.oid as usize].0.clone(),
                            hsp.hsp.score,
                            hsp.hsp.q_start + 1,
                            hsp.hsp.q_end,
                            sstart,
                            send,
                            hsp.hsp.frame,
                            payload.align_length,
                        ));
                    }
                }
            }
            assert_eq!(
                actual_order, expected_order,
                "{case} result order/report-ready fields"
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:198-250,717-719;
    // core/blast_kappa.c:3670-3733; core/blast_hits.c:3243-3297:
    // every retained local Stage C FASTA fixture must preserve its final HSP
    // values and order through the ordinary or composition result branch.
    #[test]
    fn remaining_five_local_fixture_families_match_ncbi_stage_d() {
        let croot = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/remaining_local_20260925/run_20260925/"
        );
        for case in [
            "ambiguity_20260923",
            "lowercase_20260923",
            "multi_hsp_20260924",
            "multi_hsp_no_fence_20260924",
            "six_frame_merged_20260924",
        ] {
            let read_fasta = |file: &str| -> Vec<(String, Vec<u8>)> {
                let mut rows = Vec::new();
                for line in std::fs::read_to_string(format!("{croot}{case}/{file}"))
                    .unwrap()
                    .lines()
                {
                    if let Some(id) = line.strip_prefix('>') {
                        rows.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                    } else {
                        rows.last_mut().unwrap().1.extend(line.bytes());
                    }
                }
                rows
            };
            let queries = read_fasta("query.faa");
            let subjects = read_fasta("subjects.fna");
            let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
            let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
            let mut scoring = LocalStageDScoring::default();
            if case == "multi_hsp_no_fence_20260924" {
                scoring.gap_xdrop_bits = 5.0;
                scoring.final_xdrop_bits = 5.0;
            }
            for mode in [0, 2] {
                let name = format!("{case}.mode{mode}");
                let profile = LocalStageDProfile {
                    seg: None,
                    soft_masking: false,
                    mask_lowercase: case == "lowercase_20260923",
                    genetic_code: 1,
                    expect_value: 10000.0,
                    max_target_seqs: 500,
                };
                let mut boundary = StageDBoundaryTrace::default();
                let results = run_local_search(
                    &query_sequences,
                    &subject_sequences,
                    profile,
                    mode == 2,
                    mode == 2,
                    scoring,
                    Some(&mut boundary),
                )
                .unwrap();
                assert_eq!(results.len(), 1, "{name} query count");
                let mut actual = Vec::new();
                let mut ordered = Vec::new();
                for list in results[0].lists() {
                    for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                        actual.push((
                            hsp.hsp.frame,
                            hsp.hsp.q_start,
                            hsp.hsp.q_end,
                            hsp.hsp.s_start,
                            hsp.hsp.s_end,
                            hsp.hsp.score,
                            payload.num_ident,
                            hsp.num,
                            hsp.evalue.to_bits(),
                            payload.bit_score.to_bits(),
                        ));
                        let frame = i32::from(hsp.hsp.frame);
                        let slen = subjects[list.oid as usize].1.len() as i32;
                        let (sstart, send) = if frame > 0 {
                            (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                        } else {
                            (
                                slen - 3 * hsp.hsp.s_start + frame + 1,
                                slen - 3 * hsp.hsp.s_end + frame + 2,
                            )
                        };
                        ordered.push((
                            (queries[0].0.clone(), subjects[list.oid as usize].0.clone()),
                            (
                                hsp.hsp.score,
                                payload.report_num_ident,
                                payload.report_num_positives,
                                payload.align_length,
                                payload.report_mismatches,
                                payload.gap_letters,
                                payload.gap_opens,
                            ),
                            (
                                hsp.hsp.q_start + 1,
                                hsp.hsp.q_end,
                                sstart,
                                send,
                                hsp.hsp.frame,
                            ),
                        ));
                    }
                }
                let trace = std::fs::read_to_string(format!(
                    "{root}{name}.{}",
                    if mode == 0 {
                        "calls.trace"
                    } else {
                        "kappa.trace"
                    }
                ))
                .unwrap();
                let calls = std::fs::read_to_string(format!("{root}{name}.calls.trace")).unwrap();
                assert_local_stage_d_boundaries(
                    &boundary,
                    &calls,
                    if mode == 2 { Some(&trace) } else { None },
                    mode == 2,
                );
                let mut expected: Vec<_> = trace
                    .lines()
                    .filter(|line| {
                        if mode == 0 {
                            line.starts_with("D_HSP\t") && line.contains("\tbits_after\t")
                        } else {
                            line.starts_with("K_TRACE_HEAP_HSP\t")
                        }
                    })
                    .map(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        if mode == 0 {
                            (
                                f[5].parse().unwrap(),
                                f[6].parse().unwrap(),
                                f[7].parse().unwrap(),
                                f[8].parse().unwrap(),
                                f[9].parse().unwrap(),
                                f[10].parse().unwrap(),
                                f[11].parse().unwrap(),
                                f[12].parse().unwrap(),
                                f[13].parse::<f64>().unwrap().to_bits(),
                                f[14].parse::<f64>().unwrap().to_bits(),
                            )
                        } else {
                            (
                                f[8].parse().unwrap(),
                                f[9].parse().unwrap(),
                                f[10].parse().unwrap(),
                                f[11].parse().unwrap(),
                                f[12].parse().unwrap(),
                                f[3].parse().unwrap(),
                                f[6].parse().unwrap(),
                                f[13].parse().unwrap(),
                                f[5].parse::<f64>().unwrap().to_bits(),
                                f[4].parse::<f64>().unwrap().to_bits(),
                            )
                        }
                    })
                    .collect();
                actual.sort_unstable();
                expected.sort_unstable();
                assert_eq!(actual, expected, "{name} full-precision internal HSPs");
                let expected_order: Vec<_> =
                    std::fs::read_to_string(format!("{root}{name}.report.tsv"))
                        .unwrap()
                        .lines()
                        .map(|line| {
                            let f: Vec<_> = line.split('\t').collect();
                            (
                                (f[0].to_owned(), f[1].to_owned()),
                                (
                                    f[2].parse().unwrap(),
                                    f[5].parse().unwrap(),
                                    f[6].parse().unwrap(),
                                    f[7].parse().unwrap(),
                                    f[8].parse().unwrap(),
                                    f[9].parse().unwrap(),
                                    f[10].parse().unwrap(),
                                ),
                                (
                                    f[11].parse().unwrap(),
                                    f[12].parse().unwrap(),
                                    f[13].parse().unwrap(),
                                    f[14].parse().unwrap(),
                                    f[15].parse().unwrap(),
                                ),
                            )
                        })
                        .collect();
                assert_eq!(ordered, expected_order, "{name} order/report fields");
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:283-307,870-905;
    // core/blast_traceback.c:198-250,717-719,1481-1501;
    // core/blast_kappa.c:3525-3736:
    // generated masked and unmasked translated chunk boundaries use the same
    // local Stage C search then ordinary or composition Stage D result path.
    #[test]
    fn four_generated_stage_c_chunk_families_match_ncbi_stage_d() {
        let croot = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/extended_chunks_20260925/run_20260925/"
        );
        let source = std::fs::read_to_string(format!("{croot}run_20260923/subjects.fna")).unwrap();
        let insert = source
            .split('>')
            .find(|record| record.starts_with("plus1\n"))
            .unwrap()
            .lines()
            .skip(1)
            .flat_map(str::bytes)
            .collect::<Vec<_>>();
        assert_eq!(insert.len(), 362);
        let query = std::fs::read_to_string(format!("{croot}run_20260923/query.faa"))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(str::bytes)
            .collect::<Vec<_>>();
        for case in [
            "long_chunk_20260924",
            "masked_chunk_boundary_20260924",
            "no_range_middle_20260924",
            "no_range_two_hits_20260924",
        ] {
            let mut subject = if case == "no_range_two_hits_20260924" {
                let first = 4_998_450;
                let second = 10_000_050;
                let mut seq = b"ATG".repeat(first);
                seq.extend(&insert);
                seq.push(b'A');
                seq.extend(b"ATG".repeat(second - first - 121));
                seq.extend(&insert);
                seq.extend(b"ATG".repeat(250));
                seq
            } else if case == "no_range_middle_20260924" {
                let mut seq = b"ATG".repeat(10_000_050);
                seq.extend(&insert);
                seq.extend(b"ATG".repeat(250));
                seq
            } else {
                let mut seq = b"ATG".repeat(4_999_950);
                seq.extend(&insert);
                seq.extend(b"ATG".repeat(250));
                seq
            };
            if case.starts_with("no_range") {
                subject[4_999_000 * 3..10_000_000 * 3].make_ascii_lowercase();
            } else if case == "masked_chunk_boundary_20260924" {
                subject[4_999_900 * 3..5_000_000 * 3].make_ascii_lowercase();
            }
            for mode in [0, 2] {
                let name = format!("{case}.mode{mode}");
                let profile = LocalStageDProfile {
                    seg: None,
                    soft_masking: false,
                    mask_lowercase: case != "long_chunk_20260924",
                    genetic_code: 1,
                    expect_value: 10000.0,
                    max_target_seqs: 500,
                };
                let mut boundary = StageDBoundaryTrace::default();
                let results = run_local_search(
                    &[query.clone()],
                    &[subject.clone()],
                    profile,
                    mode == 2,
                    mode == 2,
                    LocalStageDScoring::default(),
                    Some(&mut boundary),
                )
                .unwrap();
                let calls = std::fs::read_to_string(format!("{root}{name}.calls.trace")).unwrap();
                let kappa = if mode == 2 {
                    Some(std::fs::read_to_string(format!("{root}{name}.kappa.trace")).unwrap())
                } else {
                    None
                };
                assert_local_stage_d_boundaries(&boundary, &calls, kappa.as_deref(), mode == 2);
                let actual: Vec<_> = results[0]
                    .lists()
                    .iter()
                    .flat_map(|list| {
                        list.hsps
                            .hsps
                            .iter()
                            .zip(&list.payloads)
                            .map(|(hsp, payload)| {
                                let frame = i32::from(hsp.hsp.frame);
                                let slen = subject.len() as i32;
                                let (sstart, send) = if frame > 0 {
                                    (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                                } else {
                                    (
                                        slen - 3 * hsp.hsp.s_start + frame + 1,
                                        slen - 3 * hsp.hsp.s_end + frame + 2,
                                    )
                                };
                                (
                                    hsp.hsp.score,
                                    payload.report_num_ident,
                                    payload.report_num_positives,
                                    payload.align_length,
                                    payload.report_mismatches,
                                    payload.gap_letters,
                                    payload.gap_opens,
                                    hsp.hsp.q_start + 1,
                                    hsp.hsp.q_end,
                                    sstart,
                                    send,
                                    hsp.hsp.frame,
                                    hsp.evalue.to_bits(),
                                    payload.bit_score.to_bits(),
                                )
                            })
                    })
                    .collect();
                let report = std::fs::read_to_string(format!("{root}{name}.report.tsv")).unwrap();
                let expected: Vec<_> = report
                    .lines()
                    .map(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        (
                            f[2].parse().unwrap(),
                            f[5].parse().unwrap(),
                            f[6].parse().unwrap(),
                            f[7].parse().unwrap(),
                            f[8].parse().unwrap(),
                            f[9].parse().unwrap(),
                            f[10].parse().unwrap(),
                            f[11].parse().unwrap(),
                            f[12].parse().unwrap(),
                            f[13].parse().unwrap(),
                            f[14].parse().unwrap(),
                            f[15].parse().unwrap(),
                        )
                    })
                    .collect();
                assert_eq!(actual.len(), expected.len(), "{name} retained count");
                for (
                    (
                        score,
                        ident,
                        pos,
                        length,
                        mismatch,
                        gaps,
                        gapopen,
                        q0,
                        q1,
                        s0,
                        s1,
                        frame,
                        e,
                        bits,
                    ),
                    expected,
                ) in actual.into_iter().zip(expected)
                {
                    assert_eq!(
                        (
                            score, ident, pos, length, mismatch, gaps, gapopen, q0, q1, s0, s1,
                            frame
                        ),
                        expected,
                        "{name} report values/order"
                    );
                    let trace = std::fs::read_to_string(format!(
                        "{root}{name}.{}",
                        if mode == 0 {
                            "calls.trace"
                        } else {
                            "kappa.trace"
                        }
                    ))
                    .unwrap();
                    let row = trace
                        .lines()
                        .find(|line| {
                            let f: Vec<_> = line.split('\t').collect();
                            if mode == 0 {
                                f.len() > 14
                                    && f[0] == "D_HSP"
                                    && f[2] == "bits_after"
                                    && f[10].parse::<i32>().unwrap() == score
                                    && f[5].parse::<i8>().unwrap() == frame
                            } else {
                                f.len() > 13
                                    && f[0] == "K_TRACE_HEAP_HSP"
                                    && f[3].parse::<i32>().unwrap() == score
                                    && f[8].parse::<i8>().unwrap() == frame
                            }
                        })
                        .unwrap();
                    let f: Vec<_> = row.split('\t').collect();
                    let (ecol, bcol) = if mode == 0 { (13, 14) } else { (5, 4) };
                    assert_eq!(
                        e,
                        f[ecol].parse::<f64>().unwrap().to_bits(),
                        "{name} E-value"
                    );
                    assert_eq!(
                        bits,
                        f[bcol].parse::<f64>().unwrap().to_bits(),
                        "{name} bits"
                    );
                }
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:747-850,1446-1467;
    // core/blast_traceback.c:198-250,717-719:
    // the long six-frame subject enters the same local posttraceback call with
    // translated statistical length and query-indexed effective search spaces.
    #[test]
    fn long_multi_query_local_stage_d_matches_ncbi() {
        let croot = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/long_subject_20260925/run_20260925/"
        );
        let query_text =
            std::fs::read_to_string(format!("{croot}long_multi_query_20260924/query.faa")).unwrap();
        let mut queries = Vec::new();
        let mut query_names = Vec::new();
        for line in query_text.lines() {
            if let Some(id) = line.strip_prefix('>') {
                query_names.push(id.split_whitespace().next().unwrap().to_owned());
                queries.push(Vec::new());
            } else {
                queries.last_mut().unwrap().extend(line.bytes());
            }
        }
        let source = std::fs::read_to_string(format!("{croot}run_20260923/subjects.fna")).unwrap();
        let plus1: Vec<u8> = source
            .split('>')
            .find_map(|record| {
                let mut lines = record.lines();
                (lines.next()?.split_whitespace().next()? == "plus1")
                    .then(|| lines.flat_map(str::bytes).collect())
            })
            .unwrap();
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(&plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 15_000_962);
        let mut boundary = StageDBoundaryTrace::default();
        let results = run_local_search(
            &queries,
            &[subject.clone()],
            LocalStageDProfile {
                seg: None,
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 500,
            },
            false,
            false,
            LocalStageDScoring::default(),
            Some(&mut boundary),
        )
        .unwrap();
        let mut actual = Vec::new();
        let mut ordered = Vec::new();
        for (qidx, hitlist) in results.iter().enumerate() {
            for list in hitlist.lists() {
                for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                    actual.push((
                        qidx,
                        hsp.hsp.frame,
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                        hsp.hsp.score,
                        payload.num_ident,
                        hsp.evalue.to_bits(),
                        payload.bit_score.to_bits(),
                    ));
                    ordered.push((
                        (query_names[qidx].clone(), "chunk_edge".to_owned()),
                        (
                            hsp.hsp.score,
                            payload.report_num_ident,
                            payload.report_num_positives,
                            payload.align_length,
                            payload.report_mismatches,
                            payload.gap_letters,
                            payload.gap_opens,
                        ),
                        (
                            hsp.hsp.q_start + 1,
                            hsp.hsp.q_end,
                            3 * hsp.hsp.s_start + 1,
                            3 * hsp.hsp.s_end,
                            hsp.hsp.frame,
                        ),
                    ));
                }
            }
        }
        let trace = std::fs::read_to_string(format!("{root}calls.trace")).unwrap();
        assert_local_stage_d_boundaries(&boundary, &trace, None, false);
        let mut expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t") && line.contains("\tbits_after\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap(),
                    f[11].parse().unwrap(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    f[14].parse::<f64>().unwrap().to_bits(),
                )
            })
            .collect();
        actual.sort_unstable();
        expected.sort_unstable();
        assert_eq!(actual, expected, "long subject full-precision final HSPs");
        let expected_order: Vec<_> = std::fs::read_to_string(format!("{root}report.tsv"))
            .unwrap()
            .lines()
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    (f[0].to_owned(), f[1].to_owned()),
                    (
                        f[2].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                    ),
                    (
                        f[11].parse().unwrap(),
                        f[12].parse().unwrap(),
                        f[13].parse().unwrap(),
                        f[14].parse().unwrap(),
                        f[15].parse().unwrap(),
                    ),
                )
            })
            .collect();
        assert_eq!(
            ordered, expected_order,
            "long subject final order/report fields"
        );
    }

    // NCBI c++/src/algo/blast/core/blast_stat.c:2547-2605,4091-4171;
    // core/blast_traceback.c:503-536,717-719: BLOSUM45 gap-14/2 uses
    // its own score block and Karlin E-values with no Spouge gbp.
    #[test]
    fn alternate_blosum45_word2_local_stage_d_matches_ncbi() {
        let input = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
        );
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/alternate_matrix_20260925/run_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
            let mut records: Vec<(String, Vec<u8>)> = Vec::new();
            for line in std::fs::read_to_string(format!("{input}{name}"))
                .unwrap()
                .lines()
            {
                if let Some(id) = line.strip_prefix('>') {
                    records.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                } else {
                    records.last_mut().unwrap().1.extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
        let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
        let mut boundary = StageDBoundaryTrace::default();
        let results = run_local_search(
            &query_sequences,
            &subject_sequences,
            LocalStageDProfile {
                seg: None,
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 500,
            },
            false,
            false,
            LocalStageDScoring {
                matrix: ScoringMatrix::Blosum45,
                gap_open: 14,
                gap_extend: 2,
                word_size: 2,
                threshold: 16,
                window: 60,
                gap_xdrop_bits: 15.0,
                final_xdrop_bits: 25.0,
            },
            Some(&mut boundary),
        )
        .unwrap();
        let mut actual = Vec::new();
        let mut ordered = Vec::new();
        for (qidx, hitlist) in results.iter().enumerate() {
            for list in hitlist.lists() {
                for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                    actual.push((
                        qidx,
                        hsp.hsp.frame,
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                        hsp.hsp.score,
                        payload.num_ident,
                        hsp.evalue.to_bits(),
                        payload.bit_score.to_bits(),
                    ));
                    let frame = i32::from(hsp.hsp.frame);
                    let slen = subjects[list.oid as usize].1.len() as i32;
                    let (sstart, send) = if frame > 0 {
                        (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                    } else {
                        (
                            slen - 3 * hsp.hsp.s_start + frame + 1,
                            slen - 3 * hsp.hsp.s_end + frame + 2,
                        )
                    };
                    ordered.push((
                        (
                            queries[qidx].0.clone(),
                            subjects[list.oid as usize].0.clone(),
                        ),
                        (
                            hsp.hsp.score,
                            payload.report_num_ident,
                            payload.report_num_positives,
                            payload.align_length,
                            payload.report_mismatches,
                            payload.gap_letters,
                            payload.gap_opens,
                        ),
                        (
                            hsp.hsp.q_start + 1,
                            hsp.hsp.q_end,
                            sstart,
                            send,
                            hsp.hsp.frame,
                        ),
                    ));
                }
            }
        }
        let trace = std::fs::read_to_string(format!("{root}calls.trace")).unwrap();
        assert_local_stage_d_boundaries(&boundary, &trace, None, false);
        let mut expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t") && line.contains("\tbits_after\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap(),
                    f[11].parse().unwrap(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    f[14].parse::<f64>().unwrap().to_bits(),
                )
            })
            .collect();
        actual.sort_unstable();
        expected.sort_unstable();
        assert_eq!(actual, expected, "BLOSUM45 full-precision retained HSPs");
        let expected_order: Vec<_> = std::fs::read_to_string(format!("{root}report.tsv"))
            .unwrap()
            .lines()
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    (f[0].to_owned(), f[1].to_owned()),
                    (
                        f[2].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                    ),
                    (
                        f[11].parse().unwrap(),
                        f[12].parse().unwrap(),
                        f[13].parse().unwrap(),
                        f[14].parse().unwrap(),
                        f[15].parse().unwrap(),
                    ),
                )
            })
            .collect();
        assert_eq!(
            ordered, expected_order,
            "BLOSUM45 final order and report fields"
        );
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:198-250,1486-1501;
    // c++/src/algo/blast/core/blast_kappa.c:411-427,3670-3701:
    // compositionBasedStats chooses redo while do_sum_stats independently
    // chooses BLAST_LinkHsps or direct Blast_HSPListGetEvalues.
    #[test]
    fn independent_composition_and_sum_stats_options_match_ncbi() {
        let stage_c = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/option_cross_20260925/run_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
            let mut records: Vec<(String, Vec<u8>)> = Vec::new();
            for line in std::fs::read_to_string(format!("{stage_c}{name}"))
                .unwrap()
                .lines()
            {
                if let Some(id) = line.strip_prefix('>') {
                    records.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                } else {
                    records.last_mut().unwrap().1.extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
        let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
        let seg = SegParams::default();
        for (name, mode) in [("mode0_sumtrue", 0), ("mode2_sumfalse", 2)] {
            let profile = LocalStageDProfile {
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 500,
            };
            let mut boundary = StageDBoundaryTrace::default();
            let result = run_local_search(
                &query_sequences,
                &subject_sequences,
                profile,
                mode == 2,
                mode == 0,
                LocalStageDScoring::default(),
                Some(&mut boundary),
            )
            .unwrap();
            let mut actual = Vec::new();
            let mut ordered = Vec::new();
            for (query_index, hitlist) in result.iter().enumerate() {
                for list in hitlist.lists() {
                    for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                        let frame = i32::from(hsp.hsp.frame);
                        let subject_len = subjects[list.oid as usize].1.len() as i32;
                        let (sstart, send) = if frame > 0 {
                            (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                        } else {
                            (
                                subject_len - 3 * hsp.hsp.s_start + frame + 1,
                                subject_len - 3 * hsp.hsp.s_end + frame + 2,
                            )
                        };
                        actual.push((
                            query_index,
                            hsp.hsp.frame,
                            hsp.hsp.q_start,
                            hsp.hsp.q_end,
                            hsp.hsp.s_start,
                            hsp.hsp.s_end,
                            hsp.hsp.score,
                            hsp.evalue.to_bits(),
                        ));
                        ordered.push((
                            (
                                queries[query_index].0.clone(),
                                subjects[list.oid as usize].0.clone(),
                            ),
                            (
                                hsp.hsp.score,
                                payload.report_num_ident,
                                payload.report_num_positives,
                                payload.align_length,
                                payload.report_mismatches,
                                payload.gap_letters,
                                payload.gap_opens,
                            ),
                            (
                                hsp.hsp.q_start + 1,
                                hsp.hsp.q_end,
                                sstart,
                                send,
                                hsp.hsp.frame,
                            ),
                        ));
                    }
                }
            }
            let trace = std::fs::read_to_string(format!(
                "{root}{name}.{}",
                if mode == 0 {
                    "calls.trace"
                } else {
                    "kappa.trace"
                }
            ))
            .unwrap();
            let calls = std::fs::read_to_string(format!("{root}{name}.calls.trace")).unwrap();
            assert_local_stage_d_boundaries(
                &boundary,
                &calls,
                if mode == 2 { Some(&trace) } else { None },
                mode == 0,
            );
            let mut expected: Vec<_> = trace
                .lines()
                .filter(|line| {
                    if mode == 0 {
                        line.starts_with("D_HSP\t") && line.contains("\tbits_after\t")
                    } else {
                        line.starts_with("K_TRACE_RESULT_HSP\t")
                    }
                })
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    if mode == 0 {
                        (
                            f[4].parse().unwrap(),
                            f[5].parse().unwrap(),
                            f[6].parse().unwrap(),
                            f[7].parse().unwrap(),
                            f[8].parse().unwrap(),
                            f[9].parse().unwrap(),
                            f[10].parse().unwrap(),
                            f[13].parse::<f64>().unwrap().to_bits(),
                        )
                    } else {
                        (
                            usize::from(f[7] == "70"),
                            f[10].parse().unwrap(),
                            f[6].parse().unwrap(),
                            f[7].parse().unwrap(),
                            f[8].parse().unwrap(),
                            f[9].parse().unwrap(),
                            f[4].parse().unwrap(),
                            f[5].parse::<f64>().unwrap().to_bits(),
                        )
                    }
                })
                .collect();
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(
                actual, expected,
                "{name} exact retained score/coords/E-value"
            );
            let expected_order: Vec<_> =
                std::fs::read_to_string(format!("{root}{name}.report.tsv"))
                    .unwrap()
                    .lines()
                    .map(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        (
                            (f[0].to_owned(), f[1].to_owned()),
                            (
                                f[2].parse().unwrap(),
                                f[5].parse().unwrap(),
                                f[6].parse().unwrap(),
                                f[7].parse().unwrap(),
                                f[8].parse().unwrap(),
                                f[9].parse().unwrap(),
                                f[10].parse().unwrap(),
                            ),
                            (
                                f[11].parse().unwrap(),
                                f[12].parse().unwrap(),
                                f[13].parse().unwrap(),
                                f[14].parse().unwrap(),
                                f[15].parse().unwrap(),
                            ),
                        )
                    })
                    .collect();
            assert_eq!(
                ordered, expected_order,
                "{name} final order/report numeric fields"
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_setup.c:614-654;
    // c++/src/algo/blast/core/blast_traceback.c:380-391,583-605,
    // 1481-1501; c++/src/algo/blast/core/blast_kappa.c:3690-3736:
    // hard/soft SEG and query lowercase flags choose the working-query bytes
    // before Stage C, then each composition option follows its source branch.
    #[test]
    fn six_masking_fixtures_match_both_local_stage_d_option_paths() {
        let cases = [
            ("seg_cross_traceback_20260924", true, false, false),
            ("seg_soft_traceback_20260924", true, true, false),
            ("lcase_query_20260924", false, false, true),
            ("lcase_soft_query_20260924", false, true, true),
            ("seg_lcase_overlap_20260924", true, false, true),
            ("seg_lcase_overlap_soft_20260924", true, true, true),
        ];
        let read = |path: &str| -> Vec<u8> {
            std::fs::read_to_string(path)
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('>'))
                .flat_map(str::bytes)
                .collect()
        };
        let seg = SegParams::default();
        for (case, use_seg, soft, lowercase) in cases {
            let input = format!(
                "{}/../docs/evidence/tlosan_stage_c/{case}",
                env!("CARGO_MANIFEST_DIR")
            );
            let queries = vec![read(&format!("{input}/query.faa"))];
            let subjects = vec![read(&format!("{input}/subjects.fna"))];
            for mode in [0, 2] {
                let profile = LocalStageDProfile {
                    seg: if use_seg { Some(&seg) } else { None },
                    soft_masking: soft,
                    mask_lowercase: lowercase,
                    genetic_code: 1,
                    expect_value: 10000.0,
                    max_target_seqs: 500,
                };
                let mut boundary = StageDBoundaryTrace::default();
                let results = run_local_search(
                    &queries,
                    &subjects,
                    profile,
                    mode == 2,
                    mode == 2,
                    LocalStageDScoring::default(),
                    Some(&mut boundary),
                )
                .unwrap();
                assert_eq!(results.len(), 1);
                assert_eq!(results[0].lists().len(), 1, "{case} mode {mode}");
                let list = &results[0].lists()[0];
                assert_eq!(list.hsps.hsps.len(), 1, "{case} mode {mode}");
                let hsp = &list.hsps.hsps[0];
                let payload = &list.payloads[0];
                let root = concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../docs/evidence/tlosan_stage_d/masking_options_20260925/run_20260925/"
                );
                let stem = format!("{case}.mode{mode}");
                let report = std::fs::read_to_string(format!("{root}{stem}.report.tsv")).unwrap();
                let f: Vec<_> = report.trim_end().split('\t').collect();
                assert_eq!(hsp.hsp.score, f[2].parse().unwrap(), "{stem} score");
                assert_eq!(
                    (
                        payload.report_num_ident,
                        payload.report_num_positives,
                        payload.align_length,
                        payload.report_mismatches,
                        payload.gap_letters,
                        payload.gap_opens
                    ),
                    (
                        f[5].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap()
                    ),
                    "{stem} report fields"
                );
                let frame = i32::from(hsp.hsp.frame);
                let subject_len = subjects[0].len() as i32;
                let (sstart, send) = if frame > 0 {
                    (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                } else {
                    (
                        subject_len - 3 * hsp.hsp.s_start + frame + 1,
                        subject_len - 3 * hsp.hsp.s_end + frame + 2,
                    )
                };
                assert_eq!(
                    (
                        hsp.hsp.q_start + 1,
                        hsp.hsp.q_end,
                        sstart,
                        send,
                        hsp.hsp.frame
                    ),
                    (
                        f[11].parse().unwrap(),
                        f[12].parse().unwrap(),
                        f[13].parse().unwrap(),
                        f[14].parse().unwrap(),
                        f[15].parse().unwrap()
                    ),
                    "{stem} coordinates"
                );
                let trace = std::fs::read_to_string(format!(
                    "{root}{stem}.{}",
                    if mode == 0 {
                        "calls.trace"
                    } else {
                        "kappa.trace"
                    }
                ))
                .unwrap();
                let calls = std::fs::read_to_string(format!("{root}{stem}.calls.trace")).unwrap();
                assert_local_stage_d_boundaries(
                    &boundary,
                    &calls,
                    if mode == 2 { Some(&trace) } else { None },
                    mode == 2,
                );
                let row: Vec<_> = trace
                    .lines()
                    .find(|line| {
                        if mode == 0 {
                            line.starts_with("D_HSP\t") && line.contains("\tbits_after\t")
                        } else {
                            line.starts_with("K_TRACE_HEAP_HSP\t")
                        }
                    })
                    .unwrap()
                    .split('\t')
                    .collect();
                let (bits_col, e_col) = if mode == 0 { (14, 13) } else { (4, 5) };
                assert_eq!(
                    payload.bit_score.to_bits(),
                    row[bits_col].parse::<f64>().unwrap().to_bits(),
                    "{stem} full-precision bits"
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    row[e_col].parse::<f64>().unwrap().to_bits(),
                    "{stem} full-precision E-value"
                );
                assert!(!payload.edit_script.is_empty(), "{stem} owned edit script");
            }
        }
    }

    // NCBI c++/src/objects/seqfeat/gc.prt:105-357;
    // c++/src/algo/blast/core/blast_engine.c:762-774,1460-1466;
    // c++/src/algo/blast/core/blast_kappa.c:3690-3736,2494-2515:
    // the selected FindGeneticCode table enters subject translation before
    // preliminary HSPs, redo, heap insertion and final hitlist updates.
    // Each saved comparison uses the CLI-calibrated local API oracle, including
    // code 32, and independently runs the identical input under code 1.
    #[test]
    fn all_27_codes_and_code1_controls_match_local_ncbi_result_traces() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/all_codes_20260925/"
        );
        let fixtures = std::fs::read_to_string(format!("{root}fixtures.tsv")).unwrap();
        let read = |name: &str| -> Vec<u8> {
            std::fs::read_to_string(format!("{root}fixtures/{name}"))
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('>'))
                .flat_map(str::bytes)
                .collect()
        };
        let seg = SegParams::default();
        let mut compared = 0;
        for row in fixtures.lines().skip(1) {
            let code: u8 = row.split('\t').next().unwrap().parse().unwrap();
            let queries = vec![read(&format!("code{code}.faa"))];
            let subjects = vec![read(&format!("code{code}.fna"))];
            for (selected, label) in [(1, "control"), (code, "selected")] {
                let mut boundary = StageDBoundaryTrace::default();
                let results = run_local_search(
                    &queries,
                    &subjects,
                    LocalStageDProfile {
                        seg: Some(&seg),
                        soft_masking: false,
                        mask_lowercase: false,
                        genetic_code: selected,
                        expect_value: 10.0,
                        max_target_seqs: 500,
                    },
                    true,
                    true,
                    LocalStageDScoring::default(),
                    Some(&mut boundary),
                )
                .unwrap();
                let trace = std::fs::read_to_string(format!(
                    "{root}run_20260925/code{code}.{label}.kappa.trace"
                ))
                .unwrap();
                let calls = std::fs::read_to_string(format!(
                    "{root}run_20260925/code{code}.{label}.d.trace"
                ))
                .unwrap();
                assert_local_stage_d_boundaries(&boundary, &calls, Some(&trace), true);
                let expected: Vec<_> = trace
                    .lines()
                    .filter(|line| line.starts_with("K_TRACE_RESULT_HSP\t"))
                    .collect();
                let actual: Vec<_> = results[0]
                    .lists()
                    .iter()
                    .flat_map(|list| list.hsps.hsps.iter().zip(&list.payloads))
                    .collect();
                let report = std::fs::read_to_string(format!(
                    "{root}run_20260925/code{code}.{label}.report.tsv"
                ))
                .unwrap();
                let report_rows: Vec<_> = report.lines().collect();
                let fields = std::fs::read_to_string(format!(
                    "{root}fields_rules_20260925/code{code}.{label}.fields.tsv"
                ))
                .unwrap();
                let field_rows: Vec<_> = fields.lines().collect();
                let rule_trace = std::fs::read_to_string(format!(
                    "{root}fields_rules_20260925/code{code}.{label}.rules.trace"
                ))
                .unwrap();
                let rule_rows: Vec<Vec<_>> = rule_trace
                    .lines()
                    .map(|line| line.split('\t').collect())
                    .collect();
                assert_eq!(
                    actual.len(),
                    expected.len(),
                    "code {code} {label} result count"
                );
                assert_eq!(
                    actual.len(),
                    report_rows.len(),
                    "code {code} {label} report count"
                );
                assert_eq!(
                    actual.len(),
                    field_rows.len(),
                    "code {code} {label} custom report count"
                );
                let heap: Vec<_> = trace
                    .lines()
                    .filter(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
                    .collect();
                assert_eq!(heap.len(), expected.len(), "code {code} {label} heap count");
                let returns: Vec<_> = trace
                    .lines()
                    .filter(|line| line.starts_with("K_TRACE_RETURN\t"))
                    .collect();
                let edits: Vec<_> = trace
                    .lines()
                    .filter(|line| line.starts_with("K_TRACE_EDIT\t"))
                    .collect();
                assert_eq!(
                    returns.len(),
                    edits.len(),
                    "code {code} {label} redo scripts"
                );
                for ((((hsp, payload), expected), reported), custom_reported) in actual
                    .into_iter()
                    .zip(expected)
                    .zip(report_rows)
                    .zip(field_rows)
                {
                    let f: Vec<_> = expected.split('\t').collect();
                    let report_fields: Vec<_> = reported.split('\t').collect();
                    let custom_fields: Vec<_> = custom_reported.split('\t').collect();
                    // NCBI c++/src/objtools/align_format/tabular.cpp:971-1021,
                    // 1090-1095: the diagnostic field list exposes every
                    // report count on the same retained HSP ordering.
                    assert_eq!(custom_fields.len(), 16, "code {code} {label} custom width");
                    assert_eq!(
                        hsp.hsp.score,
                        custom_fields[2].parse().unwrap(),
                        "code {code} {label} custom score"
                    );
                    assert_eq!(
                        payload.report_num_ident,
                        custom_fields[5].parse().unwrap(),
                        "code {code} {label} report identities"
                    );
                    assert_eq!(
                        payload.report_num_positives,
                        custom_fields[6].parse().unwrap(),
                        "code {code} {label} positives"
                    );
                    assert_eq!(
                        payload.align_length,
                        custom_fields[7].parse().unwrap(),
                        "code {code} {label} custom length"
                    );
                    assert_eq!(
                        payload.report_mismatches,
                        custom_fields[8].parse().unwrap(),
                        "code {code} {label} custom mismatches"
                    );
                    assert_eq!(
                        payload.gap_letters,
                        custom_fields[9].parse().unwrap(),
                        "code {code} {label} gap letters"
                    );
                    assert_eq!(
                        payload.gap_opens,
                        custom_fields[10].parse().unwrap(),
                        "code {code} {label} custom gap opens"
                    );
                    // NCBI c++/src/objtools/align_format/tabular.cpp:971-1021,
                    // 1090-1095: report counts are recomputed from the aligned
                    // sequence strings after the result HSPs are retained.
                    assert_eq!(report_fields.len(), 12, "code {code} {label} report width");
                    assert_eq!(
                        payload.align_length,
                        report_fields[3].parse().unwrap(),
                        "code {code} {label} alignment length"
                    );
                    assert_eq!(
                        payload.report_mismatches,
                        report_fields[4].parse().unwrap(),
                        "code {code} {label} mismatches"
                    );
                    assert_eq!(
                        payload.gap_opens,
                        report_fields[5].parse().unwrap(),
                        "code {code} {label} gap opens"
                    );
                    let reported_identity: f64 = report_fields[2].parse().unwrap();
                    let actual_identity =
                        100.0 * payload.report_num_ident as f64 / payload.align_length as f64;
                    assert!(
                        (actual_identity - reported_identity).abs() <= 0.0005,
                        "code {code} {label} percent identity"
                    );
                    let coords = (
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                    );
                    assert_eq!(
                        hsp.hsp.score,
                        f[4].parse().unwrap(),
                        "code {code} {label} score"
                    );
                    assert_eq!(
                        hsp.evalue.to_bits(),
                        f[5].parse::<f64>().unwrap().to_bits(),
                        "code {code} {label} E-value"
                    );
                    assert_eq!(
                        coords,
                        (
                            f[6].parse().unwrap(),
                            f[7].parse().unwrap(),
                            f[8].parse().unwrap(),
                            f[9].parse().unwrap()
                        ),
                        "code {code} {label} coordinates/order"
                    );
                    assert_eq!(
                        (hsp.hsp.frame, hsp.num),
                        (f[10].parse().unwrap(), f[11].parse().unwrap()),
                        "code {code} {label} frame/link count"
                    );
                    let heap_row: Vec<_> = heap
                        .iter()
                        .map(|line| line.split('\t').collect::<Vec<_>>())
                        .find(|row| {
                            row[3].parse::<i32>().unwrap() == hsp.hsp.score
                                && (
                                    row[9].parse::<i32>().unwrap(),
                                    row[10].parse::<i32>().unwrap(),
                                    row[11].parse::<i32>().unwrap(),
                                    row[12].parse::<i32>().unwrap(),
                                ) == coords
                        })
                        .unwrap_or_else(|| panic!("code {code} {label} missing heap HSP"));
                    assert_eq!(
                        payload.bit_score.to_bits(),
                        heap_row[4].parse::<f64>().unwrap().to_bits(),
                        "code {code} {label} bits"
                    );
                    assert_eq!(
                        payload.num_ident,
                        heap_row[6].parse().unwrap(),
                        "code {code} {label} identities"
                    );
                    let matching_scripts: Vec<_> = returns
                        .iter()
                        .zip(&edits)
                        .filter(|(returned, edit)| {
                            let r: Vec<_> = returned.split('\t').collect();
                            let e: Vec<_> = edit.split('\t').collect();
                            assert_eq!(r[2], e[1], "code {code} {label} redo event identity");
                            if e[2] == "-1"
                                || (
                                    r[5].parse::<i32>().unwrap(),
                                    r[6].parse::<i32>().unwrap(),
                                    r[7].parse::<i32>().unwrap(),
                                    r[8].parse::<i32>().unwrap(),
                                ) != coords
                            {
                                return false;
                            }
                            let script: Vec<_> = e[3..]
                                .iter()
                                .map(|op| {
                                    let (kind, len) = op.split_once(':').unwrap();
                                    let len = len.parse().unwrap();
                                    match kind {
                                        "3" => crate::common::GapEditOp::Sub(len),
                                        "0" => crate::common::GapEditOp::Del(len),
                                        "6" => crate::common::GapEditOp::Ins(len),
                                        _ => panic!("unknown NCBI edit operation"),
                                    }
                                })
                                .collect();
                            payload.edit_script == script
                        })
                        .collect();
                    assert_eq!(
                        matching_scripts.len(),
                        1,
                        "code {code} {label} unique redo script owner"
                    );
                    let returned: Vec<_> = matching_scripts[0].0.split('\t').collect();
                    // NCBI c++/src/algo/blast/composition_adjustment/redo_alignment.c:103-124;
                    // core/blast_kappa.c:325-342: BlastCompo_AlignmentNew stores
                    // the raw matrix rule, then conversion moves it with the HSP.
                    let rules: Vec<_> = rule_rows
                        .iter()
                        .filter(|row| {
                            row[1] == returned[4]
                                && row[4] == returned[5]
                                && row[5] == returned[6]
                                && row[6] == returned[7]
                                && row[7] == returned[8]
                                && row[8].parse::<i32>().unwrap() == i32::from(hsp.hsp.frame)
                        })
                        .collect();
                    assert_eq!(rules.len(), 1, "code {code} {label} unique raw matrix rule");
                    assert_eq!(
                        payload.matrix_adjust_rule as i32,
                        rules[0][2].parse().unwrap(),
                        "code {code} {label} matrix rule owner"
                    );
                }
                compared += 1;
            }
        }
        assert_eq!(compared, 54);
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:305-358,3687-3736;
    // c++/src/algo/blast/core/blast_hits.c:3243-3297:
    // A broad 112-target local run retains every significant redone HSP,
    // including both distinct scripts in the two-HSP OID 5 list.
    #[test]
    fn broad_112_subject_results_keep_every_ncbi_report_payload() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<Vec<u8>> {
            let mut records: Vec<Vec<u8>> = Vec::new();
            for line in std::fs::read_to_string(format!("{root}{name}"))
                .unwrap()
                .lines()
            {
                if line.starts_with('>') {
                    records.push(Vec::new());
                } else {
                    records.last_mut().unwrap().extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let seg = SegParams::default();
        let results = run_local_mode2_sum_stats(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 112,
            },
        )
        .unwrap();
        let expected = std::fs::read_to_string(format!(
            "{root}report_payload_20260925/report_fields_all.tsv"
        ))
        .unwrap();
        let actual: Vec<_> = results[0]
            .lists()
            .iter()
            .flat_map(|list| {
                list.hsps
                    .hsps
                    .iter()
                    .zip(&list.payloads)
                    .map(move |(hsp, payload)| {
                        (
                            list.oid,
                            hsp.hsp.score,
                            payload.num_ident,
                            payload.num_positives,
                            payload.align_length,
                            payload.mismatches,
                            payload.gap_letters,
                            payload.gap_opens,
                            hsp.hsp.frame,
                        )
                    })
            })
            .collect();
        let expected: Vec<_> = expected
            .lines()
            .map(|row| {
                let f: Vec<_> = row.split('\t').collect();
                (
                    if let Some(id) = f[1].strip_prefix("weak_") {
                        id.parse::<i32>().unwrap() - 7
                    } else {
                        f[1].strip_prefix("background_")
                            .unwrap()
                            .parse::<i32>()
                            .unwrap()
                            + 12
                    },
                    f[2].parse().unwrap(),
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[13].parse().unwrap(),
                )
            })
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(
            results[0]
                .lists()
                .iter()
                .map(|list| list.payloads.len())
                .sum::<usize>(),
            19
        );
        assert_eq!(
            results[0]
                .lists()
                .iter()
                .find(|list| list.oid == 5)
                .unwrap()
                .payloads
                .len(),
            2
        );
        // NCBI c++/src/algo/blast/core/blast_kappa.c:305-358,3690-3706;
        // c++/src/algo/blast/core/blast_hits.c:1437-1456,3243-3297:
        // GapEditScript *editScript = align->context; align->context = NULL;
        // Blast_HSPInit(..., &editScript, &new_hsp);
        // Blast_HitListUpdate(hitlist, hsp_list);
        // The max-112 run retains each of the 12 alignments recorded by the
        // independent max-2 redo trace, including the one max-2 heap rejected.
        let redo_trace =
            std::fs::read_to_string(format!("{root}result_order_20260925/ncbi.trace")).unwrap();
        let returns: Vec<_> = redo_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_RETURN\t"))
            .collect();
        let edits: Vec<_> = redo_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_EDIT\t"))
            .collect();
        assert_eq!((returns.len(), edits.len()), (12, 12));
        let final_scripts: Vec<_> = results[0]
            .lists()
            .iter()
            .flat_map(|list| {
                list.hsps
                    .hsps
                    .iter()
                    .zip(&list.payloads)
                    .map(|(hsp, payload)| {
                        (
                            hsp.hsp.q_start,
                            hsp.hsp.q_end,
                            hsp.hsp.s_start,
                            hsp.hsp.s_end,
                            payload.edit_script.as_slice(),
                        )
                    })
            })
            .collect();
        for (returned, edit) in returns.iter().zip(edits) {
            let r: Vec<_> = returned.split('\t').collect();
            let e: Vec<_> = edit.split('\t').collect();
            assert_eq!(r[2], e[1], "redo event identity");
            let script: Vec<_> = e[3..]
                .iter()
                .map(|op| {
                    let (kind, len) = op.split_once(':').unwrap();
                    let len = len.parse().unwrap();
                    match kind {
                        "3" => crate::common::GapEditOp::Sub(len),
                        "0" => crate::common::GapEditOp::Del(len),
                        "6" => crate::common::GapEditOp::Ins(len),
                        _ => panic!("unexpected NCBI edit operation"),
                    }
                })
                .collect();
            let coords = (
                r[5].parse::<i32>().unwrap(),
                r[6].parse::<i32>().unwrap(),
                r[7].parse::<i32>().unwrap(),
                r[8].parse::<i32>().unwrap(),
            );
            assert!(
                final_scripts
                    .iter()
                    .any(|&(q0, q1, s0, s1, ops)| (q0, q1, s0, s1) == coords && ops == script),
                "NCBI redo event {} lost its owned script",
                r[2]
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:870-905;
    // c++/src/algo/blast/core/blast_kappa.c:3525-3736,2494-2515;
    // c++/src/algo/blast/core/blast_traceback.c:1763-1776:
    // These local -subject traces cover thirteen subjects and hard SEG with
    // the same setup, Kappa and result-post-pipe order.
    #[test]
    fn thirteen_subject_and_hard_seg_local_results_match_exact_ncbi_heap_fields() {
        for (case, expected_count) in [("run_20260923", 11), ("seg_hard_query_20260924", 1)] {
            let input_root = format!(
                "{}/../docs/evidence/tlosan_stage_c/{case}",
                env!("CARGO_MANIFEST_DIR")
            );
            let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
                let mut records: Vec<(String, Vec<u8>)> = Vec::new();
                for line in std::fs::read_to_string(format!("{input_root}/{name}"))
                    .unwrap()
                    .lines()
                {
                    if let Some(id) = line.strip_prefix('>') {
                        records
                            .push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                    } else {
                        records.last_mut().unwrap().1.extend(line.bytes());
                    }
                }
                records
            };
            let query = read_fasta("query.faa");
            let subjects = read_fasta("subjects.fna");
            let query_sequences: Vec<_> = query.iter().map(|(_, seq)| seq.clone()).collect();
            let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
            let seg = SegParams::default();
            let results = run_local_mode2_sum_stats(
                &query_sequences,
                &subject_sequences,
                LocalStageDProfile {
                    seg: Some(&seg),
                    soft_masking: false,
                    mask_lowercase: false,
                    genetic_code: 1,
                    expect_value: 10.0,
                    max_target_seqs: 500,
                },
            )
            .unwrap();
            let trace_path = format!(
                "{}/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/{case}_default.tsv",
                env!("CARGO_MANIFEST_DIR")
            );
            let trace = std::fs::read_to_string(trace_path).unwrap();
            let expected: std::collections::HashMap<(i32, usize), Vec<&str>> = trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    ((f[1].parse().unwrap(), f[2].parse().unwrap()), f)
                })
                .collect();
            assert_eq!(expected.len(), expected_count);
            let lists = results[0].lists();
            let actual_count: usize = lists.iter().map(|list| list.hsps.hsps.len()).sum();
            assert_eq!(actual_count, expected_count, "{case} HSP count");
            for list in lists {
                for (index, (hsp, payload)) in list.hsps.hsps.iter().zip(&list.payloads).enumerate()
                {
                    let row = &expected[&(list.oid, index)];
                    assert_eq!(hsp.hsp.score, row[3].parse().unwrap(), "{case} score");
                    assert_eq!(
                        payload.bit_score.to_bits(),
                        row[4].parse::<f64>().unwrap().to_bits(),
                        "{case} bits"
                    );
                    assert_eq!(
                        hsp.evalue.to_bits(),
                        row[5].parse::<f64>().unwrap().to_bits(),
                        "{case} E-value"
                    );
                    assert_eq!(
                        payload.num_ident,
                        row[6].parse().unwrap(),
                        "{case} identity"
                    );
                    assert_eq!(hsp.hsp.frame, row[8].parse().unwrap(), "{case} frame");
                    assert_eq!(
                        (
                            hsp.hsp.q_start,
                            hsp.hsp.q_end,
                            hsp.hsp.s_start,
                            hsp.hsp.s_end
                        ),
                        (
                            row[9].parse().unwrap(),
                            row[10].parse().unwrap(),
                            row[11].parse().unwrap(),
                            row[12].parse().unwrap()
                        ),
                        "{case} coordinates"
                    );
                }
            }
            let output_path = format!(
                "{}/../docs/evidence/tlosan_stage_d/run_20260924/{case}_default.out",
                env!("CARGO_MANIFEST_DIR")
            );
            let output = std::fs::read_to_string(output_path).unwrap();
            let actual_order: Vec<_> = lists
                .iter()
                .map(|list| subjects[list.oid as usize].0.as_str())
                .collect();
            let expected_order: Vec<_> = output
                .lines()
                .map(|line| line.split('\t').nth(1).unwrap())
                .collect();
            assert_eq!(actual_order, expected_order, "{case} result order");
        }
    }

    // NCBI c++/src/algo/blast/core/blast_hits.c:1437-1456,3243-3297;
    // c++/src/algo/blast/core/blast_kappa.c:305-358,3687-3736:
    // Natural cross-frame sum-statistics HSPs leave Kappa in score order
    // 238,189,135 but hitlist replacement E-value-sorts them 238,135,189.
    #[test]
    fn natural_evalue_sort_moves_the_owned_kappa_payload() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/natural_positive_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<Vec<u8>> {
            let mut records: Vec<Vec<u8>> = Vec::new();
            for line in std::fs::read_to_string(format!("{root}{name}"))
                .unwrap()
                .lines()
            {
                if line.starts_with('>') {
                    records.push(Vec::new());
                } else {
                    records.last_mut().unwrap().extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("evalue_subjects.fna");
        let results = run_local_mode2_sum_stats(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: None,
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 1,
            },
        )
        .unwrap();
        let lists = results[0].lists();
        assert_eq!(lists.len(), 1);
        assert_eq!(lists[0].oid, 0);
        assert_eq!(
            lists[0]
                .hsps
                .hsps
                .iter()
                .map(|h| h.hsp.score)
                .collect::<Vec<_>>(),
            [238, 135, 189]
        );
        let trace =
            std::fs::read_to_string(format!("{root}run_20260925/evalue_sort.kappa.trace")).unwrap();
        let heap: std::collections::HashMap<i32, Vec<&str>> = trace
            .lines()
            .filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t0\t"))
            .map(|row| {
                let f: Vec<_> = row.split('\t').collect();
                (f[3].parse().unwrap(), f)
            })
            .collect();
        for (hsp, payload) in lists[0].hsps.hsps.iter().zip(&lists[0].payloads) {
            let row = &heap[&hsp.hsp.score];
            assert_eq!(
                hsp.evalue.to_bits(),
                row[5].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(
                payload.bit_score.to_bits(),
                row[4].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(payload.num_ident, row[6].parse().unwrap());
            assert_eq!(
                (
                    hsp.hsp.q_start,
                    hsp.hsp.q_end,
                    hsp.hsp.s_start,
                    hsp.hsp.s_end
                ),
                (
                    row[9].parse().unwrap(),
                    row[10].parse().unwrap(),
                    row[11].parse().unwrap(),
                    row[12].parse().unwrap()
                )
            );
            let script_len = match hsp.hsp.score {
                238 => 45,
                135 | 189 => 35,
                _ => unreachable!(),
            };
            assert_eq!(
                payload.edit_script,
                vec![crate::common::GapEditOp::Sub(script_len)]
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:224-273,3658-3689:
    // if (hsp_list->hspcnt > 1)
    //     s_HitlistReapContained(hsp_list->hsp_array, &hsp_list->hspcnt);
    // The pinned local fixture naturally redos two same-frame alignments;
    // the smaller one is contained, leaving one postredo link input.
    #[test]
    fn natural_postredo_containment_keeps_surviving_script_and_report_fields() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/natural_positive_20260925/"
        );
        let read = |name: &str| -> Vec<u8> {
            std::fs::read_to_string(format!("{root}{name}"))
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('>'))
                .flat_map(str::bytes)
                .collect()
        };
        let queries = vec![read("containment_query.faa")];
        let subjects = vec![read("containment_subject.fna")];
        let results = run_local_mode2_sum_stats(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: None,
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 1,
            },
        )
        .unwrap();
        let lists = results[0].lists();
        assert_eq!(lists.len(), 1);
        assert_eq!(lists[0].hsps.hsps.len(), 1);
        let hsp = &lists[0].hsps.hsps[0];
        let payload = &lists[0].payloads[0];
        let trace =
            std::fs::read_to_string(format!("{root}containment_run_20260925/kappa.trace")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .find(|row| row.starts_with("K_TRACE_HEAP_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(hsp.hsp.score, expected[3].parse().unwrap());
        assert_eq!(
            payload.bit_score.to_bits(),
            expected[4].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            hsp.evalue.to_bits(),
            expected[5].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(payload.num_ident, expected[6].parse().unwrap());
        assert_eq!(
            (
                hsp.hsp.frame,
                hsp.hsp.q_start,
                hsp.hsp.q_end,
                hsp.hsp.s_start,
                hsp.hsp.s_end
            ),
            (
                expected[8].parse().unwrap(),
                expected[9].parse().unwrap(),
                expected[10].parse().unwrap(),
                expected[11].parse().unwrap(),
                expected[12].parse().unwrap()
            )
        );
        assert_eq!(
            payload.edit_script,
            vec![
                crate::common::GapEditOp::Sub(18),
                crate::common::GapEditOp::Ins(5),
                crate::common::GapEditOp::Sub(97),
            ]
        );
        assert_eq!(payload.matrix_adjust_rule as i32, 4);
        let report =
            std::fs::read_to_string(format!("{root}containment_run_20260925/report.tsv")).unwrap();
        let f: Vec<_> = report.trim_end().split('\t').collect();
        assert_eq!(
            (
                payload.num_ident,
                payload.num_positives,
                payload.align_length,
                payload.mismatches,
                payload.gap_letters,
                payload.gap_opens
            ),
            (
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
                f[10].parse().unwrap()
            )
        );
    }

    // NCBI c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391;
    // c++/src/algo/blast/core/blast_kappa.c:3719-3755:
    // Two naturally generated weak subjects make the second Kappa record
    // replace the first at heap threshold one, both above inclusion 0.002.
    #[test]
    fn natural_composition_heap_replacement_keeps_new_hsp_payloads() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/natural_positive_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<Vec<u8>> {
            let mut records: Vec<Vec<u8>> = Vec::new();
            for line in std::fs::read_to_string(format!("{root}{name}"))
                .unwrap()
                .lines()
            {
                if line.starts_with('>') {
                    records.push(Vec::new());
                } else {
                    records.last_mut().unwrap().extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("heap_subjects.fna");
        let results = run_local_mode2_sum_stats(
            &queries,
            &subjects,
            LocalStageDProfile {
                seg: None,
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 1,
            },
        )
        .unwrap();
        let lists = results[0].lists();
        assert_eq!(lists.len(), 1);
        assert_eq!(lists[0].oid, 0);
        let trace =
            std::fs::read_to_string(format!("{root}run_20260925/heap_replacement.kappa.trace"))
                .unwrap();
        assert!(trace
            .lines()
            .any(|row| row.starts_with("K_TRACE_HEAP_INSERT_RETURN\t0\t") && row.ends_with("\t1")));
        let rows: Vec<_> = trace
            .lines()
            .filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t0\t"))
            .map(|row| row.split('\t').collect::<Vec<_>>())
            .collect();
        assert_eq!(lists[0].hsps.hsps.len(), rows.len());
        for ((hsp, payload), row) in lists[0].hsps.hsps.iter().zip(&lists[0].payloads).zip(rows) {
            assert_eq!(hsp.hsp.score, row[3].parse().unwrap());
            assert_eq!(
                hsp.evalue.to_bits(),
                row[5].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(
                payload.bit_score.to_bits(),
                row[4].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(payload.num_ident, row[6].parse().unwrap());
            assert!(!payload.edit_script.is_empty());
        }
    }
}

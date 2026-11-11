//! The values of a batch's subject search that `search_query_batch` sets before its subjects: the
//! progress bar and the hit list sizes, the diagnostic flags, the scoring and X-dropoffs of the
//! extensions, the schedule of the subjects on the search pool (`subject_schedule`), and the
//! hit-saving options. `process_subject` and the tail of `search_query_batch` read them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1607
//! (`BLAST_PreliminarySearchEngine`, the set-up before its subject loop);
//! c++/src/algo/blast/core/blast_hits.c:43-71 (`GetPrelimHitlistSize`);
//! c++/src/algo/blast/core/blast_parameters.c:422-490 (`BlastExtensionParametersNew`, the gapped
//! X-dropoffs) and :158-277 (`BlastInitialWordParametersNew`, the ungapped X-dropoff and the
//! score table of the word extension); c++/src/algo/blast/core/blast_stat.c:1052-1130
//! (`BlastScoreBlkNuclMatrixCreate`); c++/src/algo/blast/api/prelim_stage.cpp:141-203
//! (`CBlastPrelimSearch::x_LaunchMultiThreadedSearch`, the threads of the subjects).
//!
//! This module is a child of `run` for the types and modules it keeps private (`SubjectSource`,
//! `x_few_subj`, `x_pair_seed`, `x_qgroup_par`). The functions are statements of
//! `search_query_batch` moved here in their order (the borrows that `search_query_batch` makes for
//! `process_subject` between them stay there); what they read is a struct of references or `Copy`
//! values, and what they make is what `search_query_batch` then reads under the same names.

use super::{
    blastn_trace, build_blastna_matrix, build_nucl_score_table, get_prelim_hitlist_size,
    x_few_subj, x_pair_seed, x_qgroup_par, BlastnArgs, ContextKarlin, SequenceData, SubjectSource,
    NCBIMATH_LN2,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::coordination::TaskConfig;
use crate::algorithm::blastn::lookup::TwoStageLookup;
use crate::algorithm::blastn::scoring::gap_x_dropoffs;
use crate::utils::threading::SearchPool;
use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};

/// What `subject_search_env` reads.
pub(super) struct SubjectEnvInputs<'a> {
    pub(super) args: &'a BlastnArgs,
    pub(super) config: &'a TaskConfig,
    pub(super) seq_data: &'a SequenceData<'a>,
    pub(super) subjects: SubjectSource<'a>,
    /// The subjects as tasks of the pool, before LOSAT_X_BNFEWSUBJ's rule.
    pub(super) subjects_parallel: bool,
    pub(super) parallel_pool: &'a SearchPool<'a>,
    pub(super) num_threads: usize,
    /// The batch's two-stage lookup table (`None` for the other tables and a split batch).
    pub(super) two_stage_lookup_ref: Option<&'a TwoStageLookup>,
    pub(super) diag_container_query_length: usize,
    pub(super) window_size: usize,
    pub(super) small_na_compressed_query: &'a [u8],
    /// The Karlin-Altschul blocks of the contexts (`context_statistics`).
    pub(super) context_karlin: &'a [Option<ContextKarlin>],
    pub(super) search_karlin: &'a [ContextKarlin],
}

/// The values of the batch's subject search that `process_subject` and the tail of
/// `search_query_batch` read, which `search_query_batch` destructures into locals of the same
/// names; `process_subject` is to take this struct as its argument when it becomes a function
/// (SW7 bn-11).
pub(super) struct SubjectSearchEnv {
    pub(super) progress_bar: Option<ProgressBar>,
    pub(super) hitlist_size: usize,
    pub(super) prelim_hitlist_size: usize,
    // The diagnostic flags (`BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`, the trace,
    // `BlastnArgs::verbose`, `BLEMIR_DEBUG_WINDOW`).
    pub(super) debug_mode: bool,
    pub(super) blastn_debug: bool,
    pub(super) blastn_trace_enabled: bool,
    pub(super) verbose: bool,
    pub(super) debug_window: Option<(usize, usize, usize, usize)>,
    pub(super) lcase_masking: bool,
    // The task's values (`TaskConfig`) and the X-dropoffs in raw score.
    pub(super) effective_word_size: usize,
    pub(super) min_ungapped_score: i32,
    pub(super) use_dp: bool,
    pub(super) use_direct_lookup: bool,
    pub(super) reward: i32,
    pub(super) penalty: i32,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) x_drop_gapped: i32,
    pub(super) x_drop_final: i32,
    pub(super) scan_range: usize,
    pub(super) min_diag_separation: i32,
    // The schedule of the subjects (`subject_schedule`).
    pub(super) x_strips: bool,
    pub(super) x_kind: x_few_subj::Kind,
    pub(super) use_parallel: bool,
    // The database size, the e-value threshold, `-subject_besthit` and LOSAT_X_QGROUPPAR.
    pub(super) db_len_total: usize,
    pub(super) db_num_seqs: usize,
    pub(super) evalue_threshold: f64,
    pub(super) subject_besthit: bool,
    pub(super) x_qgroup_batch: bool,
    // The hit-saving options, the score matrix and table, and the ungapped X-dropoff of each
    // context.
    pub(super) percent_identity: f64,
    pub(super) min_hit_length: usize,
    pub(super) score_matrix: BlastnaMatrix,
    pub(super) nucl_score_table: [i32; 256],
    pub(super) x_dropoff_init: Vec<i32>,
}

/// The values of the batch's subject search (`SubjectSearchEnv`); the error of a preliminary hit
/// list size that overflows.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1607
/// (`BLAST_PreliminarySearchEngine`, before its subject loop),
/// c++/src/algo/blast/core/blast_hits.c:43-71 (`GetPrelimHitlistSize`),
/// c++/src/algo/blast/core/blast_parameters.c:422-490 (`BlastExtensionParametersNew`) and
/// :158-277 (`BlastInitialWordParametersNew`).
pub(super) fn subject_search_env<const UNGAPPED: bool>(
    inputs: &SubjectEnvInputs<'_>,
) -> Result<SubjectSearchEnv> {
    let &SubjectEnvInputs {
        args,
        config,
        seq_data,
        subjects,
        subjects_parallel,
        parallel_pool,
        num_threads,
        two_stage_lookup_ref,
        diag_container_query_length,
        window_size,
        small_na_compressed_query,
        context_karlin,
        search_karlin,
    } = inputs;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    let subject_count = seq_data.db_num_seqs as u64;
    // The bar is the command line's (`BlastnArgs::progress_bar`): a search from memory
    // (`run_local`: the library API, the web entries) writes nothing to standard error.
    let progress_bar = if subjects_parallel && args.progress_bar {
        let bar = ProgressBar::new(subject_count);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
                .unwrap(),
        );
        Some(bar)
    } else {
        None
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    // ```c
    // if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //    hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    // }
    // m_NumDescriptions = hitlist_size;
    // m_NumAlignments = hitlist_size;
    // ```
    let hitlist_size = args.report_counts.hitlist_size;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70 (GetPrelimHitlistSize)
    // ```c
    // Int4
    // GetPrelimHitlistSize(Int4 hitlist_size, Int4 compositionBasedStats,
    //                      Boolean gapped_calculation)
    // {
    //     ...
    //     else if (gapped_calculation) {
    //          prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10),
    //                                   prelim_hitlist_size + 50);
    //     }
    //     return prelim_hitlist_size;
    // }
    // ```
    let prelim_hitlist_size =
        usize::try_from(get_prelim_hitlist_size(hitlist_size, false, !UNGAPPED))
            .map_err(|_| anyhow::anyhow!("the preliminary hit list size overflows"))?;

    // Debug mode: set BLEMIR_DEBUG=1 to enable, BLEMIR_DEBUG_WINDOW="q_start-q_end,s_start-s_end" to focus on a region
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    // BLASTN-specific debug mode: set LOSAT_DEBUG_BLASTN=1 to enable detailed hit loss diagnostics
    let blastn_debug = std::env::var("LOSAT_DEBUG_BLASTN").is_ok();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:350-692
    // ```c
    // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    // ...
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
    // ...
    // Blast_IntervalTreeReset(tree);
    // ```
    let blastn_trace_enabled = blastn_trace::enabled();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3267-3282
    // ```c
    // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
    //                  true);
    // ...
    // m_DebugOutput = static_cast<bool>(args["verbose"]);
    // ```
    let verbose = args.verbose;
    let debug_window: Option<(usize, usize, usize, usize)> =
        std::env::var("BLEMIR_DEBUG_WINDOW").ok().and_then(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                let q_parts: Vec<usize> =
                    parts[0].split('-').filter_map(|x| x.parse().ok()).collect();
                let s_parts: Vec<usize> =
                    parts[1].split('-').filter_map(|x| x.parse().ok()).collect();
                if q_parts.len() == 2 && s_parts.len() == 2 {
                    return Some((q_parts[0], q_parts[1], s_parts[0], s_parts[1]));
                }
            }
            None
        });

    // NCBI BLAST does NOT use mask_array/mask_hash for diagonal suppression
    // REMOVED: disable_mask variable (no longer needed)

    if debug_mode {
        // Build marker to verify correct code is running
        eprintln!(
            "[DEBUG] BLEMIR build: 2024-12-24-v7 (adaptive banding with MAX_WINDOW_SIZE=50000)"
        );
        eprintln!(
            "[DEBUG] Task: {}, Scoring: reward={}, penalty={}, gap_open={}, gap_extend={}",
            args.task, config.reward, config.penalty, config.gap_open, config.gap_extend
        );
        if let Some((q_start, q_end, s_start, s_end)) = debug_window {
            eprintln!(
                "[DEBUG] Focusing on window: query {}-{}, subject {}-{}",
                q_start, q_end, s_start, s_end
            );
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
    // ```c
    // const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking)
    //     ? bool(args[kArgUseLCaseMasking])
    //     : kDfltArgUseLCaseMasking;
    // m_Scope = ReadSequencesToBlast(... use_lcase_masks, subjects, ...);
    // ```
    let lcase_masking = args.lcase_masking;

    // Capture config values for use in closure
    let effective_word_size = config.effective_word_size;
    let min_ungapped_score = config.min_ungapped_score;
    let use_dp = config.use_dp;
    let use_direct_lookup = config.use_direct_lookup;
    let reward = config.reward;
    let penalty = config.penalty;
    let gap_open = config.gap_open;
    let gap_extend = config.gap_extend;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:455-463
    // ```c
    //       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
    //       params->gap_x_dropoff = (Int4)
    //           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    //     ...
    //       params->gap_x_dropoff_final = (Int4)
    //           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
    // ```
    let (x_drop_gapped, x_drop_final) =
        gap_x_dropoffs(context_karlin, config.x_drop_gapped, config.x_drop_final);
    let scan_range = config.scan_range; // For off-diagonal hit detection
    let min_diag_separation = config.min_diag_separation; // For MB_HSP_CLOSE containment check

    let SubjectSchedule {
        x_strips,
        x_kind,
        use_parallel,
    } = subject_schedule(&ScheduleInputs {
        args,
        seq_data,
        subjects,
        subjects_parallel,
        num_threads,
        two_stage_lookup_ref,
        diag_container_query_length,
        window_size,
        small_na_compressed_query,
        scan_range,
        lcase_masking,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        debug_window,
    });
    let db_len_total = seq_data.db_len_total;
    let db_num_seqs = seq_data.db_num_seqs;
    let evalue_threshold = crate::algorithm::blastn::coordination::determine_evalue(args);
    let subject_besthit = args.subject_besthit;
    // EXPERIMENT (LOSAT_X_QGROUPPAR): whether the gapped stages of a subject searched on this
    // thread may run per query on the pool (x_qgroup_par.rs). Not on the diagnostic paths, which
    // print in the loop's order, and not while the pool searches subjects in parallel
    // (`use_parallel`; several subjects one at a time under LOSAT_X_BNFEWSUBJ take it).
    let x_qgroup_batch = !UNGAPPED
        && x_qgroup_par::on()
        && !use_parallel
        && parallel_pool.enabled()
        && num_threads > 1
        && !verbose
        && !debug_mode
        && !blastn_debug
        && !blastn_trace_enabled
        && debug_window.is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS").is_none()
        && std::env::var_os("LOSAT_DEBUG_COORDS_START").is_none();
    if x_qgroup_par::on() && !x_qgroup_batch {
        crate::utils::xreach::x_reach!("LOSAT_X_QGROUPPAR", false);
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    // ```c
    // return ((hsp->num_ident * 100.0 <
    //         align_length * hit_options->percent_identity) ||
    //         align_length < hit_options->min_hit_length) ;
    // ```
    let percent_identity = args.percent_identity_value();
    let min_hit_length = args.min_hit_length;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1052-1130 (BlastScoreBlkNuclMatrixCreate)
    let score_matrix = build_blastna_matrix(reward, penalty);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:236-259
    // ```c
    // for (i = 0; i < 256; i++) {
    //    Int4 score = 0;
    //    if (i & 3) score += penalty; else score += reward;
    //    if ((i >> 2) & 3) score += penalty; else score += reward;
    //    if ((i >> 4) & 3) score += penalty; else score += reward;
    //    if (i >> 6) score += penalty; else score += reward;
    //    table[i] = score;
    // }
    // ```
    let nucl_score_table = build_nucl_score_table(reward, penalty);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:218-221
    // ```c
    // p->cutoffs[context].x_dropoff_init =
    //     (Int4)(sbp->scale_factor *
    //            ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // ```
    // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.cpp:49-54
    // m_InitWordOpts.Reset((BlastInitialWordOptions*)calloc(1, sizeof(BlastInitialWordOptions)));
    // NCBI reference: c++/src/algo/blast/api/blast_nucl_options.cpp:163-174
    // SetInitialWordOptionsDefaults() { SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL); ... }
    // SetMBInitialWordOptionsDefaults() { SetWindowSize(BLAST_WINDOW_SIZE_NUCL); }
    // NCBI reference: c++/src/algo/blast/api/disc_nucl_options.cpp:66-73
    // SetMBInitialWordOptionsDefaults() { SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_NUCL); ... }
    // Traditional megablast retains the zero-initialized X-drop option (0 bits, so 0 here):
    // its raw X-drop is the subject's word cutoff in ParametersUpdate
    // (`batch_setup::batch_cutoffs`). `-xdrop_ungap` replaces the task's bits (`configure_task`).
    let x_drop_ungapped_bits = config.x_drop_ungapped;
    let x_dropoff_init: Vec<i32> = search_karlin
        .iter()
        .map(|blocks| {
            if x_drop_ungapped_bits == 0.0 {
                0
            } else {
                ((x_drop_ungapped_bits * NCBIMATH_LN2) / blocks.ungapped.lambda).ceil() as i32
            }
        })
        .collect();
    Ok(SubjectSearchEnv {
        progress_bar,
        hitlist_size,
        prelim_hitlist_size,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        verbose,
        debug_window,
        lcase_masking,
        effective_word_size,
        min_ungapped_score,
        use_dp,
        use_direct_lookup,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_gapped,
        x_drop_final,
        scan_range,
        min_diag_separation,
        x_strips,
        x_kind,
        use_parallel,
        db_len_total,
        db_num_seqs,
        evalue_threshold,
        subject_besthit,
        x_qgroup_batch,
        percent_identity,
        min_hit_length,
        score_matrix,
        nucl_score_table,
        x_dropoff_init,
    })
}

/// What `subject_schedule` reads.
struct ScheduleInputs<'a> {
    args: &'a BlastnArgs,
    seq_data: &'a SequenceData<'a>,
    subjects: SubjectSource<'a>,
    subjects_parallel: bool,
    num_threads: usize,
    two_stage_lookup_ref: Option<&'a TwoStageLookup>,
    diag_container_query_length: usize,
    window_size: usize,
    small_na_compressed_query: &'a [u8],
    scan_range: usize,
    lcase_masking: bool,
    /// The diagnostic flags: a search with a diagnostic path on has no strips.
    debug_mode: bool,
    blastn_debug: bool,
    blastn_trace_enabled: bool,
    debug_window: Option<(usize, usize, usize, usize)>,
}

/// How the batch searches its subjects on the pool: whether the seed stage scans in strips
/// (`x_strips`) with the strips' kind of table (`x_kind`), and whether the subjects are tasks of
/// the pool (`use_parallel`).
struct SubjectSchedule {
    x_strips: bool,
    x_kind: x_few_subj::Kind,
    use_parallel: bool,
}

/// The schedule of the batch's subjects on the search pool: LOSAT_X_BNFEWSUBJ's rule, its
/// `x_reach!` record, and the stage record of the threads (`report_stage`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:141-203
/// (`CBlastPrelimSearch::x_LaunchMultiThreadedSearch`, the threads of the preliminary search).
fn subject_schedule(inputs: &ScheduleInputs<'_>) -> SubjectSchedule {
    let &ScheduleInputs {
        args,
        seq_data,
        subjects,
        subjects_parallel,
        num_threads,
        two_stage_lookup_ref,
        diag_container_query_length,
        window_size,
        small_na_compressed_query,
        scan_range,
        lcase_masking,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        debug_window,
    } = inputs;
    // EXPERIMENT (LOSAT_X_BNFEWSUBJ, x_few_subj.rs): whether the seed stage's strips would take an
    // unmasked chunk of this `-subject` search that one subject's owner searches on this thread's
    // pool, and with it whether the search takes the shape of a one-subject search (`use_parallel`
    // false: the subjects one at a time). A batch split into query chunks has no lookup table
    // here (its traceback only), so it keeps its shape; each query chunk's call decides for itself.
    // `-subject` only: a database (`-db`) keeps its shapes (plan v4 line D).
    let x_strips = x_few_subj::on()
        && subjects.prepared().is_some()
        && x_pair_seed::takes_owner_chunks(
            two_stage_lookup_ref,
            diag_container_query_length <= 8000,
            window_size,
            scan_range,
            !small_na_compressed_query.is_empty(),
            debug_mode || blastn_debug || blastn_trace_enabled || debug_window.is_some(),
        );
    let x_kind = x_few_subj::Kind::of(
        two_stage_lookup_ref.is_some_and(|two_stage| two_stage.disc().is_some()),
        &args.task,
    );
    let x_serial_subjects = subjects_parallel
        && x_strips
        && !lcase_masking
        && subjects.prepared().is_some_and(|prepared| {
            let letters = prepared.letters();
            x_few_subj::serial_subjects(
                (0..prepared.records.len()).map(|index| letters.len(index)),
                num_threads,
                x_kind,
            )
        });
    if x_few_subj::on() && subjects_parallel && subjects.prepared().is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNFEWSUBJ", x_serial_subjects);
    }
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let use_parallel = subjects_parallel && !x_serial_subjects;
    crate::utils::threading::report_stage("blastn", "subjects", seq_data.db_num_seqs, use_parallel);
    SubjectSchedule {
        x_strips,
        x_kind,
        use_parallel,
    }
}

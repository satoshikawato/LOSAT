//! The values of a batch's subject search that `search_query_batch` sets before its subjects: the
//! progress bar and the hit list sizes, the diagnostic flags, the scoring and X-dropoffs of the
//! extensions, the schedule of the subjects on the search pool (`subject_schedule`), and the
//! hit-saving options. `SubjectPasses` (`subject_passes.rs`) and `search_query_batch` read them.
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
//! `SubjectPasses` between them stay there); what they read is a struct of references or `Copy`
//! values, and what they make is what `search_query_batch` then reads under the same names.

use super::{
    blastn_trace, build_blastna_matrix, build_nucl_score_table, get_prelim_hitlist_size,
    x_few_subj, x_pair_seed, x_qgroup_par, BlastnArgs, ContextKarlin, SubjectSource,
};
use crate::algorithm::blastn::alignment::gapped::BlastnaMatrix;
use crate::algorithm::blastn::coordination::TaskConfig;
use crate::algorithm::blastn::lookup::TwoStageLookup;
use crate::algorithm::blastn::scoring::{check_gapped_x_dropoffs, gap_x_dropoffs};
use crate::core::blast_util::ncbi_int4_from_double;
use crate::utils::threading::SearchPool;
use anyhow::Result;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressFinish, ProgressStyle};

/// The command line's progress bar over the subjects of a batch's parallel search
/// (`BlastnArgs::progress_bar`), drawn to `target`: the standard error, which indicatif draws only
/// when it is a terminal and, outside Windows, `TERM` is set and is not `dumb`. Every end of the bar
/// clears it (`ProgressFinish::AndClear`): `search_query_batch` finishes it when the traceback pass,
/// the only stage that moves it, ends or fails (`finish_using_style`), and a bar dropped unfinished
/// is cleared too, so nothing of it stays on the terminal.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
/// No NCBI counterpart: NCBI's blastn writes nothing to its standard error while it searches its
/// query batches. The bar writes only to a terminal and is cleared when the batch's search ends,
/// so the terminal then shows what NCBI leaves on it.
pub(super) fn batch_progress_bar(subject_count: u64, target: ProgressDrawTarget) -> ProgressBar {
    ProgressBar::with_draw_target(Some(subject_count), target)
        .with_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
                .unwrap(),
        )
        .with_finish(ProgressFinish::AndClear)
}

/// What `subject_search_env` reads.
pub(super) struct SubjectEnvInputs<'a> {
    pub(super) args: &'a BlastnArgs,
    pub(super) config: &'a TaskConfig,
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

/// The values of the batch's subject search that `SubjectPasses` and
/// `search_query_batch` read, which `search_query_batch` destructures into locals of the same
/// names; `search_subject`'s two stages read the ones they need through
/// `prelim_subject::PrelimBatch` and `traceback_subject::TraceBatch` (or their copies per subject).
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
    // The e-value threshold, `-subject_besthit` and LOSAT_X_QGROUPPAR.
    pub(super) evalue_threshold: f64,
    pub(super) subject_besthit: bool,
    pub(super) x_qgroup_batch: bool,
    // The hit-saving options, the score matrix and table, and the ungapped X-dropoff of each
    // context.
    pub(super) percent_identity: f64,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1409
    let subject_count = subjects.num_subjects() as u64;
    // The bar is the command line's (`BlastnArgs::progress_bar`): a search from memory
    // (`run_local`: the library API, the web entries) writes nothing to standard error.
    let progress_bar = (subjects_parallel && args.progress_bar)
        .then(|| batch_progress_bar(subject_count, ProgressDrawTarget::stderr()));

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    let hitlist_size = args.report_counts.hitlist_size;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70 (GetPrelimHitlistSize)
    let prelim_hitlist_size = get_prelim_hitlist_size(hitlist_size, false, !UNGAPPED);

    // Debug mode: set BLEMIR_DEBUG=1 to enable, BLEMIR_DEBUG_WINDOW="q_start-q_end,s_start-s_end" to focus on a region
    let debug_mode = crate::utils::xdiag::var("BLEMIR_DEBUG").is_ok();
    // BLASTN-specific debug mode: set LOSAT_DEBUG_BLASTN=1 to enable detailed hit loss diagnostics
    let blastn_debug = std::env::var("LOSAT_DEBUG_BLASTN").is_ok();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:350-692
    let blastn_trace_enabled = blastn_trace::enabled();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3267-3282
    let verbose = args.verbose;
    let debug_window: Option<(usize, usize, usize, usize)> =
        crate::utils::xdiag::var("BLEMIR_DEBUG_WINDOW")
            .ok()
            .and_then(|s| {
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-463
    //
    // Under matrix-only scoring the options' values are raw scores, and the final X-drop is not
    // raised to the preliminary one.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:472-482
    // The assignment converts the double option as a C cast (`ncbi_int4_from_double`).
    let (x_drop_gapped, x_drop_final) = if config.matrix_only {
        (
            ncbi_int4_from_double(config.x_drop_gapped),
            ncbi_int4_from_double(config.x_drop_final),
        )
    } else {
        gap_x_dropoffs(context_karlin, config.x_drop_gapped, config.x_drop_final)
    };
    if !UNGAPPED {
        check_gapped_x_dropoffs(
            config.use_dp,
            config.matrix_only,
            (config.x_drop_gapped, config.x_drop_final),
            (x_drop_gapped, x_drop_final),
            (config.reward, config.penalty),
            (config.gap_open, config.gap_extend),
            args.app().name(),
        )?;
    }
    let scan_range = config.scan_range; // For off-diagonal hit detection
    let min_diag_separation = config.min_diag_separation; // For MB_HSP_CLOSE containment check

    let SubjectSchedule {
        x_strips,
        x_kind,
        use_parallel,
    } = subject_schedule(&ScheduleInputs {
        args,
        subjects,
        subjects_parallel,
        num_threads,
        two_stage_lookup_ref,
        diag_container_query_length,
        window_size,
        small_na_compressed_query,
        scan_range,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        debug_window,
    });
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    let percent_identity = args.percent_identity_value();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1052-1130 (BlastScoreBlkNuclMatrixCreate)
    // The matrix and the score table below score with `sbp->reward`/`sbp->penalty`, the engine
    // pair (1/-3 under matrix-only scoring, `coordination::engine_scoring_pair`); the option pair
    // (`reward`, `penalty`) stays for the traceback's re-evaluation (blast_hits.c:627,637).
    let score_matrix = build_blastna_matrix(config.score_reward, config.score_penalty);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:236-259
    let nucl_score_table = build_nucl_score_table(config.score_reward, config.score_penalty);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:218-221
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.cpp:49-54
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:163-174
    // `SetInitialWordOptionsDefaults` sets ungapped X-drop to `BLAST_UNGAPPED_X_DROPOFF_NUCL`,
    // whereas `SetMBInitialWordOptionsDefaults` sets `BLAST_WINDOW_SIZE_NUCL` without setting X-drop.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/disc_nucl_options.cpp:66-73
    // Discontiguous megablast overrides `SetMBInitialWordOptionsDefaults` to set `BLAST_UNGAPPED_X_DROPOFF_NUCL`.
    // Traditional megablast retains the zero-initialized X-drop option (0 bits, so 0 here):
    // its raw X-drop is the subject's word cutoff in ParametersUpdate
    // (`batch_setup::batch_cutoffs`). `-xdrop_ungap` replaces the task's bits (`configure_task`).
    // The conversion is the protein programs' (`core::x_drop::ungapped_x_dropoff`, scale factor
    // 1): 0 bits and every value whose ceiling is 0 give 0, a value past `Int4` gives `INT_MIN`.
    //
    // Under matrix-only scoring the option's value is the raw X-drop (and `batch_cutoffs` keeps
    // a 0 as 0).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:208-217
    let x_drop_ungapped_bits = config.x_drop_ungapped;
    let x_dropoff_init: Vec<i32> = if config.matrix_only {
        vec![ncbi_int4_from_double(x_drop_ungapped_bits); search_karlin.len()]
    } else {
        search_karlin
            .iter()
            .map(|blocks| {
                crate::core::x_drop::ungapped_x_dropoff(
                    x_drop_ungapped_bits,
                    blocks.ungapped.lambda,
                    1.0,
                )
            })
            .collect()
    };
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
        evalue_threshold,
        subject_besthit,
        x_qgroup_batch,
        percent_identity,
        score_matrix,
        nucl_score_table,
        x_dropoff_init,
    })
}

/// What `subject_schedule` reads.
struct ScheduleInputs<'a> {
    args: &'a BlastnArgs,
    subjects: SubjectSource<'a>,
    subjects_parallel: bool,
    num_threads: usize,
    two_stage_lookup_ref: Option<&'a TwoStageLookup>,
    diag_container_query_length: usize,
    window_size: usize,
    small_na_compressed_query: &'a [u8],
    scan_range: usize,
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
        subjects,
        subjects_parallel,
        num_threads,
        two_stage_lookup_ref,
        diag_container_query_length,
        window_size,
        small_na_compressed_query,
        scan_range,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        debug_window,
    } = inputs;
    // EXPERIMENT (LOSAT_X_BNFEWSUBJ, x_few_subj.rs): whether the seed stage's strips would take a
    // chunk (masked or not) of this `-subject` search that one subject's owner searches on this thread's
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    let use_parallel = subjects_parallel && !x_serial_subjects;
    crate::utils::threading::report_stage(
        "blastn",
        "subjects",
        subjects.num_subjects(),
        use_parallel,
    );
    SubjectSchedule {
        x_strips,
        x_kind,
        use_parallel,
    }
}

#[cfg(test)]
#[path = "subject_env_tests.rs"]
mod tests;

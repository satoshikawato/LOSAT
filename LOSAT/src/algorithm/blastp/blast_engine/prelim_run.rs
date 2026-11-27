//! The preliminary stage of a BLASTP call (`run_call`, after its kernel closures): whether a
//! subject set is searched on the pool (`use_parallel`), then the one branch that the call takes:
//! RPS-BLAST's engine (`rps_hooks::rps_search`), the PSSM search's (`pssm_hooks::pssm_prelim`),
//! none for a database without a valid query, the query blocks of LOSAT_X_BPQBLOCK (a call of
//! `run_resolved_in_pool` per block), a database's blocks (`database::search_blocks`), a subject
//! set on the pool (in LOSAT_X_BPSUBJRUN's runs or one task per subject), or the one-thread loop,
//! with the writer of `-culling_limit` or `-best_hit_*` or with the collector (`run_preliminary`). The
//! function is the statements of `run_call` at its place, in their order: `run_call` calls it
//! once, with the kernels and the hit lists, and its `?` passes the stage's errors on. The close
//! of the call's stream (`close_stream`) and the sort of its lists (`sort_preliminary_lists`),
//! which `run_call` calls after the stage, are here too.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1336
//!
//! This module is a child of `blast_engine` for the private items of that file it reads
//! (`preliminary_search_error`, `run_resolved_in_pool`, `x_block_input`, `PrelimCollector`,
//! `update_blastp_preliminary_hit_lists`, the subjects' scratch, reference and result, the
//! clocks, and for the parallel code `blastp_should_parallelize_stage`,
//! `blastp_subject_diag_offsets` and the two merges of the subjects' results, `PrelimFilter`) and
//! its sibling modules `database`, `pssm_hooks`, `rps_hooks`, `x_lookup_ahead` and
//! `x_subject_runs`.

use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::{
    blastp_timing_start, database, preliminary_search_error, pssm_hooks, rps_hooks,
    run_resolved_in_pool, update_blastp_preliminary_hit_lists, x_block_input, x_lookup_ahead,
    BlastpPrelimValues, BlastpRanges, BlastpReports, BlastpRuntimeLookup, BlastpSubjectParameters,
    BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, BlastpSubjectSource, BlastpTiming,
    PrelimCollector, PrelimFilter, PrelimOnly,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::{query_split, x_query_blocks};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::stats::search_space::SearchSpace;
use crate::stats::setup::{SearchSpaceOverflow, SubjectSpaceOverflow};
use crate::stats::{BlastGumbelBlk, KarlinParams};
use crate::utils::threading::SearchPool;

// The imports of the parallel code (the choice of the pool and the subject set on it).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    super::{
        blastp_should_parallelize_stage, blastp_subject_diag_offsets,
        merge_blastp_subject_results_filtered, merge_blastp_subject_results_in_order,
        x_subject_runs,
    },
    rayon::prelude::*,
};

/// What the preliminary stage reads of `run_call` (`run_preliminary`): each field is `run_call`'s
/// value of the same name (`Copy` values by value, the rest by reference, so that the moved
/// statements read them as they did in `run_call`), but `x_block_call`, which is
/// `x_block_inputs.is_some()`: the call is one block of LOSAT_X_BPQBLOCK's (`PrelimOnly::Block`),
/// whose read errors go to the search of all queries. The values that only the choice of the pool
/// reads are cfg'd as the parallel code.
pub(super) struct PrelimRun<'a> {
    // The call's values that its errors are made of (`preliminary_search_error`).
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) subject_source: BlastpSubjectSource<'a>,
    pub(super) subject_label: &'a str,
    pub(super) blastp_parallel_pool: &'a SearchPool<'a>,
    pub(super) ranges: &'a BlastpRanges<'a>,
    pub(super) reports: &'a BlastpReports<'a>,
    pub(super) query_batches: &'a [RangedBatch],
    pub(super) query_batch_skipped: &'a [bool],
    pub(super) query_chunk: bool,
    /// Whether a batch's search reads the database (`split_steps::CallSteps::Search`): none does
    /// where every searched batch is split into chunks without a valid context.
    pub(super) reads: bool,
    /// rpsblast: whether the engine counts each HSP's identities again for the reports
    /// (`rps_hooks::RpsSearchRun::recount_identities`), decided from the reports that the call's
    /// formats are (`call_setup::requested_reports`): none is a custom table.
    pub(super) rps_recount_identities: bool,
    // The subjects, the query blocks and the queries' contexts.
    pub(super) database: Option<&'a database::BlastpDatabase>,
    pub(super) subjects: &'a [EncodedProtein],
    pub(super) x_soft_subjects: Option<&'a [Vec<(i32, i32)>]>,
    pub(super) query_frames: &'a [Vec<QueryFrame>],
    pub(super) x_blocks: &'a Option<Vec<Range<usize>>>,
    pub(super) x_block_call: bool,
    pub(super) lookup: &'a BlastpRuntimeLookup,
    pub(super) contexts: &'a [QueryContext],
    /// psiblast `-phi_pattern`: the call's PHI-BLAST searches (`phi_hooks::phi_prelim`).
    pub(super) phi: Option<&'a super::phi_hooks::PhiCall<'a>>,
    // The statistics and the cutoffs.
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gapped_gumbel: &'a Option<BlastGumbelBlk>,
    pub(super) exact_search_spaces: &'a [Option<i64>],
    pub(super) preliminary_exact_spaces: &'a [Option<i64>],
    pub(super) preliminary_search_spaces: &'a [SearchSpace],
    pub(super) preliminary_option_spaces: &'a Option<Vec<i64>>,
    pub(super) per_subject: bool,
    pub(super) num_subjects: i32,
    pub(super) cutoff_subject_length: i32,
    /// rpsblast: the search units of the database (`call_setup::BlastpDbExtent::rps_units`).
    pub(super) rps_units: Option<&'a crate::algorithm::rps_search::volumes::RpsUnits>,
    pub(super) prelim_hitlist_size: usize,
    pub(super) prelim_evalue: f64,
    pub(super) subject_space_overflow: &'a SubjectSpaceOverflow,
    // `-db_soft_mask` and `-db_hard_mask`: the masking whose first read error of the masks ends a
    // database's blocks as a read error of its sequences (`db_masked.rs`).
    pub(super) db_masking: Option<&'a super::db_masked::BlastpDbMasking<'a>>,
    // The shape of the scan's tables.
    pub(super) window: i32,
    pub(super) diag_mask: i32,
    pub(super) offset_array_size: i32,
    pub(super) diag_array_size: i32,
    pub(super) total_query_span: i32,
    // The clocks (LOSAT_TIMING).
    pub(super) timing_enabled: bool,
    pub(super) timing: &'a Option<Arc<BlastpTiming>>,
    // The choice of the pool (`use_parallel`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) blastp_num_threads: usize,
}

/// What `run_call` reads of the preliminary stage after it (`run_preliminary`): RPS-BLAST's gapped
/// K of each query, for outfmt 0 (`rps_hooks::rps_search`; `None` for the other searches), and
/// whether the call's lists are the merged lists of several volumes of an RPS-BLAST database,
/// which are not sorted or cut and print no query warnings, and the HSPs of the profiles that
/// left those lists past the formatter's cut (`rps_hooks::RpsSearched::unlisted`).
pub(super) struct PrelimOutcome {
    pub(super) rps_gapped_k: Option<Vec<f64>>,
    pub(super) rps_multi_volume: bool,
    pub(super) rps_unlisted: Vec<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443
// A search without a database (`db_length` 0) recomputes the search spaces, cutoffs and links of
// every query from the length of each subject before its scan.
/// The preliminary stage with the search of one subject `search` (`BlastpPrelimKernel::
/// search_subject`, or with `-db_soft_mask` or `-db_hard_mask` the masked kernels of
/// `db_masked.rs`, chosen once per call so that a search without them keeps its kernel): with the
/// setup's statistics `batch_values` (`process_subject`), and where NCBI computes the statistics
/// of each subject (`per_subject`) with the subject's own, `subject_parameters` of its length with
/// the options' E-value (`process_own_subject`; a search space that overflows is noted in the
/// stage's `subject_space_overflow` and the subject gives no hits). A subject without letters takes
/// the setup's.
#[allow(clippy::too_many_arguments)]
pub(super) fn run_preliminary_with<S, P>(
    stage: &PrelimRun<'_>,
    search: &S,
    batch_values: &BlastpPrelimValues<'_>,
    subject_parameters: &P,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    prelim_filter: &mut Option<PrelimFilter>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<PrelimOutcome>
where
    S: Fn(
            &mut BlastpSubjectScratch,
            usize,
            BlastpSubjectRef<'_>,
            &mut BlastpSubjectResult,
            &BlastpPrelimValues<'_>,
        ) + Sync,
    P: Fn(usize) -> Result<BlastpSubjectParameters, SearchSpaceOverflow> + Sync,
{
    let evalue = stage.args.evalue;
    let overflow = stage.subject_space_overflow;
    let process_subject = |scratch: &mut BlastpSubjectScratch,
                           s_idx: usize,
                           subject: BlastpSubjectRef<'_>,
                           subject_result: &mut BlastpSubjectResult| {
        search(scratch, s_idx, subject, subject_result, batch_values)
    };
    let process_own_subject =
        |scratch: &mut BlastpSubjectScratch,
         s_idx: usize,
         subject: BlastpSubjectRef<'_>,
         subject_result: &mut BlastpSubjectResult| {
            if subject.aa_len == 0 {
                return search(scratch, s_idx, subject, subject_result, batch_values);
            }
            match subject_parameters(subject.aa_len) {
                Ok(parameters) => search(
                    scratch,
                    s_idx,
                    subject,
                    subject_result,
                    &parameters.values(evalue),
                ),
                Err(error) => {
                    subject_result.reset(s_idx);
                    overflow.note(error);
                }
            }
        };
    run_preliminary(
        stage,
        &process_subject,
        &process_own_subject,
        preliminary_hit_lists,
        x_compact_hit_lists,
        prelim_filter,
        outputs,
    )
}

/// The preliminary stage of the call (module header): the statements of `run_call` from the choice
/// of the pool to the end of the call's branch, in their order. The kernels are `run_call`'s
/// closures, `process_subject` with the setup's statistics and `process_own_subject` with each
/// subject's own, which only a database's blocks read (`per_subject`); each branch calls the one
/// it called in `run_call`. `preliminary_hit_lists` receive RPS-BLAST's final lists,
/// `x_compact_hit_lists` the other searches' compact lists, and `prelim_filter` is the writer of
/// `-culling_limit` or `-best_hit_*`, which `run_call` finishes. An error is the call's error: where a sequence file
/// cannot be read or RPS-BLAST's files fail, `preliminary_search_error` first writes to `outputs`
/// what NCBI writes before that batch fails (the outfmt 0 prolog, the batch's reader messages, the
/// reports of the unsearched batches before it), but a block of LOSAT_X_BPQBLOCK hands a read error
/// to the search of all queries; a database whose subject's search space overflows is refused
/// after its blocks (`subject_space_overflow`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:307-310
/// A database's threads are chosen in its driver (`database::search_blocks`). For a subject set
/// NCBI uses one thread (`-subject` resets `-num_threads`):
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3224-3227
/// No NCBI counterpart: `use_parallel`, the pool for a subject set, whose output is the one
/// thread's.
pub(super) fn run_preliminary<K, O>(
    stage: &PrelimRun<'_>,
    process_subject: &K,
    process_own_subject: &O,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    prelim_filter: &mut Option<PrelimFilter>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<PrelimOutcome>
where
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    O: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let PrelimRun {
        args,
        query_records,
        subject_source,
        subject_label,
        blastp_parallel_pool,
        ranges,
        reports,
        query_batches,
        query_batch_skipped,
        query_chunk,
        reads,
        rps_recount_identities,
        database,
        subjects,
        x_soft_subjects,
        query_frames,
        x_blocks,
        x_block_call,
        lookup,
        contexts,
        gapped_params,
        gapped_gumbel,
        exact_search_spaces,
        preliminary_exact_spaces,
        preliminary_search_spaces,
        preliminary_option_spaces,
        per_subject,
        num_subjects,
        cutoff_subject_length,
        rps_units,
        prelim_hitlist_size,
        prelim_evalue,
        subject_space_overflow,
        db_masking,
        window,
        diag_mask,
        offset_array_size,
        diag_array_size,
        total_query_span,
        timing_enabled,
        timing,
        phi,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        blastp_num_threads,
    } = *stage;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let use_parallel = {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1475
        // (A database has its own driver below, `database::search_blocks`.)
        database.is_none()
            && blastp_should_parallelize_stage(
                "prelim",
                args.num_threads,
                blastp_num_threads,
                subjects.len(),
                subjects.len(),
            )
    };
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let use_parallel = false;

    // rpsblast: the RPS engine's whole search (`rps_hooks::rps_search`): the final hit lists,
    // which the ordinary traceback in `run_call` after this stage hands on, and each query's
    // outfmt 0 gapped K.
    let mut rps_gapped_k: Option<Vec<f64>> = None;
    let mut rps_unlisted: Vec<usize> = Vec::new();
    // A database of several volumes: its merged lists are not sorted or cut after this stage
    // (`traceback::finish_hit_lists`), and it prints no query warnings; so even when no batch is searched (every query
    // invalid), as NCBI merges the volumes' empty results (review A of PS12a, S2's rpsblast half).
    let mut rps_multi_volume = args.rps_search.is_some()
        && database.is_some_and(|database| rps_hooks::volume_count(database) > 1);
    if let Some(rps) = args.rps_search.as_deref() {
        let run = rps_hooks::RpsSearchRun {
            contexts,
            gapped_params,
            gumbel: gapped_gumbel.is_some(),
            units: rps_units,
            search_spaces: exact_search_spaces,
            prelim_hitlist_size,
            batches: query_batches
                .iter()
                .map(|batch| batch.searched.clone())
                .collect(),
            recount_identities: rps_recount_identities,
        };
        let searched = rps_hooks::rps_search(
            rps,
            args,
            &run,
            database,
            blastp_parallel_pool,
            preliminary_hit_lists,
        )
        .map_err(|error| {
            match (database, error.downcast::<rps_hooks::RpsEngineError>()) {
                // NCBI's errors of the RPS files come from the first batch's set-up, after the
                // outfmt 0 prolog, whatever its queries' validity.
                (Some(database), Ok(engine_error)) => preliminary_search_error(
                    args,
                    query_records,
                    subject_source,
                    subject_label,
                    outputs,
                    blastp_parallel_pool,
                    ranges,
                    reports,
                    query_batches,
                    query_batch_skipped,
                    true,
                    database.name(),
                    &engine_error.0,
                ),
                (_, Ok(engine_error)) => anyhow::Error::new(engine_error),
                (_, Err(other)) => other,
            }
        })?;
        rps_gapped_k = Some(searched.gapped_k);
        rps_multi_volume = searched.volumes > 1;
        rps_unlisted = searched.unlisted;
    } else if let Some(psi) = args.psi_search.as_deref() {
        pssm_hooks::pssm_prelim(
            psi,
            args,
            lookup,
            &pssm_hooks::PssmPrelimRun {
                contexts,
                gapped_params,
                gumbel: gapped_gumbel.as_ref(),
                exact_search_spaces: preliminary_exact_spaces,
                search_spaces: preliminary_search_spaces,
                cutoff_subject_length,
                prelim_evalue,
                per_subject: per_subject.then(|| pssm_hooks::PsiSubjectStatistics {
                    option_space: preliminary_option_spaces
                        .as_ref()
                        .map_or(0, |spaces| spaces[0]),
                    setup_count: num_subjects,
                    dbsize: args.dbsize,
                }),
                window,
                diag_mask,
                offset_array_size,
                diag_array_size,
                total_query_span,
                prelim_hitlist_size,
            },
            database,
            subjects,
            blastp_parallel_pool,
            prelim_filter.as_mut(),
            x_compact_hit_lists,
        )
        .map_err(|error| match database {
            // A sequence file that cannot be read: NCBI's error of the batch's preliminary
            // search, as for blastp's database search below.
            Some(database) => match error.downcast::<database::PreliminaryReadError>() {
                Ok(read_error) => preliminary_search_error(
                    args,
                    query_records,
                    subject_source,
                    subject_label,
                    outputs,
                    blastp_parallel_pool,
                    ranges,
                    reports,
                    query_batches,
                    query_batch_skipped,
                    query_chunk,
                    database.name(),
                    read_error.as_bytes(),
                ),
                Err(other) => other,
            },
            None => error,
        })?;
    } else if let Some(phi) = phi {
        // psiblast `-phi_pattern`: the pattern scan in place of the word lookup
        // (`phi_hooks::phi_prelim`).
        let new_scratch = || {
            BlastpSubjectScratch::new(offset_array_size, diag_array_size, window, total_query_span)
        };
        super::phi_hooks::phi_prelim(
            phi,
            contexts,
            window,
            &new_scratch,
            database,
            subjects,
            blastp_parallel_pool,
            prelim_filter.as_mut(),
            prelim_hitlist_size,
            x_compact_hit_lists,
        )
        .map_err(|error| match database {
            // A sequence file that cannot be read: NCBI's error of the batch's preliminary
            // search, as for blastp's database search below.
            Some(database) => match error.downcast::<database::PreliminaryReadError>() {
                Ok(read_error) => preliminary_search_error(
                    args,
                    query_records,
                    subject_source,
                    subject_label,
                    outputs,
                    blastp_parallel_pool,
                    ranges,
                    reports,
                    query_batches,
                    query_batch_skipped,
                    query_chunk,
                    database.name(),
                    &read_error.0,
                ),
                Err(other) => other,
            },
            None => error,
        })?;
    } else if database.is_some() && (!reads || !contexts.iter().any(|ctx| ctx.is_valid)) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
        // No query of the call has Karlin-Altschul parameters: NCBI searches none of its
        // batches and reads no sequence of the database (a sequence file that cannot be read
        // fails no batch: NCBI 2.17.0 blastp, an all-X query, exit 0; line PS DECISIONS 30). Nor
        // does a call whose searched batches are all split into chunks without a valid context
        // (`reads`; NCBI 2.17.0 probe `g6_gapvalid_nopsq` of line RG item X-bpgap: exit 0).
    } else if let Some(blocks) = x_blocks.as_deref() {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the preliminary stage in blocks of whole batches, NCBI's
        // batches in order: each block a call of its own (`PrelimOnly::Block`: its lookup table,
        // diagonal table and scratch, freed when it returns), whose compact hit lists move to their
        // queries here. A query's preliminary hit list does not depend on the other queries of the
        // search (`x_query_blocks.rs`), and its HSPs' query coordinates are relative to its
        // context (`blastp_adjust_initial_hsp_offsets`), so only the query's index moves
        // (`query_split::rebase_hit_list`). Nothing is written before the redo and the reports
        // in `run_call` after this stage, which run over all queries as after one table. A
        // sequence file that cannot be read fails a block's search
        // (`database::PreliminaryReadError`, handed on by the block): it is this call's error,
        // made with this call's batches and reports (`preliminary_search_error`), as where one
        // table finds it.
        // The blocks run (`LOSAT_X_REACH`): recorded here, not where they are planned, so that a
        // call whose blocks do not run (no valid query of a database, or an error of the set-up)
        // does not record them.
        x_query_blocks::x_reach_block_search(true);
        for block in blocks {
            let input = x_block_input(query_batches, block);
            let block_input = ranges.query_input.slice(input.clone());
            let mut block_lists = Vec::new();
            let searched = run_resolved_in_pool(
                args.clone(),
                &query_records[block.clone()],
                subject_source,
                "",
                subject_label,
                outputs,
                blastp_parallel_pool,
                Some(PrelimOnly::Block {
                    lists: &mut block_lists,
                    subjects,
                    soft_subjects: x_soft_subjects,
                    frames: &query_frames[block.clone()],
                    first_query: block.start,
                }),
                &BlastpRanges {
                    input_records: &ranges.input_records[input.clone()],
                    query_input: &block_input,
                    subject_placements: ranges.subject_placements,
                    x_ordinal_base: ranges.x_ordinal_base + input.start,
                    // As the block search's blocks; a call of the preliminary stage only makes no
                    // search strategy requests (`run_resolved_in_pool`), the only reader.
                    searched_base: ranges.searched_base + input.start,
                    split_search_spaces: None,
                },
                &BlastpReports::prelim_only(false),
                x_lookup_ahead::XAheadBlock::none(),
            );
            if let Err(error) = searched {
                return Err(match error.downcast::<database::PreliminaryReadError>() {
                    Ok(read_error) => preliminary_search_error(
                        args,
                        query_records,
                        subject_source,
                        subject_label,
                        outputs,
                        blastp_parallel_pool,
                        ranges,
                        reports,
                        query_batches,
                        query_batch_skipped,
                        query_chunk,
                        database.map_or(b"".as_slice(), |database| database.name()),
                        read_error.as_bytes(),
                    ),
                    Err(other) => other,
                });
            }
            for (q_idx, list) in block.clone().zip(block_lists) {
                if let Some(mut list) = list {
                    query_split::rebase_hit_list(
                        &mut list,
                        q_idx as u32,
                        0,
                        contexts[q_idx].aa_len,
                    );
                    x_compact_hit_lists[q_idx] = Some(list);
                }
            }
        }
    } else if let Some(database) = database {
        // The preliminary stage of a database: its sequences are read in blocks and
        // searched in the order of their ordinal ids (`database::search_blocks`).
        let merge_sort_start = blastp_timing_start(timing_enabled);
        let new_scratch = || {
            BlastpSubjectScratch::new(offset_array_size, diag_array_size, window, total_query_span)
        };
        // The blocks of the database with the kernel of the setup's statistics or of each subject's
        // (`process_own_subject`), chosen once here.
        macro_rules! search_database_blocks {
            ($process_subject:expr) => {
                match prelim_filter.as_mut() {
                    Some(filter) => database::search_blocks(
                        database,
                        blastp_parallel_pool,
                        window,
                        &new_scratch,
                        $process_subject,
                        filter,
                    ),
                    None => database::search_blocks(
                        database,
                        blastp_parallel_pool,
                        window,
                        &new_scratch,
                        $process_subject,
                        &mut PrelimCollector {
                            x_compact_hit_lists,
                            prelim_hitlist_size,
                        },
                    ),
                }
            };
        }
        let searched = if per_subject {
            search_database_blocks!(process_own_subject)
        } else {
            search_database_blocks!(process_subject)
        };
        // `-db_soft_mask` and `-db_hard_mask`: the first read error of the masks by ordinal id
        // (`db_masked.rs`), a read error of the preliminary stage as the sequences' are.
        searched
            .and_then(
                |()| match db_masking.and_then(|masking| masking.take_error()) {
                    Some(error) => Err(error),
                    None => Ok(()),
                },
            )
            .map_err(|error| {
                // A sequence file that cannot be read: NCBI's error of the batch's
                // preliminary search (`database::PreliminaryReadError`). EXPERIMENT
                // (LOSAT_X_BPQBLOCK): a block of the preliminary stage hands it to the search of all
                // queries, whose batches and reports it is (`x_blocks`).
                match error.downcast::<database::PreliminaryReadError>() {
                    Ok(read_error) if x_block_call => anyhow::Error::new(read_error),
                    Ok(read_error) => preliminary_search_error(
                        args,
                        query_records,
                        subject_source,
                        subject_label,
                        outputs,
                        blastp_parallel_pool,
                        ranges,
                        reports,
                        query_batches,
                        query_batch_skipped,
                        query_chunk,
                        database.name(),
                        read_error.as_bytes(),
                    ),
                    Err(other) => other,
                }
            })?;
        subject_space_overflow.check()?;
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.prelim_merge_sort_ns,
                &timing.prelim_merge_sort_calls,
                merge_sort_start,
            );
        }
    } else if use_parallel {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1688
            let subject_diag_offsets = blastp_subject_diag_offsets(subjects, window);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1744-1748
            let pool = blastp_parallel_pool;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1554
            //
            // Parallel workers write one subject result into its oid-indexed
            // slot; the reducer below replays those slots in NCBI subject
            // iteration order.
            let mut subject_results: Vec<Option<BlastpSubjectResult>> =
                std::iter::repeat_with(|| None)
                    .take(subjects.len())
                    .collect();
            let new_scratch = || {
                BlastpSubjectScratch::new(
                    offset_array_size,
                    diag_array_size,
                    window,
                    total_query_span,
                )
            };
            // EXPERIMENT (LOSAT_X_BPSUBJRUN): runs of consecutive subjects, one scratch per
            // thread and one diagonal table fill per run (`x_subject_runs.rs`).
            if x_subject_runs::x_bp_subj_run() {
                crate::utils::xreach::x_reach!("LOSAT_X_BPSUBJRUN", true);
                x_subject_runs::x_search_subject_runs(
                    pool,
                    subjects,
                    &subject_diag_offsets,
                    &new_scratch,
                    process_subject,
                    &mut subject_results,
                );
            } else {
                pool.install(|| {
                    subject_results.par_iter_mut().enumerate().for_each_init(
                        new_scratch,
                        |scratch, (s_idx, subject_result_slot)| {
                            scratch.prepare_independent_subject(subject_diag_offsets[s_idx]);
                            let mut subject_result = BlastpSubjectResult::new(s_idx);
                            process_subject(
                                scratch,
                                s_idx,
                                BlastpSubjectRef::of(&subjects[s_idx]),
                                &mut subject_result,
                            );
                            *subject_result_slot = Some(subject_result);
                        },
                    );
                });
            }
            let merge_sort_start = blastp_timing_start(timing_enabled);
            match prelim_filter.as_mut() {
                Some(filter) => merge_blastp_subject_results_filtered(filter, subject_results),
                None => merge_blastp_subject_results_in_order(
                    x_compact_hit_lists,
                    subject_results,
                    prelim_hitlist_size,
                ),
            }
            if let Some(timing) = timing.as_ref() {
                BlastpTiming::record_call(
                    &timing.prelim_merge_sort_ns,
                    &timing.prelim_merge_sort_calls,
                    merge_sort_start,
                );
            }
        }
    } else if let Some(filter) = prelim_filter.as_mut() {
        // `-culling_limit` and `-best_hit_*`: the serial loop below with the option's writer in
        // place of the collector, chosen once here so that the default loop is unchanged.
        let mut scratch =
            BlastpSubjectScratch::new(offset_array_size, diag_array_size, window, total_query_span);
        let mut subject_result = BlastpSubjectResult::new(0);
        for s_idx in 0..subjects.len() {
            process_subject(
                &mut scratch,
                s_idx,
                BlastpSubjectRef::of(&subjects[s_idx]),
                &mut subject_result,
            );
            filter.add(&mut subject_result);
        }
    } else {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if x_subject_runs::x_bp_subj_run() {
            crate::utils::xreach::x_reach!("LOSAT_X_BPSUBJRUN", false);
        }
        let mut scratch =
            BlastpSubjectScratch::new(offset_array_size, diag_array_size, window, total_query_span);
        let mut subject_result = BlastpSubjectResult::new(0);
        for s_idx in 0..subjects.len() {
            process_subject(
                &mut scratch,
                s_idx,
                BlastpSubjectRef::of(&subjects[s_idx]),
                &mut subject_result,
            );
            let merge_sort_start = blastp_timing_start(timing_enabled);
            update_blastp_preliminary_hit_lists(
                x_compact_hit_lists,
                &mut subject_result,
                prelim_hitlist_size,
            );
            if let Some(timing) = timing.as_ref() {
                BlastpTiming::record_call(
                    &timing.prelim_merge_sort_ns,
                    &timing.prelim_merge_sort_calls,
                    merge_sort_start,
                );
            }
        }
    }
    Ok(PrelimOutcome {
        rps_gapped_k,
        rps_multi_volume,
        rps_unlisted,
    })
}

/// The close of the call's stream after its preliminary stage (`run_call`, after a split batch's
/// chunks merged their lists into the batch's).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1839-1846
/// `ADAPTIVE_CBS`: the stream's close trims each query's preliminary hit list
/// (`core::blast_hspstream`) before the subjects of the lists are read. RPSBLAST's engine
/// trims each volume's lists itself (`rps_search::engine`); PSIBLAST's traceback size is its
/// preliminary size, 1500 or more here, so it never trims.
pub(super) fn close_stream(
    args: &ResolvedBlastpArgs,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) {
    if args.rps_search.is_none()
        && crate::core::blast_hspstream::adaptive_cbs_close(
            args.comp_based_stats.is_enabled(),
            args.traceback_hitlist_size(),
        )
    {
        let hitlist_size = args.traceback_hitlist_size();
        for hit_list in preliminary_hit_lists.iter_mut().flatten() {
            hit_list.cbs_close_trim(hitlist_size);
        }
        for hit_list in x_compact_hit_lists.iter_mut().flatten() {
            hit_list.cbs_close_trim(hitlist_size);
        }
    }
}

/// The e-value sort of the preliminary lists, but rpsblast's (`run_call`, after `close_stream`):
/// the lists in the order in which the redo or the traceback reads them, before a database's
/// subjects are listed: the lists that its fetch excludes keep their place in that order as marks
/// (`database::exclude_fetched`, line RG item X4-excl). The sort's last key, the subject, gives
/// the same order by ordinal id as by rank (ranks are the listed ordinal ids in order).
pub(super) fn sort_preliminary_lists(
    args: &ResolvedBlastpArgs,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    timing: Option<&BlastpTiming>,
) {
    let merge_sort_start = blastp_timing_start(timing.is_some());
    // (rpsblast: the RPS engine sorted its lists, `rps_search::engine::finish_volume_list`; over
    // several volumes the merged lists keep `s_MergeAlignSet`'s order, which NCBI sorts nowhere.)
    if args.rps_search.is_none() {
        for hit_list_opt in preliminary_hit_lists.iter_mut() {
            if let Some(hit_list) = hit_list_opt.as_mut() {
                hit_list.sort_by_evalue();
            }
        }
        // The same sort of the compact lists.
        for hit_list in x_compact_hit_lists.iter_mut().flatten() {
            hit_list.sort_by_evalue();
        }
    }
    if let Some(timing) = timing {
        BlastpTiming::record_call(
            &timing.prelim_merge_sort_ns,
            &timing.prelim_merge_sort_calls,
            merge_sort_start,
        );
    }
}

//! The statistics of a batch's query contexts, set up before its subjects are searched: the
//! Karlin-Altschul blocks of every context (`karlin_blocks`), which contexts are valid, the checks
//! of the set-up on the subjects and the gap costs, and the effective search space of every
//! context (`effective_search_spaces`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
//! (`BlastSetup_ScoreBlkInit`), which calls c++/src/algo/blast/core/blast_stat.c:2736-2829
//! (`Blast_ScoreBlkKbpUngappedCalc`) and blast_setup.c:40-128 (`Blast_ScoreBlkKbpGappedCalc`);
//! blast_setup.c:887-999 (`BLAST_GapAlignSetUp`), which calls blast_setup.c:699-850
//! (`BLAST_CalcEffLengths`).
//!
//! This module is a child of `run` for the types it keeps private (`QueryContext`,
//! `SubjectSource`, `BatchStage`, `ShownSubjects`, `BatchOutputs`). The functions are statements
//! of `search_query_batch` moved here in their order; what they read is a struct of references or
//! `Copy` values, and what they make is what `search_query_batch` then reads under the same names.
//! An error of the set-up writes the outfmt 0 prolog of a file sink of the first batch through the
//! batch's outputs before it is returned, as it did in `search_query_batch`.

use super::{
    after_batch_reports, context_lengths, database, database_search_error, prolog_with_report,
    reborrow_outputs, write_format_prolog, write_pairwise_prologs, BatchOutputs, BatchStage,
    BlastnArgs, ContextKarlin, KarlinParams, NuclScoringSpec, QueryBatching, QueryContext,
    SequenceData, SetupOnly, ShownSubjects, SubjectSource, FIRST_BATCH_OUTPUTS,
};
use crate::algorithm::blastn::query_split::{ReadFails, SetupCheckFailure, UnmetChecks};
use crate::algorithm::blastn::scoring::{
    check_greedy_matrix_only, context_blocks, context_blocks_ungapped, karlin_error,
};
use crate::report::pairwise::BlastnIntro;
use crate::wording::ENGINE_ERROR;
use anyhow::Result;

/// What `context_statistics` reads: the batch's contexts with the ungapped blocks of their
/// compositions, the subjects, the stage of the batch, and what the outfmt 0 prolog of an error
/// shows.
pub(super) struct ContextInputs<'a> {
    pub(super) args: &'a BlastnArgs,
    /// The ungapped block of each context's composition (`None`: none could be computed).
    pub(super) context_ungapped: &'a [Option<KarlinParams>],
    pub(super) scoring_spec: NuclScoringSpec,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) query_lengths: &'a [usize],
    pub(super) batching: QueryBatching,
    pub(super) first_batch: bool,
    pub(super) stage: BatchStage<'a>,
    /// Whether NCBI searches the batch chunk by chunk (`search_query_chunks`): the subjects and
    /// the greedy set-up are then checked by the set-up of each chunk (`BLAST_GapAlignSetUp` in
    /// the chunk's `CPrelimSearchRunner`), after the chunk's masks, or, if no chunk has a valid
    /// context, by the batch's traceback set-up (`SThreadLocalDataArraySetup`), after the chunks.
    pub(super) split: bool,
    pub(super) subjects: SubjectSource<'a>,
    pub(super) seq_data: &'a SequenceData<'a>,
    /// What the outfmt 0 prolog of an error introduces.
    pub(super) intro: BlastnIntro,
    pub(super) shown: &'a ShownSubjects,
}

/// The statistics of the batch's contexts, which `search_query_batch` reads under the same names.
pub(super) struct ContextStatistics {
    /// The blocks of each context, `None` for an invalid context.
    pub(super) context_karlin: Vec<Option<ContextKarlin>>,
    pub(super) round_down_evalue_score: bool,
    pub(super) context_valid: Vec<bool>,
    pub(super) any_valid_context: bool,
    /// The blocks the search reads for each context (an invalid context's slot holds those of a
    /// valid one).
    pub(super) search_karlin: Vec<ContextKarlin>,
    /// The effective search space of each context.
    pub(super) query_eff_searchsp: Vec<i64>,
}

/// The Karlin-Altschul blocks, the validity and the effective search space of every context of the
/// batch, after the checks of the set-up; an error writes what NCBI writes before it (the outfmt 0
/// prolog of a file sink of the first batch) through `outputs`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
/// (`BlastSetup_ScoreBlkInit`, which `BLAST_MainSetUp` calls at :652 for
/// `CSetupFactory::CreateScoreBlock`, c++/src/algo/blast/api/setup_factory.cpp:121-158) and
/// blast_setup.c:887-999 (`BLAST_GapAlignSetUp`, which c++/src/algo/blast/core/blast_engine.c:1665
/// calls first in `Blast_RunPreliminarySearchWithInterrupt`, before the subjects).
pub(super) fn context_statistics<const UNGAPPED: bool>(
    inputs: &ContextInputs<'_>,
    mut outputs: BatchOutputs<'_, '_>,
) -> Result<ContextStatistics> {
    let (context_karlin, round_down_evalue_score) =
        karlin_blocks::<UNGAPPED>(inputs, reborrow_outputs(&mut outputs))?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:324-332
    // An invalid context has an infinite cutoff (`batch_setup::batch_cutoffs`), so it gets no hit.
    // Its slot holds the blocks of a valid context (or placeholder blocks when none is valid),
    // which no hit reads.
    let context_valid: Vec<bool> = context_karlin.iter().map(Option::is_some).collect();
    let any_valid_context = context_valid.iter().any(|&valid| valid);
    let (search_karlin, query_eff_searchsp) = effective_search_spaces::<UNGAPPED>(
        inputs,
        outputs,
        &context_karlin,
        &context_valid,
        any_valid_context,
    )?;
    Ok(ContextStatistics {
        context_karlin,
        round_down_evalue_score,
        context_valid,
        any_valid_context,
        search_karlin,
        query_eff_searchsp,
    })
}

/// The Karlin-Altschul blocks of every context: the gapped blocks made over the given ungapped
/// ones (`None` for an invalid context), and whether the e-values round the score down; the error
/// of scoring options without Karlin-Altschul values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
/// (`BlastSetup_ScoreBlkInit`), which calls c++/src/algo/blast/core/blast_stat.c:2736-2829
/// (`Blast_ScoreBlkKbpUngappedCalc`) and, for a gapped search, blast_setup.c:40-128
/// (`Blast_ScoreBlkKbpGappedCalc`).
fn karlin_blocks<const UNGAPPED: bool>(
    inputs: &ContextInputs<'_>,
    outputs: BatchOutputs<'_, '_>,
) -> Result<(Vec<Option<ContextKarlin>>, bool)> {
    let &ContextInputs {
        args,
        context_ungapped,
        scoring_spec,
        query_contexts,
        query_lengths,
        batching,
        first_batch,
        intro,
        shown,
        ..
    } = inputs;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2758-2792
    // Every query context (strand) has the ungapped block of its composition (kbp_std:
    // the gap trigger and the ungapped X-drop) and a gapped block (kbp_gap: e-values,
    // bit scores, cutoffs and the length adjustment); `scoring.rs` has the NCBI
    // references. NCBI raises an unsupported scoring system when it sets up the first
    // batch with a valid query, after the outfmt 0 prolog.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:513-529
    // An ungapped (`-ungapped`) search has the ungapped blocks only (`context_blocks_ungapped`).
    let context_blocks_result = if UNGAPPED {
        Ok(context_blocks_ungapped(context_ungapped, &scoring_spec))
    } else {
        context_blocks(context_ungapped, &scoring_spec)
    };
    let (context_karlin, round_down_evalue_score) = match context_blocks_result {
        Ok(blocks) => blocks,
        Err(message) => {
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:170-172
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2784-2791
            // NCBI raises the error only when the first message of the batch is it: an
            // invalid context puts its warning first, and NCBI then goes on without the
            // gapped blocks (it crashes; with a first batch of invalid queries only, a later
            // batch raises it). LOSAT reproduces the error only for a first batch of valid
            // queries (the error does not depend on the queries, so the first batch with a valid
            // query meets it; NCBI writes the results and warnings of the batches before it).
            // The contexts of a query without letters are invalid from their set-up and put
            // no message here (`context_ungapped_blocks`). LOSAT's refusal comes where NCBI's
            // search of the batch would crash: after the reports of the batches before (whose
            // queries were all invalid) in both paths of LOSAT_X_RPTSTREAM, or in the first
            // batch after the prologs (`setup_failure`).
            if query_contexts
                .iter()
                .zip(context_ungapped)
                .any(|(context, block)| context.len != 0 && block.is_none())
            {
                let refusal = anyhow::anyhow!(
                    "these scoring options have no Karlin-Altschul values, and the query batch (up to {} residues) has an invalid query; NCBI BLAST+ goes on without the gapped values and crashes, which is not supported by LOSAT's {}",
                    batching.batch_size,
                    args.app().name()
                );
                return Err(setup_failure(
                    refusal,
                    first_batch,
                    false,
                    outputs,
                    intro,
                    shown,
                ));
            }
            // A query without letters (also a `-query_loc` interval that holds none) keeps its
            // warning in the error (`karlin_error`).
            let empty_queries = query_lengths.iter().map(|&length| length == 0);
            // A batch after the first: every batch before it had only invalid queries (one
            // with a valid query would have raised the error), and NCBI has written their
            // reports (`run_in_pool`).
            if !first_batch {
                return Err(after_batch_reports(karlin_error(&message, empty_queries)));
            }
            // The prolog of a file sink (the others were written before the search).
            let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
            for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
                if prolog_with_report(format, output_format) {
                    let posts = crate::report::query_warnings::html_posts_sink(
                        outputs.html_posts.as_deref_mut(),
                        &mut *outputs.diagnostics,
                    );
                    write_format_prolog(format, output_format, intro, shown, posts)?;
                }
            }
            return Err(karlin_error(&message, empty_queries));
        }
    };
    Ok((context_karlin, round_down_evalue_score))
}

/// LOSAT's refusal of the greedy set-up of `BLAST_GapAlignSetUp` (`s_BlastGreedyAlignMemAlloc`),
/// which NCBI reaches with a valid context only (with no valid context, NCBI extends nothing):
/// matrix-only scoring with the greedy extension, where NCBI divides by zero or does not end
/// (`check_greedy_matrix_only`; an ungapped search sets it up too).
pub(super) fn greedy_setup_refusal(args: &BlastnArgs, ungapped: bool) -> Option<anyhow::Error> {
    check_greedy_matrix_only(args, ungapped).err()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-996
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-219
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:310-314
// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:231-253
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1831-1837
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback_mt_priv.c:218-228
/// The checks of `BLAST_GapAlignSetUp` after the Karlin-Altschul block check: first a statistics
/// count of 0 of a positive total (`common::database::check_setup_count`, where NCBI divides by
/// zero and crashes), then `stats::setup::setup_subject_length` (the total divided by `num_seqs`,
/// cast to `Int4` and then `Uint4`, is 0: NCBI's engine error `subject_error`,
/// `SUBJECT_LENGTH_TOO_SHORT` in the preliminary stage and `TRACEBACK_FAILED` in the traceback
/// stage), then LOSAT's refusals of the greedy set-up (`greedy_setup_refusal`).
#[cfg(test)]
pub(super) fn gap_align_setup_error(
    args: &BlastnArgs,
    subjects: SubjectSource<'_>,
    seq_data: &SequenceData<'_>,
    ungapped: bool,
    subject_error: &str,
) -> Option<anyhow::Error> {
    let total_length = seq_data.db_len_total;
    let (total_length, num_seqs) = if total_length > 0 {
        (total_length, seq_data.db_num_seqs)
    } else {
        // A database with a total length of 0 holds sequences without letters only:
        // the length of its first sequence (`BlastSeqSrcGetSeqLen` of ordinal id 0) is 0.
        (
            subjects.prepared().map_or(0, |prepared| {
                if prepared.records.is_empty() {
                    0
                } else {
                    prepared.letters().len(0)
                }
            }) as i64,
            1,
        )
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // A count of 0 (of a positive total) is the division by 0 where NCBI crashes, LOSAT's
    // refusal (`common::database::check_setup_count`, D1a-5 of line RG) at this batch, after
    // the reports of the batches before (NCBI's output before its crash in the tabular
    // formats), as the error below; a run whose batches have no valid context does not reach
    // it and is reported as NCBI's. Else `stats::setup::setup_subject_length`: the total
    // divided by `num_seqs`, cast to `Int4` and then to `Uint4`, is 0 (a negative `NSEQ`
    // passes; one of 2^31 fails).
    crate::algorithm::common::database::check_setup_count(
        total_length,
        num_seqs,
        false,
        args.app().name(),
    )
    .err()
    .or_else(|| {
        crate::stats::setup::setup_subject_length(total_length, num_seqs, None)
            .is_none()
            .then(|| {
                anyhow::Error::from(crate::cli::NativeError {
                    exit: 3,
                    message: format!("{ENGINE_ERROR} {subject_error}\n").into(),
                })
            })
    })
    // The greedy set-up after the subjects: NCBI 2.17.0 reports `The average subject
    // length is too short` (or `Traceback failed` in the batch traceback set-up when no
    // chunk had a valid context, exit status 3) with scores that LOSAT refuses
    // for the greedy extension, also in a split query batch (`CHUNK_SIZE=200`; line RG D1a-6).
    .or_else(|| greedy_setup_refusal(args, ungapped))
}

/// A failure of a batch with a valid context — NCBI's own error or LOSAT's refusal where NCBI
/// crashes or LOSAT does not follow it — at the point where NCBI meets it: in a batch after the
/// first, wrapped in `ErrorAfterBatchReports` so `run_in_pool` writes the reports of the batches
/// before (`write_batch_reports`, hits included, as for `database_search_error`); in the first
/// batch, after writing the prologs that come with the reports (`write_pairwise_prologs`).
///
/// The average subject length, statistics count of 0 (D1a-5), and greedy set-up refusals (D1a-6)
/// do not depend on the queries, so in a batch after the first every earlier batch had no valid
/// query. A 64-bit effective search space overflow (`effective_search_spaces`,
/// `batch_setup::batch_spaces`, or mid-search at `subject_space_overflow.check()` in `run.rs`;
/// D1a-6 review A F2) depends on the query length (and subject length for per-subject spaces), so
/// valid earlier batches with hits may precede it; the handler writes their reports alike. So does
/// the refusal of the ungapped word cutoff (`batch_setup::batch_cutoff_set`: a context of 2^30
/// letters or more, or an average subject length whose search space overflows), which a batch with a
/// long query meets after batches with shorter ones.
pub(super) fn setup_failure(
    error: anyhow::Error,
    first_batch: bool,
    setup_only: bool,
    outputs: BatchOutputs<'_, '_>,
    intro: BlastnIntro,
    shown: &ShownSubjects,
) -> anyhow::Error {
    if !first_batch {
        return after_batch_reports(error);
    }
    // In `refused_batch_setup_error` (`setup_only`), a refusal that is not NCBI's `NativeError`
    // gives way to the scope refusal, whose `fail_batch!` writes the prologs itself.
    if setup_only && !error.is::<crate::cli::NativeError>() {
        return error;
    }
    let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
    match write_pairwise_prologs(outputs, output_formats, intro, shown, false) {
        Ok(()) => error,
        Err(write_error) => write_error,
    }
}

/// Where a `-db` search first fails to read the database in a split query batch (`SplitBatchSteps`):
/// `ReadFails::Setup` when the volume has no sequences (`SeqDb::first_oid_read_fails`: OID 0 in
/// `BLAST_GapAlignSetUp`, blast_setup.c:925-937), or `ReadFails::Search` when a sequence file
/// cannot be opened in the preliminary stage.
pub(super) fn split_read_fails(subjects: SubjectSource<'_>) -> Option<ReadFails> {
    let SubjectSource::Database(database) = subjects else {
        return None;
    };
    if database.first_oid_read_fails() {
        Some(ReadFails::Setup)
    } else if database.sequence_file_read_error().is_some() {
        Some(ReadFails::Search)
    } else {
        None
    }
}

/// The read error corresponding to `split_read_fails(subjects)` when `SplitBatchSteps::before_read`
/// returns `ChunkStep::End`.
pub(super) fn split_read_error(subjects: SubjectSource<'_>) -> Option<anyhow::Error> {
    let SubjectSource::Database(database) = subjects else {
        return None;
    };
    if database.first_oid_read_fails() {
        Some(database::PreliminaryReadError::from_static("Database memory map file error").into())
    } else {
        database.sequence_file_read_error()
    }
}

/// The checks of `BLAST_GapAlignSetUp` on the subjects and greedy extension parameters after the
/// OID 0 read: a statistics count of 0 (refusal, rule 3), the average subject length (NCBI's
/// engine error `SUBJECT_LENGTH_TOO_SHORT` in the preliminary stage, `Traceback failed` via
/// `UnmetChecks::TracebackFailed` in traceback setup when no chunk of a split batch had a valid
/// context), and the greedy set-up refusals.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:475-480
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-939
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-993
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-219
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:310-314
///
/// The subject set has no BLAST database statistics, and BLASTN has no Gumbel block
/// (blast_setup.c:475-480): a batch with a valid context stops NCBI when the subjects have fewer
/// letters than subjects (subjects without letters count with none). A batch without a valid
/// context is not searched (`CheckInternalData`), so it does not.
pub(super) fn gap_align_setup_checks(
    args: &BlastnArgs,
    subjects: SubjectSource<'_>,
    ungapped: bool,
) -> std::result::Result<(), SetupCheckFailure> {
    let metadata = subjects.metadata();
    let total_length = metadata.db_len_total;
    let (total_length, num_seqs) = if total_length > 0 {
        (total_length, metadata.db_num_seqs)
    } else {
        // A database with a total length of 0 holds sequences without letters only:
        // the length of its first sequence (`BlastSeqSrcGetSeqLen` of ordinal id 0) is 0.
        (
            subjects.prepared().map_or(0, |prepared| {
                if prepared.records.is_empty() {
                    0
                } else {
                    prepared.letters().len(0)
                }
            }) as i64,
            1,
        )
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // A count of 0 (of a positive total) is the division by 0 where NCBI crashes, LOSAT's refusal
    // (`common::database::check_setup_count`, D1a-5 of line RG, rule 3) at this batch, after the
    // reports of the batches before (NCBI's output before its crash in the tabular formats), as
    // the error below; a run whose batches have no valid context does not reach it and is reported
    // as NCBI's. Else `stats::setup::setup_subject_length`: the total divided by `num_seqs`, cast
    // to `Int4` and then to `Uint4`, is 0 (a negative `NSEQ` passes; one of 2^31 fails).
    if let Err(error) = crate::algorithm::common::database::check_setup_count(
        total_length,
        num_seqs,
        false,
        args.app().name(),
    ) {
        return Err(SetupCheckFailure::Refusal(error));
    }
    if crate::stats::setup::setup_subject_length(total_length, num_seqs, None).is_none() {
        return Err(SetupCheckFailure::Engine(anyhow::Error::from(
            crate::cli::NativeError {
                exit: 3,
                message: format!(
                    "{ENGINE_ERROR} {}\n",
                    crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT
                )
                .into(),
            },
        )));
    }
    // `BLAST_GapAlignSetUp` of this batch or chunk: the subjects first, then the greedy set-up
    // below (blast_setup.c:969-979 before :993, `BLAST_GapAlignStructNew`).
    // The greedy set-up after the subjects: NCBI 2.17.0 reports `The average subject length is too
    // short` (exit status 3) with gap costs or scores that LOSAT refuses for the greedy extension,
    // also in a split query batch (each chunk checks them, `CHUNK_SIZE=200`; line RG D1a-6).
    if let Some(error) = greedy_setup_refusal(args, ungapped) {
        return Err(SetupCheckFailure::Refusal(error));
    }
    Ok(())
}

/// Formats `SplitBatchSteps::unmet_checks()` for BLASTN's traceback setup after every chunk of a
/// split batch had no valid context: NCBI's `Traceback failed` (exit 3) for an engine setup error,
/// or the preserved refusal (count of 0, greedy setup).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:231-253
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1831-1837
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback_mt_priv.c:219-226
pub(super) fn traceback_setup_error(
    unmet: UnmetChecks,
    first_batch: bool,
    outputs: BatchOutputs<'_, '_>,
    intro: BlastnIntro,
    shown: &ShownSubjects,
) -> anyhow::Error {
    let error = match unmet {
        UnmetChecks::TracebackFailed => anyhow::Error::from(crate::cli::NativeError {
            exit: 3,
            message: format!(
                "{ENGINE_ERROR} {}\n",
                crate::blastinput::app::TRACEBACK_FAILED
            )
            .into(),
        }),
        UnmetChecks::Refusal(refusal) => refusal,
    };
    setup_failure(error, first_batch, false, outputs, intro, shown)
}

/// Checks of `BLAST_GapAlignSetUp` on the subjects and greedy extension parameters for a batch
/// searched whole (`!split`):
/// - The subject set has no BLAST database statistics, and BLASTN has no Gumbel block
///   (blast_setup.c:475-480): a batch with a valid context stops NCBI when the subjects have fewer
///   letters than subjects (subjects without letters count with none). A batch without a valid
///   context is not searched (`CheckInternalData`), so it does not.
/// - A database without sequences (makeblastdb builds none; a volume whose index holds none) has a
///   total length of 0 (unless an alias file's LENGTH or STATS_TOTLEN gives one:
///   `SeqDb::first_oid_read_fails`), so the setup asks its reader for the length of OID 0, which
///   it does not have: NCBI's search stage reports the reader's error as a database file that
///   cannot be read (local_blast.cpp:240-251). Shown with the oracle (campaign dbsearch-20261012,
///   `$BUILD_ROOT/dbsearch/db4/smoke/empty`: the nucleotide fixture's index with no sequence): the
///   outfmt 0 prolog, then `BLAST engine error: Database memory map file error`, exit status 3,
///   for outfmt 0, 6 and 7.
/// - A split batch checks these in its chunks (`ContextInputs::split`, via `SplitBatchSteps`): the
///   first chunk with a valid context stops NCBI there, after the masks of the chunks up to it
///   (NCBI BLAST+ 2.17.0 on 2026-11-12, line NP M3c review `c3`: `CHUNK_SIZE=600` and subjects
///   with fewer letters than subjects post the registry of `-window_masker_taxid` for the batch
///   and for chunk 0, then `The average subject length is too short`), while OID 0 reads stop at
///   chunk 0 before its validity (`blast_setup.c:933`), and a split batch with no valid chunk
///   meets the checks in traceback setup (`traceback_setup_error`). The search spaces of the batch
///   (`SplitQuery_SetEffectiveSearchSpace`) check nothing.
/// - Count of 0 (rule 3, `check_setup_count`, line RG D1a-5) and greedy setup after the subjects
///   (line RG D1a-6): see `gap_align_setup_checks`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:475-480
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-993
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-219
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:310-314
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:448-458
#[allow(clippy::too_many_arguments)]
pub(super) fn check_gap_align_setup(
    args: &BlastnArgs,
    subjects: SubjectSource<'_>,
    ungapped: bool,
    first_batch: bool,
    setup_only: bool,
    outputs: BatchOutputs<'_, '_>,
    intro: BlastnIntro,
    shown: &ShownSubjects,
) -> Result<()> {
    if let SubjectSource::Database(database) = subjects {
        if database.first_oid_read_fails() {
            return Err(database_search_error(
                database::PreliminaryReadError::from_static("Database memory map file error")
                    .into(),
                first_batch,
                outputs,
                intro,
                shown,
            ));
        }
    }
    if let Err(failure) = gap_align_setup_checks(args, subjects, ungapped) {
        return Err(setup_failure(
            failure.into_error(),
            first_batch,
            setup_only,
            outputs,
            intro,
            shown,
        ));
    }
    Ok(())
}

/// The checks of the set-up on the subjects (a database without sequences, the average subject
/// length) and on the greedy set-up, then the blocks the search reads for each context and the
/// effective search space of each context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:887-999
/// (`BLAST_GapAlignSetUp`), which calls c++/src/algo/blast/core/blast_setup.c:699-850
/// (`BLAST_CalcEffLengths`).
fn effective_search_spaces<const UNGAPPED: bool>(
    inputs: &ContextInputs<'_>,
    mut outputs: BatchOutputs<'_, '_>,
    context_karlin: &[Option<ContextKarlin>],
    context_valid: &[bool],
    any_valid_context: bool,
) -> Result<(Vec<ContextKarlin>, Vec<i64>)> {
    let &ContextInputs {
        args,
        query_contexts,
        first_batch,
        stage,
        split,
        subjects,
        seq_data,
        intro,
        shown,
        ..
    } = inputs;
    let setup_only = matches!(
        stage,
        BatchStage::Search {
            setup_only: SetupOnly::Whole | SetupOnly::Split,
            ..
        }
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-939
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-219
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:310-314
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:448-458
    // A split batch checks these through `SplitBatchSteps` in `search_query_chunks` /
    // `refused_batch_setup_error` (`BatchStage::ChunkPrelim` and `ContextInputs::split`).
    if any_valid_context && !split && !matches!(stage, BatchStage::ChunkPrelim { .. }) {
        check_gap_align_setup(
            args,
            subjects,
            UNGAPPED,
            first_batch,
            setup_only,
            reborrow_outputs(&mut outputs),
            intro,
            shown,
        )?;
    }
    let search_karlin: Vec<ContextKarlin> = {
        let fallback = context_karlin
            .iter()
            .flatten()
            .next()
            .copied()
            .unwrap_or(ContextKarlin {
                ungapped: KarlinParams::default(),
                gapped: KarlinParams::default(),
            });
        context_karlin
            .iter()
            .map(|blocks| blocks.unwrap_or(fallback))
            .collect()
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-846
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-721
    // A `-dbsize` above 0 replaces the total length in the effective lengths (with `-db`
    // and with `-subject`); the number of sequences, the average subject length above, the
    // batch mixer and the reports keep the real numbers.
    let db_len_total_i64 = if args.db_length() > 0 {
        args.db_length()
    } else {
        seq_data.db_len_total
    };
    let db_num_seqs = seq_data.db_num_seqs;
    let program = args.app().name();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:778-786
    // An invalid context keeps the search space of the options, which is 0, except in a query
    // chunk: NCBI sets the options to the search spaces of the batch's contexts
    // (`SplitQuery_SetEffectiveSearchSpace`), and a chunk's context takes the one at its own
    // index in the chunk, not its context in the batch (0, for an invalid context of the
    // batch, has the chunk compute its own).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:177-182
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:690-693
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:826-846
    let query_eff_searchsp: Vec<i64> = match query_contexts
        .iter()
        .zip(&search_karlin)
        .zip(context_valid)
        .enumerate()
        .map(|(index, ((ctx, blocks), &valid))| {
            let options_searchsp = match stage {
                BatchStage::Search { .. } => args.options_searchsp(),
                BatchStage::ChunkPrelim {
                    batch_eff_searchsp, ..
                } => batch_eff_searchsp[index],
            };
            if !valid || options_searchsp != 0 {
                return Ok(options_searchsp);
            }
            context_lengths(
                ctx.len as i64,
                db_len_total_i64,
                db_num_seqs,
                &blocks.gapped,
                program,
            )
            .map(|(_, effective_search_space)| effective_search_space)
        })
        .collect::<Result<_>>()
    {
        Ok(spaces) => spaces,
        Err(error) => {
            return Err(setup_failure(
                error,
                first_batch,
                setup_only,
                outputs,
                intro,
                shown,
            ));
        }
    };
    Ok((search_karlin, query_eff_searchsp))
}

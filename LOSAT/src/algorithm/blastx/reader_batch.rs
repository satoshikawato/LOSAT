//! The values that NCBI computes over a whole query batch (LOSAT's reader batch, `BATCH_SIZE`
//! letters) and LOSAT's search units of at most 10002 letters must share.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    parameters::{score_block_for, subject_parameters, ContextParameters},
    query_setup::{context_frame, prepare_queries, prepare_unit_queries, PreparedQueryBatch},
    subjects::Subjects,
};
use crate::stats::tables::KarlinParams;
use anyhow::Result;
use std::{ops::Range, sync::Arc};

/// What the units of one reader batch share (phase 1 of the search of a reader batch of several
/// units; `ReaderBatchShape::reduce`). A reader batch of one unit has none: its unit computes
/// each value itself, from its own contexts, which are then the batch's.
#[derive(Debug)]
pub(crate) struct ReaderBatchShape {
    /// NCBI's `num_queries` of the batch: with more than one the batch is not split.
    pub(crate) query_count: usize,
    /// `contexts[last_context].query_offset` plus `query_length`, of the concatenated batch.
    pub(crate) concatenated_length: usize,
    /// One more than `last_context`, of the batch.
    pub(crate) contexts: usize,
    /// Whether a context of the batch has residues, and the batch's empty queries: the setup of
    /// a batch without a valid context fails with one message per empty query
    /// (`query_setup::prepare_unit_queries`).
    pub(crate) valid: bool,
    pub(crate) empty_queries: usize,
    /// The Karlin block of the first valid context with the smallest ungapped lambda, over the
    /// batch's contexts.
    pub(crate) smallest_lambda: Option<KarlinParams>,
    /// `word_params->cutoff_score_min` of the batch (an ungapped search with sum statistics).
    pub(crate) word_cutoff_min: Option<i32>,
    /// `hit_params->cutoff_score_min` of the batch (the composition redo).
    pub(crate) redo_cutoff_min: Option<i32>,
    /// Each context of the batch, in order: its length, whether it is valid and its ungapped
    /// Karlin block, of which the statistics of each subject take the batch's word cutoff minimum
    /// (`parameters::subject_word_cutoff_min`; empty before the units are prepared).
    pub(crate) batch_contexts: Vec<(usize, bool, KarlinParams)>,
    /// NCBI's error of the batch's set-up when the score block cannot load the `-matrix` of an
    /// ungapped search: the message of each query of the batch (`query_setup::unknown_matrix_message`),
    /// which every unit's preparation raises (`query_setup::prepare_unit_queries`).
    pub(crate) matrix_error: Option<Arc<str>>,
}

/// One unit of a reader batch, prepared in phase 1: its queries, translated and filtered, and
/// the statistical blocks of its contexts. Phase 2 searches the unit from these, so that
/// nothing is translated or filtered twice.
pub(crate) struct UnitPrep {
    pub(crate) prepared: PreparedQueryBatch,
    pub(crate) params: Vec<ContextParameters>,
    /// The unit's part of the batch's `word_cutoff_min` and `redo_cutoff_min` (its valid
    /// contexts, with the batch's average length), when the search reads them.
    word_cutoff_min: Option<i32>,
    redo_cutoff_min: Option<i32>,
}

impl PreparedQueryBatch {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:954-957
    // The query info is the query batch's, not a part of it: `concat_qlen` and the number of
    // contexts are those of the reader batch when this unit is a part of one.
    /// `last.query_offset + last.query_length` of the batch NCBI searches.
    pub(crate) fn batch_concatenated_length(&self) -> usize {
        match &self.reader_shape {
            Some(shape) => shape.concatenated_length,
            None => self
                .contexts
                .last()
                .map_or(0, |last| last.offset + last.length),
        }
    }

    /// One more than `last_context`, of the batch NCBI searches.
    pub(crate) fn batch_context_count(&self) -> usize {
        match &self.reader_shape {
            Some(shape) => shape.contexts,
            None => self.contexts.len(),
        }
    }
}

impl ReaderBatchShape {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:91-117
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:954-957
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-88
    // The reader batch is one query info: its contexts follow one another through the units, a
    // context of positive length taking its length and a separator. The layout of the contexts
    // follows from the records alone (`context_frame`, as `prepare_queries` lays them out), so
    // it is known before any unit is prepared. The smallest lambda is the first strict minimum
    // over the contexts in order, as above; the two cutoff minima are the minima over the units'.
    /// The values of a reader batch that follow from its records: the query count, the
    /// concatenated length, the context count and whether one is valid (phase 1 starts with
    /// these; the units' cutoff minima are taken with the average length they give).
    pub(crate) fn layout(records: &[FastaRecord], options: &ResolvedOptions) -> Self {
        let (mut offset, mut concatenated_length, mut contexts) = (0, 0, 0);
        for record in records {
            for i in 0..6 {
                let (len, enabled) = context_frame(record.sequence.len(), i, &options.strand);
                let length = if enabled { len } else { 0 };
                concatenated_length = offset + length;
                contexts += 1;
                if length > 0 {
                    offset += length + 1;
                }
            }
        }
        ReaderBatchShape {
            query_count: records.len(),
            concatenated_length,
            contexts,
            valid: offset > 0,
            empty_queries: records.iter().filter(|q| q.sequence.is_empty()).count(),
            smallest_lambda: None,
            word_cutoff_min: None,
            redo_cutoff_min: None,
            batch_contexts: Vec::new(),
            matrix_error: options
                .unknown_matrix
                .then(|| super::query_setup::matrix_error_message(records, options).into()),
        }
    }

    /// Whether NCBI's set-up fails the batch: no context is valid and a query has no letters
    /// (`query_setup::prepare_unit_queries`, whose error has one message per empty query).
    pub(crate) fn set_up_fails(&self) -> bool {
        !self.valid && self.empty_queries > 0
    }

    /// The shape of a reader batch of `layout` from its units' preparations, in order (phase 1
    /// ends with it).
    pub(crate) fn reduce<'u>(
        layout: &ReaderBatchShape,
        units: impl IntoIterator<Item = &'u UnitPrep>,
    ) -> Self {
        let mut shape = ReaderBatchShape {
            query_count: layout.query_count,
            concatenated_length: layout.concatenated_length,
            contexts: layout.contexts,
            valid: layout.valid,
            empty_queries: layout.empty_queries,
            smallest_lambda: None,
            word_cutoff_min: None,
            redo_cutoff_min: None,
            batch_contexts: Vec::new(),
            matrix_error: layout.matrix_error.clone(),
        };
        for unit in units {
            shape.batch_contexts.extend(
                unit.prepared
                    .contexts
                    .iter()
                    .zip(&unit.params)
                    .map(|(context, p)| (context.length, p.valid, p.ungapped)),
            );
            for p in unit.params.iter().filter(|p| p.valid) {
                if p.ungapped.lambda > 0.0
                    && shape
                        .smallest_lambda
                        .is_none_or(|smallest| p.ungapped.lambda < smallest.lambda)
                {
                    shape.smallest_lambda = Some(p.ungapped);
                }
            }
            shape.word_cutoff_min = [shape.word_cutoff_min, unit.word_cutoff_min]
                .into_iter()
                .flatten()
                .min();
            shape.redo_cutoff_min = [shape.redo_cutoff_min, unit.redo_cutoff_min]
                .into_iter()
                .flatten()
                .min();
        }
        shape
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980 (`BLAST_GapAlignSetUp`,
// cited at `parameters::subject_parameters`)
/// The set-up of the records `queries` of a reader batch alone (`native::add_batch_scope`): the
/// error of the set-up of the subjects where a context of theirs is valid (a database whose
/// sequences could not be read, `UnreadableSubjects`, as the search fails its first batch with a
/// valid context; the average subject length, the hit-saving parameters: `subject_parameters`),
/// then of the search's first read of the database (a sequence file that cannot be opened,
/// `Subjects::sequence_file_read_error`; line WQA D10 round 6, review 4 M2), for the caller's
/// `search_error`; `None` where the set-up passes or no context is valid. NCBI
/// splits no batch of two queries or more (`split::calculate_num_chunks`), and a refused batch
/// that is set up holds two (line WQA D10 round 5, 2026-11-20).
pub(crate) fn set_up_error(
    queries: &[FastaRecord],
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<anyhow::Error> {
    let prepared = prepare_unit_queries(queries, None, options).ok()?;
    let mut params = score_block_for(&prepared, options.scoring);
    if !params.iter().any(|p| p.valid) {
        return None;
    }
    if subjects.is_unreadable() {
        return Some(super::subjects::UnreadableSubjects.into());
    }
    let numbers = subjects.numbers().ok()?;
    subject_parameters(&prepared, &mut params, options, numbers, None)
        .err()
        .or_else(|| subjects.sequence_file_read_error())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:416-416,990-990
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2450-2453
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:954-957
// NCBI sets up and searches a reader batch (one `BATCH_SIZE` query batch) as a whole, so these
// values are taken over all its contexts. LOSAT searches the batch in units of at most 10002
// letters. Phase 1 prepares every unit (translation, filtering, the statistical blocks of the
// contexts, and the unit's cutoff minima with the batch's average length, which the layout of the
// records gives), reduces the batch's values, and hands each unit its preparation with the shape
// attached; phase 2 searches the units as before, from these. A batch with no valid context is
// not searched, and a database refused before its sequences are read (`Subjects::Unreadable`)
// ends the search before they are used, so the minima are not taken then.
/// Phase 1 for one unit (`queries`) of a reader batch of `layout`. The unit's preparation has
/// the layout attached; `attach` replaces it with the reduced shape.
pub(crate) fn prepare_unit(
    queries: &[FastaRecord],
    layout: &Arc<ReaderBatchShape>,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Result<UnitPrep> {
    let mut prepared = prepare_unit_queries(queries, Some(layout), options)?;
    let params = score_block_for(&prepared, options.scoring);
    prepared.reader_shape = Some(Arc::clone(layout));
    let kappa = options.composition != 0 && options.gapped;
    let ungapped_link = super::linking::computes_link_cutoffs(options);
    let (mut word_cutoff_min, mut redo_cutoff_min) = (None, None);
    if (kappa || ungapped_link)
        && params.iter().any(|p| p.valid)
        && !subjects.is_unreadable()
        && !subjects.is_empty()
    {
        let mut minima = params.clone();
        // The word cutoffs are the preliminary stage's, with its option spaces where they differ
        // from the setup's (`parameters::preliminary_spaces`); the redo reads the setup's.
        let mut word_minima = None;
        let computed = (|| {
            let numbers = subjects.numbers()?;
            subject_parameters(&prepared, &mut minima, options, numbers, None)?;
            if ungapped_link {
                if let Some(spaces) = super::parameters::preliminary_spaces(
                    &prepared,
                    &params,
                    options,
                    numbers,
                    matches!(subjects, Subjects::Stream(_)),
                )? {
                    let mut preliminary = params.clone();
                    subject_parameters(
                        &prepared,
                        &mut preliminary,
                        options,
                        numbers,
                        Some(&spaces),
                    )?;
                    word_minima = Some(preliminary);
                }
            }
            anyhow::Ok(())
        })();
        match computed {
            Ok(()) => {
                let valid = || minima.iter().filter(|p| p.valid);
                word_cutoff_min = ungapped_link
                    .then(|| {
                        word_minima
                            .as_ref()
                            .unwrap_or(&minima)
                            .iter()
                            .filter(|p| p.valid)
                            .map(|p| p.word_cutoff)
                            .min()
                    })
                    .flatten();
                redo_cutoff_min = kappa.then(|| valid().map(|p| p.hit_cutoff).min()).flatten();
            }
            // A database that fails the whole check takes no minima, as when it was refused
            // before the search: its error comes where the unit is searched
            // (`search::search_core`, `XPastHead`).
            Err(_) if subjects.x_damaged() => {}
            Err(error) => return Err(error),
        }
    }
    Ok(UnitPrep {
        prepared,
        params,
        word_cutoff_min,
        redo_cutoff_min,
    })
}

/// A reader batch of one unit prepared: its queries translated and filtered, and the statistical
/// blocks of its contexts. `search::search_core` prepares such a batch with it when it is searched,
/// and the capped batch flow (`LOSAT_X_BXCAPFILL`, `x_bx_mem.rs`) ahead of its search. Such a batch
/// has no shape (its unit's values are the batch's) and no cutoff minima to reduce.
pub(crate) fn prepare_lone_unit(
    queries: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<UnitPrep> {
    let prepared = prepare_queries(queries, options)?;
    let params = score_block_for(&prepared, options.scoring);
    Ok(UnitPrep {
        prepared,
        params,
        word_cutoff_min: None,
        redo_cutoff_min: None,
    })
}

/// The end of phase 1 for the units of a reader batch: the batch's shape replaces the layout.
pub(crate) fn attach<'u>(
    units: impl IntoIterator<Item = &'u mut UnitPrep>,
    shape: ReaderBatchShape,
) {
    let shape = Arc::new(shape);
    for unit in units {
        unit.prepared.reader_shape = Some(Arc::clone(&shape));
    }
}

/// Phase 1 of the search of a reader batch of `ranges.len()` > 1 units, all at once (with a pool
/// of 2+ threads the units are prepared side by side). The flow of `native.rs` runs the same
/// steps as tasks of its own.
pub(crate) fn prepare_reader_batch(
    records: &[FastaRecord],
    ranges: &[Range<usize>],
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<Vec<UnitPrep>> {
    let x_t = super::x_timing::lap();
    let layout = Arc::new(ReaderBatchShape::layout(records, options));
    let prepare = |range: &Range<usize>| -> Result<UnitPrep> {
        prepare_unit(&records[range.clone()], &layout, options, subjects)
    };
    #[cfg(feature = "parallel")]
    let mut units: Vec<UnitPrep> = if pool.enabled() {
        use rayon::prelude::*;
        pool.install(|| {
            ranges
                .par_iter()
                .with_max_len(1)
                .map(prepare)
                .collect::<Result<Vec<_>>>()
        })?
    } else {
        ranges.iter().map(prepare).collect::<Result<Vec<_>>>()?
    };
    #[cfg(not(feature = "parallel"))]
    let mut units: Vec<UnitPrep> = {
        let _ = pool;
        ranges.iter().map(prepare).collect::<Result<Vec<_>>>()?
    };
    super::x_timing::batch(super::x_timing::Batch::QuerySetup, x_t);
    let shape = ReaderBatchShape::reduce(&layout, &units);
    attach(units.iter_mut(), shape);
    Ok(units)
}

#[cfg(test)]
#[path = "reader_batch_tests.rs"]
mod tests;

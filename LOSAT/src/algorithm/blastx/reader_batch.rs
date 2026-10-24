//! The values that NCBI computes over a whole query batch (LOSAT's reader batch, `BATCH_SIZE`
//! letters) and LOSAT's search units of at most 10002 letters must share.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    parameters::{score_block, subject_parameters, ContextParameters},
    query_setup::{context_frame, prepare_queries, PreparedQueryBatch},
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
    /// `contexts[last_context].query_offset + query_length` of the concatenated batch.
    pub(crate) concatenated_length: usize,
    /// `last_context + 1` of the batch.
    pub(crate) contexts: usize,
    /// The Karlin block of the first valid context with the smallest ungapped lambda, over the
    /// batch's contexts.
    pub(crate) smallest_lambda: Option<KarlinParams>,
    /// `word_params->cutoff_score_min` of the batch (an ungapped search with sum statistics).
    pub(crate) word_cutoff_min: Option<i32>,
    /// `hit_params->cutoff_score_min` of the batch (the composition redo).
    pub(crate) redo_cutoff_min: Option<i32>,
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
    // ```c
    //          Int4 concat_qlen =
    //              query_info->contexts[query_info->last_context].query_offset +
    //              query_info->contexts[query_info->last_context].query_length;
    //          Int4 avg_qlen = concat_qlen / (query_info->last_context + 1);
    // ```
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

    /// `last_context + 1` of the batch NCBI searches.
    pub(crate) fn batch_context_count(&self) -> usize {
        match &self.reader_shape {
            Some(shape) => shape.contexts,
            None => self.contexts.len(),
        }
    }
}

impl ReaderBatchShape {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:92-117
    // ```c
    // static double
    // s_BlastFindSmallestLambda(Blast_KarlinBlk** kbp_in,
    //                           const BlastQueryInfo* query_info,
    //                           Blast_KarlinBlk** kbp_out)
    // {
    //     Int4 i;
    //     double min_lambda = (double) INT4_MAX;
    //
    //     ASSERT(kbp_in && query_info);
    //
    //     for (i=query_info->first_context; i<=query_info->last_context; i++) {
    //         ASSERT(s_BlastKarlinBlkIsValid(kbp_in[i]) ==
    //                query_info->contexts[i].is_valid);
    //         if (s_BlastKarlinBlkIsValid(kbp_in[i])) {
    //             if (min_lambda > kbp_in[i]->Lambda)
    //             {
    //                 min_lambda = kbp_in[i]->Lambda;
    //                 if (kbp_out)
    //                   *kbp_out = kbp_in[i];
    //             }
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:954-957
    // ```c
    //          Int4 concat_qlen =
    //              query_info->contexts[query_info->last_context].query_offset +
    //              query_info->contexts[query_info->last_context].query_length;
    //          Int4 avg_qlen = concat_qlen / (query_info->last_context + 1);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-88
    // ```c++
    //     if (index) {
    //         Uint4 prev_loc = qinfo->contexts[index-1].query_offset;
    //         Uint4 prev_len = qinfo->contexts[index-1].query_length;
    //
    //         Uint4 shift = prev_len ? prev_len + 1 : 0;
    //
    //         qinfo->contexts[index].query_offset = prev_loc + shift;
    // ```
    // The reader batch is one query info: its contexts follow one another through the units, a
    // context of positive length taking its length and a separator. The layout of the contexts
    // follows from the records alone (`context_frame`, as `prepare_queries` lays them out), so
    // it is known before any unit is prepared. The smallest lambda is the first strict minimum
    // over the contexts in order, as above; the two cutoff minima are the minima over the units'.
    /// The values of a reader batch that follow from its records: the query count, the
    /// concatenated length and the context count (phase 1 starts with these; the units' cutoff
    /// minima are taken with the average length they give).
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
            smallest_lambda: None,
            word_cutoff_min: None,
            redo_cutoff_min: None,
        }
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
            smallest_lambda: None,
            word_cutoff_min: None,
            redo_cutoff_min: None,
        };
        for unit in units {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:416-416,990-990
// ```c
//    parameters->cutoff_score_min = cutoff_min;
// ...
//       params->cutoff_score_min = cutoff_min;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2450-2453
// ```c
//     if (do_link_hsps) {
//         ASSERT(hitParams->link_hsp_params != NULL);
//         cutoff_s =
//             (int) (hitParams->cutoff_score_min * context->localScalingFactor);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:954-957
// ```c
//          Int4 concat_qlen =
//              query_info->contexts[query_info->last_context].query_offset +
//              query_info->contexts[query_info->last_context].query_length;
//          Int4 avg_qlen = concat_qlen / (query_info->last_context + 1);
// ```
// NCBI sets up and searches a reader batch (one `BATCH_SIZE` query batch) as a whole, so these
// values are taken over all its contexts. LOSAT searches the batch in units of at most 10002
// letters. Phase 1 prepares every unit (translation, filtering, the statistical blocks of the
// contexts, and the unit's cutoff minima with the batch's average length, which the layout of the
// records gives), reduces the batch's values, and hands each unit its preparation with the shape
// attached; phase 2 searches the units as before, from these. A batch with no valid context is
// not searched, and subjects that cannot be read end the search before they are used, so the
// minima are not taken then.
/// Phase 1 for one unit (`queries`) of a reader batch of `layout`. The unit's preparation has
/// the layout attached; `attach` replaces it with the reduced shape.
pub(crate) fn prepare_unit(
    queries: &[FastaRecord],
    layout: &Arc<ReaderBatchShape>,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Result<UnitPrep> {
    let mut prepared = prepare_queries(queries, options)?;
    let params = score_block(&prepared);
    prepared.reader_shape = Some(Arc::clone(layout));
    let kappa = options.composition == 2 && options.gapped;
    let ungapped_link = !options.gapped && options.sum_stats;
    let (mut word_cutoff_min, mut redo_cutoff_min) = (None, None);
    if (kappa || ungapped_link)
        && params.iter().any(|p| p.valid)
        && !subjects.is_unreadable()
        && !subjects.is_empty()
    {
        let mut minima = params.clone();
        subject_parameters(
            &prepared,
            &mut minima,
            options,
            subjects.total_length()?,
            i64::try_from(subjects.len())?,
            subjects.shortest_length(),
            None,
        )?;
        let valid = || minima.iter().filter(|p| p.valid);
        word_cutoff_min = ungapped_link
            .then(|| valid().map(|p| p.word_cutoff).min())
            .flatten();
        redo_cutoff_min = kappa.then(|| valid().map(|p| p.hit_cutoff).min()).flatten();
    }
    Ok(UnitPrep {
        prepared,
        params,
        word_cutoff_min,
        redo_cutoff_min,
    })
}

/// The end of phase 1 for one unit: the batch's shape replaces the layout.
pub(crate) fn attach(unit: &mut UnitPrep, shape: &Arc<ReaderBatchShape>) {
    unit.prepared.reader_shape = Some(Arc::clone(shape));
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
    let shape = Arc::new(ReaderBatchShape::reduce(&layout, &units));
    for unit in &mut units {
        attach(unit, &shape);
    }
    Ok(units)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(letters: usize) -> FastaRecord {
        // Deterministic pseudo-random nucleotides: the contexts have valid Karlin blocks.
        let mut state = letters as u64 + 12345;
        let sequence = (0..letters)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 60) as usize & 3]
            })
            .collect();
        FastaRecord {
            internal_id: String::new(),
            title: String::new(),
            sequence,
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
        }
    }

    /// The shape of a reader batch of several units is the shape of the batch prepared as one
    /// unit: the concatenation and the context count (from the records' layout), the smallest
    /// lambda and the composition redo's minimum cutoff do not depend on where the batch is cut,
    /// with either strand or both.
    #[test]
    fn the_shape_of_the_units_is_the_shape_of_the_batch_as_one_unit() {
        let records: Vec<_> = [3000, 301, 7000, 2, 5000, 1, 4999, 90].map(record).to_vec();
        let ranges = super::super::query_setup::batch_ranges(&records, 6000);
        assert!(ranges.len() > 2);
        let pool = crate::utils::threading::SearchPool::x_serial();
        for strand in ["both", "plus", "minus"] {
            let mut options = super::super::stage_c_tests::options();
            options.strand = strand.into();
            options.composition = 2;
            let subjects: Subjects<'_> = (&records[..1]).into();
            let units = prepare_reader_batch(&records, &ranges, &options, subjects, &pool).unwrap();
            let layout = Arc::new(ReaderBatchShape::layout(&records, &options));
            let whole = prepare_unit(&records, &layout, &options, subjects).unwrap();
            let shape = units[0].prepared.reader_shape.as_ref().unwrap();
            assert_eq!(shape.query_count, records.len());
            assert_eq!(shape.contexts, whole.prepared.contexts.len());
            let last = whole.prepared.contexts.last().unwrap();
            assert_eq!(shape.concatenated_length, last.offset + last.length);
            let smallest = whole
                .params
                .iter()
                .filter(|p| p.valid)
                .map(|p| p.ungapped.lambda)
                .fold(f64::INFINITY, f64::min);
            assert_eq!(shape.smallest_lambda.unwrap().lambda, smallest);
            assert!(shape.redo_cutoff_min.is_some());
            assert_eq!(shape.redo_cutoff_min, whole.redo_cutoff_min);
            // Every unit sees the one shape.
            for unit in &units {
                assert_eq!(
                    unit.prepared.batch_concatenated_length(),
                    shape.concatenated_length
                );
                assert_eq!(unit.prepared.batch_context_count(), shape.contexts);
            }
        }
    }

    /// A unit of 19,410 nucleotides or more that is alone in its search is split into chunks;
    /// the same unit as a part of a reader batch of several queries is not (NCBI splits no
    /// batch of several queries).
    #[test]
    fn a_unit_of_a_reader_batch_of_several_queries_is_not_split() {
        use super::super::split::split_queries;
        let (records, mut full, options) =
            super::super::stage_c_tests::fixture("S12_split19411.fna");
        assert_eq!(split_queries(&records, &full, &options).unwrap().len(), 2);
        full.reader_shape = Some(Arc::new(ReaderBatchShape {
            query_count: 2,
            concatenated_length: 2 * full.batch_concatenated_length(),
            contexts: 2 * full.contexts.len(),
            smallest_lambda: None,
            word_cutoff_min: None,
            redo_cutoff_min: None,
        }));
        assert_eq!(split_queries(&records, &full, &options).unwrap().len(), 1);
    }
}

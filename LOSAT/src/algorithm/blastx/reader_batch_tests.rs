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
        title: Vec::new(),
        sequence,
        lowercase_masks: Vec::new(),
        source: None,
        around: None,
        skipped: None,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:949-976
/// NCBI's sum statistics carry their e-value over the contexts of the query batch: the units of a
/// reader batch take the cutoffs of the batch set up as one unit, also a unit whose first valid
/// context follows valid contexts of earlier units. The subjects' length of 5 residues makes the
/// carried cutoff (2) differ from the first context's (1).
#[test]
fn the_units_carry_the_sum_statistics_of_their_reader_batch() {
    use super::super::parameters::SubjectNumbers;
    let records: Vec<_> = [3000, 301, 7000, 2, 5000, 1, 4999, 90].map(record).to_vec();
    let ranges = super::super::query_setup::batch_ranges(&records, 6000);
    assert!(ranges.len() > 2);
    let pool = crate::utils::threading::SearchPool::x_serial();
    let mut options = super::super::stage_c_tests::options();
    options.evalue = 1e-5;
    let numbers = SubjectNumbers {
        total: 5,
        count: 1,
        gumbel_total: 5,
        shortest: 5,
        average: 5,
        request_count: 1,
    };
    let subjects: Subjects<'_> = (&records[..1]).into();
    let units = prepare_reader_batch(&records, &ranges, &options, subjects, &pool).unwrap();
    let layout = Arc::new(ReaderBatchShape::layout(&records, &options));
    let whole = prepare_unit(&records, &layout, &options, subjects).unwrap();
    let cutoffs = |prep: &UnitPrep| {
        let mut params = prep.params.clone();
        subject_parameters(&prep.prepared, &mut params, &options, numbers, None).unwrap();
        params
            .iter()
            .map(|p| (p.valid, p.hit_cutoff))
            .collect::<Vec<_>>()
    };
    let batch = cutoffs(&whole);
    let first = batch.iter().position(|&(valid, _)| valid).unwrap();
    assert_eq!(batch[first].1, 1);
    assert!(batch[first + 1..]
        .iter()
        .filter(|&&(valid, _)| valid)
        .all(|&(_, cutoff)| cutoff == 2));
    assert!(!units[0].prepared.batch_valid_before);
    assert!(units[1..]
        .iter()
        .all(|unit| unit.prepared.batch_valid_before));
    let by_units: Vec<_> = units.iter().flat_map(cutoffs).collect();
    assert_eq!(by_units, batch);
}

/// NCBI validates the query batch as a whole: a unit of only empty records is prepared when
/// another unit of its reader batch has a valid context, and in a reader batch without one
/// every unit fails with one message per empty record of the batch.
#[test]
fn whether_a_unit_has_no_data_is_its_reader_batchs() {
    let options = super::super::stage_c_tests::options();
    let pool = crate::utils::threading::SearchPool::x_serial();
    let valid: Vec<_> = [12000, 0, 0].map(record).to_vec();
    let ranges = super::super::query_setup::batch_ranges(&valid, 10002);
    assert_eq!(ranges, [0..1, 1..3]);
    let subjects: Subjects<'_> = (&valid[..1]).into();
    let units = prepare_reader_batch(&valid, &ranges, &options, subjects, &pool).unwrap();
    assert!(!units[1].params.iter().any(|p| p.valid));
    assert!(prepare_queries(&valid[1..], &options).is_err());
    let none: Vec<_> = [0]
        .into_iter()
        .chain([2; 5001])
        .chain([0])
        .map(record)
        .collect();
    let ranges = super::super::query_setup::batch_ranges(&none, 10002);
    assert_eq!(ranges, [0..5002, 5002..5003]);
    let layout = Arc::new(ReaderBatchShape::layout(&none, &options));
    for range in ranges {
        let error = prepare_unit(&none[range], &layout, &options, subjects)
            .err()
            .expect("a reader batch without a valid context");
        assert_eq!(
            error.to_string(),
            "Warning: Sequence contains no data ".repeat(2)
        );
    }
}

/// A unit of 19,410 nucleotides or more that is alone in its search is split into chunks;
/// the same unit as a part of a reader batch of several queries is not (NCBI splits no
/// batch of several queries).
#[test]
fn a_unit_of_a_reader_batch_of_several_queries_is_not_split() {
    use super::super::split::{split_queries, DEFAULT_SIZES};
    let (records, mut full, options) = super::super::stage_c_tests::fixture("S12_split19411.fna");
    assert_eq!(
        split_queries(&records, &full, &options, DEFAULT_SIZES)
            .unwrap()
            .len(),
        2
    );
    full.reader_shape = Some(Arc::new(ReaderBatchShape {
        query_count: 2,
        concatenated_length: 2 * full.batch_concatenated_length(),
        contexts: 2 * full.contexts.len(),
        valid: true,
        empty_queries: 0,
        smallest_lambda: None,
        word_cutoff_min: None,
        redo_cutoff_min: None,
        batch_contexts: Vec::new(),
        matrix_error: None,
    }));
    assert_eq!(
        split_queries(&records, &full, &options, DEFAULT_SIZES)
            .unwrap()
            .len(),
        1
    );
}

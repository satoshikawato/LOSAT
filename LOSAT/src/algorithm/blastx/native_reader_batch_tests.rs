use super::*;

fn record(letters: usize) -> FastaRecord {
    FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence: vec![b'A'; letters],
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    }
}

/// One unit's result for `queries` queries, skipped or not.
fn unit_results(queries: usize, skipped: bool) -> Vec<BatchResults> {
    let records: Vec<_> = (0..queries).map(|_| record(30)).collect();
    vec![BatchResults {
        query_ordinal: 0,
        queries: vec![Vec::new(); queries],
        prepared: super::super::query_setup::prepare_queries(
            &records,
            &super::super::stage_c_tests::options(),
        )
        .unwrap(),
        parameters: Vec::new(),
        search_skipped: skipped,
        request_space: None,
        unlisted: Vec::new(),
        shown: Vec::new(),
        gapped_k: None,
    }]
}

/// A reader batch (any `BATCH_SIZE`) is split into units of at most 10002 letters, as a run
/// without `BATCH_SIZE` batches its queries; the last unit is marked, and the units cover the
/// batch in order.
#[test]
fn a_reader_batch_is_split_into_units_of_10002_letters() {
    let batch = [4000, 4000, 4000, 100, 10002, 5, 20000, 1]
        .map(record)
        .to_vec();
    let options = super::super::stage_c_tests::options();
    let units = XUnit::split(&batch, &options, false);
    let ranges: Vec<_> = units.iter().map(|unit| unit.range.clone()).collect();
    // 4000 + 4000 + 4000 reaches 10002 at the third record; 100 + 10002 at the fifth;
    // 5 + 20000 at the seventh; the rest is the last unit.
    assert_eq!(ranges, vec![0..3, 3..5, 5..7, 7..8]);
    assert_eq!(
        units.iter().map(|unit| unit.last).collect::<Vec<_>>(),
        [false, false, false, true]
    );
    assert_eq!(units[1].queries().len(), 2);
    // A reader batch of at most 10002 letters (the default, or a smaller `BATCH_SIZE`) is
    // one unit.
    let small = XUnit::split(&batch[..3], &options, false);
    assert_eq!(small.len(), 1);
    assert!(small[0].last);
}

/// A unit whose phase 1 did not end with a preparation (it failed) is prepared again with
/// the reader batch's layout, so its error is the batch's: here a batch without a valid
/// context, whose two empty records are in different units.
#[test]
fn a_unit_without_its_phase_1_fails_as_its_reader_batch() {
    let batch: Vec<_> = [0]
        .into_iter()
        .chain([2; 5001])
        .chain([0])
        .map(record)
        .collect();
    let options = super::super::stage_c_tests::options();
    let units = XUnit::split(&batch, &options, false);
    assert_eq!(units.len(), 2);
    let subjects: Subjects<'_> = (&batch[1..2]).into();
    for unit in &units {
        let error = unit
            .take_prep(&options, subjects)
            .err()
            .expect("a reader batch without a valid context");
        assert_eq!(
            error.to_string(),
            "Warning: Sequence contains no data ".repeat(2)
        );
    }
}

/// The results of the units become the reader batch's at its last unit: ordinals count from
/// the batch's first query, and a unit with no valid context is skipped only when no unit
/// had one.
#[test]
fn the_units_results_are_the_reader_batchs() {
    let mut pending = ReaderBatchResults::default();
    assert!(pending.add(0, false, unit_results(3, false)).is_none());
    assert!(pending.add(3, false, unit_results(1, true)).is_none());
    let results = pending.add(4, true, unit_results(2, false)).unwrap();
    let ordinals: Vec<_> = results.iter().map(|r| r.query_ordinal).collect();
    assert_eq!(ordinals, [0, 3, 4]);
    assert!(results.iter().all(|r| !r.search_skipped));
    // Nothing is left for the next reader batch, and a batch of skipped units is skipped.
    let results = pending
        .add(0, false, unit_results(1, true))
        .or_else(|| pending.add(1, true, unit_results(1, true)))
        .unwrap();
    assert_eq!(results.len(), 2);
    assert!(results.iter().all(|r| r.search_skipped));
}

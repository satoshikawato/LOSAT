use super::*;
#[test]
// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blastfilter_unit_test.cpp:1003-1027 BlastxLowerCaseMaskProteinLocations
fn ncbi_active_lowercase_frame_assertions() {
    let expected = [
        (0, 25),
        (0, 24),
        (0, 24),
        (3034, 3059),
        (3034, 3058),
        (3034, 3058),
    ];
    for (frame, expected) in [1, 2, 3, -1, -2, -3].into_iter().zip(expected) {
        assert_eq!(
            dna_to_protein_masks(&[(0, 75)], 9180, frame),
            vec![expected]
        );
    }
}
#[test]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:995-1003
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1089-1096
fn touching_and_duplicate_masks_preserve_lookup_intervals() {
    assert_eq!(
        combine_masks(vec![(0, 0), (0, 0), (1, 2), (2, 3)]),
        vec![(0, 0), (0, 0), (1, 2), (2, 3)]
    );
    let c = QueryContext {
        query_index: 0,
        frame: 1,
        offset: 0,
        length: 2,
        is_valid: true,
        lowercase_masks: vec![],
        masks: vec![(0, 0), (0, 0)],
        dna_masks: vec![],
    };
    assert_eq!(complement_masks(&[c]), vec![(1, -1), (1, 1)]);
}
#[test]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:139-165
fn ncbi_batch_crossing_record_is_retained() {
    let record = |n| FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence: vec![b'A'; n],
        lowercase_masks: vec![],
        warnings: vec![],
        source: None,
        around: None,
        skipped: None,
    };
    assert_eq!(
        batch_ranges(&[record(10001), record(1), record(3)], BATCH_SIZE),
        vec![0..2, 2..3]
    );
    assert_eq!(
        batch_ranges(&[record(10002), record(3)], BATCH_SIZE),
        vec![0..1, 1..2]
    );
    assert_eq!(
        batch_ranges(&[record(10003), record(3)], BATCH_SIZE),
        vec![0..1, 1..2]
    );
}

// `BATCH_SIZE` sets the limit (`GetQueryBatchSize`): the record that crosses it ends its
// batch, and a limit above the default keeps the longer batch whole.
#[test]
fn the_batch_size_is_the_limit_of_a_batch() {
    let record = |n| FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence: vec![b'A'; n],
        lowercase_masks: vec![],
        warnings: vec![],
        source: None,
        around: None,
        skipped: None,
    };
    let records = [
        record(100),
        record(100),
        record(100),
        record(100),
        record(100),
    ];
    assert_eq!(batch_ranges(&records, 200), vec![0..2, 2..4, 4..5]);
    assert_eq!(batch_ranges(&records, 201), vec![0..3, 3..5]);
    assert_eq!(
        batch_ranges(&records, 1),
        vec![0..1, 1..2, 2..3, 3..4, 4..5]
    );
    assert_eq!(batch_ranges(&records, 30000), vec![0..5]);
}

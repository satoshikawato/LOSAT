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

// Line NP item M4b-3b2 (decision 57 of series M): a reader batch is one search unit only where a
// unit cannot give the batch's rows alone, with `-use_sw_tback` and where the `-db_hard_mask`
// walk (the BLASTX word 3 and window 40) keeps a sequence of the database apart. Of the databases
// of `tests/db_search/dbmaskp`: ov_th10's B starts on a new table (its rows of a query depend on
// the query's place in NCBI's batch: rows `db.pmask.bx.hard.ov_th10.qfc.*`), ov_e3's chunk at 70
// reads back to the entries of the one word at 60, and the touching masks of ov_e3c read nothing
// back.
#[test]
fn a_reader_batch_is_one_unit_only_where_a_unit_cannot_give_its_rows() {
    use crate::blastinput::db_mask::{MaskKind, ResolvedMask};
    use crate::seqdb::{SeqDb, SeqType};
    let open = |name: &str| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/dbmaskp")
            .join(name);
        let name = path.to_str().unwrap();
        let mut db = SeqDb::open(name, SeqType::Protein).unwrap();
        db.read_name_records(name.as_bytes(), true).unwrap();
        db
    };
    // A new masking for each database: the walk is kept with the masking, as one search keeps it.
    let masked = |kind: MaskKind| {
        let mut options = super::super::stage_c_tests::options();
        options.db_mask = Some(ResolvedMask {
            kind,
            id: 21,
            key: None,
            shown: Vec::new(),
            starts: Default::default(),
        });
        options
    };
    for (name, apart) in [("ov_th10", true), ("ov_e3", true), ("ov_e3c", false)] {
        let db = open(name);
        assert_eq!(
            keeps_sequences_apart(&masked(MaskKind::Hard), Some(&db)),
            apart,
            "{name}"
        );
        assert!(!keeps_sequences_apart(&masked(MaskKind::Soft), Some(&db)));
    }
    assert!(!keeps_sequences_apart(&masked(MaskKind::Hard), None));
    assert!(!keeps_sequences_apart(
        &super::super::stage_c_tests::options(),
        Some(&open("ov_th10"))
    ));

    let record = |n| FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence: vec![b'A'; n],
        lowercase_masks: vec![],
        source: None,
        around: None,
        skipped: None,
    };
    let records = [record(4000), record(4000), record(4000), record(100)];
    // The units of `batch_ranges` (the third query brings the first to 10002 letters), and the
    // whole batch where a sequence is kept apart or with `-use_sw_tback`.
    for (name, whole) in [("ov_th10", true), ("ov_e3", true), ("ov_e3c", false)] {
        let db = open(name);
        let options = masked(MaskKind::Hard);
        let want = if whole { vec![0..4] } else { vec![0..3, 3..4] };
        assert_eq!(unit_ranges(&records, &options, Some(&db)), want, "{name}");
        assert_eq!(searches_batch_whole(&options, Some(&db)), whole, "{name}");
    }
    // `-subject` (no database): the units.
    assert_eq!(
        unit_ranges(&records, &masked(MaskKind::Hard), None),
        vec![0..3, 3..4]
    );
    let mut sw = super::super::stage_c_tests::options();
    sw.use_sw_tback = true;
    assert_eq!(unit_ranges(&records, &sw, None), vec![0..4]);
    assert!(searches_batch_whole(&sw, None));
    assert!(unit_ranges(&[], &sw, None).is_empty());
}

use super::*;

fn message(error: anyhow::Error) -> String {
    format!("{error:#}")
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/blastinput/unit_test/blastinput_unit_test.cpp:2760-2763 ParseSequenceRange_EmptyRange
// NCBI unit test (598d8ae6): c++/src/algo/blast/blastinput/unit_test/blastinput_unit_test.cpp:2767-2770 ParseSequenceRange_0BasedRange
// NCBI unit test (598d8ae6): c++/src/algo/blast/blastinput/unit_test/blastinput_unit_test.cpp:2774-2777 ParseSequenceRange_InvalidDelimiter
// NCBI unit test (598d8ae6): c++/src/algo/blast/blastinput/unit_test/blastinput_unit_test.cpp:2781-2788 ParseSequenceRange_IncompleteRange
// NCBI unit test (598d8ae6): c++/src/algo/blast/blastinput/unit_test/blastinput_unit_test.cpp:2792-2798 ParseSequenceRange_InvalidRange
// NCBI unit test (598d8ae6): c++/src/algo/blast/blastinput/unit_test/blastinput_unit_test.cpp:2802-2805 ParseSequenceRange_1BasedRange
#[test]
fn spellings_follow_ncbis_parser() {
    let parse = |text: &str| parse_sequence_range(text, RangeRole::Query);
    assert_eq!(parse("10-20").unwrap(), SequenceRange { from: 9, to: 19 });
    assert_eq!(parse("+1-5").unwrap(), SequenceRange { from: 0, to: 4 });
    assert_eq!(parse("01-05").unwrap(), SequenceRange { from: 0, to: 4 });
    assert_eq!(parse("1-+5").unwrap(), SequenceRange { from: 0, to: 4 });
    assert_eq!(
        parse("1-2147483647").unwrap(),
        SequenceRange {
            from: 0,
            to: 2147483646
        }
    );
    for text in [
        "", "-", "--", "-10", "10-", "10", "1--5", "1-2-3", "-5-10", "3,4",
    ] {
        assert_eq!(
            message(parse(text).unwrap_err()),
            "LOSAT engine error: Invalid specification of query location (Format: start-stop)\n",
            "{text:?}"
        );
    }
    for text in ["0-10", "10-0"] {
        assert!(message(parse(text).unwrap_err())
            .contains("(range elements cannot be less than or equal to 0)"));
    }
    for text in ["1-1", "5-5"] {
        assert!(message(parse(text).unwrap_err()).contains("(range cannot be empty)"));
    }
    assert!(message(parse("20-10").unwrap_err()).contains("(start cannot be larger than stop)"));
    // NCBI BLAST+ 2.17.0 blastn and blastp `-query_loc` / `-subject_loc` (env parity probe
    // h_probe of 2026-11-12, cases `LOC.*`: `a-10`, `1-b`, `1-1.5`): the part that
    // `NStr::StringToInt` cannot convert stops the program with its `CStringException`, the
    // start before the stop.
    let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
    for (text, part, line, kind, pos) in [
        ("a-10", "a", 861, "Int8", 0),
        ("1-b", "b", 861, "Int8", 0),
        ("1-1.5", "1.5", 870, "Int8", 1),
        (" 1-5", " 1", 861, "Int8", 0),
        ("1-5 ", "5 ", 870, "Int8", 1),
        ("1-5a", "5a", 870, "Int8", 1),
        ("1.0-5", "1.0", 870, "Int8", 1),
        ("0x10-0x20", "0x10", 870, "Int8", 1),
        ("0-a", "a", 861, "Int8", 0),
        ("1-2147483648", "2147483648", 640, "int, overflow", 0),
        (
            "1-99999999999999999999",
            "99999999999999999999",
            852,
            "Int8, overflow",
            18,
        ),
    ] {
        let function = if line == 640 {
            "StringToInt"
        } else {
            "StringToInt8"
        };
        let error = parse(text).unwrap_err();
        let native = error
            .downcast_ref::<crate::cli::NativeError>()
            .expect("NCBI's exception");
        assert_eq!(native.exit, 255, "{text:?}");
        assert_eq!(
            native.message,
            // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
            crate::wording::losat_text(&format!(
                "Error: NCBI C++ Exception:\n    T0 \"{root}/src/corelib/ncbistr.cpp\", line {line}: Error: (CStringException::eConvert) ncbi::NStr::{function}() - Cannot convert string '{part}' to {kind} (m_Pos = {pos})\n\n"
            ))
            .as_bytes(),
            "{text:?}"
        );
    }
    assert!(
        message(parse_sequence_range("5-3", RangeRole::Subject).unwrap_err())
            .starts_with("LOSAT engine error: Invalid specification of subject location")
    );
}

#[test]
fn intervals_follow_the_record_length() {
    let range = SequenceRange { from: 9, to: 19 };
    assert_eq!(
        record_interval(&range, 30),
        RecordInterval::Letters { from: 9, to: 19 }
    );
    assert_eq!(
        record_interval(&range, 15),
        RecordInterval::Letters { from: 9, to: 14 }
    );
    assert_eq!(
        record_interval(&range, 10),
        RecordInterval::Letters { from: 9, to: 9 }
    );
    assert_eq!(record_interval(&range, 9), RecordInterval::Empty);
    assert_eq!(record_interval(&range, 8), RecordInterval::PastEnd);
}

#[test]
fn skipped_queries_keep_their_ordinals_and_batches_count_whole_records() {
    let records: Vec<crate::blastinput::fasta_reader::FastaRecord> = [8000usize, 400, 3000, 200]
        .iter()
        .enumerate()
        .map(|(index, &length)| {
            crate::blastinput::fasta_reader::FastaRecord::new(
                format!("Query_{}", index + 1),
                format!("q{index}").as_bytes(),
                &vec![b'A'; length],
            )
        })
        .collect();
    let ranged = cut_queries(&records, &SequenceRange { from: 500, to: 999 });
    let input = &ranged.input;
    assert_eq!(input.ordinals, vec![0, 2]);
    assert_eq!(input.skipped, vec![false, true, false, true]);
    assert_eq!(input.range_length, Some(500));
    assert_eq!(ranged.records[0].seq().len(), 500);
    assert_eq!(input.placements.offset(1), 500);
    assert_eq!(input.placements.length(1, 500), 3000);
    // 8000 letters reach a batch of 5000: the next batch holds q1 (skipped), q2, q3.
    assert_eq!(
        input.batches(5000),
        vec![
            RangedBatch {
                input: 0..1,
                searched: 0..1
            },
            RangedBatch {
                input: 1..4,
                searched: 1..2
            }
        ]
    );
    // A batch of 3000 letters closes after q2; q3 alone is skipped: an empty batch.
    assert_eq!(
        input.batches(3000),
        vec![
            RangedBatch {
                input: 0..1,
                searched: 0..1
            },
            RangedBatch {
                input: 1..3,
                searched: 1..2
            },
            RangedBatch {
                input: 3..4,
                searched: 2..2
            }
        ]
    );
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the input records from the start of a batch to the end of a
// later one have the whole input's batches of those records (each moved by the slice's first
// record and first searched query), with their placements and positions; with `-query_loc`
// (skipped records inside a batch) and without it.
#[test]
fn a_slice_of_whole_batches_has_the_batches_of_the_input() {
    let lengths = [
        8000usize, 400, 3000, 200, 7000, 50, 9000, 600, 600, 1200, 300,
    ];
    let records: Vec<crate::blastinput::fasta_reader::FastaRecord> = lengths
        .iter()
        .enumerate()
        .map(|(index, &length)| {
            crate::blastinput::fasta_reader::FastaRecord::new(
                format!("Query_{}", index + 1),
                format!("q{index}").as_bytes(),
                &vec![b'A'; length],
            )
        })
        .collect();
    let ranged = cut_queries(&records, &SequenceRange { from: 500, to: 999 });
    assert!(ranged.input.skipped.iter().any(|&skipped| skipped));
    for input in [QueryInput::whole(&records), ranged.input] {
        for batch_size in [1u32, 700, 3000, 5000, 10000] {
            let batches = input.batches(batch_size);
            for first in 0..batches.len() {
                for last in first..batches.len() {
                    let records = batches[first].input.start..batches[last].input.end;
                    let searched = batches[first].searched.start..batches[last].searched.end;
                    let slice = input.slice(records.clone());
                    let expected: Vec<RangedBatch> = batches[first..=last]
                        .iter()
                        .map(|batch| RangedBatch {
                            input: batch.input.start - records.start
                                ..batch.input.end - records.start,
                            searched: batch.searched.start - searched.start
                                ..batch.searched.end - searched.start,
                        })
                        .collect();
                    let what = format!("batch size {batch_size}, records {records:?}");
                    assert_eq!(slice.batches(batch_size), expected, "{what}");
                    assert_eq!(slice.ordinals.len(), searched.len(), "{what}");
                    for (index, query) in searched.clone().enumerate() {
                        assert_eq!(
                            slice.ordinal(index) + records.start,
                            input.ordinal(query),
                            "{what}"
                        );
                        assert_eq!(
                            slice.placements.get(index),
                            input.placements.get(query),
                            "{what}"
                        );
                    }
                    assert_eq!(slice.range_length, input.range_length);
                }
            }
        }
    }
}

#[test]
fn frames_are_those_of_the_whole_record() {
    assert_eq!(record_frame(true, 0, 100), 1);
    assert_eq!(record_frame(true, 4, 100), 2);
    assert_eq!(record_frame(false, 99, 100), -1);
    assert_eq!(record_frame(false, 97, 100), -3);
}

#[test]
fn masks_move_into_the_record_and_stop_at_the_interval() {
    // Without a range: as found, other than the whole record.
    assert_eq!(record_mask(10, 20, 100, 0), Some((10, 20)));
    assert_eq!(record_mask(0, 99, 100, 0), None);
    assert_eq!(record_mask(5, 4, 100, 0), None);
    // The interval 501-600 of a record: moved by its start, cut at its end.
    assert_eq!(record_mask(10, 20, 100, 500), Some((510, 520)));
    assert_eq!(record_mask(90, 150, 100, 500), Some((590, 599)));
    assert_eq!(record_mask(0, 99, 100, 500), None);
    assert_eq!(record_mask(100, 120, 100, 500), None);
}

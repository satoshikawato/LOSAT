use super::*;
use std::io::Cursor;

/// A delivered batch: the ids, letters and placements of its records.
type Batch = Vec<(String, Vec<u8>, RecordPlacement)>;

fn summary(records: &[FastaRecord]) -> Batch {
    records
        .iter()
        .map(|record| {
            (
                record.internal_id.clone(),
                record.sequence.clone(),
                record.placement(),
            )
        })
        .collect()
}

/// The batches a query file of `data` gives with `range` and `batch_size`, read by
/// `parse_fasta_stream_keeping`; and the records the reader keeps.
fn batches_kept(
    data: &[u8],
    range: Option<&SequenceRange>,
    batch_size: usize,
) -> (Vec<Batch>, Vec<FastaRecord>) {
    let mut delivered = Vec::new();
    let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
    stream.keep_skipped = true;
    let kept = parse_fasta_stream_keeping(
        &mut stream,
        false,
        true,
        batch_size,
        false,
        range,
        &mut |records| {
            delivered.push(summary(records));
            Ok(())
        },
        &mut |_| Ok(()),
    )
    .unwrap();
    (delivered, kept)
}

/// `batches_kept` with a range.
fn batches(
    data: &[u8],
    range: &SequenceRange,
    batch_size: usize,
) -> (Vec<Batch>, Vec<FastaRecord>) {
    batches_kept(data, Some(range), batch_size)
}

/// The batches of `batches_kept` read by `read_fasta_stream_batches`, which keeps no record.
fn batches_unkept(data: &[u8], range: Option<&SequenceRange>, batch_size: usize) -> Vec<Batch> {
    let mut delivered = Vec::new();
    read_fasta_stream_batches(
        &mut FastaStream::new(Cursor::new(data.to_vec())),
        false,
        true,
        batch_size,
        false,
        range,
        &mut |records| {
            delivered.push(summary(records));
            Ok(())
        },
        &mut |_| Ok(()),
    )
    .unwrap();
    delivered
}

// The read that drops each batch once delivered (the search's queries and the reads of them
// before the search) delivers the batches of the read that keeps the records: the same records,
// cut the same way, with the same numbers after the dropped batches and the skipped records.
#[test]
fn a_read_that_keeps_no_record_delivers_the_same_batches() {
    let range = SequenceRange { from: 4, to: 7 };
    let data =
        b">a\nACGTACGTAC\n>b\nACG\n>c\nTTTTGGGG\n>d\nACGT\n>e\nACGTA\n>f\nAC\n>g\nACGTACGT\n";
    for batch_size in [1, 4, 10, 15, 21, usize::MAX] {
        for range in [Some(&range), None] {
            let (kept_batches, kept) = batches_kept(data, range, batch_size);
            assert_eq!(
                batches_unkept(data, range, batch_size),
                kept_batches,
                "batch size {batch_size}, range {}",
                range.is_some()
            );
            assert_eq!(kept_batches.iter().map(Vec::len).sum::<usize>(), kept.len());
        }
    }
    let (delivered, kept) = batches_kept(data, None, 10);
    assert_eq!(kept.len(), 7);
    assert_eq!(delivered.len(), 4);
    assert_eq!(delivered[3][0].0, "Query_7");
}

// The read that keeps no record frees each batch's records once it has delivered the batch,
// so that it holds one reader batch, not the file; the read that keeps them holds every
// record to its end. Under `-parse_deflines` (without the search's ids) a record is the only
// owner of its Seq-ids, so a weak reference to them shows whether the record is still held.
#[test]
fn a_read_that_keeps_no_record_frees_each_batch_once_delivered() {
    use crate::blastinput::fasta_reader::believed::BelievedIds;
    use std::sync::{Arc, Weak};
    let data = b">lcl|q1\nACGT\n>lcl|q2\nACGT\n>lcl|q3\nACGT\n";
    for keep in [false, true] {
        let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
        stream.parse_deflines = true;
        let mut ids: Vec<Weak<BelievedIds>> = Vec::new();
        // At each delivery, the records of the earlier batches still held.
        let mut held_before = Vec::new();
        let mut deliver = |records: &[FastaRecord]| {
            held_before.push(ids.iter().filter(|id| id.upgrade().is_some()).count());
            for record in records {
                let believed = record
                    .source
                    .as_ref()
                    .and_then(|source| source.believed.as_ref())
                    .expect("parsed Seq-ids");
                ids.push(Arc::downgrade(believed));
            }
            Ok(())
        };
        let kept = if keep {
            parse_fasta_stream_keeping(
                &mut stream,
                false,
                false,
                1,
                false,
                None,
                &mut deliver,
                &mut |_| Ok(()),
            )
            .unwrap()
        } else {
            read_fasta_stream_batches(
                &mut stream,
                false,
                false,
                1,
                false,
                None,
                &mut deliver,
                &mut |_| Ok(()),
            )
            .unwrap();
            Vec::new()
        };
        let held_after = ids.iter().filter(|id| id.upgrade().is_some()).count();
        if keep {
            assert_eq!(held_before, vec![0, 1, 2]);
            assert_eq!((kept.len(), held_after), (3, 3));
        } else {
            assert_eq!(held_before, vec![0, 0, 0]);
            assert_eq!(held_after, 0);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:144-168
// A batch counts the whole length of each searched record; a skipped record (its interval
// starts more than one letter past its end) counts nothing, is not delivered and keeps its
// number; the letters delivered are the interval's.
#[test]
fn a_range_cuts_the_records_and_skips_the_short_ones() {
    let range = SequenceRange { from: 4, to: 7 };
    let data = b">a\nACGTACGTAC\n>b\nACG\n>c\nTTTTGGGG\n>d\nACGT\n>e\nACGTA\n";
    let (delivered, kept) = batches(data, &range, 15);
    let place = |offset, length| RecordPlacement { offset, length };
    assert_eq!(
        delivered,
        vec![
            vec![
                ("Query_1".to_string(), b"ACGT".to_vec(), place(4, 10)),
                ("Query_3".to_string(), b"GGGG".to_vec(), place(4, 8)),
            ],
            // `d` starts just past its end: an interval without letters; `e` keeps one letter.
            vec![
                ("Query_4".to_string(), Vec::new(), place(4, 4)),
                ("Query_5".to_string(), b"A".to_vec(), place(4, 5)),
            ],
        ]
    );
    assert_eq!(kept.len(), 4);
    // The skipped `b` stays in its batch's scope, on the batch's first record (review 2 of
    // line DBID loadord, H2; blast_fasta_input.cpp:395 adds it to the scope first), in both
    // readers; the second batch skipped none.
    let ids = |records: &[FastaRecord]| -> Vec<String> {
        shown_ids(records).map(|(id, _)| id.to_owned()).collect()
    };
    assert_eq!(ids(&kept[..2]), ["Query_1", "Query_3", "Query_2"]);
    assert_eq!(ids(&kept[2..]), ["Query_4", "Query_5"]);
    let unkept = |keep_skipped: bool| {
        let mut unkept = Vec::new();
        let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
        stream.keep_skipped = keep_skipped;
        read_fasta_stream_batches(
            &mut stream,
            false,
            true,
            15,
            false,
            Some(&range),
            &mut |records| {
                unkept.push(ids(records));
                Ok(())
            },
            &mut |_| Ok(()),
        )
        .unwrap();
        unkept
    };
    assert_eq!(
        unkept(true),
        [
            vec!["Query_1", "Query_3", "Query_2"],
            vec!["Query_4", "Query_5"]
        ]
    );
    // A run whose prefetch cannot stop keeps none of them (review 5 of line DBID loadord,
    // L2); the numbers stay.
    assert_eq!(
        unkept(false),
        [vec!["Query_1", "Query_3"], vec!["Query_4", "Query_5"]]
    );
    // The scopes of the queries take the records in the order read, the skipped `b` with its
    // whole letters between `a` and `c` (`scope_records`; line WQA D10 round 6).
    let order: Vec<(usize, Option<usize>)> = scope_records(&kept[..2])
        .into_iter()
        .map(|(_, length, place)| (length, place))
        .collect();
    assert_eq!(order, [(10, Some(0)), (3, None), (8, Some(1))]);
}

// A batch of skipped records alone is delivered empty (NCBI's query factory refuses it), also
// at the end of the file.
#[test]
fn a_last_batch_of_skipped_records_is_delivered_empty() {
    let range = SequenceRange { from: 4, to: 7 };
    let (delivered, _) = batches(b">a\nACGTACGTAC\n>b\nACG\n>c\nA\n", &range, 10);
    assert_eq!(delivered.len(), 2);
    assert_eq!(delivered[0].len(), 1);
    assert!(delivered[1].is_empty());
    let (delivered, kept) = batches(b">a\nACG\n", &range, 10);
    assert!(kept.is_empty());
    assert_eq!(delivered, vec![Vec::new()]);
}

// A record whose letter before the interval is lower case keeps its whole letters after the
// cut (that letter's run, which a split query's chunks restrict), whether the run ends there
// or goes on into the interval; a record with runs only elsewhere, or none, keeps only its
// interval.
#[test]
fn a_cut_record_keeps_its_whole_letters_for_a_lower_case_letter_before_it() {
    let range = SequenceRange { from: 4, to: 7 };
    let data = b">a\nACGtACGTAC\n>b\nAcgtaCGTAC\n>c\naCGTaCGTac\n>d\nACGTACGTAC\n";
    let (_, kept) = batches(data, &range, 100);
    assert_eq!(kept[0].sequence, b"ACGT");
    assert!(kept[0].lowercase_masks.is_empty());
    assert_eq!(kept[0].archive_letters(), b"ACGtACGTAC");
    assert_eq!(kept[0].archive_lowercase_runs(), vec![(3, 3)]);
    assert_eq!(kept[1].lowercase_masks, vec![(0, 0)]);
    assert_eq!(kept[1].archive_letters(), b"AcgtaCGTAC");
    assert_eq!(kept[2].lowercase_masks, vec![(0, 0)]);
    assert_eq!(kept[2].archive_letters(), b"ACGT");
    assert_eq!(kept[3].archive_letters(), b"ACGT");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1041-1046
#[test]
fn a_cut_clips_and_moves_the_lower_case_masks() {
    let mut record = parse_fasta(b">q\nacgTTTTaaaCCgg\n", false, true)
        .unwrap()
        .remove(0);
    assert_eq!(record.lowercase_masks, vec![(0, 2), (7, 9), (12, 13)]);
    let mut kept = record.clone();
    record.cut_in_place(2, 9, false);
    assert_eq!(record.sequence, b"GTTTTAA");
    assert_eq!(record.lowercase_masks, vec![(0, 0), (5, 6)]);
    assert_eq!(
        record.placement(),
        RecordPlacement {
            offset: 2,
            length: 14
        }
    );
    // Only a record that keeps its whole letters has them (and its own runs) after the cut.
    assert_eq!(record.archive_letters(), b"GTTTTAA");
    assert_eq!(record.archive_lowercase_runs(), vec![(0, 0), (5, 6)]);
    kept.cut_in_place(2, 9, true);
    assert_eq!(kept.sequence, record.sequence);
    assert_eq!(kept.lowercase_masks, record.lowercase_masks);
    assert_eq!(kept.placement(), record.placement());
    assert_eq!(kept.archive_letters(), b"acgTTTTaaaCCgg");
    assert_eq!(
        kept.archive_lowercase_runs(),
        vec![(0, 2), (7, 9), (12, 13)]
    );
    let whole = parse_fasta(b">q\nACGT\n", false, true).unwrap().remove(0);
    assert_eq!(
        whole.placement(),
        RecordPlacement {
            offset: 0,
            length: 4
        }
    );
}

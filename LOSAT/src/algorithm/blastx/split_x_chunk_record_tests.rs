use super::{
    prepare_queries, split_queries, x_chunk_record, x_chunk_record_from_copy, FastaRecord,
    DEFAULT_SIZES,
};
use crate::algorithm::blastx::input::InputWarning;
use crate::blastinput::fasta_reader::RecordSource;

/// A query of 30,000 letters of no repeat in its frames, with masks at the chunks' borders.
fn whole() -> FastaRecord {
    FastaRecord {
        internal_id: "Query_1".to_string(),
        title: b"long query".to_vec(),
        sequence: (0..30_000u32)
            .map(|k| b"ACGT"[(k.wrapping_mul(2_654_435_761) >> 30) as usize])
            .collect(),
        lowercase_masks: vec![(0, 9), (100, 299), (9_990, 10_010), (29_990, 29_999)],
        warnings: vec![InputWarning {
            line: 3,
            kind: "defline_id",
            positions: vec![7],
            text: b"a warning".to_vec(),
        }],
        source: Some(Box::new(RecordSource {
            defline: b">long query".to_vec(),
            ..RecordSource::default()
        })),
        around: None,
        skipped: None,
    }
}

// LOSAT_X_BXCHUNKSLICE: the record without the copy equals the one from the copy in every
// field, at ranges that hold, cross and miss the masks.
#[test]
fn a_chunks_record_equals_the_one_cut_from_a_copy() {
    let whole = whole();
    for (start, end) in [
        (0, 10_002),
        (9_705, 19_707),
        (10, 100),
        (299, 300),
        (300, 9_000),
        (29_000, 30_000),
    ] {
        let new = x_chunk_record(&whole, start, end);
        let old = x_chunk_record_from_copy(&whole, start, end);
        assert_eq!(format!("{new:?}"), format!("{old:?}"), "{start}..{end}");
        assert_eq!(new.sequence.len(), end - start);
    }
}

// LOSAT_X_BXCHUNKSLICE: `split_queries` does split this query, and each chunk holds what a
// record cut from a copy prepares, masks included.
#[test]
fn a_split_querys_chunks_are_prepared_as_from_a_copy() {
    let whole = whole();
    let options = super::super::stage_c_tests::options();
    let records = std::slice::from_ref(&whole);
    let full = prepare_queries(records, &options).unwrap();
    let chunks = split_queries(records, &full, &options, DEFAULT_SIZES).unwrap();
    assert!(chunks.len() > 1, "{} chunks", chunks.len());
    for chunk in &chunks {
        let (start, end) = (chunk.range.start, chunk.range.end);
        let old = x_chunk_record_from_copy(&whole, start, end);
        let old = prepare_queries(&[old], &options).unwrap();
        let new = &chunk.prepared;
        assert_eq!(new.original_lengths, old.original_lengths, "{start}..{end}");
        assert_eq!(new.sequence_start, old.sequence_start, "{start}..{end}");
        assert_eq!(
            new.sequence_start_nomask, old.sequence_start_nomask,
            "{start}..{end}"
        );
        assert_eq!(new.lookup_segments, old.lookup_segments, "{start}..{end}");
    }
}

/// `whole` cut `offset` letters into its record, keeping the record's letters: `before` lower
/// case letters end on the letter before the interval.
fn cut_whole(mut whole: FastaRecord, offset: usize, before: usize) -> FastaRecord {
    let mut uncut = vec![b'A'; offset];
    uncut[offset - before..].make_ascii_lowercase();
    uncut.extend_from_slice(&whole.sequence);
    for &(left, right) in &whole.lowercase_masks {
        uncut[offset + left as usize..=offset + right as usize].make_ascii_lowercase();
    }
    let source = whole.source.as_deref_mut().unwrap();
    source.cut_from = offset;
    source.uncut_length = uncut.len();
    source.uncut = Some(uncut);
    whole
}

// A chunk restricts the runs it reads from its letters in the record and the one before
// them as it would restrict all the runs of the record, also runs that cross the ends of those
// letters, at offsets 0, 1 and 2,500.
#[test]
fn a_chunks_runs_read_near_it_restrict_as_all_the_records_runs() {
    use crate::algorithm::blastn::query_split::{restrict_record_masks, ChunkQuery};
    use crate::utils::dust::MaskedInterval;
    let mut whole = whole();
    whole
        .lowercase_masks
        .extend([(9_700, 9_710), (19_405, 19_415), (19_700, 19_720)]);
    whole.lowercase_masks.sort_unstable();
    for (offset, before) in [(0, 0), (1, 1), (2_500, 0), (2_500, 1), (2_500, 3)] {
        let cut = cut_whole(whole.clone(), offset, before);
        let letters = cut.source.as_deref().unwrap().uncut.as_deref().unwrap();
        let runs: Vec<_> = crate::algorithm::common::masks::lowercase_runs(letters)
            .into_iter()
            .map(|(a, b)| MaskedInterval::new(a as usize, b as usize + 1))
            .collect();
        for (start, end) in [
            (0, 10_002),
            (9_705, 19_707),
            (19_410, 30_000),
            (99, 101),
            (299, 300),
            (9_701, 9_702),
        ] {
            let part = ChunkQuery {
                query: 0,
                from: start,
                to: end,
            };
            let all: Vec<_> = restrict_record_masks(&runs, &part, offset)
                .into_iter()
                .map(|mask| (mask.start as i32, mask.end as i32 - 1))
                .collect();
            assert_eq!(
                x_chunk_record(&cut, start, end).lowercase_masks,
                all,
                "offset {offset} before {before} {start}..{end}"
            );
        }
    }
}

// `-query_loc`: a record cut 2,500 letters into its record keeps a run in a chunk only in
// `[start + 2,500, end - 1]` of the record (the run in the last chunk's last letters is
// dropped, those inside stay), the run that ends on the letter before the interval masks the
// first chunk's first letter (a record without its whole letters cannot), and no chunk copies
// the record's whole letters.
#[test]
fn a_cut_records_chunks_restrict_the_runs_in_the_record() {
    let mut whole = whole();
    whole.lowercase_masks.remove(0);
    let plain: Vec<_> = [(9_705, 19_707), (29_000, 30_000)]
        .iter()
        .map(|&(start, end)| x_chunk_record(&whole, start, end).lowercase_masks)
        .collect();
    let whole = cut_whole(whole, 2_500, 3);
    let mut runs_cut_only = whole.clone();
    runs_cut_only.source.as_deref_mut().unwrap().uncut = None;
    let first = x_chunk_record(&whole, 0, 10_002);
    assert_eq!(first.lowercase_masks, vec![(0, 0), (100, 300)]);
    assert_eq!(
        x_chunk_record(&runs_cut_only, 0, 10_002).lowercase_masks,
        vec![(100, 300)]
    );
    let inside = x_chunk_record(&whole, 9_705, 19_707);
    assert_eq!(inside.lowercase_masks, plain[0]);
    assert!(!plain[1].is_empty());
    let last = x_chunk_record(&whole, 29_000, 30_000);
    assert!(last.lowercase_masks.is_empty());
    for chunk in [first, inside, last] {
        let source = chunk.source.as_deref().unwrap();
        assert!(source.uncut.is_none());
        assert_eq!(
            (
                source.cut_from,
                source.uncut_length,
                chunk.placement().offset
            ),
            (2_500, 32_500, 2_500)
        );
    }
}

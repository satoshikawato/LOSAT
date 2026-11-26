use super::*;

fn collection(ranges: &[(u32, u32)]) -> RangeCollection {
    let mut out = RangeCollection::default();
    for &(from, to) in ranges {
        out.combine_with(SeqRange::new(from, to));
    }
    out
}

fn pairs(collection: &RangeCollection) -> Vec<(u32, u32)> {
    collection.0.iter().map(|r| (r.from, r.to())).collect()
}

#[test]
fn range_collection_merges_overlapping_and_touching_ranges() {
    let merged = collection(&[(10, 19), (30, 39), (20, 25), (0, 3)]);
    assert_eq!(pairs(&merged), [(0, 3), (10, 25), (30, 39)]);
    assert_eq!(merged.covered_length(), 4 + 16 + 10);
    let spanning = collection(&[(10, 19), (30, 39), (5, 50)]);
    assert_eq!(pairs(&spanning), [(5, 50)]);
    let zero = collection(&[(5, 9), (0, 4)]);
    assert_eq!(pairs(&zero), [(0, 9)]);
}

// Ranges out of the order of `to_open` (a range whose `to_open` wrapped to 0 after larger
// ones): the positions that `std::lower_bound` returned for them (`observe_lower_bound one`,
// `LOSAT/tests/unit/lower_bound_observed/README.md`). The first five are not the first range
// not less than the value (at 1, 1, 1, 0 and 2). The last two also tell the observed search
// from other searches that give the first five: `slice::partition_point`, for one, returns 3
// for 15.
#[test]
fn lower_bound_of_unordered_ranges_is_the_observed_one() {
    let ranges = |to_opens: &[u32]| {
        RangeCollection(
            to_opens
                .iter()
                .map(|&to_open| SeqRange {
                    from: to_open,
                    to_open,
                })
                .collect(),
        )
    };
    assert_eq!(ranges(&[40, 50, 0, 60]).lower_bound(0, 45), 3);
    assert_eq!(ranges(&[10, 30, 0, 40, 50]).lower_bound(0, 25), 3);
    assert_eq!(ranges(&[10, 30, 0, 40, 50]).lower_bound(1, 25), 3);
    assert_eq!(ranges(&[100, 0, 0, 0, 0, 0, 0, 0]).lower_bound(0, 50), 8);
    assert_eq!(ranges(&[50, 0, 120, 4, 90, 200]).lower_bound(1, 100), 5);
    assert_eq!(ranges(&[0, 0, 10, 20, 0]).lower_bound(0, 5), 2);
    assert_eq!(ranges(&[0, 0, 10, 20, 0]).lower_bound(0, 15), 5);
}

#[test]
fn range_collection_subtracts_as_ncbi() {
    let mut parts = RangeCollection(vec![SeqRange::new(0, 99)]);
    parts.subtract(SeqRange::new(10, 19));
    assert_eq!(pairs(&parts), [(0, 9), (20, 99)]);
    parts.subtract(SeqRange::new(90, 120));
    assert_eq!(pairs(&parts), [(0, 9), (20, 89)]);
    parts.subtract(SeqRange::new(0, 30));
    assert_eq!(pairs(&parts), [(31, 89)]);
    parts.subtract(SeqRange::new(0, 200));
    assert!(parts.0.is_empty());
}

fn hsp(
    q: (usize, usize),
    s: (usize, usize),
    minus: bool,
    ops: Option<Vec<GapEditOp>>,
) -> BlastnHsp {
    BlastnHsp {
        identity: 0.0,
        length: 0,
        mismatch: 0,
        gapopen: 0,
        q_start: q.0,
        q_end: q.1,
        s_start: s.0,
        s_end: s.1,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: if minus { -1 } else { 1 },
        query_length: 1000,
        q_idx: 0,
        s_idx: 0,
        raw_score: 0,
        internal_q_offset_0: 0,
        internal_q_end_0: (q.1 - q.0 + 1) as i32,
        internal_s_offset_0: 0,
        internal_s_end_0: s.0.max(s.1) - s.0.min(s.1) + 1,
        internal_query_context_offset: 0,
        gap_info: ops,
        num_positives: 0,
    }
}

#[test]
fn query_positions_follow_the_columns() {
    // Query 11..=20 against subject 101..=111 with a query gap after 4 columns.
    let gapped = hsp(
        (11, 20),
        (101, 111),
        false,
        Some(vec![
            GapEditOp::Sub(4),
            GapEditOp::Del(1),
            GapEditOp::Sub(6),
        ]),
    );
    assert_eq!(query_position(&gapped, 100, 0, false), 10);
    assert_eq!(query_position(&gapped, 103, 0, false), 13);
    assert_eq!(query_position(&gapped, 104, 0, false), -1);
    assert_eq!(query_position(&gapped, 105, 0, false), 14);
    assert_eq!(query_position(&gapped, 110, 0, false), 19);
    // A minus-strand HSP: the query descends along the subject.
    let minus = hsp((11, 20), (110, 101), true, None);
    assert_eq!(query_position(&minus, 100, 0, true), 19);
    assert_eq!(query_position(&minus, 109, 0, true), 10);
}

#[test]
fn uniq_coverage_skips_subject_letters_covered_before() {
    // Two HSPs on disjoint subject ranges: both query ranges count.
    let a = hsp((1, 100), (1, 100), false, None);
    let b = hsp((201, 300), (201, 300), false, None);
    assert_eq!(uniq_seq_coverage(&[a.clone(), b], 0, true).unwrap(), 200);
    // A repeat of the query on the same subject letters counts once.
    let repeat = hsp((401, 500), (1, 100), false, None);
    assert_eq!(
        uniq_seq_coverage(&[a.clone(), repeat], 0, true).unwrap(),
        100
    );
    // A partly covered HSP: only its uncovered subject part maps to query letters.
    let partial = hsp((501, 600), (51, 150), false, None);
    assert_eq!(uniq_seq_coverage(&[a, partial], 0, true).unwrap(), 100 + 50);
    assert_eq!(uniq_subject_coverage(&[], 100, 0, true), -1);
}

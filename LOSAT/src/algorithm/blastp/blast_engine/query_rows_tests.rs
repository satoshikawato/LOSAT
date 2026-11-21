use super::*;
use crate::common::Hit;

/// A row of query `q_idx` against subject `s_idx`.
fn row(q_idx: u32, s_idx: u32) -> PairwiseHit {
    PairwiseHit::from(Hit {
        identity: 100.0,
        length: 10,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 10,
        s_start: 1,
        s_end: 10,
        e_value: 1e-5,
        bit_score: 25.0,
        num_ident: 10,
        query_frame: 0,
        query_length: 30,
        q_idx,
        s_idx,
        raw_score: 50,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 10,
    })
}

/// The (query, subject) pairs of `rows`.
fn pairs(rows: &[PairwiseHit]) -> Vec<(u32, u32)> {
    rows.iter()
        .map(|row| (row.hit.q_idx, row.hit.s_idx))
        .collect()
}

/// Queries 0 and 3 without rows; 1 with two, 2 with one, 4 with three.
fn all_rows() -> Vec<PairwiseHit> {
    vec![
        row(1, 5),
        row(1, 2),
        row(2, 7),
        row(4, 0),
        row(4, 9),
        row(4, 1),
    ]
}

/// The rows of the queries of `window` among `all_rows`, as the call's producer makes them.
fn rows_of(window: Range<usize>) -> Result<Vec<PairwiseHit>> {
    Ok(all_rows()
        .into_iter()
        .filter(|row| window.contains(&(row.hit.q_idx as usize)))
        .collect())
}

#[test]
fn each_query_gets_its_rows_from_windows_or_from_the_rows_made_up_front() {
    let all = all_rows();
    let expected: Vec<Vec<(u32, u32)>> = (0..5)
        .map(|q_idx| {
            pairs(
                &all.iter()
                    .filter(|row| row.hit.q_idx as usize == q_idx)
                    .cloned()
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    for windows in [
        vec![0..5],
        vec![0..2, 2..5],
        vec![0..1, 1..2, 2..3, 3..4, 4..5],
        vec![0..3, 3..5],
    ] {
        for up_front in [None, Some(all.as_slice())] {
            // Every query asked.
            let mut rows = QueryRows::new(windows.clone(), up_front, &rows_of);
            for (q_idx, expected) in expected.iter().enumerate() {
                assert_eq!(&pairs(rows.query(q_idx).unwrap()), expected, "{windows:?}");
            }
            // Only the queries with rows (BLAST XML, the tabular report).
            let mut rows = QueryRows::new(windows.clone(), up_front, &rows_of);
            for q_idx in [1, 2, 4] {
                assert_eq!(pairs(rows.query(q_idx).unwrap()), expected[q_idx]);
            }
        }
    }
}

#[test]
fn a_query_asked_again_or_past_the_windows_is_an_error() {
    let mut rows = QueryRows::new(vec![0..2, 2..5], None, &rows_of);
    rows.query(2).unwrap();
    assert!(rows.query(1).is_err());
    let mut rows = QueryRows::new(vec![0..2], None, &rows_of);
    rows.query(1).unwrap();
    assert!(rows.query(1).is_err());
    assert!(rows.query(4).is_err());
}

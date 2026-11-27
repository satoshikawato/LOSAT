use super::*;
use crate::seqid::model::Dbtag;

fn domain(tag: i32) -> ParsedSeqId {
    ParsedSeqId::General(Dbtag {
        db: b"CDD".to_vec(),
        tag: ObjectId::Id(tag),
    })
}

/// A hit of the domain `tag` with segments `(query_start, subject_start, len)`.
fn hit(tag: i32, e_value: f64, segments: &[(i32, i32, i32)]) -> Hit {
    hit_of_segments(domain(tag), tag as u32, e_value, segments.iter().copied())
}

fn spans(hit: &Hit) -> Vec<(i32, i32, i32, i32)> {
    hit.segments
        .iter()
        .map(|s| (s.query.from, s.query.to, s.subject.from, s.subject.to))
        .collect()
}

#[test]
fn the_range_of_a_segment_ends_at_start_plus_length() {
    // NCBI's range of a segment ends one past its last column: a segment of 10 columns at 5 is 5..=15.
    let h = hit(1, 1e-5, &[(5, 20, 10)]);
    assert_eq!(spans(&h), [(5, 15, 20, 30)]);
}

#[test]
fn the_inclusion_threshold_keeps_the_hits_from_minus_one_below_it() {
    let hits = vec![
        hit(1, -2.0, &[(0, 0, 5)]),
        hit(2, -1.0, &[(0, 0, 5)]),
        hit(3, 0.01, &[(0, 0, 5)]),
        hit(4, 0.05, &[(0, 0, 5)]),
        hit(5, 0.5, &[(0, 0, 5)]),
    ];
    let kept: Vec<u32> = included_hits(hits, 0.05).iter().map(|h| h.oid).collect();
    assert_eq!(kept, [2, 3]);
}

#[test]
fn seq_ids_order_local_before_general_and_integers_before_strings() {
    let local = ParsedSeqId::Local(ObjectId::Id(9));
    assert_eq!(compare_ordered(&local, &domain(1)), Ordering::Less);
    assert_eq!(compare_ordered(&domain(2), &domain(10)), Ordering::Less);
    let text = ParsedSeqId::General(Dbtag {
        db: b"cdd".to_vec(),
        tag: ObjectId::Str(b"abc".to_vec()),
    });
    // The database name compares without case; an integer tag sorts before a string.
    assert_eq!(compare_ordered(&domain(5), &text), Ordering::Less);
    // A string that reads as an integer is the integer.
    let numeric = ParsedSeqId::General(Dbtag {
        db: b"CDD".to_vec(),
        tag: ObjectId::Str(b"7".to_vec()),
    });
    assert_eq!(compare_ordered(&numeric, &domain(7)), Ordering::Equal);
}

#[test]
fn a_hit_to_another_part_of_a_domain_is_kept_whole() {
    // The same domain, different subject ranges: both stay.
    let hits = vec![hit(1, 1e-9, &[(0, 0, 10)]), hit(1, 1e-5, &[(30, 40, 10)])];
    let kept = remove_multiple_cd_hits(hits).expect("no tie");
    assert_eq!(kept.len(), 2);
    assert_eq!(spans(&kept[1]), [(30, 40, 40, 50)]);
}

#[test]
fn a_later_hit_loses_the_query_range_of_the_better_one_on_the_same_domain_positions() {
    // The better hit covers domain columns 0..=10 at query 0..=10; the worse covers the same
    // domain columns at query 20..=30: the subject ranges intersect, but their query ranges
    // do not, so nothing is subtracted.
    let hits = vec![hit(1, 1e-9, &[(0, 0, 10)]), hit(1, 1e-5, &[(20, 0, 10)])];
    let kept = remove_multiple_cd_hits(hits).expect("no tie");
    assert_eq!(kept.len(), 2);
    assert_eq!(spans(&kept[1]), [(20, 30, 0, 10)]);
    // A worse hit inside the better one's query range, on the same domain columns, is dropped.
    let hits = vec![hit(1, 1e-9, &[(0, 0, 20)]), hit(1, 1e-5, &[(5, 5, 10)])];
    let kept = remove_multiple_cd_hits(hits).expect("no tie");
    assert_eq!(kept.len(), 1);
    assert_eq!(kept[0].e_value, 1e-9);
}

#[test]
fn hits_of_different_domains_do_not_trim_each_other() {
    let hits = vec![hit(2, 1e-5, &[(0, 0, 20)]), hit(1, 1e-9, &[(0, 0, 20)])];
    let kept = remove_multiple_cd_hits(hits).expect("no tie");
    // Sorted by domain.
    assert_eq!(kept.iter().map(|h| h.oid).collect::<Vec<_>>(), [1, 2]);
    assert_eq!(spans(&kept[0]), spans(&kept[1]));
}

#[test]
fn hits_of_one_domain_with_one_evalue_keep_the_order_of_the_search() {
    // Two tied hits (the same domain and E-value, disjoint query ranges) and a worse hit of the
    // domain that each of them trims: the cuts do not commute, so the order of the tie decides
    // the worse hit's segments. LOSAT keeps the order of the domain search's results (Rust's
    // stable sort, the Owner's choice of 2026-11-24); NCBI's C++ library's sort may not.
    let t1 = hit(1, 1e-5, &[(5, 5, 7)]);
    let t2 = hit(1, 1e-5, &[(20, 12, 10)]);
    let worse = hit(1, 1e-3, &[(0, 0, 10), (20, 12, 10), (40, 24, 10)]);
    let kept = remove_multiple_cd_hits(vec![t1.clone(), t2.clone(), worse.clone()]).expect("kept");
    assert_eq!(kept.len(), 3);
    assert_eq!(spans(&kept[0]), spans(&t1));
    assert_eq!(spans(&kept[2]), [(0, 5, 0, 5)]);
    let kept = remove_multiple_cd_hits(vec![t2.clone(), worse.clone(), t1.clone()]).expect("kept");
    assert_eq!(spans(&kept[0]), spans(&t2));
    assert_eq!(spans(&kept[2]), [(0, 5, 0, 5), (40, 50, 24, 34)]);
    // Tied hits with overlapping query ranges are not refused either; they too keep the search's
    // order.
    let a = hit(1, 1e-5, &[(0, 0, 20)]);
    let b = hit(1, 1e-5, &[(10, 50, 20)]);
    let kept = remove_multiple_cd_hits(vec![a.clone(), b.clone()]).expect("kept");
    assert_eq!(kept, [a.clone(), b.clone()]);
    let kept = remove_multiple_cd_hits(vec![b.clone(), a.clone()]).expect("kept");
    assert_eq!(kept, [b, a]);
}

#[test]
fn an_id_that_is_not_local_or_general_is_refused() {
    let genbank = ParsedSeqId::Genbank(Default::default());
    let hits = vec![
        hit_of_segments(genbank.clone(), 0, 1e-5, [(0, 0, 5)]),
        hit_of_segments(genbank, 1, 1e-6, [(0, 0, 5)]),
    ];
    assert!(matches!(
        remove_multiple_cd_hits(hits),
        Err(CddRefusal::IdKind(_))
    ));
}

#[test]
fn a_single_hit_is_not_trimmed() {
    let hits = vec![hit(1, 1e-5, &[(0, 0, 10), (12, 11, 5)])];
    let kept = remove_multiple_cd_hits(hits.clone()).expect("one hit");
    assert_eq!(kept, hits);
}

#[test]
fn intersected_segments_with_equal_subject_start_keep_stable_order() {
    // Adjacent subject ranges [0, 10] and [10, 20] (a gap in the query of the kept hit) intersect
    // with a segment covering [0, 20], producing pieces [0, 10] and [10, 10] / [10, 20] that tie
    // on `subject.from == 10`. Rust's stable sort keeps their intersection order (Owner
    // Q-DL-SEGTIE=A, 2026-11-26; PD-LOSAT-NCBI-DEFECTS 1.123).
    let mut target = hit(1, 1e-3, &[(0, 0, 20)]);
    let ranges = vec![Range { from: 0, to: 10 }, Range { from: 10, to: 20 }];
    target.intersect_with_ranges(&ranges, true);
    assert_eq!(spans(&target), [(0, 10, 0, 10), (10, 20, 10, 20)]);
    // Explicit point-segment tie at subject.from = 10: [10, 10] followed by [10, 15].
    let mut target2 = hit(1, 1e-3, &[(0, 0, 15)]);
    let ranges2 = vec![
        Range { from: 0, to: 10 },
        Range { from: 10, to: 10 },
        Range { from: 10, to: 15 },
    ];
    target2.intersect_with_ranges(&ranges2, true);
    assert_eq!(
        spans(&target2),
        [(0, 10, 0, 10), (10, 10, 10, 10), (10, 15, 10, 15)]
    );
}

/// The rows of `place_cells` as text: `.` an unaligned cell, else the column's index.
fn row_text(rows: &[Vec<CdCell>]) -> Vec<String> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|&cell| {
                    if cell == CdCell::UNALIGNED {
                        ".".to_string()
                    } else {
                        (0..64)
                            .find(|&index| cell == CdCell::aligned(index))
                            .expect("a column below 64")
                            .to_string()
                    }
                })
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect()
}

#[test]
fn cells_of_a_range_from_1_sit_in_their_rows() {
    // Two hits on a range of 5 letters that starts at 1 (`range_start` 0): every record position
    // is inside its row, as without a range.
    let segments = [(0, 1, 3, 0), (1, 0, 2, 3)];
    let rows = place_cells(&segments, 2, 5, 0).expect("in its rows");
    assert_eq!(row_text(&rows), [". 0 1 2 .", "3 4 . . ."]);
}

#[test]
fn cells_past_a_row_spill_into_the_next_hits_rows_and_later_hits_win() {
    // L = 5, a range at record positions 3..8 (`range_start` 3): hit 0's segment at 6..8 writes
    // flat cells 6 and 7 (NCBI's index h * L + q), row 1's positions 1 and 2; hit 1's at 3 is flat
    // 8, row 1's position 3.
    let segments = [(0, 6, 2, 0), (1, 3, 1, 2)];
    let rows = place_cells(&segments, 2, 5, 3).expect("inside NCBI's array");
    assert_eq!(row_text(&rows), [". . . . .", ". 0 1 2 ."]);
    // NCBI writes in the hits' order: a later hit's cell at a flat index that an earlier hit
    // wrote holds the later hit's column (record position 1, flat 6, over hit 0's spilled cell).
    let segments = [(0, 6, 2, 0), (1, 1, 1, 2)];
    let rows = place_cells(&segments, 2, 5, 3).expect("inside NCBI's array");
    assert_eq!(row_text(&rows), [". . . . .", ". 2 1 . ."]);
}

#[test]
fn a_cell_at_the_arrays_end_places_every_hit_at_the_ranges_coordinates() {
    // L = 5, N = 2: hit 1's segment at record positions 5..7 reaches flat index 1 * 5 + 5 = 10 =
    // N * L, past NCBI's array (the range 2-364 of q_serpin: its kept hit ends at the record's
    // last letter). Every hit is then placed at the range's coordinates (record position minus
    // `range_start` 1), the alignment of the range's letters as a query of their own.
    let segments = [(0, 1, 2, 0), (1, 5, 1, 2)];
    let rows = place_cells(&segments, 2, 5, 1).expect("inside the range");
    assert_eq!(row_text(&rows), ["0 1 . . .", ". . . . 2"]);
    // A cell just before the array's end keeps NCBI's flat placement.
    let segments = [(0, 1, 2, 0), (1, 4, 1, 2)];
    let rows = place_cells(&segments, 2, 5, 1).expect("inside NCBI's array");
    assert_eq!(row_text(&rows), [". 0 1 . .", ". . . . 2"]);
}

#[test]
fn a_segment_of_negative_length_is_still_refused() {
    let base =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/delta/deltatest");
    let files = CdProfileFiles::open(&base, "DELTABLAST").expect("the fixture's profile files");
    let mut negative = hit(1, 1e-5, &[(0, 0, 3)]);
    negative.segments[0].query.to = negative.segments[0].query.from - 1;
    let error = build_msa(&[1; 10], &[negative], &files, 0).expect_err("refused");
    assert!(
        error
            .to_string()
            .contains("has a segment of negative length or start"),
        "{error}"
    );
    // A segment past the searched letters is placed at NCBI's index (here in the next hit's row),
    // not refused.
    let hits = [hit(1, 1e-5, &[(8, 0, 3)]), hit(2, 1e-5, &[(5, 0, 1)])];
    let msa = build_msa(&[1; 10], &hits, &files, 0).expect("placed");
    assert_eq!(msa.rows.len(), 2);
    assert_eq!(msa.rows[1][0], CdCell::aligned(2));
    assert_eq!(msa.rows[1][1], CdCell::UNALIGNED);
}

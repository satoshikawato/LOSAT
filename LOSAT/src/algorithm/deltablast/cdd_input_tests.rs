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
fn equal_evalues_of_one_domain_with_overlapping_query_ranges_are_refused() {
    let hits = vec![hit(1, 1e-5, &[(0, 0, 20)]), hit(1, 1e-5, &[(10, 50, 20)])];
    assert!(matches!(
        remove_multiple_cd_hits(hits),
        Err(CddRefusal::Tie(_))
    ));
    // Disjoint query ranges are no tie.
    let hits = vec![hit(1, 1e-5, &[(0, 0, 10)]), hit(1, 1e-5, &[(50, 50, 10)])];
    assert!(remove_multiple_cd_hits(hits).is_ok());
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

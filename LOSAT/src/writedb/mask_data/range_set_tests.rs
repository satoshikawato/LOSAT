//! Tests of `CSeq_loc::Add` as ported (`range_set`): libstdc++'s sort, empty and whole ranges,
//! the handles that the address orders, `giim` handles. The answers are NCBI's: the review probes
//! of M4a-1 (`p1`, `p1b`, `p4`, makeblastdb 2.17.0), `tests/seqdb/mask_ids/cases.tsv` and the
//! `seqloc_*` files of `tests/seqdb/mask_parse` (run end to end in `tests_e2e.rs`).

use super::*;
use crate::seqid::model::{ObjectId, TextseqId};
use crate::writedb::mask_info::MaskInterval;

fn lid(n: i32) -> ParsedSeqId {
    ParsedSeqId::Local(ObjectId::Id(n))
}

fn iv(from: u32, to: u32, strand: Option<u8>) -> MaskPart {
    MaskPart::Interval(MaskInterval { from, to, strand })
}

fn plain(ranges: &[(u32, u32)]) -> Vec<MaskPart> {
    ranges.iter().map(|&(a, b)| iv(a, b, None)).collect()
}

/// The offsets of the sequence with `ids` for algorithm 30 after `masks` are inserted in order,
/// and the Seq-ids named when NCBI's order of them is the addresses' (`address_ordered`).
fn looked_up(
    masks: &[(ParsedSeqId, Vec<MaskPart>)],
    ids: &[ParsedSeqId],
) -> (Vec<(u32, u32)>, Option<String>) {
    let mut set = MaskedRangeSet::new();
    for (id, parts) in masks {
        set.insert(30, id.clone(), parts.clone()).unwrap();
    }
    let found = set.get_ranges(ids);
    let offsets = found
        .ranges
        .first()
        .map(|r| r.offsets.clone())
        .unwrap_or_default();
    (offsets, found.address_ordered)
}

/// `looked_up`'s offsets where NCBI's order is decided, else the Seq-ids it names.
fn offsets(
    masks: &[(ParsedSeqId, Vec<MaskPart>)],
    ids: &[ParsedSeqId],
) -> Result<Vec<(u32, u32)>, String> {
    match looked_up(masks, ids) {
        (offsets, None) => Ok(offsets),
        (_, Some(seq_ids)) => Err(seq_ids),
    }
}

fn parse_ranges(text: &str) -> Vec<Range> {
    text.split_whitespace()
        .map(|item| {
            let (from, to_open) = item.split_once(':').unwrap();
            Range {
                from: from.parse().unwrap(),
                to_open: to_open.parse().unwrap(),
            }
        })
        .collect()
}

/// `stl_sort` and its heap fallback leave the ranges in the order of libstdc++ (`std::sort`,
/// `std::partial_sort` over all of them), equivalent empty ranges included
/// (`tests/seqdb/stl_sort/README.md`). The last two cases reach introsort's depth limit, so
/// `stl_sort` heap-sorts 97 and 123 of their ranges, empty ones among them.
#[test]
fn stl_sort_gives_libstdcxx_s_permutations() {
    let text = include_str!("../../../tests/seqdb/stl_sort/cases.txt");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len() % 3, 0);
    let mut cases = 0;
    for triple in lines.chunks(3) {
        let mut head = triple[0].splitn(4, ' ');
        assert_eq!(head.next(), Some("CASE"));
        let reverse = head.next() == Some("1");
        let n: usize = head.next().unwrap().parse().unwrap();
        let input = parse_ranges(head.next().unwrap());
        assert_eq!(input.len(), n);
        let less = |a: &Range, b: &Range| range_less(a, b, reverse);
        let mut sorted = input.clone();
        stl_sort(&mut sorted, &less);
        assert_eq!(
            sorted,
            parse_ranges(triple[1].strip_prefix("SORTED").unwrap())
        );
        let mut heap = input.clone();
        heap_sort(&mut heap, &less);
        assert_eq!(heap, parse_ranges(triple[2].strip_prefix("HEAP").unwrap()));
        cases += 1;
    }
    assert_eq!(cases, 8);
}

/// Empty ranges (`to < from`, `to` = 4294967295) are merged into the next range or written as an
/// `empty` Seq-loc that is skipped; NCBI's answers of the review probes p1 and p1b and of the
/// `seqloc_*` fixtures.
#[test]
fn empty_ranges_combine_as_ncbi_does() {
    // p1: {10-20 lcl 1}, {60-50 lcl 1}, {30-40, 0-5 lcl 2} -> (10,21) and (30,41),(0,6).
    let p1 = [
        (lid(1), plain(&[(10, 20)])),
        (lid(1), plain(&[(60, 50)])),
        (lid(2), plain(&[(30, 40), (0, 5)])),
    ];
    assert_eq!(offsets(&p1, &[lid(1)]).unwrap(), [(10, 21)]);
    assert_eq!(offsets(&p1, &[lid(2)]).unwrap(), [(30, 41), (0, 6)]);
    // p1b: {51-70} + {60-50} -> (51,71).
    let p1b = [(lid(1), plain(&[(51, 70)])), (lid(1), plain(&[(60, 50)]))];
    assert_eq!(offsets(&p1b, &[lid(1)]).unwrap(), [(51, 71)]);
    // One Seq-loc is kept as given, empty ranges included (the writer checks them).
    assert_eq!(
        offsets(&[(lid(1), plain(&[(6, 5), (20, 25)]))], &[lid(1)]).unwrap(),
        [(6, 6), (20, 26)]
    );
    // An open end that wraps to 0 contains every later range from its start on
    // (`seqloc_to_max_absorbs`).
    let absorbs = [
        (lid(1), plain(&[(5, u32::MAX)])),
        (lid(1), plain(&[(10, 20)])),
        (lid(1), plain(&[(2, 3)])),
    ];
    assert_eq!(offsets(&absorbs, &[lid(1)]).unwrap(), [(2, 4)]);
    // The order of two empty ranges that the sort leaves as they came decides it
    // (`seqloc_empties_order_a` / `_b`).
    let order = |first: (u32, u32), second: (u32, u32)| {
        offsets(
            &[
                (lid(1), plain(&[first, second])),
                (lid(1), plain(&[(200, 300)])),
            ],
            &[lid(1)],
        )
        .unwrap()
    };
    assert_eq!(order((100, u32::MAX), (60, 50)), [(200, 301)]);
    assert_eq!(order((60, 50), (100, u32::MAX)), Vec::<(u32, u32)>::new());
    // `empty` Seq-locs are skipped.
    assert_eq!(
        offsets(
            &[(
                lid(1),
                vec![MaskPart::Empty, iv(5, 14, None), MaskPart::Empty]
            )],
            &[lid(1)]
        )
        .unwrap(),
        [(5, 15)]
    );
}

/// A whole range (a `whole` Seq-loc, or an interval from 0 to 4294967294) is written as
/// `0..4294967295`, which the writer's bounds check refuses as NCBI's does.
#[test]
fn whole_ranges_are_written_whole() {
    let whole = || iv(0, u32::MAX - 1, None);
    assert_eq!(
        offsets(&[(lid(1), vec![whole()])], &[lid(1)]).unwrap(),
        [(0, u32::MAX)]
    );
    // Sorted first, it contains the other ranges.
    assert_eq!(
        offsets(
            &[(lid(1), vec![whole()]), (lid(1), plain(&[(5, 14)]))],
            &[lid(1)]
        )
        .unwrap(),
        [(0, u32::MAX)]
    );
}

/// Two handles with the same `m_Packed` are refused only where their order changes what is
/// written: the review probe p4 (`gi 5` plus, `gb|AY000005.1|` minus) is NCBI's (10,21),(30,41).
#[test]
fn handles_that_the_address_orders_are_refused_only_on_one_strand() {
    let gi = ParsedSeqId::Gi(5);
    let gb = ParsedSeqId::Genbank(TextseqId {
        accession: Some(b"AY000005".to_vec()),
        version: Some(1),
        ..TextseqId::default()
    });
    let ids = [gi.clone(), gb.clone()];
    let p4 = [
        (gi.clone(), vec![iv(10, 20, Some(1))]),
        (gb.clone(), vec![iv(30, 40, Some(2))]),
    ];
    assert_eq!(offsets(&p4, &ids).unwrap(), [(10, 21), (30, 41)]);
    // Both on the plus strand (cases.tsv n59/n60: the order follows the addresses): noted, and
    // written in LOSAT's order (the refusal comes at the close of the build).
    let tie = [
        (gi.clone(), plain(&[(10, 40)])),
        (gb.clone(), plain(&[(50, 60)])),
    ];
    assert_eq!(
        looked_up(&tie, &ids),
        (
            vec![(10, 41), (50, 61)],
            Some("gi|5 and gb|AY000005.1|".to_string())
        )
    );
    // The same ranges from both: either order writes the same.
    let same = [
        (gi.clone(), plain(&[(10, 40)])),
        (gb.clone(), plain(&[(10, 40)])),
    ];
    assert_eq!(offsets(&same, &ids).unwrap(), [(10, 41), (10, 41)]);
    // One of them has only an empty Seq-loc there: nothing of it is written.
    let hidden = [
        (gi.clone(), plain(&[(10, 40)])),
        (gb, vec![MaskPart::Empty]),
    ];
    assert_eq!(offsets(&hidden, &ids).unwrap(), [(10, 41)]);
}

fn giim(db: Option<&str>) -> ParsedSeqId {
    ParsedSeqId::Giim(GiimportId {
        id: 7,
        db: db.map(|d| d.as_bytes().to_vec()),
        release: None,
    })
}

/// `giim` Seq-ids find the first mask handle of their number that is not `Equals` to them
/// (cases.tsv n36-n45, the sequence `gim|7`).
#[test]
fn giim_handles_follow_the_order_of_the_masks() {
    let seq = [giim(None)];
    let one = |id: ParsedSeqId, ranges: &[(u32, u32)]| (id, plain(ranges));
    // n36: an equal giim does not find it.
    assert_eq!(offsets(&[one(giim(None), &[(10, 40)])], &seq).unwrap(), []);
    // n37: one with a db does.
    assert_eq!(
        offsets(&[one(giim(Some("x")), &[(10, 40)])], &seq).unwrap(),
        [(10, 41)]
    );
    // n40: the second mask lands on the first handle (whose ranges NCBI orders by address); the
    // sequence does not find it.
    let n40 = [
        one(giim(None), &[(10, 40)]),
        one(giim(Some("x")), &[(50, 60)]),
    ];
    assert_eq!(offsets(&n40, &seq).unwrap(), []);
    // n41: it finds the combination whose order is the addresses'.
    let n41 = [
        one(giim(Some("x")), &[(10, 40)]),
        one(giim(None), &[(50, 60)]),
    ];
    assert_eq!(
        looked_up(&n41, &seq),
        (
            vec![(10, 41), (50, 61)],
            Some("gim|7 and gim|7".to_string())
        )
    );
    // n65: the same range in both: either address order writes NCBI's 10-41;10-41;.
    let n65 = [
        one(giim(Some("x")), &[(10, 40)]),
        one(giim(None), &[(10, 40)]),
    ];
    assert_eq!(offsets(&n65, &seq).unwrap(), [(10, 41), (10, 41)]);
    // n66: a third mask combined with it: its range (first handle) and the two others are
    // ordered by address (NCBI wrote 70-81;10-41;10-41;).
    let n66 = [
        one(giim(Some("x")), &[(10, 40)]),
        one(giim(None), &[(10, 40)]),
        one(giim(Some("y")), &[(70, 80)]),
    ];
    assert!(offsets(&n66, &seq).is_err());
    // n44: two algorithms: each finds its own mask under the first handle.
    let mut set = MaskedRangeSet::new();
    set.insert(11, giim(Some("x")), plain(&[(10, 40)])).unwrap();
    set.insert(12, giim(None), plain(&[(50, 60)])).unwrap();
    let found = set.get_ranges(&seq);
    assert_eq!(found.address_ordered, None);
    assert_eq!(found.ranges[0].offsets, [(10, 41)]);
    assert_eq!(found.ranges[1].offsets, [(50, 61)]);
}

/// A negative algorithm id is LOSAT's refusal (NCBI's vector cannot take it).
#[test]
fn a_negative_algorithm_id_is_refused() {
    let mut set = MaskedRangeSet::new();
    let error = set.insert(-1, lid(1), plain(&[(1, 2)])).unwrap_err();
    assert!(matches!(error, WriteDbError::Unsupported(_)), "{error:?}");
}

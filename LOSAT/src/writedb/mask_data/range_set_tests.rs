//! Tests of `CSeq_loc::Add` as ported (`range_set`): the sort as observed, empty and whole ranges,
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

/// One line of `tests/unit/sort_observed/sort_calls.txt` (fields in its README): `sort` (`s`) or
/// `partial_sort` over the whole range (`p`), the length, the answers, the calls (the positions
/// handed as first and second argument), the order after the sort, and the ranges of a line made
/// with NCBI's comparators (`G`: whether the comparator was the reverse one, and the ranges).
struct ObservedSort {
    whole_sort: bool,
    len: usize,
    answers: Vec<bool>,
    calls: Vec<(usize, usize)>,
    order: Vec<usize>,
    ranges: Option<(bool, Vec<Range>)>,
}

/// The comma-separated items of a field, none for `-`.
fn items(field: &str) -> Vec<&str> {
    if field == "-" {
        Vec::new()
    } else {
        field.split(',').collect()
    }
}

fn observed_sort(line: &str) -> ObservedSort {
    let fields: Vec<&str> = line.split(' ').collect();
    let answers = if fields[3] == "-" {
        Vec::new()
    } else {
        fields[3].bytes().map(|answer| answer == b'1').collect()
    };
    let calls = items(fields[4])
        .iter()
        .map(|call| {
            let (first, second) = call.split_once(':').unwrap();
            (first.parse().unwrap(), second.parse().unwrap())
        })
        .collect();
    let order = items(fields[5])
        .iter()
        .map(|k| k.parse().unwrap())
        .collect();
    assert_eq!(fields[6], "0", "a read outside the range: {line}");
    let ranges = (fields[0] == "G").then(|| {
        let (reverse, ranges) = fields[7].split_once('/').unwrap();
        let ranges = items(ranges).join(" ");
        (reverse == "1", parse_ranges(&ranges))
    });
    ObservedSort {
        whole_sort: fields[1] == "s",
        len: fields[2].parse().unwrap(),
        answers,
        calls,
        order,
        ranges,
    }
}

// For each recorded sort (every answer sequence of a strict weak order for 0 to 5 elements, random
// strict weak orders with ties, NCBI's comparators on NCBI-shaped ranges, comparators that decide
// lazily and so reach the heap, and `partial_sort` over whole ranges), `sort_as_observed` (or
// `heap_order` for `partial_sort`) hands the comparator the positions that g++ 13 handed it, in the
// same argument order and sequence, and leaves the same order. The lines made with NCBI's
// comparators also give the same answers when the comparator is `range_less` itself.
#[test]
fn the_sort_makes_the_observed_calls() {
    let records = include_str!("../../../tests/unit/sort_observed/sort_calls.txt");
    let (mut sorts, mut all_calls, mut by_ranges) = (0, 0, 0);
    for line in records.lines() {
        let observed = observed_sort(line);
        let mut order: Vec<usize> = (0..observed.len).collect();
        let mut asked = Vec::new();
        let mut less = |&first: &usize, &second: &usize| {
            asked.push((first, second));
            observed
                .answers
                .get(asked.len() - 1)
                .copied()
                .unwrap_or(false)
        };
        if observed.whole_sort {
            sort_as_observed(&mut order, &mut less);
        } else {
            heap_order(&mut order, &mut less);
        }
        assert_eq!(asked, observed.calls, "{line}");
        assert_eq!(order, observed.order, "{line}");
        if let Some((reverse, ranges)) = &observed.ranges {
            let mut sorted = ranges.clone();
            let mut answers = Vec::new();
            sort_as_observed(&mut sorted, &mut |a: &Range, b: &Range| {
                answers.push(range_less(a, b, *reverse));
                answers[answers.len() - 1]
            });
            assert_eq!(answers, observed.answers, "{line}");
            let expected: Vec<Range> = observed.order.iter().map(|&k| ranges[k]).collect();
            assert_eq!(sorted, expected, "{line}");
            by_ranges += 1;
        }
        sorts += 1;
        all_calls += observed.calls.len();
    }
    assert_eq!((sorts, all_calls, by_ranges), (391, 24798, 20));
}

// `tests/seqdb/stl_sort/cases.txt`: the order in which g++ 13's `sort` (`SORTED`) and
// `partial_sort` over all of them (`HEAP`) left eight lists of ranges under NCBI's comparators,
// with many different empty ranges; the last two lists make the sort order a stretch through the
// heap.
#[test]
fn the_ranges_take_the_observed_orders() {
    let text = include_str!("../../../tests/seqdb/stl_sort/cases.txt");
    let lines: Vec<&str> = text.lines().collect();
    let mut cases = 0;
    for (k, line) in lines.iter().enumerate() {
        let Some(case) = line.strip_prefix("CASE ") else {
            continue;
        };
        let mut words = case.splitn(3, ' ');
        let reverse = words.next() == Some("1");
        let len: usize = words.next().unwrap().parse().unwrap();
        let input = parse_ranges(words.next().unwrap());
        assert_eq!(input.len(), len);
        let sorted = parse_ranges(lines[k + 1].strip_prefix("SORTED").unwrap());
        let heap = parse_ranges(lines[k + 2].strip_prefix("HEAP").unwrap());
        let mut less = |a: &Range, b: &Range| range_less(a, b, reverse);
        let mut by_sort = input.clone();
        sort_as_observed(&mut by_sort, &mut less);
        assert_eq!(by_sort, sorted, "{line}");
        let mut by_heap = input.clone();
        heap_order(&mut by_heap, &mut less);
        assert_eq!(by_heap, heap, "{line}");
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

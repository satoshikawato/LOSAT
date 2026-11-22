use super::super::pattern::AMINOACID_TO_NCBISTDAA;
use super::*;

fn encode(text: &str) -> Vec<u8> {
    text.bytes()
        .map(|byte| AMINOACID_TO_NCBISTDAA[byte as usize])
        .collect()
}

fn pattern(text: &str) -> PhiPattern {
    PhiPattern::compile(text.as_bytes()).unwrap()
}

fn hit(end: i32, start: i32) -> PatternHit {
    PatternHit { end, start }
}

#[test]
fn right_one_bits_and_len() {
    assert_eq!(right_one_bits(0b1000, 0b1010), (3, 1));
    assert_eq!(right_one_bits(0, 0b1010), (0, 3));
    assert_eq!(len_of(0b1000, 0b1000), 4);
}

#[test]
fn one_word_hits_and_cap() {
    let ploop = pattern("[AG]-x(4)-G-K-[ST]");
    let seq = encode("WALLLLGKSW");
    assert_eq!(subject_hits(&ploop, &seq).unwrap(), vec![hit(8, 1)]);
    // `G-x(1,3)-K`: one hit per end, the start of the first alternative in NCBI's order that ends there.
    let variable = pattern("G-x(1,3)-K");
    assert_eq!(
        subject_hits(&variable, &encode("GAAKGAKK")).unwrap(),
        vec![hit(3, 0), hit(6, 4), hit(7, 4)]
    );
    let one = pattern("A");
    let many = encode(&"A".repeat(10_001));
    let hits = subject_hits(&one, &many).unwrap();
    assert_eq!(hits.len(), 10_000);
    assert_eq!(hits[9_999], hit(9_999, 9_999));
}

#[test]
fn multi_word_hits_and_refusal() {
    let long = pattern("A-x(30)");
    assert_eq!(long.kind, PatternKind::MultiWord);
    let hits = subject_hits(&long, &encode(&"A".repeat(40))).unwrap();
    assert_eq!(hits.len(), 10);
    assert_eq!(hits[0], hit(30, 0));
    let all = pattern("x(31)");
    assert_eq!(
        subject_hits(&all, &encode(&"L".repeat(10_030)))
            .unwrap()
            .len(),
        10_000
    );
    assert!(matches!(
        subject_hits(&all, &encode(&"L".repeat(10_031))),
        Err(PhiError::Refusal(_))
    ));
}

#[test]
fn very_long_hits_extend_right_and_left() {
    let very = pattern("[AG]-x(4)-G-K-[ST]-x(0,25)-L");
    assert_eq!(very.kind, PatternKind::VeryLong);
    assert_eq!(
        subject_hits(&very, &encode("ALLLLGKSLL")).unwrap(),
        vec![hit(8, 0), hit(9, 0)]
    );
    assert_eq!(subject_hits(&very, &encode("ALLLLGKSW")).unwrap(), vec![]);
    // A most specific word right of the first: the words before it are found leftwards.
    let left = pattern("L-x(0,40)-W-W-W-W-W-W");
    assert_eq!(left.kind, PatternKind::VeryLong);
    assert_eq!(left.most_specific, 1);
    assert_eq!(
        subject_hits(&left, &encode("LALWWWWWW")).unwrap(),
        vec![hit(8, 0), hit(8, 2)]
    );
}

#[test]
fn query_occurrences_by_segment() {
    let ll = pattern("L-L");
    assert_eq!(
        query_occurrences(&ll, &encode("LLL"), &[(0, 2)]).unwrap(),
        Occurrences::Found(vec![
            Occurrence {
                offset: 0,
                length: 2
            },
            Occurrence {
                offset: 1,
                length: 2
            }
        ])
    );
    let gk = pattern("G-K");
    assert_eq!(
        query_occurrences(&gk, &encode("GK"), &[(0, 1)]).unwrap(),
        Occurrences::WholeQuery
    );
    // A hit across two segments is not found; each segment's hits are offset by its start.
    assert_eq!(
        query_occurrences(&gk, &encode("GKAGKG"), &[(0, 3), (4, 5)]).unwrap(),
        Occurrences::Found(vec![Occurrence {
            offset: 0,
            length: 2
        }])
    );
    assert_eq!(
        query_occurrences(&gk, &encode("AGKAGK"), &[(0, 5)]).unwrap(),
        Occurrences::Found(vec![
            Occurrence {
                offset: 1,
                length: 2
            },
            Occurrence {
                offset: 4,
                length: 2
            }
        ])
    );
}

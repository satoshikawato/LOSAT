use super::*;

#[test]
fn an_empty_edit_script_is_one_segment_like_no_script() {
    let mut segments = Vec::new();
    for ops in [None, Some(&[][..])] {
        fill_segments(&mut segments, ops, 7);
        assert_eq!(segments, [(0, 6)]);
    }
    fill_segments(&mut segments, Some(&[][..]), 0);
    assert!(segments.is_empty());
    let ops = [GapEditOp::Sub(3), GapEditOp::Ins(2), GapEditOp::Sub(4)];
    fill_segments(&mut segments, Some(&ops), 9);
    assert_eq!(segments, [(0, 2), (3, 4), (5, 8)]);
}

/// The XML and JSON rows show a row's end gaps as spaces (`end_gaps_as_spaces`, NCBI's
/// `GetWholeAlnSeqString`): the gap columns of the row before its first letter and after its
/// last; the inner gaps and the other row stay.
#[test]
fn a_rows_end_gaps_are_written_as_spaces() {
    use GapEditOp::{Del, Ins, Sub};
    let ops = [Del(2), Sub(3), Ins(1), Sub(2), Del(1)];
    let mut query = b"--ABCQDE-".to_vec();
    end_gaps_as_spaces(&mut query, &ops, AlignRow::Query);
    assert_eq!(query, b"  ABCQDE ");
    let mut subject = b"XYABC-DEZ".to_vec();
    end_gaps_as_spaces(&mut subject, &ops, AlignRow::Subject);
    assert_eq!(subject, b"XYABC-DEZ");
    // A gap start in the subject's row.
    let mut subject = b"---ABC".to_vec();
    end_gaps_as_spaces(&mut subject, &[Ins(3), Sub(3)], AlignRow::Subject);
    assert_eq!(subject, b"   ABC");
    // An X-drop script keeps every `-`.
    let mut query = b"AB-CD".to_vec();
    end_gaps_as_spaces(&mut query, &[Sub(2), Del(1), Sub(2)], AlignRow::Query);
    assert_eq!(query, b"AB-CD");
}

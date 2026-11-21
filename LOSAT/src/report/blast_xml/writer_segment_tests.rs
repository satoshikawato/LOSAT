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

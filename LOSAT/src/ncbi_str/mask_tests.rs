use super::matches_mask;

#[test]
fn question_marks_and_stars() {
    assert!(matches_mask(b"p5.pin", b"p5.???"));
    assert!(!matches_mask(b"p5.ab", b"p5.???"));
    assert!(!matches_mask(b"p5.pinx", b"p5.???"));
    assert!(!matches_mask(b"pmv.00.pin", b"pmv.???"));
    assert!(matches_mask(b".ptm", b"*tm"));
    assert!(matches_mask(b"tm", b"*tm"));
    assert!(!matches_mask(b".tmx", b"*tm"));
    assert!(matches_mask(b"abc", b"**"));
    assert!(matches_mask(b"", b"*"));
    assert!(!matches_mask(b"", b"?"));
}

#[test]
fn sets_ranges_and_escapes() {
    assert!(matches_mask(b"a1", b"a[0-9]"));
    assert!(!matches_mask(b"ax", b"a[0-9]"));
    assert!(matches_mask(b"ax", b"a[!0-9]"));
    assert!(matches_mask(b"a-", b"a[x-]"));
    assert!(!matches_mask(b"a1", b"a[0-9"));
    assert!(matches_mask(b"a*", b"a\\*"));
    assert!(!matches_mask(b"ab", b"a\\*"));
    // A NUL byte ends the text as its end does.
    assert!(matches_mask(b"ab\0cd", b"ab"));
}

use super::*;

#[test]
fn has_mods_is_the_first_bracket_pair_with_an_equal_sign() {
    assert!(has_mods(b"a [organism=x] b"));
    assert!(has_mods(b"[a=b]"));
    assert!(has_mods(b"[[a=b]]"));
    assert!(!has_mods(b"a [organism] b [x=y]"));
    assert!(!has_mods(b"a ] [x=y]"));
    assert!(!has_mods(b"a [x=y"));
    assert!(!has_mods(b""));
}

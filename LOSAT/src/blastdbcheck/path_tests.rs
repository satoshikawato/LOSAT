use super::*;

#[test]
fn a_path_splits_into_its_directory_base_and_extension() {
    let file = FileParts::of(b"/run/cases/a/gi.list.gil");
    assert_eq!(file.dir, b"/run/cases/a/");
    assert_eq!(file.base, b"gi.list");
    assert_eq!(file.ext, b".gil");
    let bare = FileParts::of(b"p5");
    assert_eq!(bare.dir, [b'.', RULES.separator()]);
    assert_eq!(bare.base, b"p5");
    assert_eq!(bare.ext, b"");
    let dotted = FileParts::of(b"dir/.gil/");
    assert_eq!(dotted.path, b"dir/.gil");
    assert_eq!(dotted.base, b"");
    assert_eq!(dotted.ext, b".gil");
}

#[test]
fn make_path_joins_an_absolute_name_as_any_other() {
    assert_eq!(
        make_path(b"/run/a/", b"/run/a/gi.txt"),
        b"/run/a//run/a/gi.txt"
    );
    assert_eq!(make_path(b"", b"x"), b"x");
    assert_eq!(
        make_path(b"d", b"x"),
        [&b"d"[..], &[RULES.separator()], b"x"].concat()
    );
}

#[test]
fn the_entries_of_a_directory_are_named_after_its_reset_path() {
    assert_eq!(entries_base(b"/run/in//"), b"/run/in/");
    assert_eq!(entries_base(b"/"), b"/");
    assert_eq!(entries_base(b""), [b'.', RULES.separator()]);
}

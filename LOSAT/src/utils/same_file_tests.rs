use super::*;

#[test]
fn the_spellings_of_one_path_share_a_form() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-same-file-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(dir.join("sub")).unwrap();
    std::fs::write(dir.join("a"), b"a").unwrap();
    std::fs::write(dir.join("b"), b"b").unwrap();
    let a = path_identity(&dir.join("a")).unwrap();
    assert_eq!(path_identity(&dir.join("sub/../a")), Some(a.clone()));
    assert_eq!(path_identity(&dir.join("./sub/./../a")), Some(a.clone()));
    assert_ne!(path_identity(&dir.join("b")), Some(a));
    // A file that is not there has no form (NCBI creates `-out` later).
    assert_eq!(path_identity(&dir.join("sub/../none")), None);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_lexical_form_resolves_dots_by_name() {
    let root = Path::new(if cfg!(windows) { "C:\\" } else { "/" });
    let at = |tail: &str| lexical_absolute(&root.join(tail)).unwrap();
    assert_eq!(at("a/./b/../c"), root.join("a").join("c"));
    // `..` at the root stays there.
    assert_eq!(at("../../a"), root.join("a"));
    assert_eq!(at("a/b/../../.."), root.to_path_buf());
    // A relative name is taken from the current directory.
    let cwd = std::env::current_dir().unwrap();
    assert_eq!(
        lexical_absolute(Path::new("x/./y/../z")),
        lexical_absolute(&cwd.join("x/z"))
    );
}

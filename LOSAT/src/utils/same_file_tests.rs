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

#[test]
fn same_file_target_matches_aliases_and_missing_targets() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-same-file-target-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(dir.join("sub")).unwrap();
    let a = dir.join("a.fa");
    let b = dir.join("b.fa");
    std::fs::write(&a, b">a\nACGT\n").unwrap();
    std::fs::write(&b, b">b\nTGCA\n").unwrap();

    assert!(same_file_target(&a, &a));
    assert!(same_file_target(&dir.join("sub/../a.fa"), &a));
    assert!(!same_file_target(&a, &b));
    // Directories are never matched (opening a directory as -logfile falls back to stderr).
    assert!(!same_file_target(&dir.join("sub"), &dir.join("sub")));
    assert!(!same_file_target(&dir.join("sub"), &a));
    // Existing vs non-existent are distinct.
    let missing = dir.join("missing.fa");
    assert!(!same_file_target(&a, &missing));
    assert!(!same_file_target(&missing, &a));
    // Two spellings of the same missing file match.
    assert!(same_file_target(&missing, &dir.join("sub/../missing.fa")));
    assert!(!same_file_target(&missing, &dir.join("other_missing.fa")));

    #[cfg(unix)]
    {
        use std::os::unix::fs::symlink;
        let sym_a = dir.join("link_a.fa");
        symlink(&a, &sym_a).unwrap();
        assert!(same_file_target(&sym_a, &a));
        assert!(same_file_target(&a, &sym_a));

        let hard_a = dir.join("hard_a.fa");
        std::fs::hard_link(&a, &hard_a).unwrap();
        assert!(same_file_target(&hard_a, &a));

        let dir_link = dir.join("dir_link");
        symlink(dir.join("sub"), &dir_link).unwrap();
        assert!(same_file_target(
            &dir_link.join("../a.fa"),
            &dir.join("a.fa")
        ));
        assert!(same_file_target(
            &dir_link.join("../missing.fa"),
            &dir.join("missing.fa")
        ));

        let dangling = dir.join("dangling.fa");
        symlink(&missing, &dangling).unwrap();
        assert!(same_file_target(&dangling, &missing));
        assert!(same_file_target(&missing, &dangling));

        // Character devices such as /dev/null are never matched.
        assert!(!same_file_target(
            Path::new("/dev/null"),
            Path::new("/dev/null")
        ));
    }

    let _ = std::fs::remove_dir_all(&dir);
}

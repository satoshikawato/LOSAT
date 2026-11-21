use super::*;

/// `CDirEntry::SplitPath` ends the directory part at `:`, `/` or `\` on Windows and at `/` alone
/// elsewhere (ncbifile.cpp:105-116): checked here for both on every platform, since the Windows
/// job of the CI would be the only one to see a wrong `cfg!(windows)`.
#[test]
fn the_separators_of_windows_end_the_directory_part_there_only() {
    let base = |path: &str, xml: bool, with_dir: bool, windows: bool| {
        String::from_utf8(base_name_of(path.as_bytes(), xml, with_dir, windows)).unwrap()
    };
    // (path, xml, with_dir, base on Windows, base elsewhere)
    let cases = [
        ("d\\X.XML", true, true, "d\\X", "d\\X"),
        ("d\\X.XML", true, false, "X", "d\\X"),
        ("C:o.xml", true, true, "C:o", "C:o"),
        ("C:o.xml", true, false, "o", "C:o"),
        ("a/b\\c.json", false, true, "a/b\\c", "a/b\\c"),
        ("a/b\\c.json", false, false, "c", "b\\c"),
        // the dot before a separator is not the file's: `d.xml\o` has no extension on Windows
        ("d.xml\\o", true, true, "d.xml\\o", "d.xml\\o"),
        ("d.xml\\o", true, false, "o", "d.xml\\o"),
        ("C:\\d\\.json", false, true, "C:\\d\\", "C:\\d\\"),
    ];
    for (path, xml, with_dir, on_windows, elsewhere) in cases {
        assert_eq!(
            base(path, xml, with_dir, true),
            on_windows,
            "{path} windows"
        );
        assert_eq!(base(path, xml, with_dir, false), elsewhere, "{path} unix");
    }
}

/// The names of the side files and the master's entries on the platform of the build.
#[cfg(windows)]
#[test]
fn windows_paths_name_the_side_files_and_the_master_entries() {
    let files = FsSideFiles::new(Path::new("d\\X.XML"), Xml2Kind::Files14);
    assert_eq!(files.path(1), PathBuf::from("d\\X_1.xml"));
    assert_eq!(files.listed(1), b"X_1.xml");
    let files = FsSideFiles::new(Path::new("C:o.json"), Xml2Kind::Files13);
    assert_eq!(files.path(2), PathBuf::from("C:o_2.json"));
    assert_eq!(files.listed(2), b"o.json_2.json");
}

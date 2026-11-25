use super::*;
use crate::utils::test_dir::TempDir;

#[test]
fn fnmatch_wildcards_and_classes() {
    assert!(fnmatch_bytes(b"db.*.pin", b"db.00.pin"));
    assert!(fnmatch_bytes(b"db.*.pin", b"db..pin"));
    assert!(fnmatch_bytes(b"db.*.pin", b"db.01.02.pin"));
    assert!(!fnmatch_bytes(b"db.*.pin", b"db.00.nin"));
    assert!(fnmatch_bytes(b"db.01.p??", b"db.01.pin"));
    assert!(fnmatch_bytes(b"db.01.p??", b"db.01.psq"));
    assert!(!fnmatch_bytes(b"db.01.p??", b"db.01.pa"));
    assert!(!fnmatch_bytes(b"db.01.p??", b"db.01.pinx"));
    assert!(fnmatch_bytes(b"a[0-2].pin", b"a1.pin"));
    assert!(!fnmatch_bytes(b"a[0-2].pin", b"a3.pin"));
    assert!(fnmatch_bytes(b"a[!0-2].pin", b"a3.pin"));
    assert!(fnmatch_bytes(b"a[0.pin", b"a[0.pin"));
}

#[test]
fn surrogateescape_ordering_and_bracket_ranges() {
    // Python surrogateescape places invalid Unix bytes in U+DC80..=U+DCFF:
    // U+00E9 ('é') < U+DC80 (0x80) < U+1F600 ('😀').
    let e_acute = "db.\u{00e9}.pin".as_bytes();
    let raw_80 = b"db.\x80.pin".as_slice();
    let emoji = "db.\u{1f600}.pin".as_bytes();
    assert_eq!(cmp_py_path(e_acute, raw_80), std::cmp::Ordering::Less);
    assert_eq!(cmp_py_path(raw_80, emoji), std::cmp::Ordering::Less);
    assert_eq!(cmp_py_path(e_acute, emoji), std::cmp::Ordering::Less);

    // Bracket range `[a-\x80]` includes U+00E9 ('é') and excludes U+1F600 ('😀').
    assert!(fnmatch_bytes_case(
        b"[a-\x80]",
        "\u{00e9}".as_bytes(),
        false
    ));
    assert!(!fnmatch_bytes_case(
        b"[a-\x80]",
        "\u{1f600}".as_bytes(),
        false
    ));
}

#[test]
fn windows_path_joining_and_unicode_case_folding() {
    // `ntpath.split` + `ntpath.join` joins non-root prefixes with `\`, even when the input used `/`.
    for (pat, want_joined) in [
        (b"in/p5.*.pin".as_slice(), b"in\\p5.09.pin".as_slice()),
        (b"in//p5.*.pin", b"in\\p5.09.pin"),
        (b"in\\p5.*.pin", b"in\\p5.09.pin"),
        (b"a/b/c.*.pin", b"a/b\\p5.09.pin"),
        (b"/p5.*.pin", b"/p5.09.pin"),
        (b"C:/p5.*.pin", b"C:/p5.09.pin"),
        (b"C://p5.*.pin", b"C:/p5.09.pin"),
        (b"C:\\p5.*.pin", b"C:\\p5.09.pin"),
        (b"C:p5.*.pin", b"C:p5.09.pin"),
        (b"\\\\srv\\sh\\p5.*.pin", b"\\\\srv\\sh\\p5.09.pin"),
        (b"//srv/sh/p5.*.pin", b"//srv/sh/p5.09.pin"),
        (
            b"\\\\?\\UNC\\srv\\sh\\p5.*.pin",
            b"\\\\?\\UNC\\srv\\sh\\p5.09.pin",
        ),
        (b"p5.*.pin", b"p5.09.pin"),
    ] {
        let (dir, _base) = path_split_os(pat, true);
        let got = path_join_os(dir, b"p5.09.pin", true);
        assert_eq!(got, want_joined, "pat={}", String::from_utf8_lossy(pat));
    }

    // Windows `ntpath.normcase` uses 1-to-1 invariant lowercasing (`LCMapStringEx`):
    // U+0130 ('İ') maps to single 'i' (not 'i' + U+0307) and U+03A3 ('Σ') maps to 'σ' (not final 'ς').
    assert!(fnmatch_bytes_case(
        "\u{00c4}.pin".as_bytes(),
        "\u{00e4}.PIN".as_bytes(),
        true
    ));
    assert!(fnmatch_bytes_case(
        "\u{0130}.pin".as_bytes(),
        "i.pin".as_bytes(),
        true
    ));
    assert!(!fnmatch_bytes_case(
        "\u{0130}.pin".as_bytes(),
        "i\u{0307}.pin".as_bytes(),
        true
    ));
    assert!(fnmatch_bytes_case(
        "\u{0391}\u{039f}\u{03a3}.pin".as_bytes(),
        "\u{03b1}\u{03bf}\u{03c3}.pin".as_bytes(),
        true
    ));
    assert!(!fnmatch_bytes_case(
        "\u{0391}\u{039f}\u{03a3}.01.pin".as_bytes(),
        "\u{03b1}\u{03bf}\u{03c2}.01.pin".as_bytes(),
        true
    ));
    assert!(!fnmatch_bytes_case(
        "\u{00c4}.pin".as_bytes(),
        "\u{00e4}.pin".as_bytes(),
        false
    ));
}

#[test]
fn glob_and_glob_sorted_handle_subdirs_hidden_and_metachars() {
    let dir = TempDir::unique("cleanup-glob");
    let sub = dir.join("sub");
    std::fs::create_dir_all(&sub).unwrap();
    for name in [
        "db.00.pin",
        "db.01.pin",
        "db.01.psq",
        "db.01.phr",
        "db.01.pinx",
        "..01.pin",
    ] {
        std::fs::write(sub.join(name), b"").unwrap();
    }
    let sep = if cfg!(windows) { '\\' } else { '/' };
    let prefix = sub.to_string_lossy();
    let entry = |name: &str| format!("{prefix}{sep}{name}").into_bytes();
    let vols = glob_sorted(format!("{prefix}/db.*.pin").as_bytes());
    assert_eq!(vols, vec![entry("db.00.pin"), entry("db.01.pin")]);
    // Unsorted `glob` preserves `std::fs::read_dir` order.
    let expected_readdir: Vec<Vec<u8>> = std::fs::read_dir(&sub)
        .unwrap()
        .flatten()
        .map(|e| os_str_bytes(&e.file_name()).into_owned())
        .filter(|name| fnmatch_bytes(b"db.01.p??", name))
        .map(|name| {
            let mut full = prefix.as_bytes().to_vec();
            full.push(sep as u8);
            full.extend_from_slice(&name);
            full
        })
        .collect();
    let files = glob(format!("{prefix}/db.01.p??").as_bytes());
    assert_eq!(files, expected_readdir);
    // Pattern starting with `.` matches dot-prefixed entries.
    let hidden = glob_sorted(format!("{prefix}/.*.pin").as_bytes());
    assert_eq!(hidden, vec![entry("..01.pin")]);
    // Redundant slashes before the basename collapse like Python `os.path.split` + `os.path.join`.
    let collapsed = glob_sorted(format!("{prefix}//db.*.pin").as_bytes());
    assert_eq!(collapsed, vec![entry("db.00.pin"), entry("db.01.pin")]);
}

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
fn glob_sorted_handles_subdirs_hidden_and_metachars() {
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
    let prefix = sub.to_string_lossy();
    let entry = |name: &str| format!("{prefix}/{name}").into_bytes();
    let vols = glob_sorted(format!("{prefix}/db.*.pin").as_bytes());
    assert_eq!(vols, vec![entry("db.00.pin"), entry("db.01.pin")]);
    let files = glob_sorted(format!("{prefix}/db.01.p??").as_bytes());
    assert_eq!(
        files,
        vec![entry("db.01.phr"), entry("db.01.pin"), entry("db.01.psq")]
    );
    // Pattern starting with `.` matches dot-prefixed entries.
    let hidden = glob_sorted(format!("{prefix}/.*.pin").as_bytes());
    assert_eq!(hidden, vec![entry("..01.pin")]);
    // Redundant slashes before the basename collapse like Python `os.path.split` + `os.path.join`.
    let collapsed = glob_sorted(format!("{prefix}//db.*.pin").as_bytes());
    assert_eq!(collapsed, vec![entry("db.00.pin"), entry("db.01.pin")]);
}

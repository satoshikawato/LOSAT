use super::{remove, rename_overwrite};
use std::path::{Path, PathBuf};

/// A directory of its own under the system's temporary directory.
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir = crate::utils::test_dir::temp_dir().join(format!(
            "losat-dir-entry-{name}-{}-{:?}",
            crate::utils::test_dir::process_id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        #[cfg(unix)]
        make_writable(&self.0);
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[cfg(unix)]
fn make_writable(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    if std::fs::symlink_metadata(path).is_ok_and(|meta| meta.is_dir()) {
        let _ = std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755));
        if let Ok(entries) = std::fs::read_dir(path) {
            for entry in entries.map_while(Result::ok) {
                make_writable(&entry.path());
            }
        }
    }
}

/// A directory with a nested directory and files.
fn full_dir(path: &Path) {
    std::fs::create_dir_all(path.join("inner").join("deeper")).unwrap();
    std::fs::write(path.join("inner").join("deeper").join("f"), b"x").unwrap();
    std::fs::write(path.join("a"), b"a").unwrap();
}

#[test]
fn removes_a_file_and_a_directory_with_its_contents() {
    let dir = TempDir::new("full");
    let file = dir.0.join("db.pin");
    std::fs::write(&file, b"x").unwrap();
    assert!(remove(&file));
    assert!(!file.exists());
    let sub = dir.0.join("db.pdb");
    full_dir(&sub);
    assert!(remove(&sub));
    assert!(std::fs::symlink_metadata(&sub).is_err());
    // A missing name is not removed (`eRecursive` has no `fIgnoreMissing`).
    assert!(!remove(&sub));
}

#[cfg(unix)]
#[test]
fn makes_each_directory_writable_and_stops_at_the_first_failure() {
    use std::os::unix::fs::PermissionsExt;
    let dir = TempDir::new("modes");
    // Read-only directories and files are removed: each directory is made writable for
    // the user before its entries are removed.
    let ro = dir.0.join("db.pdb-lock");
    full_dir(&ro);
    std::fs::write(ro.join("ro.txt"), b"r").unwrap();
    std::fs::set_permissions(ro.join("ro.txt"), std::fs::Permissions::from_mode(0o444)).unwrap();
    std::fs::set_permissions(ro.join("inner"), std::fs::Permissions::from_mode(0o555)).unwrap();
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
    assert!(remove(&ro));
    assert!(std::fs::symlink_metadata(&ro).is_err());

    // A directory that cannot be read stops the removal: it gets the user's write
    // permission only, the entries listed before it are gone and the others stay.
    let stuck = dir.0.join("db.ptf");
    std::fs::create_dir_all(stuck.join("locked")).unwrap();
    std::fs::write(stuck.join("locked").join("h"), b"h").unwrap();
    std::fs::set_permissions(stuck.join("locked"), std::fs::Permissions::from_mode(0o000)).unwrap();
    std::fs::set_permissions(&stuck, std::fs::Permissions::from_mode(0o3755)).unwrap();
    // Root reads any directory; the case then says nothing.
    if std::fs::read_dir(stuck.join("locked")).is_ok() {
        eprintln!("skipped (running as root)");
        return;
    }
    assert!(!remove(&stuck));
    let mode = |path: &Path| std::fs::metadata(path).unwrap().permissions().mode() & 0o7777;
    assert_eq!(mode(&stuck.join("locked")), 0o200);
    // The set-group-ID and sticky bits are cleared with the mode change.
    assert_eq!(mode(&stuck), 0o755);
}

#[cfg(unix)]
#[test]
fn removes_a_symbolic_link_and_keeps_its_target() {
    let dir = TempDir::new("link");
    let target = dir.0.join("target");
    full_dir(&target);
    let link = dir.0.join("db.pdb");
    std::os::unix::fs::symlink("target", &link).unwrap();
    assert!(remove(&link));
    assert!(std::fs::symlink_metadata(&link).is_err());
    assert!(target.join("a").exists() && target.join("inner/deeper/f").exists());
    // A link inside a removed directory is removed, not followed.
    let sub = dir.0.join("db.pos");
    std::fs::create_dir(&sub).unwrap();
    std::os::unix::fs::symlink("../target", sub.join("to-target")).unwrap();
    assert!(remove(&sub));
    assert!(!sub.exists() && target.join("a").exists());
}

// `CDirEntry::Rename(fRF_Overwrite)` as makeblastdb 2.17.0 calls it for a volume file
// (probes of 2026-11-11; the same cases run against the built program in
// `tests/makeblastdb_rename_lock.rs`).
#[test]
fn a_file_gets_its_new_name_and_replaces_what_is_there() {
    let dir = TempDir::new("rename");
    let (src, dst) = (dir.0.join("db.00.pin"), dir.0.join("db.pin"));
    std::fs::write(&src, b"new").unwrap();
    assert!(rename_overwrite(&src, &dst));
    assert!(!src.exists());
    assert_eq!(std::fs::read(&dst).unwrap(), b"new");
    // A file at the new name goes first.
    std::fs::write(&src, b"newer").unwrap();
    assert!(rename_overwrite(&src, &dst));
    assert_eq!(std::fs::read(&dst).unwrap(), b"newer");
    // So does a directory, with its contents.
    std::fs::write(&src, b"newest").unwrap();
    std::fs::remove_file(&dst).unwrap();
    full_dir(&dst);
    assert!(rename_overwrite(&src, &dst));
    assert_eq!(std::fs::read(&dst).unwrap(), b"newest");
}

#[test]
fn a_source_that_is_not_there_renames_nothing() {
    let dir = TempDir::new("rename-missing");
    let (src, dst) = (dir.0.join("db.00.pog"), dir.0.join("db.pog"));
    std::fs::write(&dst, b"old").unwrap();
    assert!(!rename_overwrite(&src, &dst));
    // The destination is not touched when the source is missing.
    assert_eq!(std::fs::read(&dst).unwrap(), b"old");
}

#[cfg(unix)]
#[test]
fn a_destination_that_cannot_be_removed_stays_and_the_source_keeps_its_name() {
    use std::os::unix::fs::PermissionsExt;
    let dir = TempDir::new("rename-stuck");
    let (src, dst) = (dir.0.join("db.00.pin"), dir.0.join("db.pin"));
    std::fs::write(&src, b"new").unwrap();
    std::fs::create_dir_all(dst.join("locked")).unwrap();
    std::fs::write(dst.join("locked").join("h"), b"h").unwrap();
    std::fs::set_permissions(dst.join("locked"), std::fs::Permissions::from_mode(0o000)).unwrap();
    // Root reads any directory; the case then says nothing.
    if std::fs::read_dir(dst.join("locked")).is_ok() {
        eprintln!("skipped (running as root)");
        return;
    }
    assert!(!rename_overwrite(&src, &dst));
    assert_eq!(std::fs::read(&src).unwrap(), b"new");
    assert!(dst.is_dir());
}

#[cfg(unix)]
#[test]
fn links_are_entries_of_their_own() {
    let dir = TempDir::new("rename-link");
    // A source that is a link to nothing is renamed as the link (makeblastdb 2.17.0 renamed
    // `db.00.pni -> /nonexistent/x` to `db.pni`).
    let (src, dst) = (dir.0.join("db.00.pog"), dir.0.join("db.pog"));
    std::os::unix::fs::symlink("/nonexistent/target", &src).unwrap();
    assert!(rename_overwrite(&src, &dst));
    assert!(std::fs::symlink_metadata(&src).is_err());
    assert_eq!(
        std::fs::read_link(&dst).unwrap(),
        Path::new("/nonexistent/target")
    );
    // A destination that is a link to nothing exists: it is removed, then the rename goes on.
    std::fs::write(&src, b"s").unwrap();
    assert!(rename_overwrite(&src, &dst));
    assert_eq!(std::fs::read(&dst).unwrap(), b"s");
    // A source that is a link to a device moves as the link.
    let (full, to) = (dir.0.join("db.00.pin"), dir.0.join("db.pin"));
    std::os::unix::fs::symlink("/dev/null", &full).unwrap();
    assert!(rename_overwrite(&full, &to));
    assert_eq!(std::fs::read_link(&to).unwrap(), Path::new("/dev/null"));
    // A source that is a link to a file moves as the link.
    let (link, new) = (dir.0.join("db.00.nog"), dir.0.join("db.nog"));
    std::os::unix::fs::symlink("db.pog", &link).unwrap();
    assert!(rename_overwrite(&link, &new));
    assert!(std::fs::symlink_metadata(&link).is_err());
    assert_eq!(std::fs::read_link(&new).unwrap(), Path::new("db.pog"));
}

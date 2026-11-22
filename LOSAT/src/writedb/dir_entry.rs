//! `CDirEntry::Remove` of NCBI's corelib (`c++/src/corelib/ncbifile.cpp`) with its default flags
//! `eRecursive`: how makeblastdb removes a name of an existing database (`Remove` of a `CFile`
//! in `DeleteBlastDb`) and `CDirEntry::Rename` with `fRF_Overwrite`, with which it gives a volume
//! file its final name (`rename_overwrite`). A directory is removed with its contents.

use std::path::Path;

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:707-741,772
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2595-2603
// `CFile` does not override `Remove`, so `Remove` of a `CFile` is this function with
// `eRecursive`. A failure is not reported: the file API logs only with `FileAPILogging`, whose
// default is false (ncbifile.cpp:176-180, 221-230), and `DeleteBlastDb` only counts the names
// removed.
/// `Remove` with `eRecursive` of a `CDirEntry` of `path`: removes `path`, a directory with its
/// contents; `true` when it is removed.
pub(crate) fn remove(path: &Path) -> bool {
    if is_dir_ignore_links(path) {
        return remove_dir(path);
    }
    remove_entry(path)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2203-2241
// `IsDir` with `eIgnoreLinks`: a symbolic link is not a directory, whatever it points to, so the
// link is removed and its target is kept. On Windows NCBI's `GetType` calls `stat` for both values
// of `follow` (ncbifile.cpp:2208-2209); LOSAT does not follow a link there either (not compared
// with a makeblastdb run on Windows).
fn is_dir_ignore_links(path: &Path) -> bool {
    std::fs::symlink_metadata(path).is_ok_and(|meta| meta.is_dir())
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2606-2630
// `eRecursive` has no `fIgnoreMissing`, so a missing entry is a failure; `DeleteReadOnlyFiles`
// is false by default (ncbifile.cpp:146-148), so a read-only file on Windows is one too.
/// `CDirEntry::RemoveEntry`: `remove()` of an entry that is not a directory.
fn remove_entry(path: &Path) -> bool {
    std::fs::remove_file(path).is_ok()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:4342-4456
// `eRecursive` has no `fProcessAll`: the first entry that cannot be removed stops the removal,
// and the entries before it are gone (makeblastdb 2.17.0 on Linux x64, 2026-11-08: a directory
// holding a file and an unreadable directory loses the file and keeps the other).
/// `CDir::Remove` with `eRecursive` of the directory `path`.
fn remove_dir(path: &Path) -> bool {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:1061-1084 (CDirEntry::SetModeEntry)
    // `SetMode` that adds `fWrite` to the user's mode and leaves the group's and the others' modes
    // unchanged, with the default special bits
    // 0 (`CDir::SetMode` of the top entry only, ncbifile.hpp:2040-2044): the permissions of the
    // user, the group and the others are kept, the user's write permission is added, and the
    // set-user-ID, set-group-ID and sticky bits are cleared. Its result is not checked.
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        if let Ok(meta) = std::fs::metadata(path) {
            let mode = (meta.permissions().mode() & 0o777) | 0o200;
            let _ = std::fs::set_permissions(path, std::fs::Permissions::from_mode(mode));
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3996-4010 (CDir::GetEntriesPtr)
    // Every entry is listed, in the order the directory gives, before any is removed; an
    // error of `readdir` ends the list.
    let Ok(entries) = std::fs::read_dir(path) else {
        return false;
    };
    let names: Vec<_> = entries
        .map_while(Result::ok)
        .map(|entry| entry.file_name())
        .collect();
    for name in names {
        let item = path.join(name);
        let removed = if is_dir_ignore_links(&item) {
            remove_dir(&item)
        } else {
            remove(&item)
        };
        if !removed {
            return false;
        }
    }
    std::fs::remove_dir(path).is_ok()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2456-2586
// (CDirEntry::Rename)
// `CWriteDB_File::RenameSingle` and `RenameFileIndex` call it with `fRF_Overwrite` alone (no
// `fRF_FollowLinks`, `fRF_EqualTypes`, `fRF_Update` or `fRF_Backup`) and ignore the result.
// The `LOG_ERROR` macros print only with `FileAPILogging`, whose default is false
// (ncbifile.cpp:176-180, 193-198).
// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:851,3958-3964
// The type of an entry (`GetType()`, `Exists()`) is that of `lstat` on Unix (`GetType`,
// ncbifile.cpp:2203-2241, cited at `is_dir_ignore_links`): a link is an entry of its own, so a
// source that is a link (to a file, to a device, to nothing) is not a file, skips the hard
// link and is renamed as the link, and a destination that is a link is removed (makeblastdb
// 2.17.0 renamed `db.00.pni -> /nonexistent/x` to `db.pni`; probe of 2026-11-12). On Windows
// `GetType` calls `stat`, and `Exists` opens the entry; LOSAT asks `stat` for both there (not
// compared with a makeblastdb run on Windows).
// A source that is a file reaches the hard-link branch first (Unix); its failure, because the
// destination exists, is what sends an existing destination to the removal below. The two
// names of a volume file are in one directory, so a cross-device `rename()` is not expected;
// its `Copy()` is ported for a file only (the mode is copied, the owner and the times of
// `fCF_PreserveAll` are not), and a source that is a directory fails there.
/// `Rename` of a `CDirEntry` of `src` to `dst` with `CDirEntry::fRF_Overwrite`: `true` when `src`
/// got the name `dst`. An existing destination is removed first (`remove`, a directory with its
/// contents); one that stays makes the rename fail, and `src` keeps its name.
pub(crate) fn rename_overwrite(src: &Path, dst: &Path) -> bool {
    let Some(src_type) = entry_type(src) else {
        return false;
    };
    #[cfg(unix)]
    if src_type.is_file() && std::fs::hard_link(src, dst).is_ok() && remove_entry(src) {
        return true;
    }
    if entry_type(dst).is_some() {
        // `dst.Exists()` asks what the `GetType()` above asked.
        let _ = remove(dst);
    }
    if entry_type(dst).is_some() {
        return false;
    }
    let Err(error) = std::fs::rename(src, dst) else {
        return true;
    };
    #[cfg(unix)]
    let cross_device = error.kind() == std::io::ErrorKind::CrossesDevices;
    #[cfg(not(unix))]
    let cross_device = error.kind() == std::io::ErrorKind::PermissionDenied;
    if !cross_device {
        return false;
    }
    if !src_type.is_file() || std::fs::copy(src, dst).is_err() {
        let _ = remove(dst);
        return false;
    }
    remove(src)
}

/// `GetType` with `eIgnoreLinks` of a `CDirEntry` of `path`, when it is not `eUnknown`: `lstat` on
/// Unix, `stat` elsewhere.
fn entry_type(path: &Path) -> Option<std::fs::Metadata> {
    #[cfg(unix)]
    let meta = std::fs::symlink_metadata(path);
    #[cfg(not(unix))]
    let meta = std::fs::metadata(path);
    meta.ok()
}

#[cfg(test)]
#[path = "dir_entry_tests.rs"]
mod tests;

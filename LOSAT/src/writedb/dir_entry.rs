//! `CDirEntry::Remove` of NCBI's corelib (`c++/src/corelib/ncbifile.cpp`) with its default flags
//! `eRecursive`: how makeblastdb removes a name of an existing database (`CFile(name).Remove()`
//! in `DeleteBlastDb`) and `CDirEntry::Rename(..., fRF_Overwrite)`, with which it gives a volume
//! file its final name (`rename_overwrite`). A directory is removed with its contents.

use std::path::Path;

// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:707-741,772
// ```c++
//         fEntry         = (1 << 0),
//         ...
//         fDir_Self      = fEntry,    ///< Top directory entry
//         fDir_Files     = (1 << 1),  ///< Non-directory entries
//         fDir_Subdirs   = (1 << 2),  ///< Subdirectory entries (non recursive)
//         fDir_Recursive = (1 << 3),  ///< Add a recursion
//         ...
//         fDir_All = fDir_Self + fDir_Files + fDir_Subdirs,
//         ...
//         eRecursive    = fDir_All | fDir_Recursive,
//     ...
//     virtual bool Remove(TRemoveFlags flags = eRecursive) const;
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2595-2603
// ```c++
// bool CDirEntry::Remove(TRemoveFlags flags) const
// {
//     // Is this a directory ? (and processing not entry only)
//     if ( (flags & (fDir_All | fDir_Recursive)) != eEntryOnly  &&  IsDir(eIgnoreLinks) ) {
//         return CDir(GetPath()).Remove(flags);
//     }
//     // Other entries
//     return RemoveEntry(flags);
// }
// ```
// `CFile` does not override `Remove`, so `CFile(name).Remove()` is this function with
// `eRecursive`. A failure is not reported: the file API logs only with `FileAPILogging`, whose
// default is false (ncbifile.cpp:176-180, 221-230), and `DeleteBlastDb` only counts the names
// removed.
/// `CDirEntry(path).Remove(eRecursive)`: removes `path`, a directory with its contents; `true`
/// when it is removed.
pub(crate) fn remove(path: &Path) -> bool {
    if is_dir_ignore_links(path) {
        return remove_dir(path);
    }
    remove_entry(path)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2203-2241
// ```c++
// CDirEntry::EType CDirEntry::GetType(EFollowLinks follow) const
// {
//     ...
// #else // NCBI_OS_UNIX
//     if (follow == eFollowLinks) {
//         errcode = stat(GetPath().c_str(), &st);
//     } else {
//         errcode = lstat(GetPath().c_str(), &st);
//     }
// #endif
//     if (errcode != 0) {
//         ...
//         return eUnknown;
//     }
//     return GetType(st);
// }
// ```
// `IsDir(eIgnoreLinks)`: a symbolic link is not a directory, whatever it points to, so the link
// is removed and its target is kept. On Windows NCBI's `GetType` calls `stat` for both values of
// `follow` (ncbifile.cpp:2208-2209); LOSAT does not follow a link there either (not compared
// with a makeblastdb run on Windows).
fn is_dir_ignore_links(path: &Path) -> bool {
    std::fs::symlink_metadata(path).is_ok_and(|meta| meta.is_dir())
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2606-2630
// ```c++
// bool CDirEntry::RemoveEntry(TRemoveFlags flags) const
// {
//
//     if ( NcbiSys_remove(_T_XCSTRING(GetPath())) != 0 ) {
//         switch (errno) {
//         case ENOENT:
//             if ( flags & fIgnoreMissing )
//                 return true;
//             break;
//
// #if defined(NCBI_OS_MSWIN)
//         case EACCES:
//             if ( NCBI_PARAM_TYPE(NCBI, DeleteReadOnlyFiles)::GetDefault() ) {
// ```
// `eRecursive` has no `fIgnoreMissing`, so a missing entry is a failure; `DeleteReadOnlyFiles`
// is false by default (ncbifile.cpp:146-148), so a read-only file on Windows is one too.
/// `CDirEntry::RemoveEntry`: `remove()` of an entry that is not a directory.
fn remove_entry(path: &Path) -> bool {
    std::fs::remove_file(path).is_ok()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:4342-4456
// ```c++
// bool CDir::Remove(TRemoveFlags flags) const
// {
//     ...
// #if !defined(NCBI_OS_MSWIN)
//     // Make directory writable for user to remove any entry inside
//     SetMode(CDirEntry::fWrite | static_cast<TMode>(CDirEntry::fModeAdd),
//             CDirEntry::fModeNoChange,
//             CDirEntry::fModeNoChange);
// #endif
//
//     // Read all entries in directory
//     unique_ptr<TEntries> contents(GetEntriesPtr());
//     if (!contents.get()) {
//         LOG_ERROR(72, "CDir::Remove(): Cannot get content of: " + GetPath());
//         return false;
//     }
//
//     bool success = true;
//     try {
//         // Remove each entry
//         ITERATE(TEntries, entry, *contents.get()) {
//             string name = (*entry)->GetName();
//             if (name == "." || name == ".." || name == string(1, GetPathSeparator())) {
//                 continue;
//             }
//             // Get entry item with full pathname
//             CDirEntry item(GetPath() + GetPathSeparator() + name);
//
//             if (flags & fDir_Recursive) {
//                 ...
//                 int f = (flags & fDir_Subdirs) ? (flags | fDir_Self) : flags;
//                 if (item.IsDir(eIgnoreLinks)) {
//                     if (!CDir(item.GetPath()).Remove(f)) {
//                         if (flags & fProcessAll) {
//                             success = false;
//                         } else {
//                             throw "Removing subdirectory failed";
//                         }
//                     }
//                 }
//                 else if (flags & fDir_Files) {
//                     if (!item.Remove(f)) {
//                         if (flags & fProcessAll) {
//                             success = false;
//                         } else {
//                             throw "Removing directory entry failed";
//                         }
//                     }
//                 }
//             }
//             ...
//         }
//         // Remove top directory
//         if ((flags & fDir_Self) && NcbiSys_rmdir(_T_XCSTRING(GetPath())) != 0) {
//             ...
//             if (flags & fProcessAll) {
//                 success = false;
//             } else {
//                 throw "Cannot remove directory entry";
//             }
//         }
//     }
//     catch (const char* what) {
//         ...
//         return false;
//     }
//     return success;
// }
// ```
// `eRecursive` has no `fProcessAll`: the first entry that cannot be removed stops the removal,
// and the entries before it are gone (makeblastdb 2.17.0 on Linux x64, 2026-11-08: a directory
// holding a file and an unreadable directory loses the file and keeps the other).
/// `CDir::Remove(eRecursive)` of the directory `path`.
fn remove_dir(path: &Path) -> bool {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:1061-1084 (CDirEntry::SetModeEntry)
    // ```c++
    //         TNcbiSys_stat st;
    //         if (NcbiSys_stat(_T_XCSTRING(GetPath()), &st) != 0) {
    //             ...
    //             return false;
    //         }
    //         ModeFromModeT(st.st_mode, &user, &group, &other);
    //     }
    //
    //     UPDATE_PERMS(user,    user_mode);
    //     UPDATE_PERMS(group,   group_mode);
    //     UPDATE_PERMS(other,   other_mode);
    //     UPDATE_PERMS(special, special_mode);
    //
    //     // change permissions
    //
    //     mode_t mode = MakeModeT(user, group, other, special);
    //
    //     if ( NcbiSys_chmod(_T_XCSTRING(GetPath()), mode) != 0 ) {
    // ```
    // `SetMode(fWrite | fModeAdd, fModeNoChange, fModeNoChange)` with the default special bits
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
    // ```c++
    //     DIR* dir = opendir(base_path.c_str());
    //     if ( !dir ) {
    //         CNcbiError::SetFromErrno();
    //         delete contents;
    //         ...
    //         return NULL;
    //     }
    //     while (struct dirent* entry = readdir(dir)) {
    //         if ( !IS_RECURSIVE_ENTRY  &&
    //              masks.Match(entry->d_name, use_case) ) {
    //             s_AddEntry(contents, base_path, entry, flags);
    //         }
    //     }
    // ```
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2456-2472,2475-2482,2510-2559
// (CDirEntry::Rename)
// ```c++
// bool CDirEntry::Rename(const string& newname, TRenameFlags flags)
// {
//     CDirEntry src(*this);
//     CDirEntry dst(newname);
//
//     // Dereference links
//     if ( F_ISSET(flags, fRF_FollowLinks) ) {
//         ...
//     }
//     // The source entry must exists
//     EType src_type = src.GetType();
//     if ( src_type == eUnknown ) {
//         LOG_ERROR_NCBI(17, "CDirEntry::Rename():"
//                            " Source path does not exist: " + src.GetPath(),
//                            CNcbiError::eNoSuchFileOrDirectory);
//         return false;
//     }
//
//     // Try to "move" in one atomic operation if possible to avoid race
//     // conditions between check on destination existence and renaming.
//
// #ifdef NCBI_OS_MSWIN
//     // On Windows we can try to move file or whole directory, even across volumes
//     if ( ::MoveFileEx(_T_XCSTRING(src.GetPath()),
//                       _T_XCSTRING(dst.GetPath()), MOVEFILE_COPY_ALLOWED) ) {
//         Reset(newname);
//         return true;
//     }
// #else
//     // On Unix we can use "link" technique for files.
//     ...
//     if ( src_type == eFile  &&
//          link(_T_XCSTRING(src.GetPath()),
//               _T_XCSTRING(dst.GetPath())) == 0 ) {
//         // Hard link successfully created, so we can just remove source file
//         if ( src.RemoveEntry() ) {
//             Reset(newname);
//             return true;
//         }
//     }
// #endif
//     // On error, symlink, or destination existence --
//     // continue on regular processing below
//
//     EType dst_type = dst.GetType();
//
//     // If destination exists...
//
//     if ( dst_type != eUnknown ) {
//         ...
//         // Overwrite destination entry
//         if ( F_ISSET(flags, fRF_Overwrite) ) {
//             if ( dst.Exists() ) {
//                 dst.Remove();
//             }
//         }
//     }
//
//     // On some platform rename() fails if destination entry exists,
//     // on others it can overwrite destination.
//     // For consistency return FALSE if destination already exists.
//     if ( dst.Exists() ) {
//         // this means Remove() has failed, and error is set already
//         LOG_ERROR(21, "CDirEntry::Rename(): Destination path exists: " + GetPath());
//         return false;
//     }
//
//     // Rename
//
//     if ( NcbiSys_rename(_T_XCSTRING(src.GetPath()),
//                         _T_XCSTRING(dst.GetPath())) != 0 ) {
// #ifdef NCBI_OS_MSWIN
//         if ( errno != EACCES ) {
// #else
//         if ( errno != EXDEV ) {
// #endif
//             LOG_ERROR_ERRNO(21, "CDirEntry::Rename(): rename() failed for " + GetPath());
//             return false;
//         }
//         // Note that rename() fails in the case of cross-device renaming.
//         // So, try to make a copy and remove the original later.
//         unique_ptr<CDirEntry> e(CDirEntry::CreateObject(src_type, src.GetPath()));
//         if ( !e->Copy(dst.GetPath(), fCF_Recursive | fCF_PreserveAll ) ) {
//             LOG_ERROR(102, "CDirEntry::Rename(): Renaming via Copy() failed for " + GetPath());
//             unique_ptr<CDirEntry> tmp(CDirEntry::CreateObject(src_type, dst.GetPath()));
//             tmp->Remove(eRecursive);
//             return false;
//         }
//         // Remove 'src'
//         if ( !e->Remove(eRecursive) ) {
//             LOG_ERROR(102, "CDirEntry::Rename(): Renaming via Copy() failed for " + GetPath());
//             // Do not delete 'dst' here because in case of directories the
//             // source may be already partially removed, so we can lost data.
//             return false;
//         }
//     }
//     Reset(newname);
//     return true;
// }
// ```
// `CWriteDB_File::RenameSingle` and `RenameFileIndex` call it with `fRF_Overwrite` alone (no
// `fRF_FollowLinks`, `fRF_EqualTypes`, `fRF_Update` or `fRF_Backup`) and ignore the result.
// The `LOG_ERROR` macros print only with `FileAPILogging`, whose default is false
// (ncbifile.cpp:176-180, 193-198). The type of an entry (`GetType`, `Exists`, ncbifile.hpp:3961)
// is that of `stat`, which follows a link: a source that is a dangling link does not exist, and
// a destination that is one does not exist either.
// A source that is a file reaches the hard-link branch first (Unix); its failure, because the
// destination exists, is what sends an existing destination to the removal below. The two
// names of a volume file are in one directory, so a cross-device `rename()` is not expected;
// its `Copy()` is ported for a file only (the mode is copied, the owner and the times of
// `fCF_PreserveAll` are not), and a source that is a directory fails there.
/// `CDirEntry(src).Rename(dst, CDirEntry::fRF_Overwrite)`: `true` when `src` got the name
/// `dst`. An existing destination is removed first (`remove`, a directory with its contents);
/// one that stays makes the rename fail, and `src` keeps its name.
pub(crate) fn rename_overwrite(src: &Path, dst: &Path) -> bool {
    let Ok(src_meta) = std::fs::metadata(src) else {
        return false;
    };
    #[cfg(unix)]
    if src_meta.is_file() && std::fs::hard_link(src, dst).is_ok() && remove_entry(src) {
        return true;
    }
    if std::fs::metadata(dst).is_ok() {
        // `dst.Exists()` asks what the `GetType()` above asked.
        let _ = remove(dst);
    }
    if std::fs::metadata(dst).is_ok() {
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
    if !src_meta.is_file() || std::fs::copy(src, dst).is_err() {
        let _ = remove(dst);
        return false;
    }
    remove(src)
}

#[cfg(test)]
mod tests {
    use super::{remove, rename_overwrite};
    use std::path::{Path, PathBuf};

    /// A directory of its own under the system's temporary directory.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(name: &str) -> Self {
            let dir = std::env::temp_dir().join(format!(
                "losat-dir-entry-{name}-{}-{:?}",
                std::process::id(),
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
        std::fs::set_permissions(ro.join("ro.txt"), std::fs::Permissions::from_mode(0o444))
            .unwrap();
        std::fs::set_permissions(ro.join("inner"), std::fs::Permissions::from_mode(0o555)).unwrap();
        std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
        assert!(remove(&ro));
        assert!(std::fs::symlink_metadata(&ro).is_err());

        // A directory that cannot be read stops the removal: it gets the user's write
        // permission only, the entries listed before it are gone and the others stay.
        let stuck = dir.0.join("db.ptf");
        std::fs::create_dir_all(stuck.join("locked")).unwrap();
        std::fs::write(stuck.join("locked").join("h"), b"h").unwrap();
        std::fs::set_permissions(stuck.join("locked"), std::fs::Permissions::from_mode(0o000))
            .unwrap();
        std::fs::set_permissions(&stuck, std::fs::Permissions::from_mode(0o3755)).unwrap();
        // Root reads any directory; the case then says nothing.
        if std::fs::read_dir(stuck.join("locked")).is_ok() {
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
        std::fs::set_permissions(dst.join("locked"), std::fs::Permissions::from_mode(0o000))
            .unwrap();
        // Root reads any directory; the case then says nothing.
        if std::fs::read_dir(dst.join("locked")).is_ok() {
            return;
        }
        assert!(!rename_overwrite(&src, &dst));
        assert_eq!(std::fs::read(&src).unwrap(), b"new");
        assert!(dst.is_dir());
    }

    #[cfg(unix)]
    #[test]
    fn links_are_judged_by_what_they_point_to() {
        let dir = TempDir::new("rename-link");
        // A source that is a link to nothing does not exist.
        let (src, dst) = (dir.0.join("db.00.pog"), dir.0.join("db.pog"));
        std::os::unix::fs::symlink("/nonexistent/target", &src).unwrap();
        assert!(!rename_overwrite(&src, &dst));
        assert!(
            std::fs::symlink_metadata(&src).is_ok() && std::fs::symlink_metadata(&dst).is_err()
        );
        // A destination that is a link to nothing does not exist either: the rename replaces it.
        std::fs::remove_file(&src).unwrap();
        std::fs::write(&src, b"s").unwrap();
        std::os::unix::fs::symlink("/nonexistent/target", &dst).unwrap();
        assert!(rename_overwrite(&src, &dst));
        assert_eq!(std::fs::read(&dst).unwrap(), b"s");
        // A source that is a link to a file moves as the link.
        let (link, new) = (dir.0.join("db.00.nog"), dir.0.join("db.nog"));
        std::os::unix::fs::symlink("db.pog", &link).unwrap();
        assert!(rename_overwrite(&link, &new));
        assert!(std::fs::symlink_metadata(&link).is_err());
        assert_eq!(std::fs::read_link(&new).unwrap(), Path::new("db.pog"));
    }
}

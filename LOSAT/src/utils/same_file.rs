//! Whether two names are one file, for the hosts that give no device and inode numbers: the
//! check that an `-export_search_strategy` file (`report::archive::strategy`) and a psiblast
//! `-out_pssm` file (`algorithm::psiblast::pssm`) are not the `-out` file.
//!
//! No NCBI counterpart: NCBI opens each of these files through a stream of its own and never
//! asks whether two of them are one file; LOSAT refuses the pairs whose bytes NCBI leaves in the
//! order of its buffers' flushes.
//!
//! Unix compares the device and inode of the files, and does not use this module.
//! Elsewhere (Windows, WASI) two names are the same file when they have the same canonical
//! path ([`canonical_path`]: [`std::fs::canonicalize`]). WASI has no `canonicalize` in Rust's
//! standard library (`Unsupported`) and its device and inode numbers are an unstable interface
//! (`std::os::wasi::fs::MetadataExt`, `wasi_ext`), so there the path is made absolute against
//! the current directory (`/` of the preopened file system until the program changes it) and
//! resolved name by name as `realpath(3)` resolves it, through symbolic links (line NP WW4b1,
//! 2026-11-27). That form sees every spelling of one path and its symbolic links, but not a hard
//! link to the file (a second name of the same file): such an alias of the `-out` file is not
//! refused on WASI.

use std::io::ErrorKind;
use std::path::{Component, Path, PathBuf};

/// The form of the existing file `path` that two names of one path share, or `None` for a path
/// that does not exist (`metadata` fails) or that cannot be put in that form.
#[cfg_attr(unix, allow(dead_code))] // unix compares device and inode numbers
pub(crate) fn path_identity(path: &Path) -> Option<PathBuf> {
    std::fs::metadata(path).ok()?;
    canonical_path(path)
}

/// The existing entry `path` resolved through symbolic links, as `realpath(3)` resolves it
/// ([`std::fs::canonicalize`]), or `None` where that fails. Where the platform has no
/// `canonicalize` (WASI preview 1, whose files a host gives: LOSAT Web's engine), the same
/// resolution made name by name ([`resolved_by_names`]), so that two names of one file or
/// directory are one path there too.
pub(crate) fn canonical_path(path: &Path) -> Option<PathBuf> {
    match std::fs::canonicalize(path) {
        Ok(canonical) => Some(canonical),
        Err(error) if error.kind() == ErrorKind::Unsupported => resolved_by_names(path),
        Err(_) => None,
    }
}

/// The most symbolic links that one resolution follows (Linux's `MAXSYMLINKS`; `ELOOP` past it).
const MAX_LINKS: usize = 40;

/// `realpath(3)` of `path` made from its names: absolute from the current directory, each name
/// looked up without following it, a symbolic link replaced by its target (from the root when the
/// target is absolute), `..` taking away the name before it once the links before it are
/// resolved; `None` where a name is missing, where a name or a final `/` follows a file
/// (`ENOTDIR`), or after `MAX_LINKS` links (`ELOOP`).
pub(crate) fn resolved_by_names(path: &Path) -> Option<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir().ok()?.join(path)
    };
    // The names still to resolve, the next one last.
    let mut pending: Vec<std::ffi::OsString> = Vec::new();
    let push = |pending: &mut Vec<std::ffi::OsString>, path: &Path| {
        // A final `/` asks for a directory: kept as a last `.`, which a file may not precede.
        let ends_in_separator = path
            .as_os_str()
            .as_encoded_bytes()
            .last()
            .is_some_and(|&byte| std::path::is_separator(byte.into()));
        if ends_in_separator {
            pending.push(".".into());
        }
        let names: Vec<_> = path
            .components()
            .filter_map(|component| match component {
                Component::Normal(name) => Some(name.to_os_string()),
                Component::ParentDir => Some("..".into()),
                Component::CurDir | Component::RootDir | Component::Prefix(_) => None,
            })
            .collect();
        pending.extend(names.into_iter().rev());
    };
    push(&mut pending, &absolute);
    let mut resolved = PathBuf::from("/");
    let mut links = 0;
    while let Some(name) = pending.pop() {
        if name == "." {
            continue;
        }
        if name == ".." {
            resolved.pop();
            continue;
        }
        let candidate = resolved.join(&name);
        let meta = std::fs::symlink_metadata(&candidate).ok()?;
        if !meta.file_type().is_symlink() {
            if !meta.is_dir() && !pending.is_empty() {
                return None;
            }
            resolved = candidate;
            continue;
        }
        links += 1;
        if links > MAX_LINKS {
            return None;
        }
        let target = std::fs::read_link(&candidate).ok()?;
        if target.is_absolute() {
            resolved = PathBuf::from("/");
        }
        push(&mut pending, &target);
    }
    Some(resolved)
}

/// Whether `log_path` (a file path that `-logfile` would open with `create(true).append(true)`)
/// and `input_path` (an input file path read by the run) refer to the same file or would create
/// and read the same missing file (line WQA D23, `Q-D23=A`).
///
/// - Directories and character devices (such as `/dev/null` or `/dev/stderr`) return `false`.
/// - When both paths exist: Unix compares `(dev, ino)` (handling symlinks, hard links, relative
///   paths, and case-insensitive file systems); Windows compares `GetFileInformationByHandle`
///   volume/index IDs with a `path_identity` fallback; WASI compares `path_identity`.
/// - When neither path exists yet: both paths resolve dangling symlinks (up to 40 hops) and
///   canonicalize their existing parent directory (falling back to `lexical_absolute`), so
///   `-logfile missing.fa` matching a missing `-query missing.fa` (or via a symlink directory or
///   dangling symlink) is also caught before `-logfile` creates the file.
pub(crate) fn same_file_target(log_path: &Path, input_path: &Path) -> bool {
    match (std::fs::metadata(log_path), std::fs::metadata(input_path)) {
        (Ok(log_meta), Ok(input_meta)) => {
            if log_meta.is_dir() || input_meta.is_dir() {
                return false;
            }
            #[cfg(unix)]
            {
                use std::os::unix::fs::{FileTypeExt, MetadataExt};
                if log_meta.file_type().is_char_device() {
                    return false;
                }
                (log_meta.dev(), log_meta.ino()) == (input_meta.dev(), input_meta.ino())
            }
            #[cfg(windows)]
            {
                if let (Some(a), Some(b)) = (windows_file_id(log_path), windows_file_id(input_path))
                {
                    return a == b;
                }
                path_identity(log_path)
                    .zip(path_identity(input_path))
                    .is_some_and(|(a, b)| a == b)
            }
            #[cfg(not(any(unix, windows)))]
            {
                path_identity(log_path)
                    .zip(path_identity(input_path))
                    .is_some_and(|(a, b)| a == b)
            }
        }
        (Err(_), Err(_)) => {
            let (Some(a), Some(b)) = (
                resolve_missing_target(log_path),
                resolve_missing_target(input_path),
            ) else {
                return false;
            };
            paths_equal(&a, &b)
        }
        _ => false,
    }
}

/// On Windows, the `(volume_serial, index_high, index_low)` file identifier of an existing file
/// via `GetFileInformationByHandle`, catching hard links as well as symlinks and 8.3 names.
#[cfg(windows)]
fn windows_file_id(path: &Path) -> Option<(u32, u32, u32)> {
    use std::os::windows::fs::OpenOptionsExt;
    use std::os::windows::io::AsRawHandle;
    use windows_sys::Win32::Foundation::HANDLE;
    use windows_sys::Win32::Storage::FileSystem::{
        GetFileInformationByHandle, BY_HANDLE_FILE_INFORMATION, FILE_SHARE_DELETE, FILE_SHARE_READ,
        FILE_SHARE_WRITE,
    };
    let file = std::fs::OpenOptions::new()
        .read(true)
        .share_mode(FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE)
        .open(path)
        .ok()?;
    let mut info = BY_HANDLE_FILE_INFORMATION::default();
    // SAFETY: `file` holds a valid open handle and `info` is a writable struct of the right type.
    let ok = unsafe { GetFileInformationByHandle(file.as_raw_handle() as HANDLE, &mut info) };
    (ok != 0).then_some((
        info.dwVolumeSerialNumber,
        info.nFileIndexHigh,
        info.nFileIndexLow,
    ))
}

/// Resolves a non-existent file path (possibly through dangling symlinks or a symlinked parent
/// directory) to its canonical parent directory joined with the final file name, or `None` if the
/// path ends in a directory component (`.`, `..`, or root) or exceeds the symlink hop limit.
fn resolve_missing_target(path: &Path) -> Option<PathBuf> {
    const MAX_HOPS: usize = 40;
    let mut current = path.to_path_buf();
    for _ in 0..MAX_HOPS {
        match std::fs::symlink_metadata(&current) {
            Ok(meta) if meta.file_type().is_symlink() => {
                let target = std::fs::read_link(&current).ok()?;
                current = if target.is_absolute() {
                    target
                } else {
                    match current.parent() {
                        Some(parent) if !parent.as_os_str().is_empty() => parent.join(target),
                        _ => target,
                    }
                };
            }
            _ => break,
        }
    }
    if std::fs::symlink_metadata(&current)
        .is_ok_and(|meta| meta.file_type().is_symlink() || meta.is_dir())
    {
        return None;
    }
    let file_name = current.file_name()?;
    let parent = match current.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    let resolved_parent = match canonical_path(parent) {
        Some(canonical) => canonical,
        None => lexical_absolute(parent)?,
    };
    Some(resolved_parent.join(file_name))
}

/// Compares two resolved non-existent target paths, folding ASCII case on Windows.
fn paths_equal(a: &Path, b: &Path) -> bool {
    #[cfg(windows)]
    {
        let (Some(a_parent), Some(b_parent)) = (a.parent(), b.parent()) else {
            return a == b;
        };
        let (Some(a_name), Some(b_name)) = (a.file_name(), b.file_name()) else {
            return a == b;
        };
        a_parent == b_parent
            && a_name
                .as_encoded_bytes()
                .eq_ignore_ascii_case(b_name.as_encoded_bytes())
    }
    #[cfg(not(windows))]
    {
        a == b
    }
}

/// `path` made absolute against the current directory, with `.` dropped and each `..` removing
/// the component before it (none at the root), by name only.
fn lexical_absolute(path: &Path) -> Option<PathBuf> {
    let joined = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir().ok()?.join(path)
    };
    let mut resolved = PathBuf::new();
    for component in joined.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                resolved.pop();
            }
            other => resolved.push(other),
        }
    }
    Some(resolved)
}

#[cfg(test)]
#[path = "same_file_tests.rs"]
mod tests;

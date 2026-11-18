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
//! path ([`std::fs::canonicalize`]). WASI has no `canonicalize` in Rust's standard library
//! (`Unsupported`) and its device and inode numbers are an unstable interface
//! (`std::os::wasi::fs::MetadataExt`, `wasi_ext`), so there the path is made absolute against
//! the current directory (`/` of the preopened file system until the program changes it) and
//! its `.` and `..` are resolved by name. That form sees every spelling of one path, but not a
//! symbolic link or a hard link to the file (a name of a second path to the same file): such an
//! alias of the `-out` file is not refused on WASI.

use std::io::ErrorKind;
use std::path::{Component, Path, PathBuf};

/// The form of the existing file `path` that two names of one path share, or `None` for a path
/// that does not exist (`metadata` fails) or that cannot be put in that form.
#[cfg_attr(unix, allow(dead_code))] // unix compares device and inode numbers
pub(crate) fn path_identity(path: &Path) -> Option<PathBuf> {
    std::fs::metadata(path).ok()?;
    match std::fs::canonicalize(path) {
        Ok(canonical) => Some(canonical),
        Err(error) if error.kind() == ErrorKind::Unsupported => lexical_absolute(path),
        Err(_) => None,
    }
}

/// `path` made absolute against the current directory, with `.` dropped and each `..` removing
/// the component before it (none at the root), by name only.
#[cfg_attr(unix, allow(dead_code))] // unix compares device and inode numbers
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
mod tests {
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
}

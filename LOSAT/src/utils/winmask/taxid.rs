//! blastn `-window_masker_taxid`: the counts file of a taxonomy ID (`WindowMaskerTaxidToDb`).

use std::ffi::OsString;
use std::path::{PathBuf, MAIN_SEPARATOR};

use super::CountsSource;

/// The file under `<directory>/<taxid>/`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:73
/// ```c++
/// static string s_WINDOW_MASKER_STAT_FILE_NAME("wmasker.obinary");
/// ```
const STAT_FILE_NAME: &str = "wmasker.obinary";

/// `WindowMaskerTaxidToDb(taxid)`: the counts file of `taxid` under `directory` (the registry
/// entry `[WINDOW_MASKER] WINDOW_MASKER_PATH`,
/// `ncbi_environment::WindowMaskerRegistry::path`; `None` or empty: the current
/// directory), or `None` when `source` has no file (or link to one) at that path (the file
/// system for the command line, the registered files for LOSAT Web). Only the name decides:
/// NCBI does not check the format of the file it finds (an ascii counts file under that name
/// is read as one), unlike `-window_masker_db`. `Err` when the current directory is needed
/// and cannot be read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:389-408
/// ```c++
///     if (retval == kEmptyStr) {
///         retval = CDir::GetCwd();
///     }
/// #if defined(NCBI_OS_MSWIN)
///     // We address this here otherwise CDirEntry::IsAbsolutePath() fails
///     if (NStr::StartsWith(retval, "//")) {
///         NStr::ReplaceInPlace(retval, "//", "\\\\");
///     }
/// #endif
///     return retval;
/// }
///
/// string WindowMaskerTaxidToDb(const string& window_masker_path, int taxid)
/// {
///     string path = window_masker_path;
///     path += CFile::GetPathSeparator() + NStr::IntToString(taxid)
///         + CFile::GetPathSeparator();
///     const string binpath = path + s_WINDOW_MASKER_STAT_FILE_NAME;
///     return (CFile(binpath).Exists() ? binpath : kEmptyStr);
/// }
/// ```
/// `CFile::Exists` is `IsFile()`, which follows links (ncbifile.hpp:4039-4042).
pub fn taxid_stat_file(
    directory: Option<&str>,
    taxid: i32,
    source: &dyn CountsSource,
) -> std::io::Result<Option<PathBuf>> {
    let directory: OsString = match directory.filter(|value| !value.is_empty()) {
        Some(value) => value.into(),
        None => std::env::current_dir()?.into_os_string(),
    };
    let path = stat_file_path(directory, taxid);
    Ok(source.is_file(&path).then_some(path))
}

/// `<directory>/<taxid>/wmasker.obinary` as NCBI joins the strings (a separator is added
/// after a directory that ends with one).
fn stat_file_path(mut directory: OsString, taxid: i32) -> PathBuf {
    if cfg!(windows) {
        if let Some(text) = directory.to_str() {
            if text.starts_with("//") {
                directory = text.replace("//", "\\\\").into();
            }
        }
    }
    directory.push(format!(
        "{MAIN_SEPARATOR}{taxid}{MAIN_SEPARATOR}{STAT_FILE_NAME}"
    ));
    PathBuf::from(directory)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::winmask::{FileCounts, MaskerParams, StatFile, UnitCounts, WinMaskError};
    use std::path::Path;
    use std::sync::Arc;

    // NCBI 2.17.0 blastn, 2026-10-19 (line NP M3 notes): `-window_masker_taxid -1` read
    // `<dir>/-1/wmasker.obinary`, `+7227`, `07227` and `0x1C3B` read `<dir>/7227/` (the path
    // holds the number in decimal), a trailing `/` and a link to the file worked, and a
    // directory named `wmasker.obinary` was no file (the "could not open" error).
    #[test]
    fn the_file_is_found_under_the_decimal_taxid() {
        let dir = std::env::temp_dir().join(format!("losat-winmask-taxid-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        for sub in ["7227", "-1", "9606/wmasker.obinary"] {
            std::fs::create_dir_all(dir.join(sub)).unwrap();
        }
        std::fs::write(dir.join("7227").join(STAT_FILE_NAME), b"x").unwrap();
        std::fs::write(dir.join("-1").join(STAT_FILE_NAME), b"x").unwrap();
        let root = dir.to_str().unwrap();
        let found =
            |directory: &str, taxid| taxid_stat_file(Some(directory), taxid, &FileCounts).unwrap();
        assert_eq!(
            found(root, 7227),
            Some(PathBuf::from(format!(
                "{root}{MAIN_SEPARATOR}7227{MAIN_SEPARATOR}wmasker.obinary"
            )))
        );
        assert_eq!(
            found(root, -1),
            Some(PathBuf::from(format!(
                "{root}{MAIN_SEPARATOR}-1{MAIN_SEPARATOR}wmasker.obinary"
            )))
        );
        assert_eq!(found(root, 9606), None);
        assert_eq!(found(root, 10090), None);
        assert_eq!(found(&format!("{root}{MAIN_SEPARATOR}missing"), 7227), None);
        let slash = format!("{root}{MAIN_SEPARATOR}");
        assert_eq!(
            found(&slash, 7227),
            Some(PathBuf::from(format!(
                "{slash}{MAIN_SEPARATOR}7227{MAIN_SEPARATOR}wmasker.obinary"
            )))
        );
        #[cfg(unix)]
        {
            std::fs::create_dir_all(dir.join("link/7227")).unwrap();
            std::os::unix::fs::symlink(
                dir.join("7227").join(STAT_FILE_NAME),
                dir.join("link/7227").join(STAT_FILE_NAME),
            )
            .unwrap();
            assert!(found(dir.join("link").to_str().unwrap(), 7227).is_some());
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn no_directory_is_the_current_one() {
        let cwd = std::env::current_dir().unwrap().into_os_string();
        let expected = stat_file_path(cwd, 7227);
        assert!(!expected.is_file());
        assert_eq!(taxid_stat_file(None, 7227, &FileCounts).unwrap(), None);
        assert_eq!(taxid_stat_file(Some(""), 7227, &FileCounts).unwrap(), None);
        assert_eq!(
            stat_file_path("/w".into(), 7227),
            PathBuf::from(format!(
                "/w{MAIN_SEPARATOR}7227{MAIN_SEPARATOR}wmasker.obinary"
            ))
        );
    }

    /// A source with files at some paths only (LOSAT Web's registered files).
    #[derive(Debug)]
    struct Paths(&'static [&'static str]);

    impl CountsSource for Paths {
        fn discover(&self, _: &str) -> Result<StatFile, WinMaskError> {
            unreachable!("only the paths are asked")
        }

        fn counts(&self, _: &str, _: &MaskerParams) -> Result<Arc<UnitCounts>, WinMaskError> {
            unreachable!("only the paths are asked")
        }

        fn is_file(&self, path: &Path) -> bool {
            self.0.iter().any(|known| Path::new(known) == path)
        }
    }

    #[test]
    fn the_source_says_which_files_are_there() {
        let source = Paths(&["windowmasker/9606/wmasker.obinary"]);
        assert_eq!(
            taxid_stat_file(Some("windowmasker"), 9606, &source).unwrap(),
            Some(PathBuf::from(format!(
                "windowmasker{MAIN_SEPARATOR}9606{MAIN_SEPARATOR}wmasker.obinary"
            )))
        );
        assert_eq!(
            taxid_stat_file(Some("windowmasker"), 10090, &source).unwrap(),
            None
        );
        assert_eq!(taxid_stat_file(Some("other"), 9606, &source).unwrap(), None);
    }
}

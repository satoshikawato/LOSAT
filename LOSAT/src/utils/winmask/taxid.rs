//! blastn `-window_masker_taxid`: the counts file of a taxonomy ID (`WindowMaskerTaxidToDb`).

use std::ffi::OsString;
use std::path::{PathBuf, MAIN_SEPARATOR};

use super::CountsSource;

/// The file under `<directory>/<taxid>/`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:73
const STAT_FILE_NAME: &str = "wmasker.obinary";

/// `WindowMaskerTaxidToDb` of `taxid`: the counts file of `taxid` under `directory` (the registry
/// entry `[WINDOW_MASKER] WINDOW_MASKER_PATH`,
/// `ncbi_environment::WindowMaskerRegistry::path`; `None` or empty: the current
/// directory), or `None` when `source` has no file (or link to one) at that path (the file
/// system for the command line, the registered files for LOSAT Web). Only the name decides:
/// NCBI does not check the format of the file it finds (an ascii counts file under that name
/// is read as one), unlike `-window_masker_db`. `Err` when the current directory is needed
/// and cannot be read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:389-408
/// `CFile::Exists` is `IsFile()`, which follows links (ncbifile.hpp:4039-4042).
pub fn taxid_stat_file(
    directory: Option<&[u8]>,
    taxid: i32,
    source: &dyn CountsSource,
) -> std::io::Result<Option<PathBuf>> {
    let directory: OsString = match directory.filter(|value| !value.is_empty()) {
        Some(value) => crate::utils::os_bytes::os_string_from_bytes(value.to_vec()),
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
#[path = "taxid_tests.rs"]
mod tests;

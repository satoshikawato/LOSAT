use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

static PATHS: Mutex<Vec<(PathBuf, PathBuf)>> = Mutex::new(Vec::new());

/// Lists `path` (a path of the test's own directory), to be read through `with`.
pub(crate) fn refuse(path: &Path, with: &Path) {
    PATHS
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .push((path.to_path_buf(), with.to_path_buf()));
}

/// The handle `DbFile::open` reads `path` through: `with`'s for a listed path.
pub(super) fn handle(path: &Path, file: File) -> File {
    let with = PATHS
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .iter()
        .find(|(listed, _)| listed == path)
        .map(|(_, with)| with.clone());
    match with {
        Some(with) => File::open(with).expect("the test's handle opens"),
        None => file,
    }
}

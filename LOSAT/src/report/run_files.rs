//! The files that a run writes besides its reports: the file of `-logfile` (line RG item E1w-3a).
//! The command line and the library write them to the file system (`RunFiles::Disk`); LOSAT Web,
//! which has none, hands each to its host (`RunFiles::Receiver`), as it hands over the side files
//! of outfmt 13 and 14 (`report::blast_xml2::SideFileSink`). No NCBI counterpart: NCBI writes
//! files.

use std::ffi::OsStr;
use std::io;
use std::sync::Arc;

use crate::blastinput::ncbi_environment::{diag_target, special_log_target, DiagTarget};

/// The receiver of the files of a run that a caller without files takes in place of the files
/// (LOSAT Web's host): the side files of outfmt 13 and 14 (`report::blast_xml2::SideFileSink`) and
/// the files of `RunFiles::Receiver`.
pub trait FileReceiver: Send + Sync {
    /// One file, whole: `name` is the path that the command line writes it to, `bytes` its
    /// contents.
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()>;
}

/// Where the files of a run go.
#[derive(Clone)]
pub enum RunFiles {
    /// The file system, at the names that the arguments give (the command line, the library).
    Disk,
    /// The caller's receiver, which gets each file when the run ends (LOSAT Web).
    Receiver(Arc<dyn FileReceiver>),
}

impl RunFiles {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608 (cited at
    // `ncbi_environment::diag_target`): no file for `""` and `/dev/null`, `cerr` for `-`, else the
    // file, or `cerr` when its directory does not exist.
    /// The target of the `-logfile` value `name` (`None` without one). On `Disk` the file
    /// system's (`diag_target`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), `cerr` for `.` and `..` (directories, which NCBI cannot open) and
    /// for a name with a directory (its host has none, and NCBI writes to `cerr` where the
    /// directory does not exist), and a file in memory for any other name, which `deliver_log`
    /// hands over. A run that searches refuses a name with a directory after the command line's
    /// steps before its search (`log_refusal`), so the help and the version of an argv with one
    /// are printed, as on the command line (line WQA's D14, 2026-11-21).
    pub fn log_target(&self, name: Option<&OsStr>) -> DiagTarget {
        match self {
            RunFiles::Disk => diag_target(name),
            RunFiles::Receiver(_) => match (special_log_target(name), name) {
                (Some(target), _) => target,
                (None, None) => DiagTarget::Stderr,
                // The working directory and its parent, which NCBI cannot open as a file, and a
                // file in a directory: `cerr`.
                (None, Some(name)) if name == "." || name == ".." || has_directory(name) => {
                    DiagTarget::Stderr
                }
                (None, Some(_)) => DiagTarget::Memory(Default::default()),
            },
        }
    }

    /// LOSAT Web's refusal of the `-logfile` value `name` of a run that searches: a name with a
    /// directory, which its host does not have (a `Receiver`'s; `None` on `Disk`, for NCBI's
    /// special names, `.` and `..`).
    pub fn log_refusal(&self, name: Option<&OsStr>) -> Option<String> {
        let name = name.filter(|name| special_log_target(Some(name)).is_none())?;
        (matches!(self, RunFiles::Receiver(_)) && has_directory(name)).then(|| {
            format!(
                "the -logfile name {} has a directory, which LOSAT Web does not have: give a file name without one",
                name.to_string_lossy()
            )
        })
    }

    /// Hands the file of the `-logfile` value `name`, which `target` holds (`log_target`), to the
    /// receiver. Nothing on `Disk` (the file is written) and for a target without a file.
    pub fn deliver_log(&self, name: Option<&OsStr>, target: DiagTarget) -> io::Result<()> {
        match (self, name, target) {
            (RunFiles::Receiver(receiver), Some(name), DiagTarget::Memory(bytes)) => receiver
                .receive(
                    name.as_encoded_bytes(),
                    &bytes
                        .into_inner()
                        .unwrap_or_else(std::sync::PoisonError::into_inner),
                ),
            _ => Ok(()),
        }
    }
}

/// A `-logfile` name that names a file in a directory (`.` and `..` are directories themselves).
fn has_directory(name: &OsStr) -> bool {
    name != "." && name != ".." && name.as_encoded_bytes().contains(&b'/')
}

#[cfg(test)]
#[path = "run_files_tests.rs"]
mod tests;

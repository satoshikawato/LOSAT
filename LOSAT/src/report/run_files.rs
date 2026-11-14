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
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608 (quoted at
    // `ncbi_environment::diag_target`): no file for `""` and `/dev/null`, `cerr` for `-`, else the
    // file, or `cerr` when its directory does not exist.
    /// The target of the `-logfile` value `name` (`None` without one). On `Disk` the file
    /// system's (`diag_target`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), `cerr` for `.` and `..` (directories, which NCBI cannot open), and a
    /// file in memory for any other name without a directory, which `deliver_log` hands over; a
    /// name with a directory is LOSAT Web's refusal (its host has no directories, where NCBI would
    /// write the messages to `cerr` instead).
    pub fn log_target(&self, name: Option<&OsStr>) -> Result<DiagTarget, String> {
        match self {
            RunFiles::Disk => Ok(diag_target(name)),
            RunFiles::Receiver(_) => match (special_log_target(name), name) {
                (Some(target), _) => Ok(target),
                (None, None) => Ok(DiagTarget::Stderr),
                // The working directory and its parent, which NCBI cannot open as a file: `cerr`.
                (None, Some(name)) if name == "." || name == ".." => Ok(DiagTarget::Stderr),
                (None, Some(name)) if name.as_encoded_bytes().contains(&b'/') => Err(format!(
                    "the -logfile name {} has a directory, which LOSAT Web does not have: give a file name without one",
                    name.to_string_lossy()
                )),
                (None, Some(_)) => Ok(DiagTarget::Memory(Vec::new())),
            },
        }
    }

    /// Hands the file of the `-logfile` value `name`, which `target` holds (`log_target`), to the
    /// receiver. Nothing on `Disk` (the file is written) and for a target without a file.
    pub fn deliver_log(&self, name: Option<&OsStr>, target: DiagTarget) -> io::Result<()> {
        match (self, name, target) {
            (RunFiles::Receiver(receiver), Some(name), DiagTarget::Memory(bytes)) => {
                receiver.receive(name.as_encoded_bytes(), &bytes)
            }
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    #[derive(Default)]
    struct Kept(Mutex<Vec<(Vec<u8>, Vec<u8>)>>);

    impl FileReceiver for Kept {
        fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
            self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
            Ok(())
        }
    }

    /// The names of `-logfile` of a receiver: NCBI's special names, a file in memory handed over
    /// under its name, and the refusal of a name with a directory.
    #[test]
    fn a_receiver_keeps_the_log_in_memory_and_refuses_directories() {
        let kept = Arc::new(Kept::default());
        let files = RunFiles::Receiver(kept.clone());
        let target = |name: Option<&str>| files.log_target(name.map(OsStr::new));
        assert!(matches!(target(None), Ok(DiagTarget::Stderr)));
        assert!(matches!(target(Some("-")), Ok(DiagTarget::Stderr)));
        assert!(matches!(target(Some("")), Ok(DiagTarget::Discard)));
        assert!(matches!(target(Some("/dev/null")), Ok(DiagTarget::Discard)));
        assert!(matches!(target(Some(".")), Ok(DiagTarget::Stderr)));
        assert!(matches!(target(Some("..")), Ok(DiagTarget::Stderr)));
        for name in ["d/run.log", "./run.log", "/tmp/run.log"] {
            assert_eq!(
                target(Some(name)).err().as_deref(),
                Some(
                    format!(
                        "the -logfile name {name} has a directory, which LOSAT Web does not have: give a file name without one"
                    )
                    .as_str()
                )
            );
        }
        let mut log = target(Some("run.log")).unwrap();
        log.write_text(b"posted\n", &mut io::sink());
        files.deliver_log(Some(OsStr::new("run.log")), log).unwrap();
        // A target without a file hands nothing over.
        files
            .deliver_log(Some(OsStr::new("-")), DiagTarget::Stderr)
            .unwrap();
        assert_eq!(
            *kept.0.lock().unwrap(),
            vec![(b"run.log".to_vec(), b"posted\n".to_vec())]
        );
    }
}

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
    // file, or `cerr` when its directory does not exist or it cannot be opened (a directory).
    /// The target of the `-logfile` value `name` (`None` without one). On `Disk` the file
    /// system's (`diag_target`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), and any other name as the command line on Linux reads it in a
    /// working directory where the name's directories exist (`log_name`), as LOSAT Web reads a
    /// `-conffile` name as a path of the command line's working directory (line WQA's sweep, D11,
    /// 2026-11-20; it refused a name with a directory before): `cerr` for a directory, which NCBI
    /// cannot open as a file; the file that is the standard error (`DiagTarget::StderrFile`) for
    /// its names; no output for another spelling of `/dev/null`; and for any other name a file in
    /// memory under the whole name, which `deliver_log` hands over (`logs/run.log`, `/tmp/run.log`,
    /// and `/dev/stdout`, whose bytes the command line writes into its standard output among the
    /// report's, where their place depends on the report's buffering).
    pub fn log_target(&self, name: Option<&OsStr>) -> DiagTarget {
        match self {
            RunFiles::Disk => diag_target(name),
            RunFiles::Receiver(_) => match (special_log_target(name), name) {
                (Some(target), _) => target,
                (None, None) => DiagTarget::Stderr,
                (None, Some(name)) => match log_name(name.as_encoded_bytes()) {
                    LogName::Directory => DiagTarget::Stderr,
                    LogName::StandardError => DiagTarget::StderrFile,
                    LogName::NullDevice => DiagTarget::Discard,
                    LogName::File => DiagTarget::Memory(Default::default()),
                },
            },
        }
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

/// What a `-logfile` name that is not one of NCBI's special names names, for a `Receiver`
/// (`RunFiles::log_target`).
#[derive(Debug, PartialEq, Eq)]
enum LogName {
    /// A directory, which NCBI cannot open as a file (`cerr`).
    Directory,
    /// The process's standard error.
    StandardError,
    /// The null device, written otherwise than NCBI's special `/dev/null`.
    NullDevice,
    /// A file of that name.
    File,
}

/// What the `-logfile` name `name` names on the command line on Linux, in a working directory
/// where the name's directories exist (No NCBI counterpart: LOSAT Web has no file system; line
/// WQA's sweep, D11, 2026-11-20). A name whose last component is empty (it ends with `/`), `.` or
/// `..` is a directory: NCBI's open of it fails, so its messages go to `cerr` (NCBI BLAST+ 2.17.0
/// and LOSAT's command line: `-logfile .`, `logs/`, `logs/.`). An absolute name is read
/// lexically (`//` and `/./` are `/`, `..` takes the component before it away): `/dev/stderr`,
/// `/dev/fd/2` and `/proc/self/fd/2` are the standard error, where LOSAT's command line writes the
/// file's bytes (NCBI's `-logfile /dev/stderr` too), and `/dev/null` is the null device. Every
/// other name is a file, made in its directory (`logs/run.log`, `./run.log`, `/tmp/run.log`);
/// so are the standard output's names (`/dev/stdout`).
fn log_name(name: &[u8]) -> LogName {
    let last = name.rsplit(|byte| *byte == b'/').next().unwrap_or_default();
    if matches!(last, b"" | b"." | b"..") {
        return LogName::Directory;
    }
    if name.first() != Some(&b'/') {
        return LogName::File;
    }
    let mut parts: Vec<&[u8]> = Vec::new();
    for part in name.split(|byte| *byte == b'/') {
        match part {
            b"" | b"." => {}
            b".." => {
                parts.pop();
            }
            part => parts.push(part),
        }
    }
    match parts.as_slice() {
        [b"dev", b"stderr"] | [b"dev", b"fd", b"2"] | [b"proc", b"self", b"fd", b"2"] => {
            LogName::StandardError
        }
        [b"dev", b"null"] => LogName::NullDevice,
        _ => LogName::File,
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

    /// The names of `-logfile` of a receiver: NCBI's special names, the names of directories and
    /// of the standard error, other spellings of `/dev/null`, and a file in memory, handed over
    /// under its whole name, also with a directory (line WQA's sweep, D11, 2026-11-20: refused
    /// before).
    #[test]
    fn a_receiver_reads_a_log_name_as_the_command_line_on_linux() {
        let kept = Arc::new(Kept::default());
        let files = RunFiles::Receiver(kept.clone());
        let target = |name: Option<&str>| files.log_target(name.map(OsStr::new));
        assert!(matches!(target(None), DiagTarget::Stderr));
        assert!(matches!(target(Some("-")), DiagTarget::Stderr));
        assert!(matches!(target(Some("")), DiagTarget::Discard));
        assert!(matches!(target(Some("/dev/null")), DiagTarget::Discard));
        for name in [
            ".", "..", "/", "logs/", "logs/.", "a/b/..", "/tmp/", "./", "//",
        ] {
            assert!(matches!(target(Some(name)), DiagTarget::Stderr), "{name}");
        }
        for name in [
            "/dev/stderr",
            "/dev/fd/2",
            "/proc/self/fd/2",
            "//dev/stderr",
            "/dev/./stderr",
            "/dev/fd/../fd/2",
        ] {
            assert!(
                matches!(target(Some(name)), DiagTarget::StderrFile),
                "{name}"
            );
        }
        for name in ["//dev/null", "/dev/./null", "/tmp/../dev/null"] {
            assert!(matches!(target(Some(name)), DiagTarget::Discard), "{name}");
        }
        let names = [
            "run.log",
            "logs/run.log",
            "./run.log",
            "a/../run.log",
            "/tmp/run.log",
            "/dev/stdout",
            "dev/stderr",
            "../run.log",
        ];
        for name in names {
            let log = target(Some(name));
            assert!(matches!(log, DiagTarget::Memory(_)), "{name}");
            log.write_text(name.as_bytes(), &mut io::sink());
            files.deliver_log(Some(OsStr::new(name)), log).unwrap();
        }
        // A target without a file hands nothing over.
        files
            .deliver_log(Some(OsStr::new("-")), DiagTarget::Stderr)
            .unwrap();
        files
            .deliver_log(Some(OsStr::new("/dev/stderr")), DiagTarget::StderrFile)
            .unwrap();
        assert_eq!(
            *kept.0.lock().unwrap(),
            names
                .iter()
                .map(|name| (name.as_bytes().to_vec(), name.as_bytes().to_vec()))
                .collect::<Vec<_>>()
        );
    }

    /// The file that is the standard error writes where `cerr` does, and gets the messages
    /// collected before it was opened, which `cerr` has already; a file in memory keeps them, and
    /// `cerr` gets none (line WQA's sweep, D11, 2026-11-20).
    #[test]
    fn the_standard_error_file_writes_to_the_standard_error() {
        use crate::blastinput::ncbi_environment::TextPlace;
        let mut stderr = Vec::new();
        let file = DiagTarget::StderrFile;
        assert_eq!(file.place(), TextPlace::StandardError);
        file.write_collected(b"collected\n", &mut stderr);
        file.write_text(b"posted\n", &mut stderr);
        file.write_to_file(b"error\n");
        assert_eq!(stderr, b"collected\nposted\n");
        let mut stderr = Vec::new();
        DiagTarget::Stderr.write_collected(b"collected\n", &mut stderr);
        assert!(stderr.is_empty());
        let memory = DiagTarget::Memory(Default::default());
        memory.write_collected(b"collected\n", &mut stderr);
        assert!(stderr.is_empty());
        let DiagTarget::Memory(bytes) = memory else {
            unreachable!()
        };
        assert_eq!(bytes.into_inner().unwrap(), b"collected\n");
    }
}

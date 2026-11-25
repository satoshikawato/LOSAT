//! Removal of the temporary files and directories a test makes (SW8): a `TempPath` removes its
//! path when it is dropped, so a test that passes leaves nothing in `temp_dir()`. A path that a
//! test keeps for inspection when it fails (a worker log, a directory named `..._kept`) is made
//! with `TempPath::keep_on_failure`, which keeps it when the test panics. A test that spawns a
//! child process drops its `TempPath` after the child has exited (`output()` and `status()` wait
//! for it).
#![allow(
    dead_code,
    clippy::duplicate_mod,
    reason = "each test file that declares this module uses a different subset of it, and the suites load it through several #[path] files"
)]

use std::ffi::OsStr;
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// A file or a directory (with what is in it) that is removed when the guard is dropped.
#[derive(Debug)]
pub struct TempPath {
    path: PathBuf,
    keep_on_failure: bool,
}

impl TempPath {
    /// `path` is removed when the guard is dropped, whether the test passed or not.
    pub fn new(path: impl Into<PathBuf>) -> TempPath {
        TempPath {
            path: path.into(),
            keep_on_failure: false,
        }
    }

    /// `path` is removed when the guard is dropped, unless the test panics: a failed test leaves
    /// it to be looked at.
    pub fn keep_on_failure(path: impl Into<PathBuf>) -> TempPath {
        TempPath {
            path: path.into(),
            keep_on_failure: true,
        }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Deref for TempPath {
    type Target = Path;

    fn deref(&self) -> &Path {
        &self.path
    }
}

impl AsRef<Path> for TempPath {
    fn as_ref(&self) -> &Path {
        &self.path
    }
}

impl AsRef<OsStr> for TempPath {
    fn as_ref(&self) -> &OsStr {
        self.path.as_os_str()
    }
}

impl Drop for TempPath {
    fn drop(&mut self) {
        if self.keep_on_failure && std::thread::panicking() {
            return;
        }
        remove(&self.path);
    }
}

/// Removes a file, or a directory with what is in it; a path that is not there is left as it is.
pub fn remove(path: &Path) {
    if path.is_dir() {
        let _ = std::fs::remove_dir_all(path);
    } else {
        let _ = std::fs::remove_file(path);
    }
}

/// A new empty directory `temp_dir()/<label>-<module path>-<pid>-<n>`, removed when it is dropped.
/// The tests of one binary run at the same time, so no two calls share a directory; a suite binary
/// holds one copy of this module per file that declares it, each with its own counter, so the
/// module path (`::` written as `-`) is part of the name.
pub fn unique_dir(label: &str) -> TempPath {
    new_dir(label, false)
}

fn new_dir(label: &str, keep_on_failure: bool) -> TempPath {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let base = std::env::temp_dir();
    // `create_dir` fails on a name that exists: a directory that an earlier test kept for
    // inspection (same process id and counter after a pid is reused) is never reused, and so is
    // never removed by a test that passes. The next number is taken instead.
    let path = loop {
        let path = base.join(format!(
            "{label}-{}-{}-{}",
            module_path!().replace("::", "-"),
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        match std::fs::create_dir(&path) {
            Ok(()) => break path,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => panic!("temporary directory {}: {error}", path.display()),
        }
    };
    if keep_on_failure {
        TempPath::keep_on_failure(path)
    } else {
        TempPath::new(path)
    }
}

/// A file named `name` in a new directory of its own (`unique_dir(label)`): the directory, with
/// the file, is removed when the value is dropped. The file is not made: the test writes it.
#[derive(Debug)]
pub struct TempFile {
    path: PathBuf,
    _dir: TempPath,
}

impl TempFile {
    pub fn new(label: &str, name: &str) -> TempFile {
        let dir = new_dir(label, false);
        TempFile {
            path: dir.join(name),
            _dir: dir,
        }
    }

    /// As `new`, but the directory is kept when the test panics (a worker log).
    pub fn keep_on_failure(label: &str, name: &str) -> TempFile {
        let dir = new_dir(label, true);
        TempFile {
            path: dir.join(name),
            _dir: dir,
        }
    }
}

impl Deref for TempFile {
    type Target = Path;

    fn deref(&self) -> &Path {
        &self.path
    }
}

impl AsRef<Path> for TempFile {
    fn as_ref(&self) -> &Path {
        &self.path
    }
}

impl AsRef<OsStr> for TempFile {
    fn as_ref(&self) -> &OsStr {
        self.path.as_os_str()
    }
}

/// The files of one run of the binary that the test makes only for it: an empty `HOME` and a
/// worker log. They are removed when the value is dropped; the log is kept when the test panics.
pub struct RunFiles {
    _home: TempPath,
    _log: TempFile,
}

/// Gives `command` an empty `HOME` and `USERPROFILE` (no `~/.ncbirc` names a taxonomy database or
/// a BLASTDB) and a worker log (`LOSAT_BLASTX_WORKER_LOG`: a binary built with the
/// `blastx-worker-probe` feature needs that file name as soon as a search has more than one
/// thread, and it selects nothing), both in directories of their own named after `label`.
pub fn home_and_worker_log(command: &mut Command, label: &str) -> RunFiles {
    let home = unique_dir(&format!("{label}_home"));
    command.env("HOME", &home).env("USERPROFILE", &home);
    let log = TempFile::keep_on_failure(&format!("{label}_worker"), "worker.log");
    command.env("LOSAT_BLASTX_WORKER_LOG", &log);
    RunFiles {
        _home: home,
        _log: log,
    }
}

/// What a run gave (`value` derefs to it) together with the run's files `files` (`RunFiles`, or a
/// worker log alone): the files stay until the caller drops this, so the caller's assertions on the
/// result are done before they are removed, and the log is kept when one of them fails.
pub struct Finished<T, F = RunFiles> {
    value: T,
    _files: F,
}

impl<T, F> Finished<T, F> {
    pub fn new(value: T, files: F) -> Finished<T, F> {
        Finished {
            value,
            _files: files,
        }
    }
}

impl<T, F> Deref for Finished<T, F> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.value
    }
}

/// When the test's environment names no `LOSAT_BLASTX_WORKER_LOG`, a worker log of this run alone
/// (a binary built with the `blastx-worker-probe` feature needs that file name as soon as a search
/// has more than one thread, and it selects nothing), set on `command`: the caller keeps the
/// returned guard until its assertions on the run are done; the log is removed when it is dropped
/// and kept when the test panics.
pub fn worker_log_if_unset(command: &mut Command, label: &str) -> Option<TempFile> {
    std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
        .is_none()
        .then(|| {
            let log = TempFile::keep_on_failure(label, "worker.log");
            command.env("LOSAT_BLASTX_WORKER_LOG", &log);
            log
        })
}

/// LOSAT's bytes of NCBI's frozen output or message: the one mapping of
/// `losat::wording::losat_output_bytes` (AGENTS.md rule 8; `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
/// exception 12).
pub fn expected_losat_bytes(frozen: Vec<u8>) -> Vec<u8> {
    match losat::wording::losat_output_bytes(&frozen) {
        std::borrow::Cow::Borrowed(_) => frozen,
        std::borrow::Cow::Owned(mapped) => mapped,
    }
}

/// The bytes of a binary seqid list (`.bsl`, as `blastdb_aliastool -seqid_file_in` writes it: the
/// byte 0, the 8-byte size of the file, the 8-byte number of ids, the 4-byte length of the title,
/// the title, the length byte of the date, the date, then the rest) with its date and the date's
/// length byte as `<DATE>`, and its size field zeroed when it holds the length of the file: the
/// date of NCBI's writer is 21 or 22 bytes by the day of the month, and the size counts it. A size
/// field that does not hold the length stays, so a wrong one changes the masked bytes. `None` for
/// bytes of another form. The frozen NCBI files of `non_utf8_names` and of
/// `blastdb_aliastool_fixtures` are compared so (`make_fixtures.py` of `non_utf8_names` the same).
pub fn undated_seqid_list(data: &[u8]) -> Option<Vec<u8>> {
    if data.len() <= 21 || data[0] != 0 {
        return None;
    }
    let date = 21 + u32::from_le_bytes(data[17..21].try_into().unwrap()) as usize;
    let length = *data.get(date)?;
    let end = (date + 1 + usize::from(length)).min(data.len());
    let size = u64::from_le_bytes(data[1..9].try_into().unwrap());
    let size: &[u8] = if size == data.len() as u64 {
        &[0; 8]
    } else {
        &data[1..9]
    };
    Some([&data[..1], size, &data[9..date], b"<DATE>", &data[end..]].concat())
}

/// [`expected_losat_bytes`] of text.
pub fn expected_losat_text(frozen: String) -> String {
    String::from_utf8(expected_losat_bytes(frozen.into_bytes())).expect("UTF-8 in, UTF-8 out")
}

/// [`expected_losat_text`] at the end of a chain: `read_to_string(..).unwrap().in_losat_words()`.
pub trait InLosatWords {
    fn in_losat_words(self) -> String;
}

impl InLosatWords for String {
    fn in_losat_words(self) -> String {
        expected_losat_text(self)
    }
}

//! The files of outfmt 13 and 14: the names of the side files and of the entries of the master
//! file, and the side files on the file system (DESIGN_OF5 1.7).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2316-2330
//! The side file of the n-th report is named by `s_GetBaseName` of `-out` (for the report's
//! format, with the directory), then `_`, n, and ".xml" or ".json"; the master file lists
//! `s_GetBaseName` of `-out` (the XML rule, no directory), then `_`, n, and the extension, for both
//! formats (blast_format.cpp:2432, :2449). So the master of an
//! outfmt 13 run lists `Y.Json_1.json` for `-out d/Y.Json` (the side file is `d/Y_1.json`) and
//! `k_1.json` for `-out d/k.xml` (the side file is `d/k.xml_1.json`): names of files that do not
//! exist. DECISIONS 8 keeps this as NCBI's defined behaviour.

use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::framing::{SideFiles, Xml2Kind};
use crate::report::run_files::FileReceiver;
use crate::utils::os_bytes::{path_bytes, path_from_bytes};

/// The characters `CDirEntry::SplitPath` takes as the end of the directory part: `/` on Unix;
/// `:`, `/` and `\` on Windows.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:358-375, :358-377
fn is_separator(byte: u8, windows: bool) -> bool {
    byte == b'/' || (windows && (byte == b'\\' || byte == b':'))
}

/// `s_GetBaseName`: the `-out` value without its extension when that is `.xml` (`xml`) or
/// `.json` (not `xml`) in any case, with its directory part (`with_dir`) or without. The
/// extension starts at the last `.` of the file name, so `.xml` alone is a base of nothing with
/// the extension `.xml`, and `a.b.json` has the extension `.json`.
pub fn base_name(out: &Path, xml: bool, with_dir: bool) -> Vec<u8> {
    base_name_of(&path_bytes(out), xml, with_dir, cfg!(windows))
}

/// `base_name` of the bytes of a path, with the separators of Windows (`windows`) or of Unix (the
/// unit tests check both on every platform).
fn base_name_of(path: &[u8], xml: bool, with_dir: bool, windows: bool) -> Vec<u8> {
    let name_at = path
        .iter()
        .rposition(|&byte| is_separator(byte, windows))
        .map_or(0, |at| at + 1);
    let (dir, file) = path.split_at(name_at);
    let (mut base, ext) = match file.iter().rposition(|&byte| byte == b'.') {
        Some(dot) => (file[..dot].to_vec(), &file[dot..]),
        None => (file.to_vec(), &file[file.len()..]),
    };
    let dropped: &[u8] = if xml { b".xml" } else { b".json" };
    if !ext.eq_ignore_ascii_case(dropped) {
        base.extend_from_slice(ext);
    }
    if with_dir {
        let mut name = dir.to_vec();
        name.extend_from_slice(&base);
        name
    } else {
        base
    }
}

/// The side files of a run written to the file system next to `-out`: the n-th report goes to
/// `<base>_<n>.xml` (outfmt 14) or `<base>_<n>.json` (outfmt 13), opened for writing and
/// truncated if it exists (`IOS_BASE::out`; an older file of another name stays as it was).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-775
pub struct FsSideFiles {
    /// `s_GetBaseName` of `out` for the format, with the directory: the side files' path up to the
    /// `_<n>`.
    path_base: Vec<u8>,
    /// `s_GetBaseName` of `out` by the XML rule, without the directory: the name the master lists,
    /// up to the `_<n>`.
    listed_base: Vec<u8>,
    extension: &'static str,
}

impl FsSideFiles {
    /// The side files of the run with `-out` `out` of the kind `kind` (`Files14` or `Files13`).
    pub fn new(out: &Path, kind: Xml2Kind) -> Self {
        let xml = kind != Xml2Kind::Files13;
        Self {
            path_base: base_name(out, xml, true),
            listed_base: base_name(out, true, false),
            extension: if xml { ".xml" } else { ".json" },
        }
    }

    fn name(&self, base: &[u8], n: u32) -> Vec<u8> {
        let mut name = base.to_vec();
        name.push(b'_');
        name.extend_from_slice(n.to_string().as_bytes());
        name.extend_from_slice(self.extension.as_bytes());
        name
    }

    /// The path of the n-th side file.
    pub fn path(&self, n: u32) -> PathBuf {
        path_from_bytes(self.name(&self.path_base, n))
    }

    /// The path of the n-th side file as bytes.
    fn path_bytes(&self, n: u32) -> Vec<u8> {
        self.name(&self.path_base, n)
    }
}

impl SideFiles for FsSideFiles {
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>> {
        Ok(Box::new(BufWriter::new(std::fs::File::create(
            self.path(n),
        )?)))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        self.name(&self.listed_base, n)
    }

    /// Writes the n-th report to its file as NCBI's object stream and `CNcbiOfstream` write it
    /// (`ostream::write_side_file`): where a write fails, NCBI's pieces of the document
    /// (`out_buffer::xml_pieces`, `json_pieces`) through the file's buffer decide NCBI's end.
    fn write_file(&mut self, n: u32, report: &[u8]) -> anyhow::Result<()> {
        use crate::report::serial::out_buffer::{json_pieces, xml_pieces, Stream};
        let file = std::fs::File::create(self.path(n))
            .map_err(|_| super::framing::side_file_open_error())?;
        let xml = self.extension == ".xml";
        let stream = if xml { Stream::Xml } else { Stream::Json };
        crate::report::ostream::write_side_file(
            file,
            report,
            || {
                if xml {
                    xml_pieces(report, super::writer::is_string_element)
                } else {
                    json_pieces(report)
                }
            },
            stream.back_limit(),
        )
    }

    /// Opens the side file of the n-th report as NCBI does when it writes the report, without
    /// leaving a mark: a path that exists is opened for writing and not truncated (a directory or a
    /// file without write permission fails), a path that does not exist is created and removed. The
    /// run that holds its side files makes this check when the report is made, so that an
    /// unopenable file stops it where NCBI's run stops, not at the epilog after every later batch
    /// was searched (DECISIONS 35 S1). The bytes are written when the hold ends, which creates the
    /// file again.
    fn check_openable(&self, n: u32) -> io::Result<()> {
        let path = self.path(n);
        match std::fs::OpenOptions::new().write(true).open(&path) {
            Ok(_) => Ok(()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                // A dangling link is created through (`File::create` follows it); only a path that
                // did not exist at all is removed again.
                let existed = path.symlink_metadata().is_ok();
                std::fs::File::create(&path)?;
                if !existed {
                    std::fs::remove_file(&path)?;
                }
                Ok(())
            }
            Err(error) => Err(error),
        }
    }
}

/// The side files of a run given to the caller's receiver, named as `FsSideFiles` names them for a
/// `-out` of the caller's choosing; the master file lists them by the same rule
/// (`FsSideFiles::listed`). A `-out` without a directory and without the extension `.xml` or
/// `.json` gives the side files the names that the master lists, in both formats.
#[derive(Clone)]
pub struct SideFileSink {
    out: PathBuf,
    receiver: Arc<dyn FileReceiver>,
}

impl SideFileSink {
    pub fn new(out: impl Into<PathBuf>, receiver: Arc<dyn FileReceiver>) -> Self {
        Self {
            out: out.into(),
            receiver,
        }
    }

    /// The `-out` whose name the side files and the master's entries follow.
    pub fn out(&self) -> &Path {
        &self.out
    }

    /// The side files of a run of the kind `kind` (`Files14` or `Files13`).
    pub fn side_files(&self, kind: Xml2Kind) -> ReceivedSideFiles {
        ReceivedSideFiles {
            names: FsSideFiles::new(&self.out, kind),
            receiver: self.receiver.clone(),
        }
    }
}

/// `SideFiles` that hand each side file to a `FileReceiver` (`SideFileSink::side_files`).
pub struct ReceivedSideFiles {
    names: FsSideFiles,
    receiver: Arc<dyn FileReceiver>,
}

impl SideFiles for ReceivedSideFiles {
    /// A writer that hands the side file over when it is flushed or dropped (the trait's other
    /// callers; a run writes through `write_file`).
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>> {
        Ok(Box::new(Received {
            receiver: &*self.receiver,
            name: self.names.path_bytes(n),
            bytes: Vec::new(),
            delivered: false,
        }))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        self.names.listed(n)
    }

    /// The report goes to the receiver as it is, without a copy; a receiver that fails is the
    /// side-file write error of every program (`side_file_write_error`).
    fn write_file(&mut self, n: u32, report: &[u8]) -> anyhow::Result<()> {
        self.receiver
            .receive(&self.names.path_bytes(n), report)
            .map_err(|_| super::framing::side_file_write_error())
    }
}

/// One side file for a receiver from `create`: the bytes are kept and handed over at the first
/// flush, or when it is dropped without one.
struct Received<'a> {
    receiver: &'a dyn FileReceiver,
    name: Vec<u8>,
    bytes: Vec<u8>,
    delivered: bool,
}

impl Received<'_> {
    fn deliver(&mut self) -> io::Result<()> {
        if std::mem::replace(&mut self.delivered, true) {
            return Ok(());
        }
        let bytes = std::mem::take(&mut self.bytes);
        self.receiver.receive(&self.name, &bytes)
    }
}

impl Write for Received<'_> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.delivered {
            return Err(io::Error::other("the side file was handed over"));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.deliver()
    }
}

impl Drop for Received<'_> {
    fn drop(&mut self) {
        let _ = self.deliver();
    }
}

/// The side files of a run of `kind`: the caller's receiver (`sink`) when it gives one, else the
/// files next to `out` (the run's `-out`), else none (outfmt 15 and 16, and a caller with neither,
/// which then fails at its first report).
pub fn run_side_files(
    kind: Xml2Kind,
    out: Option<&Path>,
    sink: Option<&SideFileSink>,
) -> super::framing::RunSideFiles {
    match (sink, out) {
        _ if !kind.is_files() => Box::new(super::framing::NoSideFiles),
        (Some(sink), _) => Box::new(sink.side_files(kind)),
        (None, Some(out)) => Box::new(FsSideFiles::new(out, kind)),
        (None, None) => Box::new(super::framing::NoSideFiles),
    }
}

#[cfg(test)]
#[path = "files_tests.rs"]
mod tests;

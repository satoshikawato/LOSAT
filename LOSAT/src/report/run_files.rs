//! The files that a run writes besides its reports: the file of `-logfile` (line RG item E1w-3a)
//! and the output files that NCBI's tools open by name (`-out FILE` of windowmasker, line NP WW3a).
//! The command line and the library write them to the file system (`RunFiles::Disk`); LOSAT Web,
//! which has none, hands each to its host (`RunFiles::Receiver`), as it hands over the side files
//! of outfmt 13 and 14 (`report::blast_xml2::SideFileSink`). No NCBI counterpart: NCBI writes
//! files.

use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, Write};
use std::sync::{Arc, Mutex, PoisonError};

use crate::blastinput::ncbi_environment::{diag_target, special_log_target, DiagTarget};

/// The receiver of the files of a run that a caller without files takes in place of the files
/// (LOSAT Web's host): the side files of outfmt 13 and 14 (`report::blast_xml2::SideFileSink`) and
/// the files of `RunFiles::Receiver`.
pub trait FileReceiver: Send + Sync {
    /// One file, whole: `name` is the path that the command line writes it to, `bytes` its
    /// contents.
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()>;

    /// Opens the file `name` of a tool's run (`OpenFiles`: `-out FILE`, `-logfile FILE`) to take
    /// the run's writes as the run makes them (`FileStream`); the file is complete when the stream
    /// is dropped, when the run ends. `None`, the default: the receiver takes each file whole, which
    /// the run keeps in memory and hands to `receive` when it ends (the files of a search, whose
    /// bytes and order this keeps). LOSAT Web's run of a tool streams its files, which may be larger
    /// than the browser's memory (line NP WW3a, D-WW3-6, 2026-11-22).
    fn open(&self, _name: &[u8]) -> Option<Box<dyn FileStream>> {
        None
    }
}

/// A file of a tool's run that the receiver takes as the run writes it (`FileReceiver::open`): the
/// writes and the changes of length of the one file of that path, in the order the run makes them,
/// as the system applies them to the file (line NP WW3a, D-WW3-9, 2026-11-23). Most writes come
/// at the end of what came before; a run that opens the file twice (`-out x -logfile x`) writes
/// through each opening at its own position, which may be before the end.
pub trait FileStream: Send {
    /// Writes `bytes` at `offset` of the file; past its end, the bytes between are zero (a hole).
    fn write_at(&mut self, offset: u64, bytes: &[u8]) -> io::Result<()>;

    /// Makes the file `length` bytes long: cut, or with zero bytes added. The run sends it when it
    /// opens the file again to write it, which empties it (length 0).
    fn set_len(&mut self, length: u64) -> io::Result<()>;
}

/// Where the files of a run go.
#[derive(Clone)]
pub enum RunFiles {
    /// The file system, at the names that the arguments give (the command line, the library).
    Disk,
    /// The caller's receiver, which gets each file whole when it is complete, or as it is written
    /// (`FileReceiver::open`, a tool's run; LOSAT Web).
    Receiver(Arc<dyn FileReceiver>),
}

impl std::fmt::Debug for RunFiles {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            RunFiles::Disk => "Disk",
            RunFiles::Receiver(_) => "Receiver",
        })
    }
}

impl RunFiles {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608 (cited at
    // `ncbi_environment::diag_target`): no file for `""` and `/dev/null`, `cerr` for `-`, else the
    // file, or `cerr` when its directory does not exist or it cannot be opened (a directory).
    /// The target of the `-logfile` value `name` (`None` without one). On `Disk` the file
    /// system's (`diag_target`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), and any other name as the command line of a user on Linux reads it
    /// in a working directory of the user's own (`log_name`), as LOSAT Web reads a `-conffile`
    /// name as a path of the command line's working directory (line WQA's sweep, D11, 2026-11-20;
    /// it refused a name with a directory before): `cerr` for a directory, which NCBI cannot open
    /// as a file, and for a name that a user's command line cannot open (`//run.log`,
    /// `/etc/run.log`, `/data/run.log`, `/dev/kmsg`, a component of 256 bytes); the file that is
    /// the standard error (`DiagTarget::StderrFile`) for `/dev/stderr`; no output for another
    /// spelling of `/dev/null`, for a device whose writes are lost (`/dev/ptmx`) and for a name in
    /// a file system that reports no free space (`/dev/fd/2`); and for any other name, one in the
    /// working directory, in the user's home or in a directory where every user makes files, a
    /// file in memory under the whole name, which `deliver_log` hands over (`logs/run.log`,
    /// `logs/`, `/tmp/run.log`, `/var/lock/run.log`, and `/dev/stdout`, whose bytes the command
    /// line writes into its standard output among the report's, where their place depends on the
    /// report's buffering). A tool's run opens its files through `OpenFiles`, which reads the
    /// names the same way (`received_log_target`).
    pub fn log_target(&self, name: Option<&OsStr>) -> DiagTarget {
        match self {
            RunFiles::Disk => diag_target(name),
            RunFiles::Receiver(_) => {
                received_log_target(name, |_| DiagTarget::Memory(Default::default()))
            }
        }
    }

    /// Hands the file of the `-logfile` value `name`, which `target` holds (`log_target`), to the
    /// receiver. Nothing on `Disk` (the file is written) and for a target without a file.
    pub fn deliver_log(&self, name: Option<&OsStr>, target: DiagTarget) -> io::Result<()> {
        match (self, name, target) {
            (RunFiles::Receiver(receiver), Some(name), DiagTarget::Memory(bytes)) => receiver
                .receive(
                    name.as_encoded_bytes(),
                    &bytes.into_inner().unwrap_or_else(PoisonError::into_inner),
                ),
            _ => Ok(()),
        }
    }
}

/// The target of the `-logfile` value `name` for a `Receiver` (`RunFiles::log_target`,
/// `OpenFiles::log_target`): NCBI's special names, then `log_name`'s reading; `file` gives the
/// target of a name that is a file (`LogName::File`), from the name.
fn received_log_target(name: Option<&OsStr>, file: impl FnOnce(&[u8]) -> DiagTarget) -> DiagTarget {
    match (special_log_target(name), name) {
        (Some(target), _) => target,
        (None, None) => DiagTarget::Stderr,
        (None, Some(name)) => match log_name(name.as_encoded_bytes()) {
            LogName::Unopened => DiagTarget::Stderr,
            LogName::StandardError => DiagTarget::StderrFile,
            LogName::Nowhere => DiagTarget::Discard,
            LogName::File => file(name.as_encoded_bytes()),
        },
    }
}

/// The files that one run of a tool opens by name for writing (`ncbi_app::run_app`: the file of
/// `-logfile` and the outputs, `-out FILE`), as the system opens them (line NP WW3a, D-WW3-9,
/// 2026-11-23). On `Disk` the file system's: `-logfile` as `diag_target` opens it, an output made
/// or emptied by `File::create`, as before. For a `Receiver`, each name as the command line of a
/// user on Linux opens it in the run's working directory (`log_name`, `output_name`), and every
/// name of one path is one file (`file_key`: `x`, `./x` and `sub/../x` are one), shared as the
/// system shares a file that a process opens twice (`-out x -logfile x`): each opening writes at
/// its own position, the output's opening empties the file, the log's appends to its end. A file
/// goes to the receiver as it is written (`FileReceiver::open`) or whole when the run ends; the
/// run's files end when the run ends (when the last clone of this and every opening of a file is
/// dropped), in the order the run opened them.
#[derive(Clone)]
pub(crate) struct OpenFiles(Option<Arc<FileTable>>);

impl std::fmt::Debug for OpenFiles {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(if self.0.is_some() { "Receiver" } else { "Disk" })
    }
}

/// The files that a `Receiver`'s run opened (`OpenFiles`), in the order it opened them, by their
/// paths (`file_key`).
struct FileTable {
    receiver: Arc<dyn FileReceiver>,
    files: Mutex<Vec<(Vec<u8>, Arc<Mutex<SharedFile>>)>>,
}

/// How an opening writes to its file (`FileTable::open`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Opening {
    /// To its end, for every write: the file of `-logfile`, which NCBI opens to append to it.
    Append,
    /// From its start, at its own position, after it emptied the file: an output (`create`).
    Create,
}

impl OpenFiles {
    /// The files of a run whose files go to `files`.
    pub(crate) fn new(files: &RunFiles) -> Self {
        OpenFiles(match files {
            RunFiles::Disk => None,
            RunFiles::Receiver(receiver) => Some(Arc::new(FileTable {
                receiver: receiver.clone(),
                files: Mutex::new(Vec::new()),
            })),
        })
    }

    /// The files of a run on the file system (`RunFiles::Disk`).
    pub(crate) fn disk() -> Self {
        OpenFiles(None)
    }

    /// The target of the `-logfile` value `name` (`RunFiles::log_target`), whose file, for a
    /// `Receiver`, is an opening of the run's file of that path that appends to it
    /// (`DiagTarget::Stream`).
    pub(crate) fn log_target(&self, name: Option<&OsStr>) -> DiagTarget {
        match &self.0 {
            None => diag_target(name),
            Some(table) => received_log_target(name, |bytes| {
                let key = log_key(bytes).expect("a file of log_name has a path");
                DiagTarget::Stream(Mutex::new(Box::new(table.open(
                    bytes,
                    key,
                    Opening::Append,
                ))))
            }),
        }
    }

    /// The output file `name` that a tool opens by name for writing (NCBI's `CNcbiOfstream` of the
    /// name: `-out FILE` of windowmasker, which made or emptied the file; the callers take `-` and
    /// the empty name first). On `Disk` the file of the file system, made or emptied
    /// (`File::create`, as before); `Err` where it does not open. For a `Receiver`, the name as the
    /// command line of a user on Linux opens it in the run's working directory (`output_name`): the
    /// run's file of that path, emptied (a file that the run opened before under that path or
    /// another name of it, `file_key`, is that file); the run's standard output for a name of it
    /// (`/dev/stdout`); writes that are lost (the null device) or that fail (`/dev/full`); `Err`
    /// where that command line cannot open it (a directory, a folder that the run's directory does
    /// not have, a directory of root's); and LOSAT's refusal of the names whose writes LOSAT Web
    /// does not reproduce (`Created::Refused`).
    pub(crate) fn create(&self, name: &str) -> Result<Created, Unopened> {
        let Some(table) = &self.0 else {
            return File::create(name)
                .map(|file| Created::File(OpenedFile::Disk(file)))
                .map_err(|_| Unopened);
        };
        match output_name(name.as_bytes(), &run_directory_entry) {
            OutputName::Unopened => Err(Unopened),
            OutputName::File => {
                let key = output_key(name.as_bytes()).expect("a file of output_name has a path");
                Ok(Created::File(OpenedFile::Shared(table.open(
                    name.as_bytes(),
                    key,
                    Opening::Create,
                ))))
            }
            OutputName::Stdout => Ok(Created::Stdout),
            OutputName::Discard => Ok(Created::File(OpenedFile::Discard)),
            OutputName::Full => Ok(Created::File(OpenedFile::Full)),
            OutputName::Special => Ok(Created::Refused(format!(
                "an output file '{name}' that is the standard error, a device or a file of \
                 procfs (LOSAT Web writes no such file)"
            ))),
            OutputName::Existing => Ok(Created::Refused(format!(
                "an output file '{name}' that is a file of the run's directory (the command \
                 line empties it; LOSAT Web reads the run's files and writes none of them)"
            ))),
        }
    }
}

impl FileTable {
    /// An opening of the run's file at the path `key` (`file_key`), named `name` where the run
    /// opens it first: the file that the run opened before at that path, which an output's opening
    /// empties (`O_TRUNC`), or a new file, which the receiver takes as it is written or whole.
    fn open(&self, name: &[u8], key: Vec<u8>, opening: Opening) -> SharedOpening {
        let mut files = self.files.lock().unwrap_or_else(PoisonError::into_inner);
        let file = match files.iter().find(|(path, _)| *path == key) {
            Some((_, file)) => {
                if opening == Opening::Create {
                    // A receiver's stream takes every change of length (LOSAT Web's emits it).
                    let _ = file
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner)
                        .set_len(0);
                }
                file.clone()
            }
            None => {
                let held = match self.receiver.open(name) {
                    Some(stream) => Held::Stream(stream),
                    None => Held::Memory {
                        name: name.to_vec(),
                        bytes: Vec::new(),
                        receiver: self.receiver.clone(),
                    },
                };
                let file = Arc::new(Mutex::new(SharedFile { length: 0, held }));
                files.push((key, file.clone()));
                file
            }
        };
        SharedOpening {
            file,
            position: 0,
            opening,
        }
    }
}

/// A file of a `Receiver`'s run (`FileTable`): its length and where its bytes go. It ends when its
/// last holder drops it (the run's table and its openings), at the end of the run: the receiver's
/// stream is dropped, or the file in memory goes to `receive`.
struct SharedFile {
    length: u64,
    held: Held,
}

/// Where the bytes of a `SharedFile` go.
enum Held {
    /// To the receiver as they are written (`FileReceiver::open`).
    Stream(Box<dyn FileStream>),
    /// Into memory, for the receiver's `receive` when the file ends.
    Memory {
        name: Vec<u8>,
        bytes: Vec<u8>,
        receiver: Arc<dyn FileReceiver>,
    },
}

impl SharedFile {
    /// Writes `bytes` at `offset`, as the system writes them to the file.
    fn write_at(&mut self, offset: u64, bytes: &[u8]) -> io::Result<()> {
        let end = offset + bytes.len() as u64;
        match &mut self.held {
            Held::Stream(stream) => stream.write_at(offset, bytes)?,
            Held::Memory { bytes: kept, .. } => {
                let (start, end) = (to_index(offset)?, to_index(end)?);
                if kept.len() < end {
                    kept.resize(end, 0);
                }
                kept[start..end].copy_from_slice(bytes);
            }
        }
        self.length = self.length.max(end);
        Ok(())
    }

    /// Makes the file `length` bytes long (`FileStream::set_len`).
    fn set_len(&mut self, length: u64) -> io::Result<()> {
        match &mut self.held {
            Held::Stream(stream) => stream.set_len(length)?,
            Held::Memory { bytes, .. } => bytes.resize(to_index(length)?, 0),
        }
        self.length = length;
        Ok(())
    }
}

/// An offset of a file in memory, which fits the address space.
fn to_index(offset: u64) -> io::Result<usize> {
    usize::try_from(offset).map_err(|_| io::Error::from(io::ErrorKind::OutOfMemory))
}

impl Drop for SharedFile {
    fn drop(&mut self) {
        if let Held::Memory {
            name,
            bytes,
            receiver,
        } = &self.held
        {
            // LOSAT Web's receiver keeps every file it gets (it fails for none).
            let _ = receiver.receive(name, bytes);
        }
    }
}

/// An opening of a `Receiver`'s file by the run (`FileTable::open`), which writes as the system
/// writes through a descriptor of the file: at the end of the file for an `Opening::Append`, else
/// at its own position, which each write moves on.
pub(crate) struct SharedOpening {
    file: Arc<Mutex<SharedFile>>,
    position: u64,
    opening: Opening,
}

impl Write for SharedOpening {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut file = self.file.lock().unwrap_or_else(PoisonError::into_inner);
        let at = match self.opening {
            Opening::Append => file.length,
            Opening::Create => self.position,
        };
        file.write_at(at, buf)?;
        self.position = at + buf.len() as u64;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// An output file that `OpenFiles::create` could not open: the command line's open fails.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Unopened;

/// What `OpenFiles::create` opened.
pub(crate) enum Created {
    /// A file, written by the tool.
    File(OpenedFile),
    /// The run's standard output (a `Receiver`'s `/dev/stdout`): the tool writes the file's bytes
    /// there, as the command line's file of that name does.
    Stdout,
    /// LOSAT's refusal of the name, the reason (a `Receiver`'s names that LOSAT Web does not write).
    Refused(String),
}

/// An output file that `OpenFiles::create` opened.
pub(crate) enum OpenedFile {
    /// A file of the file system (`RunFiles::Disk`).
    Disk(File),
    /// A `Receiver`'s file, through this opening of it.
    Shared(SharedOpening),
    /// A `Receiver`'s file whose writes are lost (the null device).
    Discard,
    /// A `Receiver`'s file whose writes fail for want of space (`/dev/full`).
    Full,
}

impl Write for OpenedFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OpenedFile::Disk(file) => file.write(buf),
            OpenedFile::Shared(file) => file.write(buf),
            OpenedFile::Discard => Ok(buf.len()),
            OpenedFile::Full => Err(io::Error::from(io::ErrorKind::StorageFull)),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            OpenedFile::Disk(file) => file.write_all(buf),
            _ => self.write(buf).map(|_| ()),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OpenedFile::Disk(file) => file.flush(),
            _ => Ok(()),
        }
    }
}

/// The path of a file of a run's working directory `cwd` (`working_directory`), without links,
/// as one name of the run's files (`OpenFiles`): every name that the system opens at that path
/// gives it (`x`, `./x`, `sub/../x`, `/proc/self/cwd/x`). A path in the working directory's tree
/// (`WORKING_TREE`) is told by its place from the working directory (how many directories above
/// it, then the rest), as each name has a working directory of its own depth; any other path is
/// itself.
fn file_key(path: &[u8], cwd: &[u8]) -> Vec<u8> {
    let depth = cwd.len() / WORKING_TREE.len();
    let (mut rest, mut inside) = (path, 0);
    while inside < depth {
        match rest.strip_prefix(WORKING_TREE) {
            Some(after) if after.first().is_none_or(|byte| *byte == b'/') => {
                rest = after;
                inside += 1;
            }
            _ => break,
        }
    }
    if inside == 0 {
        return path.to_vec();
    }
    let mut key = format!("\0{}", depth - inside).into_bytes();
    key.extend_from_slice(rest);
    key
}

/// The working directory in which a name is read (`WORKING_TREE`): deeper than the name's `..`
/// reach, so that no name climbs out of it.
fn working_directory(name: &[u8]) -> Vec<u8> {
    let climbs = name
        .split(|byte| *byte == b'/')
        .filter(|part| *part == b"..")
        .count();
    WORKING_TREE.repeat(climbs + 1)
}

/// The path of the file that the output name `name` opens (`file_key`), where `output_name` gives
/// a file; `None` where the system refuses the name.
fn output_key(name: &[u8]) -> Option<Vec<u8>> {
    let cwd = working_directory(name);
    Walk::resolve(name, &cwd).map(|path| file_key(&path, &cwd))
}

/// The path of the file that the `-logfile` name `name` opens (`file_key`), where `log_name` gives
/// a file: the path that NCBI opens (`log_file_name`); `None` where the system refuses it.
fn log_key(name: &[u8]) -> Option<Vec<u8>> {
    let (cwd, file) = log_file(name);
    Walk::resolve(&file.path, &cwd).map(|path| file_key(&path, &cwd))
}

/// The working directory of the `-logfile` name `name` (`working_directory`) and the paths that
/// NCBI reads from it there (`log_file_name`, Unix).
fn log_file(name: &[u8]) -> (Vec<u8>, crate::blastinput::ncbi_environment::LogFileName) {
    let cwd = working_directory(crate::ncbi_str::truncate_spaces(name));
    let file = crate::blastinput::ncbi_environment::log_file_name(
        crate::blastinput::ncbi_environment::PathRules::Unix,
        name,
        Some(&cwd),
    )
    .expect("NCBI's Unix reading reads every name");
    (cwd, file)
}

/// What an entry of the run's working directory is (`output_name`): LOSAT Web's host gives the
/// engine the files of the run's directory, whose folders are those that hold them (docs/web/
/// abi_v2.md §3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Entry {
    Directory,
    File,
    Missing,
}

/// The entry of the run's working directory at `path`, a path relative to it without links, `.`
/// and `..`: the file system's, whose current directory is the run's (LOSAT Web: the host's preopen
/// `/work`, `web/adapter/src/run.rs`).
fn run_directory_entry(path: &[u8]) -> Entry {
    let Ok(path) = std::str::from_utf8(path) else {
        return Entry::Missing;
    };
    match std::fs::metadata(path) {
        Ok(metadata) if metadata.is_dir() => Entry::Directory,
        Ok(_) => Entry::File,
        Err(_) => Entry::Missing,
    }
}

/// What an output file name gives for a `Receiver` (`output_name`).
#[derive(Debug, PartialEq, Eq)]
enum OutputName {
    /// The system does not open it: the command line's open fails.
    Unopened,
    /// A file of that name, made.
    File,
    /// The process's standard output.
    Stdout,
    /// Writes that are lost: the null device and the devices that take writes, and the standard
    /// input, the null device in LOSAT Web's run.
    Discard,
    /// Writes that fail: `/dev/full`.
    Full,
    /// The standard error, another device or a file of procfs, whose writes LOSAT Web does not
    /// reproduce (they are among the run's messages, or change the process).
    Special,
    /// A file of the run's working directory, which the command line empties.
    Existing,
}

/// The devices of `LOST_DEVICES` whose writes a user's system takes and loses (`output_name`):
/// the null device, `/dev/zero` and the random devices. `/dev/full` fails them; `/dev/ptmx`,
/// `/dev/fuse` and `/dev/net/tun` give writes to no reader or fail them each in its own way.
const DISCARDING_DEVICES: [&[u8]; 4] =
    [b"/dev/null", b"/dev/zero", b"/dev/random", b"/dev/urandom"];

/// What the output file name `name` gives on the command line on Linux (No NCBI counterpart: the
/// system's open of a stream of that name, `CNcbiOfstream`, which writes to the name as typed, with
/// no reading of NCBI's own; LOSAT Web has no file system; line NP WW3a, 2026-11-22), on the model
/// of `log_name` (the model's host and its user, `USER_DIRECTORIES`, the links of `LINKS`, at most
/// `MAX_LINKS`, the directories that the user cannot search and the files that the system opens to
/// a user, `opens`), in a working directory that is the run's: `entry` tells what the run's
/// directory holds at a path inside it. Each folder of the name is walked as the system walks it
/// (`Walk`): a folder inside the working directory must be one of the run's directory (`nodir/c`,
/// `in.fa/c` do not open), one outside it a directory that the model knows (`directory_at`: a
/// directory of root's that the model does not list, passed on the way back out, as `log_name`
/// takes it); a name of `.`, `..` or with a final `/`, and a directory that the walk passed, the
/// working directory or one of their ancestors, are directories, which do not open for writing;
/// a name in a directory of root's does not open (`opens`). The opened file is then the standard
/// output (`/dev/stdout`, `/dev/fd/1`, `/proc/self/fd/1`), the standard input (the null device of
/// the run: writes lost), the standard error, a file of procfs or a device of `LOST_DEVICES`
/// (`Special`, but the null device, `/dev/zero` and the random devices, whose writes are lost, and
/// `/dev/full`, whose writes fail), a queue of `/dev/mqueue` (`Special`), a file of the run's
/// directory (`Existing`; a folder of it does not open), or a file of that name in the working
/// directory, in a home or in a directory where every user makes files (`File`).
fn output_name(name: &[u8], entry: &dyn Fn(&[u8]) -> Entry) -> OutputName {
    if name.is_empty() || !fits(name) {
        return OutputName::Unopened;
    }
    let cwd = working_directory(name);
    // The path inside the run's directory of a path inside the working directory.
    let inside = |path: &[u8]| -> Option<Vec<u8>> {
        path.strip_prefix(cwd.as_slice())?
            .strip_prefix(b"/")
            .map(<[u8]>::to_vec)
    };
    let (folders, last) = match name.iter().rposition(|byte| *byte == b'/') {
        Some(at) => (&name[..at], &name[at + 1..]),
        None => (&name[..0], name),
    };
    let mut directories: Vec<Vec<u8>> = vec![cwd.clone()];
    let mut walk = Walk::at(name, &cwd);
    for part in folders.split(|byte| *byte == b'/') {
        walk.step(part);
        if walk.refused || directory_at(walk.path()) == Directory::None {
            return OutputName::Unopened;
        }
        if inside(walk.path()).is_some_and(|path| entry(&path) != Entry::Directory) {
            return OutputName::Unopened;
        }
        if directories.last().map(Vec::as_slice) != Some(walk.path()) {
            directories.push(walk.path().to_vec());
        }
    }
    if directory_at(walk.path()) != Directory::Known || matches!(last, b"" | b"." | b"..") {
        return OutputName::Unopened;
    }
    walk.step(last);
    if walk.refused {
        return OutputName::Unopened;
    }
    let path = walk.path().to_vec();
    if directories.iter().any(|directory| within(directory, &path)) {
        return OutputName::Unopened;
    }
    let stream = |descriptor: &[u8]| {
        [SELF, THREAD]
            .iter()
            .any(|of| path.strip_prefix(*of) == Some(descriptor))
    };
    if stream(b"/fd/1".as_slice()) {
        return OutputName::Stdout;
    }
    if stream(b"/fd/0".as_slice()) {
        return OutputName::Discard;
    }
    if within(&path, SELF) {
        return if is_procfs_file(&path) {
            OutputName::Special
        } else {
            OutputName::Unopened
        };
    }
    if path == b"/dev/full" {
        return OutputName::Full;
    }
    if DISCARDING_DEVICES.contains(&path.as_slice()) {
        return OutputName::Discard;
    }
    if LOST_DEVICES.contains(&path.as_slice()) {
        return OutputName::Special;
    }
    if !opens(&path) {
        return OutputName::Unopened;
    }
    if path.starts_with(b"/dev/mqueue/") {
        return OutputName::Special;
    }
    match inside(&path).map(|path| entry(&path)) {
        Some(Entry::Directory) => OutputName::Unopened,
        Some(Entry::File) => OutputName::Existing,
        Some(Entry::Missing) | None => OutputName::File,
    }
}

/// What a `-logfile` name that is not one of NCBI's special names gives, for a `Receiver`
/// (`RunFiles::log_target`).
#[derive(Debug, PartialEq, Eq)]
enum LogName {
    /// A name that NCBI cannot open, a directory: the messages stay on `cerr`.
    Unopened,
    /// The process's standard error, opened as a file.
    StandardError,
    /// No message: the null device written otherwise than NCBI's special `/dev/null`, or a file
    /// system that reports no free space.
    Nowhere,
    /// A file of that name.
    File,
}

/// The working directory's tree on LOSAT Web's Linux (`log_name`): the working directory is a
/// chain of `"\0"` components in it (`/\0/\0/...`), deeper than a name's `..` reach, so that no
/// name climbs out of it, as a user's command line runs in a directory of the user's own
/// (`~/work`). A command line's words have no NUL byte, so no name gives another path in it.
const WORKING_TREE: &[u8] = b"/\0";

/// The model's host, whose system LOSAT Web's `-logfile` names are read on (`log_name`), and the
/// reference of every measured list here: Ubuntu 24.04.3, kernel 5.15 (WSL2), and its first user
/// (uid 1000, in that user's groups), not root (line WQA's D11 round 7, 2026-11-21: one host and
/// one user, Owner-delegated).
///
/// The directories where the model's user makes files, besides the working directory's tree and a
/// home (`user_makes`): those that the user writes on the model's host, measured as the user
/// (`find -type d -writable`, line WQA's D11 round 7): where every user makes files (mode 1777,
/// `/var/crash` 3777: `/tmp`, `/var/tmp`, `/dev/shm` and `/run/shm`, `/run/lock`, also as
/// `/var/lock`, `/var/crash`, WSL's `/mnt/wsl` and `/mnt/wslg`), the user's own `/run/user/1000`,
/// and WSL's drives of Windows (`/mnt/c`, `/mnt/d`, mode 0777). `/dev/mqueue` (1777 too) holds
/// message queues only (`opens`).
const USER_DIRECTORIES: [&[u8]; 11] = [
    b"/tmp",
    b"/var/tmp",
    b"/dev/shm",
    b"/run/shm",
    b"/run/lock",
    b"/var/crash",
    b"/run/user/1000",
    b"/mnt/wsl",
    b"/mnt/wslg",
    b"/mnt/c",
    b"/mnt/d",
];

/// Whether the user of LOSAT Web's Linux makes entries in `directory`, a path without links: in
/// the working directory's tree (`WORKING_TREE`), in a user's home (`/home/<user>`, which LOSAT
/// Web takes as the user's own, not knowing the user's name) and in `USER_DIRECTORIES`, each with
/// the directories that the user makes in it. Every other directory is root's, as only root makes
/// entries in `/` (line WQA's D11 round 4, 2026-11-20, Owner-delegated: one rule of owners in
/// place of lists of root's directories): `/`, `/data`, `/home`, `/etc`, `/var/log`.
fn user_makes(directory: &[u8]) -> bool {
    within(directory, WORKING_TREE)
        || directory
            .strip_prefix(b"/home/")
            .is_some_and(|user| !user.is_empty())
        || USER_DIRECTORIES
            .iter()
            .any(|place| within(directory, place))
}

/// The directories at the root of LOSAT Web's Linux (`directory_at`): the model's host's (`ls /`,
/// line WQA's D11 round 7), an Ubuntu of merged `/usr` (`/bin`, `/lib`, `/lib32`, `/lib64` and
/// `/sbin` are links, `LINKS`) with its marks `*.usr-is-merged`, Snap's `/snap` and the file
/// system's `/lost+found`. WSL's `/init` is a file; another name at the root is not there
/// (`/data`).
const ROOT_DIRECTORIES: [&[u8]; 20] = [
    b"/bin.usr-is-merged",
    b"/boot",
    b"/dev",
    b"/etc",
    b"/home",
    b"/lib.usr-is-merged",
    b"/lost+found",
    b"/media",
    b"/mnt",
    b"/opt",
    b"/proc",
    b"/root",
    b"/run",
    b"/sbin.usr-is-merged",
    b"/snap",
    b"/srv",
    b"/sys",
    b"/tmp",
    b"/usr",
    b"/var",
];

/// LOSAT Web's process in procfs: `/proc/self` is a link to it (`LINKS`), as to a process id, here
/// a NUL byte, which no name gives (`WORKING_TREE`).
const SELF: &[u8] = b"/proc/\0";

/// LOSAT Web's thread in procfs: `/proc/thread-self` is a link to it (`LINKS`), a task of `SELF`
/// whose id is a NUL byte.
const THREAD: &[u8] = b"/proc/\0/task/\0";

/// The links of LOSAT Web's Linux and their targets (`Walk`), besides `cwd` of `SELF` and of
/// `THREAD`, the working directory: merged `/usr`'s, the Filesystem Hierarchy Standard's
/// `/var/lock` and `/var/run`, the standard streams' (`/dev/fd`, `/dev/stdin`, `/dev/stdout`,
/// `/dev/stderr`) and procfs's (`/proc/self`, `/proc/thread-self` and their `root`).
const LINKS: [(&[u8], &[u8]); 15] = [
    (b"/bin", b"/usr/bin"),
    (b"/lib", b"/usr/lib"),
    (b"/lib32", b"/usr/lib32"),
    (b"/lib64", b"/usr/lib64"),
    (b"/sbin", b"/usr/sbin"),
    (b"/var/lock", b"/run/lock"),
    (b"/var/run", b"/run"),
    (b"/dev/fd", b"/proc/self/fd"),
    (b"/dev/stdin", b"/proc/self/fd/0"),
    (b"/dev/stdout", b"/proc/self/fd/1"),
    (b"/dev/stderr", b"/proc/self/fd/2"),
    (b"/proc/self", SELF),
    (b"/proc/thread-self", THREAD),
    (b"/proc/\0/root", b"/"),
    (b"/proc/\0/task/\0/root", b"/"),
];

/// The links that the system follows in one lookup (`Walk`), each link of the path and of the
/// links' targets: a name through 40 writes its file on the model's host's command line, and one
/// through 41 does not (`ELOOP`; observed in line WQA's D11 review r6, 2026-11-21:
/// `/var/lock/../..` 40 and 41 times, `/proc/self/root`, two links, 20 and 21 times).
const MAX_LINKS: usize = 40;

/// The directories of LOSAT Web's Linux in which its user looks up no name (`searchable`): those
/// that the model's user cannot search on the model's host, measured as that user at every depth
/// but in the user's places and other processes (`find / -type d ! -executable`, line WQA's D11
/// round 7, 2026-11-21). Some are the host's state then: Snap's FUSE connections and the
/// directories of systemd's units in `/run/systemd/propagate`.
const CLOSED_DIRECTORIES: [&[u8]; 39] = [
    b"/etc/credstore",
    b"/etc/credstore.encrypted",
    b"/etc/polkit-1/rules.d",
    b"/etc/ssl/private",
    b"/lost+found",
    b"/proc/tty/driver",
    b"/root",
    b"/run/docker",
    b"/run/systemd/inaccessible/dir",
    b"/run/systemd/propagate/systemd-logind.service",
    b"/run/systemd/propagate/systemd-resolved.service",
    b"/run/systemd/propagate/systemd-timesyncd.service",
    b"/run/systemd/propagate/systemd-udevd.service",
    b"/sys/fs/fuse/connections/86",
    b"/sys/fs/fuse/connections/87",
    b"/sys/fs/fuse/connections/88",
    b"/sys/fs/fuse/connections/89",
    b"/sys/fs/fuse/connections/90",
    b"/sys/fs/fuse/connections/91",
    b"/sys/fs/fuse/connections/92",
    b"/sys/fs/fuse/connections/93",
    b"/sys/fs/fuse/connections/94",
    b"/sys/fs/fuse/connections/95",
    b"/sys/kernel/debug",
    b"/sys/kernel/tracing",
    b"/var/cache/apt/archives/partial",
    b"/var/cache/ldconfig",
    b"/var/cache/private",
    b"/var/lib/apt/lists/partial",
    b"/var/lib/docker",
    b"/var/lib/polkit-1",
    b"/var/lib/private",
    b"/var/lib/snapd/cache",
    b"/var/lib/snapd/cookie",
    b"/var/lib/ubuntu-advantage/apt-esm/var/lib/apt/lists/partial",
    b"/var/lib/update-notifier/package-data-downloads/partial",
    b"/var/log/private",
    b"/var/spool/cron/crontabs",
    b"/var/spool/rsyslog",
];

/// The entries of another process (`/proc/<pid>` and its threads' `/proc/<pid>/task/<tid>`) in
/// which the model's user looks up no name (`searchable`): its descriptors (`fd`, `map_files`) and
/// its links to its executable and its working and root directories (`exe`, `cwd`, `root`), which
/// only the process's user may follow. LOSAT Web's own process is `SELF`; a process of
/// `/proc/<pid>` is root's, as `/proc/1` (measured as the user on the model's host: `stat` of
/// `<entry>/.` is refused there, and allowed in `fdinfo`, `ns`, `attr` and `net`; line WQA's D11
/// round 7, 2026-11-21).
const CLOSED_PROCESS_ENTRIES: [&[u8]; 5] = [b"cwd", b"exe", b"fd", b"map_files", b"root"];

/// Whether a user looks up names in `directory`, a path without links (`Walk`): not in
/// `CLOSED_DIRECTORIES` nor in another process's `CLOSED_PROCESS_ENTRIES`, where the system
/// refuses (`EACCES`), whatever the name (`.` and `..` too).
fn searchable(directory: &[u8]) -> bool {
    let number = |part: &[u8]| !part.is_empty() && part.iter().all(u8::is_ascii_digit);
    let closed_entry = directory.strip_prefix(b"/proc/").is_some_and(|entry| {
        match entry
            .split(|byte| *byte == b'/')
            .collect::<Vec<_>>()
            .as_slice()
        {
            [pid, closed] => number(pid) && CLOSED_PROCESS_ENTRIES.contains(closed),
            [pid, task, tid, closed] => {
                number(pid)
                    && *task == b"task"
                    && number(tid)
                    && CLOSED_PROCESS_ENTRIES.contains(closed)
            }
            _ => false,
        }
    });
    !closed_entry && !CLOSED_DIRECTORIES.contains(&directory)
}

/// The file systems of Linux that report no free space (`statfs`'s `f_bavail` is 0): procfs,
/// sysfs, devpts and mqueue, where they are mounted.
const NO_FREE_SPACE: [&[u8]; 4] = [b"/proc", b"/sys", b"/dev/pts", b"/dev/mqueue"];

/// The devices of LOSAT Web's Linux that a user opens for writing (mode 0666), whose writes are
/// lost: the null device written otherwise than NCBI's special `/dev/null` (`//dev/null`),
/// `/dev/full` (whose writes fail, which NCBI does not check), `/dev/zero`, `/dev/random`,
/// `/dev/urandom`, and a new terminal's master (`/dev/ptmx`), FUSE's device (`/dev/fuse`) and the
/// TUN device (`/dev/net/tun`), whose writes fail or no one reads (line WQA's D11 round 3,
/// 2026-11-20: NCBI BLAST+ 2.17.0 and LOSAT's command line write no message there). `/dev/tty`
/// is a device too, which a run without a terminal does not open.
const LOST_DEVICES: [&[u8]; 8] = [
    b"/dev/null",
    b"/dev/full",
    b"/dev/zero",
    b"/dev/random",
    b"/dev/urandom",
    b"/dev/ptmx",
    b"/dev/fuse",
    b"/dev/net/tun",
];

/// The files of procfs that a process may open for writing as a user, of the process (`SELF`) and
/// of its thread (`THREAD`): those of the model's host's Linux 5.15 that `open` opens with
/// `O_WRONLY`, `O_APPEND` and `O_CREAT` (line WQA's sweep, D11 round 2, 2026-11-20). The file
/// descriptors 0 to 2 (`fd/0`...) are the standard streams.
const PROCESS_FILES: [&[u8]; 18] = [
    b"attr/current",
    b"attr/exec",
    b"attr/fscreate",
    b"attr/keycreate",
    b"attr/sockcreate",
    b"clear_refs",
    b"comm",
    b"coredump_filter",
    b"gid_map",
    b"loginuid",
    b"mem",
    b"oom_adj",
    b"oom_score_adj",
    b"projid_map",
    b"sched",
    b"timens_offsets",
    b"timerslack_ns",
    b"uid_map",
];
const THREAD_FILES: [&[u8]; 15] = [
    b"attr/current",
    b"attr/exec",
    b"attr/fscreate",
    b"attr/keycreate",
    b"attr/sockcreate",
    b"clear_refs",
    b"comm",
    b"gid_map",
    b"loginuid",
    b"mem",
    b"oom_adj",
    b"oom_score_adj",
    b"projid_map",
    b"sched",
    b"uid_map",
];

/// Linux's limits of the paths given to the system (Linux 6.8: include/uapi/linux/limits.h:12-13,
/// `NAME_MAX` and `PATH_MAX`): a component of at most 255 bytes, a path of fewer than 4096 (with
/// its closing NUL). The system refuses a longer one (`ENAMETOOLONG`), which NCBI then cannot open.
const NAME_MAX: usize = 255;
const PATH_MAX: usize = 4096;

/// Whether the system takes `path`, a path given to it as it is (`NAME_MAX`, `PATH_MAX`).
fn fits(path: &[u8]) -> bool {
    path.len() < PATH_MAX
        && path
            .split(|byte| *byte == b'/')
            .all(|part| part.len() <= NAME_MAX)
}

/// Whether `path` is `directory` or a directory in it, component by component.
fn within(path: &[u8], directory: &[u8]) -> bool {
    path.strip_prefix(directory)
        .is_some_and(|rest| rest.is_empty() || rest[0] == b'/' || directory.ends_with(b"/"))
}

/// Whether `path`, a path without links, is a file of procfs that the model knows in `SELF` or
/// `THREAD`: a standard stream (`fd/0` to `fd/2`) or a file of `PROCESS_FILES` or `THREAD_FILES`.
fn is_procfs_file(path: &[u8]) -> bool {
    let stream = |file: &[u8]| matches!(file, b"fd/0" | b"fd/1" | b"fd/2");
    let file = |directory: &[u8]| path.strip_prefix(directory)?.strip_prefix(b"/");
    match (file(THREAD), file(SELF)) {
        (Some(file), _) => stream(file) || THREAD_FILES.contains(&file),
        (None, Some(file)) => stream(file) || PROCESS_FILES.contains(&file),
        (None, None) => false,
    }
}

/// Whether `path`, a path without links, is one of the files that the model knows (`opens`): a
/// device of `LOST_DEVICES` or a file of procfs (`is_procfs_file`). Such a file is no directory
/// (`/dev/null/x`, `/dev/stdout/../stderr`, `/dev/fd/2/..`).
fn is_known_file(path: &[u8]) -> bool {
    LOST_DEVICES.contains(&path) || is_procfs_file(path)
}

/// Whether `directory` holds a file that the model knows (`is_known_file`): `/dev`, `/dev/net`,
/// `SELF` and its `fd` and `attr`, its `task`, `THREAD` and its `fd` and `attr`.
fn holds_known_file(directory: &[u8]) -> bool {
    let holds = |path: &[u8]| path.len() > directory.len() && within(path, directory);
    let in_procfs = |process: &[u8], files: &[&[u8]]| {
        [b"fd/0".as_slice()]
            .iter()
            .chain(files)
            .any(|file| holds(&[process, b"/".as_slice(), *file].concat()))
    };
    LOST_DEVICES.iter().any(|device| holds(device))
        || in_procfs(SELF, PROCESS_FILES.as_slice())
        || in_procfs(THREAD, THREAD_FILES.as_slice())
}

/// What LOSAT Web's Linux has at `directory`, a path without links (`log_name`).
#[derive(Debug, PartialEq, Eq)]
enum Directory {
    /// A directory that the model knows: `/`, one where the user makes entries (`user_makes`),
    /// which the model takes as made where a name gives it (`logs`, `/tmp/a`; the working
    /// directory and its ancestors), a directory at the root (`ROOT_DIRECTORIES`), a mount of
    /// `NO_FREE_SPACE` (`/dev/pts`) and a directory of a file that the model knows (`/dev/net`,
    /// `SELF`'s `fd`).
    Known,
    /// Any other path in root's directories, which may be a directory there or not
    /// (`/etc/default`, `/usr/tmp`): the user cannot change them, and the model does not list
    /// them.
    Unknown,
    /// No directory: a file that the model knows (`/dev/null`, `SELF`'s `fd/1`) or a name at the
    /// root that is not there (`/data`).
    None,
}

/// What the model knows of `directory`, a path without links (`Directory`).
fn directory_at(directory: &[u8]) -> Directory {
    if directory == b"/"
        || user_makes(directory)
        || ROOT_DIRECTORIES.contains(&directory)
        || NO_FREE_SPACE.contains(&directory)
        || holds_known_file(directory)
    {
        Directory::Known
    } else if is_known_file(directory)
        || !ROOT_DIRECTORIES.iter().any(|root| within(directory, root))
    {
        Directory::None
    } else {
        Directory::Unknown
    }
}

/// Whether a user's system opens `path` for writing, a path without links in a directory that
/// exists: a file that a name may make or write, also a device (No NCBI counterpart). In procfs,
/// only the standard streams and the files of `PROCESS_FILES` and `THREAD_FILES` (a name that
/// procfs does not have is not made there); sysfs and devpts open none of a user's names; an
/// mqueue file is a name in `/dev/mqueue`; the devices of `LOST_DEVICES` are open to every user.
/// Any other name opens where the user makes entries (`user_makes`): a new name in a directory of
/// root's is not made, and another device or a file of the system's (`/dev/kmsg`, `/dev/console`,
/// `/dev/tty` in a run without a terminal, `/etc/passwd`) is not open to a user.
fn opens(path: &[u8]) -> bool {
    if within(path, SELF) {
        return is_procfs_file(path);
    }
    if let Some(queue) = path.strip_prefix(b"/dev/mqueue/") {
        return !queue.contains(&b'/');
    }
    if [b"/proc".as_slice(), b"/sys", b"/dev/pts"]
        .iter()
        .any(|system| within(path, system))
    {
        return false;
    }
    LOST_DEVICES.contains(&path)
        || user_makes(&path[..path.iter().rposition(|byte| *byte == b'/').unwrap_or(0)])
}

/// A path of LOSAT Web's Linux walked component by component as the system walks it (No NCBI
/// counterpart: the system's resolution, `log_name`): `.` and empty components stay, `..` goes to
/// the parent of the path walked so far (`/` stays), and a link (`LINKS`; `cwd` of `SELF` and
/// `THREAD`, the working directory `cwd`) is replaced by its target, which is walked from `/`. The
/// system refuses the path where a component is looked up in a directory that the user cannot
/// search (`searchable`: `/root/..`, `/root/.`), and at the link past `MAX_LINKS` (also a
/// standard stream of procfs, `fd/0` to `fd/2`, which is a link there). A name is walked once,
/// each component once and each link's target once more; a step's work is at most the length of
/// the path walked so far, under `PATH_MAX` (line WQA's D11 round 5, 2026-11-21: resolving each
/// prefix again was cubic, a minute for a name of 4 KB in LOSAT Web's serial engine).
struct Walk<'a> {
    /// The path walked so far, without links; empty for `/`.
    resolved: Vec<u8>,
    /// The working directory, the target of `cwd` in `SELF`.
    cwd: &'a [u8],
    /// The links that the walk followed.
    links: usize,
    /// Whether the system refused the path.
    refused: bool,
}

impl<'a> Walk<'a> {
    /// The walk of `path` before its first component: at `/` for an absolute path, else at the
    /// working directory `cwd`.
    fn at(path: &[u8], cwd: &'a [u8]) -> Self {
        let mut walk = Walk {
            resolved: Vec::new(),
            cwd,
            links: 0,
            refused: false,
        };
        if path.first() != Some(&b'/') {
            walk.walk(cwd);
        }
        walk
    }

    /// `path` with its links followed and its `.` and `..` taken; `None` where the system refuses
    /// it (`refused`).
    fn resolve(path: &[u8], cwd: &'a [u8]) -> Option<Vec<u8>> {
        let mut walk = Walk::at(path, cwd);
        walk.walk(path);
        (!walk.refused).then(|| walk.path().to_vec())
    }

    /// Takes every component of `path`.
    fn walk(&mut self, path: &[u8]) {
        for part in path.split(|byte| *byte == b'/') {
            self.step(part);
        }
    }

    /// Takes the component `part`.
    fn step(&mut self, part: &[u8]) {
        #[cfg(test)]
        WALK_STEPS.with(|steps| steps.set(steps.get() + 1));
        if part.is_empty() || self.refused {
            return;
        }
        if !searchable(self.path()) {
            self.refused = true;
            return;
        }
        match part {
            b"." => {}
            b".." => {
                let parent = self.resolved.iter().rposition(|byte| *byte == b'/');
                self.resolved.truncate(parent.unwrap_or(0));
            }
            _ => {
                self.resolved.push(b'/');
                self.resolved.extend_from_slice(part);
                let path = self.resolved.as_slice();
                let link = |of: &[u8], name: &[u8]| path.strip_prefix(of) == Some(name);
                let target = if link(SELF, b"/cwd".as_slice()) || link(THREAD, b"/cwd".as_slice()) {
                    Some(self.cwd)
                } else {
                    LINKS
                        .iter()
                        .find(|(link, _)| *link == path)
                        .map(|(_, target)| *target)
                };
                let stream = [SELF, THREAD].iter().any(|of| {
                    [b"/fd/0".as_slice(), b"/fd/1", b"/fd/2"]
                        .iter()
                        .any(|fd| link(of, fd))
                });
                if target.is_some() || stream {
                    self.links += 1;
                    if self.links > MAX_LINKS {
                        self.refused = true;
                        return;
                    }
                }
                if let Some(target) = target {
                    self.resolved.clear();
                    self.walk(target);
                }
            }
        }
    }

    /// The path walked so far (`/` at the root).
    fn path(&self) -> &[u8] {
        if self.resolved.is_empty() {
            b"/"
        } else {
            &self.resolved
        }
    }
}

// The components that `Walk::step` takes on this thread (the test of `log_name`'s cost).
#[cfg(test)]
thread_local! {
    static WALK_STEPS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// What the `-logfile` name `name` gives on the command line on Linux (No NCBI counterpart: LOSAT
/// Web has no file system; line WQA's sweep, D11, 2026-11-20 and its rounds 2 to 7), on the
/// model's host as its first user (`USER_DIRECTORIES`). The name is
/// read as NCBI reads it on Unix (`ncbi_environment::log_file_name`, the command line's reading
/// too) in a working directory of the user's own (`WORKING_TREE`) on a file system with free
/// space: `logs/run.log`, ` run.log` (`run.log`), `logs/` (the file `logs`) and `/tmp/run.log`
/// are files. A directory of '/'s only (`//run.log`) is none, as `CDirEntry::Reset` leaves it
/// empty. The system then follows the links of Linux (`LINKS`, at most `MAX_LINKS`;
/// `/proc/self/cwd` is the working directory) in the directory that NCBI checks, in the path that
/// it opens and in the directory whose free space it reads (`Walk`, once over the name), and
/// refuses a component of more than
/// 255 bytes and a path of 4096 or more (`fits`: `aaa...`, a directory `x/../x/../...` of 4096
/// bytes). The directory that NCBI checks (the name's as written) is walked as the system walks
/// it (`directory_at`): each directory that it passes may be one of root's that the model does
/// not know, which a name passes on its way back out (`/etc/default/../../tmp`, as on the command
/// line), but none that is not there (`/data`, a file: `/dev/null/x`, `/dev/stdout/..`) nor one
/// that the user cannot search (`searchable`: `/root/..`, `/var/cache/ldconfig/..`); the
/// directory where it ends must be one that the model knows, root's or one where the user makes
/// entries, which the user made (`logs`, `/tmp/a`, `/home/me/logs`); not `/bin/..` (`/usr`, then
/// `/usr/tmp`) or `/var/lock/..` (`/run`, then `/run/tmp`). The opened path is a directory, which
/// NCBI cannot open (`cerr`), where it is the working directory or one of its ancestors, or a
/// directory that the walk passed or one of its ancestors (`.`, `logs/.`, `a/b/..`, `/tmp/`,
/// `/proc/self/cwd`); a directory of root's opens to no user (`opens`: `/tmp`, `/proc/self/fd`);
/// NCBI cannot open a path that the system does not open to a user either (`opens`: `/dev/fd/3`,
/// `/proc/self/status`, `/sys/...`, `/dev/pts/0`, `/dev/tty`, a name in a directory of root's:
/// `/run.log`, `/etc/run.log`, `/var/log/run.log`, `/dev/run.log`, `/dev/kmsg`). A file opened in
/// a file system that reports no free space (`NO_FREE_SPACE`) gets no message, as NCBI checks
/// the space before it writes (`/dev/fd/2`, `/proc/self/fd/1`, `/proc/self/comm`,
/// `/dev/mqueue/q`, `/dev/fd/../stderr`, whose directory is `/proc/self`). Otherwise
/// `/dev/stderr` is the standard error, also written otherwise (`//dev/stderr`); the standard
/// input is the null device (`/dev/stdin`), as for a command line run without one; the devices
/// whose writes are lost (`LOST_DEVICES`: `/dev/ptmx`) give no message; `/dev/stdout` is a file
/// of that name, whose bytes the command line writes into its standard output among the
/// report's, `/proc/self/cwd/run.log` the file `run.log` of the working directory, and
/// `/var/lock/run.log` a file in `/run/lock`, each under the name as typed. The model's command
/// line is a user's (not root), with no terminal, its standard output and error being pipes, as a
/// shell's `|` gives.
fn log_name(name: &[u8]) -> LogName {
    let (cwd, file) = log_file(name);
    // The directories that the walk of the checked directory passes, each a directory of the
    // system: the opened path is none of them nor an ancestor of one.
    let mut directories: Vec<Vec<u8>> = Vec::new();
    if let Some(directory) = &file.directory {
        if directory.is_empty() || !fits(directory) {
            return LogName::Unopened;
        }
        let mut walk = Walk::at(directory, &cwd);
        for part in directory.split(|byte| *byte == b'/') {
            walk.step(part);
            if walk.refused || directory_at(walk.path()) == Directory::None {
                return LogName::Unopened;
            }
            if directories.last().map(Vec::as_slice) != Some(walk.path()) {
                directories.push(walk.path().to_vec());
            }
        }
        if directory_at(walk.path()) != Directory::Known {
            return LogName::Unopened;
        }
    }
    directories.push(cwd.clone());
    let Some(path) = Walk::resolve(&file.path, &cwd) else {
        return LogName::Unopened;
    };
    if file.path.is_empty()
        || !fits(&file.path)
        || directories.iter().any(|directory| within(directory, &path))
        || !opens(&path)
    {
        return LogName::Unopened;
    }
    // A directory whose free space the system does not read keeps the file (NCBI ignores the
    // error); its walk is a part of the checked directory's, which the system took.
    let space = Walk::resolve(&file.space, &cwd).unwrap_or_default();
    if NO_FREE_SPACE.iter().any(|system| within(&space, system)) {
        return LogName::Nowhere;
    }
    let stream = |descriptor: &[u8]| path.strip_prefix(SELF) == Some(descriptor);
    if stream(b"/fd/2".as_slice()) {
        LogName::StandardError
    } else if stream(b"/fd/0".as_slice()) || LOST_DEVICES.contains(&path.as_slice()) {
        LogName::Nowhere
    } else {
        LogName::File
    }
}

#[cfg(test)]
#[path = "run_files_tests.rs"]
mod tests;

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

    /// One file of a search's run besides its outputs, whole (`receive`), with the option that
    /// names it (`RunFileOption`), which the run knows where it hands the file over: the file of
    /// `-logfile` (`RunFiles::deliver_log`), the strategy file (`strategy::end_run`) and the PSSM
    /// files (`psiblast::pssm`, `psiblast::rounds`). LOSAT Web's host labels the file by it
    /// (docs/web/abi_v2.md §5, stream 8), as two options' names may name one file where the run
    /// opens only one of them (line RG item F1w-r2, review 3 S1). Default: `receive`, for a
    /// receiver that keeps no option.
    fn receive_run_file(&self, option: RunFileOption, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        let _ = option;
        self.receive(name, bytes)
    }

    /// Opens the file `name` of a tool's run (`OpenFiles`: `-out FILE`, `-logfile FILE`) to take
    /// the run's writes as the run makes them (`FileStream`); the file is complete when the stream
    /// is dropped, when the run ends. `None`, the default: the receiver takes each file whole, which
    /// the run keeps in memory and hands to `receive` when it ends (the files of a search, whose
    /// bytes and order this keeps). LOSAT Web's run of a tool streams its files, which may be larger
    /// than the browser's memory (line NP WW3a, D-WW3-6, 2026-11-22).
    fn open(&self, _name: &[u8]) -> Option<Box<dyn FileStream>> {
        None
    }

    /// The names of the run's files that the host holds outside the run's working directory
    /// (LOSAT Web's registered inputs of a search: `-query`, `-subject`, the registry files that
    /// the run's start reads, and the counts file of `-window_masker_db`). Default: none (a
    /// tool's run opens its files in the working directory, which `run_directory_entry` reads).
    fn registered_files(&self) -> &[Vec<u8>] {
        &[]
    }

    /// The names of files that the run opens for reading (whether present in `/work` or missing),
    /// besides `registered_files` (line WQA D19 round R2): e.g. a `-db` search's list files
    /// (`-gilist`, `-ipglist`, `-seqidlist`, `-taxidlist` and their negatives) and a tool's
    /// input files (`-in`, `-ustat`, `-ids`, `-exclude_ids`). Kept apart from `registered_files`
    /// so a missing input does not make its parent folders present; only matched when `log_name`
    /// would otherwise make a new file (`LogName::File { key }`).
    fn read_files(&self) -> &[Vec<u8>] {
        &[]
    }

    /// The subset of `read_files` that a tool opens and reads after opening `-out` in a way that
    /// changes the run when the file is missing (line WQA D19 round R2): `-in` of `windowmasker
    /// -mk_counts` (without `-fa_list` or `-infmt blastdb`) and `-ustat` of masking. Only matched
    /// when `output_name` would otherwise make a new file (`OutputName::File { key }`).
    fn output_read_files(&self) -> &[Vec<u8>] {
        &[]
    }

    /// Whether the run writes the files of its working directory through the file system, which
    /// is then the command line's directory of the run's files (LOSAT Web's `makeblastdb`, whose
    /// working directory the host makes writable, docs/web/abi_v2.md §3; line NP WW4b1 R2,
    /// 2026-11-28): a `-logfile` name of a file in that directory (`log_name`, `output_name`:
    /// a new file there or one of the run's files) opens there, through the file system, as the
    /// command line opens it in its directory, and the command line's refusal of a `-logfile`
    /// naming an input (`report::log_inputs::check_disk_logfile`) applies, as on `Disk`; only the
    /// names outside it and the devices take the model's reading. Default: `false`, the run reads
    /// the files of its working directory and writes none of them (LOSAT's refusal of a name of
    /// one of them, `LogName::Existing`, `OutputName::Existing`).
    fn writes_working_directory(&self) -> bool {
        false
    }
}

/// The option of a search's argv that names a file of its run besides its outputs
/// (`FileReceiver::receive_run_file`). No NCBI counterpart: NCBI writes files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunFileOption {
    Logfile,
    ExportSearchStrategy,
    OutPssm,
    OutAsciiPssm,
}

impl RunFileOption {
    /// The option's name without its hyphen, as LOSAT Web's host names it (`logfile`,
    /// `export_search_strategy`, `out_pssm`, `out_ascii_pssm`).
    pub fn name(self) -> &'static str {
        match self {
            RunFileOption::Logfile => "logfile",
            RunFileOption::ExportSearchStrategy => "export_search_strategy",
            RunFileOption::OutPssm => "out_pssm",
            RunFileOption::OutAsciiPssm => "out_ascii_pssm",
        }
    }

    /// The option as the argv writes it (`-out_pssm`).
    pub fn flag(self) -> &'static str {
        match self {
            RunFileOption::Logfile => "-logfile",
            RunFileOption::ExportSearchStrategy => "-export_search_strategy",
            RunFileOption::OutPssm => "-out_pssm",
            RunFileOption::OutAsciiPssm => "-out_ascii_pssm",
        }
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
#[derive(Clone, Default)]
pub enum RunFiles {
    /// The file system, at the names that the arguments give (the command line, the library).
    #[default]
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
    /// The target of the `-logfile` value `name` (`None` without one), and the file's host name
    /// when a `Receiver` keeps it in memory (`deliver_log`). On `Disk` the file system's
    /// (`diag_target`, `None`). For a `Receiver`: NCBI's special names as there
    /// (`special_log_target`), and any other name as the command line of a user on Linux reads it
    /// in a working directory holding exactly the run's files (`log_name`, `run_entries`): `cerr`
    /// for a folder of the run's directory or a name whose folder inside the working directory is
    /// missing or is a file, and for a name that a user's command line cannot open (`//run.log`,
    /// `/etc/run.log`, `/data/run.log`, `/dev/kmsg`, a component of 256 bytes); `Err` with LOSAT's
    /// refusal reason when the name is a file of the run's directory (`LogName::Existing`); the
    /// file that is the standard error (`DiagTarget::StderrFile`) for `/dev/stderr`; no output for
    /// another spelling of `/dev/null`, for a device whose writes are lost (`/dev/ptmx`) and for a
    /// name in a file system that reports no free space (`/dev/fd/2`, `/dev/fd/1`,
    /// `/proc/self/fd/1`); a file in memory named `/dev/stdout` for `/dev/stdout` (a search's
    /// reports are split by format); and for any other file name in the working directory, in a
    /// home or in a directory where every user makes files, a file in memory under the command
    /// line's file name (`host_file_name`: `run.log` for ` run.log` and `./run.log`, `logs` for
    /// `logs/`, `/tmp/run.log` outside the working directory). A tool's run opens its files
    /// through `OpenFiles`, which reads the names the same way (`received_log_target`).
    pub fn log_target(
        &self,
        name: Option<&OsStr>,
    ) -> Result<(DiagTarget, Option<Vec<u8>>), String> {
        match self {
            RunFiles::Disk => Ok((diag_target(name), None)),
            RunFiles::Receiver(receiver) => received_log_target(
                name,
                receiver.as_ref(),
                |target| (target, None),
                |host_name, _key| (DiagTarget::Memory(Default::default()), Some(host_name)),
                || {
                    (
                        DiagTarget::Memory(Default::default()),
                        Some(b"/dev/stdout".to_vec()),
                    )
                },
                // A file of a working directory that the run writes (`writes_working_directory`):
                // in memory too, so that a check of the words (`losat_web2_validate` of
                // makeblastdb) writes nothing; the run opens it through `OpenFiles::log_target`.
                |_relative| (DiagTarget::Memory(Default::default()), None),
            ),
        }
    }

    /// Like [`RunFiles::log_target`], also enforcing the command-line / library self-reference
    /// refusal on [`RunFiles::Disk`] before opening or creating `-logfile` (line WQA D23,
    /// `Q-D23=A`): refuses any `-logfile` pointing to an input file read by `start.program()`
    /// with `words` and `start` (`report::log_inputs::check_disk_logfile`).
    pub(crate) fn log_target_for<S: AsRef<OsStr>>(
        &self,
        words: &[S],
        start: &crate::blastinput::ncbi_environment::ApplicationStart,
        name: Option<&OsStr>,
    ) -> Result<(DiagTarget, Option<Vec<u8>>), String> {
        if self.writes_on_disk() {
            crate::report::log_inputs::check_disk_logfile(start.program(), words, start, name)?;
        }
        self.log_target(name)
    }

    /// Whether the run writes its files through the file system: `Disk`, or a `Receiver` whose run
    /// writes its working directory (`FileReceiver::writes_working_directory`), where the command
    /// line's refusal of a `-logfile` that names an input applies as on `Disk`.
    fn writes_on_disk(&self) -> bool {
        match self {
            RunFiles::Disk => true,
            RunFiles::Receiver(receiver) => receiver.writes_working_directory(),
        }
    }

    /// The host file name of the `-logfile` value `name` (what a `Receiver` names the stream 8
    /// entry if the file is made: `run.log` for `./run.log` and ` run.log`, `logs` for `logs/`,
    /// `/tmp/run.log` outside the working directory, `/dev/stdout` for a search's `/dev/stdout`),
    /// without the run's files (`losat_web2_start`, which has no run directory yet): read in the
    /// run's files that make it a file if any do (`made_file`, as `output_host_name`), so that
    /// every file that a run of the name makes has this name (`/tmp/run.log` for
    /// `/data/../tmp/run.log`, which a run makes beside a registered `/data/q.fa`). `None` for
    /// NCBI's special names (`-`, `""`, `/dev/null`) and names that never make a stream 8 file
    /// (`.`, `..`, `/dev/stderr`, `/dev/fd/2`, `/etc/run.log`).
    pub fn log_host_name(name: &OsStr) -> Option<Vec<u8>> {
        if special_log_target(Some(name)).is_some() {
            return None;
        }
        let name = name.as_encoded_bytes();
        let (cwd, file) = log_file(name);
        let made = made_file(Walk::resolve(&file.path, &cwd));
        match log_name(name, Some(&made)) {
            LogName::File { name, .. } => Some(name),
            LogName::Stdout => Some(b"/dev/stdout".to_vec()),
            LogName::Unopened | LogName::StandardError | LogName::Nowhere | LogName::Existing => {
                None
            }
        }
    }

    /// Hands the file of `-logfile`, named `name` (`log_target`), which `target` holds, to the
    /// receiver. Nothing on `Disk` (the file is written) and for a target without a file.
    pub fn deliver_log(&self, name: Option<&[u8]>, target: DiagTarget) -> io::Result<()> {
        match (self, name, target) {
            (RunFiles::Receiver(receiver), Some(name), DiagTarget::Memory(bytes)) => receiver
                .receive_run_file(
                    RunFileOption::Logfile,
                    name,
                    &bytes.into_inner().unwrap_or_else(PoisonError::into_inner),
                ),
            _ => Ok(()),
        }
    }

    /// What a psiblast PSSM output file (`-out_pssm` or `-out_ascii_pssm`) named `name` resolves
    /// to for a `Receiver` (`PssmTarget`): uses the same working-directory output-name rules as
    /// `OpenFiles::create` (`receiver_output_name`). No NCBI counterpart: NCBI writes files.
    pub(crate) fn pssm_target(&self, option: &str, name: &OsStr) -> PssmTarget {
        let RunFiles::Receiver(receiver) = self else {
            return PssmTarget::Unopened;
        };
        match receiver_output_name(receiver.as_ref(), name.as_encoded_bytes()) {
            OutputName::File {
                name: host_name, ..
            } => PssmTarget::File(host_name),
            OutputName::Discard => PssmTarget::Discard,
            OutputName::Full => PssmTarget::Full,
            OutputName::Unopened => PssmTarget::Unopened,
            OutputName::Stdout => PssmTarget::Refused(format!(
                "an {option} file that is also the standard output ({}) is not supported by LOSAT's PSIBLAST",
                name.display()
            )),
            OutputName::StandardError | OutputName::Special => PssmTarget::Refused(format!(
                "an {option} file '{}' that is the standard error, a device or a file of \
                 procfs (LOSAT Web writes no such file)",
                name.display()
            )),
            OutputName::Existing => PssmTarget::Refused(format!(
                "an {option} file '{}' that is a file of the run's directory (the command \
                 line replaces it; LOSAT Web reads the run's files and writes none of them)",
                name.display()
            )),
        }
    }

    /// The host file name of an output file name `name` that a search writes besides its outputs
    /// (`-out_pssm`, `-out_ascii_pssm`, `-export_search_strategy`): what a `Receiver` names the
    /// stream 8 entry if the run makes the file (`output_name`'s `File`: `s.asn` for `./s.asn`,
    /// `in/../s.asn` and `/proc/self/cwd/s.asn`, `/tmp/s.asn` outside the working directory),
    /// without the run's files (`losat_web2_start`, which has no run directory yet). The name is
    /// read in the run's files that make it a file if any do (`made_file`), so that every file
    /// that a run of the name makes has this name (`/tmp/s.asn` for `/data/../tmp/s.asn`, which a
    /// run makes beside a registered `/data/q.fa`). `None` for a name that makes no file in any
    /// run's files (`/dev/null`, `/dev/stdout`, `/etc/s.asn`, `.`, `..`, `s.asn/`).
    pub fn output_host_name(name: &OsStr) -> Option<Vec<u8>> {
        let name = name.as_encoded_bytes();
        if name.is_empty() || !fits(name) {
            return None;
        }
        let cwd = working_directory(name);
        let made = made_file(Some(Walk::resolve(name, &cwd)?));
        match output_name(name, &made) {
            OutputName::File { name, .. } => Some(name),
            OutputName::Unopened
            | OutputName::Stdout
            | OutputName::StandardError
            | OutputName::Discard
            | OutputName::Full
            | OutputName::Special
            | OutputName::Existing => None,
        }
    }
}

/// What a psiblast PSSM output file name (`-out_pssm` or `-out_ascii_pssm`) gives for a
/// `Receiver` (`RunFiles::pssm_target`). No NCBI counterpart: NCBI writes files.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum PssmTarget {
    /// A file of the run under `host_name` (`host_file_name`), which the receiver gets.
    File(Vec<u8>),
    /// No file: writes are lost (the null device and discarding devices).
    Discard,
    /// Writes fail (`/dev/full`).
    Full,
    /// Cannot be opened by the command line (a directory, a missing folder inside the working
    /// directory, or a path a user cannot open).
    Unopened,
    /// The standard output, standard error, a special device/procfs file, or an existing file in
    /// the run's directory (refused by LOSAT Web).
    Refused(String),
}

/// LOSAT's refusal reason for a `-logfile` naming an existing file of the run's directory or a
/// file the run opens for reading.
fn existing_log_refusal(name: &OsStr) -> String {
    format!(
        "a -logfile file '{}' that is a file of the run's directory (LOSAT Web reads the run's \
         files and writes none of them)",
        name.display()
    )
}

/// LOSAT's refusal reason for a tool's `-out` naming an existing file of the run's directory or
/// an input the command line reads after opening `-out`.
fn existing_output_refusal(name: &OsStr) -> String {
    format!(
        "an output file '{}' that is a file of the run's directory (the command \
         line empties it; LOSAT Web reads the run's files and writes none of them)",
        name.display()
    )
}

/// LOSAT's refusal reason for a tool's `-out` naming the standard error, a device or a file of
/// procfs (`OutputName::StandardError`, `OutputName::Special`), whose writes LOSAT Web does not
/// reproduce.
fn special_output_refusal(name: &OsStr) -> String {
    format!(
        "an output file '{}' that is the standard error, a device or a file of \
         procfs (LOSAT Web writes no such file)",
        name.display()
    )
}

/// The target of the `-logfile` value `name` for a `Receiver` (`RunFiles::log_target`,
/// `OpenFiles::log_target`): NCBI's special names, then `log_name`'s reading against the run's
/// real files (`RunView`), plus the case-2 refusal when a new file's `file_key` matches a name in
/// `FileReceiver::read_files`; `direct` wraps a `DiagTarget` without a new file,
/// `file` opens a new file given its host name and key (`LogName::File`), `stdout` handles
/// `/dev/stdout` (`LogName::Stdout`), and `here` a file of the run's working directory (a new one
/// or one of the run's files) where the run writes that directory
/// (`FileReceiver::writes_working_directory`), given its path relative to the directory
/// (`run_directory_path`). `Err` with LOSAT's refusal reason for `LogName::Existing` elsewhere.
fn received_log_target<T>(
    name: Option<&OsStr>,
    receiver: &dyn FileReceiver,
    direct: impl FnOnce(DiagTarget) -> T,
    file: impl FnOnce(Vec<u8>, Vec<u8>) -> T,
    stdout: impl FnOnce() -> T,
    here: impl FnOnce(Vec<u8>) -> T,
) -> Result<T, String> {
    match (special_log_target(name), name) {
        (Some(target), _) => Ok(direct(target)),
        (None, None) => Ok(direct(DiagTarget::Stderr)),
        (None, Some(name)) => {
            let view = RunView::new(receiver, receiver.read_files());
            let bytes = name.as_encoded_bytes();
            // The file's path in the run's working directory, where the run writes it.
            let written_here = |key: &[u8]| {
                receiver
                    .writes_working_directory()
                    .then(|| run_directory_path(key))
                    .flatten()
            };
            match view.log_name(bytes) {
                LogName::Unopened => Ok(direct(DiagTarget::Stderr)),
                LogName::StandardError => Ok(direct(DiagTarget::StderrFile)),
                LogName::Nowhere => Ok(direct(DiagTarget::Discard)),
                LogName::Stdout => Ok(stdout()),
                LogName::Existing => match log_key(bytes).as_deref().and_then(written_here) {
                    Some(path) => Ok(here(path)),
                    None => Err(existing_log_refusal(name)),
                },
                LogName::File {
                    name: host_name,
                    key,
                } => {
                    if let Some(path) = written_here(&key) {
                        Ok(here(path))
                    } else if view.reads(&key) {
                        Err(existing_log_refusal(name))
                    } else {
                        Ok(file(host_name, key))
                    }
                }
            }
        }
    }
}

/// The `-logfile` target opened for a tool run (`OpenFiles::log_target`).
pub(crate) enum ToolLog {
    /// A diagnostic target (`DiagTarget`).
    Target(DiagTarget),
    /// The run's standard output (`LogName::Stdout` on a `Receiver`), shared with the tool's
    /// `out` writer (`ncbi_app::run_app`).
    Stdout,
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
///
/// The tools' output streams that the run opened (`keep_stream`: NCBI's `-out` and `cout`,
/// `ncbi_app::files::tool_stream`) are ended after the run and its texts (`finish_streams`), as
/// NCBI's arguments close their streams after `Run`.
#[derive(Clone)]
pub(crate) struct OpenFiles {
    table: Option<Arc<FileTable>>,
    streams: Arc<Mutex<Vec<crate::report::ostream::ReportStream>>>,
    /// The run is the command line's process (`ncbi_app::run`, `of_command_line`): its output
    /// streams are NCBI's and tied to its standard error (`ncbi_app::files::tool_stream`).
    command_line: bool,
    /// The run's `-logfile` shares standard output on a `Receiver` (`ToolLog::Stdout`):
    /// `OutputFile::Caller` forwards NCBI's cuts and flushes to the shared buffered stdout stream
    /// so interleaving matches the CLI's `Cout` stream.
    tool_log_stdout: bool,
}

impl std::fmt::Debug for OpenFiles {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(if self.table.is_some() {
            "Receiver"
        } else {
            "Disk"
        })
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
        OpenFiles {
            table: match files {
                RunFiles::Disk => None,
                RunFiles::Receiver(receiver) => Some(Arc::new(FileTable {
                    receiver: receiver.clone(),
                    files: Mutex::new(Vec::new()),
                })),
            },
            streams: Arc::default(),
            command_line: false,
            tool_log_stdout: false,
        }
    }

    /// The run's files as the command line's process's (`ncbi_app::run`): its standard output
    /// and its `-out` are NCBI's streams, tied to its standard error (`command_line`).
    pub(crate) fn of_command_line(mut self) -> Self {
        self.command_line = true;
        self
    }

    /// Marks the run's files as sharing standard output with `-logfile /dev/stdout`
    /// (`ToolLog::Stdout` in `ncbi_app::run_app`).
    pub(crate) fn with_tool_log_stdout(mut self) -> Self {
        self.tool_log_stdout = true;
        self
    }

    /// Whether the run is the command line's process (`of_command_line`); a library run or a run
    /// of LOSAT Web writes to its caller's writers.
    pub(crate) fn command_line(&self) -> bool {
        self.command_line
    }

    /// Whether `-logfile` writes to the run's shared standard output (`with_tool_log_stdout`).
    pub(crate) fn tool_log_stdout(&self) -> bool {
        self.tool_log_stdout
    }

    /// The files of a run on the file system (`RunFiles::Disk`).
    pub(crate) fn disk() -> Self {
        OpenFiles {
            table: None,
            streams: Arc::default(),
            command_line: false,
            tool_log_stdout: false,
        }
    }

    /// Keeps `stream`, a tool's output stream of the run, to end it after the run
    /// (`finish_streams`).
    pub(crate) fn keep_stream(&self, stream: crate::report::ostream::ReportStream) {
        self.streams
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .push(stream);
    }

    /// Ends the tools' output streams of the run, in the order it opened them
    /// (`ReportStream::finish_tool`): whether a write met a closed pipe, which ended NCBI's
    /// process (`report::ostream::ends_like_sigpipe`).
    pub(crate) fn finish_streams(&self) -> bool {
        let streams = std::mem::take(
            &mut *self
                .streams
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner),
        );
        let mut closed = false;
        for stream in streams {
            closed |= stream.finish_tool().is_err() || stream.closed();
        }
        closed
    }

    /// The target of the `-logfile` value `name` (`RunFiles::log_target`), whose file, for a
    /// `Receiver`, is an opening of the run's file of that path that appends to it
    /// (`DiagTarget::Stream`), or `ToolLog::Stdout` for `/dev/stdout`. `Err` with LOSAT's refusal
    /// reason when the name is a file of the run's directory.
    pub(crate) fn log_target(&self, name: Option<&OsStr>) -> Result<ToolLog, String> {
        match &self.table {
            None => Ok(ToolLog::Target(diag_target(name))),
            Some(table) => received_log_target(
                name,
                table.receiver.as_ref(),
                ToolLog::Target,
                |host_name, key| {
                    ToolLog::Target(DiagTarget::Stream(Mutex::new(Box::new(table.open(
                        &host_name,
                        key,
                        Opening::Append,
                    )))))
                },
                || ToolLog::Stdout,
                // A file of the working directory that the run writes: the file system's, as
                // `diag_target` opens it there (`cerr` where it does not open).
                |path| {
                    ToolLog::Target(
                        append_in_run_directory(&path).map_or(DiagTarget::Stderr, DiagTarget::File),
                    )
                },
            ),
        }
    }

    /// Like [`OpenFiles::log_target`], also enforcing the `Disk` self-reference refusal before
    /// opening or creating `-logfile` (`report::log_inputs::check_disk_logfile`, line WQA D23).
    pub(crate) fn log_target_for<S: AsRef<OsStr>>(
        &self,
        words: &[S],
        start: &crate::blastinput::ncbi_environment::ApplicationStart,
        name: Option<&OsStr>,
    ) -> Result<ToolLog, String> {
        if self
            .table
            .as_ref()
            .is_none_or(|table| table.receiver.writes_working_directory())
        {
            crate::report::log_inputs::check_disk_logfile(start.program(), words, start, name)?;
        }
        self.log_target(name)
    }

    /// The output file `name` that a tool opens by name for writing (NCBI's `CNcbiOfstream` of the
    /// name: `-out FILE` of windowmasker, which made or emptied the file; the callers take `-` and
    /// the empty name first). On `Disk` the file of the file system, made or emptied
    /// (`File::create`, as before); `Err` where it does not open. For a `Receiver`, the name as the
    /// command line of a user on Linux opens it in the run's working directory (`output_name`): the
    /// run's file of that path, emptied (a file that the run opened before under that path or
    /// another name of it, `file_key`, is that file, named by `host_file_name`); the run's standard
    /// output for a name of it (`/dev/stdout`); writes that are lost (the null device) or that
    /// fail (`/dev/full`); `Err` where that command line cannot open it (a directory, a folder
    /// that the run's directory does not have, a directory of root's); and LOSAT's refusal of the
    /// names whose writes LOSAT Web does not reproduce (`Created::Refused`).
    pub(crate) fn create(&self, name: &(impl AsRef<OsStr> + ?Sized)) -> Result<Created, Unopened> {
        let name = name.as_ref();
        let Some(table) = &self.table else {
            return File::create(name)
                .map(|file| Created::File(OpenedFile::Disk(file)))
                .map_err(|_| Unopened);
        };
        let bytes = name.as_encoded_bytes();
        match receiver_output_name(table.receiver.as_ref(), bytes) {
            OutputName::Unopened => Err(Unopened),
            OutputName::File {
                name: host_name,
                key,
            } => Ok(Created::File(OpenedFile::Shared(table.open(
                &host_name,
                key,
                Opening::Create,
            )))),
            OutputName::Stdout => Ok(Created::Stdout),
            OutputName::Discard => Ok(Created::File(OpenedFile::Discard)),
            OutputName::Full => Ok(Created::File(OpenedFile::Full)),
            OutputName::StandardError | OutputName::Special => {
                Ok(Created::Refused(special_output_refusal(name)))
            }
            OutputName::Existing => Ok(Created::Refused(existing_output_refusal(name))),
        }
    }

    /// LOSAT Web's refusal reason for a tool's `-out` `name` (without opening or creating the
    /// file, for the checks of `validate` where they know that the run reaches `-out`:
    /// `windowmasker::app::check`, `seqmasks_io::run::check_output`): `Some(reason)` when `self`
    /// is a `Receiver` and `create` would refuse `name` (`receiver_output_name` resolves it to
    /// `OutputName::Existing`, the case-2 refusal, or to `OutputName::StandardError` or
    /// `OutputName::Special`, the standard error, a device or a file of procfs; line NP WW4c,
    /// 2026-11-28), with `create`'s reason; else `None`. Always `None` on `Disk`.
    pub(crate) fn refused_output(&self, name: &OsStr) -> Option<String> {
        let table = self.table.as_ref()?;
        match receiver_output_name(table.receiver.as_ref(), name.as_encoded_bytes()) {
            OutputName::Existing => Some(existing_output_refusal(name)),
            OutputName::StandardError | OutputName::Special => Some(special_output_refusal(name)),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1199-1201 (cited at
    // `writedb::makeblastdb::run_with`)
    /// The file of the `-logfile` value `name` that a tool opens as an output to append its log to,
    /// besides its diagnostics' opening (`log_target`): makeblastdb's `AsOutputFile` with
    /// `fAppend`, which the system opens under the name as written (line WQA's D11 round 2). On
    /// `Disk` the file of the file system, made where it is missing (`Err` where it does not open,
    /// as for the empty name and `logs/`). For a `Receiver`, the name as the command line of a user
    /// on Linux opens it in the run's working directory (`output_name`): a file of a working
    /// directory that the run writes (`FileReceiver::writes_working_directory`) through the file
    /// system, appended to, as there; another file, an opening that appends to the run's file of
    /// that path (`file_key`), which the diagnostics' opening of the same path shares
    /// (`-logfile /tmp/mk.log`); the run's standard output (`/dev/stdout`) and standard error
    /// (`/dev/stderr`), where the log's lines go among the run's other bytes as the system writes
    /// them; lost writes (the null device) and failing ones (`/dev/full`), which makeblastdb does
    /// not check; `Err` where that command line cannot open it; and LOSAT's refusal of another
    /// device, a file of procfs and a file of the run's directory where the run writes none
    /// (`log_output_refusal`).
    pub(crate) fn log_output(&self, name: &OsStr) -> Result<LogOutput, Unopened> {
        let bytes = name.as_encoded_bytes();
        let Some(table) = &self.table else {
            if bytes.is_empty() {
                return Err(Unopened);
            }
            return std::fs::OpenOptions::new()
                .append(true)
                .create(true)
                .open(name)
                .map(|file| LogOutput::File(OpenedFile::Disk(file)))
                .map_err(|_| Unopened);
        };
        let receiver = table.receiver.as_ref();
        let here = |key: &[u8]| {
            receiver
                .writes_working_directory()
                .then(|| run_directory_path(key))
                .flatten()
        };
        let in_run_directory = |path: Vec<u8>| {
            append_in_run_directory(&path)
                .map(|file| LogOutput::File(OpenedFile::Disk(file)))
                .ok_or(Unopened)
        };
        let named = receiver_output_name(receiver, bytes);
        if let Some(reason) = log_output_refusal(name, &named) {
            return Ok(LogOutput::Refused(reason));
        }
        match named {
            OutputName::Unopened => Err(Unopened),
            OutputName::File {
                name: host_name,
                key,
            } => match here(&key) {
                Some(path) => in_run_directory(path),
                None => Ok(LogOutput::File(OpenedFile::Shared(table.open(
                    &host_name,
                    key,
                    Opening::Append,
                )))),
            },
            OutputName::Existing => match output_key(bytes).as_deref().and_then(here) {
                Some(path) => in_run_directory(path),
                None => Ok(LogOutput::Refused(existing_log_refusal(name))),
            },
            OutputName::Stdout => Ok(LogOutput::Stdout),
            OutputName::StandardError => Ok(LogOutput::Stderr),
            OutputName::Discard => Ok(LogOutput::File(OpenedFile::Discard)),
            OutputName::Full => Ok(LogOutput::File(OpenedFile::Full)),
            OutputName::Special => unreachable!("refused above (`log_output_refusal`)"),
        }
    }
}

impl RunFiles {
    /// LOSAT's refusal reason for the opening of the `-logfile` value `name` as an output that a
    /// tool appends its log to (`OpenFiles::log_output`), without opening anything (the checks of
    /// a tool's words, which refuse what its run refuses, before the run opens a file): `Some` for
    /// a `Receiver` where the name is a device that LOSAT Web does not write, a file of procfs
    /// (`OutputName::Special`) or, where the run writes none of its directory's files, one of
    /// them. Always `None` on `Disk`.
    pub(crate) fn refused_log_output(&self, name: &OsStr) -> Option<String> {
        let RunFiles::Receiver(receiver) = self else {
            return None;
        };
        let named = receiver_output_name(receiver.as_ref(), name.as_encoded_bytes());
        log_output_refusal(name, &named).or_else(|| {
            (named == OutputName::Existing && !receiver.writes_working_directory())
                .then(|| existing_log_refusal(name))
        })
    }
}

/// LOSAT's refusal reason of the `-logfile` value `name` as an output that a tool appends its log
/// to (`OpenFiles::log_output`), which `named` reads: a device whose writes LOSAT Web does not
/// reproduce or a file of procfs, whose writes change the process (`OutputName::Special`; the
/// standard error is the run's, `OutputName::StandardError`).
fn log_output_refusal(name: &OsStr, named: &OutputName) -> Option<String> {
    (*named == OutputName::Special).then(|| {
        format!(
            "a -logfile file '{}' that is a device or a file of procfs (LOSAT Web writes no \
             such file)",
            name.display()
        )
    })
}

/// What `OpenFiles::log_output` opened: the file of `-logfile` that a tool appends its log to.
pub(crate) enum LogOutput {
    /// A file, appended to (`OpenedFile::Disk` on the file system).
    File(OpenedFile),
    /// The run's standard output (`/dev/stdout` of a `Receiver`).
    Stdout,
    /// The run's standard error (`/dev/stderr` of a `Receiver`).
    Stderr,
    /// LOSAT's refusal of the name, the reason.
    Refused(String),
}

/// Resolves `bytes` as an output file name of `receiver` (`output_name` in the run's view,
/// `RunView`), turning a new file whose `file_key` matches `FileReceiver::output_read_files` into
/// `OutputName::Existing`.
fn receiver_output_name(receiver: &dyn FileReceiver, bytes: &[u8]) -> OutputName {
    let view = RunView::new(receiver, receiver.output_read_files());
    match view.output_name(bytes) {
        OutputName::File { key, .. } if view.reads(&key) => OutputName::Existing,
        other => other,
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

    fn write_vectored(&mut self, bufs: &[io::IoSlice<'_>]) -> io::Result<usize> {
        match self {
            OpenedFile::Disk(file) => file.write_vectored(bufs),
            OpenedFile::Shared(file) => file.write_vectored(bufs),
            OpenedFile::Discard => Ok(bufs.iter().map(|b| b.len()).sum()),
            OpenedFile::Full => Err(io::Error::from(io::ErrorKind::StorageFull)),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OpenedFile::Disk(file) => file.flush(),
            _ => Ok(()),
        }
    }
}

/// Where `path` (without links) sits relative to the working directory `cwd`: `Some((up, rest))`
/// when `path` is in `WORKING_TREE` (`up` directories above `cwd`, then `rest` without a leading
/// `/`), else `None`.
fn in_working_tree<'a>(path: &'a [u8], cwd: &[u8]) -> Option<(usize, &'a [u8])> {
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
    (inside > 0).then(|| (depth - inside, rest.strip_prefix(b"/").unwrap_or(rest)))
}

/// The path of a file of a run's working directory `cwd` (`working_directory`), without links,
/// as one name of the run's files (`OpenFiles`): every name that the system opens at that path
/// gives it (`x`, `./x`, `sub/../x`, `/proc/self/cwd/x`). A path in the working directory's tree
/// (`WORKING_TREE`) is told by its place from the working directory (how many directories above
/// it, then the rest), as each name has a working directory of its own depth; any other path is
/// itself.
fn file_key(path: &[u8], cwd: &[u8]) -> Vec<u8> {
    match in_working_tree(path, cwd) {
        Some((up, rest)) => {
            let mut key = format!("\0{up}").into_bytes();
            if !rest.is_empty() {
                key.push(b'/');
                key.extend_from_slice(rest);
            }
            key
        }
        None => path.to_vec(),
    }
}

/// The file's name that the host gets (`log_name`, `output_name`): for a file in the working
/// directory's tree (`in_working_tree`), its path relative to the run's directory (`run.log`,
/// `logs/run.log`, `logs`, `../x.log`); for a file outside it, `opened`: for `-logfile` the path
/// that NCBI opens, which NCBI normalizes lexically (`/tmp/run.log`, `/var/lock/run.log`,
/// `/var/run/../tmp/run.log` as `/var/tmp/run.log`), and for a file that the system opens as
/// typed, the name read lexically where it walks to the file, else the file's path
/// (`output_name`: `/var/lock/s.asn`, `/var/run/../tmp/s.asn` as `/tmp/s.asn`).
fn host_file_name(path: &[u8], opened: &[u8], cwd: &[u8]) -> Vec<u8> {
    match in_working_tree(path, cwd) {
        Some((up, rest)) => {
            let mut relative = b"../".repeat(up);
            relative.extend_from_slice(rest);
            relative
        }
        None => opened.to_vec(),
    }
}

/// The run's files in which a name whose opened path is `file` makes it a new file, if any run's
/// files do (`RunFiles::output_host_name`, `RunFiles::log_host_name`): every other path that the
/// name looks up is a folder of them, and `file` is none (a registered file makes its folders
/// present outside the working directory too, `run_entry`, as `/data` of `/data/q.fa`). Such a
/// lookup only makes folders present, and a folder that the walk passes is never the opened path.
fn made_file(file: Option<Vec<u8>>) -> impl Fn(&[u8], &[u8]) -> Option<Entry> {
    move |path, _| (Some(path) != file.as_deref()).then_some(Entry::Directory)
}

/// The relative path of `path` strictly inside the working directory `cwd`.
fn inside_cwd<'a>(path: &'a [u8], cwd: &[u8]) -> Option<&'a [u8]> {
    path.strip_prefix(cwd)?.strip_prefix(b"/")
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
/// a file; `None` where the system refuses the name. A name that a run reads (`-query`, compared
/// with the strategy file's in `RunFiles::strategy_target`) is at the same path.
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

/// The path relative to the run's working directory of the file whose path is `key` (`file_key`),
/// where that file is in the working directory (`file_key`'s `up` 0: `run.log`, `logs/run.log`);
/// `None` for a file elsewhere (`../x.log`, `/tmp/run.log`) and for the directory itself.
fn run_directory_path(key: &[u8]) -> Option<Vec<u8>> {
    key.strip_prefix(b"\x000/")
        .filter(|rest| !rest.is_empty())
        .map(<[u8]>::to_vec)
}

/// Opens the file `path` of the run's working directory (`run_directory_path`) through the file
/// system, which is the run's working directory, to append to it, made where it is missing, as
/// NCBI opens the file of `-logfile` (`ncbi_environment::diag_target`) and makeblastdb its log
/// (`AsOutputFile` with `fAppend`); `None` where the system does not open it.
fn append_in_run_directory(path: &[u8]) -> Option<File> {
    std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(crate::utils::os_bytes::path_of_bytes(path))
        .ok()
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

/// What an entry of the run's working directory is (`log_name`, `output_name`): LOSAT Web's host
/// gives the engine the files of the run's directory, whose folders are those that hold them
/// (docs/web/abi_v2.md §3), and for a search also the host's registered input names
/// (`FileReceiver::registered_files`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Entry {
    Directory,
    File,
    Missing,
}

/// Whether `path` (resolved without links) is a directory that the model knows can never be a
/// regular file of a run (`run_entries`): `cwd` or its ancestors (`in_working_tree` with empty
/// `rest`), `/`, or a system directory (`ROOT_DIRECTORIES`, `USER_DIRECTORIES`, `NO_FREE_SPACE`,
/// `holds_known_file`).
fn is_known_directory(path: &[u8], cwd: &[u8]) -> bool {
    path == b"/"
        || in_working_tree(path, cwd).is_some_and(|(_, rest)| rest.is_empty())
        || ROOT_DIRECTORIES.contains(&path)
        || USER_DIRECTORIES.contains(&path)
        || NO_FREE_SPACE.contains(&path)
        || holds_known_file(path)
}

/// The keys (`file_key`) of `raw_names`, inside or outside the working directory, resolved with
/// `Walk` in each name's own `working_directory` so that `./q.fa`, `in/../q.fa`, and
/// `/proc/self/cwd/q.fa` share one key (`\00/q.fa`), `../q.fa` carries its climb depth
/// (`\01/q.fa`), and `/tmp/q.fa` or `/data/q.fa` keeps its resolved path. Non-file directory
/// names (`.`, `..`, `/`, `/tmp`, trailing `/`) are skipped.
fn path_keys(raw_names: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let mut files = Vec::new();
    for raw in raw_names {
        if raw.is_empty() || raw.ends_with(b"/") || !fits(raw) {
            continue;
        }
        let cwd = working_directory(raw);
        if let Some(path) = Walk::resolve(raw, &cwd) {
            if !is_known_directory(&path, &cwd) {
                files.push(file_key(&path, &cwd));
            }
        }
    }
    files
}

/// The keys (`file_key`) of `receiver`'s registered files (`path_keys`, `run_entry`).
fn run_entries(receiver: &dyn FileReceiver) -> Vec<Vec<u8>> {
    path_keys(receiver.registered_files())
}

/// The one view of a `Receiver`'s run files through which the run reads a name that it writes
/// (line WQA D19; line RG item F1w-r2): `-logfile` (`received_log_target`), a tool's output
/// (`receiver_output_name`) and the search strategy file (`RunFiles::strategy_target`). It holds
/// the keys of the run's registered files (`run_entries`), looked up with the run's directory
/// (`run_entry`), and the names of the files that the run reads, which a new file of the name may
/// not be (`read`: `FileReceiver::read_files`, or `FileReceiver::output_read_files` for a tool's
/// output; keyed by `path_keys` where a name makes a new file).
#[derive(Default)]
struct RunView<'a> {
    registered: Vec<Vec<u8>>,
    read: &'a [Vec<u8>],
}

impl<'a> RunView<'a> {
    /// The view of `receiver`'s files, `read` being the names of the files that its run reads.
    fn new(receiver: &dyn FileReceiver, read: &'a [Vec<u8>]) -> Self {
        RunView {
            registered: run_entries(receiver),
            read,
        }
    }

    /// What the run's files hold at `path` in the working directory `cwd` (`run_entry`).
    fn entry(&self, path: &[u8], cwd: &[u8]) -> Option<Entry> {
        run_entry(&self.registered, path, cwd)
    }

    /// The output file name `name` in the run's files (`output_name`).
    fn output_name(&self, name: &[u8]) -> OutputName {
        output_name(name, &|path: &[u8], cwd: &[u8]| self.entry(path, cwd))
    }

    /// The `-logfile` name `name` in the run's files (`log_name`).
    fn log_name(&self, name: &[u8]) -> LogName {
        log_name(name, Some(&|path: &[u8], cwd: &[u8]| self.entry(path, cwd)))
    }

    /// Whether the run reads the file of the path `key` (`file_key`), present or missing.
    fn reads(&self, key: &[u8]) -> bool {
        path_keys(self.read)
            .iter()
            .any(|read| read.as_slice() == key)
    }
}

/// What the run's files hold at `path` (resolved without links) in working directory `cwd`:
/// `Some(Entry::File)` when a registered file's `file_key` equals `key`,
/// `Some(Entry::Directory)` when a registered file sits below `key`, else
/// `Some(run_directory_entry(rel))` when `path` is strictly inside `cwd`, else `None` (outside
/// `cwd` and not named by a registered file, where the caller keeps the model's `directory_at`
/// and `opens` rules).
fn run_entry(registered: &[Vec<u8>], path: &[u8], cwd: &[u8]) -> Option<Entry> {
    let key = file_key(path, cwd);
    if registered
        .iter()
        .any(|file| file.as_slice() == key.as_slice())
    {
        return Some(Entry::File);
    }
    if registered
        .iter()
        .any(|file| file.len() > key.len() && within(file, &key))
    {
        return Some(Entry::Directory);
    }
    inside_cwd(path, cwd).map(run_directory_entry)
}

/// Whether `path` (resolved without links) lives in a virtual kernel file system or device tree
/// rather than a regular file directory: procfs, sysfs, or `/dev` outside `/dev/shm`.
fn is_virtual_or_device(path: &[u8]) -> bool {
    within(path, b"/proc")
        || within(path, b"/sys")
        || (within(path, b"/dev") && !within(path, b"/dev/shm"))
}

/// Lookup of an entry of the run's files (`run_entry`) at a resolved path in `cwd`.
type EntryLookup<'a> = &'a dyn Fn(&[u8], &[u8]) -> Option<Entry>;

/// Whether `path` (resolved without links) in working directory `cwd` is a usable folder
/// (`walks_to_directory`, `log_name`, `output_name`): an entry of the run's files must be
/// `Entry::Directory` (which also makes a registered file's ancestor outside `cwd` present even
/// where `directory_at` returns `Directory::None`, such as `/data` for `/data/q.fa`), else the
/// model's `directory_at` must be `Directory::Known` when `final_dir` is true or not
/// `Directory::None` for an intermediate step.
fn folder_present(
    path: &[u8],
    cwd: &[u8],
    entry: Option<EntryLookup<'_>>,
    final_dir: bool,
) -> bool {
    match entry.and_then(|entry| entry(path, cwd)) {
        Some(Entry::Directory) => true,
        Some(Entry::File | Entry::Missing) => false,
        None if final_dir => directory_at(path) == Directory::Known,
        None => directory_at(path) != Directory::None,
    }
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
    /// A file made at this host name (`host_file_name`) and file-table key (`file_key`).
    File { name: Vec<u8>, key: Vec<u8> },
    /// The process's standard output.
    Stdout,
    /// The process's standard error (`/dev/stderr`, `/dev/fd/2`, `/proc/self/fd/2`), whose writes
    /// are among the run's messages: LOSAT Web reproduces only a log written there a line at a
    /// time (`OpenFiles::log_output`) and refuses an output (`OpenFiles::create`).
    StandardError,
    /// Writes that are lost: the null device and the devices that take writes, and the standard
    /// input, the null device in LOSAT Web's run.
    Discard,
    /// Writes that fail: `/dev/full`.
    Full,
    /// Another device or a file of procfs, whose writes LOSAT Web does not reproduce (they change
    /// the process, or no one reads them).
    Special,
    /// A file of the run's working directory, which the command line empties.
    Existing,
}

/// The devices of `LOST_DEVICES` whose writes a user's system takes and loses (`output_name`,
/// `RunFiles::strategy_target`): the null device, `/dev/zero` and the random devices. `/dev/full`
/// fails them; `/dev/ptmx`, `/dev/fuse` and `/dev/net/tun` give writes to no reader or fail them
/// each in its own way.
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
/// the run: writes lost), the standard error (`StandardError`: `/dev/stderr`, `/dev/fd/2`,
/// `/proc/self/fd/2`), a file of procfs or a device of `LOST_DEVICES`
/// (`Special`, but the null device, `/dev/zero` and the random devices, whose writes are lost, and
/// `/dev/full`, whose writes fail), a queue of `/dev/mqueue` (`Special`), a file of the run's
/// directory (`Existing`; a folder of it does not open), or a file of that name in the working
/// directory, in a home or in a directory where every user makes files (`File`), named as the
/// command line's file (`host_file_name`: outside the working directory, the name read lexically
/// where it walks to the file, else the file's path, as a `..` after a link steps back from the
/// link's target: `/var/run/../tmp/s.asn` is `/tmp/s.asn`).
fn output_name(name: &[u8], entry: EntryLookup<'_>) -> OutputName {
    if name.is_empty() || !fits(name) {
        return OutputName::Unopened;
    }
    let cwd = working_directory(name);
    let (folders, last) = match name.iter().rposition(|byte| *byte == b'/') {
        Some(at) => (&name[..at], &name[at + 1..]),
        None => (&name[..0], name),
    };
    let mut directories: Vec<Vec<u8>> = vec![cwd.clone()];
    let mut walk = Walk::at(name, &cwd);
    for part in folders.split(|byte| *byte == b'/') {
        walk.step(part);
        if walk.refused || !folder_present(walk.path(), &cwd, Some(entry), false) {
            return OutputName::Unopened;
        }
        if directories.last().map(Vec::as_slice) != Some(walk.path()) {
            directories.push(walk.path().to_vec());
        }
    }
    if !folder_present(walk.path(), &cwd, Some(entry), true) || matches!(last, b"" | b"." | b"..") {
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
    if stream(b"/fd/2".as_slice()) {
        return OutputName::StandardError;
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
    let found = entry(&path, &cwd);
    if found == Some(Entry::Directory) {
        return OutputName::Unopened;
    }
    let existing = found == Some(Entry::File) && !is_virtual_or_device(&path);
    if !existing && !opens(&path) {
        return OutputName::Unopened;
    }
    if path.starts_with(b"/dev/mqueue/") {
        return OutputName::Special;
    }
    if existing {
        return OutputName::Existing;
    }
    // The name as typed, read lexically, names the file where it walks to it (`/var/lock/s.asn`);
    // where a `..` follows a link, the system steps back from the link's target, so the file is
    // its path (`/var/run/../tmp/s.asn` makes `/tmp/s.asn`, not `/var/tmp/s.asn`).
    let typed = crate::blastinput::ncbi_environment::PathRules::Unix
        .normalize(name, None)
        .expect("a path normalized without links has no link limit");
    let opened = match Walk::resolve(&typed, &cwd) {
        Some(file) if file == path => typed,
        _ => path.clone(),
    };
    OutputName::File {
        name: host_file_name(&path, &opened, &cwd),
        key: file_key(&path, &cwd),
    }
}

/// What a `-logfile` name that is not one of NCBI's special names gives, for a `Receiver`
/// (`RunFiles::log_target`, `OpenFiles::log_target`).
#[derive(Debug, PartialEq, Eq)]
enum LogName {
    /// A name that NCBI cannot open (a directory, a missing folder or a file as a folder inside
    /// the working directory, or a path a user cannot open): the messages stay on `cerr`.
    Unopened,
    /// The process's standard error, opened as a file.
    StandardError,
    /// The process's standard output (`/dev/stdout`, where the free-space check does not give
    /// `Nowhere`): in a tool's run, stream 11; in a search's run, a stream 8 file named
    /// `/dev/stdout`.
    Stdout,
    /// No message: the null device written otherwise than NCBI's special `/dev/null`, or a file
    /// system that reports no free space.
    Nowhere,
    /// A file of the run's working directory (refused by LOSAT Web).
    Existing,
    /// A new file made at this host name (`host_file_name`) and file-table key (`file_key`).
    File { name: Vec<u8>, key: Vec<u8> },
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
    /// A directory that the model knows: `/`, one where the user makes entries (`user_makes`:
    /// outside the working directory the model takes a user directory as present, such as
    /// `/tmp/a`; inside the working directory `entry` tells the real folders of the run), a
    /// directory at the root (`ROOT_DIRECTORIES`), a mount of `NO_FREE_SPACE` (`/dev/pts`) and a
    /// directory of a file that the model knows (`/dev/net`, `SELF`'s `fd`).
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
/// Web has no file system; line WQA's sweep, D11 and D19), on the model's host as its first user
/// (`USER_DIRECTORIES`), in a working directory holding the entries that `entry` reports
/// (`run_entries`, `run_directory_entry`). The name is read as NCBI reads it on Unix
/// (`ncbi_environment::log_file_name`, the command line's reading too) in a working directory of
/// the user's own (`WORKING_TREE`) on a file system with free space: `run.log`, ` run.log`
/// (`run.log`), `logs/` (the file `logs` when `logs` is not a folder or file of the run),
/// `logs/run.log` when `logs` is a folder of the run, and `/tmp/run.log` are files. A directory of
/// '/'s only (`//run.log`) is none, as `CDirEntry::Reset` leaves it empty. The system then follows
/// the links of Linux (`LINKS`, at most `MAX_LINKS`; `/proc/self/cwd` is the working directory) in
/// the directory that NCBI checks (`LogFileName::directory`), in the path that it opens
/// (`LogFileName::path`) and in the directory whose free space it reads (`Walk`, once over each),
/// and refuses a component of more than 255 bytes and a path of 4096 or more (`fits`). Every
/// folder looked up inside the working directory must be `Entry::Directory`, else
/// `LogName::Unopened`; outside it, folders keep the model's `directory_at` rules. The opened path
/// cannot be the working directory, an ancestor, a directory passed by the walk, or a path that a
/// user cannot open (`opens`). A file opened in a file system that reports no free space
/// (`NO_FREE_SPACE`) gets no message (`/dev/fd/2`, `/dev/fd/1`, `/proc/self/fd/1`,
/// `/proc/self/fd/2`, `/proc/self/comm`, `/dev/mqueue/q`, `/dev/fd/../stderr`). Otherwise
/// `/dev/stderr` (also `//dev/stderr`) is `StandardError`; `/dev/stdout` (also `//dev/stdout`) is
/// `Stdout`; `/dev/stdin` and `LOST_DEVICES` are `Nowhere`; inside the working directory,
/// `Entry::Directory` is `Unopened`, `Entry::File` is `Existing` (refused by LOSAT Web), and
/// `Entry::Missing` is `File` named relative to the working directory (`host_file_name`); outside
/// it, `File` is named by the path that NCBI opens (`LogFileName::path`).
fn log_name(name: &[u8], entry: Option<EntryLookup<'_>>) -> LogName {
    let (cwd, file) = log_file(name);
    // The directories that the walk of the checked directory passes, each a directory of the
    // system: the opened path is none of them nor an ancestor of one.
    let mut directories: Vec<Vec<u8>> = Vec::new();
    if let Some(directory) = &file.directory {
        if !walks_to_directory(directory, &cwd, entry, &mut directories) {
            return LogName::Unopened;
        }
    }
    directories.push(cwd.clone());
    if file.path.is_empty() || !fits(&file.path) {
        return LogName::Unopened;
    }
    let (folders, last) = match file.path.iter().rposition(|byte| *byte == b'/') {
        Some(at) => (&file.path[..at], &file.path[at + 1..]),
        None => (&file.path[..0], file.path.as_slice()),
    };
    let mut walk = Walk::at(&file.path, &cwd);
    for part in folders.split(|byte| *byte == b'/') {
        walk.step(part);
        if walk.refused || !folder_present(walk.path(), &cwd, entry, false) {
            return LogName::Unopened;
        }
    }
    walk.step(last);
    if walk.refused {
        return LogName::Unopened;
    }
    let path = walk.path().to_vec();
    if directories.iter().any(|directory| within(directory, &path)) {
        return LogName::Unopened;
    }
    let found = entry.and_then(|entry| entry(&path, &cwd));
    if found == Some(Entry::Directory) {
        return LogName::Unopened;
    }
    let existing = found == Some(Entry::File) && !is_virtual_or_device(&path);
    if !existing && !opens(&path) {
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
        return LogName::StandardError;
    }
    if stream(b"/fd/1".as_slice()) {
        return LogName::Stdout;
    }
    if stream(b"/fd/0".as_slice()) || LOST_DEVICES.contains(&path.as_slice()) {
        return LogName::Nowhere;
    }
    if existing {
        return LogName::Existing;
    }
    LogName::File {
        name: host_file_name(&path, &file.path, &cwd),
        key: file_key(&path, &cwd),
    }
}

/// Whether `directory`, a directory of a name in the working directory `cwd` (`log_name`), is one
/// that the system walks to and the model knows: each directory that its walk passes may be one
/// of root's that the model does not know, but none that is not there nor one that the user
/// cannot search, and the directory where it ends is known (`directory_at`). When `entry` is
/// given, every folder looked up inside `cwd` (or named by a registered file) must be
/// `Entry::Directory`, and an ancestor of a registered file outside `cwd` is present even where
/// `directory_at` returns `Directory::None`. Each directory passed is added to `directories` (the
/// opened path is none of them nor an ancestor of one).
fn walks_to_directory(
    directory: &[u8],
    cwd: &[u8],
    entry: Option<EntryLookup<'_>>,
    directories: &mut Vec<Vec<u8>>,
) -> bool {
    if directory.is_empty() || !fits(directory) {
        return false;
    }
    let mut walk = Walk::at(directory, cwd);
    for part in directory.split(|byte| *byte == b'/') {
        walk.step(part);
        if walk.refused || !folder_present(walk.path(), cwd, entry, false) {
            return false;
        }
        if directories.last().map(Vec::as_slice) != Some(walk.path()) {
            directories.push(walk.path().to_vec());
        }
    }
    folder_present(walk.path(), cwd, entry, true)
}

/// What an `-export_search_strategy` name gives for a `Receiver` (`RunFiles::strategy_target`).
#[derive(Debug, PartialEq, Eq)]
pub enum StrategyTarget {
    /// A file of the run under this host name, the command line's file (`output_name`'s `File`:
    /// `s.asn` for `./s.asn`, `in/../s.asn` and `/proc/self/cwd/s.asn`; the *start*'s name of it,
    /// `RunFiles::output_host_name`), which the receiver gets when the run ends.
    File(Vec<u8>),
    /// No file: the writes are lost (the null device; the standard input, which is the null
    /// device of a command line run without one).
    Nowhere,
    /// A name that NCBI cannot open for writing: NCBI's argument error
    /// (`cli::inaccessible`).
    Inaccessible,
    /// The standard output, which is the report's stream (LOSAT Web has no `-out`): LOSAT Web's
    /// refusal of a strategy file that is also the output (`strategy::open_strategy_file_in`,
    /// `strategy::create_strategy_file_in`).
    StandardOutput,
    /// A file that is not a file of the run: the standard error, a device whose writes fail, a
    /// file of procfs or of a message queue (LOSAT Web's refusal).
    System,
    /// The file of another option of the run, whose bytes the command line's strategy file
    /// replaces or mixes with its own (LOSAT Web's refusal): `-query`, `-subject`, a PSSM file
    /// (`-in_pssm`, `-in_msa`, `-out_pssm`, `-out_ascii_pssm`), `-logfile` or `-window_masker_db`.
    OtherFile(&'static str),
    /// Any other file of the run's files (`RunView`: `run_entries`, `run_directory_entry`), or a
    /// file that the search reads (`FileReceiver::read_files`), which the command line's strategy
    /// file empties or makes before the search reads it (LOSAT Web's refusal, as `-logfile`'s
    /// `LogName::Existing` and a tool's `OutputName::Existing`: LOSAT Web reads the run's files and
    /// writes none of them).
    Existing,
}

/// A file that the run reads or writes besides its strategy file (`RunFiles::strategy_target`):
/// the option that names it and its name as the command line takes it.
pub struct RunFileName<'a> {
    pub option: &'static str,
    pub name: &'a OsStr,
}

impl RunFiles {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535 (cited
    // at `strategy::create_strategy_file`): an output file (`CArgDescriptions::eOutputFile`), which
    // NCBI's arguments open as a stream of that name, truncated.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:623-630
    /// What the `-export_search_strategy` name `name` gives for a `Receiver`, as the command line
    /// on Linux opens it for writing in the run's directory, in a run whose other files are
    /// `others` (`StrategyTarget`). LOSAT's own names (`-`) are the strategy module's. NCBI opens
    /// the name as it is written (no reading of NCBI's, unlike `-logfile`), for writing, creating
    /// it and none of its folders (the option has no `fCreatePath`), as a tool opens its output
    /// (`output_name`), against the run's files through the one view of the run's directory that
    /// `-logfile` and a tool's output read (`RunView`; line WQA D19): a folder
    /// inside the working directory must be one of the run's files (`out/s.asn` without a folder
    /// `out`, `q.fa/s.asn` and `a/b/../s.asn` without `a/b` are NCBI's error), a folder of the
    /// run's files is a directory, and a file of them is LOSAT Web's refusal (`Existing`). Then the
    /// standard output (`/dev/stdout`, `/dev/fd/1`) is the report's stream; the standard input and
    /// the devices of `DISCARDING_DEVICES` take the bytes and keep none; the standard error, the
    /// other devices that a user opens, the files of procfs (`/proc/self/comm`) and message queues
    /// (`/dev/mqueue/q`) are no file of the run; the file of another option of the run is that
    /// option's (`OtherFile`: compared by path, `file_key`; the `-logfile` name read as
    /// `received_log_target` reads it, so a `-logfile` that stays unopened does not match), then a
    /// file that the search reads (`FileReceiver::read_files`) is `Existing`; every other name is a
    /// file where the user makes files (`/tmp/s.asn`, `out/s.asn` with a folder `out`,
    /// `/proc/self/cwd/s.asn`), named as the command line's file (`File`).
    pub fn strategy_target(&self, name: &OsStr, others: &[RunFileName<'_>]) -> StrategyTarget {
        let view = match self {
            RunFiles::Disk => RunView::default(),
            RunFiles::Receiver(receiver) => RunView::new(receiver.as_ref(), receiver.read_files()),
        };
        let name = name.as_encoded_bytes();
        let (made, key) = match view.output_name(name) {
            OutputName::Unopened => return StrategyTarget::Inaccessible,
            OutputName::Stdout => return StrategyTarget::StandardOutput,
            OutputName::Discard => return StrategyTarget::Nowhere,
            OutputName::Full | OutputName::StandardError | OutputName::Special => {
                return StrategyTarget::System;
            }
            // Each name is read in a working directory of its own depth; their paths compare as
            // the run's files' (`file_key`). An existing file of the run is at the path of the
            // same walk (`output_key`).
            OutputName::Existing => (None, output_key(name).unwrap_or_default()),
            OutputName::File { name, key } => (Some(name), key),
        };
        for other in others {
            let bytes = other.name.as_encoded_bytes();
            let same = match other.option != "logfile" {
                // The standard input (`-`) is no file of a name.
                true => bytes != b"-" && output_key(bytes).is_some_and(|other| other == key),
                // NCBI's special names (`-`, `/dev/null`) are no file of the run.
                false => {
                    special_log_target(Some(other.name)).is_none()
                        && matches!(
                            view.log_name(bytes),
                            LogName::File { key: other, .. } if other == key
                        )
                }
            };
            if same {
                return StrategyTarget::OtherFile(other.option);
            }
        }
        match made {
            Some(host_name) if !view.reads(&key) => StrategyTarget::File(host_name),
            _ => StrategyTarget::Existing,
        }
    }
}

#[cfg(test)]
#[path = "run_files_tests.rs"]
mod tests;

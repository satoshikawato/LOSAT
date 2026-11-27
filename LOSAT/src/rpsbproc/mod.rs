//! `losat rpsbproc`: NCBI's post-RPSBLAST processing utility rpsbproc 0.5.1 (the package of
//! 2025-05-13 on `ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/rpsbproc/`), which reads the ASN.1 text
//! archives of `rpsblast -outfmt 11` and the CDD annotation data and writes the conserved domain
//! hits of each query. LOSAT gives the domain hits (`-t doms`, with `-m rep`, `std` or `full`
//! and `-f`) and the sites and motifs (`-t feats`, `-t both`, the default) of protein queries
//! (rpsblast's archives) and of nucleotide queries (rpstblastn's archives, per reading frame),
//! and reads the BlastOutput XML of `-outfmt 5` too, which NCBI's rpsbproc reads as deprecated
//! (`blastxml`, `offline`).
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1853-1871
//!
//! NCBI's rpsbproc is not part of BLAST+ and its source is not in the pinned BLAST tree; the
//! references `RpsbProc 0.5.1, 2025-05-13` name the files of its source package
//! (`RpsbProc-src.tar.gz` of 2025-05-13, sha256 a20f4088...6086), the same source that bioconda
//! builds on the BLAST+ 2.17.0 toolkit.

mod archive;
mod args;
mod blastxml;
mod data;
mod domains;
mod frames;
mod nucleotide;
mod offline;
mod prosite;
pub(crate) use prosite::{HOST_DEPENDENT_BYTES, USER_SPACE_BYTES};
mod query;
mod segset;
mod sites;

pub use args::RpsbprocArgs;

use std::ffi::{OsStr, OsString};
use std::io::{BufReader, Read, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_exception::{err_post_x, err_post_x_bytes, NcbiException};
use crate::report::serial::asn_types::read::{Post, PostSink};

use archive::{ArchiveStream, ParamValue, Queries, RequestBody, SeqEntry, Subject};
use args::TargetData;
use blastxml::XmlStream;
use data::{CddData, DataError, DataFiles};
use domains::{
    domain_records, feature_records, map_seg_set, translated_pos_map, DataMode, FeatureRecord,
    Record,
};
use prosite::{MapEnd, MatchError};
use query::{
    check_masks, choice_selection, process_blast_results, unassigned_member, DomQuery, ProcessError,
};
use segset::SegSet;
use sites::{MotifError, MotifState};

/// `EErrorSubCode`: the exit codes.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:91-101
const E_OK: i32 = 0;
const E_DATA_LOAD_ERROR: i32 = 2;
const E_INPUT_FILE_OPEN_ERROR: i32 = 3;
pub(crate) const E_OUTPUT_FILE_OPEN_ERROR: i32 = 4;
const E_INPUT_DATA_ERROR: i32 = 5;

/// `CRpsbProcApp::COORDSBASE`: the printed coordinates count from 1.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:335
const COORDSBASE: i64 = 1;

/// `position + COORDSBASE`, the `int` sum of every coordinate that rpsbproc prints in `int`
/// (`m_iFrom` and `m_iTo` plus `COORDSBASE`, of the DOMAINS and SUPERFAMILIES lines, the
/// residues of the SITES and MOTIFS lines, the first base of a nucleotide query's codon and motif).
/// NCBI's sum overflows at `INT_MAX` (undefined in C++; both 0.5.1 builds print -2147483648, PS16c
/// review oracle `plusmax`): LOSAT adds in `i64`, so it prints 2147483648, as NCBI prints the
/// position of a hit moved away from the limit (PD-LOSAT-NCBI-DEFECTS exception 12, the Owner's
/// rule of 2026-10-10 and Q1 = B of 2026-11-14). `position` holds an `int`'s value, past `int`'s
/// range where NCBI's arithmetic took it there (exception 12 too).
fn printed_position(position: impl Into<i64>) -> i64 {
    position.into() + COORDSBASE
}

/// `OVERLAPLEADING` as rpsbproc stores it (`unsigned int` of 10000000000).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:338
const OVERLAPLEADING: u32 = (10_000_000_000u64 % (1u64 << 32)) as u32;

/// How `losat rpsbproc` ends: an exit status, or a refusal (its message, which may quote data
/// bytes; exit 1).
pub enum Outcome {
    Exit(i32),
    Refused(Vec<u8>),
}

impl Outcome {
    fn refused(message: impl Into<Vec<u8>>) -> Self {
        Outcome::Refused(message.into())
    }
}

/// glibc's `%g` (an `ostream`'s default output of a double), `-0` for negative zero.
fn format_g(value: f64) -> String {
    crate::report::pairwise::cpp_default_double(value)
}

/// The writer of the report (`m_ostream`): standard output or the `-o` file, with the run's
/// diagnostics (`diag`). LOSAT writes the report in blocks; before each line it writes to
/// standard error or the diagnostics, the report's file holds what NCBI's stream has written by
/// then (`cli::NcbiBuffer`): every line, which NCBI ends with `endl`, and of the line being written
/// what the stream's buffer did not keep. `cerr` is tied to `cout`, so before a write to `cerr`
/// a report on standard output is written whole. The streams then interleave as NCBI's do, also
/// in one file (`-o F -logfile F`, `-logfile /dev/stdout`, `2>&1`): NCBI rpsbproc 0.5.1 with an
/// exception in the middle of the `SESSION` line keeps that part in its stream, writes the
/// report of the exception at the end of the `-logfile` file and then the rest at the stream's
/// offset, over the start of the report (line RG E1c review 1, `m_subj_same`). The buffers are
/// those of NCBI's builds of rpsbproc, which are Linux builds (glibc, GCC 13's libstdc++): LOSAT
/// keeps their counts on every platform (NCBI's rpsbproc on macOS or Windows is not measured).
///
/// A failed write stops LOSAT (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 9): NCBI's
/// `ostream` sets its bad bit and the run goes on to exit 0 (a closed pipe ends it by SIGPIPE).
struct Report {
    out: Box<dyn Write>,
    /// The report's bytes not yet written to `out`.
    pending: Vec<u8>,
    /// NCBI's buffer of the stream, which keeps the last bytes of `pending`
    /// (`cli::NcbiBuffer::held`): what of a line that is being written has reached the file.
    buffer: crate::cli::NcbiBuffer,
    quiet: bool,
    /// The first write error; nothing is written after it.
    failed: Option<std::io::Error>,
    diag: Diagnostics,
}

/// The size of the blocks in which LOSAT writes the report.
const REPORT_BLOCK: usize = 1 << 16;

/// `write_all` of `first` and then `second`, the two in one write where the device takes them so
/// (`writev`): a block of the report is one write wherever its bytes are, so a write that fails
/// cuts the report where it cuts the block (`standard_output`: a writer's write sizes decide it).
/// As `write_all`, a write that is interrupted is made again, and one that takes no byte fails.
fn write_all2(out: &mut dyn Write, mut first: &[u8], mut second: &[u8]) -> std::io::Result<()> {
    while !first.is_empty() || !second.is_empty() {
        let written = if first.is_empty() {
            out.write(second)
        } else if second.is_empty() {
            out.write(first)
        } else {
            out.write_vectored(&[std::io::IoSlice::new(first), std::io::IoSlice::new(second)])
        };
        match written {
            Ok(0) => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::WriteZero,
                    "failed to write whole buffer",
                ))
            }
            Ok(n) => {
                let from_first = n.min(first.len());
                first = &first[from_first..];
                second = &second[n - from_first..];
            }
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(())
}

impl Report {
    /// A report written to `out` through NCBI's `stream`.
    fn new(
        out: Box<dyn Write>,
        stream: crate::cli::NcbiStdio,
        quiet: bool,
        diag: Diagnostics,
    ) -> Self {
        Report {
            out,
            pending: Vec::with_capacity(REPORT_BLOCK),
            buffer: crate::cli::NcbiBuffer::new(stream),
            quiet,
            failed: None,
            diag,
        }
    }

    /// An insertion of a string or a number (`sputn`).
    fn write(&mut self, bytes: &[u8]) {
        if self.failed.is_some() {
            return;
        }
        // (`sputn` reads no more than the insertion's own bytes of the stream.)
        self.buffer.sputn(bytes.len(), bytes);
        // An insertion that ends a line is `endl`, which flushes the stream. LOSAT writes the
        // report in the pieces of NCBI's insertions where a piece can reach a buffer's limit: the
        // strings of the archive and the data. NCBI rpsbproc 0.5.1, both builds, on a `SESSION`
        // line cut by an exception after a database name of 1,100, 5,000 and 9,000 bytes, with
        // `-o F -logfile F` and with `-logfile /dev/stdout` on a pipe: line RG `E/e1c/bufprobe`.
        if bytes.last() == Some(&b'\n') {
            self.buffer.flush();
        }
        let held = self.buffer.held();
        if self.pending.len() + bytes.len() - held < REPORT_BLOCK {
            self.pending.extend_from_slice(bytes);
        } else {
            self.write_through(bytes, held);
        }
    }

    /// The write of a full block or more that NCBI's stream has written, ending in the insertion
    /// `bytes`: the bytes before the last `held` ones (which NCBI's buffer holds) go to the file
    /// from `pending` and from `bytes` where they are, in one write as the block would
    /// (`write_all2`), and only the held ones are kept, so a long string of the input (a title of
    /// many megabytes) is not copied into the report as well.
    fn write_through(&mut self, bytes: &[u8], held: usize) {
        let held_of_bytes = held.min(bytes.len());
        let (written, kept) = bytes.split_at(bytes.len() - held_of_bytes);
        let n = self.pending.len() - (held - held_of_bytes);
        if let Err(error) =
            write_all2(&mut *self.out, &self.pending[..n], written).and_then(|()| self.out.flush())
        {
            self.failed = Some(error);
        }
        self.pending.drain(..n);
        self.pending.extend_from_slice(kept);
    }

    /// An insertion of a `char` (a stream insertion of `DELIMIT`: `put`), where a line can be cut
    /// by an exception (the `SESSION` line; in the other lines it holds what `write` would).
    fn put(&mut self, byte: u8) {
        if self.failed.is_some() {
            return;
        }
        self.pending.push(byte);
        self.buffer.sputc(&self.pending);
        self.write_blocks();
    }

    /// Writes the report's full blocks that NCBI's stream has written.
    fn write_blocks(&mut self) {
        if self.pending.len() - self.buffer.held() >= REPORT_BLOCK {
            self.write_out(self.pending.len() - self.buffer.held());
        }
    }

    /// Writes the first `n` bytes of `pending` to the file.
    fn write_out(&mut self, n: usize) {
        if self.failed.is_none() && n > 0 {
            if let Err(error) =
                write_all2(&mut *self.out, &self.pending[..n], &[]).and_then(|()| self.out.flush())
            {
                self.failed = Some(error);
            }
        }
        self.pending.drain(..n);
    }

    /// The file holds what NCBI's stream has written by now.
    fn sync(&mut self) {
        self.write_out(self.pending.len() - self.buffer.held());
    }

    /// Before a write to `cerr`, which is tied to `cout`: a report on standard output is written
    /// whole (its stream flushed), one in the `-o` file as far as NCBI's stream has written it.
    fn before_cerr(&mut self) {
        match self.buffer.stream() {
            crate::cli::NcbiStdio::Cout { .. } => self.flush(),
            crate::cli::NcbiStdio::File => self.sync(),
        }
    }

    /// Before a post of the diagnostics: as before a write to `cerr` when they go there, else
    /// as far as NCBI's stream has written the report (the `-logfile` file is not tied to it).
    fn before_post(&mut self) {
        if self.diag.is_cerr() {
            self.before_cerr();
        } else {
            self.sync();
        }
    }

    fn line(&mut self, parts: &[&[u8]]) {
        for part in parts {
            self.write(part);
        }
        self.write(b"\n");
    }

    /// Writes the whole report (NCBI's stream flushed, as at the end of the program).
    fn flush(&mut self) {
        self.write_out(self.pending.len());
        self.buffer.flush();
    }

    /// A progress line on standard error (`cerr`; not with `-q`, nor after a failed write).
    fn progress(&mut self, text: &str) {
        if !self.quiet {
            self.before_cerr();
            if self.failed.is_none() {
                cerr_line(text);
            }
        }
    }

    /// A post of the run's diagnostics (`Diagnostics::post`).
    fn post(&mut self, text: &str) {
        self.before_post();
        self.diag.post(text);
    }

    /// Whether a write failed: the run stops (`run` reports it). What NCBI's stream has written
    /// by now is written first.
    fn stopped(&mut self) -> bool {
        self.sync();
        self.failed.is_some()
    }
}

/// What NCBI's reader posts while it reads an archive, written when it is posted (the stream's
/// posts, `ArchiveStream::posts`): before rpsbproc uses the object, also when the object then
/// fails to read. rpsbproc sets no prefix for the posts (`Post::line`).
impl PostSink for Report {
    fn post(&mut self, post: Post) {
        self.before_post();
        self.diag.post_bytes(&post.line(None));
    }
}

/// A line that rpsbproc writes to `cerr` (`text` streamed to `cerr`, then `endl`): standard error,
/// whatever the `-logfile`.
fn cerr_line(text: &str) {
    write_stderr(&format!("{text}\n"));
}

/// Where rpsbproc's diagnostics go: its `LOG_POST_X` messages and the application's posts and
/// reports of an exception (`ERR_POST_X`, `NCBI_REPORT_EXCEPTION_X`), to standard error or, with
/// `-logfile`, to that file or nowhere (`ncbi_environment::DiagTarget`, set up by the command
/// line's entry before the arguments, as NCBI's `SetupDiag`). rpsbproc's `cerr` lines (the
/// program information, the progress, `DryRun`'s line: `cerr_line`) and LOSAT's own messages stay
/// on standard error. NCBI's rpsbproc 0.5.1 with `-logfile log.txt` (line RG `E/SURVEY_E1c.md`,
/// probes `c_log_*`): `Input file does not exist`, `Cannot find data files`, the reports of an
/// exception in `Init` or `Run` go into the file, the banner and `RPSBlast session 1 done...` to
/// standard error.
pub struct Diagnostics(crate::blastinput::ncbi_environment::DiagTarget);

impl Diagnostics {
    /// The diagnostics of `target`.
    pub fn new(target: crate::blastinput::ncbi_environment::DiagTarget) -> Self {
        Diagnostics(target)
    }

    /// Whether the diagnostics go to `cerr` (no `-logfile`, or one that NCBI could not open).
    fn is_cerr(&self) -> bool {
        matches!(
            self.0,
            crate::blastinput::ncbi_environment::DiagTarget::Stderr
        )
    }

    /// `LOG_POST_X`: rpsbproc's messages are its text alone (the old post format).
    fn post(&mut self, text: &str) {
        self.post_line_bytes(text.as_bytes());
    }

    /// `LOG_POST_X` of raw bytes (adds the trailing newline).
    fn post_line_bytes(&mut self, bytes: &[u8]) {
        let mut line = Vec::with_capacity(bytes.len() + 1);
        line.extend_from_slice(bytes);
        line.push(b'\n');
        self.post_bytes(&line);
    }

    /// Writes `bytes` where the diagnostics go; on standard error as `write_stderr_bytes`
    /// writes (a closed one ends LOSAT, exception 9).
    fn post_bytes(&mut self, bytes: &[u8]) {
        self.0.write_text(bytes, &mut Cerr);
    }
}

/// rpsbproc's standard error as a stream (`write_stderr_bytes`).
struct Cerr;

impl Write for Cerr {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        write_stderr_bytes(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Writes `text` to standard error. A failed write (a reader that has closed its pipe) ends
/// LOSAT at once with exit status 4 and no message, where NCBI's rpsbproc dies of SIGPIPE
/// (exit status 141 from the shell; PS16a review oracle M1): PD-LOSAT-CLI-NONSEARCH-DIFFERENCES
/// exception 9. The report is flushed before every line that rpsbproc writes to standard error
/// (`Report::progress`, `Report::stopped`), so it holds no bytes then.
pub(crate) fn write_stderr(text: &str) {
    write_stderr_bytes(text.as_bytes());
}

/// `write_stderr` of bytes (a message that quotes data bytes as NCBI writes them).
pub(crate) fn write_stderr_bytes(bytes: &[u8]) {
    let mut stderr = std::io::stderr().lock();
    if stderr.write_all(bytes).is_err() {
        std::process::exit(E_OUTPUT_FILE_OPEN_ERROR);
    }
}

/// Writes `text` (help or version) to standard output; a failed write ends LOSAT as a failed
/// write of the report does (exception 9).
pub(crate) fn write_stdout(text: &str) {
    let mut stdout = crate::cli::standard_output();
    if let Err(error) = stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.flush())
    {
        write_stderr(&format!("rpsbproc failed to write output: {error}\n"));
        std::process::exit(E_OUTPUT_FILE_OPEN_ERROR);
    }
}

/// `NCBI_USE_ERRCODE_X` of ncbiapp.cpp (`Corelib_App`), the error code of the application's
/// reports.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:65
// NCBI reference (598d8ae6): c++/include/corelib/error_codes.hpp:49
const CORELIB_APP: u32 = 106;

/// The report of an exception that leaves rpsbproc's `Run` (`CNcbiApplication::x_TryMain`):
/// the exception's own line, then the application's, exit 3. For the exceptions of `Init`
/// (`initialization`) the application's line has subcode 15 and the exit is 2.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:669-674,711-728
// (rpsbproc sets no diagnostic flags, so the report is `NcbiException::report_default_ex`'s; the
// lines as NCBI's rpsbproc printed them: `Error: (CInvalidChoiceSelection::eFail) <message>` and
// `Error: (106.16) Application's execution failed (CInvalidChoiceSelection::eFail) <message>`
// (probe 2, a database with `lcl|` ids); `Error: (CStringException::eConvert) <message>` and
// `Error: (106.15) Application's initialization failed (CStringException::eConvert) <message>`
// (probe 4, data files with numbers that do not convert).)
fn report_exception(
    diag: &mut Diagnostics,
    exception: &NcbiException,
    initialization: bool,
) -> i32 {
    let (subcode, title, status) = if initialization {
        (15, "Application's initialization failed", 2)
    } else {
        (16, "Application's execution failed", 3)
    };
    diag.post_bytes(&exception.report_default_ex(CORELIB_APP, subcode, title));
    status
}

/// `CDir`'s `Exists` on the path.
fn dir_exists(path: &OsStr) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_dir())
}

/// `CFile`'s `Exists` on the path.
fn file_exists(path: &OsStr) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_file())
}

/// `a + "/" + b` of `std::string`s.
fn join(a: &OsStr, b: &str) -> OsString {
    let mut out = a.to_os_string();
    out.push("/");
    out.push(b);
    out
}

/// `CDirEntry::NormalizePath` of the path (`eIgnoreLinks`): the path's `.` and empty components
/// dropped and `..` taken back where possible, without its trailing `/`
/// (`ncbi_environment::normalize_path`, the Unix form, on every target).
// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:501-502
fn normalize_path(path: &[u8]) -> Vec<u8> {
    crate::blastinput::ncbi_environment::normalize_path(path, None)
        .expect("a path normalized without links has no link limit")
}

/// The paths where rpsbproc looks for its data without `-d`: the launch directory of the
/// program (`argv[0]` normalized) and its real directory (links followed, made absolute).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1855-1868
// The binary is `losat` (or a link `rpsbproc` to it, `--install-shims`). The launch path is
// ported as NCBI makes it. For the real path NCBI follows the links of each component of the
// path (`s_DereferencePath`, ncbifile.cpp) relative to the current directory, normalizes it
// lexically and puts the current directory before a relative one; LOSAT takes the canonical
// path of the launch path (`fs::canonicalize`), else the current directory with the path. They
// differ only where a `..` follows a linked directory (the canonical path goes up from the
// link's target, NCBI's lexical normalization from the link).
fn program_paths() -> (OsString, OsString) {
    let argv0 = std::env::args_os().next().unwrap_or_default();
    let alias = normalize_path(argv0.as_encoded_bytes());
    let dir_of = |path: &[u8]| -> Vec<u8> {
        match path.iter().rposition(|&b| b == b'/') {
            Some(at) => path[..at].to_vec(),
            None => Vec::new(),
        }
    };
    // SAFETY: the bytes are an `OsStr`'s (`argv[0]`, or a path made of its parts and `/`).
    let to_os = |bytes: Vec<u8>| unsafe { OsString::from_encoded_bytes_unchecked(bytes) };
    let launch = dir_of(&alias);
    let alias_path = std::path::PathBuf::from(to_os(alias.clone()));
    let real_path = match std::fs::canonicalize(&alias_path) {
        Ok(path) => path.into_os_string().into_encoded_bytes(),
        Err(_) => {
            let mut absolute = std::env::current_dir()
                .map(|dir| dir.into_os_string().into_encoded_bytes())
                .unwrap_or_default();
            if alias.first() == Some(&b'/') {
                absolute = alias.clone();
            } else {
                absolute.push(b'/');
                absolute.extend_from_slice(&alias);
            }
            normalize_path(&absolute)
        }
    };
    (to_os(launch), to_os(dir_of(&real_path)))
}

/// The data directory (`#Annotation data at:`), searched as rpsbproc searches it.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:572-643
// (The registry's `[datafiles]` names, main:535-566, are rejected with every other registry
// entry that changes the output: `ApplicationStart::settings`.)
fn find_data_path(diag: &mut Diagnostics, data_path: Option<&OsStr>) -> Option<OsString> {
    let found = |path: &OsStr| dir_exists(path) && file_exists(&join(path, "cddid.tbl"));
    match data_path {
        Some(datapath) => {
            if datapath.is_empty() {
                return None;
            }
            let candidate = join(datapath, "data");
            if found(&candidate) {
                return Some(candidate);
            }
            if found(datapath) {
                return Some(datapath.to_os_string());
            }
            let mut msg = b"Specified data location (".to_vec();
            msg.extend_from_slice(&crate::utils::os_bytes::os_str_bytes(datapath));
            msg.extend_from_slice(b") does not exist or does not contain data files");
            diag.post_line_bytes(&msg);
            None
        }
        None => {
            for candidate in [OsString::from("./data"), OsString::from("./")] {
                if found(&candidate) {
                    return Some(candidate);
                }
            }
            let (launch, real) = program_paths();
            let candidate = join(&launch, "data");
            if found(&candidate) {
                return Some(candidate);
            }
            if real != launch {
                let candidate = join(&real, "data");
                if found(&candidate) {
                    return Some(candidate);
                }
                if found(&real) {
                    return Some(real);
                }
            }
            None
        }
    }
}

/// The 9 lines that rpsbproc writes to standard error first (also with `-q`); the first names
/// LOSAT's version (`Post-RPSBLAST Processing Utility (LOSAT 0.2.0)`, `wording::RPSBPROC_TITLE`
/// and `REPORT_VERSION`) where NCBI's names its own (`v0.5.1`, `VERDEF`; the Owner 2026-11-21
/// 「Q1=B」).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:74-87 (`VERDEF`)
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:397-406
// (`m_silent` is still false here: `-q` is read at main:677.)
fn banner() {
    write_stderr(&format!(
        "{} {}\n\
         -------------------------------------------------------\n\
         Download:\n\
         \tftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/\n\
         \n\
         Blast applications download:\n\
         \tftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n\
         \n\
         For more information please refer to the README file\n",
        crate::wording::RPSBPROC_TITLE,
        crate::wording::REPORT_VERSION
    ));
}

/// The header of the report.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:679-716
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:699-711
fn header_tail(evalue: f64, mode: DataMode, show_families: bool, target: TargetData) -> String {
    let mode = match mode {
        DataMode::Rep => "Concise",
        DataMode::Std => "Standard",
        DataMode::Full => "Full",
    };
    let (requested, domains, sites) = match target {
        TargetData::Doms => ("Domain hits", true, false),
        TargetData::Feats => ("Site annotations", false, true),
        TargetData::Both => ("Domain hits and site annotations", true, true),
    };
    let mut text = format!(
        "#E-Value cutoff:\t{}\n\
         #Redundancy:\t{mode}\n\
         #Data requested:\t{requested}\n\
         #Output format -- tab-delimited table\n\
         #Show superfamilies: {}\n\
         #DATA\n\
         #SESSION\t<session-ordinal>\t<program>\t<database>\t<score-matrix>\t<evalue-threshold>\n\
         #QUERY\t<query-id>\t<seq-type>\t<seq-length>\t<definition-line>\n",
        format_g(evalue),
        if show_families { "YES" } else { "NO" }
    );
    if domains {
        text.push_str(
            "#DOMAINS\n\
             #<session-ordinal>\t<query-id[readingframe]>\t<hit-type>\t<PSSM-ID>\t<from>\t<to>\t<E-Value>\t<bitscore>\t<accession>\t<short-name>\t<incomplete>\t<superfamily PSSM-ID>\n\
             #more such lines......\n\
             #ENDDOMAINS\n",
        );
    }
    if sites {
        text.push_str(
            "#SITES\n\
             #<session-ordinal>\t<query-id[readingframe]>\t<annot-type>\t<title>\t<residue(coordinates)>\t<complete-size>\t<mapped-size>\t<source-domain>\n\
             #more such lines......\n\
             #ENDSITES\n\
             #MOTIFS\n\
             #<session-ordinal>\t<query-id[readingframe]>\t<title>\t<from>\t<to>\t<source-domain>\n\
             #more such lines......\n\
             #ENDMOTIFS\n",
        );
    }
    text.push_str(
        "#ENDQUERY\t<query-id>\n\
         #more query sections..\n\
         #ENDSESSION\t<session-ordinal>\n\
         #more session sections..\n\
         #ENDDATA\n\
         #=====================================================================\n",
    );
    text
}

/// The first variable of `variables` that changes how NCBI's object stream reads an archive or
/// what its objects give, with a value that changes it; `None` when there is none.
// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:57-68
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:257-287
// (`CSerialObject::x_GetVerifyData`, serialobject.cpp:180-221, reads `SERIAL_VERIFY_DATA_GET`
// the same way; `NO`, `NEVER`, `DEFVALUE` and `DEFVALUE_ALWAYS` turn the checks off
// (`ncbi_environment::verification_off`): a missing member is no longer a read error, and the
// getter of an unset member gives its default instead of throwing `CUnassignedMember`. `YES`,
// `ALWAYS` and other values keep NCBI's default.)
// rpsbproc refuses the whole run where the other programs refuse only an input that fails a check
// (`ncbi_environment::ReadVerification`): its readers of the archive and of the BLAST XML report
// read past members that they do not keep without checking them, where NCBI's check of a
// missing mandatory member stops the read by default and only posts under these values; a run
// whose input would pass the check is not told apart from one that would not.
// `OBJECTS_SEQ_DESCR_ALLOW_EMPTY` lets an empty `Seq-descr` read; LOSAT rejects it with any value
// when the input is an archive (`xml` false: a BlastOutput document has no `Seq-descr`).
// The other parameters of the readers are registry entries (`NCBI_CONFIG__SERIAL__...`, which
// the application start rejects, `ApplicationStart::settings`) or change no result
// (`SERIAL_READ_MMAPBYTESOURCE`, `SERIAL_DISABLE_DELAY_BUFFERS` and the binary ASN.1 ones); the
// gi parameters of the Seq-ids are refused in `run_sessions`, where a query has a gi.
fn reader_variable(
    variables: impl Iterator<Item = (OsString, OsString)>,
    xml: bool,
) -> Option<String> {
    for (name, value) in variables {
        let (Some(name), value) = (name.to_str(), value.to_string_lossy()) else {
            continue;
        };
        let changes = match name {
            "OBJECTS_SEQ_DESCR_ALLOW_EMPTY" => !xml,
            "SERIAL_VERIFY_DATA_READ" | "SERIAL_VERIFY_DATA_GET" => {
                crate::blastinput::ncbi_environment::verification_off(&value)
            }
            _ => false,
        };
        if changes {
            return Some(format!(
                "the environment variable {name}={value}, which changes how NCBI's object stream reads its input, is not supported by LOSAT's rpsbproc"
            ));
        }
    }
    None
}

/// What rpsbproc's `Init` leaves for `Run` and `DryRun` (`init`): its status (`m_app_status`),
/// and when that is `e_OK`, what the run reads and writes.
enum Init {
    /// A status other than `e_OK`: `Run` returns it, `DryRun` writes it.
    Failed(i32),
    /// `e_OK`.
    Ready(Box<Ready>),
}

/// What `Run` takes from a successful `Init`.
struct Ready {
    /// The report's header (`m_init_msg`).
    header: Vec<u8>,
    input: BufReader<Box<dyn Read>>,
    /// The first byte of the input that is not white space, read by `Init`.
    first: u8,
    out: Box<dyn Write>,
    data: CddData,
    /// Whether the input is XML (`first` is `<`): `Run` reads BlastOutput documents.
    xml: bool,
    seq_id_rank_variable: Option<&'static str>,
}

/// Runs `losat rpsbproc` with `args` (NCBI's `Init` and `Run`), its diagnostics going to `diag`.
/// `settings` are those of rpsbproc's application start, or LOSAT's refusal of them
/// (`ApplicationStart::settings`).
pub fn run(
    args: RpsbprocArgs,
    settings: Result<ApplicationSettings, String>,
    mut diag: Diagnostics,
) -> Outcome {
    match init(&args, settings, &mut diag) {
        Err(outcome) => outcome,
        Ok(Init::Failed(status)) => Outcome::Exit(status),
        Ok(Init::Ready(ready)) => run_ready(&args, *ready, diag),
    }
}

// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:924-932
/// `losat rpsbproc -dryrun` with `args`: rpsbproc's own `DryRun`, which NCBI's application runs
/// in place of `Run` after `Init` (ncbiapp.cpp:705-707, cited at `cli::main_entry`). `Init`
/// runs as for a run (the program information, its messages, the `-i` file opened and read up
/// to its first byte that is not white space, the `-o` file made empty, the data files loaded),
/// then `DryRun`'s line on standard error, and its status is the exit status. An exception that
/// leaves `Init` and LOSAT's refusals in `Init` (a file name, a setting that NCBI applies before
/// `DryRun`) end the program as in a run. NCBI rpsbproc 0.5.1 (both builds; line RG
/// `E/SURVEY_E1c.md`, probes `c_dry_*`): an archive that `Run` would reject or LOSAT's `Run`
/// would refuse (`arch/lcl.asn`, a truncated one, the binary archive) passes, and so does the
/// XML input;
/// `-o out.txt` leaves an empty file; a missing `-i` file with an empty standard input fails
/// with 5 (the input's), with an archive on it with 3.
pub fn dry_run(
    args: RpsbprocArgs,
    settings: Result<ApplicationSettings, String>,
    mut diag: Diagnostics,
) -> Outcome {
    let status = match init(&args, settings, &mut diag) {
        Err(outcome) => return outcome,
        Ok(Init::Failed(status)) => status,
        Ok(Init::Ready(_)) => E_OK,
    };
    if status == E_OK {
        cerr_line("Dryrun passed, all resources are present and loaded successfully.");
    } else {
        cerr_line(&format!("Dryrun failed, error subcode = {status}"));
    }
    Outcome::Exit(status)
}

/// rpsbproc's `Init` with `args`: LOSAT's checks of the names and the settings that NCBI's
/// application reads before `DryRun` too (`settings`, those of the application start), the
/// program information, the input and output files, the data files and the report's header.
/// `Err` with
/// how the program ends there, before `Run` or `DryRun`: LOSAT's refusal, or the report of an
/// exception that leaves `Init` (exit status 2).
fn init(
    args: &RpsbprocArgs,
    settings: Result<ApplicationSettings, String>,
    diag: &mut Diagnostics,
) -> Result<Init, Outcome> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044 (cited at the BLAST
    // programs in `cli::main_entry`): the application layer reads its environment and its
    // registry files (`rpsbproc.ini`, `.ncbirc`) before `Init`.
    let settings = match settings {
        Ok(settings) => settings,
        Err(reason) => return Err(Outcome::refused(reason)),
    };
    banner();
    let mut status = E_OK;
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:447-498
    let mut header: Vec<u8> = b"#Post-RPSBLAST Processing Utility\n\n#Input data file:".to_vec();
    let mut input: Option<Box<dyn Read>> = None;
    if let Some(infile) = &args.infile {
        let name = crate::utils::os_bytes::os_str_bytes(infile);
        if file_exists(infile) {
            match std::fs::File::open(infile) {
                Ok(file) => {
                    input = Some(Box::new(file));
                    header.push(b'\t');
                    header.extend_from_slice(&name);
                }
                Err(_) => {
                    diag.post_line_bytes(&[b"Unable to open input file: ", &name[..]].concat());
                    status = E_INPUT_FILE_OPEN_ERROR;
                }
            }
        } else {
            diag.post_line_bytes(&[b"Input file does not exist: ", &name[..]].concat());
            status = E_INPUT_FILE_OPEN_ERROR;
        }
    }
    let input: Box<dyn Read> = match input {
        Some(file) => file,
        None => {
            header.extend_from_slice(b"\tstdin");
            Box::new(std::io::stdin())
        }
    };
    header.extend_from_slice(b"\n#Output data file:");
    let mut out: Box<dyn Write> = Box::new(crate::cli::standard_output());
    if let Some(outfile) = &args.outfile {
        let name = crate::utils::os_bytes::os_str_bytes(outfile);
        match std::fs::File::create(outfile) {
            Ok(file) => {
                out = Box::new(file);
                header.push(b'\t');
                header.extend_from_slice(&name);
            }
            Err(_) => {
                diag.post_line_bytes(&[b"Unable to open output file: ", &name[..]].concat());
                status = E_OUTPUT_FILE_OPEN_ERROR;
            }
        }
    } else {
        header.extend_from_slice(b"\t(stdout)");
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:501-522
    // (`c` is a signed `char`: the bytes from 0x80 are skipped as white space.)
    let mut input = BufReader::with_capacity(1 << 16, input);
    let first = loop {
        let mut byte = [0u8];
        match input.read(&mut byte) {
            Ok(1) => {
                if (byte[0] as i8) > b' ' as i8 {
                    break Some(byte[0]);
                }
            }
            Ok(_) => break None,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(_) => break None,
        }
    };
    let Some(first) = first else {
        return Ok(Init::Failed(E_INPUT_DATA_ERROR));
    };
    // The warning comes from `Init`: after a failed `-i`, `-o` or data load NCBI posts it and
    // ends with the status (PS16a review oracle M4), and its `DryRun` posts it and passes.
    let xml = first == b'<';
    if xml {
        diag.post("Warning: XML format support is deprecated!");
    }
    header.push(b'\n');
    header.extend_from_slice(b"#Annotation data at:");
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:645-669
    let mut data = CddData::default();
    match find_data_path(diag, args.data_path.as_deref()) {
        Some(path) => {
            let base = crate::utils::os_bytes::os_str_bytes(&path);
            header.push(b'\t');
            header.extend_from_slice(&base);
            let file = |name: &str| [&base[..], b"/", name.as_bytes()].concat();
            let files = DataFiles {
                cddids: file("cddid.tbl"),
                spfeats: file("cddannot.dat"),
                genfeats: file("cddannot_generic.dat"),
                bscores: file("bitscore_specific.txt"),
                cdtrack: file("cdtrack.txt"),
                famlinks: file("family_superfamily_links"),
            };
            match CddData::load(&files) {
                Ok(loaded) => data = loaded,
                Err(DataError::Missing(path)) => {
                    diag.post_line_bytes(
                        &[
                            b"Annotation data loading unsuccessful: Data file ".as_slice(),
                            &path,
                            b" does not exist",
                        ]
                        .concat(),
                    );
                    status = E_DATA_LOAD_ERROR;
                }
                Err(DataError::Conversion(error)) => {
                    return Err(Outcome::Exit(report_exception(diag, &error, true)));
                }
                // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:675-679
                // (A `std::exception` of the loader is no `CSimpleException`, so it leaves
                // `Init`; NCBI's rpsbproc printed `Error: (106.6) ...` with the `what()` of
                // libstdc++, PS16a review oracle A4-1 cases `annot10`, `gen10` and `ct3`.)
                Err(DataError::Std(what)) => {
                    diag.post_bytes(
                        err_post_x(
                            CORELIB_APP,
                            6,
                            &format!("Application's initialization failed: {what}"),
                        )
                        .as_bytes(),
                    );
                    return Err(Outcome::Exit(2));
                }
                Err(DataError::Unsupported(message)) => return Err(Outcome::refused(message)),
            }
        }
        None => {
            diag.post("Cannot find data files");
            status = E_DATA_LOAD_ERROR;
        }
    }
    header.push(b'\n');
    if status != E_OK {
        return Ok(Init::Failed(status));
    }
    header.extend_from_slice(
        header_tail(args.evalue, args.mode(), args.show_families, args.target()).as_bytes(),
    );
    Ok(Init::Ready(Box::new(Ready {
        header,
        input,
        first,
        out,
        data,
        xml,
        seq_id_rank_variable: settings.seq_id_rank_variable,
    })))
}

/// rpsbproc's `Run` after a successful `Init` (`ready`), its diagnostics going to `diag`: the
/// archives of the input (`run_sessions`), or its BlastOutput XML documents
/// (`run_xml_sessions`). LOSAT refuses there, before NCBI's `Run` reads the input, the
/// environment variables that change how NCBI's object stream reads it (`reader_variable`).
/// `DryRun` reads no input (`dry_run`).
fn run_ready(args: &RpsbprocArgs, ready: Ready, diag: Diagnostics) -> Outcome {
    let Ready {
        header,
        input,
        first,
        out,
        data,
        xml,
        seq_id_rank_variable,
    } = ready;
    if let Some(reason) = reader_variable(std::env::vars_os(), xml) {
        return Outcome::refused(reason);
    }
    // `Ready` has the `-o` file when the run names one (`init` failed otherwise).
    let stream = if args.outfile.is_some() {
        crate::cli::NcbiStdio::File
    } else {
        crate::cli::NcbiStdio::cout()
    };
    let mut report = Report::new(out, stream, args.quiet, diag);
    let outcome = if xml {
        run_xml_sessions(
            &mut report,
            &header,
            XmlStream::new(input, &[first]),
            &data,
            args.mode(),
            args.show_families,
            args.target(),
            args.evalue,
        )
    } else {
        run_sessions(
            &header,
            ArchiveStream::new(input, &[first], &mut report),
            &data,
            args.mode(),
            args.show_families,
            args.target(),
            seq_id_rank_variable,
        )
    };
    // PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 9 (version 1.10): the first failed write of
    // the report stops LOSAT with this message and exit 4, rpsbproc's code of the output file
    // (`e_OutputFileOpenError`); NCBI prints nothing for it.
    report.flush();
    if let Some(error) = report.failed {
        write_stderr(&format!("rpsbproc failed to write output: {error}\n"));
        return Outcome::Exit(E_OUTPUT_FILE_OPEN_ERROR);
    }
    outcome
}

/// The name of a `Blast4-value` variant.
fn value_variant(value: &ParamValue) -> &'static str {
    match value {
        ParamValue::String(_) => "string",
        ParamValue::Cutoff(_) => "cutoff",
        ParamValue::Other(name) => name,
    }
}

/// The sessions: one per archive (`Run`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:727-739
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:793-921
// The `request` and `results` members are mandatory in the module, so an object that reads has
// both (main:886-896 are not reached). The refusals that a session's request decides come
// before anything of that session.
fn run_sessions<R: Read>(
    header: &[u8],
    mut stream: ArchiveStream<R, &mut Report>,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    seq_id_rank_variable: Option<&'static str>,
) -> Outcome {
    // NCBI's function-static motif positions (`MotifCheck`): one list for the whole run.
    let mut motif_state = MotifState::default();
    // The report is where the stream's posts go (`PostSink for Report`): the stream writes them
    // while it reads an object, rpsbproc between the objects.
    let report: &mut Report = stream.posts();
    report.write(header);
    report.line(&[b"DATA"]);
    let mut session = 0usize;
    let mut queries_done = 0usize;
    loop {
        if stream.posts().stopped() {
            return Outcome::Exit(E_OK);
        }
        let archive = stream.next_archive();
        let report: &mut Report = stream.posts();
        let mut archive = match archive {
            Ok(Some(archive)) => archive,
            Ok(None) => break,
            Err(refusal) => return Outcome::refused(refusal.message),
        };
        if let RequestBody::QueueSearch(search) = &archive.body {
            if let Queries::BioseqSet(entries) = &search.queries {
                // `SEQ_ID_PREFER_ACCESSION_OVER_GI` and `SEQ_ID_AVOID_GI` change the rank of a gi
                // (`ApplicationSettings::seq_id_rank_variable`, Seq_id.cpp:3064-3174), so the
                // label of a query with a gi among its ids (`query_accession`).
                if let Some(name) = seq_id_rank_variable {
                    if entries.iter().any(|entry| {
                        matches!(entry, SeqEntry::Seq(query) if query.ids.iter().any(|id| matches!(id, crate::seqid::model::ParsedSeqId::Gi(_))))
                    }) {
                        return Outcome::refused(format!(
                            "the environment variable {name}, which changes the Seq-id that NCBI's rpsbproc prints for a query with a gi, is not supported by LOSAT's rpsbproc"
                        ));
                    }
                }
            }
        }
        report.progress(&format!("RPSBlast session {session} start..."));
        let RequestBody::QueueSearch(search) = &mut archive.body else {
            report
                .post("Wrong type of request data: CBlast4_request_body::e_Queue_search expected");
            return Outcome::Exit(E_INPUT_DATA_ERROR);
        };
        session += 1;
        let session_text = session.to_string();
        report.write(b"SESSION\t");
        report.write(session_text.as_bytes());
        report.put(b'\t');
        report.write(&search.program);
        if let Some(ident) = &archive.ident {
            report.put(b'\t');
            report.write(ident);
        }
        // Each choice's tab is written before its getter, which may throw (C++17: the left
        // operand of `<<` is evaluated first; NCBI rpsbproc 0.5.1 on an archive whose subject is
        // `sequences`: `SESSION\t1\tblastp\t2.17.0+\t`, then the exception; line RG E1c review 1).
        report.put(b'\t');
        match &search.subject {
            Subject::Database(database) => {
                report.write(database);
            }
            Subject::Other(name) => {
                return process_failure(
                    report,
                    choice_selection("NCBI-Blast4", "Blast4-subject", "database", name),
                );
            }
        }
        if search.has_algorithm_options {
            if search.matrix_name.is_some() {
                report.put(b'\t');
            }
            match &search.matrix_name {
                None => {}
                Some(ParamValue::String(name)) => {
                    report.write(name);
                }
                Some(other) => {
                    return process_failure(
                        report,
                        choice_selection(
                            "NCBI-Blast4",
                            "Blast4-value",
                            "string",
                            value_variant(other),
                        ),
                    );
                }
            }
            if search.evalue_threshold.is_some() {
                report.put(b'\t');
            }
            match &search.evalue_threshold {
                None => {}
                Some(ParamValue::Cutoff(archive::Cutoff::EValue(evalue))) => {
                    report.write(format_g(*evalue).as_bytes());
                }
                Some(ParamValue::Cutoff(archive::Cutoff::RawScore)) => {
                    return process_failure(
                        report,
                        choice_selection("NCBI-Blast4", "Blast4-cutoff", "e-value", "raw-score"),
                    );
                }
                Some(other) => {
                    return process_failure(
                        report,
                        choice_selection(
                            "NCBI-Blast4",
                            "Blast4-value",
                            "cutoff",
                            value_variant(other),
                        ),
                    );
                }
            }
        }
        report.write(b"\n");
        let Queries::BioseqSet(entries) = &mut search.queries else {
            if report.stopped() {
                return Outcome::Exit(E_OK);
            }
            report.post("Wrong type of request data: CBlast4_queries::e_Bioseq_set expected");
            return Outcome::Exit(E_INPUT_DATA_ERROR);
        };
        if let Some(masks) = &archive.results.masks {
            if let Err(error) = check_masks(masks) {
                return process_failure(report, error);
            }
        }
        match &archive.results.alignments {
            None => {
                if report.stopped() {
                    return Outcome::Exit(E_OK);
                }
                report.post(
                    "Wrong type of results data: CBlast4_get_search_results_reply::e_alignments expected",
                );
            }
            Some(aligns) => {
                report.progress("Calculating....");
                let queries = match process_blast_results(entries, aligns, data) {
                    Ok(queries) => queries,
                    Err(error) => return process_failure(report, error),
                };
                for query in &queries {
                    if let Err(error) = print_dom_query(
                        report,
                        query,
                        &session_text,
                        data,
                        mode,
                        show_families,
                        target,
                        &mut motif_state,
                    ) {
                        return print_failure(report, error);
                    }
                    queries_done += 1;
                    if report.stopped() {
                        return Outcome::Exit(E_OK);
                    }
                }
            }
        }
        report.line(&[b"ENDSESSION\t", session_text.as_bytes()]);
        report.progress(&format!("RPSBlast session {session} done..."));
    }
    let report: &mut Report = stream.posts();
    report.line(&[b"ENDDATA"]);
    report.progress(&format!("All done, {queries_done} queries processed\n"));
    Outcome::Exit(E_OK)
}

/// The sessions of an XML input: one per BlastOutput document (`Run`'s XML loop), until one does
/// not read. A document is read whole before its session line (NCBI reads it into one object):
/// what it prints is its version, database, matrix and E-value threshold, and its queries.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:741-778
// (The `DATA` and `ENDDATA` lines and the last line on standard error are those of every input,
// cited at `run_sessions`. C++17 evaluates the left operand of `<<` first: the session line is
// written up to the tab before the matrix when `GetMatrix` throws for a document without
// `Parameters_matrix` (`CUnassignedMember`). `GetParam`, `GetVersion`, `GetDb`, `GetExpect` and
// `GetIterations` get mandatory members, which every document read has. The motif positions of
// `MotifCheck` are the run's one list.)
#[allow(clippy::too_many_arguments)]
fn run_xml_sessions<R: Read>(
    report: &mut Report,
    header: &[u8],
    mut stream: XmlStream<R>,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    evcut: f64,
) -> Outcome {
    let mut motif_state = MotifState::default();
    report.write(header);
    report.line(&[b"DATA"]);
    let mut session = 0usize;
    let mut queries_done = 0usize;
    loop {
        if report.stopped() {
            return Outcome::Exit(E_OK);
        }
        let document = match stream.next_document(evcut, data) {
            Ok(Some(document)) => document,
            Ok(None) => break,
            Err(refusal) => return Outcome::refused(refusal),
        };
        report.progress(&format!("RPSBlast session {session} start..."));
        session += 1;
        let session_text = session.to_string();
        report.write(b"SESSION");
        report.put(b'\t');
        report.write(session_text.as_bytes());
        report.put(b'\t');
        report.write(&document.version);
        report.put(b'\t');
        report.write(&document.db);
        report.put(b'\t');
        let Some(matrix) = &document.matrix else {
            return process_failure(
                report,
                unassigned_member("NCBI-BlastOutput", "Parameters", "matrix"),
            );
        };
        report.write(matrix);
        report.put(b'\t');
        report.write(format_g(document.expect).as_bytes());
        report.write(b"\n");
        report.progress(&format!("RPSBlast session {session} start..."));
        for query in document.queries {
            report.progress("Calculating....");
            let query = match query {
                Ok(query) => query,
                Err(error) => return process_failure(report, error),
            };
            if let Err(error) = print_dom_query(
                report,
                &query,
                &session_text,
                data,
                mode,
                show_families,
                target,
                &mut motif_state,
            ) {
                return print_failure(report, error);
            }
            queries_done += 1;
            if report.stopped() {
                return Outcome::Exit(E_OK);
            }
        }
        report.line(&[b"ENDSESSION\t", session_text.as_bytes()]);
        report.progress(&format!("RPSBlast session {session} done..."));
    }
    report.line(&[b"ENDDATA"]);
    report.progress(&format!("All done, {queries_done} queries processed\n"));
    Outcome::Exit(E_OK)
}

/// How the run ends at an exception of `PrintDomQuery`: a `std::exception` that leaves `Run`
/// (`CNcbiApplication` reports it, ncbiapp.cpp:730-734; probe 3 of PS16b and review oracles R6,
/// R20, R21: the line below with `what()` as bytes, exit 3, the report ends before the site's
/// line), or LOSAT's refusal.
fn print_failure(report: &mut Report, error: PrintError) -> Outcome {
    match error {
        PrintError::Exception(what) => std_exception(report, &what),
        PrintError::Unsupported(message) => {
            if report.stopped() {
                return Outcome::Exit(E_OK);
            }
            Outcome::refused(message)
        }
    }
}

/// How the run ends at an exception of a session (`ProcessError`): NCBI's report of a
/// `CException` or of a `std::exception` that leaves `Run`, or LOSAT's refusal (after what the
/// report has written; nothing more after a failed write).
fn process_failure(report: &mut Report, error: ProcessError) -> Outcome {
    match error {
        ProcessError::Exception(exception) => {
            report.before_post();
            if report.failed.is_some() {
                return Outcome::Exit(E_OK);
            }
            Outcome::Exit(report_exception(&mut report.diag, &exception, false))
        }
        ProcessError::Std(what) => std_exception(report, what.as_bytes()),
        ProcessError::Unsupported(message) => {
            if report.stopped() {
                return Outcome::Exit(E_OK);
            }
            Outcome::refused(message)
        }
    }
}

/// A `std::exception` that leaves `Run`: `CNcbiApplication` reports its `what()` (ncbiapp.cpp:
/// 730-734; PS16b probe 3 and review oracles R6, R20, R21: the line below with `what()` as
/// bytes, exit 3, the report ends where the exception was thrown).
fn std_exception(report: &mut Report, what: &[u8]) -> Outcome {
    report.before_post();
    if report.failed.is_some() {
        return Outcome::Exit(E_OK);
    }
    let mut message = b"Application's execution failed: ".to_vec();
    message.extend_from_slice(what);
    report
        .diag
        .post_bytes(&err_post_x_bytes(CORELIB_APP, 7, &message));
    Outcome::Exit(3)
}

/// `PrintDomQuery` of a protein query with `-t doms`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:947-1031
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1801-1802
// `m_uiClusterPssmId` is never set by the offline data (`TDomain` keeps its 0): the key of a
// superfamily is its region times `OVERLAPLEADING`, an `unsigned int` product that wraps, and
// `emplace` keeps the first record of a key. `Multidom` (`m_iRepClass` 1) is never made
// (`domains::calculate`), so it is never printed.
#[allow(clippy::too_many_arguments)]
fn print_dom_query(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let length = query.length.to_string();
    let kind: &[u8] = if query.is_na {
        b"\tNucleotide\t"
    } else {
        b"\tPeptide\t"
    };
    report.line(&[
        b"QUERY\t",
        &query.accession,
        kind,
        length.as_bytes(),
        b"\t",
        &query.title,
    ]);
    if query.is_na {
        nucleotide::print_frames(
            report,
            query,
            session,
            data,
            mode,
            show_families,
            target,
            motif_state,
        )?;
    } else {
        print_protein(
            report,
            query,
            session,
            data,
            mode,
            show_families,
            target,
            motif_state,
        )?;
    }
    report.line(&[b"ENDQUERY"]);
    if report.quiet {
        return Ok(());
    }
    report.before_cerr();
    if report.failed.is_some() {
        return Ok(());
    }
    let mut stderr = std::io::stderr().lock();
    let _ = stderr.write_all(b"End processing Query ");
    let _ = stderr.write_all(&query.title);
    let _ = stderr.write_all(b"\n");
    Ok(())
}

/// The DOMAINS, SUPERFAMILIES, SITES and MOTIFS parts of a protein query.
#[allow(clippy::too_many_arguments)]
fn print_protein(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let records = if target == TargetData::Feats {
        Vec::new()
    } else {
        domain_records(&query.frames[0], &query.aligns, data, mode)
    };
    if !records.is_empty() {
        let mut families: std::collections::BTreeMap<u32, Record> = Default::default();
        report.line(&[b"DOMAINS"]);
        for record in &records {
            let align = &query.aligns[record.align];
            report.write(session.as_bytes());
            report.write(b"\t");
            report.write(&query.accession);
            report.write(b"\t");
            let key = (align.region as u32).wrapping_mul(OVERLAPLEADING);
            match mode {
                DataMode::Rep => {
                    if align.spec_qualified {
                        print_dom_line(report, b"Specific", record, query, data);
                        if record.cluster.is_some() {
                            families.entry(key).or_insert(*record);
                        }
                    } else if record.cluster.is_some() {
                        print_cluster_line(report, record, query, data);
                    } else {
                        print_dom_line(report, b"Non-specific", record, query, data);
                    }
                }
                DataMode::Std | DataMode::Full => {
                    let kind: &[u8] = if align.rep && align.spec_qualified {
                        b"Specific"
                    } else {
                        b"Non-specific"
                    };
                    print_dom_line(report, kind, record, query, data);
                    if record.cluster.is_some() {
                        families.entry(key).or_insert(*record);
                    }
                }
            }
            report.write(b"\n");
        }
        report.line(&[b"ENDDOMAINS"]);
        if show_families {
            report.line(&[b"SUPERFAMILIES"]);
            for record in families.values() {
                report.write(session.as_bytes());
                report.write(b"\t");
                report.write(&query.accession);
                report.write(b"\t");
                print_cluster_line(report, record, query, data);
                report.write(b"\n");
            }
            report.line(&[b"ENDSUPERFAMILIES"]);
        }
    }
    if target != TargetData::Doms {
        print_features(report, query, session, data, mode, motif_state)?;
    }
    Ok(())
}

/// Why the report of a query stops.
enum PrintError {
    /// The `what()` of a `std::exception` that ends NCBI's rpsbproc: a motif that its parser
    /// refuses (`THROW_SIMPLE`), or a minimal map that its string cannot hold.
    Exception(Vec<u8>),
    /// A case LOSAT does not reproduce (the message may quote data bytes).
    Unsupported(Vec<u8>),
}

/// How the motif check of `site` stops the report. The refusals name the site and quote its
/// motif's bytes as the data file has them, as NCBI's parse error quotes the motif.
fn motif_error(error: MotifError, data: &CddData, site: &data::Site) -> PrintError {
    let refusal = |what: &str| {
        let mut message = what.as_bytes().to_vec();
        message.extend_from_slice(b" is not supported by LOSAT's rpsbproc: site '");
        message.extend_from_slice(data.text(site.title));
        message.extend_from_slice(b"', motif '");
        message.extend_from_slice(data.text(site.motif));
        message.push(b'\'');
        PrintError::Unsupported(message)
    };
    match error {
        MotifError::Parse { motif, position } => {
            // `CSimpleException`'s text: the motif's bytes as the data file has them.
            let mut what = b"Motif string parse error -- Motif = ".to_vec();
            what.extend_from_slice(&motif);
            what.extend_from_slice(format!(", error position: {position}").as_bytes());
            PrintError::Exception(what)
        }
        MotifError::Map(MapEnd::Throw(what)) => PrintError::Exception(what.as_bytes().to_vec()),
        MotifError::Map(MapEnd::HostDependent(bytes)) => refusal(&format!(
            "a motif whose minimal map makes NCBI's rpsbproc hold {bytes} bytes at once (1 GiB or more in allocations below 2^47 bytes, where NCBI's result depends on the host's memory)"
        )),
        MotifError::Match(MatchError::ReadPast) => refusal(
            "a motif match that reads beyond the residues of a site where NCBI's result depends on the bytes it reads there (or on its memory for some of them), or that may read beyond NCBI's buffer of them (undefined),",
        ),
        MotifError::Match(MatchError::Exhausted) => refusal(
            "a motif match for which NCBI's rpsbproc holds 1 GiB or more at once (its match records with their copies on its stack and the string of the residues; a repeat that does not advance, such as `<*`, grows the records until the memory runs out; NCBI's result depends on the host's memory)",
        ),
    }
}

/// `toupper` of `q.m_strSeqData[res]`, cast to `char`: the residue at `res`.
/// `std::string::operator[]` gives the terminating 0 at its size; beyond that, and before 0, NCBI
/// reads bytes outside the string (undefined), which LOSAT refuses.
fn upper_residue(query: &DomQuery, res: i64) -> Result<u8, PrintError> {
    upper_at(&query.residues, res)
}

/// `toupper` of `text[res]`, cast to `char`, for a `std::string` (the query's residues or the
/// translation of one of its frames).
fn upper_at(text: &[u8], res: i64) -> Result<u8, PrintError> {
    let len = text.len() as i64;
    if (0..len).contains(&res) {
        Ok(text[res as usize].to_ascii_uppercase())
    } else if res == len {
        Ok(0)
    } else {
        Err(PrintError::Unsupported(
            format!(
                "a site or motif at position {} of a query with {len} residues (NCBI reads outside the sequence) is not supported by LOSAT's rpsbproc",
                printed_position(res)
            )
            .into_bytes(),
        ))
    }
}

/// A site of a structural motif type, collected for the MOTIFS part
/// (`__MotifType_base`: the site, the record, the source PSSM-ID).
struct MotifRecord<'a> {
    site: &'a data::Site,
    record: usize,
    source: u32,
}

/// The source column: the PSSM-ID, or the accession of a model with a PSSM-ID of
/// `INVALIDPSSMID` or more that the data has.
fn write_source(report: &mut Report, data: &CddData, source: u32) {
    if source < data::INVALID_PSSM_ID {
        report.write(source.to_string().as_bytes());
    } else if let Some(domain) = data.find_domain(source) {
        report.write(data.text(data.domains[domain].accession));
    } else {
        report.write(source.to_string().as_bytes());
    }
}

/// A site's positions as a set (`featsegs` is a `CSegSet` made from `*iterFeat`: `ori_from` is `from`).
fn site_segset(data: &CddData, site: &data::Site) -> SegSet {
    let mut set = SegSet::new();
    for seg in data.site_segs(site) {
        set.add(i64::from(seg[0]), i64::from(seg[1]));
    }
    set
}

/// Whether `featsegs.GetTotalResidues()` divided by `GetCompleteSize()` (as a `double`) is above
/// 0.8, for a mapped site.
fn mostly_mapped(featsegs: &SegSet, site: &data::Site) -> bool {
    !featsegs.is_empty() && (featsegs.total_residues() as f64 / f64::from(site.complete_size)) > 0.8
}

/// The SITES and MOTIFS parts of a protein query (`PrintDomQuery`, `-t feats` and `both`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1033-1266
// (An `sd` record first in the list prints no SITES part and no motif is collected then.)
fn print_features(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let records: Vec<FeatureRecord> = feature_records(&query.frames[0], &query.aligns, data, mode);
    let Some(first) = records.first() else {
        return Ok(());
    };
    let mut motifs: Vec<MotifRecord> = Vec::new();
    if !data.domains[first.domain].is_struct_dom {
        report.line(&[b"SITES"]);
        let mut at = 0usize;
        while at < records.len() {
            let record = &records[at];
            if data.domains[record.domain].is_struct_dom {
                break;
            }
            let align = &query.aligns[record.align];
            let is_specific = align.rep && align.spec_qualified;
            let (kind, specific, source): (&[u8], bool, u32) =
                if !is_specific && record.root.is_some() {
                    (
                        b"Generic",
                        false,
                        data.domains[record.domain].hierarchy_root,
                    )
                } else {
                    (b"Specific", true, align.pssm_id)
                };
            for site in data.sites(record.domain, specific) {
                if site.site_type == data::TYPE_STRUCT_MOTIF {
                    motifs.push(MotifRecord {
                        site,
                        record: at,
                        source,
                    });
                    continue;
                }
                print_site(
                    report,
                    query,
                    session,
                    data,
                    align,
                    site,
                    kind,
                    Some(source),
                    motif_state,
                )?;
            }
            at += 1;
        }
        // Sites of the structural domains: source 0.
        while at < records.len() {
            let record = &records[at];
            let align = &query.aligns[record.align];
            for site in data.sites(record.domain, true) {
                if site.site_type == data::TYPE_STRUCT_MOTIF {
                    motifs.push(MotifRecord {
                        site,
                        record: at,
                        source: 0,
                    });
                    continue;
                }
                print_site(
                    report,
                    query,
                    session,
                    data,
                    align,
                    site,
                    b"Specific",
                    None,
                    motif_state,
                )?;
            }
            at += 1;
        }
        report.line(&[b"ENDSITES"]);
    }
    if !motifs.is_empty() {
        report.line(&[b"MOTIFS"]);
        // The first loop stops at the first `sd` record; the second takes the rest.
        let mut structural = false;
        for motif in &motifs {
            let record = &records[motif.record];
            structural = structural || data.domains[record.domain].is_struct_dom;
            let align = &query.aligns[record.align];
            let mut featsegs = site_segset(data, motif.site);
            if structural && !align.clip_set.is_empty() {
                featsegs.cross(&align.clip_set);
            }
            map_seg_set(&mut featsegs, align);
            if mostly_mapped(&featsegs, motif.site) {
                let (res0, res1) = (featsegs.left(), featsegs.right());
                // The residues first: a refusal leaves no part of the line.
                let (letter0, letter1) = (upper_residue(query, res0)?, upper_residue(query, res1)?);
                report.write(session.as_bytes());
                report.write(b"\t");
                report.write(&query.accession);
                report.write(b"\t");
                report.write(data.text(motif.site.title));
                report.write(b"\t");
                report.write(&[letter0]);
                report.write(printed_position(res0).to_string().as_bytes());
                report.write(b"\t");
                report.write(printed_position(res1).to_string().as_bytes());
                report.write(&[letter1]);
                report.write(b"\t");
                write_source(report, data, motif.source);
                report.write(b"\n");
            }
        }
        report.line(&[b"ENDMOTIFS"]);
    }
    Ok(())
}

/// One site line of the SITES part. `source`: the PSSM-ID column (`None`: the structural
/// domains' `0`).
#[allow(clippy::too_many_arguments)]
fn print_site(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    align: &domains::DomAlign,
    site: &data::Site,
    kind: &[u8],
    source: Option<u32>,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let mut featsegs = site_segset(data, site);
    map_seg_set(&mut featsegs, align);
    if !mostly_mapped(&featsegs, site) {
        return Ok(());
    }
    let mapped_positions = translated_pos_map(&featsegs, align, query.length)
        .map_err(|what| PrintError::Exception(what.into_bytes()))?;
    if motif_state
        .motif_check(data, site, &mapped_positions, &query.residues)
        .map_err(|error| motif_error(error, data, site))?
        > 0
    {
        return Ok(());
    }
    // The residues first: a refusal leaves no part of the line.
    let mut residues = Vec::new();
    for seg in featsegs.segs() {
        for res in seg.from..=seg.to {
            residues.push((upper_residue(query, res)?, res));
        }
    }
    report.write(session.as_bytes());
    report.write(b"\t");
    report.write(&query.accession);
    report.write(b"\t");
    report.write(kind);
    report.write(b"\t");
    report.write(data.text(site.title));
    report.write(b"\t");
    let mapped = residues.len();
    let mut delimit: &[u8] = b"";
    for (letter, res) in residues {
        report.write(delimit);
        report.write(&[letter]);
        report.write(printed_position(res).to_string().as_bytes());
        delimit = b",";
    }
    let complete: i64 = data
        .site_segs(site)
        .iter()
        .map(|seg| i64::from(seg[1]) - i64::from(seg[0]) + 1)
        .sum();
    report.write(format!("\t{complete}\t{mapped}\t").as_bytes());
    match source {
        Some(source) => write_source(report, data, source),
        None => report.write(b"0"),
    }
    report.write(b"\n");
    Ok(())
}

/// `N`, `C`, both or `-`: more than 20 % of the model missing at that end.
fn print_incomplete(report: &mut Report, n_missing: f64, c_missing: f64) {
    let mut missing = false;
    if n_missing > 0.2 {
        report.write(b"N");
        missing = true;
    }
    if c_missing > 0.2 {
        report.write(b"C");
        missing = true;
    }
    if !missing {
        report.write(b"-");
    }
}

/// `PrintDomLine`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1807-1832
// (`m_iFrom` plus `COORDSBASE` is `int` arithmetic, which overflows at `INT_MAX` (undefined in C++):
// `printed_position`, for protein and nucleotide queries alike.)
fn print_dom_line(
    report: &mut Report,
    kind: &[u8],
    record: &Record,
    query: &DomQuery,
    data: &CddData,
) {
    let align = &query.aligns[record.align];
    let domain = &data.domains[record.domain];
    let pssm_id = if align.pssm_id >= data::INVALID_PSSM_ID {
        0
    } else {
        align.pssm_id
    };
    let numbers = format!(
        "\t{pssm_id}\t{}\t{}\t{}\t{}\t",
        printed_position(align.from),
        printed_position(align.to),
        format_g(align.e_value),
        format_g(align.bit_score)
    );
    report.write(kind);
    report.write(numbers.as_bytes());
    report.write(data.text(domain.accession));
    report.write(b"\t");
    report.write(data.text(domain.short_name));
    report.write(b"\t");
    print_incomplete(report, align.n_missing, align.c_missing);
    report.write(b"\t");
    match record.cluster {
        None => report.write(b"-"),
        Some(cluster) => {
            let cluster = &data.clusters[cluster];
            if cluster.pssm_id >= data::INVALID_PSSM_ID {
                report.write(data.text(cluster.accession));
            } else {
                report.write(cluster.pssm_id.to_string().as_bytes());
            }
        }
    }
}

/// `PrintClstLine`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1834-1851
// (`m_uiClusterPssmId` is 0: the offline loader never sets it.)
fn print_cluster_line(report: &mut Report, record: &Record, query: &DomQuery, data: &CddData) {
    let align = &query.aligns[record.align];
    let cluster = &data.clusters[record.cluster.expect("a superfamily line has its cluster")];
    let numbers = format!(
        "Superfamily\t0\t{}\t{}\t{}\t{}\t",
        printed_position(align.from),
        printed_position(align.to),
        format_g(align.e_value),
        format_g(align.bit_score)
    );
    report.write(numbers.as_bytes());
    report.write(data.text(cluster.accession));
    report.write(b"\t");
    report.write(data.text(cluster.short_name));
    report.write(b"\t");
    print_incomplete(report, align.n_missing, align.c_missing);
    report.write(b"\t-");
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;

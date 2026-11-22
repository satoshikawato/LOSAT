//! The streams of the file arguments of NCBI's argument reader (`CArg_InputFile`,
//! `CArg_OutputFile`, c++/src/corelib/ncbiargs.cpp) and the lines that NCBI's line reader
//! (`NcbiGetlineEOL`, c++/src/corelib/ncbistre.cpp) reads from them, on each platform NCBI
//! builds for.

use std::fs::File;
use std::io::{Read, Write};

use super::args::{arg_message, ArgException};
use crate::ncbi_exception::{NcbiException, ThrowSite};

/// The exception of a file argument whose stream does not open.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-621
fn not_accessible(name: &str, value: &str) -> ArgException {
    ArgException::new(
        "eNoFile",
        arg_message(name, "File is not accessible", value),
    )
}

/// The place of `CArg_Ios::x_Open`'s throw (the `NCBI_THROW` of ncbiargs.cpp:618), which
/// `e.what()` prints without a function: segmasker 2.17.0 catches the `CArgException` of
/// `AsInputFile` and `AsOutputFile` in its `Run` and writes `e.what()` (fixtures `noin_file`,
/// `out_bad`).
const X_OPEN: ThrowSite = ThrowSite {
    file: "src/corelib/ncbiargs.cpp",
    line: 618,
    kind: "CArgException::eNoFile",
    function: None,
};

/// The exception that `error` is when a program catches it (`e.what()`): the `eNoFile` of a
/// file argument that does not open (`open_input_file`, `open_output_file`); `None` for the
/// others, which those functions do not raise.
pub(crate) fn file_exception(error: &ArgException) -> Option<NcbiException> {
    (error.code == "eNoFile" && error.messages.len() == 1)
        .then(|| NcbiException::new(X_OPEN, error.message()))
}

/// The platforms on which NCBI reads the lines of an input file differently: the text mode
/// of the stream (`InputFile`) and `NcbiGetlineEOL` (`input_lines`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LineEnds {
    /// Linux and the other systems of NCBI's `#else` (LOSAT's WASI build too).
    Unix,
    /// `NCBI_OS_MSWIN`.
    Windows,
    /// `NCBI_OS_DARWIN` (macOS).
    Darwin,
}

impl LineEnds {
    /// The platform LOSAT is built for.
    pub(crate) const HOST: LineEnds = if cfg!(windows) {
        LineEnds::Windows
    } else if cfg!(target_os = "macos") {
        LineEnds::Darwin
    } else {
        LineEnds::Unix
    };
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:580-585
/// The text mode of an input stream opened without `fBinary` on `platform`, as a filter of the
/// file's bytes in the order they are read: on Windows (the C runtime's text mode, under
/// NCBI's `CNcbiIfstream` and the `cin` that `NcbiSys_setmode` sets to `O_TEXT`,
/// ncbiargs.cpp:717-720) the bytes before the first Ctrl-Z, with each `\r\n` read as `\n`
/// (probe O of line BDC: O04, O09-O12); elsewhere the bytes of the file.
#[derive(Debug)]
struct TextMode {
    platform: LineEnds,
    /// A `\r` read last, which is kept unless a `\n` follows.
    carriage_return: bool,
    /// A Ctrl-Z was read: the stream reads nothing more.
    ended: bool,
}

impl TextMode {
    fn new(platform: LineEnds) -> Self {
        TextMode {
            platform,
            carriage_return: false,
            ended: false,
        }
    }

    /// Appends to `out` what the stream reads of `input`, the next bytes of the file.
    fn push(&mut self, input: &[u8], out: &mut Vec<u8>) {
        if self.platform != LineEnds::Windows {
            out.extend_from_slice(input);
            return;
        }
        for &byte in input {
            if self.ended {
                return;
            }
            if self.carriage_return {
                self.carriage_return = false;
                if byte != b'\n' {
                    out.push(b'\r');
                }
            }
            match byte {
                0x1a => self.ended = true,
                b'\r' => self.carriage_return = true,
                _ => out.push(byte),
            }
        }
    }

    /// The end of the file: a last `\r` is read.
    fn finish(&mut self, out: &mut Vec<u8>) {
        if self.carriage_return {
            self.carriage_return = false;
            out.push(b'\r');
        }
    }
}

/// A stream of `source` read in the text mode of `platform` (`TextMode`), which reads nothing
/// more once it has ended (the `eofbit` of NCBI's stream stays set: a later read gets nothing).
struct TextReader<R> {
    source: R,
    text: TextMode,
    /// Bytes of the text mode not handed out yet (Windows only), from `at`.
    pending: Vec<u8>,
    at: usize,
    /// The source has ended.
    ended: bool,
}

impl<R: Read> TextReader<R> {
    fn new(source: R, platform: LineEnds) -> Self {
        TextReader {
            source,
            text: TextMode::new(platform),
            pending: Vec::new(),
            at: 0,
            ended: false,
        }
    }
}

impl<R: Read> Read for TextReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.text.platform != LineEnds::Windows {
            if self.ended || buf.is_empty() {
                return Ok(0);
            }
            let read = self.source.read(buf)?;
            self.ended = read == 0;
            return Ok(read);
        }
        let mut raw = [0u8; 8192];
        while self.at == self.pending.len() {
            self.pending.clear();
            self.at = 0;
            if self.text.ended || self.ended {
                return Ok(0);
            }
            let read = self.source.read(&mut raw)?;
            if read == 0 {
                self.ended = true;
                self.text.finish(&mut self.pending);
            } else {
                self.text.push(&raw[..read], &mut self.pending);
            }
        }
        let count = buf.len().min(self.pending.len() - self.at);
        buf[..count].copy_from_slice(&self.pending[self.at..self.at + count]);
        self.at += count;
        Ok(count)
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-720 (`-` is `cin`, set to the
// text mode on Windows: `NcbiSys_setmode` on the file number of stdin (`NcbiSys_fileno`), with
// `O_TEXT`)
/// The standard input of a run (`cin`): the process's on the command line, the caller's reader
/// in the library API (`api::run_dustmasker`, `api::run_blastdbcmd`), read in the text mode of
/// LOSAT's platform. It is one stream for the whole run, as `cin` is: a second argument `-`
/// reads on where the first stopped, and nothing once the stream has ended (a terminal's next
/// lines, the bytes after a Ctrl-Z on Windows).
pub(crate) struct Cin<'a>(TextReader<&'a mut dyn Read>);

impl<'a> Cin<'a> {
    pub(crate) fn new(input: &'a mut dyn Read) -> Self {
        Cin(TextReader::new(input, LineEnds::HOST))
    }
}

impl Read for Cin<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.0.read(buf)
    }
}

/// Where an input stream reads from.
enum Source<'a> {
    /// The standard input (`Cin`, in its text mode already).
    Stdin(&'a mut dyn Read),
    File(TextReader<File>),
}

/// The stream of an `eInputFile` argument (`AsInputFile()`, without `fBinary`): the standard
/// input or the file, read in the text mode of LOSAT's platform. It is not buffered (wrap it in a
/// `std::io::BufReader` to read lines); a read error is the operating system's (NCBI's
/// libstdc++ stream throws its own text for it, which the program that reads writes).
pub(crate) struct InputFile<'a> {
    source: Source<'a>,
}

impl Read for InputFile<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match &mut self.source {
            Source::Stdin(stdin) => stdin.read(buf),
            Source::File(file) => file.read(buf),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-737
/// `AsInputFile()` of the argument `name` whose value is `value`: `-` is the standard input,
/// `stdin` (probes K30-K32 of line BDC), the empty name opens nothing (K39), and another name
/// opens the file; on Linux a directory opens too (its reads fail). `Err` is NCBI's `eNoFile`.
pub(crate) fn open_input_file<'a>(
    name: &str,
    value: &str,
    stdin: &'a mut Cin<'_>,
) -> Result<InputFile<'a>, ArgException> {
    let source = if value == "-" {
        Source::Stdin(stdin)
    } else if value.is_empty() {
        return Err(not_accessible(name, value));
    } else {
        let file = File::open(value).map_err(|_| not_accessible(name, value))?;
        Source::File(TextReader::new(file, LineEnds::HOST))
    };
    Ok(InputFile { source })
}

/// The file `path` opened as a program's own `CNcbiIfstream` of `path` opens it (not an argument's
/// stream): read in the text mode of LOSAT's platform; `Err` when it does not open (on Linux a
/// directory opens, and its reads fail).
pub(crate) fn open_text_file(path: &str) -> std::io::Result<Box<dyn Read>> {
    Ok(Box::new(TextReader::new(File::open(path)?, LineEnds::HOST)))
}

/// The bytes of the input file `value` of the argument `name` (`AsInputFile`, `-` is `stdin`)
/// as the loops of `NcbiGetlineEOL` over its stream read them (blastdbcmd's batch files,
/// `input_lines`): `Err` with NCBI's message (`File is not accessible`) when it does not open or
/// a read fails. On Linux the stream of a directory opens and reads nothing there (probes
/// K33-K35 of line BDC).
pub(crate) fn read_input_file(
    name: &str,
    value: &str,
    stdin: &mut Cin<'_>,
) -> Result<Vec<u8>, String> {
    if cfg!(unix) && value != "-" && !value.is_empty() && std::path::Path::new(value).is_dir() {
        return Ok(Vec::new());
    }
    let mut file =
        open_input_file(name, value, stdin).map_err(|error| error.message().to_string())?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .map_err(|_| not_accessible(name, value).message().to_string())?;
    Ok(bytes)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:196-208
// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:90-104
// The loops of blastdbcmd that read lines so skip the empty ones (a test that the line is not
// empty, blastdbcmd.cpp:333, 385, 450, 497): on macOS a line then ends at any `\r` or `\n`, and the
// `\r\n` that NCBI takes for one end leaves no other line.
/// The lines of an input file (`bytes` as `read_input_file` gives them) that the loops of
/// `NcbiGetlineEOL` read on `platform`, without the empty ones.
pub(crate) fn input_lines(bytes: &[u8], platform: LineEnds) -> Vec<&[u8]> {
    let filled = |line: &&[u8]| !line.is_empty();
    match platform {
        LineEnds::Unix => bytes.split(|&byte| byte == b'\n').filter(filled).collect(),
        LineEnds::Windows => bytes
            .split(|&byte| byte == b'\n')
            .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
            .filter(filled)
            .collect(),
        LineEnds::Darwin => bytes
            .split(|&byte| byte == b'\r' || byte == b'\n')
            .filter(filled)
            .collect(),
    }
}

/// The stream of an `eOutputFile` argument (`AsOutputFile()`): the standard output, or the file,
/// made or emptied. On Windows NCBI's stream of a text output (opened without `fBinary`) writes
/// each `\n` as `\r\n`; LOSAT's programs write their bytes as they are on every platform.
pub(crate) enum OutputFile<'a> {
    Stdout(&'a mut dyn Write),
    File(File),
}

impl Write for OutputFile<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputFile::Stdout(out) => out.write(buf),
            OutputFile::File(file) => file.write(buf),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        match self {
            OutputFile::Stdout(out) => out.write_all(buf),
            OutputFile::File(file) => file.write_all(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputFile::Stdout(out) => out.flush(),
            OutputFile::File(file) => file.flush(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:777-799 (`CArg_OutputFile` opens a
// `CNcbiOfstream`, or takes `cout`, without exceptions enabled on it: a failed write sets
// `badbit` and later writes do nothing)
/// An output as NCBI's `ostream` writes to it: after a write fails (`badbit`), the later writes
/// do nothing and nothing reports the failure, so the run goes on to its end (NCBI 2.17.0 with
/// `-out /dev/full`: dustmasker's `interval`, `acclist` and `fasta` writers, blastdbcmd's
/// retrieval, which still posts the errors of later entries). [`Ostream::bad`] tells a writer
/// that NCBI's serial stream would report the failure (`CObjectOStream`), which LOSAT refuses.
///
/// A write to a pipe whose reader has closed ends NCBI's process by SIGPIPE outside Windows
/// (exit status 141, no message; PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10). LOSAT's
/// runtime ignores the signal, so [`Ostream::closed`] is set there instead and every later write
/// returns that error: the program ends at once, with exit status 0 and no message
/// ([`ends_like_sigpipe`]). On Windows NCBI's write fails without a signal and is a failed
/// write like the others.
pub(crate) struct Ostream<'a> {
    inner: &'a mut dyn Write,
    bad: bool,
    closed: bool,
}

/// Whether the failed write `error` is one that ends NCBI's process by SIGPIPE: a closed pipe,
/// outside Windows.
pub(crate) fn ends_like_sigpipe(error: &std::io::Error) -> bool {
    cfg!(not(windows)) && error.kind() == std::io::ErrorKind::BrokenPipe
}

impl<'a> Ostream<'a> {
    pub(crate) fn new(inner: &'a mut dyn Write) -> Self {
        Ostream {
            inner,
            bad: false,
            closed: false,
        }
    }

    /// Whether a write or a flush has failed (`badbit`).
    pub(crate) fn bad(&self) -> bool {
        self.bad
    }

    /// Whether a write found the pipe closed where that ends NCBI's process (SIGPIPE).
    pub(crate) fn closed(&self) -> bool {
        self.closed
    }

    /// The outcome of a write or a flush of the inner stream.
    fn result(&mut self, outcome: std::io::Result<()>) -> std::io::Result<()> {
        if let Err(error) = outcome {
            self.bad = true;
            if ends_like_sigpipe(&error) {
                self.closed = true;
                return Err(error);
            }
        }
        Ok(())
    }

    fn closed_pipe() -> std::io::Error {
        std::io::Error::from(std::io::ErrorKind::BrokenPipe)
    }
}

/// `cerr`'s tie to `cout` for an output buffered over an `Ostream`: before a diagnostic line
/// NCBI's stream flushes the output (to the file or the pipe), so the two come out in the order
/// they were written. Flushes `sink`; whether the pipe is closed, which ends the program
/// before the line (SIGPIPE ends NCBI's process there).
pub(crate) fn flush_before_diagnostic(sink: &mut std::io::BufWriter<&mut Ostream<'_>>) -> bool {
    let _ = sink.flush();
    sink.get_ref().closed()
}

impl Write for Ostream<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if self.closed {
            return Err(Self::closed_pipe());
        }
        if !self.bad {
            let outcome = self.inner.write_all(buf);
            self.result(outcome)?;
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if self.closed {
            return Err(Self::closed_pipe());
        }
        if !self.bad {
            let outcome = self.inner.flush();
            self.result(outcome)?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:777-799
/// `AsOutputFile()` of the argument `name` whose value is `value`: `-` is `stdout`, the empty
/// name opens nothing, and another name makes or empties the file (`IOS_BASE::out`, without
/// `fCreatePath`). `Err` is NCBI's `eNoFile`.
pub(crate) fn open_output_file<'a>(
    name: &str,
    value: &str,
    stdout: &'a mut dyn Write,
) -> Result<OutputFile<'a>, ArgException> {
    if value == "-" {
        return Ok(OutputFile::Stdout(stdout));
    }
    if value.is_empty() {
        return Err(not_accessible(name, value));
    }
    File::create(value)
        .map(OutputFile::File)
        .map_err(|_| not_accessible(name, value))
}

#[cfg(test)]
#[path = "files_tests.rs"]
mod tests;

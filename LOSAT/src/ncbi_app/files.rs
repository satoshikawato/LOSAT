//! The streams of the file arguments of NCBI's argument reader (`CArg_InputFile`,
//! `CArg_OutputFile`, c++/src/corelib/ncbiargs.cpp) and the lines that NCBI's line reader
//! (`NcbiGetlineEOL`, c++/src/corelib/ncbistre.cpp) reads from them, on each platform NCBI
//! builds for.

use std::fs::File;
use std::io::{Read, Write};

use super::args::{arg_message, ArgException};

/// The exception of a file argument whose stream does not open.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-621
/// ```c++
/// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
/// {
///     if ( !m_Ios ) {
///         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
///             "File is not accessible",AsString()));
///     }
/// }
/// ```
fn not_accessible(name: &str, value: &str) -> ArgException {
    ArgException::new(
        "eNoFile",
        arg_message(name, "File is not accessible", value),
    )
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
// ```c++
// IOS_BASE::openmode CArg_Ios::IosMode(CArgValue::TFileFlags flags)
// {
//     IOS_BASE::openmode openmode = (IOS_BASE::openmode) 0;
//     if (flags & CArgValue::fBinary) {
//         openmode |= IOS_BASE::binary;
//     }
// ```
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
// text mode on Windows: `NcbiSys_setmode(NcbiSys_fileno(stdin), ... O_TEXT)`)
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
// ```c++
//     if (AsString() == "-") {
// #if defined(NCBI_OS_MSWIN)
//         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
// #endif
//         m_Ios  = &cin;
//     } else if ( !AsString().empty() ) {
//         if (!fstrm) {
//             fstrm = new CNcbiIfstream;
//         }
//         if (fstrm) {
//             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
//                 fstrm = NULL;
//             } else {
//                 m_DeleteFlag = true;
//             }
//         }
//         m_Ios = fstrm;
//     }
//     CArg_Ios::x_Open(flags);
// ```
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
// ```c++
// CNcbiIstream& NcbiGetlineEOL(CNcbiIstream& is, string& str, SIZE_TYPE* count)
// {
// #if   defined(NCBI_OS_MSWIN)
//     NcbiGetline(is, str, '\n', count);
//     if (!str.empty()  &&  str[str.length() - 1] == '\r')
//         str.resize(str.length() - 1);
// #elif defined(NCBI_OS_DARWIN)
//     NcbiGetline(is, str, "\r\n", count);
// #else /* assume UNIX-like EOLs */
//     NcbiGetline(is, str, '\n', count);
// #endif //NCBI_OS_...
//     return is;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:90-104
// ```c++
//         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
//         if (delim_pos != NPOS) {
//             // Special case -- if two different delimiters are back to
//             // back and in the same order as in delims, treat them as
//             // a single delimiter (necessary for correct handling of
//             // DOS/MAC-style CR/LF endings).
//             ch = is.rdbuf()->sgetc();
//             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
//                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
//                 is.rdbuf()->sbumpc();
//                 delim_count = 2;
//             } else {
//                 delim_count = 1;
//             }
//             break;
//         }
// ```
// The loops of blastdbcmd that read lines so skip the empty ones (`if ( !line.empty() )`,
// blastdbcmd.cpp:333, 385, 450, 497): on macOS a line then ends at any `\r` or `\n`, and the
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
// ```c++
//     if (AsString() == "-") {
// #if defined(NCBI_OS_MSWIN)
//         NcbiSys_setmode(NcbiSys_fileno(stdout), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
// #endif
//         m_Ios = &cout;
//     } else if ( !AsString().empty() ) {
//         if (!fstrm) {
//             fstrm = new CNcbiOfstream;
//         }
//         if (fstrm) {
//             if (x_CreatePath(m_CurrentFlags)) {
//                 fstrm->open(AsString().c_str(),IOS_BASE::out | mode);
//             }
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
//                 fstrm = NULL;
//             } else {
//                 m_DeleteFlag = true;
//             }
//         }
//         m_Ios = fstrm;
//     }
//     CArg_Ios::x_Open(flags);
// ```
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
mod tests {
    use super::*;

    /// The bytes that the stream of `bytes` reads on `platform`, fed in pieces of `piece` bytes.
    fn text_mode(bytes: &[u8], platform: LineEnds, piece: usize) -> Vec<u8> {
        let mut text = TextMode::new(platform);
        let mut out = Vec::new();
        for chunk in bytes.chunks(piece.max(1)) {
            text.push(chunk, &mut out);
        }
        text.finish(&mut out);
        out
    }

    /// The lines that `NcbiGetlineEOL` reads from `bytes` on `platform`, after the text mode.
    fn lines(bytes: &str, platform: LineEnds) -> Vec<String> {
        let read = text_mode(bytes.as_bytes(), platform, bytes.len());
        // The text mode reads the same bytes whatever pieces the file is read in.
        for piece in 1..4 {
            assert_eq!(text_mode(bytes.as_bytes(), platform, piece), read);
        }
        input_lines(&read, platform)
            .into_iter()
            .map(|line| String::from_utf8(line.to_vec()).unwrap())
            .collect()
    }

    #[test]
    fn batch_lines_end_as_on_each_platform_of_ncbi() {
        use LineEnds::{Darwin, Unix, Windows};
        // The inputs of probe O of line BDC (BDC5) and the lines NCBI 2.17.0 read from them.
        let cases: &[(&str, [&[&str]; 3])] = &[
            ("a\nb\n", [&["a", "b"], &["a", "b"], &["a", "b"]]),
            ("a\r\nb\r\n", [&["a\r", "b\r"], &["a", "b"], &["a", "b"]]),
            ("a\rb\r", [&["a\rb\r"], &["a\rb"], &["a", "b"]]),
            (
                "a\r\r\nb\r\n",
                [&["a\r\r", "b\r"], &["a", "b"], &["a", "b"]],
            ),
            ("a\n\rb\n", [&["a", "\rb"], &["a", "\rb"], &["a", "b"]]),
            ("a\r\nb\r", [&["a\r", "b\r"], &["a", "b"], &["a", "b"]]),
            (
                "\r\n \r\nb\n",
                [&["\r", " \r", "b"], &[" ", "b"], &[" ", "b"]],
            ),
            ("a\n\x1ab\n", [&["a", "\x1ab"], &["a"], &["a", "\x1ab"]]),
            (
                "a\r\x1a\nb\n",
                [&["a\r\x1a", "b"], &["a"], &["a", "\x1a", "b"]],
            ),
            ("a\r\rb\r\n", [&["a\r\rb\r"], &["a\r\rb"], &["a", "b"]]),
            (
                "a\r\n\rb\r\n",
                [&["a\r", "\rb\r"], &["a", "\rb"], &["a", "b"]],
            ),
            ("", [&[], &[], &[]]),
        ];
        for &(bytes, [unix, windows, darwin]) in cases {
            assert_eq!(lines(bytes, Unix), unix, "{bytes:?} on Linux");
            assert_eq!(lines(bytes, Windows), windows, "{bytes:?} on Windows");
            assert_eq!(lines(bytes, Darwin), darwin, "{bytes:?} on macOS");
        }
    }

    // `-` is the caller's standard input (the process's on the command line, a reader of the
    // library API's caller), read in the text mode as a file is.
    #[test]
    fn the_standard_input_is_the_reader_given() {
        let mut stdin: &[u8] = b"a\nb\n";
        let mut cin = Cin::new(&mut stdin);
        let mut bytes = Vec::new();
        open_input_file("in", "-", &mut cin)
            .unwrap()
            .read_to_end(&mut bytes)
            .unwrap();
        assert_eq!(bytes, b"a\nb\n");
        assert_eq!(
            read_input_file("entry_batch", "-", &mut Cin::new(&mut &b"x\n"[..])).unwrap(),
            b"x\n"
        );
    }

    /// A terminal: a first batch, its end (Ctrl-D), then a second batch.
    struct Terminal(Vec<&'static [u8]>);

    impl Read for Terminal {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.0.is_empty() {
                return Ok(0);
            }
            let piece = self.0.remove(0);
            buf[..piece.len()].copy_from_slice(piece);
            Ok(piece.len())
        }
    }

    // `cin` is one stream: once it has ended, a second argument `-` reads nothing (NCBI 2.17.0
    // blastdbcmd `-taxidlist - -ipg_batch -`, review of line RG item K1b, F5).
    #[test]
    fn the_standard_input_reads_nothing_after_its_end() {
        let mut terminal = Terminal(vec![b"1\n", b"", b"2\n"]);
        let mut cin = Cin::new(&mut terminal);
        assert_eq!(read_input_file("taxidlist", "-", &mut cin).unwrap(), b"1\n");
        assert_eq!(read_input_file("ipg_batch", "-", &mut cin).unwrap(), b"");
    }

    /// A writer that takes `room` bytes, then fails.
    struct Room(usize, Vec<u8>);

    impl Write for Room {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            if buf.len() > self.0 {
                return Err(std::io::Error::other("no space left"));
            }
            self.0 -= buf.len();
            self.1.extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // NCBI's `ostream` after a failed write: `badbit`, and later writes, even ones that would
    // fit, do nothing; nothing reports it.
    #[test]
    fn an_ostream_drops_every_write_after_the_first_that_fails() {
        let mut room = Room(4, Vec::new());
        let mut out = Ostream::new(&mut room);
        out.write_all(b"abc").unwrap();
        assert!(!out.bad());
        out.write_all(b"de").unwrap();
        assert!(out.bad());
        out.write_all(b"f").unwrap();
        out.flush().unwrap();
        assert!(out.bad() && !out.closed());
        assert_eq!(room.1, b"abc");
    }

    /// A pipe whose reader has closed.
    struct ClosedPipe;

    impl Write for ClosedPipe {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // A closed pipe ends NCBI's process by SIGPIPE outside Windows: every write after it is an
    // error, which ends the program. On Windows it is a failed write like the others.
    #[test]
    fn a_closed_pipe_ends_the_writes_where_sigpipe_ends_ncbi() {
        let mut pipe = ClosedPipe;
        let mut out = Ostream::new(&mut pipe);
        let first = out.write_all(b"a");
        let later = out.write_all(b"b");
        assert!(out.bad());
        if cfg!(windows) {
            assert!(first.is_ok() && later.is_ok() && !out.closed());
        } else {
            assert!(out.closed() && ends_like_sigpipe(&first.unwrap_err()));
            assert!(ends_like_sigpipe(&later.unwrap_err()));
        }
    }

    #[test]
    fn a_stream_that_does_not_open_is_ncbi_s_exception() {
        let mut empty = std::io::empty();
        let mut cin = Cin::new(&mut empty);
        let error = open_input_file("in", "", &mut cin).err().unwrap();
        assert_eq!(
            error,
            ArgException::new(
                "eNoFile",
                "Argument \"in\". File is not accessible".to_string()
            )
        );
        let missing = "no/such/dir/x.fa";
        assert_eq!(
            open_input_file("in", missing, &mut cin)
                .err()
                .unwrap()
                .message(),
            "Argument \"in\". File is not accessible:  `no/such/dir/x.fa'"
        );
        let mut stdout = Vec::new();
        assert!(matches!(
            open_output_file("out", "-", &mut stdout),
            Ok(OutputFile::Stdout(_))
        ));
        assert_eq!(
            open_output_file("out", "/no/such/dir/x.out", &mut stdout)
                .err()
                .unwrap()
                .message(),
            "Argument \"out\". File is not accessible:  `/no/such/dir/x.out'"
        );
    }
}

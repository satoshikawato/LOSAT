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

/// Where an input stream reads from.
enum Source {
    Stdin(std::io::Stdin),
    File(File),
}

impl Read for Source {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Source::Stdin(stdin) => stdin.read(buf),
            Source::File(file) => file.read(buf),
        }
    }
}

/// The stream of an `eInputFile` argument (`AsInputFile()`, without `fBinary`): the standard
/// input or the file, read in the text mode of LOSAT's platform. It is not buffered (wrap it in a
/// `std::io::BufReader` to read lines); a read error is the operating system's (NCBI's
/// libstdc++ stream throws its own text for it, which the program that reads writes).
pub(crate) struct InputFile {
    source: Source,
    text: TextMode,
    /// Bytes of the text mode not handed out yet (Windows only), from `at`.
    pending: Vec<u8>,
    at: usize,
}

impl Read for InputFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        if self.text.platform != LineEnds::Windows {
            return self.source.read(buf);
        }
        let mut raw = [0u8; 8192];
        while self.at == self.pending.len() {
            self.pending.clear();
            self.at = 0;
            if self.text.ended {
                return Ok(0);
            }
            let read = self.source.read(&mut raw)?;
            if read == 0 {
                self.text.finish(&mut self.pending);
                if self.pending.is_empty() {
                    return Ok(0);
                }
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
/// `AsInputFile()` of the argument `name` whose value is `value`: `-` is the standard input
/// (probes K30-K32 of line BDC), the empty name opens nothing (K39), and another name opens the
/// file; on Linux a directory opens too (its reads fail). `Err` is NCBI's `eNoFile`.
pub(crate) fn open_input_file(name: &str, value: &str) -> Result<InputFile, ArgException> {
    let source = if value == "-" {
        Source::Stdin(std::io::stdin())
    } else if value.is_empty() {
        return Err(not_accessible(name, value));
    } else {
        Source::File(File::open(value).map_err(|_| not_accessible(name, value))?)
    };
    Ok(InputFile {
        source,
        text: TextMode::new(LineEnds::HOST),
        pending: Vec::new(),
        at: 0,
    })
}

/// The bytes of the input file `value` of the argument `name` (`AsInputFile`) as the loops of
/// `NcbiGetlineEOL` over its stream read them (blastdbcmd's batch files, `input_lines`): `Err`
/// with NCBI's message (`File is not accessible`) when it does not open or a read fails. On Linux
/// the stream of a directory opens and reads nothing there (probes K33-K35 of line BDC).
pub(crate) fn read_input_file(name: &str, value: &str) -> Result<Vec<u8>, String> {
    if cfg!(unix) && value != "-" && !value.is_empty() && std::path::Path::new(value).is_dir() {
        return Ok(Vec::new());
    }
    let mut file = open_input_file(name, value).map_err(|error| error.message().to_string())?;
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

    #[test]
    fn a_stream_that_does_not_open_is_ncbi_s_exception() {
        let error = open_input_file("in", "").err().unwrap();
        assert_eq!(
            error,
            ArgException::new(
                "eNoFile",
                "Argument \"in\". File is not accessible".to_string()
            )
        );
        let missing = "no/such/dir/x.fa";
        assert_eq!(
            open_input_file("in", missing).err().unwrap().message(),
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

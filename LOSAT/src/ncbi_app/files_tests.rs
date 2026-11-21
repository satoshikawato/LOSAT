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

/// A source that hands out `piece` bytes at a time.
struct Pieces<'a>(&'a [u8], usize);

impl Read for Pieces<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let count = self.1.min(buf.len()).min(self.0.len());
        buf[..count].copy_from_slice(&self.0[..count]);
        self.0 = &self.0[count..];
        Ok(count)
    }
}

// The stream reads what the text mode makes of the file on each platform, whatever pieces
// the source and the reader take (the read loop of Windows runs on every host here).
#[test]
fn a_text_reader_reads_the_text_mode_of_its_platform() {
    use LineEnds::{Darwin, Unix, Windows};
    let inputs: &[&[u8]] = &[b"a\r\nb\r\n", b"a\r\rb\r", b"a\n\x1ab\n", b"x\r", b""];
    for &bytes in inputs {
        for platform in [Unix, Windows, Darwin] {
            let expected = text_mode(bytes, platform, bytes.len());
            for piece in 1..4 {
                for chunk in [1, 2, 8192] {
                    let mut reader = TextReader::new(Pieces(bytes, piece), platform);
                    let mut read = Vec::new();
                    let mut buf = vec![0u8; chunk];
                    loop {
                        let count = reader.read(&mut buf).unwrap();
                        if count == 0 {
                            break;
                        }
                        read.extend_from_slice(&buf[..count]);
                    }
                    assert_eq!(read, expected, "{bytes:?} {platform:?} {piece} {chunk}");
                    // Ended: nothing more, and an empty buffer does not end it early.
                    assert_eq!(reader.read(&mut buf).unwrap(), 0);
                }
            }
        }
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

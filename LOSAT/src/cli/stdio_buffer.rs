//! The model of NCBI's buffers of a report stream: glibc's buffer of standard output and
//! libstdc++'s `basic_filebuf` of an `ofstream` (`NcbiStdio`), and what each holds of the
//! bytes inserted into it (`NcbiBuffer`). `cli` re-exports them.

// glibc 2.17 to 2.39: libio/filedoalloc.c `_IO_file_doallocate`
/// The size of glibc's buffer of standard output (`cout` writes into it, synchronized with C's
/// `stdout`): `st_blksize` of the standard output where it is below `BUFSIZ`, else `BUFSIZ` (4096
/// for a pipe and for most files, 1024 for a terminal); 4096 where LOSAT cannot ask.
pub fn stdio_buffer_size() -> usize {
    #[cfg(unix)]
    let block = rustix::fs::fstat(std::io::stdout())
        .ok()
        .and_then(|stat| usize::try_from(stat.st_blksize).ok())
        .filter(|&size| size > 0 && size < NCBI_BUFSIZ)
        .unwrap_or(NCBI_BUFSIZ);
    #[cfg(not(unix))]
    let block = 4096;
    block
}

/// glibc's and libstdc++'s `BUFSIZ`.
pub const NCBI_BUFSIZ: usize = 8192;

/// NCBI's buffer of a report stream, which decides when the bytes inserted into the stream reach
/// the file or device (`NcbiBuffer`). The buffers are those of NCBI's builds of BLAST+ and
/// rpsbproc, which are Linux builds (glibc, GCC 13's libstdc++); LOSAT keeps their counts on every
/// platform.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NcbiStdio {
    /// `cout`, synchronized with C's `stdout` (NCBI does not call `sync_with_stdio(false)`):
    /// each insertion is an `fwrite` into glibc's buffer of `block` bytes (`stdio_buffer_size`),
    /// which writes at each newline too where `lines` is set (a terminal: `_IO_LINE_BUF`).
    Cout { block: usize, lines: bool },
    /// An `ofstream` (`-out`, rpsbproc's `-o`): libstdc++'s `basic_filebuf` with a buffer of
    /// `BUFSIZ` bytes, of which an insertion that does not fit in what is left is written
    /// directly with it.
    File,
}

impl NcbiStdio {
    /// The stream of a report on standard output: glibc's buffer of it, line-buffered on a
    /// terminal (`_IO_file_doallocate`, cited at `stdio_buffer_size`: `local_isatty`, as
    /// `is_terminal`).
    pub fn cout() -> Self {
        use std::io::IsTerminal;
        NcbiStdio::Cout {
            block: stdio_buffer_size(),
            lines: std::io::stdout().is_terminal(),
        }
    }
}

/// What NCBI's buffer of a report stream (`NcbiStdio`) holds of the bytes inserted into it: the
/// last `held` bytes; the others are written to the file or device. One model for the organism
/// report's streams (`NcbiReportBuffer`) and rpsbproc's report (`rpsbproc::Report`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NcbiBuffer {
    stream: NcbiStdio,
    held: usize,
    /// Whether glibc has made `stdout`'s buffer (at the first insertion of a byte).
    allocated: bool,
}

impl NcbiBuffer {
    /// The empty buffer of `stream`.
    pub fn new(stream: NcbiStdio) -> Self {
        NcbiBuffer {
            stream,
            held: 0,
            allocated: false,
        }
    }

    /// The buffer's stream.
    pub fn stream(&self) -> NcbiStdio {
        self.stream
    }

    /// How many of the bytes inserted the buffer holds (the last ones).
    pub fn held(&self) -> usize {
        self.held
    }

    /// An insertion of `n` bytes (`sputn`: a string or a number), the last `n` bytes of `stream`
    /// (the stream's bytes up to the insertion's end; a line-buffered `cout` reads the insertion's
    /// newlines).
    // glibc 2.17 to 2.39: libio/fileops.c `_IO_new_file_xsputn`
    // Before the first insertion `stdout` has no buffer (no space available): `_IO_OVERFLOW`
    // makes it, and the first insertion keeps `n % block` (glibc 2.39, `T/t2/fix2/
    // glibc_first_write.c` on a file of `st_blksize` 4096: an `fwrite` of 4,096 bytes is written at
    // once, 4,095 and then 1 are held, and 10 more write 4,096).
    // libstdc++ (GCC 13): include/c++/13/bits/fstream.tcc:766-778
    // (`basic_filebuf::xsputn`)
    // A line-buffered `stdout` (a terminal, `_IO_LINE_BUF` from `_IO_file_doallocate`): the same
    // function's branch for it copies, of an insertion that fits, up to its last newline and then
    // writes the buffer (`must_flush`); what follows is written with it where the buffer is under
    // 128 bytes (`do_write` is all of it), else put byte by byte (`_IO_default_xsputn`), where
    // `_IO_new_file_overflow` writes at each newline, as the rest of an insertion that does not
    // fit; so the buffer keeps what follows the last newline of what it keeps (glibc 2.39,
    // `T/t2/fix3/glibc_tty_lines.c` on a pty of `st_blksize` 1024: 100 then 50 ending in a newline
    // write 150; 1000 then 30 with a newline at 27 write 1,028; 500 then 2,000 with newlines at 100
    // and 1,900 write 2,401; line RG T2 review 5's `T/t2/review5/glibc/line_small.c`, three
    // `fwrite`s of `xxxxx\nyyy` with a line buffer of 64 bytes write 9, 18 and 27 bytes, of 1,024
    // bytes 6, 15 and 24).
    // glibc 2.17 to 2.39: libio/fileops.c `_IO_new_file_xsputn`
    // glibc 2.17 to 2.39: libio/fileops.c `_IO_new_file_overflow`
    pub fn sputn(&mut self, n: usize, stream: &[u8]) {
        if n == 0 {
            return;
        }
        self.held = match self.stream {
            NcbiStdio::Cout { block, lines } => {
                let insertion = &stream[stream.len() - n..];
                // `count`: the bytes of the insertion copied into the buffer, up to its last
                // newline where a line-buffered stream's insertion fits (`must_flush`).
                let room = if self.allocated { block - self.held } else { 0 };
                // glibc looks for the newline only in its line-buffered branch, when the
                // insertion fits; a fully buffered stream never scans the bytes.
                let last_newline = || insertion.iter().rposition(|&byte| byte == b'\n');
                let (count, must_flush) = match (lines && n <= room).then(last_newline).flatten() {
                    Some(at) => (at + 1, true),
                    None => (n.min(room), false),
                };
                let to_do = n - count;
                if to_do == 0 && !must_flush {
                    self.held + n
                } else {
                    // The buffer is written (`_IO_OVERFLOW`, which makes it before the first
                    // insertion), then `do_write` bytes of the rest; the buffer keeps the others,
                    // a line-buffered one what follows their last newline.
                    let do_write = if block >= 128 {
                        to_do - to_do % block
                    } else {
                        to_do
                    };
                    let kept = &insertion[n - (to_do - do_write)..];
                    match lines
                        .then(|| kept.iter().rposition(|&byte| byte == b'\n'))
                        .flatten()
                    {
                        Some(at) => kept.len() - at - 1,
                        None => kept.len(),
                    }
                }
            }
            NcbiStdio::File => {
                if n >= NCBI_BUFSIZ - 1 - self.held {
                    0
                } else {
                    self.held + n
                }
            }
        };
        self.allocated = true;
    }

    /// An insertion of one `char` (`put`, `sputc`), the last byte of `stream` (as `sputn`'s): the
    /// buffer takes it while it has room; a full buffer is written with it
    /// (`basic_filebuf::overflow`). NCBI rpsbproc 0.5.1, both builds, with 8,190 bytes held
    /// before a tab of the `SESSION` line (line RG E1c review 2, `s_same_mat8165`,
    /// `s_same_ev8156`). On `cout` the `char` is glibc's `putc` (`stdio_sync_filebuf::overflow`),
    /// which a full or unallocated buffer passes to `__overflow`: the buffer is made or written,
    /// then it keeps the byte (an `fwrite` of one byte keeps none where the buffer is under 128
    /// bytes, `sputn`), and a line-buffered `stdout` writes it with what it holds at a newline.
    // libstdc++ (GCC 13): include/c++/13/ostream:570-576 (a `char` is put, not
    // written with `sputn`)
    // libstdc++ (GCC 13): include/ext/stdio_sync_filebuf.h `overflow`
    // glibc 2.17 to 2.39: libio/fileops.c `_IO_new_file_overflow` (`putc` calls it
    // when `_IO_write_ptr >= _IO_write_end`: a full or unallocated buffer, and every `char` of a
    // line-buffered stream, whose `_IO_write_end` stays at `_IO_write_base`)
    pub fn sputc(&mut self, stream: &[u8]) {
        match self.stream {
            NcbiStdio::Cout { block, lines } => {
                self.held = if lines && stream.last() == Some(&b'\n') {
                    0
                } else if !self.allocated || self.held == block {
                    1
                } else {
                    self.held + 1
                };
                self.allocated = true;
            }
            NcbiStdio::File => {
                self.held = if self.held < NCBI_BUFSIZ - 1 {
                    self.held + 1
                } else {
                    0
                };
            }
        }
    }

    /// The stream flushed (`flush`, `endl`, a write to `cerr`, which is tied to `cout`): the
    /// buffer writes what it holds.
    pub fn flush(&mut self) {
        self.held = 0;
    }
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;

    #[test]
    fn a_char_on_standard_output_is_kept_as_glibcs_putc_keeps_it() {
        // glibc 2.39 (`T/t2/fix4/glibc_putc.log`, and line RG T2 review 5's
        // `T/t2/review5/glibc/line_small.log`): `stdout` with a buffer of `block` bytes (`setvbuf`,
        // full or line buffering) on a file; `cN` is N `putc` of 'x', `n` a `putc` of a newline,
        // `wN` an `fwrite` of N bytes (`wN:K`: with a newline at K), `W` an `fwrite` of a newline;
        // the file's size after each step. `putc` keeps its byte after a full or unallocated
        // buffer, where an `fwrite` of one byte into a buffer under 128 bytes keeps none; a line
        // buffer under 128 bytes keeps none of an `fwrite` that fits with a newline either.
        let run = |lines: bool, block: usize, steps: &str| -> Vec<usize> {
            let mut buffer = NcbiBuffer::new(NcbiStdio::Cout { block, lines });
            let mut stream = Vec::new();
            steps
                .split(' ')
                .map(|step| {
                    let mut numbers = step[1..].split(':').map(|n| n.parse::<usize>());
                    let count = numbers.next().and_then(Result::ok).unwrap_or(1);
                    match step.as_bytes()[0] {
                        b'c' => (0..count).for_each(|_| {
                            stream.push(b'x');
                            buffer.sputc(&stream);
                        }),
                        b'n' => {
                            stream.push(b'\n');
                            buffer.sputc(&stream);
                        }
                        b'w' => {
                            let start = stream.len();
                            stream.resize(start + count, b'x');
                            numbers.for_each(|at| stream[start + at.unwrap()] = b'\n');
                            buffer.sputn(count, &stream);
                        }
                        _ => {
                            stream.push(b'\n');
                            buffer.sputn(1, &stream);
                        }
                    }
                    stream.len() - buffer.held()
                })
                .collect()
        };
        assert_eq!(run(false, 64, "c1"), [0]);
        assert_eq!(run(false, 64, "w1"), [1]);
        assert_eq!(run(false, 64, "w63 c64 c1 w1"), [63, 63, 127, 127]);
        assert_eq!(run(false, 64, "w63 c64 w1"), [63, 63, 128]);
        assert_eq!(run(false, 64, "w64 c1"), [64, 64]);
        assert_eq!(run(false, 64, "w64 w1"), [64, 64]);
        assert_eq!(run(false, 4096, "w4096 c1"), [4096, 4096]);
        assert_eq!(run(false, 4096, "w4096 w1"), [4096, 4096]);
        assert_eq!(run(false, 4096, "w4095 c1 c1"), [0, 0, 4096]);
        assert_eq!(run(false, 4096, "w4095 w1 w1"), [0, 0, 4096]);
        assert_eq!(run(true, 64, "c1"), [0]);
        assert_eq!(run(true, 64, "n"), [1]);
        assert_eq!(run(true, 64, "W"), [1]);
        assert_eq!(run(true, 64, "c10 n"), [0, 11]);
        assert_eq!(run(true, 64, "w63 c1 c1"), [63, 63, 63]);
        assert_eq!(run(true, 64, "w63 w1 w1"), [63, 63, 63]);
        assert_eq!(run(true, 64, "c64 c1"), [0, 64]);
        assert_eq!(run(true, 64, "c63 n"), [0, 64]);
        assert_eq!(run(true, 64, "c64 n"), [0, 65]);
        assert_eq!(run(true, 1024, "c100 n c5"), [0, 101, 101]);
        assert_eq!(run(true, 64, "w9:5 w9:5 w9:5"), [9, 18, 27]);
        assert_eq!(run(true, 1024, "w9:5 w9:5 w9:5"), [6, 15, 24]);
    }
}

//! NCBI's buffer of a report stream, as its Linux builds make it: where the bytes of an insertion
//! (`operator<<`) reach the file or device, from the insertions alone.
//!
//! Two buffers exist. Standard output (`cout`, the stream that NCBI's formatter and rpsbproc write
//! to) is C's `stdout`: a fully buffered stream of `block` bytes, or a line-buffered one on a
//! terminal ("Stream Buffering" of the C library's manual; the size is the block size of the
//! device, `stat`'s `st_blksize`, capped at `BUFSIZ`). The `-out` file (and rpsbproc's `-o`) is a
//! C++ `std::ofstream`, whose file buffer is of its own (`BUFSIZ` = 8192 bytes of which one is
//! not usable).
//!
//! The rules here are read from the writes that real programs make, not from any source: the
//! tables of `tests/fixtures/stdio_buffer/` (the `strace` of a program that inserts the same
//! operations into `std::cout` and `std::ofstream`, and of NCBI's rpsbproc 0.5.1 and BLAST+
//! 2.17.0) and the tests of this module, which replay every row.

/// The block size below which standard output's buffer writes whole insertions
/// (`tests/fixtures/stdio_buffer/cout_blksize_*.tsv`).
const SMALL_BLOCK: usize = 128;

/// The C library's `BUFSIZ` of NCBI's builds.
pub const NCBI_BUFSIZ: usize = 8192;

/// The buffer size that a standard output whose `st_blksize` is `blksize` gets: the block size
/// below `NCBI_BUFSIZ`, else `NCBI_BUFSIZ` (4096 for a pipe and most files, 1024 for a terminal;
/// observed for 1 to 4,194,304, `tests/fixtures/stdio_buffer/cout_blksize_*.tsv`). A `blksize` of
/// 0 is not observed (the tables start at 1; NCBI's build gets it only from a device that reports
/// no block size): it gets `NCBI_BUFSIZ`, as does a standard output whose `fstat` fails
/// (`stdio_buffer_size`), also not observed.
#[cfg_attr(not(unix), allow(dead_code))]
pub fn block_for(blksize: u64) -> usize {
    if blksize > 0 && blksize < NCBI_BUFSIZ as u64 {
        blksize as usize
    } else {
        NCBI_BUFSIZ
    }
}

/// The size of the buffer that NCBI's BLAST+ gives standard output, for the standard output of
/// this process (`block_for` of its `st_blksize`). The rule is observed on NCBI's Linux builds
/// (glibc 2.35 and 2.39). Other Unix systems (macOS, BSD) use the same rule, unobserved; Windows
/// and WASI (the Web build) have no `st_blksize` here and get 4096, with `NcbiStdio::cout`'s
/// `lines` from `is_terminal`, unobserved.
pub fn stdio_buffer_size() -> usize {
    #[cfg(unix)]
    {
        use std::os::fd::AsFd;
        use std::os::unix::fs::MetadataExt;
        std::io::stdout()
            .as_fd()
            .try_clone_to_owned()
            .and_then(|fd| std::fs::File::from(fd).metadata())
            .map_or(NCBI_BUFSIZ, |metadata| block_for(metadata.blksize()))
    }
    #[cfg(not(unix))]
    {
        4096
    }
}

/// Which buffer a report stream has: standard output with a buffer of `block` bytes, written at
/// each newline too when `lines` (a terminal); or an output file (`-out`, rpsbproc's `-o`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NcbiStdio {
    Cout { block: usize, lines: bool },
    File,
}

impl NcbiStdio {
    /// The stream of this process's standard output.
    pub fn cout() -> NcbiStdio {
        NcbiStdio::Cout {
            block: stdio_buffer_size(),
            lines: stdout_is_terminal(),
        }
    }
}

fn stdout_is_terminal() -> bool {
    use std::io::IsTerminal;
    std::io::stdout().is_terminal()
}

/// How many bytes follow the last newline among the last `scan` bytes of `stream`, if there is one.
fn bytes_after_newline(stream: &[u8], scan: usize) -> Option<usize> {
    let tail = &stream[stream.len().saturating_sub(scan)..];
    tail.iter()
        .rposition(|&byte| byte == b'\n')
        .map(|at| tail.len() - at - 1)
}

/// What one buffer holds of the bytes inserted into it, and whether it exists yet.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NcbiBuffer {
    stream: NcbiStdio,
    held: usize,
    /// Standard output's buffer is made by its first insertion.
    made: bool,
}

impl NcbiBuffer {
    /// An empty buffer that has written nothing.
    pub fn new(stream: NcbiStdio) -> NcbiBuffer {
        NcbiBuffer {
            stream,
            held: 0,
            made: false,
        }
    }

    /// The stream it was made for.
    pub fn stream(&self) -> NcbiStdio {
        self.stream
    }

    /// How many of the last bytes inserted the buffer still holds; all bytes before them have
    /// reached the file or device.
    pub fn held(&self) -> usize {
        self.held
    }

    /// The buffer holding `held` bytes, where its stream keeps the bytes and counts them itself
    /// (`report::ostream`'s report stream, after its writes to the device).
    pub(crate) fn set_held(&mut self, held: usize) {
        debug_assert!(held == 0 || self.made || self.stream == NcbiStdio::File);
        self.held = held;
    }

    /// An insertion (`sputn`) of fewer bytes than this is held whole, with no write: the room
    /// left in a made buffer of a stream that is not a terminal's, plus one (a full block is
    /// written only by the next insertion); 0 where the insertion must be read (standard output
    /// before its first insertion, and a terminal's, which writes at its newlines). The report
    /// stream skips the model for such insertions (`report::ostream`).
    pub fn room(&self) -> usize {
        match self.stream {
            NcbiStdio::File => NCBI_BUFSIZ - 1 - self.held,
            NcbiStdio::Cout {
                block,
                lines: false,
            } if self.made => block - self.held + 1,
            NcbiStdio::Cout { .. } => 0,
        }
    }

    /// An insertion of the last `n` bytes of `stream` (a string or a number; `stream` is
    /// everything inserted so far, or the insertion alone: no byte before the last `n` is read).
    /// Afterwards `held` is what NCBI's buffer holds.
    pub fn sputn(&mut self, n: usize, stream: &[u8]) {
        if n == 0 {
            return;
        }
        match self.stream {
            NcbiStdio::Cout { block, lines } if block < SMALL_BLOCK => {
                // A stream of a few bytes (a block size below 128, which no device reports but
                // which the tables observe): its first insertion is written at once, a later
                // one is held while it fits in the block and else written with all that is
                // held, and on a terminal any insertion with a newline is written whole.
                self.held = if !self.made {
                    0
                } else if self.held + n <= block {
                    self.held + n
                } else {
                    0
                };
                self.made = true;
                if lines && self.held > 0 && bytes_after_newline(stream, n).is_some() {
                    self.held = 0;
                }
            }
            NcbiStdio::Cout { block, lines } => {
                self.held = if self.made {
                    let room = block.saturating_sub(self.held);
                    if n <= room {
                        self.held + n
                    } else {
                        (n - room) % block
                    }
                } else {
                    n % block
                };
                self.made = true;
                // On a terminal what follows the last newline is held (the bytes held before
                // this insertion have none).
                if lines && self.held > 0 {
                    if let Some(after) = bytes_after_newline(stream, n.min(self.held)) {
                        self.held = after;
                    }
                }
            }
            NcbiStdio::File => {
                self.held = if self.held + n < NCBI_BUFSIZ - 1 {
                    self.held + n
                } else {
                    0
                };
            }
        }
    }

    /// An insertion of the one byte that ends `stream` (`put`).
    pub fn sputc(&mut self, stream: &[u8]) {
        match self.stream {
            NcbiStdio::Cout { block, lines } => {
                self.held = if self.made && self.held >= block {
                    1
                } else {
                    self.held + 1
                };
                self.made = true;
                if lines && stream.last() == Some(&b'\n') {
                    self.held = 0;
                }
            }
            NcbiStdio::File => {
                self.held = if self.held >= NCBI_BUFSIZ - 1 {
                    0
                } else {
                    self.held + 1
                };
            }
        }
    }

    /// A flush (`endl`, `flush`, a write to standard error); holds 0.
    pub fn flush(&mut self) {
        self.held = 0;
    }
}

#[cfg(test)]
#[path = "stdio_buffer_tests.rs"]
mod tests;

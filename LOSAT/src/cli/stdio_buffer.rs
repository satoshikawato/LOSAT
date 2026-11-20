//! NCBI's buffer of a report stream, as its Linux builds make it: where the bytes of an insertion
//! (`out << value`) reach the file or device, from the insertions alone.
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

/// The size of the buffer that NCBI's BLAST+ (a Linux build) gives standard output, for the
/// standard output of this process: the block size of what it writes to where that is below
/// `NCBI_BUFSIZ`, else `NCBI_BUFSIZ` (4096 for a pipe and most files, 1024 for a terminal); 4096
/// on a platform that is not Unix.
pub fn stdio_buffer_size() -> usize {
    #[cfg(unix)]
    {
        use std::os::fd::AsFd;
        use std::os::unix::fs::MetadataExt;
        let block = std::io::stdout()
            .as_fd()
            .try_clone_to_owned()
            .and_then(|fd| std::fs::File::from(fd).metadata())
            .map(|metadata| metadata.blksize());
        match block {
            Ok(block) if block > 0 && block < NCBI_BUFSIZ as u64 => block as usize,
            _ => NCBI_BUFSIZ,
        }
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

    /// An insertion of the last `n` bytes of `stream` (a string or a number; `stream` is
    /// everything inserted so far). Afterwards `held` is what NCBI's buffer holds.
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
                if lines && self.held > 0 {
                    let inserted = &stream[stream.len().saturating_sub(n)..];
                    if inserted.contains(&b'\n') {
                        self.held = 0;
                    }
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
                if lines && self.held > 0 {
                    let tail = &stream[stream.len().saturating_sub(self.held)..];
                    if let Some(at) = tail.iter().rposition(|&byte| byte == b'\n') {
                        self.held = tail.len() - at - 1;
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
        let _ = stream;
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

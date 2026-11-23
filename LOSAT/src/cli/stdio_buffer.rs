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

/// The process's standard output as LOSAT's command line writes it (`standard_output`): each
/// write goes to the descriptor (the handle on Windows) at once and nothing of it is held, so the
/// bytes of a write that fails are lost there, as the bytes of NCBI's buffer are (a failed write
/// leaves glibc's buffer of `cout` empty: `crate::writedb::cout`), and no later write or the exit
/// writes them. Rust's `Stdout` keeps up to 1023 bytes of a failed write in its line buffer and
/// writes them when a later write or the exit's flush succeeds (a failure that passes: `EAGAIN` on
/// a non-blocking standard output, a device that gets space back; line RG item KC-11b-2). The
/// writers that hold bytes as NCBI's buffers do (`NcbiBuffer`: `NcbiReportBuffer`, rpsbproc's
/// report, blastdbcmd's `-metadata`) and LOSAT's own `BufWriter`s are over it. The library leaves
/// the process's standard output to `Stdout`, which its caller shares.
pub struct StandardOutput(Device);

/// Where `StandardOutput` writes.
enum Device {
    /// Descriptor 1 (`STD_OUTPUT_HANDLE` on Windows) viewed as a file and never closed.
    Raw(std::mem::ManuallyDrop<std::fs::File>),
    /// A Windows console, through Rust's `Stdout` (which writes UTF-8 to it as UTF-16), flushed
    /// at each write: a console's writes do not fail and pass.
    #[cfg(windows)]
    Console(std::io::Stdout),
    /// No standard output (a Windows process started without the handle): the writes succeed and
    /// go nowhere, as Rust's `Stdout` makes them.
    #[cfg(windows)]
    Missing,
}

/// The process's standard output, written without a buffer (`StandardOutput`). The Rust standard
/// library writes its `Stdout` the same way under its line buffer (descriptor 1 in a
/// `ManuallyDrop`; on Windows the handle, a console as UTF-16), as `blastinput::input_files`
/// reads the standard input.
pub fn standard_output() -> StandardOutput {
    #[cfg(any(unix, target_os = "wasi"))]
    {
        use std::os::fd::{AsRawFd, FromRawFd};
        // SAFETY: descriptor 1 is open for the whole process (on Unix the Rust runtime opens
        // `/dev/null` there before `main` when it was closed; on WASI it is the host's standard
        // output) and LOSAT never closes it. The `File` is held in a `ManuallyDrop` that is never
        // dropped, so it only borrows the descriptor: it is never closed through this `File`.
        let file = unsafe { std::fs::File::from_raw_fd(std::io::stdout().as_raw_fd()) };
        StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)))
    }
    #[cfg(windows)]
    {
        use std::io::IsTerminal;
        use std::os::windows::io::{AsRawHandle, FromRawHandle};
        let stdout = std::io::stdout();
        let handle = stdout.as_raw_handle();
        if handle.is_null() {
            return StandardOutput(Device::Missing);
        }
        if stdout.is_terminal() {
            return StandardOutput(Device::Console(stdout));
        }
        // SAFETY: the standard output handle is open for the whole process and LOSAT never
        // closes it. The `File` is held in a `ManuallyDrop` that is never dropped, so it only
        // borrows the handle: it is never closed through this `File`.
        let file = unsafe { std::fs::File::from_raw_handle(handle) };
        StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)))
    }
}

/// Whether `error` is the one of a descriptor or handle that is not open, whose writes Rust's
/// `Stdout` reports as done (and `StandardOutput` too).
fn is_closed_descriptor(error: &std::io::Error) -> bool {
    #[cfg(unix)]
    let closed = rustix::io::Errno::BADF.raw_os_error();
    // WASI preview 1's `errno::badf`.
    #[cfg(target_os = "wasi")]
    let closed = 8;
    #[cfg(windows)]
    let closed = windows_sys::Win32::Foundation::ERROR_INVALID_HANDLE as i32;
    #[cfg(not(any(unix, windows, target_os = "wasi")))]
    let closed = i32::MIN;
    error.raw_os_error() == Some(closed)
}

impl std::io::Write for StandardOutput {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match &mut self.0 {
            Device::Raw(file) => match std::io::Write::write(&mut **file, buf) {
                Err(error) if is_closed_descriptor(&error) => Ok(buf.len()),
                written => written,
            },
            #[cfg(windows)]
            Device::Console(stdout) => {
                std::io::Write::write_all(stdout, buf)?;
                std::io::Write::flush(stdout)?;
                Ok(buf.len())
            }
            #[cfg(windows)]
            Device::Missing => Ok(buf.len()),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
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

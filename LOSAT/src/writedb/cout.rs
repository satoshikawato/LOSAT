//! NCBI's `cout` as a stream with a state: the stream (`super::ofstream::BasicStream`) over
//! standard output's buffer, for the serializer of blastdbcmd's `-metadata`
//! (`super::metadata`), whose failures NCBI reports, for makeblastdb's log
//! (`super::makeblastdb`), whose failures it does not, and for blast_formatter's `-out -`
//! (`crate::blast_formatter::archive`), whose standard-output write failures abort as NCBI's C++
//! runtime does (line RG-F3 C6; `blast_formatter_standard_output_failure_aborts`).
//!
//! When the bytes of an insertion reach standard output is `crate::cli::NcbiBuffer`'s rule (read
//! from the writes of NCBI's programs, line LIB-O O1b). What a failed write leaves is read from
//! blastdbcmd 2.17.0 (the official Linux x64 release) writing `-metadata` to standard output on
//! a full device and on a file under a size limit (`RLIMIT_FSIZE`, SIGXFSZ ignored: the write
//! that reaches the limit is short, the next fails), 161 runs over JSON documents of 401 to
//! 21,408 bytes and limits at and around every 4096-byte block: the write that fails sets
//! `badbit` on the insertion or the flush that made it, and the buffer holds nothing afterwards
//! (a flush after a failed one writes nothing and succeeds: the run whose document fits and whose
//! last line feed fails exits with status 0 and no message). This model is written from those
//! runs and the C++ standard's text on the stream state, not from the source of a C++ runtime or
//! of a C library.
//!
//! On the command line the output under it is the process's standard output with no buffer
//! of Rust's (`crate::cli::standard_output`). A failure that passes before the exit (`EAGAIN`
//! on a non-blocking standard output) loses the bytes of the failed write there too, as NCBI's
//! buffer does (line RG item KC-11b-2).

use std::io::Write;

use super::ofstream::{BasicStream, StreamBuf};
use crate::cli::{NcbiBuffer, NcbiStdio};
use crate::ncbi_app::files::ends_like_sigpipe;

/// Standard output's buffer over `out`: the bytes it holds and when it writes them
/// (`NcbiBuffer`).
pub(crate) struct CoutBuf<'a> {
    out: &'a mut dyn Write,
    buffer: NcbiBuffer,
    held: Vec<u8>,
    closed_pipe: bool,
}

/// NCBI's `cout` over `out`.
pub(crate) type Cout<'a> = BasicStream<CoutBuf<'a>>;

/// `cout` over `out`, the standard output of this process, whose buffer is `stream`'s
/// (`NcbiStdio::cout()` for the process's own standard output).
pub(crate) fn cout(out: &mut dyn Write, stream: NcbiStdio) -> Cout<'_> {
    BasicStream::with_buf(CoutBuf::new(out, stream))
}

impl<'a> CoutBuf<'a> {
    /// Standard output's buffer over `out`, as `stream`'s (`cout`).
    pub(crate) fn new(out: &'a mut dyn Write, stream: NcbiStdio) -> Self {
        CoutBuf {
            out,
            buffer: NcbiBuffer::new(stream),
            held: Vec::new(),
            closed_pipe: false,
        }
    }
}

impl CoutBuf<'_> {
    /// After `NcbiBuffer` took the bytes `added`: writes what it no longer holds, the bytes held
    /// before first; whether the write succeeded. A failed write leaves the buffer empty.
    fn settle(&mut self, added: &[u8]) -> bool {
        let keep = self.buffer.held();
        let total = self.held.len() + added.len();
        if keep == total {
            self.held.extend_from_slice(added);
            return true;
        }
        let mut all = std::mem::take(&mut self.held);
        all.extend_from_slice(added);
        self.held = all.split_off(total - keep);
        self.write_out(&all)
    }

    /// Writes `bytes` to the output and on to its device (a writer that buffers, as Rust's
    /// standard output does, is flushed). A failure empties the buffer.
    fn write_out(&mut self, bytes: &[u8]) -> bool {
        match self.out.write_all(bytes).and_then(|()| self.out.flush()) {
            Ok(()) => true,
            Err(error) => {
                self.closed_pipe |= ends_like_sigpipe(&error);
                self.held.clear();
                self.buffer.flush();
                false
            }
        }
    }
}

impl StreamBuf for CoutBuf<'_> {
    fn put_bytes(&mut self, data: &[u8]) -> bool {
        self.buffer.sputn(data.len(), data);
        self.settle(data)
    }

    fn put_byte(&mut self, byte: u8) -> bool {
        self.buffer.sputc(&[byte]);
        self.settle(&[byte])
    }

    fn sync(&mut self) -> bool {
        self.buffer.flush();
        let held = std::mem::take(&mut self.held);
        held.is_empty() || self.write_out(&held)
    }

    fn closed_pipe(&self) -> bool {
        self.closed_pipe
    }
}

#[cfg(test)]
#[path = "cout_tests.rs"]
mod tests;

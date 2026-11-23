//! makeblastdb's output file stream: an output-only `std::ofstream` as makeblastdb 2.17.0 (the
//! official Linux x64 release) has it; blast_formatter's `-out` file of the same release is one
//! too (`blast_formatter::archive`), whose formatter makes it throw on `badbit` (`exceptions`).
//!
//! makeblastdb writes every database file (`.phr`, `.psq`, `.pin`, the ISAM files, `.pos`,
//! `.pot`, `.pto`, the alias file, the metadata file ...) through an `ofstream` whose state it
//! rarely looks at. A write that succeeds puts the same bytes in the file however they are
//! buffered; the stream matters where something fails: which bytes reach the disk before a
//! failure, that a stream that failed does nothing more, the state bits a caller reads
//! (`rdstate`, `fail`), and when the failing system call happens (the position of a metadata
//! error message, the exit status of a failing flush).
//!
//! This model is written from three sources only, none of them the source of a C++ runtime:
//! the text of the C++ standard (the stream state bits `goodbit`, `badbit`, `failbit`,
//! `eofbit`; `sentry`; the unformatted and formatted output functions; `flush`, `seekp`,
//! `open`, `close`; the file buffer's `overflow`, `sync`, `close`), POSIX (`write`, `writev`,
//! `lseek`, `close`, `ENOSPC`, `EFBIG`), and observation of what makeblastdb's runtime does,
//! made with a program of our own that runs scripted operations on a `std::ofstream` while
//! a preloaded library logs every `write`, `writev`, `lseek` and `close` of the file and can
//! make `write` fail past a byte limit (`LOSAT/tests/fixtures/makeblastdb_ofstream/`, its
//! `README` names the binaries, hosts and dates). The rules below are the ones that
//! reproduce those tables; the standard leaves a file buffer's size and when it calls `write`
//! open. Where a rule is what the observed C++ runtime (GCC 13 series) does and the standard's
//! text says something else or nothing, the rule says so.
//!
//! The call pattern (`writev` of the held bytes and the new bytes, then `write` for what is
//! left) is the Linux x64 one. On other targets the bytes written are the same, because
//! every short count is continued, but the calls differ (Rust's default vectored write on
//! Windows and WASI takes the first non-empty segment only).
//!
//! # Observed rules
//!
//! The file buffer holds data; its put area is one byte less than `BUFSIZ`, 8191 bytes (`BUFSIZ` is
//! 8192). With `b` bytes held:
//!
//! * `write`/`insert` of `n` bytes: held when `n` is below 8191 less `b` (nothing is written).
//!   Otherwise the held bytes and the `n` bytes go out together in one `writev` (two segments, the
//!   held bytes first, even when none are held); a short count is followed by more calls for the
//!   rest, a `writev` while two segments remain and a `write` once the held bytes are done.
//!   When everything went out the buffer is empty. When the calls fail the held bytes stay
//!   held (the `n` bytes are dropped) and the stream sets `badbit`, except that a call of
//!   `n == 0` never sets it.
//! * `put` of one byte: held when `b < 8191`. Otherwise the byte is added to the held bytes
//!   (making `b + 1`, one past the put area) and all of them go out in one `write`; if that
//!   fails the byte stays held, so the next `put` after `clear` holds `b + 2`, and the stream
//!   sets `badbit` (observed behaviour; the standard does not say how much an `overflow` that
//!   failed leaves in the buffer).
//! * `flush`: the held bytes go out in one `write` (nothing when none are held, and no
//!   `fsync`); on failure they stay held and the stream sets `badbit`.
//! * `seekp` to `pos`: the held bytes go out as in `flush`, then `lseek` to `pos` from the start of
//!   the file (`SEEK_SET`); either failing sets `failbit` (a stream with no file fails).
//! * `close`: the held bytes go out as in `flush`, then `close(2)` is called whatever happened
//!   (the held bytes are ignored from then on); a failure of the write, a failure of
//!   `close(2)` (observed with `EIO`, `ENOSPC` and `EDQUOT`, `ENOSPC` also on a full
//!   device; no retry, the descriptor is gone either way), or a stream with no file, sets
//!   `failbit` (the standard: `close()` sets `failbit` when the file buffer's `close` fails).
//!   Nothing is written for a failed `close(2)`, and later writes see a stream with no file.
//! * The destructor does what `close` does, without a state.
//! * `open`: a stream with a file sets `failbit` and keeps it. Otherwise the file is created
//!   or truncated; success clears the state (`clear()`), failure sets `failbit`. The buffer starts
//!   empty.
//! * The `sentry` of `write`, `put`, `flush`, `endl` and `insert`: a stream whose state is not
//!   `goodbit` does nothing, and if `badbit` is set it also sets `failbit` (`eofbit` alone
//!   does nothing without a second bit). `seekp` acts with `eofbit` set and does nothing with
//!   `failbit` set (both as the standard says); with `badbit` set it also sets `failbit`
//!   (observed behaviour, not in the standard's text).
//! * `endl` is `put('\n')` and then `flush`, each with its own sentry.
//! * `exceptions(mask)` sets the mask and then calls `clear(rdstate())`; `clear` and `setstate`
//!   throw `ios_base::failure` when the new state has a bit of the mask (the standard's text). The
//!   model records the throw ([`Stream::take_thrown`]); an operation that throws ends there, as
//!   each operation here ends with its state change. makeblastdb sets no mask.
//!
//! blast_formatter's standard output (`cout`, synchronized with the C library's `stdout`) is a
//! stream of the same rules whose buffer is the C library's `stdout` ([`Stream::cout`]): what
//! reaches the device at each insertion, `put` and `flush` is `cli::NcbiBuffer`'s, the one model of
//! NCBI's standard output (`NcbiStdio::cout`: its block size, line-buffered on a terminal; line
//! LIB-O O1b). Observed with blast_formatter 2.17.0 on `/dev/full` and on a pipe (block size 4096:
//! a Seq-annot of 2390 bytes fails at its close, one of 4845 bytes in its write; line RG item F3,
//! agent C6, 2026-11-23, `remaining-gaps-20261112/F/f3/tmp/c6/pc`). That stream is never closed.
//!
//! `close(2)` is reported on Unix (rustix's `try_close`); other targets close by dropping the
//! file and never see a close error.

use std::fs::{File, OpenOptions};
use std::io::{self, IoSlice, Seek, SeekFrom, Write};
use std::ops::BitOr;
use std::path::Path;

/// `BUFSIZ` of the Linux x64 build: the file buffer's size. Its put area is one byte shorter.
pub(crate) const BUFSIZ: usize = 8192;

/// The bytes the put area holds: `BUFSIZ` less one.
const PUT_AREA: usize = BUFSIZ - 1;

/// `what()` of the `ios_base::failure` that `clear` throws, as the release's runtime words it
/// (observed: NCBI BLAST+ 2.17.0 blast_formatter, the official Linux x64 release, with `-out
/// /dev/full` and with its standard output on `/dev/full`; line RG item F3, agent C6, 2026-11-23,
/// `remaining-gaps-20261112/F/f3/tmp/c6/pb`).
pub(crate) const IOS_FAILURE_WHAT: &str = "basic_ios::clear: iostream error";

/// The stream state (`std::ios_base::iostate`): `goodbit`, `badbit` and `failbit`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct IoState(u8);

impl IoState {
    pub(crate) const GOOD: IoState = IoState(0);
    pub(crate) const BAD: IoState = IoState(1);
    pub(crate) const FAIL: IoState = IoState(4);

    fn contains(self, bit: IoState) -> bool {
        self.0 & bit.0 != 0
    }
}

impl BitOr for IoState {
    type Output = IoState;
    fn bitor(self, other: IoState) -> IoState {
        IoState(self.0 | other.0)
    }
}

/// A file as the stream uses it: `write` and `writev` (one system call each), `lseek`, `close`.
pub(crate) trait RawFile: Sized {
    /// Creates or truncates the file for writing.
    fn open(path: &Path) -> io::Result<Self>;
    /// One `write` (one segment) or `writev` (several): the number of bytes the call took.
    fn write(&mut self, segments: &[&[u8]]) -> io::Result<usize>;
    /// `lseek` to `pos` from the start of the file (`SEEK_SET`).
    fn seek(&mut self, pos: u64) -> io::Result<()>;
    /// `close(2)`: the descriptor is released whether or not it reports an error.
    fn close(self) -> io::Result<()>;
}

/// A file on disk.
pub(crate) struct DiskFile(File);

impl From<File> for DiskFile {
    /// A file already created or truncated (`ncbi_app::files::open_output_file`).
    fn from(file: File) -> Self {
        DiskFile(file)
    }
}

impl RawFile for DiskFile {
    fn open(path: &Path) -> io::Result<Self> {
        OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(path)
            .map(DiskFile)
    }

    fn write(&mut self, segments: &[&[u8]]) -> io::Result<usize> {
        match segments {
            [one] => self.0.write(one),
            [a, b] => self.0.write_vectored(&[IoSlice::new(a), IoSlice::new(b)]),
            _ => unreachable!("the stream writes one or two segments"),
        }
    }

    fn seek(&mut self, pos: u64) -> io::Result<()> {
        self.0.seek(SeekFrom::Start(pos)).map(|_| ())
    }

    #[cfg(unix)]
    fn close(self) -> io::Result<()> {
        use std::os::fd::IntoRawFd;
        let fd = self.0.into_raw_fd();
        // SAFETY: `fd` was just taken out of the `File`, which no longer owns it, and nothing
        // else closes it.
        unsafe { rustix::io::try_close(fd) }.map_err(io::Error::from)
    }

    #[cfg(not(unix))]
    fn close(self) -> io::Result<()> {
        drop(self.0);
        Ok(())
    }
}

/// What the buffer holds of the bytes it takes, and when it writes them.
#[derive(Clone, Copy, Debug)]
enum Buffering {
    /// `basic_filebuf`'s put area of `PUT_AREA` bytes: a file.
    File,
    /// The C library's `stdout` under `cout` (`cli::NcbiBuffer`): standard output.
    Cout(crate::cli::NcbiBuffer),
}

/// The file buffer: the held bytes and the file.
struct FileBuf<F: RawFile> {
    file: Option<F>,
    /// The bytes held; `len()` is `b` of the rules.
    held: Vec<u8>,
    buffering: Buffering,
}

impl<F: RawFile> FileBuf<F> {
    /// Writes all of `segments` (the first call as given, the calls after a short count with
    /// the segments the count finished dropped). `false` when a call fails; the file keeps the
    /// bytes the earlier calls took.
    fn write_all(file: &mut F, segments: &mut [&[u8]]) -> bool {
        let mut segments = segments;
        loop {
            if segments.iter().all(|s| s.is_empty()) {
                return true;
            }
            let mut taken = match file.write(segments) {
                // A call that takes nothing of bytes still to go would repeat for ever.
                Ok(0) => return false,
                Ok(n) => n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(_) => return false,
            };
            while let Some(first) = segments.first() {
                if taken < first.len() || segments.len() == 1 {
                    break;
                }
                taken -= first.len();
                segments = &mut std::mem::take(&mut segments)[1..];
            }
            if let Some(first) = segments.first_mut() {
                *first = &first[taken..];
            }
        }
    }

    /// `sputn`: whether all `data` was taken (held or written).
    fn put_bytes(&mut self, data: &[u8]) -> bool {
        let FileBuf {
            file,
            held,
            buffering,
        } = self;
        let Some(file) = file.as_mut() else {
            return data.is_empty();
        };
        match buffering {
            Buffering::File => {
                if data.len() < PUT_AREA - held.len().min(PUT_AREA) {
                    held.extend_from_slice(data);
                    return true;
                }
                if FileBuf::write_all(file, &mut [&held[..], data]) {
                    held.clear();
                    true
                } else {
                    // The held bytes stay; a call of no bytes takes all of them.
                    data.is_empty()
                }
            }
            Buffering::Cout(buffer) => {
                held.extend_from_slice(data);
                let mut next = *buffer;
                next.sputn(data.len(), held);
                if FileBuf::write_out(file, held, next.held()) {
                    *buffer = next;
                    true
                } else {
                    // The held bytes stay, as in a file's buffer.
                    held.truncate(held.len() - data.len());
                    data.is_empty()
                }
            }
        }
    }

    /// `sputc`.
    fn put_byte(&mut self, byte: u8) -> bool {
        let FileBuf {
            file,
            held,
            buffering,
        } = self;
        let Some(file) = file.as_mut() else {
            return false;
        };
        let overflow = match buffering {
            Buffering::File => held.len() >= PUT_AREA,
            Buffering::Cout(buffer) => {
                held.push(byte);
                let mut next = *buffer;
                next.sputc(held);
                return if FileBuf::write_out(file, held, next.held()) {
                    *buffer = next;
                    true
                } else {
                    held.pop();
                    false
                };
            }
        };
        held.push(byte);
        if !overflow {
            return true;
        }
        if FileBuf::write_all(file, &mut [&held[..]]) {
            held.clear();
            true
        } else {
            false
        }
    }

    /// Standard output's write of what its buffer no longer holds: all of `held` but its last
    /// `keep` bytes, which stay; `false` when the write fails (nothing is dropped then).
    fn write_out(file: &mut F, held: &mut Vec<u8>, keep: usize) -> bool {
        let written = held.len() - keep;
        if written == 0 {
            return true;
        }
        if FileBuf::write_all(file, &mut [&held[..written]]) {
            held.drain(..written);
            true
        } else {
            false
        }
    }

    /// `sync`: the held bytes go out; a stream with no file is in sync.
    fn sync(&mut self) -> bool {
        let Some(file) = self.file.as_mut() else {
            return true;
        };
        if self.held.is_empty() {
            return true;
        }
        if FileBuf::write_all(file, &mut [&self.held]) {
            self.held.clear();
            if let Buffering::Cout(buffer) = &mut self.buffering {
                buffer.flush();
            }
            true
        } else {
            false
        }
    }

    /// `seekpos`: `sync`, then `lseek`.
    fn seek(&mut self, pos: u64) -> bool {
        if self.file.is_none() || !self.sync() {
            return false;
        }
        self.file
            .as_mut()
            .is_some_and(|file| file.seek(pos).is_ok())
    }

    /// `close`: `sync`, then the file is closed whatever `sync` did. `false` for a stream with
    /// no file, a failed `sync` or a failed `close(2)`.
    fn close(&mut self) -> bool {
        if self.file.is_none() {
            return false;
        }
        let synced = self.sync();
        let closed = self.file.take().is_some_and(|file| file.close().is_ok());
        synced && closed
    }
}

/// An output-only `std::ofstream` of makeblastdb 2.17.0 (Linux x64).
pub(crate) struct Stream<F: RawFile> {
    state: IoState,
    /// `exceptions()`.
    exceptions: IoState,
    /// Whether a state change met `exceptions` since [`Stream::take_thrown`] was last asked.
    thrown: bool,
    buf: FileBuf<F>,
}

/// The stream on a file on disk, as the writer uses it.
pub(crate) type OfStream = Stream<DiskFile>;

impl<F: RawFile> Stream<F> {
    /// A stream with no file.
    pub(crate) fn new() -> Self {
        Self::with_buffering(Buffering::File)
    }

    /// A stream with no file whose buffer is the C library's `stdout` of this process
    /// (`cli::NcbiStdio::cout`): blast_formatter's standard output, whose file is the process's
    /// standard output (the library API's writer in its place).
    pub(crate) fn cout() -> Self {
        Self::cout_of(crate::cli::NcbiStdio::cout())
    }

    /// A stream with no file whose buffer is a `stdout` of `stdio` (`cli::NcbiStdio::Cout`).
    pub(crate) fn cout_of(stdio: crate::cli::NcbiStdio) -> Self {
        Self::with_buffering(Buffering::Cout(crate::cli::NcbiBuffer::new(stdio)))
    }

    fn with_buffering(buffering: Buffering) -> Self {
        Stream {
            state: IoState::GOOD,
            exceptions: IoState::GOOD,
            thrown: false,
            buf: FileBuf {
                file: None,
                held: Vec::with_capacity(BUFSIZ),
                buffering,
            },
        }
    }

    /// `ofstream(path)`: the file created or truncated; a failure leaves the stream failed.
    pub(crate) fn create(path: &Path) -> Self {
        let mut stream = Self::new();
        stream.open(path);
        stream
    }

    /// `open` of `path` for binary output.
    pub(crate) fn open(&mut self, path: &Path) {
        self.open_with(|| F::open(path));
    }

    /// `open` where `opener` makes the file; it is not called for a stream that has a file.
    pub(crate) fn open_with(&mut self, opener: impl FnOnce() -> io::Result<F>) {
        if self.buf.file.is_some() {
            self.setstate(IoState::FAIL);
            return;
        }
        match opener() {
            Ok(file) => {
                self.buf.file = Some(file);
                self.buf.held.clear();
                self.clear(IoState::GOOD);
            }
            Err(_) => self.setstate(IoState::FAIL),
        }
    }

    pub(crate) fn rdstate(&self) -> IoState {
        self.state
    }

    /// `clear(state)`.
    pub(crate) fn clear(&mut self, state: IoState) {
        self.state = state;
        self.check_exceptions();
    }

    /// `setstate(state)`.
    pub(crate) fn setstate(&mut self, state: IoState) {
        self.state = self.state | state;
        self.check_exceptions();
    }

    /// `exceptions(mask)`: the mask set, then `clear(rdstate())`, which throws when the state
    /// already has a bit of `mask`.
    pub(crate) fn set_exceptions(&mut self, mask: IoState) {
        self.exceptions = mask;
        self.clear(self.state);
    }

    /// Whether the stream threw `ios_base::failure` since the last call (a state change that met
    /// `exceptions`), which is then forgotten.
    pub(crate) fn take_thrown(&mut self) -> bool {
        std::mem::take(&mut self.thrown)
    }

    /// The throw of `clear`: the state has a bit of `exceptions`.
    fn check_exceptions(&mut self) {
        if self.state.0 & self.exceptions.0 != 0 {
            self.thrown = true;
        }
    }

    /// `fail()`: `failbit` or `badbit`.
    pub(crate) fn fail(&self) -> bool {
        self.state.contains(IoState::FAIL | IoState::BAD)
    }

    pub(crate) fn good(&self) -> bool {
        self.state == IoState::GOOD
    }

    pub(crate) fn bad(&self) -> bool {
        self.state.contains(IoState::BAD)
    }

    /// The `sentry` of an output function: whether the stream may be written.
    fn sentry(&mut self) -> bool {
        if self.good() {
            return true;
        }
        if self.bad() {
            self.setstate(IoState::FAIL);
        }
        false
    }

    /// `write(s, n)`.
    pub(crate) fn write(&mut self, s: &[u8]) {
        if !self.sentry() {
            return;
        }
        if !self.buf.put_bytes(s) {
            self.setstate(IoState::BAD);
        }
    }

    /// `<<` of a string.
    pub(crate) fn insert(&mut self, s: &[u8]) {
        self.write(s);
    }

    /// `put(c)`.
    pub(crate) fn put(&mut self, c: u8) {
        if !self.sentry() {
            return;
        }
        if !self.buf.put_byte(c) {
            self.setstate(IoState::BAD);
        }
    }

    /// `flush()`.
    pub(crate) fn flush(&mut self) {
        if !self.sentry() {
            return;
        }
        if !self.buf.sync() {
            self.setstate(IoState::BAD);
        }
    }

    /// The stream insertion of `endl`: a newline, then `flush()`.
    pub(crate) fn endl(&mut self) {
        self.put(b'\n');
        self.flush();
    }

    /// `seekp` to `pos`.
    pub(crate) fn seekp(&mut self, pos: u64) {
        if self.bad() {
            self.setstate(IoState::FAIL);
        }
        if self.fail() {
            return;
        }
        if !self.buf.seek(pos) {
            self.setstate(IoState::FAIL);
        }
    }

    /// `close()`.
    pub(crate) fn close(&mut self) {
        if !self.buf.close() {
            self.setstate(IoState::FAIL);
        }
    }
}

impl<F: RawFile> std::fmt::Debug for Stream<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OfStream")
            .field("state", &self.state)
            .field("open", &self.buf.file.is_some())
            .field("held", &self.buf.held.len())
            .finish()
    }
}

impl<F: RawFile> Drop for Stream<F> {
    fn drop(&mut self) {
        self.buf.close();
    }
}

#[cfg(test)]
#[path = "ofstream_fixture_tests.rs"]
mod fixture_tests;

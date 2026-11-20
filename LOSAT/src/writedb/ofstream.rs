//! An output file as makeblastdb writes it: `std::ofstream` of the official makeblastdb 2.17.0
//! (Linux x64), which links GCC 13.2.0's libstdc++ statically, so its streams behave as
//! libstdc++ 13.2's `basic_ofstream` over `basic_filebuf` wherever the program runs.
//!
//! Every output file of the database writer is one: the files of a volume (`CWriteDB_File`:
//! the index, header and sequence files, the ISAM files, the OID to GI file), the lookup files of
//! version 5 (`.pos`, `.pot`, `.pto`), the alias file and the metadata file. makeblastdb never
//! tests these streams. A file that cannot be opened leaves the stream failed and every later
//! output does nothing; a write that fails sets `badbit` and every later output does nothing; the
//! build goes on. The metadata file is the one stream whose state NCBI's code reads (its
//! `COStreamBuffer`, `super::metadata`), which makes the points where a write fails part of
//! makeblastdb's output: they follow the put area of 8191 bytes and the rule of `xsputn` below.
//!
//! The model (source: libstdc++-v3 of the gcc-13.2.0 release, whose line numbers are cited; the
//! cited functions are the same in 13.3.0):
//! - the open is `fopen(name, "wb")`: the file is created or truncated; a failure sets `failbit`;
//! - the buffer has `BUFSIZ` = 8192 bytes, of which the put area holds `BUFSIZ - 1`; it starts
//!   "uncommitted" (no put area) until the first output;
//! - `write` of `n` bytes copies them into the put area when `n` is less than the room left in
//!   it (all of it while uncommitted), and otherwise writes the buffered bytes and the `n` bytes
//!   with one `writev`, retried for a partial write, and empties the put area only when all were
//!   written;
//! - `flush`, `seekp` and the close write what the put area holds; a failure keeps it there;
//! - the sentry of every output: a stream that is not good does nothing (one that is bad gets
//!   `failbit` as well).
//!
//! The errors of `close(2)` (`fclose`) are not seen: Rust's `File` drops them. makeblastdb does
//! not test `close()` either, and no output of NCBI's code depends on it.

use std::fs::File;
use std::io::{IoSlice, Seek, SeekFrom, Write};
use std::path::Path;

// libstdc++ 13.2.0 reference: libstdc++-v3/include/std/fstream:49-52, /usr/include/stdio.h:100
/// `_M_buf_size`: glibc's `BUFSIZ`.
pub(crate) const BUFSIZ: usize = 8192;

/// `ios_base::iostate`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct IoState(u8);

impl IoState {
    /// `goodbit`.
    pub(crate) const GOOD: IoState = IoState(0);
    /// `badbit`.
    pub(crate) const BAD: IoState = IoState(1);
    /// `failbit`.
    pub(crate) const FAIL: IoState = IoState(4);

    fn has(self, bits: IoState) -> bool {
        self.0 & bits.0 != 0
    }
}

impl std::ops::BitOr for IoState {
    type Output = IoState;
    fn bitor(self, other: IoState) -> IoState {
        IoState(self.0 | other.0)
    }
}

/// `std::ofstream` (`basic_ofstream<char>` with its `basic_filebuf<char>`), output only.
#[derive(Debug)]
pub(crate) struct OfStream {
    /// `_M_file`: open while `Some`. `_M_mode` is `ios_base::out` exactly while the file is open
    /// (`open` sets both, `close` clears both), so `is_open()` stands for its `__testout` too.
    file: Option<File>,
    /// `_M_buf`: `BUFSIZ` bytes from the open to the close.
    buf: Box<[u8]>,
    /// `pptr() - pbase()`.
    pptr: usize,
    /// `epptr() - pbase()`: 0 without a put area (`setp(0, 0)`), `BUFSIZ - 1` with one.
    epptr: usize,
    /// `_M_writing`.
    writing: bool,
    /// `rdstate()`.
    state: IoState,
}

// libstdc++ 13.2.0 reference: libstdc++-v3/config/io/basic_file_stdio.cc:118-150
/// `xwrite`: `write(2)` until all of `s` is written or one fails; the bytes written.
fn xwrite(file: &mut File, s: &[u8]) -> usize {
    let mut done = 0;
    while done < s.len() {
        match file.write(&s[done..]) {
            Ok(n) => done += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(_) => break,
        }
    }
    done
}

// libstdc++ 13.2.0 reference: libstdc++-v3/config/io/basic_file_stdio.cc:152-197
// `xsputn_2` (basic_file_stdio.cc:356-378) is this function where `writev` exists. Rust's
// `write_vectored` is `writev` on Unix; elsewhere it writes the first slice that is not empty,
// a partial write that the loop continues as it continues one of `writev`.
/// `xsputn_2`: `s1` then `s2` with `writev(2)`; the bytes written.
fn xwritev(file: &mut File, s1: &[u8], s2: &[u8]) -> usize {
    let total = s1.len() + s2.len();
    let mut nleft = total;
    let mut s1 = s1;
    loop {
        let ret = match file.write_vectored(&[IoSlice::new(s1), IoSlice::new(s2)]) {
            Ok(ret) => ret,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(_) => break,
        };
        nleft -= ret;
        if nleft == 0 {
            break;
        }
        if ret >= s1.len() {
            let off = ret - s1.len();
            nleft -= xwrite(file, &s2[off..]);
            break;
        }
        s1 = &s1[ret..];
    }
    total - nleft
}

impl OfStream {
    /// `ofstream()`: a stream with no file (`m_RealFile` of a `CWriteDB_File` made with
    /// `always_create` false, until `Create`).
    pub(crate) fn new() -> Self {
        Self {
            file: None,
            buf: Box::new([]),
            pptr: 0,
            epptr: 0,
            writing: false,
            state: IoState::GOOD,
        }
    }

    /// `ofstream(name)` and `ofstream(name, ios::out | ios::binary)`: the two are one mode on
    /// POSIX systems.
    pub(crate) fn create(path: &Path) -> Self {
        let mut stream = Self::new();
        stream.open(path);
        stream
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/std/fstream:927-936
    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:176-209
    // `_M_file.open` is `fopen64(name, "w")` (`"wb"` with `ios::binary`, the same on POSIX):
    // `O_WRONLY | O_CREAT | O_TRUNC` with mode 0666 (basic_file_stdio.cc:70-116, 246-264), as
    // `File::create`.
    /// `ofstream::open(name)`.
    pub(crate) fn open(&mut self, path: &Path) {
        if self.file.is_none() {
            if let Ok(file) = File::create(path) {
                self.file = Some(file);
                self.buf = vec![0u8; BUFSIZ].into_boxed_slice();
                self.writing = false;
                self.set_buffer(-1);
                self.clear(IoState::GOOD);
                return;
            }
        }
        self.setstate(IoState::FAIL);
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/std/fstream:465-481
    fn set_buffer(&mut self, off: isize) {
        self.pptr = 0;
        self.epptr = if self.file.is_some() && off == 0 {
            BUFSIZ - 1
        } else {
            0
        };
    }

    /// `rdstate()`.
    pub(crate) fn rdstate(&self) -> IoState {
        self.state
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/basic_ios.tcc:39-49
    // An `ofstream` always has its `filebuf` and no exception mask.
    /// `clear(state)`.
    pub(crate) fn clear(&mut self, state: IoState) {
        self.state = state;
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/basic_ios.h:156-158
    /// `setstate(state)`.
    pub(crate) fn setstate(&mut self, state: IoState) {
        self.clear(self.state | state);
    }

    /// `fail()`, which is also `!stream`: `failbit` or `badbit`.
    pub(crate) fn fail(&self) -> bool {
        self.state.has(IoState::FAIL | IoState::BAD)
    }

    /// `good()`.
    pub(crate) fn good(&self) -> bool {
        self.state == IoState::GOOD
    }

    /// `bad()`.
    pub(crate) fn bad(&self) -> bool {
        self.state.has(IoState::BAD)
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/ostream.tcc:45-58
    // An `ofstream` is tied to no stream.
    /// `ostream::sentry`: whether the output may go on.
    fn sentry(&mut self) -> bool {
        if self.good() {
            return true;
        }
        if self.bad() {
            self.setstate(IoState::FAIL);
        }
        false
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/ostream.tcc:186-217
    // The writers call this once for each record, offset or id (millions of times for a large
    // database), so the case that only copies is decided here and inlined: a good stream whose
    // put area is committed (`_M_writing`) and holds more room than `s` needs. There the sentry
    // passes, `xsputn` takes `__streambuf_type::xsputn` (`__n < __bufavail`), which copies all
    // of `s` without an `overflow`, and nothing else changes; every other case takes the whole
    // path (`write_slow`).
    /// `ostream::write(s, n)`.
    #[inline]
    pub(crate) fn write(&mut self, s: &[u8]) {
        if self.state == IoState::GOOD && self.writing && self.pptr < self.epptr {
            let end = self.pptr + s.len();
            if end < self.epptr {
                self.buf[self.pptr..end].copy_from_slice(s);
                self.pptr = end;
                return;
            }
        }
        self.write_slow(s);
    }

    /// `ostream::write(s, n)` in full.
    #[inline(never)]
    fn write_slow(&mut self, s: &[u8]) {
        if self.sentry() && self.xsputn(s) != s.len() {
            self.setstate(IoState::BAD);
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/ostream_insert.h:45-56,76-106
    // The width of these streams stays 0: `<<` of a string is `write` of its bytes.
    /// `ostream << string` (`__ostream_insert`).
    pub(crate) fn insert(&mut self, s: &[u8]) {
        self.write(s);
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/ostream.tcc:152-183
    // libstdc++ 13.2.0 reference: libstdc++-v3/include/std/streambuf:430-443
    /// `ostream::put(c)`.
    pub(crate) fn put(&mut self, c: u8) {
        if self.sentry() {
            let put = if self.pptr < self.epptr {
                self.buf[self.pptr] = c;
                self.pptr += 1;
                true
            } else {
                self.overflow(Some(c))
            };
            if !put {
                self.setstate(IoState::BAD);
            }
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/ostream.tcc:219-252
    /// `ostream::flush()`.
    pub(crate) fn flush(&mut self) {
        if self.sentry() && !self.sync() {
            self.setstate(IoState::BAD);
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/std/ostream:733-736
    /// `stream << endl` (`NcbiEndl`).
    pub(crate) fn endl(&mut self) {
        self.put(b'\n');
        self.flush();
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/ostream.tcc:266-283
    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:885-921
    /// `ostream::seekp(pos)`.
    pub(crate) fn seekp(&mut self, pos: u64) {
        let _ = self.sentry();
        if !self.fail() {
            let mut sought = false;
            if self.file.is_some() && self.terminate_output() {
                if let Some(file) = self.file.as_mut() {
                    if file.seek(SeekFrom::Start(pos)).is_ok() {
                        self.writing = false;
                        self.set_buffer(-1);
                        sought = true;
                    }
                }
            }
            if !sought {
                self.setstate(IoState::FAIL);
            }
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/std/fstream:999-1004
    /// `ofstream::close()`.
    pub(crate) fn close(&mut self) {
        if !self.filebuf_close() {
            self.setstate(IoState::FAIL);
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:245-290
    /// `filebuf::close()`: `false` when it was not open or what it held could not be written.
    fn filebuf_close(&mut self) -> bool {
        if self.file.is_none() {
            return false;
        }
        let written = self.terminate_output();
        self.buf = Box::new([]);
        self.writing = false;
        self.pptr = 0;
        self.epptr = 0;
        self.file = None;
        written
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:751-795
    // The `char` facet never converts and an output stream never reads. `pptr()` passes
    // `epptr()` only after an `overflow` of a full put area that failed; NCBI's code puts at
    // most one byte (`NcbiEndl`) into a stream, so its put area is never full. The room is then
    // taken as none, which takes the branch the negative room of libstdc++ takes.
    /// `filebuf::xsputn` (`sputn`): the bytes of `s` taken.
    fn xsputn(&mut self, s: &[u8]) -> usize {
        let Some(file) = self.file.as_mut() else {
            return self.streambuf_xsputn(s);
        };
        let mut bufavail = self.epptr.saturating_sub(self.pptr);
        if !self.writing {
            bufavail = BUFSIZ - 1;
        }
        if s.len() >= bufavail {
            let buffill = self.pptr;
            let ret = xwritev(file, &self.buf[..buffill], s);
            if ret == buffill + s.len() {
                self.set_buffer(0);
                self.writing = true;
            }
            ret.saturating_sub(buffill)
        } else {
            self.streambuf_xsputn(s)
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/streambuf.tcc:77-109
    /// `streambuf::xsputn`.
    fn streambuf_xsputn(&mut self, s: &[u8]) -> usize {
        let mut ret = 0;
        while ret < s.len() {
            let buf_len = self.epptr - self.pptr;
            if buf_len > 0 {
                let len = buf_len.min(s.len() - ret);
                self.buf[self.pptr..self.pptr + len].copy_from_slice(&s[ret..ret + len]);
                ret += len;
                self.pptr += len;
            }
            if ret < s.len() {
                if self.overflow(Some(s[ret])) {
                    ret += 1;
                } else {
                    break;
                }
            }
        }
        ret
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:537-597
    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:602-614,663
    // The put area ends one byte before the buffer, so the appended byte always has room.
    // `_M_file.xsputn` is `xwrite` (basic_file_stdio.cc:346-353).
    /// `filebuf::overflow(c)` (`None` is `eof()`): `false` for `eof()`, a failure.
    fn overflow(&mut self, c: Option<u8>) -> bool {
        let Some(file) = self.file.as_mut() else {
            return false;
        };
        if self.pptr > 0 {
            if let Some(c) = c {
                self.buf[self.pptr] = c;
                self.pptr += 1;
            }
            let len = self.pptr;
            if xwrite(file, &self.buf[..len]) == len {
                self.set_buffer(0);
                return true;
            }
            false
        } else {
            self.set_buffer(0);
            self.writing = true;
            if let Some(c) = c {
                self.buf[0] = c;
                self.pptr = 1;
            }
            true
        }
    }

    // libstdc++ 13.2.0 reference: libstdc++-v3/include/bits/fstream.tcc:1005-1022
    // `_M_terminate_output` (fstream.tcc:944-1003) is the same test; its second part writes
    // the unshift sequence of a converting facet, which `char` does not have.
    /// `filebuf::sync()` (`pubsync`) and `_M_terminate_output()`: `false` for `-1`.
    fn sync(&mut self) -> bool {
        self.pptr == 0 || self.overflow(None)
    }

    fn terminate_output(&mut self) -> bool {
        self.sync()
    }

    /// The bytes the put area holds (for the tests).
    #[cfg(test)]
    fn pending(&self) -> usize {
        self.pptr
    }
}

// libstdc++ 13.2.0 reference: libstdc++-v3/include/std/fstream:249-259
impl Drop for OfStream {
    fn drop(&mut self) {
        let _ = self.filebuf_close();
    }
}

#[cfg(test)]
#[path = "ofstream_tests.rs"]
mod tests;

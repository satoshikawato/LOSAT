//! The report's output stream as NCBI's writes it: the buffers of NCBI's streams (`NcbiStdio`,
//! `NcbiBuffer`), the command line's report stream (`ReportStream`), the stream that NCBI's
//! formatter writes to (`NcbiOstream`) and the tie of standard error to standard output
//! (`tie_standard_output`).

// NCBI reference (glibc 2.17 to 2.39): libio/filedoalloc.c `_IO_file_doallocate`
// ```c
//   size = BUFSIZ;
//   if (fp->_fileno >= 0 && __builtin_expect (_IO_SYSSTAT (fp, &st), 0) >= 0)
//     {
//       if (S_ISCHR (st.st_mode))
// 	{
// 	  /* Possibly a tty.  */
// 	  if (
// #ifdef DEV_TTY_P
// 	      DEV_TTY_P (&st) ||
// #endif
// 	      local_isatty (fp->_fileno))
// 	    fp->_flags |= _IO_LINE_BUF;
// 	}
// #if defined _STATBUF_ST_BLKSIZE
//       if (st.st_blksize > 0 && st.st_blksize < BUFSIZ)
// 	size = st.st_blksize;
// #endif
//     }
// ```
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
    /// terminal (`_IO_file_doallocate`, quoted at `stdio_buffer_size`: `local_isatty`, as
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
/// report's streams (`ReportStream`) and rpsbproc's report (`rpsbproc::Report`).
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
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_xsputn`
    // ```c
    //   /* First figure out how much space is available in the buffer. */
    //   ...
    //   else if (f->_IO_write_end > f->_IO_write_ptr)
    //     count = f->_IO_write_end - f->_IO_write_ptr; /* Space available. */
    //   ...
    //   if (to_do + must_flush > 0)
    //     {
    //       size_t block_size, do_write;
    //       /* Next flush the (full) buffer. */
    //       if (_IO_OVERFLOW (f, EOF) == EOF)
    //     ...
    //       block_size = f->_IO_buf_end - f->_IO_buf_base;
    //       do_write = to_do - (block_size >= 128 ? to_do % block_size : 0);
    // ```
    // Before the first insertion `stdout` has no buffer (no space available): `_IO_OVERFLOW`
    // makes it, and the first insertion keeps `n % block` (glibc 2.39, `T/t2/fix2/
    // glibc_first_write.c` on a file of `st_blksize` 4096: an `fwrite` of 4,096 bytes is written at
    // once, 4,095 and then 1 are held, and 10 more write 4,096).
    // NCBI reference (GCC 13 libstdc++): include/c++/13/bits/fstream.tcc:766-778
    // (`basic_filebuf::xsputn`)
    // ```c++
    // 	  streamsize __bufavail = this->epptr() - this->pptr();
    //
    // 	  // Don't mistake 'uncommitted' mode buffered with unbuffered.
    // 	  if (!_M_writing && _M_buf_size > 1)
    // 	    __bufavail = _M_buf_size - 1;
    //
    // 	  if (__n >= __bufavail)
    // 	    {
    // 	      const streamsize __buffill = this->pptr() - this->pbase();
    // 	      const char* __buf = reinterpret_cast<const char*>(this->pbase());
    // 	      __ret = _M_file.xsputn_2(__buf, __buffill,
    // 				       reinterpret_cast<const char*>(__s),
    // 				       __n);
    // ```
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
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_xsputn`
    // ```c
    //   if ((f->_flags & _IO_LINE_BUF) && (f->_flags & _IO_CURRENTLY_PUTTING))
    //     {
    //       count = f->_IO_buf_end - f->_IO_write_ptr;
    //       if (count >= n)
    //         {
    //           const char *p;
    //           for (p = s + n; p > s; )
    //             {
    //               if (*--p == '\n')
    //                 {
    //                   count = p - s + 1;
    //                   must_flush = 1;
    //                   break;
    //                 }
    //             }
    //         }
    //     }
    //   ...
    //   /* Then fill the buffer. */
    //   if (count > 0)
    //     {
    //       if (count > to_do)
    // 	count = to_do;
    //       f->_IO_write_ptr = __mempcpy (f->_IO_write_ptr, s, count);
    //       s += count;
    //       to_do -= count;
    //     }
    //   ...
    //       /* Now write out the remainder.  Normally, this will fit in the
    // 	 buffer, but it's somewhat messier for line-buffered files,
    // 	 so we let _IO_default_xsputn handle the general case. */
    //       if (to_do)
    // 	to_do -= _IO_default_xsputn (f, s+do_write, to_do);
    // ```
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_overflow`
    // ```c
    //   if ((f->_flags & _IO_UNBUFFERED)
    //       || ((f->_flags & _IO_LINE_BUF) && ch == '\n'))
    //     if (_IO_do_write (f, f->_IO_write_base,
    // 		      f->_IO_write_ptr - f->_IO_write_base) == EOF)
    // ```
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
    // NCBI reference (GCC 13 libstdc++): include/c++/13/ostream:570-576 (a `char` is put, not
    // written with `sputn`)
    // ```c++
    //   template<typename _Traits>
    //     inline basic_ostream<char, _Traits>&
    //     operator<<(basic_ostream<char, _Traits>& __out, char __c)
    //     {
    //       if (__out.width() != 0)
    // 	return __ostream_insert(__out, &__c, 1);
    //       __out.put(__c);
    //       return __out;
    //     }
    // ```
    // NCBI reference (GCC 13 libstdc++): include/ext/stdio_sync_filebuf.h `overflow`
    // ```c++
    //       virtual int_type
    //       overflow(int_type __c = traits_type::eof())
    //       {
    // 	int_type __ret;
    // 	if (traits_type::eq_int_type(__c, traits_type::eof()))
    // 	  {
    // 	    if (std::fflush(_M_file))
    // 	      __ret = traits_type::eof();
    // 	    else
    // 	      __ret = traits_type::not_eof(__c);
    // 	  }
    // 	else
    // 	  __ret = std::putc(__c, _M_file);
    // 	return __ret;
    //       }
    // ```
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_overflow` (`putc` calls it
    // when `_IO_write_ptr >= _IO_write_end`: a full or unallocated buffer, and every `char` of a
    // line-buffered stream, whose `_IO_write_end` stays at `_IO_write_base`)
    // ```c
    //   if ((f->_flags & _IO_CURRENTLY_PUTTING) == 0 || f->_IO_write_base == NULL)
    //     {
    //       /* Allocate a buffer if needed. */
    //       if (f->_IO_write_base == NULL)
    // 	{
    // 	  _IO_doallocbuf (f);
    //     ...
    //       if (f->_mode <= 0 && f->_flags & (_IO_LINE_BUF | _IO_UNBUFFERED))
    // 	f->_IO_write_end = f->_IO_write_ptr;
    //     }
    //   ...
    //   if (f->_IO_write_ptr == f->_IO_buf_end ) /* Buffer is really full */
    //     if (_IO_do_flush (f) == EOF)
    //       return EOF;
    //   *f->_IO_write_ptr++ = ch;
    //   if ((f->_flags & _IO_UNBUFFERED)
    //       || ((f->_flags & _IO_LINE_BUF) && ch == '\n'))
    //     if (_IO_do_write (f, f->_IO_write_base,
    // 		      f->_IO_write_ptr - f->_IO_write_base) == EOF)
    //       return EOF;
    // ```
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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:116-119
// ```c++
//     m_OrigExceptionMask = m_Outfile.exceptions();
//     m_Outfile.exceptions(NcbiBadbit);
// ```
/// NCBI's region of the report where an insertion is made: what NCBI's code does where the
/// stream's exception (`exceptions(badbit)`, which `CBlastFormat` sets) is thrown. The formatter
/// declares it (`ReportWrite::site`) before it writes there; a report starts in `Report`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Site {
    /// The exception reaches the application's `CATCH_ALL`: "BLAST failed to write output", exit
    /// 6 (`blast_app_util.hpp:252-255`).
    Report,
    /// The alignments of outfmt 0-4 (`CDisplaySeqalign::DisplaySeqalign`), from the flush of
    /// `"\nAlignments" << endl` to the flush at their end: the vestigial `out2`, a
    /// `CObjectOStream` on the report's stream (`showalign.cpp:1931`), flushes the stream again
    /// in its destructor and posts 802.4 (`objostr.cpp:317-324`), then as `Report`. On standard
    /// output the post's `cerr` flushes the failed `cout` first (the tie, `ncbidiag.cpp:6476-6493`)
    /// and that throws inside `~CNcbiDiag`: `std::terminate`.
    Alignments,
    /// A `CBlastTabularInfo` lives (outfmt 6, 7, 10; `blast_format.cpp:773-782`, `:2234-2235`):
    /// its destructor flushes the failed stream again (`tabular.cpp:160-163`), the exception
    /// leaves an implicitly `noexcept` destructor: `std::terminate`.
    Tabular,
    /// SAM (outfmt 17): `CSAM_Formatter` holds the report and writes it in its destructor
    /// (`sam_formatter.cpp:313-316`), a member of `CBlastFormat` destroyed after the body of
    /// `~CBlastFormat` resets the stream's exceptions (`blast_format.hpp:415`): a failed write is
    /// not reported.
    Silent,
    /// One query's `Iteration` of BLAST XML (outfmt 5), from its first insertion to its last
    /// (`s_BlastXMLAddIteration`, inside the `try` of `BlastXML_FormatReport`,
    /// blastxml_format.cpp:745-764): the exception is caught there, `Failed
    /// s_BlastXMLAddIteration Q(0/1` (sic) is posted, and the run goes on with the next query,
    /// whose first insertion into the bad stream throws again. The error of the failed insertion
    /// is `ReportWriteFailure::caught` (`is_caught_report_write_failure`): the query's report stops
    /// there (no later `Hit` is made). On standard output the post's tie throws
    /// (`std::terminate`).
    XmlIteration,
    /// One object that an object stream of NCBI's writes on the report's stream (`m_Outfile <<
    /// MSerial_AsnText << obj`, `WriteObject`, serialobject.cpp:815-827: the Seq-annot of a query in
    /// outfmt 8, 9 and 12): its buffer (`COStreamBuffer`) inserts the object in pieces
    /// (`FlushBuffer`, at the end `EndOfWrite`), and `~CObjectOStream`'s `Close` flushes the
    /// stream (`ReportWrite::ncbi_flush`). A failed insertion throws `ios_base::failure` through
    /// the object's frames (they catch `CException`s only, objstack.hpp:201-214); `Close`, which
    /// clears the stream's state for its flush (`STemporarilyClearStreamState`), fails again on a
    /// file (`basic_filebuf` keeps the bytes) and posts 802.4, and the exception ends the run at
    /// `CATCH_ALL`. Where only `Close`'s flush fails, its destructor's catch posts 802.4 and the run
    /// goes on with a bad stream, whose next insertion throws. On standard output each post's tie
    /// throws (`std::terminate`). `whole`: LOSAT inserts the object in one insertion, and a failure
    /// of its device write is decided by NCBI's pieces (`ReportWrite::resolve_object`, which the
    /// writer calls after an `is_undecided_report_write_failure`); else it inserts NCBI's pieces.
    Serial { whole: bool },
    /// One query's report of outfmt 16 or 15 (`json`), written by an object stream of NCBI's on
    /// the report's stream inside the `try` of `BlastXML2_FormatReport` / `BlastJSON_FormatReport`
    /// (blastxml2_format.cpp:750-763, :877-890): as `Serial` (LOSAT inserts the report whole), but
    /// a failed insertion is caught there after `Close`'s 802.4: `XML format failed` (`JSON format
    /// failed`) is posted and the run goes on (`ReportWriteFailure::caught`); where only `Close`'s
    /// flush fails, the 802.4 alone, and the report ends well.
    Xml2Report { json: bool },
}

impl Site {
    /// The region of a run's report in format `number` (`-outfmt`'s): the tabular reports'
    /// insertions are all made while a `CBlastTabularInfo` lives (outfmt 6, 7, 10 and 20; a
    /// tabular report has no prolog), SAM's after `~CBlastFormat`; every other format starts in
    /// `Report` and declares its other regions (`ReportWrite::site`).
    pub fn of_format(number: i32) -> Site {
        use crate::blastinput::app::ReportFormat;
        match ReportFormat::from_number(number) {
            Some(
                ReportFormat::Tabular
                | ReportFormat::TabularWithComments
                | ReportFormat::CommaSeparatedValues
                | ReportFormat::CommaSeparatedValuesWithHeader,
            ) => Site::Tabular,
            Some(ReportFormat::Sam) => Site::Silent,
            _ => Site::Report,
        }
    }
}

/// What NCBI does when the report's stream failed (`ReportStream::finish`).
#[derive(Clone, Debug, PartialEq, Eq)]
enum Failure {
    /// The exception reached `CATCH_ALL` after the texts `posts` (each ending in a line feed).
    Thrown { posts: String },
    /// `std::terminate`: libstdc++'s verbose handler, then `abort()`.
    Terminate,
}

// NCBI reference (GCC 13 libstdc++): libsupc++/vterminate.cc `__verbose_terminate_handler`
// ```c++
// 	  fputs("terminate called after throwing an instance of '", stderr);
// 	  fputs(dem, stderr);
// 	  fputs("'\n", stderr);
//       ...
// 	      fputs("  what():  ", stderr);
// 	      fputs(w, stderr);
// 	      fputs("\n", stderr);
// ```
/// The text of `std::terminate` after the report's stream threw (`basic_ios::clear` of a stream
/// whose `exceptions(badbit)` is set and whose state has `badbit`), before `abort()`.
pub const TERMINATE_TEXT: &str = "terminate called after throwing an instance of \
'std::__ios_failure'\n  what():  basic_ios::clear: iostream error\n";

// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:317-324
// ```c++
// CObjectOStream::~CObjectOStream(void)
// {
//     try {
//         Close();
//         ResetLocalHooks();
//     }
//     NCBI_CATCH_X(4, "Cannot close serializing output stream")
// }
// ```
/// The post of `~CObjectOStream` whose `Close` flushed the failed report stream of `program`
/// (`[program]` is the application's name in NCBI's posts).
pub fn cannot_close_post(program: &str) -> String {
    format!(
        "Error: (802.4) [{program}] [Cannot close serializing output stream] Exception: \
         basic_ios::clear: iostream error\n"
    )
}

/// The error of an insertion into (or a flush of) the report's stream while NCBI's stream throws
/// (`exceptions(badbit)`) and a write of its device failed: the run stops there, as NCBI's
/// exception unwinds; `ReportStream::finish` turns it into NCBI's end (`Site`). Where NCBI's
/// formatter catches the exception (`caught`: `Site::XmlIteration`, whose post the stream has
/// made), the formatter stops that part of the report and goes on.
#[derive(Debug)]
pub struct ReportWriteFailure {
    kind: FailureKind,
}

/// What the formatter does with a `ReportWriteFailure`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FailureKind {
    /// NCBI's exception unwinds: the run stops.
    Thrown,
    /// NCBI's formatter catches it (`Site::XmlIteration`): that part of the report stops.
    Caught,
    /// The device write of a whole object (`Site::Serial { whole: true }`) failed: NCBI's end
    /// depends on the pieces of the object (`ReportWrite::resolve_object`).
    Undecided,
}

impl std::fmt::Display for ReportWriteFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("BLAST failed to write output")
    }
}

impl std::error::Error for ReportWriteFailure {}

/// Whether `error` is the report stream's failure (`ReportWriteFailure`).
pub fn is_report_write_failure(error: &std::io::Error) -> bool {
    error
        .get_ref()
        .is_some_and(|inner| inner.is::<ReportWriteFailure>())
}

/// Whether `error` is the report stream's failure in a region whose exception NCBI's formatter
/// catches (`ReportWriteFailure::caught`, `Site::XmlIteration`): the formatter stops that part
/// of the report (its post is made) and goes on.
pub fn is_caught_report_write_failure(error: &std::io::Error) -> bool {
    failure_kind(error) == Some(FailureKind::Caught)
}

/// Whether `error` is the failed device write of an object that LOSAT inserted whole
/// (`Site::Serial { whole: true }`): the writer hands NCBI's pieces of the object to
/// `ReportWrite::resolve_object`, which ends it as NCBI's does.
pub fn is_undecided_report_write_failure(error: &std::io::Error) -> bool {
    failure_kind(error) == Some(FailureKind::Undecided)
}

fn failure_kind(error: &std::io::Error) -> Option<FailureKind> {
    error
        .get_ref()
        .and_then(|inner| inner.downcast_ref::<ReportWriteFailure>())
        .map(|failure| failure.kind)
}

fn write_failure() -> std::io::Error {
    std::io::Error::other(ReportWriteFailure {
        kind: FailureKind::Thrown,
    })
}

fn caught_failure() -> std::io::Error {
    std::io::Error::other(ReportWriteFailure {
        kind: FailureKind::Caught,
    })
}

fn undecided_failure() -> std::io::Error {
    std::io::Error::other(ReportWriteFailure {
        kind: FailureKind::Undecided,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:745-764
// ```c++
// 	query_label = string("Q(")+NStr::NumericToString(index) +
// 	    string("/")+ NStr::NumericToString(data->GetNumQueries() );
// 	try{
//         s_BlastXMLAddIteration(bxmlout, data->GetAlignment(index), *seqloc,
//         ...
// 	}
// 	...
// 	catch(...){
// 	    ERR_POST(Error << "Failed s_BlastXMLAddIteration " << query_label  );
// 	    return;
// 	}
// ```
/// The post of `BlastXML_FormatReport`'s catch of a failed insertion of an `Iteration` of
/// `program` (`std::ios_base::failure` is no `CException`: no `what()` follows the label, which
/// NCBI writes without its `)`; each call formats one query, `Q(0/1`).
fn iteration_post(program: &str) -> String {
    format!("Error: [{program}] Failed s_BlastXMLAddIteration Q(0/1\n")
}

/// The command line's report stream (`-out` or standard output) as NCBI's writes it: NCBI's
/// buffer of the stream (`NcbiBuffer`: the `-out` file's `basic_filebuf`, standard output's glibc
/// buffer, line-buffered on a terminal) before the file or device, so that a write reaches the
/// device where NCBI's buffer writes (each `write` is one insertion of NCBI's, `insert`; a flush
/// of NCBI's is `ReportWrite::ncbi_flush`, LOSAT's flushes hold), and NCBI's stream state:
/// - a failed write of the device sets `badbit`; while the formatter lives (`CBlastFormat` sets
///   `exceptions(badbit)` and `~CBlastFormat` resets it, `finish`), that insertion and every later
///   one fail (`ReportWriteFailure`: libstdc++'s sentry sets `failbit` on a bad stream, which
///   throws), and what NCBI does then is the region's (`Site`, `finish`);
/// - after `~CBlastFormat` a bad stream drops the insertions, and a failed write is not reported
///   (the exit flushes C's `stdout` without a check).
///
/// `cerr` is tied to `cout`: a write to standard error writes what standard output holds first
/// (`tie_standard_output`, `CliStderr`); where that fails while the formatter lives, the post's
/// destructor throws (`std::terminate`). The `-out` file is not tied: a text on standard error
/// does not write it.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:317-323
// ```c++
// CBlastFormat::~CBlastFormat()
// {
//     try {
//         m_Outfile.exceptions(m_OrigExceptionMask);
//     } catch (...) {/*ignore exceptions*/}
//     m_Outfile.flush();
// }
// ```
pub struct ReportStream(std::sync::Arc<std::sync::Mutex<Core>>);

/// What `ReportStream` holds.
struct Core {
    device: Box<dyn std::io::Write + Send>,
    /// The bytes that NCBI's buffer holds.
    held: Vec<u8>,
    /// NCBI's buffer of the stream (`None`: none, each insertion written at once).
    buffer: Option<NcbiBuffer>,
    /// `badbit`: a write of the device failed.
    bad: bool,
    /// `exceptions(badbit)`: the formatter lives.
    throws: bool,
    /// The region of the insertions (`Site`).
    site: Site,
    /// NCBI's end after the failure (`finish`).
    failure: Option<Failure>,
    /// Whether the stream is `cout` (tied to `cerr`).
    standard_output: bool,
    /// Whether the format is 0-4 or 18 (`ClosedPipe::reported`).
    pipe_reported: bool,
    /// The application's name in NCBI's posts.
    program: &'static str,
    /// A write to a closed pipe failed (`ClosedPipe`).
    pipe: Option<ClosedPipe>,
    /// Whether the failure of the caught region being written was posted (`Site::XmlIteration`:
    /// one post per `Iteration`).
    region_posted: bool,
    /// Where NCBI's posts during the report go (`post`): the process's standard error, or a
    /// test's writer.
    diagnostics: Option<Box<dyn std::io::Write + Send>>,
    /// The bytes inserted so far (the stream's position).
    inserted: u64,
    /// The bytes that the device took.
    landed: u64,
    /// The failed device write of a whole object (`Site::Serial { whole: true }`), until
    /// `resolve_object`: NCBI's buffer before the object and the stream's position there.
    undecided: Option<(NcbiBuffer, u64)>,
    /// NCBI's buffer and the stream's position where the object being written started (the
    /// declaration of its region, `set_site`): LOSAT may insert a whole object in several writes
    /// (the archive's), all of which before the failed one reached the device.
    object_start: Option<(NcbiBuffer, u64)>,
}

/// Writes `data` to `device`, counting the bytes it takes in `landed` (a write that stops at a
/// file's size limit takes part of the bytes, the next fails).
fn write_device(
    device: &mut dyn std::io::Write,
    landed: &mut u64,
    mut data: &[u8],
) -> std::io::Result<()> {
    while !data.is_empty() {
        match device.write(data) {
            Ok(0) => return Err(std::io::ErrorKind::WriteZero.into()),
            Ok(n) => {
                *landed += n as u64;
                data = &data[n..];
            }
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/connect/ncbi_socket.c:866-873 (NCBI keeps the default action
// of SIGPIPE unless SOCK starts, which only the usage report does, after the run)
// ```c
// #ifdef NCBI_OS_UNIX
//         if (s_SigPipeIgnored == eOff) {
//             ...
//                 signal(SIGPIPE, SIG_IGN);
// ```
/// A write of the report to a pipe whose reader closed it, where NCBI is ended by SIGPIPE (exit
/// status 141, no text, every format) and LOSAT cannot be (the Rust runtime ignores SIGPIPE;
/// approved exceptions 3 and 5 of `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`). Until the Owner's answer
/// to line RG's question T6-Q1 the run ends as before this stream (item T6c): the formats 0-4 and
/// 18 with "BLAST failed to write output" (exit 6), every other format with the error of the
/// write (`Broken pipe`, exit 1), BLASTX's and RPSTBLASTN's with their own (outfmt 6, 7, 10 and 20
/// `std::terminate`'s text and an abort, every other format "BLAST failed to write output").
#[derive(Clone, Copy, Debug)]
struct ClosedPipe {
    /// Whether the format is 0-4 or 18.
    reported: bool,
    /// The error of the failed write (`EPIPE`).
    raw: Option<i32>,
}

impl ClosedPipe {
    fn error(self) -> std::io::Error {
        self.raw.map_or_else(
            || std::io::ErrorKind::BrokenPipe.into(),
            std::io::Error::from_raw_os_error,
        )
    }
}

impl Core {
    /// One insertion (`NcbiBuffer::sputn`): what NCBI's buffer writes of the bytes held and
    /// `buf` is written, the first of them, and it holds the rest (at most a buffer's bytes: the
    /// insertion is not copied before it is written).
    fn insert(&mut self, buf: &[u8]) -> std::io::Result<()> {
        // An empty write is no insertion of NCBI's (a part of a report that NCBI does not write).
        if buf.is_empty() {
            return Ok(());
        }
        if self.bad {
            // An object written to the bad stream after `~CBlastFormat` (the error archive of
            // `PrintErrorArchive`): its buffer's first write fails without an exception mask, a
            // write fault at the buffer's position then (`write_fault`).
            if !self.throws
                && self.pipe.is_none()
                && self.site == (Site::Serial { whole: true })
                && self.undecided.is_none()
            {
                if let Some(buffer) = self.buffer {
                    self.undecided = Some((buffer, self.inserted));
                    return Err(undecided_failure());
                }
            }
            return self.after_failure();
        }
        let before = self.buffer;
        let start = self.inserted;
        self.inserted += buf.len() as u64;
        let kept = match &mut self.buffer {
            Some(buffer) => {
                buffer.sputn(buf.len(), buf);
                buffer.held()
            }
            None => 0,
        };
        let written = self.held.len() + buf.len() - kept;
        if written == 0 {
            self.held.extend_from_slice(buf);
            return Ok(());
        }
        // NCBI's buffer writes all that it held, then the insertion's first bytes.
        let from_buf = written - self.held.len();
        let Core {
            device,
            landed,
            held,
            ..
        } = self;
        let result = write_device(&mut **device, landed, held)
            .and_then(|()| write_device(&mut **device, landed, &buf[..from_buf]))
            .and_then(|()| device.flush());
        self.held.clear();
        match result {
            Ok(()) => {
                self.held.extend_from_slice(&buf[from_buf..]);
                Ok(())
            }
            Err(error) => {
                if matches!(
                    self.site,
                    Site::Serial { whole: true } | Site::Xml2Report { .. }
                ) && !self.standard_output
                    && self.throws
                    && error.kind() != std::io::ErrorKind::BrokenPipe
                {
                    if let Some(before) = before {
                        self.bad = true;
                        self.undecided = Some(self.object_start.unwrap_or((before, start)));
                        return Err(undecided_failure());
                    }
                }
                self.failed_write(error, false)
            }
        }
    }

    /// A flush of NCBI's (`flush`, `endl`, `cerr`'s tie): the buffer writes what it holds.
    fn flush(&mut self) -> std::io::Result<()> {
        if self.bad {
            return self.after_failure();
        }
        if let Some(buffer) = &mut self.buffer {
            buffer.flush();
        }
        if !self.held.is_empty() {
            let Core {
                device,
                landed,
                held,
                ..
            } = self;
            let result = write_device(&mut **device, landed, held).and_then(|()| device.flush());
            self.held.clear();
            if let Err(error) = result {
                return self.failed_write(error, true);
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-383
    // ```c++
    // CObjectOStream::TFailFlags CObjectOStream::SetFailFlags(TFailFlags flags,
    //                                                         const char* message)
    // {
    //     TFailFlags old = m_Fail;
    //     m_Fail |= flags;
    //     if ( !old && flags ) {
    //         // first fail
    //         ERR_POST_X(5, "CObjectOStream: error at "<<
    //                    GetPosition()<<": "<<GetStackTrace() << ": " << message);
    //     }
    // ```
    /// The write fault of an object written to the bad stream after the formatter
    /// (`undecided` without the exception mask): `Close` posts 802.5 at `position`, the bytes that
    /// the object stream had put into its buffer at its first write, with an empty stack; the
    /// object is not written.
    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        if self.undecided.take().is_none() {
            return Ok(());
        }
        let program = self.program;
        self.post(&format!(
            "Error: (802.5) [{program}] CObjectOStream: error at byte {position}: stack is \
             empty: write fault\n"
        ))
    }

    /// NCBI's end of an object whose whole insertion failed at the device
    /// (`Site::Serial { whole: true }`, `undecided`), from NCBI's pieces of it (`pieces`): from
    /// NCBI's buffer before the object, a piece that the buffer writes past the bytes the device
    /// took fails inside the object (802.4 at `Close`, then `CATCH_ALL`); where every write inside
    /// the object succeeds, only `Close`'s flush fails (802.4 posted, the run goes on).
    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        let Some((mut buffer, mut position)) = self.undecided.take() else {
            return Ok(());
        };
        for &piece in pieces {
            position += piece as u64;
            let held = buffer.held();
            // The pieces' bytes are not read: a file's buffer counts them.
            buffer.sputn(piece, &[]);
            let writes = piece > 0 && buffer.held() < held + piece;
            if writes && position > self.landed {
                if let Site::Xml2Report { json } = self.site {
                    return self.xml2_caught(json);
                }
                self.failure.get_or_insert(Failure::Thrown {
                    posts: cannot_close_post(self.program),
                });
                return Err(write_failure());
            }
        }
        self.post(&cannot_close_post(self.program))
    }

    /// A write of the device that failed with `error`, for an insertion or for a flush
    /// (`flushing`): the stream is bad, and NCBI's end is the region's (`Site`).
    fn failed_write(&mut self, error: std::io::Error, flushing: bool) -> std::io::Result<()> {
        match error {
            // A closed pipe that only the writes after the formatter reach (`finish`) ends the
            // run as it would have (exception 5: LOSAT's completed writes did not reach it).
            error if error.kind() == std::io::ErrorKind::BrokenPipe && !self.throws => {
                self.bad = true;
                Ok(())
            }
            error if error.kind() == std::io::ErrorKind::BrokenPipe => {
                self.bad = true;
                self.pipe = Some(ClosedPipe {
                    reported: self.pipe_reported,
                    raw: error.raw_os_error(),
                });
                Err(error)
            }
            _ => {
                self.bad = true;
                if !self.throws {
                    return Ok(());
                }
                self.thrown(flushing)
            }
        }
    }

    /// The exception of a bad stream at an insertion (or a flush, `flushing`) while the formatter
    /// lives, in the region being written (`Site`): what NCBI does with it.
    fn thrown(&mut self, flushing: bool) -> std::io::Result<()> {
        let failure = match self.site {
            Site::Report => Failure::Thrown {
                posts: String::new(),
            },
            // glibc empties `stdout`'s buffer after a failed write (`new_do_write`), so `out2`'s
            // flush during the unwinding of an insertion's exception, which clears the stream's
            // state first, writes nothing and succeeds: no post. Its own flush at the end of the
            // alignments fails and posts, and the post's tie throws. `basic_filebuf` keeps the
            // bytes: the second flush fails too.
            Site::Alignments if self.standard_output && !flushing => Failure::Thrown {
                posts: String::new(),
            },
            Site::Alignments if self.standard_output => Failure::Terminate,
            Site::Alignments => Failure::Thrown {
                posts: cannot_close_post(self.program),
            },
            Site::Tabular => Failure::Terminate,
            Site::Silent => {
                self.throws = false;
                return Ok(());
            }
            // The catch's post, once for the `Iteration`; the run goes on.
            Site::XmlIteration => {
                if !std::mem::replace(&mut self.region_posted, true) {
                    self.post(&iteration_post(self.program))?;
                }
                return Err(caught_failure());
            }
            Site::Serial { .. } if self.standard_output => Failure::Terminate,
            // `Close`'s flush alone failed: its destructor posts and the run goes on.
            Site::Serial { .. } if flushing => {
                self.post(&cannot_close_post(self.program))?;
                return Ok(());
            }
            Site::Serial { .. } => Failure::Thrown {
                posts: cannot_close_post(self.program),
            },
            Site::Xml2Report { .. } if self.standard_output => Failure::Terminate,
            Site::Xml2Report { .. } if flushing => {
                self.post(&cannot_close_post(self.program))?;
                return Ok(());
            }
            Site::Xml2Report { json } => return self.xml2_caught(json),
        };
        // The first failure decides NCBI's end; a later one of the same exception (an insertion
        // after it, while the run unwinds) adds nothing.
        self.failure.get_or_insert(failure);
        Err(write_failure())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:884-890
    // ```c++
    // 	catch(CException &e){
    // 	    ERR_POST(Error << e.GetMsg() << e.what() );
    // 	    return;
    // 	}
    // 	catch(...){
    // 	    ERR_POST(Error << "JSON format failed" );
    // 	    return;
    // 	}
    // ```
    /// The catch of a report of outfmt 16 or 15 (`json`) whose insertion failed: `Close`'s 802.4
    /// while the exception unwinds, then the catch's post; the run goes on.
    fn xml2_caught(&mut self, json: bool) -> std::io::Result<()> {
        let format = if json { "JSON" } else { "XML" };
        self.post(&format!(
            "{}Error: [{}] {format} format failed\n",
            cannot_close_post(self.program),
            self.program
        ))?;
        Err(caught_failure())
    }

    /// An insertion or flush of a bad stream: it throws while the formatter lives (the region's
    /// end, `thrown`: the first insertion of a query's `Iteration` after a caught failure throws
    /// again).
    fn after_failure(&mut self) -> std::io::Result<()> {
        if let Some(pipe) = self.pipe {
            Err(pipe.error())
        } else if self.throws {
            self.thrown(false)
        } else {
            Ok(())
        }
    }

    /// A post of NCBI's on `cerr` (`ERR_POST`) while the report is written: `cerr` is tied to
    /// `cout`, so a report on standard output writes what it holds first; where that stream
    /// failed while the formatter lives, the tie throws inside the post's destructor
    /// (`std::terminate`) and the text is not written (`tie`).
    fn post(&mut self, text: &str) -> std::io::Result<()> {
        if self.standard_output {
            self.tie()?;
        }
        match &mut self.diagnostics {
            Some(diagnostics) => {
                diagnostics.write_all(text.as_bytes())?;
                diagnostics.flush()
            }
            None => {
                use std::io::Write;
                let mut stderr = std::io::stderr();
                stderr.write_all(text.as_bytes())?;
                stderr.flush()
            }
        }
    }

    /// The region of the insertions from here on (`ReportWrite::site`): each declaration starts a
    /// region (the next query's `Iteration` posts its failure again).
    fn set_site(&mut self, site: Site) {
        self.region_posted = false;
        self.site = site;
        self.object_start = match site {
            Site::Serial { whole: true } | Site::Xml2Report { .. } => {
                self.buffer.map(|buffer| (buffer, self.inserted))
            }
            _ => None,
        };
    }

    /// `cerr`'s tie: standard output writes what it holds before a text on standard error. A
    /// failed `cout` while the formatter lives throws there, inside the post's destructor
    /// (`std::terminate`); the text is not written.
    fn tie(&mut self) -> std::io::Result<()> {
        if self.bad && self.pipe.is_none() {
            // `cout.flush()` of a bad stream: its sentry throws while the formatter lives.
            if self.throws {
                self.failure = Some(Failure::Terminate);
                return Err(write_failure());
            }
            return Ok(());
        }
        let tied = self.flush();
        if self.pipe.is_some() {
            return tied;
        }
        if tied.is_err() {
            self.failure = Some(Failure::Terminate);
            return Err(write_failure());
        }
        Ok(())
    }

    /// `~CBlastFormat` and the application's `CATCH_ALL`: the stream no longer throws and
    /// writes what it holds (a failure unreported); the run's end is NCBI's for the failure, or
    /// `result`.
    fn finish(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        // A tabular report's bytes reach the device while its `CBlastTabularInfo` lives (the
        // destructor of each query's, and of the epilog's, flushes): what is held fails there.
        if self.site == Site::Tabular {
            let _ = self.flush();
        }
        self.throws = false;
        let _ = self.flush();
        if let Some(pipe) = self.pipe {
            let blastx = matches!(self.program, "blastx" | "rpstblastn");
            let failed = |exit: i32, message: &str| {
                Err(crate::cli::NativeError {
                    exit,
                    message: message.to_string(),
                }
                .into())
            };
            return match (pipe.reported, result) {
                (true, _) => failed(6, "BLAST failed to write output\n"),
                (false, Err(error)) => Err(error),
                (false, Ok(())) if blastx && self.site == Site::Tabular => {
                    failed(-6, TERMINATE_TEXT)
                }
                (false, Ok(())) if blastx => failed(6, "BLAST failed to write output\n"),
                (false, Ok(())) => Err(pipe.error().into()),
            };
        }
        match self.failure.take() {
            None => result,
            Some(Failure::Thrown { posts }) => Err(crate::cli::NativeError {
                exit: 6,
                message: posts + "BLAST failed to write output\n",
            }
            .into()),
            Some(Failure::Terminate) => Err(crate::cli::NativeError {
                exit: -6,
                message: TERMINATE_TEXT.to_string(),
            }
            .into()),
        }
    }
}

impl Drop for Core {
    fn drop(&mut self) {
        if !self.bad {
            let _ = self.device.write_all(&self.held);
            let _ = self.device.flush();
        }
    }
}

/// The run's report stream on standard output, which a write to standard error writes first
/// (`tie_standard_output`): standard output is the process's.
static TIED_STANDARD_OUTPUT: std::sync::Mutex<Option<std::sync::Weak<std::sync::Mutex<Core>>>> =
    std::sync::Mutex::new(None);

/// The application's name in NCBI's posts of the command line's run (`ReportStream::of_cli`), for
/// the posts of the streams that NCBI writes apart from the report (`write_side_file`).
static RUN_PROGRAM: std::sync::Mutex<Option<&'static str>> = std::sync::Mutex::new(None);

fn lock<T>(mutex: &std::sync::Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

impl ReportStream {
    fn new(
        device: Box<dyn std::io::Write + Send>,
        buffer: Option<NcbiBuffer>,
        standard_output: bool,
        program: &'static str,
        site: Site,
    ) -> Self {
        ReportStream(std::sync::Arc::new(std::sync::Mutex::new(Core {
            device,
            held: Vec::new(),
            buffer,
            bad: false,
            throws: true,
            site,
            failure: None,
            standard_output,
            program,
            pipe_reported: false,
            pipe: None,
            region_posted: false,
            diagnostics: None,
            inserted: 0,
            landed: 0,
            undecided: None,
            object_start: None,
        })))
    }

    /// The command line's report stream of `program` (the application's name in NCBI's posts)
    /// for the report of `-outfmt` `outfmt` (`Site::of_format`): the `-out` file `out_file`
    /// (NCBI's `basic_filebuf`, `NcbiStdio::File`: `BUFSIZ - 2` bytes held, the 8191st written
    /// with them) or, with `None`, standard output (glibc's buffer, `NcbiStdio::cout`, tied to
    /// standard error).
    pub fn of_cli(out_file: Option<std::fs::File>, program: &'static str, outfmt: &str) -> Self {
        let number = crate::blastinput::app::parse_formatting_string(outfmt)
            .map_or(0, |choice| choice.number);
        let stream = match out_file {
            Some(file) => Self::file(Box::new(file), program, Site::of_format(number)),
            None => Self::standard_output(
                Box::new(std::io::stdout()),
                program,
                Site::of_format(number),
            ),
        };
        lock(&stream.0).pipe_reported = matches!(number, 0..=4 | 18);
        *lock(&RUN_PROGRAM) = Some(program);
        stream
    }

    /// The `-out` file `inner`, its report starting in `site`.
    pub fn file(inner: Box<dyn std::io::Write + Send>, program: &'static str, site: Site) -> Self {
        Self::new(
            inner,
            Some(NcbiBuffer::new(NcbiStdio::File)),
            false,
            program,
            site,
        )
    }

    /// A report kept in memory (`inner`, TBLASTN's Web ABI v1): each insertion is written at once
    /// (no stream of NCBI's to model).
    pub fn memory(inner: Box<dyn std::io::Write + Send>, program: &'static str) -> Self {
        Self::new(inner, None, false, program, Site::Report)
    }

    /// Standard output, writing to `inner` (the process's standard output), tied to standard
    /// error (`tie_standard_output`).
    ///
    /// NCBI writes the report to `cout`; when the process starts with its standard output
    /// closed (`>&-`), its first write fails (outfmt 0: "BLAST failed to write output", exit
    /// 6; outfmt 6 and 7 abort). Rust's runtime opens `/dev/null` on a closed standard
    /// descriptor before `main`, which LOSAT cannot tell from a `/dev/null` that the caller
    /// opened (read and write, as Python's `subprocess.DEVNULL`), so the report is written
    /// there and the run succeeds (approved exception 6 of
    /// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, accepted by the maintainer in session S08b).
    pub fn standard_output(
        inner: Box<dyn std::io::Write + Send>,
        program: &'static str,
        site: Site,
    ) -> Self {
        // NCBI's Windows (MSVC) and macOS (libc++, BSD stdio) builds buffer otherwise: there
        // the exit status of a failing write is unverified (this model is their Linux builds').
        let report = Self::new(
            inner,
            Some(NcbiBuffer::new(NcbiStdio::cout())),
            true,
            program,
            site,
        );
        *lock(&TIED_STANDARD_OUTPUT) = Some(std::sync::Arc::downgrade(&report.0));
        report
    }

    /// The handle by which a formatter that writes into the stream declares NCBI's regions and
    /// flushes (`NcbiOstream::report`, `ReportWrite`).
    pub fn handle(&self) -> ReportHandle {
        ReportHandle(self.0.clone())
    }

    /// `~CBlastFormat` and the application's `CATCH_ALL`, after the run's `result`: the stream
    /// stops throwing and writes what it holds, unchecked; a failure of the stream while the
    /// formatter lived ends the run as NCBI's does (`Site`: the region's posts and "BLAST failed
    /// to write output", exit 6, or `std::terminate`'s text and an abort, `NativeError` of exit
    /// -6), whatever `result` is (the failure's error, `ReportWriteFailure`, or what the run made
    /// of it); otherwise `result`.
    pub fn finish(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        lock(&self.0).finish(result)
    }

    /// Whether a write of the stream's device failed (`badbit`).
    pub fn failed(&self) -> bool {
        lock(&self.0).bad
    }
}

impl std::io::Write for ReportStream {
    /// One insertion of NCBI's.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        lock(&self.0).insert(buf)?;
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        lock(&self.0).insert(buf)
    }

    /// LOSAT's flushes of the report hold (NCBI's are `ReportWrite::ncbi_flush`).
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// A formatter's hold on the command line's report stream (`ReportStream::handle`), by which it
/// declares NCBI's regions (`Site`) and flushes; the stream of a library caller or LOSAT Web has
/// none (`NcbiOstream::report`).
#[derive(Clone)]
pub struct ReportHandle(std::sync::Arc<std::sync::Mutex<Core>>);

impl ReportHandle {
    fn site(&self, site: Site) {
        lock(&self.0).set_site(site);
    }

    fn flush(&self) -> std::io::Result<()> {
        lock(&self.0).flush()
    }

    fn resolve_object(&self, pieces: &[usize]) -> std::io::Result<()> {
        lock(&self.0).resolve_object(pieces)
    }

    fn write_fault(&self, position: usize) -> std::io::Result<()> {
        lock(&self.0).write_fault(position)
    }

    /// The insertions from here on are made in NCBI's region `site` (`ReportWrite::site`, for a
    /// hold that keeps the regions in order with its bytes).
    pub(crate) fn declare(&self, site: Site) -> std::io::Result<()> {
        self.site(site);
        Ok(())
    }

    /// A flush of NCBI's (`ReportWrite::ncbi_flush`, for a hold).
    pub(crate) fn ncbi_flush(&self) -> std::io::Result<()> {
        self.flush()
    }

    /// `ReportWrite::resolve_object`, for a hold.
    pub(crate) fn resolve(&self, pieces: &[usize]) -> std::io::Result<()> {
        self.resolve_object(pieces)
    }

    /// `~CBlastFormat`'s reset of the stream's exceptions, before the writes of a format that
    /// NCBI makes after it (SAM's `~CSAM_Formatter`, a member of `CBlastFormat`).
    pub fn stop_throwing(&self) {
        lock(&self.0).throws = false;
    }

    /// `ReportStream::finish`.
    pub fn finish(&self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        lock(&self.0).finish(result)
    }
}

/// A writer of a report that declares NCBI's regions and flushes to the command line's report
/// stream behind it (`ReportStream`, through `NcbiOstream`); every other writer ignores them.
pub trait ReportWrite: std::io::Write {
    /// The insertions from here on are made in NCBI's region `site`.
    fn site(&mut self, _site: Site) -> std::io::Result<()> {
        Ok(())
    }

    /// A flush of NCBI's (`flush`, `endl`): NCBI's buffer writes what it holds.
    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }

    /// After `is_undecided_report_write_failure`: NCBI's pieces of the object that LOSAT inserted
    /// whole (`Site::Serial { whole: true }`), which decide NCBI's end (an error: the run stops).
    fn resolve_object(&mut self, _pieces: &[usize]) -> std::io::Result<()> {
        Ok(())
    }

    /// After `is_undecided_report_write_failure` of an object written once the formatter ended
    /// (the error archive): NCBI's write fault at `position` (`ReportStream`'s 802.5 post).
    fn write_fault(&mut self, _position: usize) -> std::io::Result<()> {
        Ok(())
    }
}

/// The command line's report stream itself (after the formatter: the error archive).
impl ReportWrite for ReportStream {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        lock(&self.0).set_site(site);
        Ok(())
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        lock(&self.0).flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        lock(&self.0).resolve_object(pieces)
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        lock(&self.0).write_fault(position)
    }
}

impl<W: ReportWrite + ?Sized> ReportWrite for &mut W {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        (**self).site(site)
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        (**self).ncbi_flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        (**self).resolve_object(pieces)
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        (**self).write_fault(position)
    }
}

/// LOSAT's buffer before the stream: it writes what it holds first.
impl<W: ReportWrite> ReportWrite for std::io::BufWriter<W> {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        std::io::Write::flush(self)?;
        self.get_mut().site(site)
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        std::io::Write::flush(self)?;
        self.get_mut().ncbi_flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        self.get_mut().resolve_object(pieces)
    }
}

/// A buffer of the report's bytes (a piece formatted before it is written).
impl ReportWrite for Vec<u8> {}

/// A writer that takes no region or flush of NCBI's (`ReportWrite`'s defaults): the destination of
/// a report part that declares none (`query_warnings::write_marked`).
pub struct Unregioned<'a>(pub &'a mut dyn std::io::Write);

impl std::io::Write for Unregioned<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write(buf)
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.0.write_all(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

impl ReportWrite for Unregioned<'_> {}

/// What NCBI's formatter does at an offset of a part of the report that LOSAT formats into a
/// buffer before it writes it (`Cuts`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cut {
    /// One insertion of NCBI's ends here and the next starts.
    Insertion,
    /// The insertions from here on are made in NCBI's region `site` (`ReportWrite::site`); the
    /// insertion before ends here.
    Site(Site),
    /// NCBI flushes the stream here (`endl`, `ReportWrite::ncbi_flush`); the insertion before ends
    /// here.
    Flush,
}

/// The places in a part of the report formatted into a buffer where NCBI's insertions end and
/// where NCBI's formatter enters a region or flushes (`Cut`), at offsets of the buffer in
/// increasing order: the buffer is written as NCBI inserts it (`CutCursor::write`), each insertion
/// one write of the stream (`insert`), so that the stream's buffer (`ReportStream`) writes and
/// fails where NCBI's does. A buffer without cuts is one insertion.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Cuts(Vec<(usize, Cut)>);

impl Cuts {
    /// `cut` at `at`, after the cuts so far.
    pub fn push(&mut self, at: usize, cut: Cut) {
        debug_assert!(self.0.last().is_none_or(|&(last, _)| last <= at));
        self.0.push((at, cut));
    }

    /// The cuts of `other`, after these.
    pub fn extend(&mut self, other: Cuts) {
        for (at, cut) in other.0 {
            self.push(at, cut);
        }
    }

    /// A cursor that writes the buffer from its start.
    pub fn cursor(&self) -> CutCursor<'_> {
        CutCursor {
            cuts: &self.0,
            next: 0,
            region: (0, Site::Report),
            skipping: false,
        }
    }
}

/// NCBI's insertions of `report`, one query's report of outfmt 16 or 15 (`site`,
/// `Site::Xml2Report`), into the report's stream: the writes of its object stream's buffer, taken
/// back from the report's bytes (`out_buffer::xml_pieces`, `json_pieces`).
pub(crate) fn report_pieces(site: Site, report: &[u8]) -> Vec<usize> {
    match site {
        Site::Xml2Report { json: true } => crate::report::serial::out_buffer::json_pieces(report),
        _ => crate::report::serial::out_buffer::xml_pieces(report, |tag| {
            crate::report::blast_xml2::is_string_element(tag)
        }),
    }
}

/// The writing of a buffer with its `Cuts`, part by part (`write`).
pub struct CutCursor<'c> {
    cuts: &'c [(usize, Cut)],
    next: usize,
    /// Where the region being written starts, and its site.
    region: (usize, Site),
    /// The rest of the region was done where NCBI's pieces decided its end (`insert`).
    skipping: bool,
}

impl CutCursor<'_> {
    /// Writes `report[from..to]` to `writer`, the part after the parts written before: an
    /// insertion up to each cut in it, and each cut's region or flush (those at `to` too, which
    /// come before what follows the part).
    pub fn write(
        &mut self,
        writer: &mut dyn ReportWrite,
        report: &[u8],
        from: usize,
        to: usize,
    ) -> std::io::Result<()> {
        let mut at = from;
        while let Some(&(cut_at, cut)) = self
            .cuts
            .get(self.next)
            .filter(|&&(cut_at, _)| cut_at <= to)
        {
            let cut_at = cut_at.max(at);
            if cut_at > at {
                self.insert(writer, report, at, cut_at)?;
                at = cut_at;
            }
            if self.skipping && !matches!(cut, Cut::Site(_)) {
                self.next += 1;
                continue;
            }
            self.skipping = false;
            match cut {
                Cut::Insertion => {}
                Cut::Site(site) => {
                    writer.site(site)?;
                    self.region = (cut_at, site);
                }
                Cut::Flush => writer.ncbi_flush()?,
            }
            self.next += 1;
        }
        if to > at {
            self.insert(writer, report, at, to)?;
        }
        Ok(())
    }

    /// One insertion, `report[from..to]`. A report of outfmt 16 or 15 is inserted whole
    /// (`Site::Xml2Report`): where its write fails at the device, NCBI's pieces of the report
    /// (`out_buffer::xml_pieces`, `json_pieces`, from its bytes) decide NCBI's end
    /// (`ReportWrite::resolve_object`); where only its `Close` fails, the rest of the region
    /// (`Close`'s flush) is done.
    fn insert(
        &mut self,
        writer: &mut dyn ReportWrite,
        report: &[u8],
        from: usize,
        to: usize,
    ) -> std::io::Result<()> {
        if self.skipping {
            return Ok(());
        }
        match writer.write_all(&report[from..to]) {
            Err(error) if is_undecided_report_write_failure(&error) => {
                let (start, site) = self.region;
                let end = self.cuts[self.next..]
                    .iter()
                    .find(|(_, cut)| matches!(cut, Cut::Site(_)))
                    .map_or(report.len(), |&(at, _)| at);
                let pieces = report_pieces(site, &report[start.min(end)..end]);
                writer.resolve_object(&pieces)?;
                self.skipping = true;
                Ok(())
            }
            written => written,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1911-1931
// ```c++
// void CDisplaySeqalign::DisplaySeqalign(CNcbiOstream& out)
// {
//     CSeq_align_set actual_aln_list;
//     CAlignFormatUtil::ExtractSeqalignSetFromDiscSegs(actual_aln_list,
//                                                      *m_SeqalignSetRef);
//     if (actual_aln_list.Get().empty()){
//         return;
//     }
//     ...
//     unique_ptr<CObjectOStream> out2(CObjectOStream::Open(eSerial_AsnText, out));
// ```
/// NCBI's `DisplaySeqalign` of a query's alignments: from its first alignment (none shown: no
/// call that writes) they are written in its region (`Site::Alignments`), and at their end the
/// destructor of its `out2` flushes the stream (`CObjectOStream::Close`, `objostr.cpp:338-354`),
/// before the query's footer.
#[derive(Default)]
pub(crate) struct DisplaySeqalign {
    open: bool,
}

impl DisplaySeqalign {
    /// Before an alignment: the region starts at the first.
    pub(crate) fn enter<W: ReportWrite + ?Sized>(&mut self, writer: &mut W) -> std::io::Result<()> {
        if !self.open {
            self.open = true;
            writer.site(Site::Alignments)?;
        }
        Ok(())
    }

    /// After the query's last alignment: `out2`'s flush, then the report's region.
    pub(crate) fn leave<W: ReportWrite + ?Sized>(&mut self, writer: &mut W) -> std::io::Result<()> {
        if self.open {
            self.open = false;
            writer.ncbi_flush()?;
            writer.site(Site::Report)?;
        }
        Ok(())
    }
}

/// `out << value` in NCBI's report: one insertion (a string or a number) of `bytes`, one write
/// call, which the command line's report stream takes as one (`ReportStream`); other writers just
/// take the bytes. A string that LOSAT streams is `insert_streamed`.
pub fn insert(out: &mut (impl std::io::Write + ?Sized), bytes: &[u8]) -> std::io::Result<()> {
    out.write_all(bytes)
}

/// One insertion of NCBI's of `len` bytes that LOSAT streams instead of holding (the organism
/// report, whose runs of about 4 GiB of spaces `tax_report` keeps as counts): `fill` gives the bytes to
/// the `Insertion`, which writes them in pieces that NCBI's buffers (`NcbiBuffer`) write as they
/// write the whole. Each piece is longer than a buffer (`basic_filebuf` writes it with what it
/// holds and keeps nothing; glibc keeps `(held + n) % block`, which adds up over the pieces as
/// over the whole), and each but the last ends in `NCBI_BUFSIZ` bytes without a newline (so a
/// line-buffered `cout` keeps its whole remainder there, as within the whole); the last is the
/// rest. A string without such a place is written in one piece.
pub fn insert_streamed<W: std::io::Write + ?Sized>(
    out: &mut W,
    len: u64,
    fill: impl FnOnce(&mut Insertion<'_, W>) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let mut insertion = Insertion {
        out,
        piece: Vec::new(),
        left: len,
        since_newline: 0,
    };
    fill(&mut insertion)?;
    debug_assert_eq!(insertion.left, insertion.piece.len() as u64);
    insertion.out.write_all(&insertion.piece)
}

/// The bytes of one insertion that LOSAT streams (`insert_streamed`).
pub struct Insertion<'w, W: std::io::Write + ?Sized> {
    out: &'w mut W,
    /// The piece being made.
    piece: Vec<u8>,
    /// The bytes of the insertion not written yet, `piece`'s among them.
    left: u64,
    /// How many bytes follow the last newline.
    since_newline: u64,
}

impl<W: std::io::Write + ?Sized> Insertion<'_, W> {
    /// The length from which a piece is written where it may end.
    const PIECE: usize = 1 << 16;

    /// The next bytes of the insertion.
    pub fn bytes(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        self.piece.extend_from_slice(bytes);
        self.since_newline = match bytes.iter().rposition(|&byte| byte == b'\n') {
            Some(at) => (bytes.len() - at - 1) as u64,
            None => self.since_newline + bytes.len() as u64,
        };
        self.end_piece()
    }

    /// The next `count` bytes of the insertion, spaces.
    pub fn spaces(&mut self, mut count: u64) -> std::io::Result<()> {
        while count > 0 {
            let n = count.min(Self::PIECE as u64);
            self.piece.resize(self.piece.len() + n as usize, b' ');
            self.since_newline += n;
            count -= n;
            self.end_piece()?;
        }
        Ok(())
    }

    /// Writes the piece where it may end (`insert_streamed`) once it is `PIECE` bytes long.
    fn end_piece(&mut self) -> std::io::Result<()> {
        let n = self.piece.len() as u64;
        if n >= Self::PIECE as u64
            && self.since_newline >= NCBI_BUFSIZ as u64
            && self.left - n > NCBI_BUFSIZ as u64
        {
            self.out.write_all(&self.piece)?;
            self.left -= n;
            self.piece.clear();
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
// ```c
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
//     : m_FormatType(format_type), ..., m_Outfile(outfile),
// ```
/// The stream of a report as NCBI's formatter writes to it (`CBlastFormat::m_Outfile`), with
/// nothing of a program's before it: each write is one insertion of NCBI's (`out << value`,
/// `insert`, `insert_streamed`), passed on as one write. A writer becomes one only through the
/// constructors below, one per destination, each saying what comes after it: they are the only
/// places where a buffer or hold could come between the report's insertions and NCBI's buffer
/// (`ReportStream`, which counts each write as one insertion). A format's sink (`sink`: the
/// command line's `ReportStream`, `OutputSink::Report`, or the library's and LOSAT Web's
/// writers), blastx's command line (`report`), and the hold of a library call (`held`). Only the
/// command line's stream takes NCBI's regions and flushes (`ReportWrite`).
pub struct NcbiOstream<'a> {
    out: Box<dyn std::io::Write + 'a>,
    report: Option<ReportHandle>,
    /// A hold before the report stream (`before_report`), which takes the regions and flushes in
    /// order with its bytes.
    hold: Option<Box<dyn ReportWrite + 'a>>,
}

impl<'a> NcbiOstream<'a> {
    fn new(out: Box<dyn std::io::Write + 'a>) -> Self {
        NcbiOstream {
            out,
            report: None,
            hold: None,
        }
    }

    /// A format's sink opened as its report's stream (`OutputSink::open`): the command line's
    /// `ReportStream` (`OutputSink::Report`); or the caller's writer as it is, which a library
    /// caller or LOSAT Web gives (no stream of NCBI's to model); or, for a library call's
    /// `Stdout` and `File` sinks, LOSAT's buffered standard output or file (no stream of NCBI's
    /// either: the report's bytes are NCBI's, the times at which they reach the device are not).
    pub(crate) fn sink(
        sink: &'a mut crate::api::local_blast::OutputSink<'_>,
    ) -> std::io::Result<Self> {
        let report = sink.report_handle();
        let out = sink.open()?;
        Ok(NcbiOstream {
            out,
            report,
            hold: None,
        })
    }

    /// The command line's report stream `stream` of blastx and rpstblastn (`native::run_inner`),
    /// every format's.
    pub(crate) fn report(stream: ReportStream) -> NcbiOstream<'static> {
        let report = stream.handle();
        NcbiOstream {
            out: Box::new(stream),
            report: Some(report),
            hold: None,
        }
    }

    /// `hold`, a hold of blastx's command line before its report stream `report`
    /// (`blast_xml::HeldOutput`): it takes the bytes, and NCBI's regions and flushes in order with
    /// them (it passes them to `report`).
    pub(crate) fn before_report<H: ReportWrite + Clone + 'a>(
        hold: H,
        report: ReportHandle,
    ) -> Self {
        NcbiOstream {
            out: Box::new(hold.clone()),
            report: Some(report),
            hold: Some(Box::new(hold)),
        }
    }

    /// The command line's report stream behind this one (`report`, `sink`).
    pub(crate) fn report_handle(&self) -> Option<&ReportHandle> {
        self.report.as_ref()
    }

    /// A format's stream (`sink`'s) behind the hold of a library call (`local::report_batches`,
    /// `blast_xml::HeldFormats`), which passes the writes and flushes on in their order (at the
    /// release, while it holds). The caller's writer has no stream of NCBI's (`sink`), so no buffer
    /// of NCBI's comes after the hold.
    pub(crate) fn held(
        stream: crate::algorithm::blastx::HeldFormatWriter<NcbiOstream<'a>>,
    ) -> Self {
        Self::new(Box::new(stream))
    }
}

impl std::io::Write for NcbiOstream<'_> {
    /// One insertion, written with one write of the stream.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.out.write_all(buf)?;
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.out.write_all(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.out.flush()
    }
}

impl ReportWrite for NcbiOstream<'_> {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        if let Some(hold) = &mut self.hold {
            return hold.site(site);
        }
        if let Some(report) = &self.report {
            report.site(site);
        }
        Ok(())
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        if let Some(hold) = &mut self.hold {
            return hold.ncbi_flush();
        }
        match &self.report {
            Some(report) => report.flush(),
            None => Ok(()),
        }
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        match &self.report {
            Some(report) => report.resolve_object(pieces),
            None => Ok(()),
        }
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        match &self.report {
            Some(report) => report.write_fault(position),
            None => Ok(()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-778
// ```c++
// BlastXML2_FormatReport(const IBlastXML2ReportData* data, string file_name)
// {
// 	blastxml2::CBlastOutput2 bxmlout;
// 		CNcbiOfstream out_stream;
// 		out_stream.open(file_name.c_str(), IOS_BASE::out);
// 		...
// 		s_FillBlastOutput(bxmlout, data);
// 		s_WriteXML2Object(bxmlout, &out_stream);
// }
// ```
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:893-896 (`COStreamBuffer::FlushBuffer`)
// ```c++
//         if ( !m_Output.write(m_Buffer, count) ) {
//             m_Error = "write fault";
//             NCBI_THROW(CIOException,eWrite,m_Error);
//         }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:339-354, :385-395
// ```c++
// void CObjectOStream::Close(void)
// {
//     if (m_Fail != fNotOpen) {
//         try {
//             DefaultFlush();
//             ...
//         catch (CException& exc) {
//             if ( InGoodState() )
//                 RethrowError(fWriteError, "cannot close output stream",exc);
//         }
// ...
// bool CObjectOStream::InGoodState(void)
// {
//     ...
//     else if ( m_Output.fail() ) {
//         // IO exception thrown without setting fail flag
//         SetFailFlags(fWriteError, m_Output.GetError());
// ```
/// Writes `document`, the report of a side file of outfmt 13 or 14, to `file` as NCBI's
/// `CNcbiOfstream` (no exception mask; `basic_filebuf`, `NcbiStdio::File`) gets it from the buffer
/// of its object stream: `pieces` are that buffer's writes (`out_buffer::xml_pieces`,
/// `json_pieces`; `back_limit` its `m_BackLimit`), each one insertion of the stream.
///
/// A piece whose write of the file fails is the buffer's "write fault" (`CIOException::eWrite`):
/// the stream's `Close` posts 802.5 at the buffer's position then (the bytes it had written and
/// those it held, the piece's and the back limit's) with an empty stack, and the exception ends
/// the run at `CATCH_ALL`, which swallows that code: exit status 0, no later report and no master
/// file (an error of exit status 0, its text posted here). Where only `Close`'s flush fails, its
/// chain of posts ends in the destructor's 802.4 and the run goes on. Without the command line's
/// run (a library call writes its side files without a stream of NCBI's) a failed write is the
/// old error, `blast_xml2::framing::side_file_write_error`.
pub fn write_side_file(
    mut file: std::fs::File,
    document: &[u8],
    pieces: &[usize],
    back_limit: usize,
) -> anyhow::Result<()> {
    use std::io::Write;
    let Some(program) = *lock(&RUN_PROGRAM) else {
        return file
            .write_all(document)
            .and_then(|()| file.flush())
            .map_err(|_| crate::report::blast_xml2::framing::side_file_write_error());
    };
    let mut buffer = NcbiBuffer::new(NcbiStdio::File);
    let mut landed = 0u64;
    // The bytes before the held ones: what `basic_filebuf` wrote.
    let mut start = 0usize;
    let mut end = 0usize;
    let last = pieces.len().saturating_sub(1);
    for (k, &piece) in pieces.iter().enumerate() {
        end += piece;
        buffer.sputn(piece, &[]);
        if buffer.held() == 0 && piece > 0 {
            if write_device(&mut file, &mut landed, &document[start..end]).is_err() {
                // `EndOfWrite`'s write leaves nothing in the object stream's buffer.
                let leave = if k == last { 0 } else { back_limit };
                let post = format!(
                    "Error: (802.5) [{program}] CObjectOStream: error at byte {}: stack is empty: \
                     write fault\n",
                    end + leave
                );
                return Err(crate::cli::NativeError {
                    exit: 0,
                    message: post,
                }
                .into());
            }
            start = end;
        }
    }
    if write_device(&mut file, &mut landed, &document[start..end])
        .and_then(|()| file.flush())
        .is_err()
    {
        // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:899-907 (`COStreamBuffer::Flush`)
        // ```c++
        //     STemporarilyClearStreamState state(m_Output);
        //     FlushBuffer();
        //     if ( !m_Output.flush() ) {
        //         NCBI_THROW(CIOException,eFlush,"COStreamBuffer::Flush: failed");
        //     }
        // ```
        let position = end;
        let chain = format!(
            "Error: (802.5) [{program}] CObjectOStream: error at byte {position}: stack is empty: \
             cannot close output stream\n\
             Error: (CIOException::eFlush) [{program}] COStreamBuffer::Flush: failed\n\
             Error: (CSerialException::eIoError) [{program}] byte {position}: cannot close output \
             stream\n\
             Error: (802.4) [{program}] Cannot close serializing output stream \
             (CSerialException::eIoError) byte {position}: cannot close output stream\n"
        );
        tie_standard_output()?;
        let mut stderr = std::io::stderr();
        stderr.write_all(chain.as_bytes())?;
        stderr.flush()?;
    }
    Ok(())
}

/// `cerr.tie()`: before a write to standard error, the run's report on standard output writes
/// what it holds (`ReportStream::standard_output`); nothing for a report in `-out`.
pub fn tie_standard_output() -> std::io::Result<()> {
    let tied = lock(&TIED_STANDARD_OUTPUT)
        .as_ref()
        .and_then(std::sync::Weak::upgrade);
    match tied {
        Some(report) => lock(&report).tie(),
        None => Ok(()),
    }
}

/// The command line's standard error, tied to the organism report's standard output
/// (`tie_standard_output`): the process's (`process`), so that what standard output writes before
/// a text is written before it where the text is made. A whole message is written with one
/// `write_all` of the stream, which `std::io::Stderr` writes under one lock of standard error, so
/// that no other thread's message comes inside it. An empty write is no write of NCBI's to `cerr`
/// (a query without warnings, whose `ERR_POST` NCBI does not make): it writes nothing of the
/// report.
pub struct CliStderr(std::io::Stderr);

impl CliStderr {
    /// The process's standard error.
    pub fn process() -> Self {
        CliStderr(std::io::stderr())
    }
}

impl std::io::Write for CliStderr {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        tie_standard_output()?;
        self.0.write(buf)
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        if buf.is_empty() {
            return Ok(());
        }
        tie_standard_output()?;
        self.0.write_all(buf)
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        self.write_all(args.to_string().as_bytes())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;

    /// A device that takes nothing (`/dev/full`), counting the bytes offered to it.
    struct Full(std::sync::Arc<std::sync::atomic::AtomicUsize>);

    impl std::io::Write for Full {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0
                .fetch_add(buf.len(), std::sync::atomic::Ordering::Relaxed);
            Err(std::io::Error::other("full"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn the_organism_reports_file_buffer_writes_at_its_8191st_byte() {
        // probe `T/t2/fix1` `f1.d5.*` and `stdout_probe.sh`: NCBI exits 0 for a report that its
        // buffer holds to the end, 6 where the buffer writes.
        use std::io::Write;
        let offered = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let mut buffer =
            ReportStream::file(Box::new(Full(offered.clone())), "blastn", Site::Report);
        buffer.write_all(&[b'x'; 8190]).unwrap();
        buffer.flush().unwrap();
        assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 0);
        assert!(buffer.write_all(b"x").is_err());
        // The bytes held are offered first (then the insertion's, which a full device does not
        // reach): the write comes at the 8191st byte.
        assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 8190);
        // What a buffer holds at the end is written, a failure unreported.
        let quiet = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let mut buffer = ReportStream::file(Box::new(Full(quiet.clone())), "blastn", Site::Report);
        buffer.write_all(b"held").unwrap();
        drop(buffer);
        assert_eq!(quiet.load(std::sync::atomic::Ordering::Relaxed), 4);
    }

    /// A device that takes everything, keeping the size of each write.
    struct Writes(std::sync::Arc<std::sync::Mutex<Vec<usize>>>);

    impl std::io::Write for Writes {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            lock(&self.0).push(buf.len());
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn the_organism_reports_standard_output_writes_where_glibc_writes() {
        // glibc 2.39 `_IO_new_file_xsputn` on a file of `st_blksize` 4096
        // (`T/t2/fix2/glibc_first_write.log`): the writes end at multiples of the block and the
        // rest stays held; the first `fwrite` has no buffer yet; a flush empties the buffer.
        use std::io::Write;
        let run = |pieces: &[usize]| {
            let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut after = Vec::new();
            let mut buffer = ReportStream::new(
                Box::new(Writes(writes.clone())),
                Some(NcbiBuffer::new(NcbiStdio::Cout {
                    block: 4096,
                    lines: false,
                })),
                false,
                "blastn",
                Site::Report,
            );
            for &piece in pieces {
                match piece {
                    0 => lock(&buffer.0).flush().unwrap(),
                    n => buffer.write_all(&vec![b'x'; n]).unwrap(),
                }
                after.push(lock(&writes).iter().sum::<usize>());
            }
            drop(buffer);
            let total = lock(&writes).iter().sum::<usize>();
            (after, total)
        };
        // `4000 200 4000`: 4,096 written at the second, 8,192 at the third, 8 held to the end.
        assert_eq!(run(&[4000, 200, 4000]), (vec![0, 4096, 8192], 8200));
        assert_eq!(run(&[4096, 10]), (vec![4096, 4096], 4106));
        assert_eq!(run(&[4095, 1, 10]), (vec![0, 0, 4096], 4106));
        assert_eq!(run(&[8192, 1]), (vec![8192, 8192], 8193));
        assert_eq!(run(&[100, 4096]), (vec![0, 4096], 4196));
        // `100 f 4096 1` (0: the flush of a write to `cerr`, which is tied to `cout`).
        assert_eq!(run(&[100, 0, 4096, 1]), (vec![0, 100, 100, 4196], 4197));
        assert_eq!(run(&[8000, 300]), (vec![4096, 8192], 8300));
    }

    /// The bytes that `buffer` has written to its device after each step of `steps`
    /// (`Ok(bytes)`: one insertion; `Err(())`: a flush, as before a write to `cerr`), and in all.
    fn written_after(
        stream: NcbiStdio,
        steps: &[Result<Vec<u8>, ()>],
    ) -> (Vec<usize>, usize, usize) {
        use std::io::Write;
        let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut after = Vec::new();
        let mut buffer = ReportStream::new(
            Box::new(Writes(writes.clone())),
            Some(NcbiBuffer::new(stream)),
            false,
            "blastn",
            Site::Report,
        );
        for step in steps {
            match step {
                Ok(bytes) => buffer.write_all(bytes).unwrap(),
                Err(()) => lock(&buffer.0).flush().unwrap(),
            }
            after.push(lock(&writes).iter().sum::<usize>());
        }
        drop(buffer);
        let writes = lock(&writes);
        (after, writes.iter().sum(), writes.len())
    }

    #[test]
    fn a_terminals_standard_output_writes_at_its_newlines_as_glibc_does() {
        // glibc 2.39 on a pty of `st_blksize` 1024 (`T/t2/fix3/glibc_tty_lines.log`): `N` is an
        // `fwrite` of N bytes, `N:K` with a newline at K, `f` a flush.
        let tty = NcbiStdio::Cout {
            block: 1024,
            lines: true,
        };
        let run = |steps: &str| {
            let steps = steps
                .split(' ')
                .map(|step| match step {
                    "f" => Err(()),
                    _ => {
                        let mut parts = step.split(':').map(|n| n.parse::<usize>().unwrap());
                        let mut bytes = vec![b'x'; parts.next().unwrap()];
                        parts.for_each(|at| bytes[at] = b'\n');
                        Ok(bytes)
                    }
                })
                .collect::<Vec<_>>();
            written_after(tty, &steps).0
        };
        assert_eq!(run("100"), [0]);
        assert_eq!(run("100 50:49"), [0, 150]);
        assert_eq!(run("100 50:10"), [0, 111]);
        assert_eq!(run("3000"), [2048]);
        assert_eq!(run("3000:2500"), [2501]);
        assert_eq!(run("1000 30"), [0, 1024]);
        assert_eq!(run("1000 24 1"), [0, 0, 1024]);
        assert_eq!(run("1000 30:5"), [0, 1024]);
        assert_eq!(run("1000 30:27"), [0, 1028]);
        assert_eq!(run("500 2000:100"), [0, 2048]);
        assert_eq!(run("500 2000:100:1900"), [0, 2401]);
        assert_eq!(run("f 100 f 50"), [0, 0, 100, 100]);
    }

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

    #[test]
    fn a_streamed_insertion_is_written_as_the_whole_insertion() {
        // `insert_streamed` writes pieces; every model of NCBI's buffers must write what it writes
        // for the whole insertion as one write, before it and after it. A row of 70,000 bytes
        // ends 99 bytes after a newline: no piece may end there (a terminal would keep only 99).
        let row = [&[b'r'; 69_900][..], b"\n", &[b'r'; 99]].concat();
        let segments: &[(&[u8], u64)] = &[
            (b"\n  Organism Report\nrow\n", 200_000),
            (&row, 3_000),
            (b"|x\ny", 5_000),
            (b"z", 70_000),
            (b"tail\nend", 0),
        ];
        let mut whole = Vec::new();
        for (bytes, spaces) in segments {
            whole.extend_from_slice(bytes);
            whole.resize(whole.len() + *spaces as usize, b' ');
        }
        let streamed = |out: &mut dyn std::io::Write| {
            insert_streamed(out, whole.len() as u64, |insertion| {
                for (bytes, spaces) in segments {
                    insertion.bytes(bytes)?;
                    insertion.spaces(*spaces)?;
                }
                Ok(())
            })
            .unwrap()
        };
        let mut flat = Vec::new();
        streamed(&mut flat);
        assert_eq!(flat, whole);
        for stream in [
            NcbiStdio::File,
            NcbiStdio::Cout {
                block: 4096,
                lines: false,
            },
            NcbiStdio::Cout {
                block: 1024,
                lines: true,
            },
            NcbiStdio::Cout {
                block: 100,
                lines: false,
            },
        ] {
            for before in [0, 1, 1000, 4095, 4096, 8000] {
                for after in [&b"ab"[..], b"ab\ncd", &[b'y'; 9000]] {
                    let lead = vec![b'w'; before];
                    // The whole as one insertion, between the same insertions.
                    let one = written_after(
                        stream,
                        &[Ok(lead.clone()), Ok(whole.clone()), Ok(after.to_vec())],
                    );
                    // The whole streamed.
                    let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
                    let mut buffer = ReportStream::new(
                        Box::new(Writes(writes.clone())),
                        Some(NcbiBuffer::new(stream)),
                        false,
                        "blastn",
                        Site::Report,
                    );
                    let mut pieces = Vec::new();
                    std::io::Write::write_all(&mut buffer, &lead).unwrap();
                    pieces.push(lock(&writes).iter().sum::<usize>());
                    streamed(&mut buffer);
                    pieces.push(lock(&writes).iter().sum::<usize>());
                    std::io::Write::write_all(&mut buffer, after).unwrap();
                    pieces.push(lock(&writes).iter().sum::<usize>());
                    drop(buffer);
                    let total = lock(&writes).iter().sum::<usize>();
                    assert_eq!((pieces, total), (one.0, one.1), "{stream:?} {before}");
                }
            }
        }
        // In pieces: the run of 200,000 spaces ends three, the rest is the last.
        let pieces = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        streamed(&mut Writes(pieces.clone()));
        assert_eq!(*lock(&pieces), [65_559, 65_536, 65_536, 151_405]);
    }

    /// A device that takes `0` more bytes and then fails (a file at its size limit,
    /// `prlimit --fsize`: the write that reaches the limit is cut there, the next fails).
    struct FailAt(usize);

    impl std::io::Write for FailAt {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            if self.0 == 0 {
                return Err(std::io::Error::other("File too large"));
            }
            let n = buf.len().min(self.0);
            self.0 -= n;
            Ok(n)
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    fn native(result: anyhow::Result<()>) -> (i32, String) {
        let error = result.expect_err("the run fails");
        let native = error
            .downcast_ref::<crate::cli::NativeError>()
            .expect("NCBI's end");
        (native.exit, native.message.clone())
    }

    #[test]
    fn a_tabular_report_that_fails_at_byte_n_terminates_as_ncbis() {
        use std::io::Write;
        // Held to the end: the flush of the last `CBlastTabularInfo` fails while it throws.
        let mut stream = ReportStream::file(Box::new(FailAt(10)), "blastn", Site::Tabular);
        stream.write_all(&[b'x'; 100]).unwrap();
        assert_eq!(
            native(stream.finish(Ok(()))),
            (-6, TERMINATE_TEXT.to_string())
        );
        // Where NCBI's buffer writes: the insertion fails, and every later one.
        let mut stream = ReportStream::file(Box::new(FailAt(10)), "blastn", Site::Tabular);
        stream.write_all(&[b'x'; 8000]).unwrap();
        let failed = stream.write_all(&[b'x'; 500]).unwrap_err();
        assert!(is_report_write_failure(&failed));
        assert!(stream.write_all(b"x").is_err());
        assert_eq!(native(stream.finish(Err(failed.into()))).0, -6);
    }

    #[test]
    fn the_alignments_post_802_4_in_a_file_and_on_standard_output_at_their_flush() {
        use std::io::Write;
        let post = "Error: (802.4) [blastp] [Cannot close serializing output stream] \
                    Exception: basic_ios::clear: iostream error\nBLAST failed to write output\n";
        let mut stream = ReportStream::file(Box::new(FailAt(5)), "blastp", Site::Report);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        ostream.write_all(b"prolog\n").unwrap();
        let mut display = DisplaySeqalign::default();
        display.enter(&mut ostream).unwrap();
        ostream.write_all(b"> subject\n").unwrap();
        assert!(display.leave(&mut ostream).is_err());
        drop(ostream);
        assert_eq!(native(stream.finish(Ok(()))), (6, post.to_string()));
        // Standard output: a failed insertion leaves glibc's buffer empty (no post); the flush
        // at the end of the alignments posts, and the post's tie terminates.
        let mut stream = ReportStream::new(
            Box::new(FailAt(5)),
            Some(NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            })),
            true,
            "blastp",
            Site::Alignments,
        );
        assert!(stream.write_all(&[b'x'; 5000]).is_err());
        let message = "BLAST failed to write output\n".to_string();
        assert_eq!(native(stream.finish(Ok(()))), (6, message));
        let mut stream = ReportStream::new(
            Box::new(FailAt(5)),
            Some(NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            })),
            true,
            "blastp",
            Site::Alignments,
        );
        stream.write_all(b"> subject\n").unwrap();
        assert!(lock(&stream.0).flush().is_err());
        assert_eq!(native(stream.finish(Ok(()))).0, -6);
    }

    /// A second stream on the device of `stream` (the test's `NcbiOstream::report` owns it).
    fn stream_of(stream: &mut ReportStream) -> ReportStream {
        ReportStream(stream.0.clone())
    }

    #[test]
    fn a_failure_after_the_formatter_or_of_sam_is_not_reported() {
        use std::io::Write;
        // Held to the end, written after `~CBlastFormat`: unreported.
        let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
        stream.write_all(b"held").unwrap();
        assert!(stream.finish(Ok(())).is_ok());
        assert!(stream.failed());
        // SAM writes after `~CBlastFormat`.
        let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Silent);
        stream.write_all(&[b'x'; 9000]).unwrap();
        stream.write_all(b"more").unwrap();
        assert!(stream.finish(Ok(())).is_ok());
    }

    /// The posts of `stream` go to a buffer, which is returned.
    fn posts_of(stream: &ReportStream) -> std::sync::Arc<std::sync::Mutex<Vec<u8>>> {
        let posts = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        lock(&stream.0).diagnostics = Some(Box::new(Writes2(posts.clone())));
        posts
    }

    /// A device that keeps what it is given.
    struct Writes2(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

    impl std::io::Write for Writes2 {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            lock(&self.0).extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn an_xml_iteration_that_fails_is_posted_once_and_the_run_goes_on_to_the_epilog() {
        use std::io::Write;
        // probe `T/t6/survey` blastn.f5.qm full2: one post per query from the failing one, then
        // the epilog's insertion throws: "BLAST failed to write output", exit 6.
        let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
        let posts = posts_of(&stream);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        let mut cuts = Cuts::default();
        cuts.push(0, Cut::Site(Site::Report));
        cuts.push(100, Cut::Site(Site::XmlIteration));
        cuts.push(8200, Cut::Insertion);
        cuts.push(8300, Cut::Site(Site::Report));
        let report = vec![b'x'; 8300];
        let failed = cuts
            .cursor()
            .write(&mut ostream, &report, 0, report.len())
            .unwrap_err();
        assert!(is_caught_report_write_failure(&failed));
        // The next query's `Iteration` fails at its first insertion.
        let mut next = Cuts::default();
        next.push(0, Cut::Site(Site::XmlIteration));
        next.push(10, Cut::Insertion);
        let failed = next
            .cursor()
            .write(&mut ostream, &report[..20], 0, 20)
            .unwrap_err();
        assert!(is_caught_report_write_failure(&failed));
        let post = "Error: [blastn] Failed s_BlastXMLAddIteration Q(0/1\n";
        assert_eq!(
            String::from_utf8(lock(&posts).clone()).unwrap(),
            post.repeat(2)
        );
        // The epilog, outside the `try`.
        ostream.site(Site::Report).unwrap();
        let failed = ostream.write_all(b"</BlastOutput>").unwrap_err();
        assert!(is_report_write_failure(&failed) && !is_caught_report_write_failure(&failed));
        drop(ostream);
        assert_eq!(
            native(stream.finish(Err(failed.into()))),
            (6, "BLAST failed to write output\n".to_string())
        );
        // On standard output the post's tie throws: `std::terminate`, no post.
        let mut stream = ReportStream::new(
            Box::new(FailAt(0)),
            Some(NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            })),
            true,
            "blastn",
            Site::XmlIteration,
        );
        let posts = posts_of(&stream);
        let failed = stream.write_all(&[b'x'; 5000]).unwrap_err();
        assert!(!is_caught_report_write_failure(&failed));
        assert!(lock(&posts).is_empty());
        assert_eq!(native(stream.finish(Err(failed.into()))).0, -6);
    }

    #[test]
    fn a_report_held_to_its_epilogs_flush_fails_there() {
        // A small BLAST XML report that `basic_filebuf` holds: `endl` at the epilog flushes it
        // while the formatter lives (exit 6); without that flush the bytes would reach the device
        // after `~CBlastFormat`, unreported.
        let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastp", Site::Report);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        let mut cuts = Cuts::default();
        cuts.push(0, Cut::Site(Site::Report));
        cuts.push(5, Cut::Insertion);
        cuts.push(6, Cut::Flush);
        let failed = cuts
            .cursor()
            .write(&mut ostream, b"<x/>\n\n", 0, 6)
            .unwrap_err();
        drop(ostream);
        assert_eq!(native(stream.finish(Err(failed.into()))).0, 6);
    }

    /// One Seq-annot written as `seqalign::write_query_to` writes it: whole, `Close`'s flush, and
    /// NCBI's pieces where the whole insertion fails.
    fn serial_object(
        ostream: &mut NcbiOstream<'_>,
        len: usize,
        pieces: &[usize],
    ) -> std::io::Result<()> {
        use std::io::Write;
        ostream.site(Site::Serial { whole: true })?;
        match ostream
            .write_all(&vec![b'x'; len])
            .and_then(|()| ostream.ncbi_flush())
        {
            Err(error) if is_undecided_report_write_failure(&error) => {
                ostream.resolve_object(pieces)?
            }
            written => written?,
        }
        ostream.site(Site::Report)
    }

    #[test]
    fn a_serial_object_that_fails_inside_ends_the_run_and_one_that_fails_at_close_goes_on() {
        let post = "Error: (802.4) [blastn] [Cannot close serializing output stream] Exception: \
                    basic_ios::clear: iostream error\n";
        let failed = "BLAST failed to write output\n";
        let pieces = [4096, 4096, 1808];
        // `basic_filebuf` holds the first piece and writes 8,192 bytes at the second: past 8,000
        // bytes that write fails inside the object (802.4 at `Close`, then `CATCH_ALL`).
        let mut stream = ReportStream::file(Box::new(FailAt(8000)), "blastn", Site::Report);
        let posts = posts_of(&stream);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        let error = serial_object(&mut ostream, 10_000, &pieces).unwrap_err();
        drop(ostream);
        assert!(lock(&posts).is_empty());
        assert_eq!(
            native(stream.finish(Err(error.into()))),
            (6, format!("{post}{failed}"))
        );
        // Up to 9,000 bytes that write succeeds and only `Close`'s flush fails: posted, the run
        // goes on, and the next object's first insertion throws.
        let mut stream = ReportStream::file(Box::new(FailAt(9000)), "blastn", Site::Report);
        let posts = posts_of(&stream);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        serial_object(&mut ostream, 10_000, &pieces).unwrap();
        assert_eq!(String::from_utf8(lock(&posts).clone()).unwrap(), post);
        let error = serial_object(&mut ostream, 50, &[50]).unwrap_err();
        drop(ostream);
        assert_eq!(
            native(stream.finish(Err(error.into()))),
            (6, format!("{post}{failed}"))
        );
        // A report that the file's buffer holds to `Close`: its flush fails, the run ends well
        // (NCBI 2.17.0, `blastp.f9.q1`: exit 0 after one 802.4).
        let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastp", Site::Report);
        let posts = posts_of(&stream);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        serial_object(&mut ostream, 3000, &[3000]).unwrap();
        drop(ostream);
        assert!(stream.finish(Ok(())).is_ok());
        assert_eq!(lock(&posts).len(), post.len());
        // On standard output the post's tie throws.
        let mut stream = ReportStream::new(
            Box::new(FailAt(0)),
            Some(NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            })),
            true,
            "blastn",
            Site::Report,
        );
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        let error = serial_object(&mut ostream, 3000, &[3000]).unwrap_err();
        drop(ostream);
        assert_eq!(native(stream.finish(Err(error.into()))).0, -6);
    }

    #[test]
    fn an_error_archive_into_the_failed_stream_is_a_write_fault() {
        use std::io::Write;
        // `PrintErrorArchive` after `CATCH_ALL` (probe `T/t6/survey` blastn.f11.q1 full2: 802.5
        // at byte 356, the error archive's length).
        let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
        let posts = posts_of(&stream);
        let failed = stream.write_all(&[b'x'; 9000]).unwrap_err();
        assert_eq!(native(stream.finish(Err(failed.into()))).0, 6);
        stream.site(Site::Serial { whole: true }).unwrap();
        let error = stream.write_all(&[b'a'; 356]).unwrap_err();
        assert!(is_undecided_report_write_failure(&error));
        stream.write_fault(356).unwrap();
        stream.site(Site::Report).unwrap();
        assert_eq!(
            String::from_utf8(lock(&posts).clone()).unwrap(),
            "Error: (802.5) [blastn] CObjectOStream: error at byte 356: stack is empty: write \
             fault\n"
        );
    }

    #[test]
    fn an_xml2_report_that_fails_is_caught_after_its_close() {
        // probe blastn.f16: 802.4 at `Close`, then the catch's post; the run goes on.
        let mut stream = ReportStream::file(Box::new(FailAt(100)), "blastn", Site::Report);
        let posts = posts_of(&stream);
        let mut ostream = NcbiOstream::report(stream_of(&mut stream));
        let report = [
            &b"<BlastOutput2>\n"[..],
            &b"  <x>12345</x>\n".repeat(700),
            b"</BlastOutput2>\n",
        ]
        .concat();
        let mut cuts = Cuts::default();
        cuts.push(0, Cut::Site(Site::Xml2Report { json: false }));
        cuts.push(report.len(), Cut::Flush);
        cuts.push(report.len(), Cut::Site(Site::Report));
        let error = cuts
            .cursor()
            .write(&mut ostream, &report, 0, report.len())
            .unwrap_err();
        assert!(is_caught_report_write_failure(&error));
        assert_eq!(
            String::from_utf8(lock(&posts).clone()).unwrap(),
            "Error: (802.4) [blastn] [Cannot close serializing output stream] Exception: \
             basic_ios::clear: iostream error\nError: [blastn] XML format failed\n"
        );
    }

    #[test]
    fn a_post_while_standard_output_failed_terminates() {
        use std::io::Write;
        let mut stream = ReportStream::new(
            Box::new(FailAt(0)),
            Some(NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            })),
            true,
            "blastn",
            Site::Report,
        );
        stream.write_all(b"held").unwrap();
        // `cerr`'s tie writes what `cout` holds and fails inside the post's destructor.
        assert!(lock(&stream.0).tie().is_err());
        assert_eq!(native(stream.finish(Ok(()))).0, -6);
    }
}

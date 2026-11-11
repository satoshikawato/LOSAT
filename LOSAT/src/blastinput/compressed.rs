//! Compressed FASTA inputs: a file whose name ends in `.gz`, `.bz2` or `.zst` (in any case) is
//! read through the matching decoder, the rule by which NCBI's makeblastdb reads its `-in` files
//! (line GZ, DECISIONS 1, 2, 4, 8 and 9). The decoders are pure Rust (flate2's zlib-rs backend,
//! a port of zlib's inflate, for deflate; bzip2's libbz2-rs-sys backend; ruzstd), as the pure-Rust
//! runtime boundary requires.
//!
//! A `.gz` file is read as makeblastdb's `CDecompressIStream(f, eGZipFile)` reads it
//! (`GzipMembers`, a port of `CZipDecompressor` with the flags `fGZip`): each member's header is
//! looked for in the next 4 KiB of the file and skipped by makeblastdb's own parser (its header
//! CRC is not checked), the deflate data is inflated raw, and the 8-byte trailer is skipped
//! without being checked (a wrong CRC-32 or ISIZE, or a trailer cut short, reads as the intact
//! file). Every member is read, as bgzip and `cat a.gz b.gz` write them; bytes after the last
//! member that are not a gzip header are inflated as deflate data too, and are accepted only
//! where makeblastdb accepts them (fewer than 4 KiB that inflate to nothing). A file that does
//! not start with a gzip header is read as a zlib stream (RFC 1950), up to that stream's end.
//! NCBI's makeblastdb reads only the first bzip2 stream and the first zstd frame and drops the
//! rest without a message; LOSAT reads them all (`PD-LOSAT-NCBI-DEFECTS`, makeblastdb exception).
//!
//! Each format has one decoder that is fed the file's bytes in pieces (`Inflater`, in the shape
//! of NCBI's `CCompressionProcessor::Process` and `Finish`): the reader of the command line
//! (`Compression::decoder`) and LOSAT Web's adapter (`PushInflater`) both run it, so a file
//! inflates to the same bytes, or fails with the same words (`Compression::not_data`), however
//! it is cut into pieces.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050
//! ```c++
//! 		unique_ptr<CDecompressIStream>  comp_seq_istream;
//! 		bool is_gz = false, is_bz2 = false, is_zst=false;
//! 		CCompressStream::EMethod method = CDecompressIStream::eGZipFile; // default gzip method if compressed
//!
//! 		if( (is_gz  = NStr::EndsWith( seq_file, ".gz", NStr::eNocase)) ||
//! 		    (is_bz2 = NStr::EndsWith( seq_file, ".bz2",NStr::eNocase)) ||
//! 		    (is_zst = NStr::EndsWith( seq_file, ".zst",NStr::eNocase))
//! 		  )
//! 		{
//! 		    if( is_bz2  ) method = CDecompressIStream::eBZip2;
//! 		    if( is_zst  ) method = CDecompressIStream::eZstd;
//! 		    comp_seq_istream.reset(new CDecompressIStream(f , method ));
//! 		    data_in =  comp_seq_istream.get();
//! 		}
//! ```

use std::io::{self, BufRead, Read};

/// The compression a file name announces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Compression {
    Gzip,
    Bzip2,
    Zstd,
}

impl Compression {
    pub const ALL: [Compression; 3] = [Compression::Gzip, Compression::Bzip2, Compression::Zstd];

    /// The suffix of the name, as makeblastdb tests it (`.gz`, `.bz2`, `.zst`).
    pub fn suffix(self) -> &'static str {
        match self {
            Compression::Gzip => ".gz",
            Compression::Bzip2 => ".bz2",
            Compression::Zstd => ".zst",
        }
    }

    /// The name of the format in LOSAT's messages.
    pub fn format(self) -> &'static str {
        match self {
            Compression::Gzip => "gzip",
            Compression::Bzip2 => "bzip2",
            Compression::Zstd => "zstd",
        }
    }

    /// A reader of the inflated bytes of `input`: every member, stream or frame, one after the
    /// other. A fault of the compressed data (not this format, cut short, an empty file) is an
    /// error of the read that meets it; for gzip, that is what makeblastdb fails on
    /// (`GzipMembers`), and a file that is neither gzip nor zlib data fails at the first read.
    /// The reader runs this format's `Inflater` on blocks of `input` (`PullInflater`,
    /// `read_block`).
    pub fn decoder<'a>(self, input: impl Read + 'a) -> Box<dyn Read + 'a> {
        Box::new(PullInflater {
            source: io::BufReader::with_capacity(self.read_block(), input),
            inflater: self.inflater(),
            ended: false,
        })
    }

    /// The size of the blocks in which the command line's reader reads a file of this format.
    /// 64 KiB for gzip and bzip2. 256 KiB for zstd: a block of the largest size (128 KiB, with
    /// its 3-byte header and the frame's 4-byte checksum) then mostly lies inside one read
    /// block, and ruzstd decodes it there instead of from a copy (`ZstdFrames::gather`); the
    /// output of each read does not hang on the size, which only says when a whole block is at
    /// hand. A gzip read does hang on it where the data is damaged (an `inflate()` call that
    /// fails gives none of its output, so the size decides how many inflated bytes reach the
    /// reader before the error), and bzip2 reads through `BZIP2_WINDOW` windows, which are the
    /// same at any multiple of 8 KiB; both keep the size they had.
    fn read_block(self) -> usize {
        match self {
            Compression::Gzip | Compression::Bzip2 => 1 << 16,
            Compression::Zstd => 1 << 18,
        }
    }

    /// The inflater of a file of this format, to be fed its bytes in pieces (`Inflater::step`).
    pub fn inflater(self) -> Inflater {
        let core = match self {
            Compression::Gzip => Core::Gzip(GzipMembers::new()),
            Compression::Bzip2 => Core::Bzip2(Bzip2Streams::new()),
            Compression::Zstd => Core::Zstd(Box::new(ZstdFrames::new())),
        };
        Inflater {
            core,
            started: false,
            failed: None,
        }
    }

    /// LOSAT's sentence for the file `name` (with this format's suffix) whose bytes are not this
    /// format, or are cut short: `fault` is the error of the inflater (line GZ, DECISIONS 5;
    /// for `.gz`: where makeblastdb's decompressor fails, `GzipMembers`). The command line puts
    /// the option before it (`-query: `, `-subject: `, `-in: `).
    pub fn not_data(self, name: &str, fault: &str) -> String {
        format!(
            "'{name}' is named {suffix} but is not {format} data ({fault}); LOSAT reads a file named {suffix} through {format}",
            suffix = self.suffix(),
            format = self.format(),
        )
    }
}

/// The compression a file name announces by its suffix, in any case (`NStr::eNocase`), or `None`
/// for a plain name. `-` and names without a suffix are plain.
pub fn compression_of_name(path: &std::path::Path) -> Option<Compression> {
    let bytes = path.as_os_str().as_encoded_bytes();
    Compression::ALL.into_iter().find(|compression| {
        let suffix = compression.suffix().as_bytes();
        bytes.len() >= suffix.len()
            && bytes[bytes.len() - suffix.len()..].eq_ignore_ascii_case(suffix)
    })
}

/// What one `Inflater::step` did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Progress {
    /// The bytes taken from the start of the input. The caller gives the rest of the input again,
    /// before what follows it.
    pub taken: usize,
    /// The inflated bytes written at the start of the output.
    pub written: usize,
    /// The data has ended (`eStatus_EndOfData`): no more output follows, and the input after
    /// the end is not read (a zlib stream's end, or the end of the file).
    pub ended: bool,
}

impl Progress {
    /// The input is all taken and more is needed before the next output.
    fn more(taken: usize) -> Self {
        Self {
            taken,
            written: 0,
            ended: false,
        }
    }

    fn written(taken: usize, written: usize) -> Self {
        Self {
            taken,
            written,
            ended: false,
        }
    }

    fn ended(taken: usize) -> Self {
        Self {
            taken,
            written: 0,
            ended: true,
        }
    }
}

/// The inflater of one compressed file, fed its bytes in pieces: the one decoder of each format,
/// which the reader of the command line (`Compression::decoder`) and the pieces of LOSAT Web
/// (`PushInflater`) both run, so that a file gives the same bytes, or fails with the same words,
/// however it is cut. Built by `Compression::inflater`.
///
/// NCBI reference (598d8ae6): c++/include/util/compress/compress.hpp:510-521, 532-540
/// ```c++
///     /// Compress/decompress as much data as possible, and stops when the input
///     /// buffer becomes empty or the output buffer becomes full. It may
///     /// introduce some output latency (reading input without producing any
///     /// output).
///     virtual EStatus Process
///     (const char* in_buf,      // [in]  input buffer
///      size_t      in_len,      // [in]  input data length
///      char*       out_buf,     // [in]  output buffer
///      size_t      out_size,    // [in]  output buffer size
///      size_t*     in_avail,    // [out] count unproc.bytes in input buffer
///      size_t*     out_avail    // [out] count bytes putted into out buffer
///      ) = 0;
/// ...
///     /// Finish the compression/decompression process.
///     /// Process pending input, flush pending output.
///     /// This function slightly like to Flush(), but it must be called only
///     /// at the end of compression process just before End().
///     virtual EStatus Finish
///     (char*       out_buf,     // [in]  output buffer
///      size_t      out_size,    // [in]  output buffer size
///      size_t*     out_avail    // [out] count bytes putted into out buffer
///      ) = 0;
/// ```
/// `step` is `Process` (`at_end` false) and `Finish` (`at_end` true, after the last input) in
/// one: `Progress::taken` is `in_len - *in_avail`, `Progress::written` is `*out_avail`,
/// `Progress::ended` is `eStatus_EndOfData`, and an error is `eStatus_Error`.
pub struct Inflater {
    core: Core,
    /// A byte of input was given: the file is not empty.
    started: bool,
    /// The error that ended the inflation, which every later step gives again.
    failed: Option<(io::ErrorKind, String)>,
}

/// The error of a step that, given input or at the end of the input, took none of it, wrote
/// nothing and did not end: a decoder that would make a reader or a push loop wait forever. No
/// decoder does so on any input known; the check keeps a fault of one from hanging the command
/// line or a browser tab.
const STALLED: &str = "the decoder neither took input nor gave output";

/// The decoder of each format, on pushed input.
enum Core {
    Gzip(GzipMembers),
    Bzip2(Bzip2Streams),
    Zstd(Box<ZstdFrames>),
}

impl Inflater {
    /// Inflates from `input` (any length, cut anywhere) into `out`. `at_end` says that no input
    /// follows `input`. A step writes no output and does not end only where it needs more input
    /// than it was given: it has then taken all of `input` (a unit that must be whole before it
    /// is decoded, such as a zstd block or the 4 KiB of a gzip header check, is kept inside),
    /// and it is never so at the end. A step that, given input or at the end, neither takes,
    /// writes nor ends is an error (`STALLED`), so that a loop of steps always ends, whatever a
    /// decoder does. An empty `out` takes and gives nothing (NCBI's `eStatus_Overflow`). After
    /// an error, every step gives the same error.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1358-1361
    /// ```c++
    ///     *out_avail = 0;
    ///     if ( !out_size ) {
    ///         return eStatus_Overflow;
    ///     }
    /// ```
    /// No NCBI counterpart for the empty file: makeblastdb's `x_VerifyInputFilesType` refuses
    /// an empty `-in` before it is opened; at the end of an input that brought no byte, the
    /// inflater fails with "the file is empty" for every format (no gzip or bzip2 header, no
    /// zstd frame), where the decoders would otherwise read it as an empty input (line GZ,
    /// DECISIONS 5).
    pub fn step(&mut self, input: &[u8], at_end: bool, out: &mut [u8]) -> io::Result<Progress> {
        if let Some((kind, message)) = &self.failed {
            return Err(io::Error::new(*kind, message.clone()));
        }
        if out.is_empty() {
            return Ok(Progress::default());
        }
        self.started |= !input.is_empty();
        let result = if at_end && !self.started {
            Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "the file is empty",
            ))
        } else {
            match &mut self.core {
                Core::Gzip(core) => core.step(input, at_end, out),
                Core::Bzip2(core) => core.step(input, at_end, out),
                Core::Zstd(core) => core.step(input, at_end, out),
            }
        };
        let result = result.and_then(|progress| {
            if progress == Progress::default() && (at_end || !input.is_empty()) {
                return Err(io::Error::new(io::ErrorKind::InvalidData, STALLED));
            }
            Ok(progress)
        });
        if let Err(error) = &result {
            self.failed = Some((error.kind(), error.to_string()));
        }
        result
    }
}

/// The reader of the command line (`Compression::decoder`): the inflater on the blocks of a
/// buffered source, as `CCompressionStreambuf` reads makeblastdb's `CDecompressIStream`.
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/streambuf.cpp:377-441
/// ```c++
///     // Put data into the (de)compressor until there is something
///     // in the output buffer
///     do {
///         ...
///             // Refill the input buffer if necessary
///             if ( m_Reader->m_Begin == m_Reader->m_End ) {
///                 n_read = m_Stream->rdbuf()->sgetn(m_Reader->m_InBuf, m_Reader->m_InBufSize);
///                 ...
///                 if ( !n_read ) {
///                     // We can't read more of data.
///                     // Automatically 'finalize' (de)compressor.
///                     m_Reader->m_State = CSP::eFinalize;
///                     return Flush(CCompressionStream::eRead) == 0;
///                 }
///                 ...
///             }
///             // Process next data portion
///             in_len = m_Reader->m_End - m_Reader->m_Begin;
///             m_Reader->m_LastStatus = m_Reader->m_Processor->Process(
///                                         m_Reader->m_Begin, in_len, egptr(), out_size,
///                                         &in_available, &out_available);
///         ...
///         if ( m_Reader->m_LastStatus == CP::eStatus_Error ) {
///             THROW_FAILURE;
///             // return false;
///         }
///         ...
///         // Update pointer to an unprocessed data
///         m_Reader->m_Begin += (in_len - in_available);
///         ...
///         if ( m_Reader->m_LastStatus == CP::eStatus_EndOfData   &&  !out_available ) {
///             return false;
///         }
///
///     } while ( !out_available );
/// ```
/// Each read gives the output of one step that writes; the end of the source is the step with
/// `at_end` (`Flush` with the state `eFinalize` calls `Finish`).
struct PullInflater<R: Read> {
    source: io::BufReader<R>,
    inflater: Inflater,
    ended: bool,
}

impl<R: Read> Read for PullInflater<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() || self.ended {
            return Ok(0);
        }
        // Each step takes input, writes, ends or fails (`Inflater::step`, `STALLED`), and an
        // empty block from the source is its end: the loop ends.
        loop {
            let input = self.source.fill_buf()?;
            let at_end = input.is_empty();
            let progress = self.inflater.step(input, at_end, buf)?;
            self.source.consume(progress.taken);
            self.ended = progress.ended;
            if progress.written > 0 || progress.ended {
                return Ok(progress.written);
            }
        }
    }
}

/// The inflater of LOSAT Web's adapter: the pieces of a compressed file go in as they are read
/// (`push`), the inflated bytes come out through a sink in pieces of up to `OUTPUT_SIZE` bytes,
/// and `finish` ends the file. It runs the same `Inflater` as the command line's reader: a file
/// that the reader inflates gives the reader's bytes, and one that the reader refuses fails
/// with the reader's error (kind and words). After an error the value gives that error again.
/// The bytes the sink got before an error are not the reader's bytes before it: how far a
/// decoder gets before it meets the fault hangs on the sizes of the pieces (a failing
/// `inflate()` call gives none of its output, for one), so the caller drops them; LOSAT Web's
/// adapter refuses the whole file.
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/streambuf.cpp:469-496
/// ```c++
///     // Loop until no data is left
///     while ( in_available ) {
///         // Process next data portion
///         size_t out_available = 0;
///         size_t out_size = m_Writer->m_OutBuf +
///                           m_Writer->m_OutBufSize - m_Writer->m_End;
///         m_Writer->m_LastStatus = m_Writer->m_Processor->Process(
///             in_buf + count - in_available, in_available, m_Writer->m_End, out_size,
///             &in_available, &out_available);
///
///         // Check on error / small output buffer
///         if ( m_Writer->m_LastStatus == CP::eStatus_Error ) {
///             THROW_FAILURE;
///             //return false;
///         }
///         // No more data -- automatically finalize stream
///         if ( m_Writer->m_LastStatus == CP::eStatus_EndOfData ) {
///             m_Writer->m_State = CSP::eFinalize;
///         }
///         // Update the output buffer pointer
///         m_Writer->m_End += out_available;
///
///         // Write data to the underlying stream only if the output buffer
///         // is full or an overflow occurs.
///         if ( !WriteOutBufToStream() ) {
///             return false;
///         }
///     }
/// ```
/// (`CCompressionStreambuf::ProcessStreamWrite`, the writer of a `CDecompressOStream`.) The
/// output buffer goes to the sink when it is full and at the end of each `push`; after the end
/// of the data the rest of the input is not read.
pub struct PushInflater {
    inflater: Inflater,
    out: Box<[u8]>,
    filled: usize,
    ended: bool,
}

impl PushInflater {
    /// The size of the output buffer: the largest piece given to a sink.
    pub const OUTPUT_SIZE: usize = 1 << 20;

    pub fn new(compression: Compression) -> Self {
        Self {
            inflater: compression.inflater(),
            out: vec![0; Self::OUTPUT_SIZE].into_boxed_slice(),
            filled: 0,
            ended: false,
        }
    }

    /// Inflates all of `chunk`, the next bytes of the file, and gives `sink` every inflated
    /// byte the inflater has released. On an error, what the sink got from this file is to be
    /// dropped (see `PushInflater`).
    pub fn push(&mut self, mut chunk: &[u8], sink: &mut dyn FnMut(&[u8])) -> io::Result<()> {
        while !chunk.is_empty() && !self.ended {
            let progress = self
                .inflater
                .step(chunk, false, &mut self.out[self.filled..])?;
            // A step given input takes, writes, ends or fails (`STALLED`): the loop ends.
            chunk = &chunk[progress.taken..];
            self.keep(progress, sink);
        }
        self.flush(sink);
        Ok(())
    }

    /// The end of the file: gives `sink` the rest of the inflated bytes, or the error of data
    /// cut short (or of an empty file); on an error, what the sink got from this file is to be
    /// dropped (see `PushInflater`).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/streambuf.cpp:238-240, 277-279
    /// ```c++
    ///     // Finish. Change state to 'finalized'.
    ///     sp->m_State = CSP::eFinalize;
    ///     return Flush(dir);
    /// ...
    ///         if ( sp->m_State == CSP::eFinalize ) {
    ///             // State is eFinalize
    ///             sp->m_LastStatus = sp->m_Processor->Finish(buf, out_size, &out_available);
    /// ```
    pub fn finish(mut self, sink: &mut dyn FnMut(&[u8])) -> io::Result<()> {
        while !self.ended {
            let progress = self
                .inflater
                .step(&[], true, &mut self.out[self.filled..])?;
            // At the end a step writes, ends or fails (`STALLED`): the loop ends.
            self.keep(progress, sink);
        }
        self.flush(sink);
        Ok(())
    }

    /// Counts a step's output, and gives a full buffer to the sink.
    fn keep(&mut self, progress: Progress, sink: &mut dyn FnMut(&[u8])) {
        self.filled += progress.written;
        self.ended = progress.ended;
        if self.filled == self.out.len() {
            self.flush(sink);
        }
    }

    fn flush(&mut self, sink: &mut dyn FnMut(&[u8])) {
        if self.filled > 0 {
            sink(&self.out[..self.filled]);
            self.filled = 0;
        }
    }
}

/// The largest gzip header that makeblastdb's decompressor looks for: the header of each member
/// is looked for in the next `MAX_HEADER_SIZE` bytes of the file.
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:99-100
/// ```c++
/// // Maximum size of gzip file header
/// const size_t kMaxHeaderSize = 1024*4;
/// ```
const MAX_HEADER_SIZE: usize = 1024 * 4;

/// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:137-146
/// ```c++
/// // gzip magic header
/// const unsigned char gz_magic[2] = {0x1f, 0x8b};
///
/// // gzip flag byte
/// #define ASCII_FLAG   0x01 // bit 0 set: file probably ascii text
/// #define HEAD_CRC     0x02 // bit 1 set: header CRC present
/// #define EXTRA_FIELD  0x04 // bit 2 set: extra field present
/// #define ORIG_NAME    0x08 // bit 3 set: original file name present
/// #define COMMENT      0x10 // bit 4 set: file comment present
/// #define RESERVED     0xE0 // bits 5..7: reserved
/// ```
const GZ_MAGIC: [u8; 2] = [0x1f, 0x8b];
const HEAD_CRC: u8 = 0x02;
const EXTRA_FIELD: u8 = 0x04;
const ORIG_NAME: u8 = 0x08;
const COMMENT: u8 = 0x10;
const RESERVED: u8 = 0xE0;
/// zlib.h `Z_DEFLATED`.
const Z_DEFLATED: u8 = 8;

/// The length of the gzip header at the start of `buf`, or 0 where makeblastdb finds none: fewer
/// than 10 bytes, another magic or method, a reserved flag set, or a header that does not end
/// within `buf`. The extra field, the name and the comment are skipped, and so are the 2 bytes of
/// the header CRC, which is not checked. A name or comment without its NUL ends at the end of
/// `buf` (so a name longer than `buf` makes a header of all of `buf`).
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:151-211
/// ```c++
/// static
/// size_t s_CheckGZipHeader(const void* src_buf, size_t src_len,
///                          CZipCompression::SFileInfo* info = 0)
/// {
///     unsigned char* buf = (unsigned char*)src_buf;
///     // .gz header cannot be less than 10 bytes
///     if (src_len < 10) {
///         return 0;
///     }
///     // Check the gzip magic header
///     if (buf[0] != gz_magic[0]  ||
///         buf[1] != gz_magic[1]) {
///         return 0;
///     }
///     int method = buf[2];
///     int flags  = buf[3];
///     if (method != Z_DEFLATED  ||  (flags & RESERVED) != 0) {
///         return 0;
///     }
///     // Header length:
///     // gz_magic (2) + method (1) + flags (1) + time, xflags and OS code (6)
///     size_t header_len = 10;
///     ...
///     // Skip the extra fields
///     if ((flags & EXTRA_FIELD) != 0) {
///         if (header_len + 2 > src_len) {
///             return 0;
///         }
///         size_t len = buf[header_len++];
///         len += ((size_t)buf[header_len++])<<8;
///         header_len += len;
///     }
///     // Skip the original file name
///     if ((flags & ORIG_NAME) != 0) {
///         size_t pos = header_len;
///         while (header_len < src_len  &&  buf[header_len++] != 0);
///         ...
///     }
///     // Skip the file comment
///     if ((flags & COMMENT) != 0) {
///         size_t pos = header_len;
///         while (header_len < src_len  &&  buf[header_len++] != 0);
///         ...
///     }
///     // Skip the header CRC
///     if ((flags & HEAD_CRC) != 0) {
///         header_len += 2;
///     }
///     if (header_len > src_len) {
///         return 0;
///     }
///     return header_len;
/// }
/// ```
fn check_gzip_header(buf: &[u8]) -> usize {
    let src_len = buf.len();
    if src_len < 10 {
        return 0;
    }
    if buf[0] != GZ_MAGIC[0] || buf[1] != GZ_MAGIC[1] {
        return 0;
    }
    let method = buf[2];
    let flags = buf[3];
    if method != Z_DEFLATED || (flags & RESERVED) != 0 {
        return 0;
    }
    let mut header_len = 10;
    if (flags & EXTRA_FIELD) != 0 {
        if header_len + 2 > src_len {
            return 0;
        }
        let len = usize::from(buf[header_len]) | (usize::from(buf[header_len + 1]) << 8);
        header_len += 2 + len;
    }
    // `while (header_len < src_len  &&  buf[header_len++] != 0);`: past the NUL, or to the end.
    let skip_string = |mut header_len: usize| {
        while header_len < src_len {
            header_len += 1;
            if buf[header_len - 1] == 0 {
                break;
            }
        }
        header_len
    };
    if (flags & ORIG_NAME) != 0 {
        header_len = skip_string(header_len);
    }
    if (flags & COMMENT) != 0 {
        header_len = skip_string(header_len);
    }
    if (flags & HEAD_CRC) != 0 {
        header_len += 2;
    }
    if header_len > src_len {
        return 0;
    }
    header_len
}

/// Where `GzipMembers` is in the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GzipState {
    /// At the start of the file or after a member's trailer: a header is looked for.
    CheckHeader,
    /// Inflating the data after a header check.
    Inflate,
    /// A member ended: so many bytes of its trailer are still to skip.
    Trailer(usize),
    /// No more data: the end of the file, or the end of a zlib stream.
    End,
}

/// The inflated bytes of a `.gz` file, read as makeblastdb reads them: a port of NCBI's
/// `CZipDecompressor` with the flags of `CDecompressIStream::eGZipFile` (`fGZip` =
/// `fCheckFileHeader | fWriteGZipFormat | fAllowConcatenatedGZip`, of which `fWriteGZipFormat`
/// concerns only the compressor; not `fAllowTransparentRead` nor `fAllowEmptyData`), as
/// `CCompressionStreambuf` drives it. Where makeblastdb's stream gives bytes, this gives the same
/// bytes; where it fails, the read fails, with LOSAT's message (line GZ, DECISIONS 5, narrowed for
/// `.gz` by swebok-improve SW8 DECISIONS 38 and 41: a damaged or cut-short `.gz` is refused where
/// makeblastdb's decompressor fails).
///
/// The rules, at each header check (the start of the file and the end of each member's data):
/// - the header is looked for in the next `MAX_HEADER_SIZE` bytes (`check_gzip_header`); where
///   one is found it is skipped and the data is inflated as raw deflate, and where none is found
///   the data is inflated as raw deflate still if an earlier member had a header (`m_IsGZ` is
///   never cleared), or as a zlib stream (RFC 1950, with its Adler-32 checked) at the start of
///   the file;
/// - at the end of a member's deflate data, its 8-byte trailer is skipped, not checked (as many
///   of its bytes as the file still has), and the next header is looked for;
/// - at the end of a zlib stream the reading ends, whatever follows;
/// - a deflate error is an error; at the end of the file inside the data, the read fails, except
///   where the rest of the file after the header check was shorter than `MAX_HEADER_SIZE` (it
///   was then kept for the check, and makeblastdb inflates it in `Finish`) and inflated to no
///   bytes at all: the reading then ends there without an error.
///
/// The deflate and zlib data is inflated by flate2's zlib-rs backend, a Rust port of zlib's
/// `inflate.c` and `inftrees.c`, where makeblastdb calls the system's zlib `inflate()` (the
/// official Linux build links `libz.so.1`) with `Z_SYNC_FLUSH`, window bits 15 (raw: -15), and
/// a new stream for each member (`inflateEnd` + `inflateInit2_`, so no match reaches into an
/// earlier member): the same checks (a distance past the start of the window, a missing
/// end-of-block code, incomplete or over-subscribed code sets), the same input consumed when the
/// data runs out. flate2's default backend, miniz_oxide, reads some data that zlib refuses and
/// refuses some that zlib waits on. A host whose `libz.so.1` is another inflater (zlib-ng in its
/// zlib-compatible build) could in principle differ from zlib where the data is damaged.
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/stream_util.cpp:47, 166-178
/// ```c++
///     const ICompression::TFlags kDefault_GZipFile            = CZipCompression::fGZip;
/// ...
///     case CCompressStream::eGZipFile:
///     case CCompressStream::eConcatenatedGZipFile:
/// #if defined(HAVE_LIBZ)
///         if (flags == CCompressStream::fDefault) {
///             flags = kDefault_GZipFile;
///         } else {
///             flags |= kDefault_GZipFile;
///         }
///         if (type == eCompress) {
///             processor = new CZipStreamCompressor(level, flags);
///         } else {
///             processor = new CZipStreamDecompressor(flags);
///         }
/// ```
/// c++/include/util/compress/zlib.hpp:103-104, 118-120
/// ```c++
///         /// Check (and skip) gzip file header on decompression stage
///         fCheckFileHeader       = (1<<2),
/// ...
///         fAllowConcatenatedGZip = (1<<4),
///         /// Set of flags for gzip file support. See each flag description above.
///         fGZip = fCheckFileHeader | fWriteGZipFormat | fAllowConcatenatedGZip,
/// ```
///
/// The input is pushed (`step`), as `CZipDecompressor::Process` takes it: the bytes of a header
/// check are kept in the cache until it holds `MAX_HEADER_SIZE` bytes or the input ends, and
/// the bytes of a trailer that a later input holds are skipped there (`m_SkipInput`).
struct GzipMembers {
    /// The bytes taken for the last header check (`m_Cache`), inflated from `cache_at` before
    /// more input is taken.
    cache: Vec<u8>,
    cache_at: usize,
    inflater: flate2::Decompress,
    state: GzipState,
    /// A gzip header was found (`m_IsGZ`).
    is_gz: bool,
    /// The rest of the file was shorter than `MAX_HEADER_SIZE` at the last header check.
    rest_in_cache: bool,
}

impl GzipMembers {
    fn new() -> Self {
        Self {
            cache: Vec::with_capacity(MAX_HEADER_SIZE),
            cache_at: 0,
            inflater: flate2::Decompress::new(false),
            state: GzipState::CheckHeader,
            is_gz: false,
            rest_in_cache: false,
        }
    }

    /// Looks for a header in the next `MAX_HEADER_SIZE` bytes and starts the inflater.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1412-1469
    /// ```c++
    ///             if ( m_NeedCheckHeader ) {
    ///                 if (!x_in_buf  &&  !m_Cache.size()) {
    ///                     // Possible Flush(), but we should refill the cache
    ///                     // to perform header check -- so, just ignore.
    ///                     *in_avail  = 0;
    ///                     *out_avail = 0;
    ///                     return eStatus_Success;
    ///                 }
    ///                 if (x_in_buf  &&  m_Cache.size() < kMaxHeaderSize) {
    ///                     size_t n = min(kMaxHeaderSize - m_Cache.size(), x_in_len);
    ///                     m_Cache.append(x_in_buf, n);
    ///                     x_in_buf += n;
    ///                     x_in_len -= n;
    ///                     if (m_Cache.size() < kMaxHeaderSize) {
    ///                         // Data block is very small and was fully cached.
    ///                         *in_avail  = 0;
    ///                         *out_avail = 0;
    ///                         return eStatus_Success;
    ///                     }
    ///                 }
    ///                 // Check gzip header in the buffer
    ///                 header_len = s_CheckGZipHeader(m_Cache.data(), m_Cache.size());
    ///                 _ASSERT(header_len < kMaxHeaderSize);
    ///
    ///                 // If gzip header found, skip it
    ///                 if ( header_len ) {
    ///                     m_Cache.erase(0, header_len);
    ///                     IncreaseProcessedSize(header_len);
    ///                     m_DecompressMode = eMode_Decompress;
    ///                     m_IsGZ = true;
    ///                 }
    ///                 // Reinit decompression stream
    ///                 Z(inflateEnd)(STREAM);
    ///                 int errcode = Z(inflateInit2_)(
    ///                                 STREAM, m_IsGZ ? -GetWindowBits() : GetWindowBits(),
    ///                                 ZLIB_VERSION, (int)sizeof(z_stream));
    ///                 ...
    ///                 // Already skipped, or we don't have header here
    ///                 m_NeedCheckHeader = false;
    ///             }
    /// ```
    /// The cache is filled from the input until it holds `kMaxHeaderSize` bytes; with fewer
    /// bytes left in the file, the input ends first and the check is made in `Finish`, on the rest
    /// of the file (`rest_in_cache`). The window bits are `MAX_WBITS` (15), flate2's. `reset`
    /// starts a new stream with an empty window, as `inflateEnd` + `inflateInit2_` do: a match of
    /// a later member cannot reach back into the member before. False while the cache waits for
    /// more input.
    fn check_header(&mut self, input: &[u8], taken: &mut usize, at_end: bool) -> bool {
        self.cache.drain(..self.cache_at);
        self.cache_at = 0;
        let count = (MAX_HEADER_SIZE - self.cache.len()).min(input.len() - *taken);
        self.cache.extend_from_slice(&input[*taken..*taken + count]);
        *taken += count;
        if self.cache.len() < MAX_HEADER_SIZE && !at_end {
            return false;
        }
        if self.cache.is_empty() {
            // The end of the file after a member's trailer (an empty file is refused by
            // `Inflater::step` before).
            self.state = GzipState::End;
            return true;
        }
        self.rest_in_cache = self.cache.len() < MAX_HEADER_SIZE;
        let header_len = check_gzip_header(&self.cache);
        if header_len > 0 {
            self.cache_at = header_len;
            self.is_gz = true;
        }
        self.inflater.reset(!self.is_gz);
        self.state = GzipState::Inflate;
        true
    }

    /// The deflate data ended (`Z_STREAM_END`): a member's trailer follows, then another header
    /// check; a zlib stream ends the reading.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1513-1526, 1555-1564
    /// ```c++
    ///             // Concatenated file? Try to process next .gz chunk, if present
    ///             if ((errcode == Z_STREAM_END)  &&  m_IsGZ) {
    ///                 // Skip .gz file footer (8 bytes)
    ///                 if (STREAM->avail_in < 8) {
    ///                     m_SkipInput = 8 - STREAM->avail_in;
    ///                     STREAM->avail_in = 0;
    ///                 } else {
    ///                     STREAM->avail_in -= 8;
    ///                 }
    ///                 if ( F_ISSET(fAllowConcatenatedGZip) ) {
    ///                     m_NeedCheckHeader = true;
    ///                     errcode = Z_OK;
    ///                 }
    ///             }
    /// ...
    ///             switch (errcode) {
    ///             case Z_OK:
    ///                 ...
    ///                 return eStatus_Success;
    ///             case Z_STREAM_END:
    ///                 return eStatus_EndOfData;
    ///             }
    /// ```
    /// The trailer's CRC-32 and ISIZE are skipped, never compared; `m_SkipInput` skips the bytes
    /// of the trailer that a later input block holds, and at the end of the file (a trailer
    /// cut short) `Finish` ends the data (`Process` returns `eStatus_Success` with no output).
    fn stream_end(&mut self) {
        self.state = if self.is_gz {
            GzipState::Trailer(8)
        } else {
            GzipState::End
        };
    }

    /// Skips up to `left` bytes of a member's trailer, from the cache, then from the input; false
    /// while more input is to be skipped.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1382-1404
    /// ```c++
    ///         if ( m_SkipInput ) {
    ///             // Skip from cache if present
    ///             if ( m_Cache.size() ) {
    ///                 size_t n = min(m_Cache.size(), m_SkipInput);
    ///                 m_Cache.erase(0, n);
    ///                 m_SkipInput -= n;
    ///                 IncreaseProcessedSize(n);
    ///             }
    ///             // And/or from input stream also
    ///             if ( m_SkipInput ) {
    ///                 size_t n = min(x_in_len, m_SkipInput);
    ///                 x_in_buf += n;
    ///                 x_in_len -= n;
    ///                 m_SkipInput -= n;
    ///                 IncreaseProcessedSize(n);
    ///                 if ( m_SkipInput ) {
    ///                     // Data block is very small... and was skipped.
    ///                     *in_avail  = x_in_len;
    ///                     *out_avail = 0;
    ///                     return eStatus_Success;
    ///                 }
    ///             }
    ///         }
    /// ```
    /// At the end of the input the trailer is cut short, and the next header check finds the end
    /// of the file.
    fn skip_trailer(
        &mut self,
        mut left: usize,
        input: &[u8],
        taken: &mut usize,
        at_end: bool,
    ) -> bool {
        let cached = (self.cache.len() - self.cache_at).min(left);
        self.cache_at += cached;
        left -= cached;
        let count = left.min(input.len() - *taken);
        *taken += count;
        left -= count;
        if left > 0 && !at_end {
            self.state = GzipState::Trailer(left);
            return false;
        }
        self.state = GzipState::CheckHeader;
        true
    }

    /// The file ended inside deflate or zlib data, and the inflater has given out all it can.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1617-1630
    /// ```c++
    ///     // Process only if we have some data -- otherwise zlib return error
    ///     CCompressionProcessor::EStatus status = eStatus_Success;
    ///     if ( GetProcessedSize() || m_Cache.size() ) {
    ///         status = Process(0, 0, out_buf, out_size, &in_avail, out_avail);
    ///     } else {
    ///         ...
    ///     }
    ///     if (status == eStatus_Success  &&  *out_avail == 0) {
    ///         // Possible case with .gz files only (check on concatenated files)
    ///         return eStatus_EndOfData;
    ///     }
    ///     return status;
    /// ```
    /// c++/src/util/compress/api/zlib.cpp:1555-1566
    /// ```c++
    ///             switch (errcode) {
    ///             case Z_OK:
    ///                 ...
    ///                 return eStatus_Success;
    ///             case Z_STREAM_END:
    ///                 return eStatus_EndOfData;
    ///             }
    ///             ERR_COMPRESS(70, FormatErrorMessage("CZipDecompressor::Process", GetProcessedSize()));
    ///             return eStatus_Error;
    /// ```
    /// At the end of the input the stream buffer calls `Finish`, which calls `Process` without
    /// input: `inflate()` makes no progress (`Z_BUF_ERROR`), an error. One `Process` in `Finish`
    /// ends the data instead: the one that inflates a cache holding the rest of the file and
    /// gives no output (`Z_OK`, input consumed). zlib's `inflate()` stops early only for want of
    /// output room, after decoding the next symbol; a later call then gives that symbol's bytes,
    /// so a call without output comes only first after the header check, when nothing of the
    /// rest inflates to a byte.
    fn end_of_input(&mut self) -> io::Result<()> {
        if self.rest_in_cache && self.inflater.total_in() > 0 && self.inflater.total_out() == 0 {
            self.state = GzipState::End;
            return Ok(());
        }
        Err(self.fault("incomplete"))
    }

    /// LOSAT's words for a fault that makes makeblastdb's stream fail (line GZ, DECISIONS 5).
    fn fault(&self, what: &str) -> io::Error {
        let message = if self.is_gz {
            format!("{what} deflate stream")
        } else if self.inflater.total_out() == 0 && what == "corrupt" {
            "invalid gzip header".to_string()
        } else {
            format!("{what} zlib stream")
        };
        let kind = if what == "incomplete" {
            io::ErrorKind::UnexpectedEof
        } else {
            io::ErrorKind::InvalidData
        };
        io::Error::new(kind, message)
    }

    /// One `Inflater::step` of a `.gz` file: the header checks, trailers and inflation up to
    /// the first output, the end of the data, or the end of `input` (all of it taken).
    fn step(&mut self, input: &[u8], at_end: bool, out: &mut [u8]) -> io::Result<Progress> {
        let mut taken = 0;
        loop {
            match self.state {
                GzipState::End => return Ok(Progress::ended(taken)),
                GzipState::CheckHeader => {
                    if !self.check_header(input, &mut taken, at_end) {
                        return Ok(Progress::more(taken));
                    }
                }
                GzipState::Trailer(left) => {
                    if !self.skip_trailer(left, input, &mut taken, at_end) {
                        return Ok(Progress::more(taken));
                    }
                }
                GzipState::Inflate => {
                    // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1472-1489
                    // ```c++
                    //         // Prepare STREAM for decompressing
                    //         if ( m_Cache.size() ) {
                    //             // Possible, we have some unprocessed data in the cache
                    //             STREAM->next_in   = (unsigned char*)(m_Cache.data());
                    //             STREAM->avail_in  = (unsigned int)m_Cache.size();
                    //             ...
                    //             from_cache        = true;
                    //         } else {
                    //             STREAM->next_in   = (unsigned char*)x_in_buf;
                    //             STREAM->avail_in  = (unsigned int)x_in_len;
                    //             ...
                    //         }
                    //
                    //         // Try to decompress data
                    //         int errcode = Z(inflate)(STREAM, Z_SYNC_FLUSH);
                    // ```
                    // The cache first, then the input; any `inflate()` error is the stream's
                    // error (`eStatus_Error`, zlib.cpp:1565-1566). The inflater is called
                    // without input only at the end of the file (`Finish`).
                    let from_cache = self.cache_at < self.cache.len();
                    let data: &[u8] = if from_cache {
                        &self.cache[self.cache_at..]
                    } else {
                        &input[taken..]
                    };
                    if data.is_empty() && !at_end {
                        return Ok(Progress::more(taken));
                    }
                    let end_of_file = data.is_empty();
                    let (read_before, given_before) =
                        (self.inflater.total_in(), self.inflater.total_out());
                    let status = self
                        .inflater
                        .decompress(data, out, flate2::FlushDecompress::Sync);
                    let consumed = (self.inflater.total_in() - read_before) as usize;
                    let produced = (self.inflater.total_out() - given_before) as usize;
                    if from_cache {
                        self.cache_at += consumed;
                    } else {
                        taken += consumed;
                    }
                    match status {
                        Err(_) => return Err(self.fault("corrupt")),
                        Ok(flate2::Status::StreamEnd) => self.stream_end(),
                        Ok(_) if consumed == 0 && produced == 0 => {
                            if end_of_file {
                                self.end_of_input()?;
                                continue;
                            }
                            // Input that the inflater neither reads nor answers with an error.
                            return Err(self.fault("corrupt"));
                        }
                        Ok(_) => {}
                    }
                    if produced > 0 {
                        return Ok(Progress::written(taken, produced));
                    }
                }
            }
        }
    }
}

/// The capacity of the `BufReader` through which `bzip2::read::MultiBzDecoder` reads its source
/// (`BufReader::new`, std's `DEFAULT_BUF_SIZE`).
const BZIP2_WINDOW: usize = 8 * 1024;

/// The streams of a `.bz2` file, one after the other, as pbzip2 and `cat a.bz2 b.bz2` write
/// them: the loop of bzip2 0.6.1's `bufread::MultiBzDecoder` (`BzDecoder` with `multi`), which
/// the command line ran before the inflaters were pushed (`read::MultiBzDecoder`, through a
/// `BufReader` of `BZIP2_WINDOW` bytes, kept here as `window`), on the same
/// `bzip2::Decompress` (libbz2-rs-sys, a Rust port of libbzip2).
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/bzip2.cpp:795-824
/// ```c++
///         STREAM->next_in   = const_cast<char*>(in_buf);
///         STREAM->avail_in  = (unsigned int)in_len;
///         STREAM->next_out  = out_buf;
///         STREAM->avail_out = (unsigned int)out_size;
///
///         int errcode = BZ2_bzDecompress(STREAM);
///         ...
///         if ( m_DecompressMode == eMode_Decompress ) {
///             SetError(errcode, GetBZip2ErrorDescription(errcode));
///             *in_avail  = STREAM->avail_in;
///             *out_avail = out_size - STREAM->avail_out;
///             ...
///             switch (errcode) {
///             case BZ_OK:
///                 return eStatus_Success;
///             case BZ_STREAM_END:
///                 return eStatus_EndOfData;
///             }
/// ```
/// makeblastdb's `CBZip2Decompressor` ends the data at the first stream's end; LOSAT reads the
/// streams after it (`PD-LOSAT-NCBI-DEFECTS`, makeblastdb exception; line GZ, DECISIONS 9).
///
/// bzip2 0.6.1 reference: src/bufread.rs:178-221 (`impl<R: BufRead> Read for BzDecoder<R>`)
/// ```text
///         loop {
///             if self.done && !self.multi {
///                 return Ok(0);
///             }
///             let (read, consumed, remaining, ret);
///             {
///                 let input = self.obj.fill_buf()?;
///                 if self.done {
///                     assert!(self.multi);
///                     if input.is_empty() {
///                         // beyond last stream in multi-stream case
///                         return Ok(0);
///                     } else {
///                         // previous stream ended, more data follows => create new decompressor
///                         self.data = Decompress::new(false);
///                         self.done = false;
///                     }
///                 }
///                 let before_out = self.data.total_out();
///                 let before_in = self.data.total_in();
///                 ret = self.data.decompress(input, buf);
///                 read = (self.data.total_out() - before_out) as usize;
///                 consumed = (self.data.total_in() - before_in) as usize;
///                 remaining = input.len() - consumed;
///             }
///             self.obj.consume(consumed);
///
///             let ret = ret.map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
///             if ret == Status::StreamEnd {
///                 self.done = true;
///             } else if consumed == 0 && remaining == 0 && read == 0 {
///                 return Err(io::Error::new(
///                     io::ErrorKind::UnexpectedEof,
///                     "decompression not finished but EOF reached",
///                 ));
///             }
///
///             if read > 0 || buf.is_empty() {
///                 return Ok(read);
///             }
///         }
/// ```
struct Bzip2Streams {
    data: bzip2::Decompress,
    /// The stream ended (`done`): a new one starts where more input follows.
    done: bool,
    /// The bytes of the `BufReader`'s buffer, given to the decompressor from `window_at`; it is
    /// filled again (with up to `BZIP2_WINDOW` bytes) only when it is empty.
    window: Vec<u8>,
    window_at: usize,
}

impl Bzip2Streams {
    fn new() -> Self {
        Self {
            data: bzip2::Decompress::new(false),
            done: false,
            window: Vec::with_capacity(BZIP2_WINDOW),
            window_at: 0,
        }
    }

    /// One `Inflater::step` of a `.bz2` file: the loop of `BzDecoder::read`, up to the first
    /// output, the end of the last stream, or the end of `input` (all of it taken). An empty
    /// window at the end of the input is the end of the source (`input.is_empty()`).
    fn step(&mut self, input: &[u8], at_end: bool, out: &mut [u8]) -> io::Result<Progress> {
        let mut taken = 0;
        loop {
            if self.window_at == self.window.len() {
                let count = BZIP2_WINDOW.min(input.len() - taken);
                if count == 0 && !at_end {
                    return Ok(Progress::more(taken));
                }
                self.window.clear();
                self.window.extend_from_slice(&input[taken..taken + count]);
                self.window_at = 0;
                taken += count;
            }
            let window = &self.window[self.window_at..];
            if self.done {
                if window.is_empty() {
                    return Ok(Progress::ended(taken));
                }
                self.data = bzip2::Decompress::new(false);
                self.done = false;
            }
            let before_out = self.data.total_out();
            let before_in = self.data.total_in();
            let ret = self.data.decompress(window, out);
            let read = (self.data.total_out() - before_out) as usize;
            let consumed = (self.data.total_in() - before_in) as usize;
            let remaining = window.len() - consumed;
            self.window_at += consumed;

            let ret = ret.map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            if ret == bzip2::Status::StreamEnd {
                self.done = true;
            } else if consumed == 0 && remaining == 0 && read == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "decompression not finished but EOF reached",
                ));
            }

            if read > 0 {
                return Ok(Progress::written(taken, read));
            }
        }
    }
}

/// Where `ZstdFrames::gather` found the next unit whole.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Whole {
    /// In the piece of input, from the bytes taken so far: so many bytes.
    Input(usize),
    /// In `ZstdFrames::unit`.
    Gathered,
}

/// Where `ZstdFrames` is in the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ZstdState {
    /// Between frames: the next frame's header is gathered.
    Between,
    /// A skippable frame: so many bytes of it are still to pass over.
    Skip(u64),
    /// Inside a frame: its bytes are given out, then its next block is gathered.
    Frame,
}

/// The frames of a zstd file, read one after the other. ruzstd's `StreamingDecoder` reads one
/// frame and stops; pzstd and `cat a.zst b.zst` write several, pzstd with a skippable frame
/// before each, which this reader passes over (`ReadFrameHeaderError::SkipFrame`). At the end of
/// each frame the bytes are checked against the frame's `Content_Checksum` and
/// `Frame_Content_Size` (`check_frame`), which ruzstd reads but does not compare.
///
/// The input is pushed (`step`); ruzstd reads from a `Read`. Each call of ruzstd is made on a
/// slice that holds its whole unit: a frame header (`FrameDecoder::reset`) or a block with,
/// after the frame's last block, the frame's `Content_Checksum` (`FrameDecoder::decode_blocks`
/// with `UptoBlocks(1)`), sized from their first bytes (`frame_header_len`, `block_len`). The
/// slice is the piece of input itself where it holds the whole unit, else the unit gathered in
/// `unit` from the pieces (`gather`, `decode_unit`). At the end of the input a unit that is
/// not whole is given as it is, so that ruzstd fails on it as it failed on a source at its end.
/// (`FrameDecoder::decode_from_to`, ruzstd's own decoder of slices, waits on a block cut short
/// instead of failing.)
///
/// No NCBI counterpart: makeblastdb's `CZstdDecompressor` inflates with libzstd's
/// `ZSTD_decompressStream` and ends the data at the first frame's end; LOSAT reads the frames
/// after it (`PD-LOSAT-NCBI-DEFECTS`, makeblastdb exception; line GZ, DECISIONS 9).
/// NCBI reference (598d8ae6): c++/src/util/compress/api/zstd.cpp:937-967
/// ```c++
///         // Set I/O buffers for zstd stream.
///         ZSTD_inBuffer  m_In  = { in_buf, in_len, 0 };
///         ZSTD_outBuffer m_Out = { out_buf, out_size, 0 };
///
///         size_t result = ZSTD_decompressStream(DCTX, &m_Out, &m_In);
///         ...
///             *in_avail = m_In.size - m_In.pos;  // count unproc.bytes in input buffer
///             *out_avail = m_Out.pos;            // count bytes putted into out buffer
///             ...
///             if ( result == 0 ) {
///                 // frame is completely decoded and fully flushed
///                 return eStatus_EndOfData;
///             }
///             if ( !ZSTD_isError(result) ) {
///                 // decompressor still have some data
///                 return eStatus_Success;
///             }
/// ```
struct ZstdFrames {
    decoder: ruzstd::decoding::FrameDecoder,
    state: ZstdState,
    /// The bytes taken for the unit that ruzstd decodes next, where no piece of input held it
    /// whole.
    unit: Vec<u8>,
    /// The current frame ends with a `Content_Checksum`.
    checksum: bool,
    /// Bytes of the current frame given out so far.
    produced: u64,
    /// The decoder has read a frame header and reserved that frame's window. LOSAT: ruzstd
    /// reserves the ring buffer for the window only when a decoder that already has a state is
    /// reset, and a decoder's first frame creates its state without reserving, so the buffer
    /// grows by doubling as the blocks arrive. A native allocator hands the smaller regions back;
    /// a Wasm allocator cannot reuse them for the next larger one, and the linear memory, which
    /// never shrinks, ends near twice the buffer (the Web app's Data worker, `--long=27` file:
    /// 270.6 MiB where the window is 91.7 MiB). The first frame's header is therefore read a
    /// second time, with the `reset` of a decoder that has its state, which reserves what every
    /// later frame reserves. Every frame and input, no size rule.
    /// ruzstd 0.9.1 reference: src/decoding/frame_decoder.rs:200-210 (`FrameDecoder::reset`:
    /// ```text
    /// let state = match &mut self.state {
    ///     Some(s) => { s.reset(source, self.max_window_size)?; s }
    ///     None => { self.state = Some(FrameDecoderState::new(source, self.max_window_size)?); ... }
    /// ```
    /// ), 103-119 (`FrameDecoderState::new`: `DecoderScratch::new(window_size)`, whose
    /// `DecodeBuffer::new` reserves nothing, src/decoding/decode_buffer.rs:35-44), 121-130
    /// (`FrameDecoderState::reset`: `decoder_scratch.reset(window_size)`), decode_buffer.rs:46-49
    /// (`DecodeBuffer::reset`: `self.buffer.reserve(self.window_size)`), ringbuffer.rs:64-100
    /// (`reserve`, `reserve_amortized`: the doubling).
    reserved: bool,
}

/// The bytes of the frame header at the start of `unit`, as far as its first bytes tell: the
/// magic number (4 bytes; a skippable frame's 4-byte length after it, a wrong magic alone), then
/// the Frame_Header_Descriptor, the Window_Descriptor (without Single_Segment_Flag), the
/// Dictionary_ID and the Frame_Content_Size, as ruzstd reads them.
///
/// ruzstd 0.9.1 reference: src/decoding/frame.rs:6-85 (`read_frame_header`), 213-240
/// (`frame_content_size_bytes`, `dictionary_id_bytes`); RFC 8878 3.1.1.1 and 3.1.2.
fn frame_header_len(unit: &[u8]) -> usize {
    const MAGIC: u32 = 0xFD2F_B528;
    let Some(magic) = unit.first_chunk::<4>() else {
        return 4;
    };
    let magic = u32::from_le_bytes(*magic);
    if (0x184D_2A50..=0x184D_2A5F).contains(&magic) {
        return 8;
    }
    if magic != MAGIC {
        return 4;
    }
    let Some(&descriptor) = unit.get(4) else {
        return 5;
    };
    let single_segment = descriptor & 0x20 != 0;
    let window = usize::from(!single_segment);
    let dictionary_id = [0, 1, 2, 4][usize::from(descriptor & 0x3)];
    let content_size = match descriptor >> 6 {
        0 => usize::from(single_segment),
        1 => 2,
        2 => 4,
        _ => 8,
    };
    5 + window + dictionary_id + content_size
}

/// ruzstd's largest Block_Size (`MAX_BLOCK_SIZE`, 128 KiB): a block header that declares more is
/// refused from its 3 bytes.
const ZSTD_MAX_BLOCK_SIZE: usize = 128 * 1024;

/// The bytes of the block at the start of `unit`, as far as its first bytes tell: the 3-byte
/// Block_Header, the Block_Content (Block_Size bytes for a raw or compressed block, 1 for an
/// RLE block), then the frame's 4-byte Content_Checksum after its last block where the frame has
/// one. A header that ruzstd refuses from its 3 bytes, of the reserved type or with a Block_Size
/// over `ZSTD_MAX_BLOCK_SIZE`, is a unit of 3 bytes: it is refused as soon as it is whole, with
/// the error of a source that held the rest, without waiting for up to 2 MiB of input.
///
/// ruzstd 0.9.1 reference: src/decoding/block_decoder.rs:214-260 (`read_block_header`:
/// `content_size` = Block_Size for raw and compressed, 1 for RLE), 224-231 (the 3 bytes read,
/// then the reserved type and `block_content_size` refused), 283-290 (`block_content_size`:
/// `val > MAX_BLOCK_SIZE`), src/common/mod.rs:21 (`MAX_BLOCK_SIZE: u32 = 128 * 1024`),
/// src/decoding/frame_decoder.rs:347-358 (the checksum after the last block, in the same
/// `decode_blocks`); RFC 8878 3.1.1.2.
fn block_len(unit: &[u8], checksum: bool) -> usize {
    let Some(header) = unit.first_chunk::<3>() else {
        return 3;
    };
    let header = u32::from(header[0]) | (u32::from(header[1]) << 8) | (u32::from(header[2]) << 16);
    let last = header & 1 == 1;
    let block_size = (header >> 3) as usize;
    let block_type = (header >> 1) & 3;
    if block_type == 3 || block_size > ZSTD_MAX_BLOCK_SIZE {
        return 3;
    }
    let content = if block_type == 1 { 1 } else { block_size };
    3 + content + if last && checksum { 4 } else { 0 }
}

impl ZstdFrames {
    fn new() -> Self {
        Self {
            decoder: ruzstd::decoding::FrameDecoder::new(),
            state: ZstdState::Between,
            unit: Vec::new(),
            checksum: false,
            produced: 0,
            reserved: false,
        }
    }

    /// LOSAT: the checks of a frame's end that the zstd CLI makes (RFC 8878 3.1.1: the
    /// `Content_Checksum` is the low 4 bytes of XXH64 of the decoded data, and
    /// 3.1.1.1.4: `Frame_Content_Size` is the size of the decoded data). `zstd -d` refuses a
    /// mismatch ("Restored data doesn't match checksum", "Data corruption detected"); without
    /// the checks a damaged frame would be read as other letters, without a message.
    /// `content_size` is 0 when the header has no size or gives 0 (ruzstd does not tell them
    /// apart), which then is not compared: a frame without a checksum whose size reads 0 but
    /// that holds data passes, where `zstd -d` refuses it. The zstd CLI writes the checksum by
    /// default, and a frame with one is checked whatever its size field says.
    fn check_frame(&self) -> io::Result<()> {
        if let Some(stored) = self.decoder.get_checksum_from_data() {
            if self.decoder.get_calculated_checksum() != Some(stored) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "a zstd frame's checksum does not match its data",
                ));
            }
        }
        let declared = self.decoder.content_size();
        if declared != 0 && declared != self.produced {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "a zstd frame holds {} bytes where its header gives {declared}",
                    self.produced
                ),
            ));
        }
        Ok(())
    }

    /// Where the next unit is whole, or `None` while it waits for more input. Where nothing of
    /// the unit is gathered yet and the piece of input holds all of it, it is there
    /// (`Whole::Input`, its length from `taken`); else the input is moved into `unit` until that
    /// holds the whole unit (`Whole::Gathered`; at the end of the input, the unit as it is).
    /// `len` reads the unit's length from its first bytes.
    fn gather(
        &mut self,
        input: &[u8],
        taken: &mut usize,
        at_end: bool,
        len: impl Fn(&[u8]) -> usize,
    ) -> Option<Whole> {
        if self.unit.is_empty() {
            let rest = &input[*taken..];
            let want = len(rest);
            if rest.len() >= want {
                return Some(Whole::Input(want));
            }
        }
        loop {
            let want = len(&self.unit);
            if self.unit.len() >= want {
                return Some(Whole::Gathered);
            }
            let count = (want - self.unit.len()).min(input.len() - *taken);
            if count == 0 {
                return at_end.then_some(Whole::Gathered);
            }
            self.unit.extend_from_slice(&input[*taken..*taken + count]);
            *taken += count;
        }
    }

    /// The one call of ruzstd (`decode`) on the next unit, on the slice where `gather` found it
    /// whole; the bytes ruzstd read are then taken from there. `None` while the unit waits for
    /// more input.
    fn decode_unit<T>(
        &mut self,
        input: &[u8],
        taken: &mut usize,
        at_end: bool,
        len: impl Fn(&[u8]) -> usize,
        decode: impl FnOnce(&mut ruzstd::decoding::FrameDecoder, &mut &[u8]) -> T,
    ) -> Option<T> {
        let whole = self.gather(input, taken, at_end, len)?;
        let bytes: &[u8] = match whole {
            Whole::Input(unit_len) => &input[*taken..*taken + unit_len],
            Whole::Gathered => &self.unit,
        };
        let mut rest = bytes;
        let result = decode(&mut self.decoder, &mut rest);
        let used = bytes.len() - rest.len();
        match whole {
            Whole::Input(_) => *taken += used,
            Whole::Gathered => {
                self.unit.drain(..used);
            }
        }
        Some(result)
    }

    /// One `Inflater::step` of a `.zst` file: the frames, skippable frames and blocks up to the
    /// first output, the end of the file, or the end of `input` (all of it taken). The calls of
    /// ruzstd and their order are those of a reader on a source: `reset` at each frame (a
    /// skippable frame passed over), then for each block the bytes it lets out (`can_collect`,
    /// `read`), `check_frame` when the frame is finished, else `decode_blocks` of one block.
    fn step(&mut self, input: &[u8], at_end: bool, out: &mut [u8]) -> io::Result<Progress> {
        let mut taken = 0;
        loop {
            match self.state {
                ZstdState::Skip(left) => {
                    let count =
                        (input.len() - taken).min(usize::try_from(left).unwrap_or(usize::MAX));
                    taken += count;
                    let left = left - count as u64;
                    if left > 0 {
                        if at_end {
                            return Err(io::Error::new(
                                io::ErrorKind::UnexpectedEof,
                                "a skippable frame is cut short",
                            ));
                        }
                        self.state = ZstdState::Skip(left);
                        return Ok(Progress::more(taken));
                    }
                    self.state = ZstdState::Between;
                }
                ZstdState::Between => {
                    if self.unit.is_empty() && taken == input.len() {
                        return Ok(if at_end {
                            Progress::ended(taken)
                        } else {
                            Progress::more(taken)
                        });
                    }
                    let first = !self.reserved;
                    let Some((result, descriptor)) = self.decode_unit(
                        input,
                        &mut taken,
                        at_end,
                        frame_header_len,
                        |decoder, bytes| {
                            let descriptor = bytes.get(4).copied().unwrap_or(0);
                            let header = *bytes;
                            let mut result = decoder.reset(bytes);
                            if first && result.is_ok() {
                                // The first frame of a decoder: see `ZstdFrames::reserved`.
                                result = decoder.reset(header);
                            }
                            (result, descriptor)
                        },
                    ) else {
                        return Ok(Progress::more(taken));
                    };
                    match result {
                        Ok(()) => {
                            self.reserved = true;
                            self.state = ZstdState::Frame;
                            self.checksum = descriptor & 0x04 != 0;
                            self.produced = 0;
                        }
                        Err(ruzstd::decoding::errors::FrameDecoderError::ReadFrameHeaderError(
                            ruzstd::decoding::errors::ReadFrameHeaderError::SkipFrame {
                                length,
                                ..
                            },
                        )) => self.state = ZstdState::Skip(u64::from(length)),
                        Err(error) => return Err(invalid_zstd(error)),
                    }
                }
                ZstdState::Frame => {
                    if self.decoder.can_collect() > 0 {
                        let count = self.decoder.read(out)?;
                        self.produced += count as u64;
                        return Ok(Progress::written(taken, count));
                    }
                    if self.decoder.is_finished() {
                        self.check_frame()?;
                        self.state = ZstdState::Between;
                        continue;
                    }
                    let checksum = self.checksum;
                    let Some(result) = self.decode_unit(
                        input,
                        &mut taken,
                        at_end,
                        |unit| block_len(unit, checksum),
                        |decoder, bytes| {
                            decoder.decode_blocks(
                                bytes,
                                ruzstd::decoding::BlockDecodingStrategy::UptoBlocks(1),
                            )
                        },
                    ) else {
                        return Ok(Progress::more(taken));
                    };
                    result.map_err(invalid_zstd)?;
                }
            }
        }
    }
}

fn invalid_zstd(error: ruzstd::decoding::errors::FrameDecoderError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

#[cfg(test)]
pub(crate) mod test_data {
    //! Compressed test inputs made without a compressor binary: flate2 and bzip2 encode; zstd,
    //! which LOSAT only decodes, is written as frames of raw (stored) blocks, which every zstd
    //! decoder reads.
    use super::{Compression, PushInflater};
    use std::io::{Read, Write};

    /// The sizes of the pieces in which the tests push a file: one byte, a few, around the 4 KiB
    /// of a gzip header check, the 64 KiB blocks of the command line's reader, and all at once.
    pub(crate) const PIECES: [usize; 7] = [1, 3, 4095, 4096, 4097, 65536, usize::MAX];

    /// What an inflation gives: the inflated bytes, or the error's kind and words.
    pub(crate) type Outcome = Result<Vec<u8>, (std::io::ErrorKind, String)>;

    fn outcome(result: std::io::Result<Vec<u8>>) -> Outcome {
        result.map_err(|error| (error.kind(), error.to_string()))
    }

    fn describe(outcome: &Outcome) -> String {
        match outcome {
            Ok(bytes) => format!("{} bytes", bytes.len()),
            Err((kind, words)) => format!("{kind:?} {words:?}"),
        }
    }

    /// What the command line's reader (`Compression::decoder`) gives for `bytes`.
    pub(crate) fn read(compression: Compression, bytes: &[u8]) -> Outcome {
        let mut out = Vec::new();
        outcome(
            compression
                .decoder(bytes)
                .read_to_end(&mut out)
                .map(|_| out),
        )
    }

    /// What `PushInflater` gives for `bytes` pushed in pieces of `piece` bytes, then finished.
    pub(crate) fn pushed(compression: Compression, bytes: &[u8], piece: usize) -> Outcome {
        let mut inflater = PushInflater::new(compression);
        let mut out = Vec::new();
        let mut sink = |inflated: &[u8]| {
            assert!(!inflated.is_empty() && inflated.len() <= PushInflater::OUTPUT_SIZE);
            out.extend_from_slice(inflated);
        };
        let result = bytes
            .chunks(piece)
            .try_for_each(|chunk| inflater.push(chunk, &mut sink))
            .and_then(|()| inflater.finish(&mut sink));
        outcome(result.map(|()| out))
    }

    /// Checks that `bytes` pushed in pieces of every size of `PIECES` give what the command
    /// line's reader gave (`read`): the same bytes, or an error of the same kind and words.
    pub(crate) fn assert_pushed_as_read(
        compression: Compression,
        bytes: &[u8],
        read: &std::io::Result<Vec<u8>>,
    ) {
        let read = match read {
            Ok(bytes) => Ok(bytes.clone()),
            Err(error) => Err((error.kind(), error.to_string())),
        };
        assert!(
            !matches!(&read, Err((_, words)) if words == super::STALLED),
            "{compression:?}: the reader stalled"
        );
        for piece in PIECES {
            let pushed = pushed(compression, bytes, piece);
            assert!(
                pushed == read,
                "{compression:?} in pieces of {piece}: pushed {}, read {}",
                describe(&pushed),
                describe(&read)
            );
        }
    }

    /// Gzip members, each with bgzip's extra field, concatenated as bgzip and `cat` write them.
    pub(crate) fn gzip_members(members: &[&[u8]]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for member in members {
            let mut encoder = flate2::GzBuilder::new()
                .extra(b"BC\x02\x00\x00\x00".to_vec())
                .write(&mut bytes, flate2::Compression::default());
            encoder.write_all(member).unwrap();
            encoder.finish().unwrap();
        }
        bytes
    }

    /// Bzip2 streams, concatenated as pbzip2 and `cat` write them.
    pub(crate) fn bzip2_streams(streams: &[&[u8]]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for stream in streams {
            let mut encoder =
                bzip2::write::BzEncoder::new(&mut bytes, bzip2::Compression::default());
            encoder.write_all(stream).unwrap();
            encoder.finish().unwrap();
        }
        bytes
    }

    /// Zstd frames of raw blocks (Single_Segment, a 4-byte Frame_Content_Size, no checksum), each
    /// preceded by a skippable frame as pzstd writes one, concatenated.
    pub(crate) fn zstd_frames(frames: &[&[u8]]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for frame in frames {
            // A skippable frame: magic 0x184D2A50, a 4-byte length, the payload.
            bytes.extend_from_slice(&0x184D_2A50u32.to_le_bytes());
            bytes.extend_from_slice(&4u32.to_le_bytes());
            bytes.extend_from_slice(b"skip");
            bytes.extend_from_slice(&0xFD2F_B528u32.to_le_bytes());
            // Frame_Header_Descriptor: Frame_Content_Size_Flag 2 (4 bytes), Single_Segment_Flag.
            bytes.push(0x80 | 0x20);
            bytes.extend_from_slice(&(frame.len() as u32).to_le_bytes());
            let mut blocks = frame.chunks(128 * 1024).peekable();
            if blocks.peek().is_none() {
                bytes.extend_from_slice(&[1, 0, 0]);
            }
            while let Some(block) = blocks.next() {
                let last = u32::from(blocks.peek().is_none());
                let header = last | ((block.len() as u32) << 3);
                bytes.extend_from_slice(&header.to_le_bytes()[..3]);
                bytes.extend_from_slice(block);
            }
        }
        bytes
    }
}

#[cfg(test)]
#[path = "compressed_push_tests.rs"]
mod push_tests;

#[cfg(test)]
mod tests {
    use super::test_data::{bzip2_streams, gzip_members, zstd_frames};
    use super::{compression_of_name, Compression};
    use std::io::Read;
    use std::path::Path;

    /// The inflated bytes, as the command line's reader gives them; the pushed inflater must
    /// give the same in pieces of any size (`assert_pushed_as_read`).
    fn inflate(compression: Compression, bytes: &[u8]) -> std::io::Result<Vec<u8>> {
        let mut out = Vec::new();
        let result = compression
            .decoder(bytes)
            .read_to_end(&mut out)
            .map(|_| out);
        super::test_data::assert_pushed_as_read(compression, bytes, &result);
        result
    }

    // NCBI reference: c++/src/app/blastdb/makeblastdb.cpp:1037-1039 (`NStr::EndsWith(seq_file,
    // ".gz" / ".bz2" / ".zst", NStr::eNocase)`): the suffix in any case, nothing else.
    #[test]
    fn a_name_announces_its_compression_by_its_suffix_in_any_case() {
        let cases = [
            (".gz", Some(Compression::Gzip)),
            ("q.fa.GZ", Some(Compression::Gzip)),
            ("dir.gz/q.fa.gz", Some(Compression::Gzip)),
            ("q.fa.bz2", Some(Compression::Bzip2)),
            ("q.fa.BZ2", Some(Compression::Bzip2)),
            ("q.fa.zst", Some(Compression::Zstd)),
            ("q.fa.Zst", Some(Compression::Zstd)),
            ("-", None),
            ("q.fa", None),
            ("q.gzip", None),
            ("q.gz.fa", None),
            ("q.zstd", None),
            ("q.bz", None),
            ("dir.gz/q.fa", None),
            ("", None),
        ];
        for (name, expected) in cases {
            assert_eq!(compression_of_name(Path::new(name)), expected, "{name}");
        }
    }

    // No NCBI counterpart for bzip2 and zstd (line GZ, DECISIONS 4 and 9): every member, stream
    // and frame is read, skippable frames passed over, an empty member included.
    #[test]
    fn every_member_stream_and_frame_is_read() {
        let parts: [&[u8]; 3] = [b">q1\nACGT\n", b"", b">q2\nTTTT\n"];
        assert_eq!(
            inflate(Compression::Gzip, &gzip_members(&parts)).unwrap(),
            b">q1\nACGT\n>q2\nTTTT\n"
        );
        assert_eq!(
            inflate(Compression::Bzip2, &bzip2_streams(&parts)).unwrap(),
            b">q1\nACGT\n>q2\nTTTT\n"
        );
        assert_eq!(
            inflate(Compression::Zstd, &zstd_frames(&parts)).unwrap(),
            b">q1\nACGT\n>q2\nTTTT\n"
        );
        let long = vec![b'A'; 300 * 1024];
        assert_eq!(
            inflate(Compression::Zstd, &zstd_frames(&[&long])).unwrap(),
            long
        );
    }

    // No NCBI counterpart for bzip2 and zstd (line GZ, DECISIONS 5): text under a compressed name,
    // data cut short, or an empty file is an error of the read. For gzip these are faults that
    // make makeblastdb's stream fail too (the tests of `GzipMembers` below).
    #[test]
    fn text_and_cut_data_are_errors_of_the_read() {
        let text: &[u8] = b">q1\nACGT\n";
        for compression in Compression::ALL {
            assert!(inflate(compression, text).is_err(), "{compression:?}: text");
            let empty = inflate(compression, b"").unwrap_err();
            assert_eq!(empty.to_string(), "the file is empty", "{compression:?}");
        }
        let whole: &[u8] = b">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let cut = |bytes: Vec<u8>| bytes[..bytes.len() - 12].to_vec();
        assert!(inflate(Compression::Gzip, &cut(gzip_members(&[whole]))).is_err());
        assert!(inflate(Compression::Bzip2, &cut(bzip2_streams(&[whole]))).is_err());
        assert!(inflate(Compression::Zstd, &cut(zstd_frames(&[whole]))).is_err());
    }

    /// The tests of `GzipMembers` read the edge cases of a `.gz` file that NCBI's makeblastdb
    /// 2.17.0 (the official Linux x64 build, zlib of the system) was run on on 2026-11-08
    /// (`makeblastdb -in <case>.fa.gz -dbtype nucl`, then `blastdbcmd -entry all`); each
    /// expectation is that run's: the stored sequences where it built the database (exit 0, no
    /// message), an error where its decompressor failed (`[CZipDecompressor::Process]  data
    /// error` or `buffer error`, exit 6), and no bytes where it found no sequence ("No sequences
    /// added"). The small cases are built from the deflate data that zlib (and `gzip -n`) wrote
    /// for the oracle's files, byte for byte.
    mod gzip {
        use super::super::{check_gzip_header, Compression, MAX_HEADER_SIZE};
        use std::io::{Read, Write};

        const SMALL: &[u8] = b">q1 test\nACGTACGTACGTTTGACCAGTACGATCGATCGACTAGCTAGCATCGACTGACTAG\n";
        const SECOND: &[u8] = b">q2 second\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n";
        /// zlib's (level 6) raw deflate data of `SMALL`, as in the oracle's `good.fa.gz`.
        const SMALL_DEFLATE: [u8; 39] = [
            0xb3, 0x2b, 0x34, 0x54, 0x28, 0x49, 0x2d, 0x2e, 0xe1, 0x72, 0x74, 0x76, 0x0f, 0x81,
            0xe1, 0x90, 0x10, 0x77, 0x47, 0x67, 0x67, 0x47, 0x30, 0xc7, 0x31, 0x04, 0x8a, 0x9d,
            0x43, 0x1c, 0xdd, 0xc1, 0x18, 0xca, 0x83, 0x88, 0x70, 0x01, 0x00,
        ];
        /// zlib's raw deflate data of `SECOND`.
        const SECOND_DEFLATE: [u8; 26] = [
            0xb3, 0x2b, 0x34, 0x52, 0x28, 0x4e, 0x4d, 0xce, 0xcf, 0x4b, 0xe1, 0x0a, 0x01, 0x02,
            0x77, 0x20, 0x70, 0x06, 0x02, 0x47, 0x20, 0x40, 0xe7, 0x73, 0x01, 0x00,
        ];
        /// Adler-32 of `SMALL` (the trailer of the oracle's `zlib.fa.gz`, big-endian).
        const SMALL_ADLER: [u8; 4] = [0x68, 0x91, 0x12, 0x2d];
        /// An empty stored block that is not the last one: no output, five bytes of input.
        const EMPTY_STORED: [u8; 5] = [0x00, 0x00, 0x00, 0xff, 0xff];

        fn crc32(data: &[u8]) -> u32 {
            let mut crc = flate2::Crc::new();
            crc.update(data);
            crc.sum()
        }

        /// A gzip header as `gzip -n` writes it (mtime 0, OS 3), with the given fields.
        fn header(flags: u8, method: u8, extra: Option<&[u8]>, name: Option<&[u8]>) -> Vec<u8> {
            let mut h = vec![0x1f, 0x8b, method, flags, 0, 0, 0, 0, 0, 3];
            if let Some(extra) = extra {
                h.extend_from_slice(&(extra.len() as u16).to_le_bytes());
                h.extend_from_slice(extra);
            }
            if let Some(name) = name {
                h.extend_from_slice(name);
                h.push(0);
            }
            h
        }

        fn trailer(data: &[u8]) -> Vec<u8> {
            let mut t = crc32(data).to_le_bytes().to_vec();
            t.extend_from_slice(&(data.len() as u32).to_le_bytes());
            t
        }

        fn small() -> Vec<u8> {
            [
                &header(0, 8, None, None)[..],
                &SMALL_DEFLATE,
                &trailer(SMALL),
            ]
            .concat()
        }

        fn second() -> Vec<u8> {
            [
                &header(0, 8, None, None)[..],
                &SECOND_DEFLATE,
                &trailer(SECOND),
            ]
            .concat()
        }

        fn with_header(h: Vec<u8>) -> Vec<u8> {
            [&h[..], &SMALL_DEFLATE, &trailer(SMALL)].concat()
        }

        /// 300 000 letters in three records of 60-letter lines, as one member (any deflate
        /// encoder: the cases read past the end of its data).
        fn big_plain() -> Vec<u8> {
            let mut state = 1u64;
            let mut text = Vec::new();
            for record in 1..=3 {
                text.extend_from_slice(format!(">big{record} record {record}\n").as_bytes());
                for i in 0..100_000 {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    text.push(b"ACGT"[(state >> 62) as usize]);
                    if i % 60 == 59 || i == 99_999 {
                        text.push(b'\n');
                    }
                }
            }
            text
        }

        fn member_of(data: &[u8]) -> Vec<u8> {
            let mut encoder =
                flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(data).unwrap();
            let deflate = encoder.finish().unwrap();
            [&header(0, 8, None, None)[..], &deflate, &trailer(data)].concat()
        }

        /// A source that gives one byte per read, so that the header window is filled from many
        /// reads of the source.
        struct ByteByByte<'a>(&'a [u8]);

        impl Read for ByteByByte<'_> {
            fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
                if buf.is_empty() || self.0.is_empty() {
                    return Ok(0);
                }
                buf[0] = self.0[0];
                self.0 = &self.0[1..];
                Ok(1)
            }
        }

        /// The inflated bytes, read whole and also a byte at a time from a source read a byte
        /// at a time, which must agree: the result does not hang on the sizes of the reads. The
        /// pushed inflater must give the same in pieces of any size (`assert_pushed_as_read`).
        fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, String> {
            let mut whole = Vec::new();
            let whole = Compression::Gzip
                .decoder(bytes)
                .read_to_end(&mut whole)
                .map(|_| whole);
            super::super::test_data::assert_pushed_as_read(Compression::Gzip, bytes, &whole);
            let whole = whole.map_err(|error| error.to_string());
            let mut decoder = Compression::Gzip.decoder(ByteByByte(bytes));
            let mut slow = Vec::new();
            let mut byte = [0u8; 1];
            let slow = loop {
                match decoder.read(&mut byte) {
                    Ok(0) => break Ok(slow),
                    Ok(_) => slow.push(byte[0]),
                    Err(error) => break Err(error.to_string()),
                }
            };
            assert_eq!(whole, slow, "the reads' sizes changed the result");
            whole
        }

        fn ok(name: &str, bytes: &[u8], expected: &[u8]) {
            match gunzip(bytes) {
                Ok(out) => assert!(out == expected, "{name}: inflated {} bytes", out.len()),
                Err(error) => panic!("{name}: makeblastdb reads it, LOSAT refused: {error}"),
            }
        }

        fn refused(name: &str, bytes: &[u8]) -> String {
            match gunzip(bytes) {
                Ok(out) => panic!("{name}: makeblastdb fails, LOSAT read {} bytes", out.len()),
                Err(error) => error,
            }
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1513-1526: the trailer is
        // skipped (`STREAM->avail_in -= 8`, or `m_SkipInput` for the bytes the input lacks), never
        // compared. Oracle cases good, badcrc, badlen, cut1-cut8 (built), cut9 and cut12 (buffer
        // error).
        #[test]
        fn a_trailer_is_skipped_unchecked_and_may_be_cut_short() {
            let good = small();
            ok("good", &good, SMALL);
            let mut badcrc = good.clone();
            let crc_at = good.len() - 8;
            badcrc[crc_at..crc_at + 4].copy_from_slice(&(crc32(SMALL) ^ 0xFFFF_FFFF).to_le_bytes());
            ok("badcrc", &badcrc, SMALL);
            let mut badlen = good.clone();
            badlen[crc_at + 4..].copy_from_slice(&(SMALL.len() as u32 + 1).to_le_bytes());
            ok("badlen", &badlen, SMALL);
            for cut in 1..=8 {
                ok(&format!("cut{cut}"), &good[..good.len() - cut], SMALL);
            }
            for cut in [9, 12] {
                let error = refused(&format!("cut{cut}"), &good[..good.len() - cut]);
                assert_eq!(error, "incomplete deflate stream", "cut{cut}");
            }
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:151-209
        // (`s_CheckGZipHeader`). Oracle cases fhcrc_bad, fhcrc_ok, fextra, fname, fcomment, fall,
        // ftext, fname4000 (built); res20, res40, res80, cm7, fname5000, fextra5000, tiny2
        // (data error), hdr_only, hdr_plus2 (buffer error); hdr_plus1 ("No sequences added").
        #[test]
        fn a_header_is_read_as_makeblastdb_reads_it() {
            let mut fhcrc_ok = header(0x02, 8, None, None);
            let crc16 = (crc32(&fhcrc_ok) & 0xFFFF) as u16;
            fhcrc_ok.extend_from_slice(&crc16.to_le_bytes());
            ok("fhcrc_ok", &with_header(fhcrc_ok), SMALL);
            let mut fhcrc_bad = header(0x02, 8, None, None);
            fhcrc_bad.extend_from_slice(&0x1234u16.to_le_bytes());
            ok("fhcrc_bad", &with_header(fhcrc_bad), SMALL);
            ok(
                "fextra",
                &with_header(header(0x04, 8, Some(b"AB\x03\x00xyz"), None)),
                SMALL,
            );
            ok(
                "fname",
                &with_header(header(0x08, 8, None, Some(b"a.fa"))),
                SMALL,
            );
            ok(
                "fcomment",
                &with_header(header(0x10, 8, None, Some(b"a comment"))),
                SMALL,
            );
            let mut fall = header(0x1e, 8, Some(b"AB\x01\x00z"), Some(b"a.fa"));
            fall.extend_from_slice(b"c\0");
            fall.extend_from_slice(&0xBEEFu16.to_le_bytes());
            ok("fall", &with_header(fall), SMALL);
            ok("ftext", &with_header(header(0x01, 8, None, None)), SMALL);
            let long = |n: usize| vec![b'n'; n];
            ok(
                "fname4000",
                &with_header(header(0x08, 8, None, Some(&long(4000)))),
                SMALL,
            );
            for (name, h) in [
                ("res20", header(0x20, 8, None, None)),
                ("res40", header(0x40, 8, None, None)),
                ("res80", header(0x80, 8, None, None)),
                ("cm7", header(0, 7, None, None)),
                ("fextra5000", header(0x04, 8, Some(&[b'E'; 5000]), None)),
            ] {
                let error = refused(name, &with_header(h));
                assert_eq!(error, "invalid gzip header", "{name}");
            }
            let error = refused(
                "fname5000",
                &with_header(header(0x08, 8, None, Some(&long(5000)))),
            );
            assert_eq!(error, "corrupt deflate stream");
            assert_eq!(refused("tiny2", &[0x1f, 0x8b]), "invalid gzip header");
            let good = small();
            assert_eq!(
                refused("hdr_only", &good[..10]),
                "incomplete deflate stream"
            );
            ok("hdr_plus1", &good[..11], b"");
            assert_eq!(
                refused("hdr_plus2", &good[..12]),
                "incomplete deflate stream"
            );
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:151-209: the header
        // parser's lengths, at the edges of its window.
        #[test]
        fn the_header_parser_stops_at_the_end_of_its_window() {
            let h = header(0x08, 8, None, Some(b"name"));
            assert_eq!(check_gzip_header(&h), 15);
            assert_eq!(
                check_gzip_header(&h[..14]),
                14,
                "a name without its NUL ends at the end"
            );
            assert_eq!(check_gzip_header(&h[..9]), 0, "fewer than 10 bytes");
            let mut h = header(0x04, 8, Some(b"xy"), None);
            assert_eq!(check_gzip_header(&h), 14);
            assert_eq!(
                check_gzip_header(&h[..13]),
                0,
                "the extra field passes the end"
            );
            assert_eq!(check_gzip_header(&h[..11]), 0, "the length passes the end");
            h[3] |= 0x02;
            assert_eq!(check_gzip_header(&h), 0, "the header CRC passes the end");
            h.extend_from_slice(&[0, 0]);
            assert_eq!(check_gzip_header(&h), 16);
            let name = [&header(0x08, 8, None, None)[..], &[b'n'; MAX_HEADER_SIZE]].concat();
            assert_eq!(check_gzip_header(&name[..MAX_HEADER_SIZE]), MAX_HEADER_SIZE);
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1437-1447 (`m_IsGZ` is
        // set by a header and never cleared: a later member whose header is not found is
        // inflated as raw deflate) and 1522-1525 (`fAllowConcatenatedGZip`). Oracle cases
        // concat2, concat_empty, bgzf_small, m2_fhcrc_bad, m2_badcrc, m1_badcrc_m2 (built);
        // m2_hdr_plus1 (built from the first member); m2_res20, m2_badmagic, m2_cm7,
        // m1_cut3_m2 (data error), m2_hdr_plus3 (buffer error).
        #[test]
        fn every_member_is_read_and_a_later_bad_header_is_read_as_deflate() {
            let both = [SMALL, SECOND].concat();
            let (first, second) = (small(), second());
            ok("concat2", &[&first[..], &second].concat(), &both);
            let empty = [&header(0, 8, None, None)[..], &[0x03, 0x00], &trailer(b"")].concat();
            ok(
                "concat_empty",
                &[&first[..], &empty, &second].concat(),
                &both,
            );
            let bc = |data: &[u8], deflate: &[u8]| {
                [
                    &header(0x04, 8, Some(b"BC\x02\x00\x00\x00"), None)[..],
                    deflate,
                    &trailer(data),
                ]
                .concat()
            };
            ok(
                "bgzf_small",
                &[bc(SMALL, &SMALL_DEFLATE), bc(b"", &[0x03, 0x00])].concat(),
                SMALL,
            );
            let mut m2_fhcrc_bad = header(0x02, 8, None, None);
            m2_fhcrc_bad.extend_from_slice(&1u16.to_le_bytes());
            let m2_fhcrc_bad = [&m2_fhcrc_bad[..], &SECOND_DEFLATE, &trailer(SECOND)].concat();
            ok("m2_fhcrc_bad", &[&first[..], &m2_fhcrc_bad].concat(), &both);
            let mut m2_badcrc = second.clone();
            let at = m2_badcrc.len() - 8;
            m2_badcrc[at..at + 4].copy_from_slice(&1u32.to_le_bytes());
            ok("m2_badcrc", &[&first[..], &m2_badcrc].concat(), &both);
            let mut m1_badcrc = first.clone();
            let at = m1_badcrc.len() - 8;
            m1_badcrc[at..at + 4].copy_from_slice(&1u32.to_le_bytes());
            ok("m1_badcrc_m2", &[&m1_badcrc[..], &second].concat(), &both);
            ok("m2_hdr_plus1", &[&first[..], &second[..11]].concat(), SMALL);
            let bad_second =
                |h: Vec<u8>| [&first[..], &h, &SECOND_DEFLATE, &trailer(SECOND)].concat();
            for (name, bytes) in [
                ("m2_res20", bad_second(header(0x20, 8, None, None))),
                ("m2_cm7", bad_second(header(0, 7, None, None))),
                (
                    "m2_badmagic",
                    [&first[..], &[0x1f, 0x8c], &second[2..]].concat(),
                ),
                ("m1_cut3_m2", [&first[..first.len() - 3], &second].concat()),
            ] {
                assert_eq!(refused(name, &bytes), "corrupt deflate stream", "{name}");
            }
            let error = refused("m2_hdr_plus3", &[&first[..], &second[..13]].concat());
            assert_eq!(error, "incomplete deflate stream");
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1617-1630 (`Finish`:
        // a `Process` without output ends the data) with 1412-1431 (the rest of a file shorter
        // than `kMaxHeaderSize` is held for the header check until the input ends). Oracle
        // cases zeros1, zeros2, zeros4, junk_nl, junk_crlf, esb819, big_zeros3, big_esb819,
        // bgzf_big_zeros3 (built); zeros5, zeros8, zeros512, zeros5000, junk_hdr3, big_zeros5,
        // big_zeros5000 (data error); junk_text, junk_hdr10, esb819z, esb820, big_junk,
        // big_esb820 (buffer error).
        #[test]
        fn bytes_after_the_last_member_are_read_as_makeblastdb_reads_them() {
            let good = small();
            let after = |tail: &[u8]| [&good[..], tail].concat();
            for zeros in [1, 2, 4] {
                ok(&format!("zeros{zeros}"), &after(&vec![0; zeros]), SMALL);
            }
            ok("junk_nl", &after(b"\n"), SMALL);
            ok("junk_crlf", &after(b"\r\n"), SMALL);
            let stored = |count: usize| EMPTY_STORED.repeat(count);
            assert_eq!(stored(819).len(), MAX_HEADER_SIZE - 1);
            ok("esb819", &after(&stored(819)), SMALL);
            for zeros in [5, 8, 512, 5000] {
                let error = refused(&format!("zeros{zeros}"), &after(&vec![0; zeros]));
                assert_eq!(error, "corrupt deflate stream", "zeros{zeros}");
            }
            assert_eq!(
                refused("junk_hdr3", &after(&[0x1f, 0x8b, 0x08])),
                "corrupt deflate stream"
            );
            for (name, tail) in [
                ("junk_text", b"junk\n".to_vec()),
                ("junk_hdr10", header(0, 8, None, None)),
                ("esb819z", [&stored(819)[..], &[0]].concat()),
                ("esb820", stored(820)),
            ] {
                assert_eq!(
                    refused(name, &after(&tail)),
                    "incomplete deflate stream",
                    "{name}"
                );
            }
            let plain = big_plain();
            let big = member_of(&plain);
            let big_after = |tail: &[u8]| [&big[..], tail].concat();
            ok("big_zeros3", &big_after(&[0; 3]), &plain);
            ok("big_esb819", &big_after(&stored(819)), &plain);
            ok("big_cut4", &big[..big.len() - 4], &plain);
            let mut big_badcrc = big.clone();
            let at = big.len() - 8;
            big_badcrc[at..at + 4].copy_from_slice(&0u32.to_le_bytes());
            ok("big_badcrc", &big_badcrc, &plain);
            assert_eq!(
                refused("big_zeros5", &big_after(&[0; 5])),
                "corrupt deflate stream"
            );
            assert_eq!(
                refused("big_zeros5000", &big_after(&[0; 5000])),
                "corrupt deflate stream"
            );
            assert_eq!(
                refused("big_junk", &big_after(b"junk\n")),
                "incomplete deflate stream"
            );
            assert_eq!(
                refused("big_esb820", &big_after(&stored(820))),
                "incomplete deflate stream"
            );
            assert_eq!(
                refused("big_cut100", &big[..big.len() - 100]),
                "incomplete deflate stream"
            );
            // bgzip's blocks of 65280 bytes, its empty end block, then three zero bytes.
            let blocks: Vec<u8> = plain
                .chunks(65280)
                .chain(std::iter::once(&b""[..]))
                .flat_map(member_of)
                .collect();
            ok("bgzf_big_zeros3", &[&blocks[..], &[0; 3]].concat(), &plain);
            // A second member's name that passes the 4 KiB window (oracle big_m2_fname5000: data
            // error; big_m2_fname4000 built).
            let named = |n: usize| {
                let h = header(0x08, 8, None, Some(&vec![b'n'; n]));
                [&big[..], &h, &SECOND_DEFLATE, &trailer(SECOND)].concat()
            };
            ok(
                "big_m2_fname4000",
                &named(4000),
                &[&plain[..], SECOND].concat(),
            );
            assert_eq!(
                refused("big_m2_fname5000", &named(5000)),
                "corrupt deflate stream"
            );
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1443-1447 (without a
        // gzip header, `inflateInit2_` gets positive window bits: a zlib stream, Adler-32
        // checked) and 1562-1563 (its `Z_STREAM_END` ends the data). Oracle cases zlib,
        // zlib_junk, zlib_then_gz (built, the zlib stream's record only); zlib_badadler,
        // zlib_cinfo8, zlib_fdict, notgz (data error, "need dictionary"); onebyte,
        // zlib_hdr_only ("No sequences added").
        #[test]
        fn a_file_without_a_gzip_header_is_read_as_a_zlib_stream() {
            let zlib = [&[0x78, 0x9c][..], &SMALL_DEFLATE, &SMALL_ADLER].concat();
            ok("zlib", &zlib, SMALL);
            ok(
                "zlib_junk",
                &[&zlib[..], b"junk junk junk\n"].concat(),
                SMALL,
            );
            ok("zlib_then_gz", &[&zlib[..], &second()].concat(), SMALL);
            let mut badadler = zlib.clone();
            *badadler.last_mut().unwrap() ^= 1;
            assert_eq!(refused("zlib_badadler", &badadler), "corrupt zlib stream");
            for (name, cmf, flg) in [("zlib_cinfo8", 0x88, 0x1c), ("zlib_fdict", 0x78, 0x20)] {
                assert_eq!(((u16::from(cmf) << 8) | u16::from(flg)) % 31, 0, "{name}");
                let bytes = [&[cmf, flg][..], &SMALL_DEFLATE, &SMALL_ADLER].concat();
                assert_eq!(refused(name, &bytes), "invalid gzip header", "{name}");
            }
            assert_eq!(refused("notgz", SMALL), "invalid gzip header");
            ok("onebyte", b"x", b"");
            ok("zlib_hdr_only", &[0x78, 0x9c], b"");
        }

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1444-1447, 1489
        // ```c++
        //                 Z(inflateEnd)(STREAM);
        //                 int errcode = Z(inflateInit2_)(
        //                                 STREAM, m_IsGZ ? -GetWindowBits() : GetWindowBits(),
        //                                 ZLIB_VERSION, (int)sizeof(z_stream));
        // ...
        //         int errcode = Z(inflate)(STREAM, Z_SYNC_FLUSH);
        // ```
        // makeblastdb inflates with zlib, a new stream for each member: zlib's checks decide what
        // is read (zlib-rs makes them as zlib does; miniz_oxide, flate2's default, did not). Oracle
        // cases of the review of 2026-11-08 (swebok-improve SW8 DECISIONS 41): far_title (a
        // distance past the bytes written so far), far_m2 (a second member's match reaching back
        // into the first member), tail_noeob_aligned (a dynamic table without an end-of-block
        // code), tail_sym16_full (code-length symbol 16 with no length before it): data error;
        // tail_hclen_empty (a code-length code without lengths, the file ending before the code
        // lengths) and tail_sym16_cut (symbol 16, the file ending before its 2 extra bits):
        // built, the first member's record.
        #[test]
        fn zlib_checks_decide_what_is_read() {
            let far_title = [
                0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xb3, 0x2b, 0x34, 0x54,
                0xa8, 0x00, 0xea, 0xe0, 0x72, 0x74, 0x76, 0x0f, 0x71, 0x74, 0x76, 0x0f, 0x71, 0x74,
                0x76, 0x0f, 0x09, 0x09, 0x71, 0x77, 0x74, 0x76, 0x76, 0x74, 0x0f, 0x71, 0x74, 0x76,
                0x77, 0x0c, 0x71, 0x76, 0xe7, 0x02, 0x00, 0xa7, 0x7d, 0x49, 0x5f, 0x26, 0x00, 0x00,
                0x00,
            ];
            assert_eq!(refused("far_title", &far_title), "corrupt deflate stream");
            // ">q2 y\n", then a match of 10 bytes at distance 40 (past the 6 bytes of this member),
            // then text; the trailer is zero.
            let m2_deflate = [
                0xb3, 0x2b, 0x34, 0x52, 0xa8, 0xe4, 0x42, 0xa8, 0x73, 0x74, 0x76, 0x0f, 0x71, 0x74,
                0x76, 0x0f, 0x71, 0x74, 0x76, 0x0f, 0x09, 0x09, 0x71, 0x77, 0x74, 0x76, 0x76, 0x74,
                0x0f, 0x71, 0x74, 0x76, 0x77, 0x0c, 0x71, 0x76, 0xe7, 0x02, 0x00,
            ];
            let far_m2 = [
                &small()[..],
                &header(0, 8, None, None),
                &m2_deflate,
                &[0; 8],
            ]
            .concat();
            assert_eq!(refused("far_m2", &far_m2), "corrupt deflate stream");
            // Three empty fixed blocks, then a dynamic block whose table has the letters A and C
            // and no end-of-block code.
            let noeob = [
                &[
                    0x02, 0x08, 0x20, 0x40, 0x01, 0x78, 0x00, 0x24, 0x49, 0x92, 0x24, 0x49, 0x92,
                ][..],
                &[0; 32],
                &[0x80, 0x80],
                &[0; 95],
            ]
            .concat();
            let after = |tail: &[u8]| [&small()[..], tail].concat();
            assert_eq!(
                refused("tail_noeob_aligned", &after(&noeob)),
                "corrupt deflate stream"
            );
            ok("tail_hclen_empty", &after(&[0x04, 0x00, 0x00, 0x00]), SMALL);
            ok("tail_sym16_cut", &after(&[0x04, 0x00, 0x24, 0x24]), SMALL);
            assert_eq!(
                refused("tail_sym16_full", &after(&[0x04, 0x00, 0x24, 0x24, 0x00])),
                "corrupt deflate stream"
            );
        }
    }
}

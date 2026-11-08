//! Compressed FASTA inputs: a file whose name ends in `.gz`, `.bz2` or `.zst` (in any case) is
//! read through the matching decoder, the rule by which NCBI's makeblastdb reads its `-in` files
//! (line GZ, DECISIONS 1, 2, 4, 8 and 9). The decoders are pure Rust (flate2's miniz_oxide
//! backend for deflate, bzip2's libbz2-rs-sys backend, ruzstd), as the pure-Rust runtime boundary
//! requires.
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
    pub fn decoder<'a>(self, input: impl Read + 'a) -> Box<dyn Read + 'a> {
        let input = io::BufReader::with_capacity(
            1 << 16,
            NonEmptySource {
                inner: input,
                started: false,
            },
        );
        match self {
            Compression::Gzip => Box::new(GzipMembers::new(input)),
            Compression::Bzip2 => Box::new(bzip2::read::MultiBzDecoder::new(input)),
            Compression::Zstd => Box::new(ZstdFrames::new(input)),
        }
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

/// A source whose first read must bring a byte: an empty file is "cut short" for every format
/// (no gzip or bzip2 header, no zstd frame), where the decoders would otherwise read it as an
/// empty input.
struct NonEmptySource<R: Read> {
    inner: R,
    started: bool,
}

impl<R: Read> Read for NonEmptySource<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let count = self.inner.read(buf)?;
        if count == 0 && !self.started && !buf.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "the file is empty",
            ));
        }
        self.started = true;
        Ok(count)
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

/// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:135-145
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
/// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:148-209
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
/// `CZipDecompressor` with the flags of `CDecompressIStream::eGZipFile` (`fGZip`:
/// `fCheckFileHeader`, `fAllowConcatenatedGZip`; not `fAllowTransparentRead` nor
/// `fAllowEmptyData`), as `CCompressionStreambuf` drives it. What makeblastdb reads, this reads
/// (the same bytes); where makeblastdb's stream fails, the read fails, with LOSAT's message.
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
/// The deflate data is inflated by flate2 (miniz_oxide) where makeblastdb calls zlib's
/// `inflate()`: the same bytes, the same consumed input at a stream's end. Where the bytes are
/// not deflate data at all (garbage after the last member), miniz_oxide may find the fault at
/// another byte than zlib or accept a damaged Huffman table that zlib refuses; both then fail or
/// accept as their inflaters do.
///
/// NCBI reference (598d8ae6): c++/src/util/compress/api/stream_util.cpp:47, 166-177
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
/// c++/include/util/compress/zlib.hpp:103-104, 117-120
/// ```c++
///         /// Check (and skip) gzip file header on decompression stage
///         fCheckFileHeader       = (1<<2),
/// ...
///         fAllowConcatenatedGZip = (1<<4),
///         /// Set of flags for gzip file support. See each flag description above.
///         fGZip = fCheckFileHeader | fWriteGZipFormat | fAllowConcatenatedGZip,
/// ```
struct GzipMembers<R: BufRead> {
    source: R,
    /// The bytes read for the last header check (`m_Cache`), inflated from `cache_at` before the
    /// source is read again.
    cache: Vec<u8>,
    cache_at: usize,
    inflater: flate2::Decompress,
    state: GzipState,
    /// A gzip header was found (`m_IsGZ`).
    is_gz: bool,
    /// The rest of the file was shorter than `MAX_HEADER_SIZE` at the last header check.
    rest_in_cache: bool,
}

impl<R: BufRead> GzipMembers<R> {
    fn new(source: R) -> Self {
        Self {
            source,
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
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1412-1468
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
    /// of the file (`rest_in_cache`). The window bits are `MAX_WBITS` (15), flate2's.
    fn check_header(&mut self) -> io::Result<()> {
        self.cache.drain(..self.cache_at);
        self.cache_at = 0;
        while self.cache.len() < MAX_HEADER_SIZE {
            let chunk = self.source.fill_buf()?;
            if chunk.is_empty() {
                break;
            }
            let count = chunk.len().min(MAX_HEADER_SIZE - self.cache.len());
            self.cache.extend_from_slice(&chunk[..count]);
            self.source.consume(count);
        }
        if self.cache.is_empty() {
            // The end of the file after a member's trailer (an empty file is refused by
            // `NonEmptySource` before).
            self.state = GzipState::End;
            return Ok(());
        }
        self.rest_in_cache = self.cache.len() < MAX_HEADER_SIZE;
        let header_len = check_gzip_header(&self.cache);
        if header_len > 0 {
            self.cache_at = header_len;
            self.is_gz = true;
        }
        self.inflater.reset(!self.is_gz);
        self.state = GzipState::Inflate;
        Ok(())
    }

    /// The deflate data ended (`Z_STREAM_END`): a member's trailer follows, then another header
    /// check; a zlib stream ends the reading.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1513-1526, 1555-1563
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

    /// Skips up to `left` bytes of a member's trailer, from the cache, then from the source.
    fn skip_trailer(&mut self, mut left: usize) -> io::Result<()> {
        let cached = (self.cache.len() - self.cache_at).min(left);
        self.cache_at += cached;
        left -= cached;
        while left > 0 {
            let available = self.source.fill_buf()?.len();
            if available == 0 {
                break;
            }
            let count = available.min(left);
            self.source.consume(count);
            left -= count;
        }
        self.state = GzipState::CheckHeader;
        Ok(())
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
    fn end_of_input(&mut self) -> io::Result<usize> {
        if self.rest_in_cache && self.inflater.total_in() > 0 && self.inflater.total_out() == 0 {
            self.state = GzipState::End;
            return Ok(0);
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
}

impl<R: BufRead> Read for GzipMembers<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            match self.state {
                GzipState::End => return Ok(0),
                GzipState::CheckHeader => self.check_header()?,
                GzipState::Trailer(left) => self.skip_trailer(left)?,
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
                    // error (`eStatus_Error`, zlib.cpp:1565-1566).
                    let from_cache = self.cache_at < self.cache.len();
                    let input: &[u8] = if from_cache {
                        &self.cache[self.cache_at..]
                    } else {
                        self.source.fill_buf()?
                    };
                    let at_end = input.is_empty();
                    let (read_before, given_before) =
                        (self.inflater.total_in(), self.inflater.total_out());
                    let status =
                        self.inflater
                            .decompress(input, buf, flate2::FlushDecompress::None);
                    let consumed = (self.inflater.total_in() - read_before) as usize;
                    let produced = (self.inflater.total_out() - given_before) as usize;
                    if from_cache {
                        self.cache_at += consumed;
                    } else {
                        self.source.consume(consumed);
                    }
                    match status {
                        Err(_) => return Err(self.fault("corrupt")),
                        Ok(flate2::Status::StreamEnd) => self.stream_end(),
                        Ok(_) if consumed == 0 && produced == 0 => {
                            if at_end {
                                return self.end_of_input();
                            }
                            // Input that the inflater neither reads nor answers with an error.
                            return Err(self.fault("corrupt"));
                        }
                        Ok(_) => {}
                    }
                    if produced > 0 {
                        return Ok(produced);
                    }
                }
            }
        }
    }
}

/// The frames of a zstd file, read one after the other. ruzstd's `StreamingDecoder` reads one
/// frame and stops; pzstd and `cat a.zst b.zst` write several, pzstd with a skippable frame
/// before each, which this reader passes over (`ReadFrameHeaderError::SkipFrame`). At the end of
/// each frame the bytes are checked against the frame's `Content_Checksum` and
/// `Frame_Content_Size` (`check_frame`), which ruzstd reads but does not compare.
pub struct ZstdFrames<R: BufRead> {
    source: R,
    decoder: ruzstd::decoding::FrameDecoder,
    in_frame: bool,
    /// Bytes of the current frame given out so far.
    produced: u64,
}

impl<R: BufRead> ZstdFrames<R> {
    pub fn new(source: R) -> Self {
        Self {
            source,
            decoder: ruzstd::decoding::FrameDecoder::new(),
            in_frame: false,
            produced: 0,
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

    /// Starts the next frame, skipping skippable frames; `false` at the end of the source.
    fn next_frame(&mut self) -> io::Result<bool> {
        loop {
            if self.source.fill_buf()?.is_empty() {
                return Ok(false);
            }
            match self.decoder.reset(&mut self.source) {
                Ok(()) => return Ok(true),
                Err(ruzstd::decoding::errors::FrameDecoderError::ReadFrameHeaderError(
                    ruzstd::decoding::errors::ReadFrameHeaderError::SkipFrame { length, .. },
                )) => {
                    let skipped = io::copy(
                        &mut (&mut self.source).take(u64::from(length)),
                        &mut io::sink(),
                    )?;
                    if skipped != u64::from(length) {
                        return Err(io::Error::new(
                            io::ErrorKind::UnexpectedEof,
                            "a skippable frame is cut short",
                        ));
                    }
                }
                Err(error) => return Err(invalid_zstd(error)),
            }
        }
    }
}

fn invalid_zstd(error: ruzstd::decoding::errors::FrameDecoderError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

impl<R: BufRead> Read for ZstdFrames<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            if !self.in_frame {
                if !self.next_frame()? {
                    return Ok(0);
                }
                self.in_frame = true;
                self.produced = 0;
            }
            if self.decoder.can_collect() > 0 {
                let count = self.decoder.read(buf)?;
                self.produced += count as u64;
                return Ok(count);
            }
            if self.decoder.is_finished() {
                self.check_frame()?;
                self.in_frame = false;
                continue;
            }
            self.decoder
                .decode_blocks(
                    &mut self.source,
                    ruzstd::decoding::BlockDecodingStrategy::UptoBlocks(1),
                )
                .map_err(invalid_zstd)?;
        }
    }
}

#[cfg(test)]
pub(crate) mod test_data {
    //! Compressed test inputs made without a compressor binary: flate2 and bzip2 encode; zstd,
    //! which LOSAT only decodes, is written as frames of raw (stored) blocks, which every zstd
    //! decoder reads.
    use std::io::Write;

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
mod tests {
    use super::test_data::{bzip2_streams, gzip_members, zstd_frames};
    use super::{compression_of_name, Compression};
    use std::io::Read;
    use std::path::Path;

    fn inflate(compression: Compression, bytes: &[u8]) -> std::io::Result<Vec<u8>> {
        let mut out = Vec::new();
        compression.decoder(bytes).read_to_end(&mut out)?;
        Ok(out)
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
        /// at a time, which must agree: the result does not hang on the sizes of the reads.
        fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, String> {
            let mut whole = Vec::new();
            let whole = Compression::Gzip
                .decoder(bytes)
                .read_to_end(&mut whole)
                .map(|_| whole)
                .map_err(|error| error.to_string());
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

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:148-209
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

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:148-209: the header
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
    }
}

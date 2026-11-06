//! Compressed FASTA inputs: a file whose name ends in `.gz`, `.bz2` or `.zst` (in any case) is
//! read through the matching decoder, the rule by which NCBI's makeblastdb reads its `-in` files
//! (line GZ, DECISIONS 1, 2, 4, 8 and 9). The decoders are pure Rust (flate2's miniz_oxide
//! backend, bzip2's libbz2-rs-sys backend, ruzstd), as the pure-Rust runtime boundary requires.
//! Every gzip member, bzip2 stream and zstd frame of a file is read: bgzip, pbzip2 and pzstd write
//! several. NCBI's makeblastdb reads every gzip member too (`fGZip` holds
//! `fAllowConcatenatedGZip`) but only the first bzip2 stream and the first zstd frame, and drops
//! the rest without a message (`PD-LOSAT-NCBI-DEFECTS`, makeblastdb exception).
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050
//! ```c++
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
    /// error of the read that meets it; the first read already fails on a file that does not
    /// start with the format's magic bytes.
    pub fn decoder<'a>(self, input: impl Read + 'a) -> Box<dyn Read + 'a> {
        let input = io::BufReader::with_capacity(
            1 << 16,
            NonEmptySource {
                inner: input,
                started: false,
            },
        );
        match self {
            Compression::Gzip => Box::new(flate2::read::MultiGzDecoder::new(input)),
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

    // No NCBI counterpart (line GZ, DECISIONS 5): text under a compressed name, data cut short,
    // or an empty file is an error of the read.
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
}

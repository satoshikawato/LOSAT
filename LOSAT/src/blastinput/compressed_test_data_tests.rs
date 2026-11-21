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
        let mut encoder = bzip2::write::BzEncoder::new(&mut bytes, bzip2::Compression::default());
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

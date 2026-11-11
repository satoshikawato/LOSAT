//! The pushed inflater (`PushInflater`, LOSAT Web's) against the command line's reader
//! (`Compression::decoder`) on zstd and bzip2 files: whole, damaged, cut at every byte, with
//! skippable frames, several frames or streams, and bytes after the last one, each pushed in
//! pieces of every size of `test_data::PIECES`. Both run the same `Inflater`; they must give the
//! same bytes, or an error of the same kind and words. The `.gz` cases are those of
//! `tests::gzip` (makeblastdb's edge cases), whose helper makes the same check.
//!
//! No NCBI counterpart: makeblastdb reads only the first bzip2 stream and the first zstd frame
//! (`PD-LOSAT-NCBI-DEFECTS`, makeblastdb exception), and the words of a fault are LOSAT's (line
//! GZ, DECISIONS 5).

use super::test_data::{assert_pushed_as_read, bzip2_streams, pushed, read, zstd_frames, PIECES};
use super::{block_len, frame_header_len, Compression, PushInflater, Whole, ZstdFrames};
use std::io::Write;

fn fixture(path: &str) -> Vec<u8> {
    std::fs::read(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join(path),
    )
    .unwrap()
}

/// The reader's outcome for `bytes`, after the check that every piece size pushes the same.
fn checked(compression: Compression, bytes: &[u8]) -> super::test_data::Outcome {
    let outcome = read(compression, bytes);
    let result = outcome
        .clone()
        .map_err(|(kind, words)| std::io::Error::new(kind, words));
    assert_pushed_as_read(compression, bytes, &result);
    outcome
}

/// The check of `checked`, where the outcome itself is not asserted.
fn alike(compression: Compression, bytes: &[u8]) {
    let _outcome = checked(compression, bytes);
}

fn words(outcome: super::test_data::Outcome) -> String {
    match outcome {
        Ok(bytes) => panic!("inflated {} bytes", bytes.len()),
        Err((_, words)) => words,
    }
}

/// Text in 60-letter lines, `records` records of `letters` letters, from a fixed generator.
fn fasta(records: usize, letters: usize) -> Vec<u8> {
    let mut state = 7u64;
    let mut text = Vec::new();
    for record in 1..=records {
        text.extend_from_slice(format!(">r{record}\n").as_bytes());
        for i in 0..letters {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            text.push(b"ACGT"[(state >> 62) as usize]);
            if i % 60 == 59 || i + 1 == letters {
                text.push(b'\n');
            }
        }
    }
    text
}

/// Every prefix of `bytes` (where `every` is 1), or every `every`-th one and the last 16.
fn cuts(bytes: &[u8], every: usize) -> impl Iterator<Item = &[u8]> {
    (0..bytes.len())
        .filter(move |&len| len % every == 0 || len + 16 >= bytes.len())
        .map(move |len| &bytes[..len])
}

// The zstd CLI's files (compressed blocks, a Content_Checksum, with and without a content
// size), frames of raw blocks after skippable frames (`zstd_frames`), and frames written here: an
// RLE block, a block of the reserved type, a block larger than 128 KiB.
#[test]
fn zstd_pieces_inflate_as_the_reader_reads() {
    let nucl = fixture("compressed/mask_nucl.fna.zst");
    let prot = fixture("compressed/prot.faa.zst");
    let plain_nucl = fixture("seqdb/mask_nucl.fna");
    assert_eq!(checked(Compression::Zstd, &nucl), Ok(plain_nucl.clone()));
    let parts: [&[u8]; 3] = [b">q1\nACGT\n", b"", b">q2\nTTTT\n"];
    assert_eq!(
        checked(Compression::Zstd, &zstd_frames(&parts)),
        Ok(b">q1\nACGT\n>q2\nTTTT\n".to_vec())
    );
    let long = fasta(5, 60_000);
    assert_eq!(
        checked(Compression::Zstd, &zstd_frames(&[&long, &long])),
        Ok([&long[..], &long].concat())
    );
    // Output past the push buffer: the sink gets it in full buffers.
    let big = vec![b'A'; 3 * PushInflater::OUTPUT_SIZE + 5];
    assert_eq!(checked(Compression::Zstd, &zstd_frames(&[&big])), Ok(big));
    let joined = [&zstd_frames(&parts[..1])[..], &nucl, &prot, &nucl].concat();
    assert!(checked(Compression::Zstd, &joined).is_ok());
    // Every cut of a frame (in its header, its blocks, its checksum) and of skippable frames.
    for bytes in [&nucl, &zstd_frames(&parts), &joined] {
        for cut in cuts(bytes, 1) {
            alike(Compression::Zstd, cut);
        }
    }
    // The checksum, then the content size, wrong.
    let mut badsum = nucl.clone();
    *badsum.last_mut().unwrap() ^= 1;
    assert_eq!(
        words(checked(Compression::Zstd, &badsum)),
        "a zstd frame's checksum does not match its data"
    );
    let mut badsize = zstd_frames(&[b">q1\nACGT\n"]);
    badsize[17..21].copy_from_slice(&10u32.to_le_bytes());
    assert_eq!(
        words(checked(Compression::Zstd, &badsize)),
        "a zstd frame holds 9 bytes where its header gives 10"
    );
    // Skippable frames alone, cut in their length, cut in their bytes.
    let skippable = |length: u32, bytes: &[u8]| {
        [
            &0x184D_2A50u32.to_le_bytes()[..],
            &length.to_le_bytes(),
            bytes,
        ]
        .concat()
    };
    assert_eq!(
        checked(Compression::Zstd, &skippable(0, b"")),
        Ok(Vec::new())
    );
    assert_eq!(
        checked(Compression::Zstd, &skippable(3, b"abc")),
        Ok(Vec::new())
    );
    alike(Compression::Zstd, &skippable(3, b"abc")[..6]);
    assert_eq!(
        words(checked(Compression::Zstd, &skippable(10, b"abc"))),
        "a skippable frame is cut short"
    );
    // Bytes after the last frame, text, an empty file.
    for tail in [&b"junk"[..], b"ju", &0xFD2F_B528u32.to_le_bytes(), b"\n"] {
        assert!(checked(Compression::Zstd, &[&nucl[..], tail].concat()).is_err());
    }
    assert!(checked(Compression::Zstd, b">q1\nACGT\n").is_err());
    assert_eq!(words(checked(Compression::Zstd, b"")), "the file is empty");
    // Frames written here: Single_Segment, a 1-byte content size, then one block.
    let frame =
        |size: u8, block: &[u8]| [&0xFD2F_B528u32.to_le_bytes()[..], &[0x20, size], block].concat();
    let block_header = |kind: u32, size: u32| (1 | (kind << 1) | (size << 3)).to_le_bytes();
    let rle = frame(10, &[&block_header(1, 10)[..3], b"A"].concat());
    assert_eq!(checked(Compression::Zstd, &rle), Ok(vec![b'A'; 10]));
    for cut in cuts(&rle, 1) {
        alike(Compression::Zstd, cut);
    }
    assert!(checked(Compression::Zstd, &frame(5, &block_header(3, 5)[..3])).is_err());
    let oversize = frame(5, &[&block_header(0, 200_000)[..3], b"ACGTA"].concat());
    assert!(checked(Compression::Zstd, &oversize).is_err());
}

// A zstd unit (a frame header, a block) that one piece of input holds whole is decoded in that
// piece; one that crosses pieces is gathered first (`ZstdFrames::gather`). Pushed in two pieces
// cut inside and between units, the frames give the reader's bytes, and so does the reader on
// frames larger than its 256 KiB read blocks.
#[test]
fn a_zstd_unit_is_decoded_in_its_piece_or_gathered_across_pieces() {
    let text: &[u8] = b">q1\nACGTACGTAC\n";
    // Each frame: a skippable frame (12 bytes), the header (9), one raw block (3 + 15).
    let frames = zstd_frames(&[text, text]);
    let (header_at, block_at) = (12, 21);
    let header = &frames[header_at..block_at];
    let block = &frames[block_at..block_at + 3 + text.len()];
    let mut core = ZstdFrames::new();
    let mut taken = 0;
    assert_eq!(
        core.gather(header, &mut taken, false, frame_header_len),
        Some(Whole::Input(9))
    );
    assert_eq!((taken, core.unit.len()), (0, 0));
    assert_eq!(
        core.gather(block, &mut taken, false, |unit| block_len(unit, false)),
        Some(Whole::Input(block.len()))
    );
    let mut taken = 0;
    assert_eq!(
        core.gather(&header[..6], &mut taken, false, frame_header_len),
        None
    );
    assert_eq!(taken, 6);
    let mut taken = 0;
    assert_eq!(
        core.gather(&header[6..], &mut taken, false, frame_header_len),
        Some(Whole::Gathered)
    );
    assert_eq!((taken, core.unit.as_slice()), (3, header));
    let whole = read(Compression::Zstd, &frames);
    assert_eq!(whole, Ok([text, text].concat()));
    for cut in [
        header_at + 4,
        block_at,
        block_at + 2,
        block_at + 10,
        block_at + 18,
        frames.len() - 1,
    ] {
        let mut inflater = PushInflater::new(Compression::Zstd);
        let mut out = Vec::new();
        let mut sink = |inflated: &[u8]| out.extend_from_slice(inflated);
        inflater.push(&frames[..cut], &mut sink).unwrap();
        inflater.push(&frames[cut..], &mut sink).unwrap();
        inflater.finish(&mut sink).unwrap();
        assert_eq!(Ok(out), whole, "cut at {cut}");
    }
    // Raw blocks of 128 KiB: some lie inside the reader's 256 KiB read blocks, some cross them.
    let big = fasta(8, 120_000);
    assert_eq!(checked(Compression::Zstd, &zstd_frames(&[&big])), Ok(big));
}

// Every single-bit flip of the zstd CLI's files: the pieces give what the reader gives (the
// intact bytes, or the same error).
#[test]
fn damaged_zstd_pieces_fail_as_the_reader_fails() {
    for path in ["compressed/mask_nucl.fna.zst", "compressed/prot.faa.zst"] {
        let bytes = fixture(path);
        for at in 0..bytes.len() {
            for bit in 0..8 {
                let mut damaged = bytes.clone();
                damaged[at] ^= 1 << bit;
                alike(Compression::Zstd, &damaged);
            }
        }
    }
}

/// Bzip2 streams of the given level, concatenated.
fn bzip2_level(streams: &[&[u8]], level: u32) -> Vec<u8> {
    let mut bytes = Vec::new();
    for stream in streams {
        let mut encoder = bzip2::write::BzEncoder::new(&mut bytes, bzip2::Compression::new(level));
        encoder.write_all(stream).unwrap();
        encoder.finish().unwrap();
    }
    bytes
}

// Streams one after the other (pbzip2, `cat`), an empty stream, a stream of three blocks (level
// 1: 100 000-byte blocks) whole, cut and damaged, bytes after the last stream, text, an empty
// file.
#[test]
fn bzip2_pieces_inflate_as_the_reader_reads() {
    let parts: [&[u8]; 3] = [b">q1\nACGT\n", b"", b">q2\nTTTT\n"];
    let small = bzip2_streams(&parts);
    assert_eq!(
        checked(Compression::Bzip2, &small),
        Ok(b">q1\nACGT\n>q2\nTTTT\n".to_vec())
    );
    for cut in cuts(&small, 1) {
        alike(Compression::Bzip2, cut);
    }
    let text = fasta(4, 60_000);
    let blocks = bzip2_level(&[&text, b">q3\nGG\n"], 1);
    assert_eq!(
        checked(Compression::Bzip2, &blocks),
        Ok([&text[..], b">q3\nGG\n"].concat())
    );
    for cut in cuts(&blocks, 4999) {
        alike(Compression::Bzip2, cut);
    }
    for at in (0..blocks.len()).step_by(3001) {
        let mut damaged = blocks.clone();
        damaged[at] ^= 0x10;
        alike(Compression::Bzip2, &damaged);
    }
    let big = vec![b'C'; 3 * PushInflater::OUTPUT_SIZE + 5];
    assert_eq!(
        checked(Compression::Bzip2, &bzip2_streams(&[&big])),
        Ok(big)
    );
    for tail in [&b"junk"[..], b"B", b"BZh9", b"\n"] {
        assert!(checked(Compression::Bzip2, &[&small[..], tail].concat()).is_err());
    }
    assert!(checked(Compression::Bzip2, b">q1\nACGT\n").is_err());
    assert_eq!(words(checked(Compression::Bzip2, b"")), "the file is empty");
}

// Output past the push buffer for gzip (the other `.gz` cases are those of `tests::gzip`), an
// empty file pushed as empty pieces, and the error that a pushed inflater gives again after it
// failed.
#[test]
fn pushed_pieces_end_and_fail_as_the_reader_does() {
    let big = fasta(2, 2_000_000);
    let gz = super::test_data::gzip_members(&[&big, b">q3\nGG\n"]);
    assert_eq!(
        checked(Compression::Gzip, &gz),
        Ok([&big[..], b">q3\nGG\n"].concat())
    );
    for compression in Compression::ALL {
        let mut inflater = PushInflater::new(compression);
        inflater.push(b"", &mut |_| panic!("no bytes")).unwrap();
        let error = inflater.finish(&mut |_| panic!("no bytes")).unwrap_err();
        assert_eq!(
            (error.kind(), error.to_string()),
            (
                std::io::ErrorKind::UnexpectedEof,
                "the file is empty".to_string()
            ),
            "{compression:?}"
        );
        // Text of more than the 4 KiB of a gzip header check, so that every format fails on
        // the push.
        let text = fasta(1, 5000);
        let mut inflater = PushInflater::new(compression);
        let first = inflater.push(&text, &mut |_| {}).unwrap_err();
        let again = inflater.push(b"more", &mut |_| {}).unwrap_err();
        assert_eq!(
            (again.kind(), again.to_string()),
            (first.kind(), first.to_string()),
            "{compression:?}"
        );
        let last = inflater.finish(&mut |_| {}).unwrap_err();
        assert_eq!(last.to_string(), first.to_string(), "{compression:?}");
        assert_eq!(
            Err((first.kind(), first.to_string())),
            read(compression, &text),
            "{compression:?}"
        );
    }
    // The piece sizes of the tests reach past the reader's 64 KiB blocks.
    assert_eq!(PIECES[5], 1 << 16);
    assert_eq!(
        pushed(Compression::Gzip, &gz, 1 << 16),
        read(Compression::Gzip, &gz)
    );
}

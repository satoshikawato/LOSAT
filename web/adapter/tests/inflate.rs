//! The inflation exports of ABI v2 (`losat_web2_inflate_begin`, `_chunk`, `_end`;
//! docs/web/abi_v2.md §4, §5, §9) on gzip, bzip2 and zstd files, whole, damaged, cut and empty,
//! pushed in pieces of one byte and whole: the bytes of stream 5 are those of the command line's
//! reader (`Compression::decoder`), and an error is the command line's sentence without its
//! `-query: ` (the sentences below are the stderr of `losat blastn -query <file>` at the commit of
//! the push core, swebok-improve SW8 D2 part A, after its `Error: -query: `).
//!
//! No NCBI counterpart: NCBI's search programs do not inflate their inputs, and the words of a
//! fault are LOSAT's (line GZ, DECISIONS 5).

use std::io::{Read, Write};
use std::sync::Mutex;

use losat_web_adapter::emit::{CHUNK, EMITTED, STREAM_INFLATED, STREAM_RESPONSE};
use losat_web_adapter::{
    losat_web2_inflate_begin, losat_web2_inflate_chunk, losat_web2_inflate_end,
    losat_web2_last_error_len, losat_web2_last_error_ptr,
};
use LOSAT::web_adapter::Compression;

/// The module has one last error: the tests take turns to call the exports and read it.
static CALLS: Mutex<()> = Mutex::new(());

fn turn() -> std::sync::MutexGuard<'static, ()> {
    CALLS
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

fn last_error() -> String {
    let (ptr, len) = (losat_web2_last_error_ptr(), losat_web2_last_error_len());
    // SAFETY: the last error stays in place until the next export call.
    String::from_utf8(unsafe { std::slice::from_raw_parts(ptr, len) }.to_vec()).unwrap()
}

fn taken() -> Vec<(u32, Vec<u8>)> {
    EMITTED.with(|emitted| std::mem::take(&mut *emitted.borrow_mut()))
}

fn begin(name: &str) -> i32 {
    // SAFETY: the name is a live buffer of its length.
    unsafe { losat_web2_inflate_begin(name.as_ptr(), name.len()) }
}

fn chunk(handle: i32, bytes: &[u8]) -> i32 {
    // SAFETY: the bytes are a live buffer of their length.
    unsafe { losat_web2_inflate_chunk(handle as u32, bytes.as_ptr(), bytes.len()) }
}

/// What the exports give for the file `name` pushed in pieces of `piece` bytes: the inflated
/// bytes of stream 5, or the last error. Checks the response of `begin`, that nothing but stream
/// 5 is emitted after it, in chunks of 1 to `CHUNK` bytes, and that the handle is gone after
/// the error or the end.
fn inflated(
    compression: Compression,
    name: &str,
    plain: &str,
    bytes: &[u8],
    piece: usize,
) -> Result<Vec<u8>, String> {
    taken();
    let handle = begin(name);
    assert!(handle >= 1, "{name}: begin gave {handle}: {}", last_error());
    let response = format!(
        "{{\"compression\":\"{}\",\"plain_name\":\"{plain}\"}}",
        compression.format()
    );
    assert_eq!(
        taken(),
        [(STREAM_RESPONSE, response.into_bytes())],
        "{name}"
    );
    let mut failed = bytes.chunks(piece).any(|piece| chunk(handle, piece) != 0);
    if !failed {
        failed = losat_web2_inflate_end(handle as u32) != 0;
    }
    let result = if failed { Err(last_error()) } else { Ok(()) };
    // The handle is dropped after an error and after the end.
    assert_eq!(losat_web2_inflate_end(handle as u32), -1, "{name}");
    assert_eq!(
        last_error(),
        format!("unknown inflation {handle}"),
        "{name}"
    );
    let mut out = Vec::new();
    for (stream, chunk) in taken() {
        assert_eq!(stream, STREAM_INFLATED, "{name}");
        assert!(
            !chunk.is_empty() && chunk.len() <= CHUNK,
            "{name}: {}",
            chunk.len()
        );
        out.extend_from_slice(&chunk);
    }
    result.map(|()| out)
}

/// What the command line's reader gives for the file `name`: the inflated bytes, or the
/// sentence of `Compression::not_data` with the reader's error.
fn read(compression: Compression, name: &str, bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    match compression.decoder(bytes).read_to_end(&mut out) {
        Ok(_) => Ok(out),
        Err(error) => Err(compression.not_data(name, &error.to_string())),
    }
}

/// Checks that the file inflates through the exports as the reader reads it, pushed whole and in
/// pieces of one byte (of 4097 bytes for a file over 64 KiB, which would take a million calls of
/// the debug build), and returns that outcome. The engine's own tests push every input of its
/// decoders in pieces of 1, 3, 4095, 4096, 4097 and 65536 bytes (`compressed::test_data`).
fn alike(
    compression: Compression,
    name: &str,
    plain: &str,
    bytes: &[u8],
) -> Result<Vec<u8>, String> {
    let read = read(compression, name, bytes);
    let small = if bytes.len() <= 1 << 16 { 1 } else { 4097 };
    for piece in [small, usize::MAX] {
        let pushed = inflated(compression, name, plain, bytes, piece);
        assert!(
            pushed == read,
            "{name} in pieces of {piece}: pushed {:?}, read {:?}",
            pushed.as_ref().map(Vec::len),
            read.as_ref().map(Vec::len)
        );
    }
    read
}

/// Text in 60-letter lines, `records` records of `letters` letters, from a fixed generator.
fn fasta(records: usize, letters: usize) -> Vec<u8> {
    let mut state = 11u64;
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

fn gzip(text: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(text).unwrap();
    encoder.finish().unwrap()
}

fn bzip2(text: &[u8]) -> Vec<u8> {
    let mut encoder = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
    encoder.write_all(text).unwrap();
    encoder.finish().unwrap()
}

/// One zstd frame of raw blocks (Single_Segment, a 4-byte Frame_Content_Size, no checksum),
/// which every zstd decoder reads; LOSAT only decodes zstd.
fn zstd_raw(text: &[u8]) -> Vec<u8> {
    let mut bytes = 0xFD2F_B528u32.to_le_bytes().to_vec();
    bytes.push(0x80 | 0x20);
    bytes.extend_from_slice(&(text.len() as u32).to_le_bytes());
    let mut blocks = text.chunks(128 * 1024).peekable();
    while let Some(block) = blocks.next() {
        let header = u32::from(blocks.peek().is_none()) | ((block.len() as u32) << 3);
        bytes.extend_from_slice(&header.to_le_bytes()[..3]);
        bytes.extend_from_slice(block);
    }
    bytes
}

/// The zstd CLI's frame of `LOSAT/tests/compressed/mask_nucl.fna.zst` (entropy-coded blocks, a
/// content size and an XXH64 checksum).
fn zstd_cli_frame() -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../LOSAT/tests/compressed/mask_nucl.fna.zst");
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

#[test]
fn a_plain_name_has_no_inflation() {
    let _turn = turn();
    taken();
    for name in ["a.fa", "a.gzip", "a.zstd", "a.bz", "gz", "-", ""] {
        assert_eq!(begin(name), 0, "{name}");
        assert_eq!(last_error(), "", "{name}");
    }
    assert_eq!(taken(), []);
}

#[test]
fn intact_files_inflate_as_the_reader_reads_them() {
    let _turn = turn();
    let small = fasta(3, 100);
    // More than one 1 MiB chunk of stream 5 and of the inflater's output buffer.
    let big = fasta(40, 60_000);
    assert!(big.len() > 2 * CHUNK);
    for (compression, suffix, pack) in [
        (Compression::Gzip, "gz", gzip as fn(&[u8]) -> Vec<u8>),
        (Compression::Bzip2, "bz2", bzip2),
        (Compression::Zstd, "zst", zstd_raw),
    ] {
        for (text, name) in [(&small, "small"), (&big, "big")] {
            let file = format!("{name}.fa.{suffix}");
            let plain = format!("{name}.fa");
            assert_eq!(
                alike(compression, &file, &plain, &pack(text)).as_ref(),
                Ok(text)
            );
        }
        // Members, streams or frames one after the other (bgzip, pbzip2, `cat`), and the suffix
        // in another case.
        let joined = [pack(&small), pack(&big)].concat();
        let file = format!("x.FA.{}", suffix.to_ascii_uppercase());
        assert_eq!(
            alike(compression, &file, "x.FA", &joined),
            Ok([&small[..], &big[..]].concat())
        );
        // Nothing left of the name without its suffix.
        assert_eq!(
            alike(compression, &format!(".{suffix}"), "input", &pack(&small)).as_ref(),
            Ok(&small)
        );
    }
    let cli = zstd_cli_frame();
    let text = alike(
        Compression::Zstd,
        "mask_nucl.fna.zst",
        "mask_nucl.fna",
        &cli,
    )
    .unwrap();
    assert!(text.starts_with(b">"));
}

#[test]
fn damaged_files_fail_with_the_command_line_sentence() {
    let _turn = turn();
    let text = fasta(200, 440);
    // The checksum of the zstd CLI's frame, wrong.
    let mut badsum = zstd_cli_frame();
    *badsum.last_mut().unwrap() ^= 1;
    // A gzip member and a bzip2 stream cut in their compressed data.
    let gz = gzip(&text);
    let bz = bzip2(&text);
    let cases: [(Compression, &str, &[u8], &str); 7] = [
        (
            Compression::Zstd,
            "badsum.fa.zst",
            &badsum,
            "a zstd frame's checksum does not match its data",
        ),
        (
            Compression::Gzip,
            "cut.fa.gz",
            &gz[..gz.len() / 2],
            "incomplete deflate stream",
        ),
        (
            Compression::Bzip2,
            "cut.fa.bz2",
            &bz[..bz.len() / 2],
            "decompression not finished but EOF reached",
        ),
        (
            Compression::Gzip,
            "text.fa.gz",
            b">q\nACGT\n",
            "invalid gzip header",
        ),
        (Compression::Gzip, "empty.fa.gz", b"", "the file is empty"),
        (Compression::Bzip2, "empty.fa.bz2", b"", "the file is empty"),
        (Compression::Zstd, "empty.fa.zst", b"", "the file is empty"),
    ];
    for (compression, name, bytes, fault) in cases {
        let (suffix, format) = (compression.suffix(), compression.format());
        let sentence = format!(
            "'{name}' is named {suffix} but is not {format} data ({fault}); LOSAT reads a file named {suffix} through {format}"
        );
        let plain = &name[..name.len() - suffix.len()];
        assert_eq!(
            alike(compression, name, plain, bytes),
            Err(sentence),
            "{name}"
        );
    }
}

#[test]
fn an_unknown_handle_is_an_error() {
    let _turn = turn();
    assert_eq!(chunk(i32::MAX, b"x"), -1);
    assert_eq!(last_error(), format!("unknown inflation {}", i32::MAX));
    assert_eq!(losat_web2_inflate_end(0), -1);
    assert_eq!(last_error(), "unknown inflation 0");
}

//! The lexical readers read the same values, errors, positions and lines from every source: the
//! text held whole (`SliceSource`), a stream read in blocks (`StreamSource`), and a stream whose
//! reads give one byte each, so that every byte is read after the buffer was filled again (a value
//! split across any two reads, a text that ends anywhere). The texts are the ASN.1 text fixtures
//! of the readers that use this module (rpsbproc archives, PSSMs, mask files, `-outfmt 11`
//! archives), whole and cut at many lengths.

use super::*;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};

/// A `Read` that gives one byte per call.
pub(crate) struct OneByte<'a>(pub(crate) &'a [u8]);

impl Read for OneByte<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match (self.0.split_first(), buf.first_mut()) {
            (Some((&byte, rest)), Some(slot)) => {
                *slot = byte;
                self.0 = rest;
                Ok(1)
            }
            _ => Ok(0),
        }
    }
}

/// One step of a walk: what was read, and the reader's position and line after it.
fn note<S: AsnSource>(r: &Reader<S>, log: &mut Vec<String>, what: String) {
    log.push(format!("{what} @{}:{}", r.pos, r.line));
}

fn show(bytes: &[u8]) -> String {
    format!("{:?}", String::from_utf8_lossy(bytes))
}

/// Values nested deeper than this (in blocks or after names) are skipped whole (`skip_value`
/// keeps no stack), so that the walk of the deeply nested fixtures fits a test thread's stack.
const WALK_DEPTH: usize = 48;

/// Reads one value of any type with the lexical readers, choosing among the readers of a kind
/// by `step` so that each is used: blocks element by element (or skipped whole), strings, octet
/// strings, numbers, REALs after the names that hold one, booleans, NULL, and the names of
/// members, variants and enumerated values. `depth` counts the blocks and names around the value.
fn value<S: AsnSource>(
    r: &mut Reader<S>,
    log: &mut Vec<String>,
    name: &[u8],
    step: &mut usize,
    depth: usize,
) -> Result<()> {
    *step += 1;
    let c = r.skip_white_space()?;
    if depth >= WALK_DEPTH {
        r.skip_value()?;
        note(r, log, "skipped deep".into());
        return Ok(());
    }
    match c {
        b'{' if matches!(name, b"real" | b"lambda" | b"k" | b"h") => {
            let v = r.real()?;
            note(r, log, format!("real {:x}", v.to_bits()));
        }
        b'{' if *step % 11 == 0 => {
            r.skip_value()?;
            note(r, log, "skipped".into());
        }
        b'{' => {
            let mut block_start = r.start_block()?;
            while r.next_element(&mut block_start)? {
                value(r, log, b"", step, depth + 1)?;
            }
            r.end_block()?;
            note(r, log, "}".into());
        }
        b'"' => match *step % 4 {
            0 => {
                let text = r.string_as_read()?;
                note(r, log, format!("string {}", show(&text)));
            }
            1 => {
                let text = r.string()?;
                note(r, log, format!("visible {}", show(&text)));
            }
            n => {
                let (mut out, mut bad) = (Vec::new(), Vec::new());
                r.visible_string(&mut out, true, n == 3, &mut |line, byte, piece| {
                    bad.push(format!("{line} {byte} {}", show(piece)))
                })?;
                note(r, log, format!("pieces {} {bad:?}", show(&out)));
            }
        },
        b'\'' => {
            let octets = r.octets()?;
            note(r, log, format!("octets {}", show(&octets)));
        }
        b'0'..=b'9' if *step % 3 == 2 => {
            let n = r.uint4()?;
            note(r, log, format!("uint4 {n}"));
        }
        b'-' | b'+' | b'0'..=b'9' => {
            if *step % 3 == 0 {
                let n = r.int8()?;
                note(r, log, format!("int8 {n}"));
            } else {
                let digits = r.number()?;
                note(r, log, format!("number {}", show(r.bytes(digits))));
            }
        }
        b'T' | b'F' => {
            let b = r.boolean()?;
            note(r, log, format!("bool {b}"));
        }
        b'N' => {
            r.null()?;
            note(r, log, "null".into());
        }
        b'[' | b'a'..=b'z' => {
            let id = match r.member_or_tag()? {
                MemberId::Name(id) => r.bytes(id).to_vec(),
                MemberId::Tag(id) => [b"tag ", r.bytes(id)].concat(),
            };
            note(r, log, format!("id {}", show(&id)));
            if !matches!(r.skip_white_space()?, b',' | b'}') {
                value(r, log, &id, step, depth + 1)?;
            }
        }
        _ if matches!(name, b"real" | b"lambda" | b"k" | b"h") => {
            let v = r.real_word()?;
            note(r, log, format!("real word {:?}", v.map(f64::to_bits)));
        }
        _ => {
            r.skip_value()?;
            note(r, log, "skipped".into());
        }
    }
    Ok(())
}

/// The objects of a text, each a type name, `::=` and a value, up to the error that ends it.
fn objects<S: AsnSource>(r: &mut Reader<S>, log: &mut Vec<String>) -> Result<()> {
    let mut step = 0usize;
    loop {
        let id = r.file_header()?;
        let name = r.bytes(id).to_vec();
        note(r, log, format!("header {}", show(&name)));
        value(r, log, &name, &mut step, 0)?;
    }
}

/// The walk of `src` with the string rule `mode` (0: any byte; 1: bytes outside the visible
/// range refused; 2: refused, deferred): every step, and the error and refusal it ended with.
fn walk<S: AsnSource>(src: S, mode: usize) -> Vec<String> {
    let mut r = Reader::new(src, "text", "test reader");
    if mode >= 1 {
        r = r.refuse_invisible_bytes();
    }
    if mode == 2 {
        r = r.defer_refusals();
    }
    let mut log = Vec::new();
    let error = objects(&mut r, &mut log).expect_err("a walk ends with an error");
    note(
        &r,
        &mut log,
        format!("error {error:?} refusal {:?}", r.refusal),
    );
    log
}

/// The walks of `text` from the three sources, which must be equal; returns a digest of the walk.
fn walk_every_source(label: &str, text: &[u8], mode: usize) -> u64 {
    let slice = walk(SliceSource::new(text), mode);
    let first = text.len().min(1);
    let blocks = walk(StreamSource::new(&text[first..], &text[..first]), mode);
    assert_eq!(slice, blocks, "{label}: slice and stream differ");
    let bytes = walk(StreamSource::new(OneByte(text), &[]), mode);
    assert_eq!(slice, bytes, "{label}: slice and one byte per read differ");
    let mut hasher = DefaultHasher::new();
    slice.hash(&mut hasher);
    hasher.finish()
}

/// The ASN.1 text fixtures of the readers of this module (binary and other files among them are
/// read too: their walks end at their first byte or two).
fn fixture_texts() -> Vec<(String, Vec<u8>)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let dirs = [
        "fixtures/rpsbproc/arch",
        "fixtures/rpsbproc/arch/bad",
        "fixtures/pssm",
        "seqdb/mask_parse",
        "fixtures/outfmt11",
    ];
    let mut texts = Vec::new();
    for dir in dirs {
        let mut paths: Vec<PathBuf> = std::fs::read_dir(root.join(dir))
            .unwrap_or_else(|e| panic!("{dir}: {e}"))
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.is_file())
            .collect();
        paths.sort();
        for path in paths {
            let data = std::fs::read(&path).unwrap();
            let data = if path.extension().is_some_and(|e| e == "gz") {
                let mut out = Vec::new();
                flate2::read::MultiGzDecoder::new(&data[..])
                    .read_to_end(&mut out)
                    .unwrap();
                out
            } else {
                data
            };
            let label = format!("{dir}/{}", path.file_name().unwrap().to_string_lossy());
            texts.push((label, data));
        }
    }
    texts
}

/// Every fixture whole in the three string modes, and the fixtures of up to 32 KiB cut at 24
/// lengths over the text and at each of its last 8 lengths (one mode per cut, in turn). The
/// digests are printed (`--nocapture`) to compare the walks of two builds.
#[test]
fn every_source_reads_the_same() {
    let texts = fixture_texts();
    assert!(texts.len() > 500, "{} fixtures", texts.len());
    for (label, text) in &texts {
        for mode in 0..3 {
            let digest = walk_every_source(label, text, mode);
            println!("walk {label} whole {mode} {digest:016x}");
        }
        if text.len() > 32 * 1024 {
            continue;
        }
        let mut cuts: Vec<usize> = (1..24).map(|k| text.len() * k / 24).collect();
        cuts.extend(text.len().saturating_sub(8)..text.len());
        cuts.sort_unstable();
        cuts.dedup();
        for (k, &cut) in cuts.iter().enumerate() {
            let mode = k % 3;
            let digest = walk_every_source(&format!("{label} cut {cut}"), &text[..cut], mode);
            println!("walk {label} {cut} {mode} {digest:016x}");
        }
    }
}

/// A value split across the end of a read at every byte: numbers, names, strings across lines
/// with `""`, comments, octet strings and REALs, whole and cut at every length.
#[test]
fn values_split_at_every_byte() {
    let text: &[u8] = b"Demo ::= {\r\n  -- a comment -- name [12] \"one\"\"two\r\nthree\" ,\n\
        x-y 123456789012 , z '0A1b'H, real { 314159, 10, -5 }, k MINUS-INFINITY,\n\
        t TRUE , f FALSE , n NULL , big 99999999999999999999 , u +42 , e enum-val -- c\n } ";
    for cut in 0..=text.len() {
        for mode in 0..3 {
            walk_every_source(&format!("demo cut {cut}"), &text[..cut], mode);
        }
    }
    let walked = walk(SliceSource::new(text), 0);
    assert!(
        walked.iter().any(|step| step.contains("123456789012")),
        "{walked:?}"
    );
}

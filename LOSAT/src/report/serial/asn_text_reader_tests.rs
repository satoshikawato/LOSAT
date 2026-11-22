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
        b'{' if step.is_multiple_of(11) => {
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
            if step.is_multiple_of(3) {
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

/// Every fixture whole in the three string modes, and the fixtures of up to 32 KiB cut at 23
/// lengths over the text (k/24 of it for k = 1..23) and at each of its last 8 lengths, at most 31
/// cuts (one mode per cut, in turn). The test checks that the three sources agree; it does not
/// compare with another build. The digests are printed (`--nocapture`) so that the walks of two
/// builds can be compared by hand; `literal_walks` holds results of dev's reader.
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

/// One operation of a literal walk, named as in `LITERAL_WALKS`, and what it read. `blocks` and
/// `members` are the blocks and SEQUENCEs open at that point.
fn apply<S: AsnSource>(
    r: &mut Reader<S>,
    op: &'static str,
    blocks: &mut Vec<bool>,
    members: &mut Vec<Members>,
) -> Result<String> {
    let (name, arg) = op.split_once(':').unwrap_or((op, ""));
    let names: Vec<&'static str> = arg.split(',').collect();
    let values: Vec<(&'static str, i32)> = names
        .iter()
        .filter_map(|value| value.split_once('='))
        .map(|(name, n)| (name, n.parse().unwrap()))
        .collect();
    Ok(match name {
        "ws" => format!("ws {:?}", r.skip_white_space()? as char),
        "int4" => format!("int4 {}", r.int4()?),
        "uint4" => format!("uint4 {}", r.uint4()?),
        "int8" => format!("int8 {}", r.int8()?),
        "string" => format!("string {}", show(&r.string_as_read()?)),
        "visible" => format!("visible {}", show(&r.string()?)),
        "id" => match r.member_or_tag()? {
            MemberId::Name(id) => format!("id {}", show(r.bytes(id))),
            MemberId::Tag(id) => format!("tag {}", show(r.bytes(id))),
        },
        "type" => match r.type_id()? {
            Some(id) => format!("type {}", show(r.bytes(id))),
            None => "type none".into(),
        },
        "ptr" => format!("ptr {:?}", r.pointer_type()?),
        "choice" => format!("choice {}", r.choice(&names)?),
        "enum" => format!("enum {:?}", r.enumerated(&values)?),
        "{" => {
            blocks.push(r.start_block()?);
            "{".into()
        }
        "next" => format!("next {}", r.next_element(blocks.last_mut().unwrap())?),
        "}" => {
            r.end_block()?;
            blocks.pop();
            "}".into()
        }
        "members" => {
            members.push(r.begin_members()?);
            "members".into()
        }
        "member" => {
            let member = r.next_member(members.last_mut().unwrap(), &names)?;
            if member.is_none() {
                members.pop();
            }
            format!("member {member:?}")
        }
        "each" => {
            let mut elements = Vec::new();
            r.each(|r| {
                elements.push(r.int4()?);
                Ok(())
            })?;
            format!("each {elements:?}")
        }
        "skip" => {
            r.skip_value()?;
            "skip".into()
        }
        other => panic!("no operation {other}"),
    })
}

/// The steps of the operations `ops` on `src` with the string rule `mode` (as `walk`), up to the
/// first error.
fn literal_walk<S: AsnSource>(src: S, mode: usize, ops: &'static str) -> Vec<String> {
    let mut r = Reader::new(src, "text", "test reader");
    if mode >= 1 {
        r = r.refuse_invisible_bytes();
    }
    if mode == 2 {
        r = r.defer_refusals();
    }
    let (mut blocks, mut members, mut log) = (Vec::new(), Vec::new(), Vec::new());
    for op in ops.split_whitespace() {
        match apply(&mut r, op, &mut blocks, &mut members) {
            Ok(what) => note(&r, &mut log, what),
            Err(error) => {
                let what = format!("error {error:?} refusal {:?}", r.refusal);
                note(&r, &mut log, what);
                return log;
            }
        }
    }
    note(&r, &mut log, format!("end refusal {:?}", r.refusal));
    log
}

/// Texts, string rules and operations, with the steps that the reader gave for them before its
/// position moved into locals (recorded from dev's reader, not from this one): line ends of every form in white space, comments and
/// strings, a comment and a `-` at the end of the text, the stop positions of the overflows,
/// errors at their positions, and the readers of members, variants, enumerated values, pointer
/// types, type names and lists. Each step is `what was read @position:line`.
const LITERAL_WALKS: &[(&[u8], usize, &str, &[&str])] = &[
    (
        b"{\n1,\r2,\r\n3,\n\r4,\r\r5,\t\n 6 }",
        0,
        "{ next int4 next int4 next int4 next int4 next int4 next int4 next }",
        &[
            "{ @1:1",
            "next true @2:2",
            "int4 1 @3:2",
            "next true @4:2",
            "int4 2 @6:3",
            "next true @7:3",
            "int4 3 @10:4",
            "next true @11:4",
            "int4 4 @14:5",
            "next true @15:5",
            "int4 5 @18:7",
            "next true @19:7",
            "int4 6 @23:8",
            "next false @24:8",
            "} @25:8",
            "end refusal None @25:8",
        ],
    ),
    (
        b"{ -- a\n1 -- b -- , --\r\n2 --c-\n\r }",
        0,
        "{ next int4 next int4 next }",
        &[
            "{ @1:1",
            "next true @7:2",
            "int4 1 @8:2",
            "next true @18:2",
            "int4 2 @24:3",
            "next false @32:4",
            "} @33:4",
            "end refusal None @33:4",
        ],
    ),
    (
        b"1 -- tail",
        0,
        "int4 ws",
        &[
            "int4 1 @1:1",
            "error AsnTextError { line: 1, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } refusal None @9:1",
        ],
    ),
    (
        b"1 \n-",
        0,
        "int4 ws",
        &[
            "int4 1 @1:1",
            "error AsnTextError { line: 2, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } refusal None @3:2",
        ],
    ),
    (
        b" \n-5 - 3",
        0,
        "int4 int4",
        &[
            "int4 -5 @4:2",
            "error AsnTextError { line: 2, message: \"bad number\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 574, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::BadNumber\") })) } refusal None @6:2",
        ],
    ),
    (
        b"\"a\nb\r\nc\n\rd\" \"e\"\"f\"",
        0,
        "string visible",
        &[
            "string \"abcd\" @11:4",
            "visible \"e\\\"f\" @18:4",
            "end refusal None @18:4",
        ],
    ),
    (
        b"\"abc\ndef",
        0,
        "string",
        &[
            "error AsnTextError { line: 2, message: \"unclosed string starts at line 1\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 869, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::UnendedString\") })) } refusal None @8:2",
        ],
    ),
    (
        b"\"a\x01b\"\n5",
        1,
        "visible int4",
        &[
            "error AsnTextError { line: 1, message: \"a string with the byte 0x1 (NCBI replaces it by '#' and posts an error) is not supported by LOSAT's test reader\", kind: Unsupported, ncbi: None } refusal None @3:1",
        ],
    ),
    (
        b"\"a\x01b\"\n5",
        2,
        "visible int4",
        &[
            "visible \"a\\u{1}b\" @5:1",
            "int4 5 @7:2",
            "end refusal Some(AsnTextError { line: 1, message: \"a string with the byte 0x1 (NCBI replaces it by '#' and posts an error) is not supported by LOSAT's test reader\", kind: Unsupported, ncbi: None }) @7:2",
        ],
    ),
    (
        b"2147483647 2147483648",
        0,
        "int4 int4",
        &[
            "int4 2147483647 @10:1",
            "error AsnTextError { line: 1, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } refusal None @21:1",
        ],
    ),
    (
        b"-2147483648\n-2147483649",
        0,
        "int4 int4",
        &[
            "int4 -2147483648 @11:1",
            "error AsnTextError { line: 2, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } refusal None @22:2",
        ],
    ),
    (
        b"4294967295 +7 4294967296",
        0,
        "uint4 uint4 uint4",
        &[
            "uint4 4294967295 @10:1",
            "uint4 7 @13:1",
            "error AsnTextError { line: 1, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } refusal None @24:1",
        ],
    ),
    (
        b"\r-1",
        0,
        "uint4",
        &[
            "error AsnTextError { line: 2, message: \"bad number\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 574, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::BadNumber\") })) } refusal None @1:2",
        ],
    ),
    (
        b"9223372036854775807 -9223372036854775808 9223372036854775808",
        0,
        "int8 int8 int8",
        &[
            "int8 9223372036854775807 @19:1",
            "int8 -9223372036854775808 @40:1",
            "error AsnTextError { line: 1, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } refusal None @60:1",
        ],
    ),
    (
        b"12",
        0,
        "int4 ws",
        &[
            "int4 12 @2:1",
            "error AsnTextError { line: 1, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } refusal None @2:1",
        ],
    ),
    (
        b"{ 1\n 2 }",
        0,
        "{ next int4 next",
        &[
            "{ @1:1",
            "next true @2:1",
            "int4 1 @3:1",
            "error AsnTextError { line: 2, message: \"',' or '}' expected\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 1220, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::NextElement\") })) } refusal None @5:2",
        ],
    ),
    (
        b"{ a 1, c 2,\n b 3 }",
        0,
        "members member:a,b,c int4 member:a,b,c int4 member:a,b,c",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "int4 1 @5:1",
            "member Some(2) @8:1",
            "int4 2 @10:1",
            "error AsnTextError { line: 2, message: \"member b out of order\", kind: Format, ncbi: None } refusal None @14:2",
        ],
    ),
    (
        b"{ a 1,\n x 2 }",
        0,
        "members member:a,b int4 member:a,b",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "int4 1 @5:1",
            "error AsnTextError { line: 2, message: \"unexpected member x\", kind: Format, ncbi: None } refusal None @9:2",
        ],
    ),
    (
        b"{\n [1] 5 }",
        0,
        "members member:a",
        &[
            "members @1:1",
            "error AsnTextError { line: 2, message: \"a member or variant named by its tag is not supported by LOSAT's test reader\", kind: Unsupported, ncbi: None } refusal None @6:2",
        ],
    ),
    (
        b"{ a 1\n}",
        0,
        "members member:a,b int4 member:a,b",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "int4 1 @5:1",
            "member None @7:2",
            "end refusal None @7:2",
        ],
    ),
    (
        b"b 5",
        0,
        "choice:a,b int4",
        &[
            "choice 1 @1:1",
            "int4 5 @3:1",
            "end refusal None @3:1",
        ],
    ),
    (
        b"\nz 5",
        0,
        "choice:a,b",
        &[
            "error AsnTextError { line: 2, message: \"unexpected variant z\", kind: Format, ncbi: None } refusal None @2:2",
        ],
    ),
    (
        b" }",
        0,
        "choice:a",
        &[
            "error AsnTextError { line: 1, message: \"choice variant id expected\", kind: Format, ncbi: None } refusal None @1:1",
        ],
    ),
    (
        b"red 2 green\n9",
        0,
        "enum:red=0,green=2 enum:red=0,green=2 enum:red=0,green=2 enum:red=0,green=2",
        &[
            "enum (\"red\", 0) @3:1",
            "enum (\"green\", 2) @5:1",
            "enum (\"green\", 2) @11:1",
            "error AsnTextError { line: 2, message: \"invalid value of enumerated type: 9\", kind: Format, ncbi: Some(Plain(ThrowSite { file: \"src/serial/enumerated.cpp\", line: 156, kind: \"CSerialException::eInvalidData\", function: None })) } refusal None @13:2",
        ],
    ),
    (
        b"NULL @ :\n5",
        0,
        "ptr ptr ptr ptr int4",
        &[
            "ptr Null @4:1",
            "ptr Object @6:1",
            "ptr Other @8:1",
            "ptr This @9:2",
            "int4 5 @10:2",
            "end refusal None @10:2",
        ],
    ),
    (
        b"Seq-id x",
        0,
        "type",
        &[
            "type \"Seq-id\" @6:1",
            "end refusal None @6:1",
        ],
    ),
    (
        b"[abc]",
        0,
        "type",
        &[
            "type \"ab\" @4:1",
            "end refusal None @4:1",
        ],
    ),
    (
        b"[]",
        0,
        "type",
        &[
            "type none @1:1",
            "end refusal None @1:1",
        ],
    ),
    (
        b"[x\ny]",
        0,
        "type",
        &[
            "error AsnTextError { line: 1, message: \"end of line: expected ']'\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 299, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::ReadTypeId\") })) } refusal None @0:1",
        ],
    ),
    (
        b"[name] [7] x-y-1 ab-",
        0,
        "id id id id",
        &[
            "id \"name\" @6:1",
            "tag \"7\" @10:1",
            "id \"x-y-1\" @16:1",
            "error AsnTextError { line: 1, message: \"unexpected end of the text\", kind: Eof, ncbi: None } refusal None @17:1",
        ],
    ),
    (
        b"{ 1, -2 ,3 }",
        0,
        "each",
        &[
            "each [1, -2, 3] @12:1",
            "end refusal None @12:1",
        ],
    ),
    (
        b"{ 1\n 2 }",
        0,
        "each",
        &[
            "error AsnTextError { line: 2, message: \"',' or '}' expected\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 1220, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::NextElement\") })) } refusal None @5:2",
        ],
    ),
    (
        b"{ a { 1, \"}\" }, -- } --\n b 'FF'H } 7",
        0,
        "skip int4",
        &[
            "error AsnTextError { line: 2, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } refusal None @36:2",
        ],
    ),
];

/// `LITERAL_WALKS`: each walk is the same from the three sources and gives dev's steps. On a
/// difference the message lists every walk as read.
#[test]
fn literal_walks() {
    let mut read = String::new();
    let mut differ = Vec::new();
    for (k, &(text, mode, ops, expected)) in LITERAL_WALKS.iter().enumerate() {
        let slice = literal_walk(SliceSource::new(text), mode, ops);
        let first = text.len().min(1);
        let blocks = literal_walk(StreamSource::new(&text[first..], &text[..first]), mode, ops);
        assert_eq!(slice, blocks, "walk {k}: slice and stream differ");
        let bytes = literal_walk(StreamSource::new(OneByte(text), &[]), mode, ops);
        assert_eq!(slice, bytes, "walk {k}: slice and one byte per read differ");
        if slice != expected {
            differ.push(k);
        }
        read += &format!("walk {k}: {slice:?}\n");
    }
    assert!(differ.is_empty(), "walks {differ:?} differ; read:\n{read}");
}

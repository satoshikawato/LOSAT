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

/// A source of `text` that reads one byte per fill, so that the bytes it keeps end at every byte
/// of the text.
pub(crate) fn one_byte_per_fill(text: &[u8]) -> StreamSource<OneByte<'_>> {
    let mut src = StreamSource::new(OneByte(text), &[]);
    src.fill_one_byte_at_a_time();
    src
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
        skip_value(r)?;
        note(r, log, "skipped deep".into());
        return Ok(());
    }
    match c {
        b'{' if matches!(name, b"real" | b"lambda" | b"k" | b"h") => {
            let v = r.real()?;
            note(r, log, format!("real {:x}", v.to_bits()));
        }
        b'{' if step.is_multiple_of(11) => {
            skip_value(r)?;
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
                let text = r.utf8_string()?;
                note(r, log, format!("string {}", show(&text)));
            }
            1 => {
                let text = r.string(&Frame::ROOT, &mut Noted)?;
                note(r, log, format!("visible {}", show(&text)));
            }
            n => {
                let (mut out, mut bad) = (Vec::new(), Vec::new());
                r.visible_string(&mut out, true, n == 3, &mut |_, line, byte, piece| {
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
            skip_value(r)?;
            note(r, log, "skipped".into());
        }
    }
    Ok(())
}

/// Skips one value of any type with the lexical readers: a block with what it holds, a string, an
/// OCTET STRING or BIT STRING, or the words and numbers up to the `,` or `}` that ends it (a
/// CHOICE is its variant's name and value), checking the blocks, strings and octet strings of the
/// value but not the grammar of its type (no reader of LOSAT skips a value so: NCBI's readers read
/// each value by its type, `asn_types::read`).
fn skip_value<S: AsnSource>(r: &mut Reader<S>) -> Result<()> {
    let mut depth = 0usize;
    let mut any = false;
    loop {
        let c = r.skip_white_space()?;
        match c {
            b'{' => {
                r.pos += 1;
                depth += 1;
            }
            b'}' | b',' if depth == 0 => {
                if !any {
                    return Err(r.error("value expected"));
                }
                return Ok(());
            }
            b'}' => {
                r.pos += 1;
                depth -= 1;
            }
            b'"' => {
                r.utf8_string()?;
            }
            b'\'' => {
                r.pos += 1;
                loop {
                    let Some(c) = r.src.byte_at(r.pos) else {
                        return Err(r.end_of_text());
                    };
                    r.pos += 1;
                    match c {
                        b'\'' => break,
                        b'\r' | b'\n' => r.skip_end_of_line(c),
                        _ => {}
                    }
                }
                if matches!(r.peek_at(0), b'H' | b'B') {
                    r.pos += 1;
                }
            }
            _ => r.pos += 1,
        }
        any = true;
    }
}

thread_local! {
    /// The posts that the strings of this thread's walk gave (`Noted`).
    static NOTED: std::cell::RefCell<Vec<ReplacedChar>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Where a walk's strings give their posts (`Reader::string`): they are kept for `replaced`.
struct Noted;

impl ReplacedSink for Noted {
    fn replaced(&mut self, replaced: ReplacedChar) {
        NOTED.with(|noted| noted.borrow_mut().push(replaced));
    }
}

/// The posts that the strings read so far gave (`Noted`): the trace, the byte, its position and
/// the piece of each.
fn replaced() -> String {
    let posts: Vec<String> = NOTED.with(|noted| {
        noted
            .borrow()
            .iter()
            .map(|post| {
                format!(
                    "{} {:x} {} {}",
                    post.trace,
                    post.byte,
                    post.position,
                    show(&post.string)
                )
            })
            .collect()
    });
    format!("{posts:?}")
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

/// The walk of `src`: every step, and the error and the posts of its strings it ended with.
fn walk<S: AsnSource>(src: S) -> Vec<String> {
    NOTED.with(|noted| noted.borrow_mut().clear());
    let mut r = Reader::new(src, "text", "test reader");
    let mut log = Vec::new();
    let error = objects(&mut r, &mut log).expect_err("a walk ends with an error");
    note(
        &r,
        &mut log,
        format!("error {error:?} replaced {}", replaced()),
    );
    log
}

/// The walks of `text` from the three sources, which must be equal; returns a digest of the walk.
fn walk_every_source(label: &str, text: &[u8]) -> u64 {
    let slice = walk(SliceSource::new(text));
    let first = text.len().min(1);
    let blocks = walk(StreamSource::new(&text[first..], &text[..first]));
    assert_eq!(slice, blocks, "{label}: slice and stream differ");
    let bytes = walk(one_byte_per_fill(text));
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

/// Every fixture whole, and the fixtures of up to 32 KiB cut at 23 lengths over the text (k/24 of
/// it for k = 1..23) and at each of its last 8 lengths, at most 31 cuts. The test checks that the three sources agree; it does not
/// compare with another build. The digests are printed (`--nocapture`) so that the walks of two
/// builds can be compared by hand; `literal_walks` holds results of dev's reader.
#[test]
fn every_source_reads_the_same() {
    let texts = fixture_texts();
    assert!(texts.len() > 500, "{} fixtures", texts.len());
    for (label, text) in &texts {
        let digest = walk_every_source(label, text);
        println!("walk {label} whole {digest:016x}");
        if text.len() > 32 * 1024 {
            continue;
        }
        let mut cuts: Vec<usize> = (1..24).map(|k| text.len() * k / 24).collect();
        cuts.extend(text.len().saturating_sub(8)..text.len());
        cuts.sort_unstable();
        cuts.dedup();
        for &cut in &cuts {
            let digest = walk_every_source(&format!("{label} cut {cut}"), &text[..cut]);
            println!("walk {label} {cut} {digest:016x}");
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
        walk_every_source(&format!("demo cut {cut}"), &text[..cut]);
    }
    let walked = walk(SliceSource::new(text));
    assert!(
        walked.iter().any(|step| step.contains("123456789012")),
        "{walked:?}"
    );
}

/// A frame open in a literal walk (`Frame`, kept by its fields between the operations): a
/// SEQUENCE being read (`members`), a list (`list`) or a CHOICE's variant (`choice`).
struct Level {
    /// Whether it is a SEQUENCE or a list (whose `member` or `item` reads the next of it).
    block: bool,
    block_start: bool,
    names: &'static [&'static str],
    next: usize,
}

/// `f` with the frame of the innermost of `levels`, each held by the one before it, the first by
/// `at`.
fn in_frames<R>(levels: &[Level], at: &Frame<'_>, f: impl FnOnce(&Frame<'_>) -> R) -> R {
    match levels.split_first() {
        None => f(at),
        Some((level, rest)) => {
            let frame = Frame {
                parent: Some(at),
                names: level.names,
                next: level.next,
            };
            in_frames(rest, &frame, f)
        }
    }
}

/// The names of an operation's argument (`a,b`), kept for the walk.
fn leaked_names(arg: &str) -> &'static [&'static str] {
    let names: Vec<&'static str> = arg
        .split(',')
        .map(|name| &*Box::leak(name.to_string().into_boxed_str()))
        .collect();
    Box::leak(names.into_boxed_slice())
}

/// One operation of a literal walk, named as in `LITERAL_WALKS`, and what it read. `blocks` are
/// the blocks open at that point, `levels` the frames (`member` and `item` read the next member
/// or element of the innermost SEQUENCE or list, whose variants end there).
fn apply<S: AsnSource>(
    r: &mut Reader<S>,
    op: &'static str,
    blocks: &mut Vec<bool>,
    levels: &mut Vec<Level>,
) -> Result<String> {
    let (name, arg) = op.split_once(':').unwrap_or((op, ""));
    let names = leaked_names(arg);
    let values: Vec<(&'static str, i32)> = names
        .iter()
        .filter_map(|value| value.split_once('='))
        .map(|(name, n)| (name, n.parse().unwrap()))
        .collect();
    if matches!(name, "member" | "item") {
        let block = levels.iter().rposition(|level| level.block);
        levels.truncate(block.map_or(0, |block| block + 1));
    }
    Ok(match name {
        "ws" => format!("ws {:?}", r.skip_white_space()? as char),
        "int4" => format!("int4 {}", r.int4()?),
        "uint4" => format!("uint4 {}", r.uint4()?),
        "int8" => format!("int8 {}", r.int8()?),
        "string" => format!("string {}", show(&r.utf8_string()?)),
        "visible" => {
            let text = in_frames(levels, &Frame::ROOT, |at| r.string(at, &mut Noted))?;
            format!("visible {}", show(&text))
        }
        "id" => match r.member_or_tag()? {
            MemberId::Name(id) => format!("id {}", show(r.bytes(id))),
            MemberId::Tag(id) => format!("tag {}", show(r.bytes(id))),
        },
        "type" => match r.type_id()? {
            Some(id) => format!("type {}", show(r.bytes(id))),
            None => "type none".into(),
        },
        "ptr" => format!("ptr {:?}", r.pointer_type()?),
        "choice" => {
            let next = in_frames(levels, &Frame::ROOT, |at| {
                r.choice(at, names).map(|variant| variant.next)
            })?;
            levels.push(Level {
                block: false,
                block_start: false,
                names,
                next,
            });
            format!("choice {}", next - 1)
        }
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
        "members" | "list" => {
            let block_start = in_frames(levels, &Frame::ROOT, |at| {
                if name == "list" {
                    r.begin_list(at).map(|list| list.block_start)
                } else {
                    r.begin_members(at, names)
                        .map(|members| members.block_start)
                }
            })?;
            levels.push(Level {
                block: true,
                block_start,
                names: if name == "list" { &[] } else { names },
                next: 0,
            });
            name.into()
        }
        "member" | "item" => {
            let (level, outer) = levels.split_last_mut().expect("a SEQUENCE or list is open");
            let read = in_frames(outer, &Frame::ROOT, |at| {
                let frame = Frame {
                    parent: Some(at),
                    names: level.names,
                    next: level.next,
                };
                if name == "item" {
                    let mut list = List {
                        block_start: level.block_start,
                        frame,
                    };
                    let more = r.next_item(&mut list);
                    level.block_start = list.block_start;
                    more.map(|more| format!("item {more}"))
                } else {
                    let mut members = Members {
                        block_start: level.block_start,
                        frame,
                    };
                    let member = r.next_member(&mut members);
                    (level.block_start, level.next) = (members.block_start, members.frame.next);
                    member.map(|member| format!("member {member:?}"))
                }
            })?;
            if read.ends_with("None") || read.ends_with("false") {
                levels.pop();
            }
            read
        }
        "each" => {
            let mut elements = Vec::new();
            r.each(&Frame::ROOT, |r, _| {
                elements.push(r.int4()?);
                Ok(())
            })?;
            format!("each {elements:?}")
        }
        "skip" => {
            skip_value(r)?;
            "skip".into()
        }
        other => panic!("no operation {other}"),
    })
}

/// The steps of the operations `ops` on `src` (as `walk`), up to the first error, with the
/// reader's hand grammar rooted at `Root`.
fn literal_walk<S: AsnSource>(src: S, ops: &'static str) -> Vec<String> {
    NOTED.with(|noted| noted.borrow_mut().clear());
    let mut r = Reader::new(src, "text", "test reader");
    r.root = "Root";
    let (mut blocks, mut levels, mut log) = (Vec::new(), Vec::new(), Vec::new());
    for op in ops.split_whitespace() {
        match apply(&mut r, op, &mut blocks, &mut levels) {
            Ok(what) => note(&r, &mut log, what),
            Err(error) => {
                let what = format!("error {error:?} replaced {}", replaced());
                note(&r, &mut log, what);
                return log;
            }
        }
    }
    note(&r, &mut log, format!("end replaced {}", replaced()));
    log
}

/// Texts, string rules and operations, with the steps that the reader gave for them before its
/// position moved into locals (recorded from dev's reader, not from this one): line ends of every form in white space, comments and
/// strings, a comment and a `-` at the end of the text, the stop positions of the overflows,
/// errors at their positions (those of the numbers since the readers of numbers became one, at
/// NCBI's: `GetInt4` and `GetUint4` take the byte after a sign, the byte that is no digit and the
/// digit that overflows before they throw), and the readers of members, variants, enumerated values, pointer
/// types, type names and lists. Each step is `what was read @position:line`; the last also shows
/// the posts that the strings gave (`Noted`: trace, byte, position, piece). Since item
/// 9(a) of line RG-KB (a VisibleString of a hand grammar reads a byte outside 0x20..0x7E as `#`
/// with its post, where it was refused, and a UTF8String is read by `visible_string`, which leaves
/// the position of an unended string at the start of its last run) the walks of those strings and
/// of their frames are this reader's.
const LITERAL_WALKS: &[(&[u8], &str, &[&str])] = &[
    (
        b"{\n1,\r2,\r\n3,\n\r4,\r\r5,\t\n 6 }",
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
            "end replaced [] @25:8",
        ],
    ),
    (
        b"{ -- a\n1 -- b -- , --\r\n2 --c-\n\r }",
        "{ next int4 next int4 next }",
        &[
            "{ @1:1",
            "next true @7:2",
            "int4 1 @8:2",
            "next true @18:2",
            "int4 2 @24:3",
            "next false @32:4",
            "} @33:4",
            "end replaced [] @33:4",
        ],
    ),
    (
        b"1 -- tail",
        "int4 ws",
        &[
            "int4 1 @1:1",
            "error AsnTextError { line: 1, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } replaced [] @9:1",
        ],
    ),
    (
        b"1 \n-",
        "int4 ws",
        &[
            "int4 1 @1:1",
            "error AsnTextError { line: 2, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } replaced [] @3:2",
        ],
    ),
    (
        b" \n-5 - 3",
        "int4 int4",
        &[
            "int4 -5 @4:2",
            "error AsnTextError { line: 2, message: \"bad number\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 574, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::BadNumber\") })) } replaced [] @7:2",
        ],
    ),
    (
        b"\"a\nb\r\nc\n\rd\" \"e\"\"f\"",
        "string visible",
        &[
            "string \"abcd\" @11:4",
            "visible \"e\\\"f\" @18:4",
            "end replaced [] @18:4",
        ],
    ),
    (
        b"\"abc\ndef",
        "string",
        &[
            "error AsnTextError { line: 2, message: \"unclosed string starts at line 1\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 869, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::UnendedString\") })) } replaced [] @5:2",
        ],
    ),
    (
        b"\"a\x01b\"\n5",
        "visible int4",
        &[
            "visible \"a#b\" @5:1",
            "int4 5 @7:2",
            "end replaced [\"Root 1 line 1 \\\"a\\\\u{1}b\\\"\"] @7:2",
        ],
    ),
    (
        b"{ a { \"x\x7f\", \"\x80\x81\" }, b c \"\x01\" } \"\x02\"",
        "members:a,b member list item visible item visible item } member choice:c,d visible member visible",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "list @5:1",
            "item true @6:1",
            "visible \"x#\" @10:1",
            "item true @11:1",
            "visible \"##\" @16:1",
            "item false @17:1",
            "} @18:1",
            "member Some(1) @21:1",
            "choice 0 @23:1",
            "visible \"#\" @27:1",
            "member None @29:1",
            "visible \"#\" @33:1",
            "end replaced [\"Root.a.E 7f line 1 \\\"x\\\\u{7f}\\\"\", \"Root.a.E 80 line 1 \\\"\u{fffd}\u{fffd}\\\"\", \"Root.a.E 81 line 1 \\\"\u{fffd}\u{fffd}\\\"\", \"Root.b.c 1 line 1 \\\"\\\\u{1}\\\"\", \"Root 2 line 1 \\\"\\\\u{2}\\\"\"] @33:1",
        ],
    ),
    (
        b"2147483647 2147483648",
        "int4 int4",
        &[
            "int4 2147483647 @10:1",
            "error AsnTextError { line: 1, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } replaced [] @21:1",
        ],
    ),
    (
        b"-2147483648\n-2147483649",
        "int4 int4",
        &[
            "int4 -2147483648 @11:1",
            "error AsnTextError { line: 2, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } replaced [] @23:2",
        ],
    ),
    (
        b"4294967295 +7 4294967296",
        "uint4 uint4 uint4",
        &[
            "uint4 4294967295 @10:1",
            "uint4 7 @13:1",
            "error AsnTextError { line: 1, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } replaced [] @24:1",
        ],
    ),
    (
        b"\r-1",
        "uint4",
        &[
            "error AsnTextError { line: 2, message: \"bad number\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 574, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::BadNumber\") })) } replaced [] @2:2",
        ],
    ),
    (
        b"9223372036854775807 -9223372036854775808 9223372036854775808",
        "int8 int8 int8",
        &[
            "int8 9223372036854775807 @19:1",
            "int8 -9223372036854775808 @40:1",
            "error AsnTextError { line: 1, message: \"number overflow\", kind: Format, ncbi: Some(Buffer(ThrowSite { file: \"src/util/strbuffer.cpp\", line: 583, kind: \"CUtilException::eWrongData\", function: Some(\"ncbi::CIStreamBuffer::NumberOverflow\") })) } replaced [] @60:1",
        ],
    ),
    (
        b"12",
        "int4 ws",
        &[
            "int4 12 @2:1",
            "error AsnTextError { line: 1, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } replaced [] @2:1",
        ],
    ),
    (
        b"{ 1\n 2 }",
        "{ next int4 next",
        &[
            "{ @1:1",
            "next true @2:1",
            "int4 1 @3:1",
            "error AsnTextError { line: 2, message: \"',' or '}' expected\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 1220, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::NextElement\") })) } replaced [] @5:2",
        ],
    ),
    (
        b"{ a 1, c 2,\n b 3 }",
        "members:a,b,c member int4 member int4 member",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "int4 1 @5:1",
            "member Some(2) @8:1",
            "int4 2 @10:1",
            "error AsnTextError { line: 2, message: \"member b out of order\", kind: Format, ncbi: None } replaced [] @14:2",
        ],
    ),
    (
        b"{ a 1,\n x 2 }",
        "members:a,b member int4 member",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "int4 1 @5:1",
            "error AsnTextError { line: 2, message: \"unexpected member x\", kind: Format, ncbi: None } replaced [] @9:2",
        ],
    ),
    (
        b"{\n [1] 5 }",
        "members:a member",
        &[
            "members @1:1",
            "error AsnTextError { line: 2, message: \"a member or variant named by its tag is not supported by LOSAT's test reader\", kind: Unsupported, ncbi: None } replaced [] @6:2",
        ],
    ),
    (
        b"{ a 1\n}",
        "members:a,b member int4 member",
        &[
            "members @1:1",
            "member Some(0) @3:1",
            "int4 1 @5:1",
            "member None @7:2",
            "end replaced [] @7:2",
        ],
    ),
    (
        b"b 5",
        "choice:a,b int4",
        &[
            "choice 1 @1:1",
            "int4 5 @3:1",
            "end replaced [] @3:1",
        ],
    ),
    (
        b"\nz 5",
        "choice:a,b",
        &[
            "error AsnTextError { line: 2, message: \"unexpected variant z\", kind: Format, ncbi: None } replaced [] @2:2",
        ],
    ),
    (
        b" }",
        "choice:a",
        &[
            "error AsnTextError { line: 1, message: \"choice variant id expected\", kind: Format, ncbi: None } replaced [] @1:1",
        ],
    ),
    (
        b"red 2 green\n9",
        "enum:red=0,green=2 enum:red=0,green=2 enum:red=0,green=2 enum:red=0,green=2",
        &[
            "enum (\"red\", 0) @3:1",
            "enum (\"green\", 2) @5:1",
            "enum (\"green\", 2) @11:1",
            "error AsnTextError { line: 2, message: \"invalid value of enumerated type: 9\", kind: Format, ncbi: Some(Plain(ThrowSite { file: \"src/serial/enumerated.cpp\", line: 156, kind: \"CSerialException::eInvalidData\", function: None })) } replaced [] @13:2",
        ],
    ),
    (
        b"NULL @ :\n5",
        "ptr ptr ptr ptr int4",
        &[
            "ptr Null @4:1",
            "ptr Object @6:1",
            "ptr Other @8:1",
            "ptr This @9:2",
            "int4 5 @10:2",
            "end replaced [] @10:2",
        ],
    ),
    (
        b"Seq-id x",
        "type",
        &[
            "type \"Seq-id\" @6:1",
            "end replaced [] @6:1",
        ],
    ),
    (
        b"[abc]",
        "type",
        &[
            "type \"ab\" @4:1",
            "end replaced [] @4:1",
        ],
    ),
    (
        b"[]",
        "type",
        &[
            "type none @1:1",
            "end replaced [] @1:1",
        ],
    ),
    (
        b"[x\ny]",
        "type",
        &[
            "error AsnTextError { line: 1, message: \"end of line: expected ']'\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 299, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::ReadTypeId\") })) } replaced [] @0:1",
        ],
    ),
    (
        b"[name] [7] x-y-1 ab-",
        "id id id id",
        &[
            "id \"name\" @6:1",
            "tag \"7\" @10:1",
            "id \"x-y-1\" @16:1",
            "error AsnTextError { line: 1, message: \"unexpected end of the text\", kind: Eof, ncbi: None } replaced [] @17:1",
        ],
    ),
    (
        b"{ 1, -2 ,3 }",
        "each",
        &[
            "each [1, -2, 3] @12:1",
            "end replaced [] @12:1",
        ],
    ),
    (
        b"{ 1\n 2 }",
        "each",
        &[
            "error AsnTextError { line: 2, message: \"',' or '}' expected\", kind: Format, ncbi: Some(Stream(ThrowSite { file: \"src/serial/objistrasn.cpp\", line: 1220, kind: \"CSerialException::eFormatError\", function: Some(\"ncbi::CObjectIStreamAsn::NextElement\") })) } replaced [] @5:2",
        ],
    ),
    (
        b"{ a { 1, \"}\" }, -- } --\n b 'FF'H } 7",
        "skip int4",
        &[
            "error AsnTextError { line: 2, message: \"unexpected end of the text\", kind: Eof, ncbi: Some(WhiteSpaceEof) } replaced [] @36:2",
        ],
    ),
];

/// `LITERAL_WALKS`: each walk is the same from the three sources and gives dev's steps. On a
/// difference the message lists every walk as read.
#[test]
fn literal_walks() {
    let mut read = String::new();
    let mut differ = Vec::new();
    for (k, &(text, ops, expected)) in LITERAL_WALKS.iter().enumerate() {
        let slice = literal_walk(SliceSource::new(text), ops);
        let first = text.len().min(1);
        let blocks = literal_walk(StreamSource::new(&text[first..], &text[..first]), ops);
        assert_eq!(slice, blocks, "walk {k}: slice and stream differ");
        let bytes = literal_walk(one_byte_per_fill(text), ops);
        assert_eq!(slice, bytes, "walk {k}: slice and one byte per read differ");
        if slice != expected {
            differ.push(k);
        }
        read += &format!("walk {k}: {slice:?}\n");
    }
    assert!(differ.is_empty(), "walks {differ:?} differ; read:\n{read}");
}

/// `StreamSource::read_into` past what the source keeps reads the bytes straight from the input
/// into `out`, reserved at their count first as NCBI reserves them (one block of their size, not
/// one that doubles as it fills): 100,000 bytes, more than a block of the source
/// (`STREAM_BLOCK`), through a reader that gives one byte per call; a count past the end of the
/// input reads what there is.
#[test]
fn a_long_read_comes_from_the_input_into_one_block_of_its_size() {
    let data: Vec<u8> = (0..100_000u32).map(|i| (i % 251) as u8).collect();
    assert!(data.len() > STREAM_BLOCK);
    let mut src = StreamSource::new(OneByte(&data), &[]);
    let mut out = Vec::new();
    assert_eq!(src.read_into(0, 100_000, &mut out), 100_000);
    assert_eq!(out, data);
    assert_eq!(out.capacity(), 100_000);
    let mut src = StreamSource::new(OneByte(&data), &data[..10]);
    let mut out = Vec::new();
    assert_eq!(src.read_into(0, 200_000, &mut out), 100_010);
    assert_eq!(out[..10], data[..10]);
    assert_eq!(out[10..], data[..]);
}

/// A long single-line quoted string (> `8 * STREAM_BLOCK`) releases completed 128-byte pieces in
/// `StreamSource` so the stream buffer stays bounded instead of holding the full run alongside
/// `out`, while matching `SliceSource` byte-for-byte on output, `bad_char` pieces, skipped runs,
/// and unended strings.
#[test]
fn long_single_line_visible_string_stays_bounded_and_matches_slice() {
    let body_len = 10 * STREAM_BLOCK + 37;
    let mut text = Vec::with_capacity(body_len + 8);
    text.push(b'"');
    for idx in 0..body_len {
        let b = match idx {
            130 | 65_530 | 65_540 | 524_290 => 0x01,
            _ => b'A' + (idx % 26) as u8,
        };
        text.push(b);
    }
    text.push(b'"');

    for skip in [false, true] {
        let mut slice_out = Vec::new();
        let mut slice_posts = Vec::new();
        let mut rs = Reader::new(SliceSource::new(&text), "text", "test");
        rs.visible_string(&mut slice_out, true, skip, &mut |_, line, byte, piece| {
            slice_posts.push((line, byte, piece.to_vec()));
        })
        .unwrap();

        let mut stream_out = Vec::new();
        let mut stream_posts = Vec::new();
        let mut r = Reader::new(StreamSource::new(text.as_slice(), &[]), "text", "test");
        r.visible_string(&mut stream_out, true, skip, &mut |_, line, byte, piece| {
            stream_posts.push((line, byte, piece.to_vec()));
        })
        .unwrap();

        assert_eq!(stream_out, slice_out);
        assert_eq!(stream_posts, slice_posts);
        assert!(
            r.src.buffered() < 8 * STREAM_BLOCK,
            "StreamSource buffered {} bytes for a {} byte single-line string",
            r.src.buffered(),
            body_len
        );
    }

    // Unended long single-line string: both sources report identical posts, error, and position.
    let unended = &text[..text.len() - 1];
    let mut slice_out = Vec::new();
    let mut slice_posts = Vec::new();
    let mut rs = Reader::new(SliceSource::new(unended), "text", "test");
    let err_slice = rs
        .visible_string(&mut slice_out, true, false, &mut |_, line, byte, piece| {
            slice_posts.push((line, byte, piece.to_vec()));
        })
        .unwrap_err();

    let mut stream_out = Vec::new();
    let mut stream_posts = Vec::new();
    let mut r = Reader::new(StreamSource::new(unended, &[]), "text", "test");
    let err_stream = r
        .visible_string(&mut stream_out, true, false, &mut |_, line, byte, piece| {
            stream_posts.push((line, byte, piece.to_vec()));
        })
        .unwrap_err();

    assert_eq!(stream_out, slice_out);
    assert_eq!(stream_posts, slice_posts);
    assert_eq!(err_stream, err_slice);
    assert_eq!(r.pos, rs.pos);
}

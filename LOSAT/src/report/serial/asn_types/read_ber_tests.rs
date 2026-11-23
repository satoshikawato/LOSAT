//! Unit tests of the generic reader over ASN.1 binary (`super::read::read_ber_document`,
//! `ObjectStream`): a text file and its binary twin written by NCBI's `datatool` give the same
//! values, the objects of one stream are each read by a stream of their own, and the errors and
//! posts are NCBI BLAST+ 2.17.0's (makeblastdb probes of line RG D3c and D3c-3a; the files are
//! `tests/fixtures/asn1`, whose README names their origin).

use super::read::{read_ber_document, read_document, ObjectStream, Post, ReadError};
use super::table::named;
use super::Value;
use crate::ncbi_exception::NCBI_BUILD_ROOT;
use crate::report::serial::asn_input::{Input, Names};
use crate::report::serial::asn_text_reader::{SliceSource, StreamSource};
use crate::report::serial::encoder::StreamPosition;

fn fixture(name: &str) -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/asn1")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// What a post says apart from where it was posted (`line N` in text, `byte N` in binary).
fn post_text(post: &Post) -> String {
    match post {
        Post::BadChar(replaced) => format!("BadChar({}, {:#x})", replaced.trace, replaced.byte),
        other => format!("{other:?}"),
    }
}

// The main test of the binary reader (DESIGN_D3c3 §4): the reader's frames, members, numbers,
// strings and hooks are those of the text reader, so a document and its twin give equal values.
#[test]
fn text_and_binary_twins_read_the_same_values() {
    for (name, root) in [
        ("bioseq", named::BIOSEQ),
        ("e1", named::SEQ_ENTRY),
        ("fr", named::SEQ_ENTRY),
        ("nucprot", named::SEQ_ENTRY),
        ("ids", named::BIOSEQ_SET),
        ("titles", named::BIOSEQ_SET),
        ("ptitles", named::BIOSEQ_SET),
        ("seqdata", named::BIOSEQ_SET),
    ] {
        let (text, text_posts) = read_document(&fixture(&format!("{name}.asn")), root, Names::TEST);
        let (binary, binary_posts) =
            read_ber_document(&fixture(&format!("{name}.asnb")), root, Names::TEST);
        let text = text.unwrap_or_else(|e| panic!("{name}.asn: {e:?}"));
        let binary = binary.unwrap_or_else(|e| panic!("{name}.asnb: {e:?}"));
        assert_eq!(text, binary, "{name}");
        assert_eq!(
            text_posts.iter().map(post_text).collect::<Vec<_>>(),
            binary_posts.iter().map(post_text).collect::<Vec<_>>(),
            "{name}"
        );
    }
}

/// The objects of `stream`, of the root types `roots` in turn, until no object follows
/// (`EndOfData`).
fn objects<I: Input>(mut stream: ObjectStream<I>, roots: &[u16]) -> Vec<Value> {
    let mut values = Vec::new();
    while !stream.end_of_data() {
        let (value, posts) = stream.read_object(roots[values.len()]);
        assert!(posts.is_empty(), "{posts:?}");
        values.push(value.unwrap_or_else(|e| panic!("object {}: {e:?}", values.len())));
    }
    values
}

// Objects one after another (makeblastdb's `-input_type asn1_txt` / `asn1_bin` of a file with
// several): each from the byte after the last one the object before took, from a stream read in
// blocks too; `EndOfData` after the last.
#[test]
fn the_objects_of_one_stream_are_read_one_after_another() {
    let roots = [named::SEQ_ENTRY, named::SEQ_ENTRY, named::BIOSEQ_SET];
    let text = fixture("multi.asn");
    let binary = fixture("multi.asnb");
    let from_text = objects(
        ObjectStream::text(SliceSource::new(&text), Names::TEST),
        &roots,
    );
    let from_binary = objects(
        ObjectStream::binary(SliceSource::new(&binary), Names::TEST),
        &roots,
    );
    assert_eq!(from_text.len(), 3);
    assert_eq!(from_text, from_binary);
    let streamed = objects(
        ObjectStream::binary(StreamSource::new(&binary[..], &[]), Names::TEST),
        &roots,
    );
    assert_eq!(streamed, from_binary);
    let streamed = objects(
        ObjectStream::text(StreamSource::new(&text[..], &[]), Names::TEST),
        &roots,
    );
    assert_eq!(streamed, from_text);
}

// Each object is read by a new stream (`s_ReadObject`): the byte count of a binary object and
// its data limits start again at its first byte, so an error in a later object names the place
// it names in that object alone; the text's line count starts again at 1.
#[test]
fn a_later_object_counts_its_own_bytes_and_lines() {
    let e1 = fixture("e1.asnb");
    let cut = &e1[..e1.len() - 3];
    let alone = read_ber_document(cut, named::SEQ_ENTRY, Names::TEST)
        .0
        .unwrap_err();
    let both = [&e1[..], cut].concat();
    let mut stream = ObjectStream::binary(SliceSource::new(&both), Names::TEST);
    assert!(stream.read_object(named::SEQ_ENTRY).0.is_ok());
    assert_eq!(
        stream.input().position(),
        StreamPosition::Byte(e1.len() as u64)
    );
    assert_eq!(stream.read_object(named::SEQ_ENTRY).0.unwrap_err(), alone);

    let first = fixture("e1.asn");
    let first = first.trim_ascii_end();
    let truncated = fixture("trunc.asn");
    let alone = read_document(&truncated, named::SEQ_ENTRY, Names::TEST)
        .0
        .unwrap_err();
    let text = [first, &truncated].concat();
    let mut stream = ObjectStream::text(SliceSource::new(&text), Names::TEST);
    assert!(stream.read_object(named::SEQ_ENTRY).0.is_ok());
    assert_eq!(stream.read_object(named::SEQ_ENTRY).0.unwrap_err(), alone);
}

/// A string longer than a block of the stream comes from the input straight into the value
/// (`StreamSource::read_into`): an Object-id `str` of 102,400 bytes (`a1 80 1a 83 01 90 00`, the
/// bytes, then the end of the CHOICE) read through a source that reads one byte per fill is the
/// value read from memory.
#[test]
fn a_long_string_comes_from_the_stream_as_from_memory() {
    use crate::report::serial::asn_text_reader::tests::one_byte_per_fill;
    let mut bytes = vec![0xa1, 0x80, 0x1a, 0x83, 0x01, 0x90, 0x00];
    bytes.extend((0..0x19000u32).map(|i| b'a' + (i % 26) as u8));
    bytes.extend([0, 0]);
    let from_memory =
        ObjectStream::binary(SliceSource::new(&bytes), Names::TEST).read_object(named::OBJECT_ID);
    let from_stream =
        ObjectStream::binary(one_byte_per_fill(&bytes), Names::TEST).read_object(named::OBJECT_ID);
    assert!(from_memory.1.is_empty(), "{:?}", from_memory.1);
    match &from_memory.0 {
        Ok(Value::Choice(1, value)) => {
            assert!(matches!(&**value, Value::Str(text) if text[..] == bytes[7..bytes.len() - 2]))
        }
        other => panic!("{:?}", other.as_ref().err()),
    }
    assert_eq!(from_stream, from_memory);
}

/// NCBI's report of an exception (`what()`, NCBI's text before LOSAT's words) of the official
/// build: one place, the message and the frames.
fn ncbi_report(file: &str, line: u32, kind: &str, function: &str, text: &str) -> String {
    format!(
        "NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/{file}\", line {line}: Error: ({kind}) \
         {function}() - {text}\n"
    )
}

fn serial_report(error: ReadError) -> String {
    match error {
        ReadError::Serial(e) => String::from_utf8(e.ncbi_what()).unwrap(),
        other => panic!("{other:?}"),
    }
}

// NCBI BLAST+ 2.17.0 makeblastdb (`D/d3/c/probe/asn2/out/<case>/stderr`, the conda build's tree
// replaced by the official one): `b_bioseq` (a binary Bioseq read as Bioseq-set: its `SET OF
// Seq-id` where the `Object-id` CHOICE of `Bioseq-set.id` has a context tag), `b_txt` (ASN.1 text
// read as binary), `b_fasta`, `t_bin` (binary read as text), `t_bioseq` (a Bioseq document read as
// Bioseq-set), `t_trunc` (a text cut inside its object: the end of the text in `SkipWhiteSpace`).
#[test]
fn errors_are_ncbis() {
    let unexpected = |text: &str| {
        ncbi_report(
            "src/serial/objistrasnb.cpp",
            177,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStreamAsnBinary::UnexpectedTagClassByte",
            text,
        )
    };
    let read = read_ber_document(&fixture("bioseq.asnb"), named::BIOSEQ_SET, Names::TEST);
    assert_eq!(
        serial_report(read.0.unwrap_err()),
        unexpected(
            "byte 4: unexpected tag: constructed/Set (49), should be \
             contextspecific/constructed/0 (160) ( at Bioseq-set.id)"
        )
    );
    let read = read_ber_document(&fixture("e1.asn"), named::SEQ_ENTRY, Names::TEST);
    assert_eq!(
        serial_report(read.0.unwrap_err()),
        unexpected(
            "byte 0: unexpected tag: application/19 (83), should be \
             contextspecific/constructed/0 (160) ( at Seq-entry)"
        )
    );
    let read = read_ber_document(&fixture("f.fa"), named::SEQ_ENTRY, Names::TEST);
    assert_eq!(
        serial_report(read.0.unwrap_err()),
        unexpected(
            "byte 0: unexpected tag: constructed/ObjectReference (62), should be \
             contextspecific/constructed/0 (160) ( at Seq-entry)"
        )
    );
    let read = read_document(&fixture("e1.asnb"), named::SEQ_ENTRY, Names::TEST);
    assert_eq!(
        serial_report(read.0.unwrap_err()),
        ncbi_report(
            "src/serial/objistrasn.cpp",
            719,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStreamAsn::ReadFileHeader",
            "line 1: '::=' expected ( at Seq-entry)"
        )
    );
    let read = read_document(&fixture("bioseq.asn"), named::BIOSEQ_SET, Names::TEST);
    assert_eq!(
        serial_report(read.0.unwrap_err()),
        ncbi_report(
            "src/serial/objistr.cpp",
            896,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStream::SkipFileHeader",
            "line 1: incompatible type Bioseq<>Bioseq-set ( at Bioseq-set)"
        )
    );
    // The end of the text quotes the `CEofException` of the stream buffer; LOSAT's report of it is
    // NCBI's in LOSAT's words.
    let read = read_document(&fixture("trunc.asn"), named::SEQ_ENTRY, Names::TEST);
    let ReadError::Serial(e) = read.0.unwrap_err() else {
        panic!("a CSerialException")
    };
    let nested = ncbi_report(
        "src/util/strbuffer.cpp",
        435,
        "CEofException::eEof",
        "ncbi::CIStreamBuffer::FillBuffer",
        "end of file",
    );
    let ncbi = ncbi_report(
        "src/serial/objistrasn.cpp",
        239,
        "CSerialException::eEOF",
        "ncbi::CObjectIStreamAsn::SkipWhiteSpace",
        &format!("line 5: {nested} ( at Seq-entry.seq.inst)"),
    );
    assert_eq!(
        String::from_utf8(e.what()).unwrap(),
        String::from_utf8(crate::wording::losat_bytes(ncbi.as_bytes()).into_owned()).unwrap()
    );
}

/// What NCBI BLAST+ 2.17.0 makeblastdb did with a probe of the binary reader (`ber/<case>.asnb`).
enum Probe {
    /// Read, with the value of `ber/base.asnb` (its encoding changed, not its value).
    Base,
    /// Read, with its own value (a byte that the reader keeps or replaces).
    Read,
    /// The `CSerialException` at `file:line` (`kind`, `function`) with `text`.
    Serial(&'static str, u32, &'static str, &'static str, &'static str),
}

/// The probes of the binary reader's rules (`D/d3/c3a/probe/ber`, DESIGN_D3c3 §2.3), with the
/// lines NCBI's makeblastdb posted (its standard error; `[makeblastdb]` is the program's prefix).
const BER_PROBES: &[(&str, Probe, &str)] = &[
    ("base", Probe::Base, ""),
    // Definite lengths in the short form, and a primitive's length in the long form, are read.
    ("def", Probe::Base, ""),
    ("def_member", Probe::Base, ""),
    ("def_inst", Probe::Base, ""),
    ("def_title_variant", Probe::Base, ""),
    ("def_list_only", Probe::Base, ""),
    ("len81", Probe::Base, ""),
    // Long tags: a member's (`bf 02`) and an application tag for INTEGER (`42`).
    ("tag_long", Probe::Base, ""),
    ("tag_app_int", Probe::Base, ""),
    // A constructed value's length in the long form: the data limit is as many octets short as
    // the length has after its first (objistrasnb.inl:296).
    (
        "def_outer",
        Probe::Serial(
            "src/serial/objistrasnb.cpp",
            197,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStreamAsnBinary::UnexpectedContinuation",
            "byte 187: EndOfContentsByte expected ( at Seq-entry.seq)",
        ),
        "",
    ),
    (
        "len81_member",
        Probe::Serial(
            "src/serial/objistrasnb.cpp",
            197,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStreamAsnBinary::UnexpectedContinuation",
            "byte 158: EndOfContentsByte expected ( at Seq-entry.seq.inst.length)",
        ),
        "",
    ),
    // The end of the input inside the root CHOICE's variant (`HandleEOF` at depth 3).
    (
        "cut_1",
        Probe::Serial(
            "src/serial/objistr.cpp",
            606,
            "CSerialException::eEOF",
            "ncbi::CObjectIStream::HandleEOF",
            "byte 1:  Frame type= eFrameChoiceVariant ( at Seq-entry)",
        ),
        "",
    ),
    (
        "cut_last",
        Probe::Serial(
            "src/serial/objistr.cpp",
            606,
            "CSerialException::eEOF",
            "ncbi::CObjectIStream::HandleEOF",
            "byte 188:  Frame type= eFrameChoiceVariant, Member name= seq ( at Seq-entry)",
        ),
        "",
    ),
    (
        "member9",
        Probe::Serial(
            "src/serial/objistrasnb.cpp",
            1066,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStreamAsnBinary::UnexpectedMember",
            "byte 187: unexpected member: [9], should be one of: id[0] descr[1] inst[2] annot[3]  \
             ( at Seq-entry.seq.inst)",
        ),
        "",
    ),
    (
        "order",
        Probe::Serial(
            "src/serial/objistrasnb.cpp",
            1066,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStreamAsnBinary::UnexpectedMember",
            "byte 77: unexpected member: [1], should be one of: id[0] descr[1] inst[2] annot[3]  \
             ( at Seq-entry.seq.inst)",
        ),
        "",
    ),
    (
        "missing_inst",
        Probe::Serial(
            "src/serial/objistr.cpp",
            1037,
            "CSerialException::eFormatError",
            "ncbi::CObjectIStream::ExpectedMember",
            "byte 134: member inst expected ( at Seq-entry.seq)",
        ),
        "",
    ),
    // A UTF8String tag for a VisibleString: read as it is, with NCBI's warning; a VisibleString
    // byte outside 0x20..0x7E: `#`, posted at the string's end; a StringStore keeps its bytes.
    (
        "title_utf8tag",
        Probe::Read,
        "Warning: (803.10) [makeblastdb] CObjectIStreamAsnBinary: UTF8String data for VisibleString \
         member Seq-entry.seq.descr.E.title, ASN.1 specification may need an update\n",
    ),
    (
        "title_ff",
        Probe::Read,
        "Error: (803.7) [makeblastdb] Seq-entry.seq.descr.E.title\nBad char [0xFF] in string at \
         byte 40\nmy \u{fffd}seq.\n",
    ),
    (
        "title_c3a9",
        Probe::Read,
        "Error: (803.7) [makeblastdb] Seq-entry.seq.descr.E.title\nBad char [0xC3] in string at \
         byte 40\nmy s\u{e9}q.\nError: (803.7) [makeblastdb] Seq-entry.seq.descr.E.title\n\
         Bad char [0xA9] in string at byte 40\nmy s#\u{fffd}q.\n",
    ),
    ("store_ff", Probe::Read, ""),
];

// NCBI BLAST+ 2.17.0 makeblastdb `-input_type asn1_bin` of each probe (`D/d3/c3a/probe/ber/<case>`:
// its exit status, standard error, and `blastdbcmd` dump of the sequence read). The posts are
// compared as text with each byte that is not UTF-8 as U+FFFD (the probes' bytes are those of the
// string).
#[test]
fn the_binary_probes_are_read_as_ncbi_reads_them() {
    let base = read_ber_document(&fixture("ber/base.asnb"), named::SEQ_ENTRY, Names::TEST)
        .0
        .unwrap();
    for (case, probe, posted) in BER_PROBES {
        let (value, posts) = read_ber_document(
            &fixture(&format!("ber/{case}.asnb")),
            named::SEQ_ENTRY,
            Names::TEST,
        );
        let lines: Vec<u8> = posts
            .iter()
            .flat_map(|post| post.line(Some("makeblastdb")))
            .collect();
        assert_eq!(String::from_utf8_lossy(&lines), *posted, "{case}");
        match probe {
            Probe::Base => assert_eq!(value.unwrap(), base, "{case}"),
            Probe::Read => assert!(value.is_ok(), "{case}: {value:?}"),
            Probe::Serial(file, line, kind, function, text) => assert_eq!(
                serial_report(value.unwrap_err()),
                ncbi_report(file, *line, kind, function, text),
                "{case}"
            ),
        }
    }
}

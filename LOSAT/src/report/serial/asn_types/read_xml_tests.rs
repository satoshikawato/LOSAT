//! Unit tests of the XML front end (`read_xml`): what NCBI's `CObjectIStreamXml` reads of the
//! table's types in the ASN.1 style (`CObjectIStream::Open(eSerial_Xml, f)`, blast_formatter's
//! archive) and where it stops, with its messages (objistrxml.cpp, quoted at each step of
//! `read_xml`), on small values of `Dbtag`, `Seq-interval`, `Seq-loc` and `Seq-data` written by
//! NCBI's XML writer of the archive (`XmlEncoder::ncbi_object`).

use super::read::{read_object, Post, Read, ReadError};
use super::read_xml::{Xml, XmlStyle};
use super::table::named;
use super::write::write_document;
use super::{read_document, Value};
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::asn_text_reader::SliceSource;
use crate::report::serial::encoder::{Encoder, ReplacedChar, StreamPosition};
use crate::report::serial::xml::XmlEncoder;

/// What the XML reader reads of `bytes` as the root type `root`, in the ASN.1 style.
fn read(bytes: &[u8], root: u16) -> Read<Value> {
    read_object(&mut Xml::new(SliceSource::new(bytes), XmlStyle::Asn), root)
}

/// The ASN.1 text of a value of the root type `root`.
fn text_of(value: &Value, root: u16) -> String {
    let mut out = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut out);
    write_document(&mut enc, root, value);
    assert_eq!(enc.error(), None);
    String::from_utf8(out).expect("ASCII")
}

/// The XML that NCBI's object stream writes of a text's value (`os << MSerial_Xml << object`).
fn xml(root: u16, module: &'static str, text: &str) -> String {
    let (value, posts) = read_document(text.as_bytes(), root);
    assert!(posts.is_empty());
    let value = value.expect("a value");
    let mut out = Vec::new();
    let mut enc = XmlEncoder::ncbi_object(&mut out, module);
    write_document(&mut enc, root, &value);
    assert_eq!(enc.error(), None);
    String::from_utf8(out).expect("ASCII")
}

/// The text of what the XML reader reads of `bytes`, or its error; no message.
fn outcome(bytes: &str, root: u16) -> Result<String, ReadError> {
    let (value, posts) = read(bytes.as_bytes(), root);
    assert!(posts.is_empty(), "{posts:?}");
    value.map(|value| text_of(&value, root))
}

fn serial(message: &str) -> Result<String, ReadError> {
    Err(ReadError::Serial(message.into()))
}

/// `text` with the first `from` replaced by `to`.
fn replace(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from} in {text}");
    text.replacen(from, to, 1)
}

const DBTAG: &str = "Dbtag ::= {\n  db \"abc\",\n  tag id 5\n}\n";
const INTERVAL: &str =
    "Seq-interval ::= {\n  from 4294967295,\n  to 2147483648,\n  strand other,\n  \
                        id local str \"q\"\n}\n";
const LOC: &str =
    "Seq-loc ::= packed-int {\n  {\n    from 0,\n    to 9,\n    id gi 5\n  },\n  {\n    \
                   from 20,\n    to 29,\n    strand minus,\n    id local id 7\n  }\n}\n";
const DATA: &str = "Seq-data ::= ncbi2na 'E4'H\n";

#[test]
fn values_written_by_ncbis_writer_read_back() {
    for (root, module, text) in [
        (named::DBTAG, "NCBI-General", DBTAG),
        (named::SEQ_INTERVAL, "NCBI-Seqloc", INTERVAL),
        (named::SEQ_LOC, "NCBI-Seqloc", LOC),
        (named::SEQ_DATA, "NCBI-Sequence", DATA),
    ] {
        let document = xml(root, module, text);
        assert_eq!(outcome(&document, root), Ok(text.into()), "{document}");
    }
    // The tags of the types and members, the list's elements of a named type (no `_E`), the
    // ENUMERATED as an attribute of its named type's tag, an OCTET STRING in hexadecimal.
    let interval = xml(named::SEQ_INTERVAL, "NCBI-Seqloc", INTERVAL);
    assert!(interval.contains(
        "<Seq-interval_strand>\n    <Na-strand value=\"other\"/>\n  </Seq-interval_strand>"
    ));
    let loc = xml(named::SEQ_LOC, "NCBI-Seqloc", LOC);
    assert!(
        loc.contains("<Packed-seqint>\n      <Seq-interval>"),
        "{loc}"
    );
    assert!(xml(named::SEQ_DATA, "NCBI-Sequence", DATA).contains("<NCBI2na>E4</NCBI2na>"));
}

/// White space, comments and processing instructions between tags, attributes NCBI does not
/// know, a self-closed tag of a string (the empty string) and the `value` attribute among others.
/// NCBI's blast_formatter read a processing instruction of the forms `<?pi?>` and `<?pi a="1"?>`
/// between tags (probe D3, 2026-11-20); a word without `=` in it (`<?pi x?>`) ended its input
/// as a serial error (`SkipQDecl` reads the words as attributes).
#[test]
fn layout_attributes_and_empty_tags() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    let spaced = replace(
        &replace(
            &document,
            "<Dbtag_tag>",
            "<!-- a comment -->\n<?pi a=\"1\"?><?pi?><Dbtag_tag  >",
        ),
        "<Object-id>",
        "<Object-id x=\"1\" y='2'>",
    );
    assert_eq!(outcome(&spaced, named::DBTAG), Ok(DBTAG.into()));
    let word = replace(&document, "<Dbtag_tag>", "<?pi x?><Dbtag_tag>");
    assert_eq!(outcome(&word, named::DBTAG), serial("'=' expected"));
    let empty = replace(&document, "<Dbtag_db>abc</Dbtag_db>", "<Dbtag_db/>");
    assert_eq!(
        outcome(&empty, named::DBTAG),
        Ok(DBTAG.replace("\"abc\"", "\"\""))
    );
    let interval = xml(named::SEQ_INTERVAL, "NCBI-Seqloc", INTERVAL);
    let attributes = replace(
        &interval,
        "<Na-strand value=\"other\"/>",
        "<Na-strand a=\"b\" value=\" other \"/>",
    );
    assert_eq!(
        outcome(&attributes, named::SEQ_INTERVAL),
        Ok(INTERVAL.into())
    );
}

/// The messages of NCBI's stream for tags that are not the ones the stack names.
#[test]
fn tags_are_checked_as_ncbi_checks_them() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    // The root's tag names another type: the walk's check of the header.
    assert_eq!(
        outcome(&document.replace("Dbtag", "Seq-id"), named::DBTAG),
        serial("incompatible type Seq-id<>Dbtag")
    );
    // A member that is not one of the class's (`UnexpectedMember`).
    assert_eq!(
        outcome(&document.replace("Dbtag_db", "Dbtag_xx"), named::DBTAG),
        serial("\"xx\": unexpected member, should be one of: \"db\" \"tag\" ")
    );
    // A member's tag not prefixed by its class (`SkipTagName`).
    assert_eq!(
        outcome(&document.replace("Dbtag_db", "Other_db"), named::DBTAG),
        serial("invalid tag name: Other_db")
    );
    // The class's name without `_` (`SkipStackTagName(tag, level, '_')`).
    assert_eq!(
        outcome(&document.replace("Dbtag_db", "Dbtagdb"), named::DBTAG),
        serial("invalid tag name: db")
    );
    // A closing tag that is longer than the member's (`CloseStackTag`: the tag and the rest).
    assert_eq!(
        outcome(
            &replace(&document, "</Dbtag_db>", "</Dbtag_dbx>"),
            named::DBTAG
        ),
        serial("unexpected tag: Dbtag_dbxx")
    );
    // The closing tag of a class (`CloseTag`).
    assert_eq!(
        outcome(&replace(&document, "</Dbtag>", "</Dbtagx>"), named::DBTAG),
        serial("tag 'Dbtag' expected: Dbtagx")
    );
    // A mandatory member missing (the walk's `ExpectedMember`).
    let start = document.find("  <Dbtag_tag>").expect("the member");
    let end = document.find("</Dbtag_tag>\n").expect("its end") + "</Dbtag_tag>\n".len();
    let missing = format!("{}{}", &document[..start], &document[end..]);
    assert_eq!(
        outcome(&missing, named::DBTAG),
        serial("member tag expected")
    );
    // A CHOICE without a variant (`ReadChoiceSimple`).
    let start = document.find("<Object-id_id>").expect("the variant");
    let end = document.find("</Object-id_id>").expect("its end") + "</Object-id_id>".len();
    let no_variant = format!("{}{}", &document[..start], &document[end..]);
    assert_eq!(
        outcome(&no_variant, named::DBTAG),
        serial("choice variant id expected")
    );
}

/// ENUMERATED values by the `value` attribute (`ReadEnum`, `FindValue`), numbers by
/// `CIStreamBuffer`'s readers (their `CUtilException` with the line).
#[test]
fn enumerated_values_and_numbers() {
    let interval = xml(named::SEQ_INTERVAL, "NCBI-Seqloc", INTERVAL);
    assert_eq!(
        outcome(
            &replace(&interval, "value=\"other\"", "value=\"up\""),
            named::SEQ_INTERVAL
        ),
        serial("invalid value of enumerated type: up")
    );
    assert_eq!(
        outcome(
            &replace(&interval, "value=\"other\"", "valu=\"other\""),
            named::SEQ_INTERVAL
        ),
        serial("attribute 'value' is missing")
    );
    assert_eq!(
        outcome(
            &replace(
                &interval,
                "<Na-strand value=\"other\"/>",
                "<Na-strand>255</Na-strand>"
            ),
            named::SEQ_INTERVAL
        ),
        serial("attribute 'value' expected")
    );
    let line = |marker: &str| {
        interval[..interval.find(marker).expect("marker")]
            .matches('\n')
            .count()
            + 1
    };
    let number = |from: &str| match read(
        replace(&interval, ">4294967295<", from).as_bytes(),
        named::SEQ_INTERVAL,
    )
    .0
    {
        Err(ReadError::Exception(e)) => (e.site().line, e.get_msg().to_string()),
        other => panic!("{other:?}"),
    };
    let at = line("<Seq-interval_from>");
    assert_eq!(number(">x<"), (574, format!("bad number in line {at}")));
    assert_eq!(
        number(">4294967296<"),
        (583, format!("number overflow in line {at}"))
    );
    // Uint4: no sign but `+`.
    assert_eq!(number(">-1<"), (574, format!("bad number in line {at}")));
    assert_eq!(
        outcome(
            &replace(&interval, ">4294967295<", "> +4294967295<"),
            named::SEQ_INTERVAL
        ),
        Ok(INTERVAL.into())
    );
}

/// The end of the input: in the root's own frame the end of the data (`CEofException`), below
/// it the walk's `CSerialException` (`HandleEOF`).
#[test]
fn the_end_of_the_input_is_an_error_below_the_root() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    assert_eq!(read(b"", named::DBTAG).0, Err(ReadError::Eof));
    let cut = document.find("abc").expect("the string");
    match read(&document.as_bytes()[..cut], named::DBTAG).0 {
        Err(ReadError::Serial(message)) => {
            assert!(
                message.starts_with("unexpected end of the data"),
                "{message}"
            )
        }
        other => panic!("{other:?}"),
    }
}

/// A control character in character data: `#`, posted with the line it was read on and no
/// string (`x_VerifyChar`, `ReplaceVisibleChar` with `kEmptyStr`); NCBI's standard error of
/// `blast_formatter` is the post's one line (probe D3: 0x01 posted, 0x7F not).
#[test]
fn control_characters_are_replaced_and_posted_with_their_line() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    let bad = replace(&document, ">abc<", ">\u{1}b\u{7}<");
    let (value, posts) = read(bad.as_bytes(), named::DBTAG);
    assert_eq!(
        text_of(&value.expect("a value"), named::DBTAG),
        DBTAG.replace("abc", "#b#")
    );
    let line = document[..document.find(">abc<").expect("the string")]
        .matches('\n')
        .count()
        + 1;
    let replaced = |byte: u8| {
        Post::BadChar(ReplacedChar {
            trace: "Dbtag.db".into(),
            byte,
            position: StreamPosition::Line(line as u64),
            string: Vec::new().into(),
        })
    };
    assert_eq!(posts, vec![replaced(0x01), replaced(0x07)]);
    assert_eq!(
        posts[0].line(None),
        format!("Error: (803.7) Dbtag.db\nBad char [0x1] in string at line {line}\n").into_bytes()
    );
}

/// The standard-XML style (`SetEnforcedStdXml(true)`, the search strategy of line RG-F2) is not
/// read yet: the stream refuses it at the header.
#[test]
fn the_standard_style_is_refused_until_it_is_read() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    let read = read_object(
        &mut Xml::new(SliceSource::new(document.as_bytes()), XmlStyle::Std),
        named::DBTAG,
    );
    assert!(
        matches!(read.0, Err(ReadError::Unsupported(_))),
        "{:?}",
        read.0
    );
}

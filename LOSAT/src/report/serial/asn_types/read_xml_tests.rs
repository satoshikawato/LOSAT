//! Unit tests of the XML stream (`read_xml`) under the object stack of `read`: what NCBI's
//! `CObjectIStreamXml` reads of the table's types in the ASN.1 style (`CObjectIStream::Open` with
//! `eSerial_Xml`: psiblast's and tblastn's `-in_pssm`, blast_formatter's archive) and where it
//! stops, with its messages (objistrxml.cpp, cited at each step of `read_xml`), on small values
//! written by LOSAT's port of NCBI's XML writer (`XmlEncoder::ncbi_object`).

use super::read::{read_delayed, ObjectStream, Post, Read, ReadError};
use super::read_xml::{find_empty, XmlStyle};
use super::table::named;
use super::write::write_document;
use super::{read_document, type_info, Encoding, Kind, Value};
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::asn_text_reader::SliceSource;
use crate::report::serial::encoder::{Encoder, ReplacedChar, SerialError, StreamPosition, Traced};
use crate::report::serial::xml::XmlEncoder;

/// What the XML stream reads of `bytes` as the root type `root`, in the ASN.1 style.
fn read(bytes: &[u8], root: u16) -> Read<Value> {
    ObjectStream::xml(SliceSource::new(bytes), XmlStyle::Asn, Names::TEST).read_object(root)
}

/// The ASN.1 text of a value of the root type `root`.
fn text_of(value: &Value, root: u16) -> String {
    let mut out = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut out);
    write_document(&mut enc, root, value);
    assert_eq!(enc.error(), None);
    drop(enc);
    String::from_utf8(out).expect("ASCII")
}

/// The XML that NCBI's object stream writes of a text's value (an object written to an output
/// stream with the manipulator `MSerial_Xml`).
fn xml(root: u16, module: &'static str, text: &str) -> String {
    let (value, posts) = read_document(text.as_bytes(), root, Names::TEST);
    assert!(posts.is_empty());
    let value = value.expect("a value");
    let mut out = Vec::new();
    let mut enc = XmlEncoder::ncbi_object(&mut out, module);
    write_document(&mut enc, root, &value);
    assert_eq!(enc.error(), None);
    drop(enc);
    String::from_utf8(out).expect("ASCII")
}

/// The text of what the XML stream reads of `bytes`, or its error: the message of a
/// `CSerialException` (without the stream's position), else the error itself; no post.
fn outcome(bytes: &str, root: u16) -> Result<String, String> {
    let (value, posts) = read(bytes.as_bytes(), root);
    assert!(posts.is_empty(), "{posts:?}");
    value
        .map(|value| text_of(&value, root))
        .map_err(|error| match error {
            ReadError::Serial(e) => {
                let message = String::from_utf8_lossy(e.get_msg()).into_owned();
                match message.split_once(": ") {
                    Some((position, rest)) if position.starts_with("line ") => rest.to_owned(),
                    _ => message,
                }
            }
            other => format!("{other:?}"),
        })
}

fn serial(message: &str) -> Result<String, String> {
    Err(message.into())
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
/// between tags (probe D3 of line RG-F3, 2026-11-20); a word without `=` in it (`<?pi x?>`)
/// ended its input as a serial error (`SkipQDecl` reads the words as attributes).
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
    // The root's tag names another type: the stack's check of the header.
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
    // The class's name without `_` (`SkipStackTagName` with the separator `_`).
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
    // A mandatory member missing (the stack's `ExpectedMember`).
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
/// `CIStreamBuffer`'s readers (their `CUtilException` with the line), one parser of the numbers
/// for the text and XML streams (`asn_text_reader::buffer_number`).
#[test]
fn enumerated_values_and_numbers() {
    let interval = xml(named::SEQ_INTERVAL, "NCBI-Seqloc", INTERVAL);
    // `strand` is OPTIONAL: `ReadWithSetFlagMember` rethrows the error of its value.
    let strand = |bytes: &str, thrown: &str| {
        assert_eq!(
            outcome(bytes, named::SEQ_INTERVAL),
            serial("error while reading strand")
        );
        match read(bytes.as_bytes(), named::SEQ_INTERVAL).0 {
            Err(ReadError::Serial(e)) => {
                let report = String::from_utf8_lossy(&e.ncbi_what()).into_owned();
                assert!(report.contains(thrown), "{report}");
            }
            other => panic!("{other:?}"),
        }
    };
    strand(
        &replace(&interval, "value=\"other\"", "value=\"up\""),
        "invalid value of enumerated type: up",
    );
    strand(
        &replace(
            &interval,
            "<Na-strand value=\"other\"/>",
            "<Na-strand>255</Na-strand>",
        ),
        "attribute 'value' expected",
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
        // The exception as thrown (the oldest of NCBI's report of the chain: each frame left
        // rethrows it).
        Err(ReadError::Exception(e)) => {
            let report = String::from_utf8_lossy(&e.ncbi_what()).into_owned();
            let first = report.lines().nth(1).expect("the throw").to_string();
            let line = first
                .split(", line ")
                .nth(1)
                .and_then(|rest| rest.split(':').next())
                .and_then(|line| line.parse::<u32>().ok())
                .expect("a line");
            let message = first.rsplit("() - ").next().expect("a message").to_string();
            (line, message)
        }
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
/// it the stack's `CSerialException` (`HandleEOF`: the stream's position and the information of
/// the frame it leaves, here the member `db` whose string the input cuts; that information starts
/// with a space, objstack.cpp:399-410).
#[test]
fn the_end_of_the_input_is_an_error_below_the_root() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    assert_eq!(read(b"", named::DBTAG).0, Err(ReadError::Eof));
    let cut = document.find("abc").expect("the string");
    let line = document[..cut].matches('\n').count() + 1;
    match read(&document.as_bytes()[..cut], named::DBTAG).0 {
        Err(ReadError::Serial(e)) => assert_eq!(
            String::from_utf8_lossy(e.get_msg()),
            format!("line {line}:  Frame type= eFrameClassMember, Member name= db")
        ),
        other => panic!("{other:?}"),
    }
}

/// A control character in character data: `#`, posted at once with the stack trace of the
/// moment, the line it was read on and no string (`x_VerifyChar`, `ReplaceVisibleChar` with
/// `kEmptyStr`); NCBI's standard error of `blast_formatter` is the post's one line (probe D3 of
/// line RG-F3: 0x01 posted, 0x7F not). The posts of an attribute that the stream skips come at
/// once too, before an error that follows them (review D of line RG-F3, F4).
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
    let replaced = |trace: &str, byte: u8, line: usize| {
        Post::BadChar(ReplacedChar {
            trace: trace.into(),
            byte,
            position: StreamPosition::Line(line as u64),
            string: Vec::new().into(),
        })
    };
    assert_eq!(
        posts,
        vec![
            replaced("Dbtag.db", 0x01, line),
            replaced("Dbtag.db", 0x07, line)
        ]
    );
    assert_eq!(
        posts[0].line("blast_formatter"),
        format!(
            "Error: (803.7) [blast_formatter] Dbtag.db\nBad char [0x1] in string at line {line}\n"
        )
        .into_bytes()
    );
    // In an attribute skipped at a member's tag, then an unknown member: the post, then the error.
    let tag_line = document[..document.find("<Dbtag_tag>").expect("the member")]
        .matches('\n')
        .count()
        + 1;
    let bad = replace(
        &replace(&document, "<Dbtag_tag>", "<Dbtag_tag a=\"\u{2}\">"),
        "<Object-id_id>",
        "<Object-id_xx>",
    );
    let bad = bad.replacen("</Object-id_id>", "</Object-id_xx>", 1);
    let (value, posts) = read(bad.as_bytes(), named::DBTAG);
    assert!(matches!(value, Err(ReadError::Serial(_))), "{value:?}");
    assert_eq!(posts, vec![replaced("Dbtag.tag", 0x02, tag_line)]);
}

/// Namespace prefixes and declarations (`ReadName`): a prefixed tag is read by the name after its
/// last `:`; `xmlns` and `xmlns:<prefix>` attributes are read with their values; a `nil`
/// attribute whose prefix is not declared or names the schema-instance namespace holds a boolean
/// (else `NStr::StringToBool`'s exception at position 0, through the frames), one of another
/// namespace is skipped (NCBI 2.17.0 probes ns_root, ns_tags, ns_member, nil_true, nil_bad,
/// xsnil_bad, othernil_bad of line RG-F3 agent D6, 2026-11-23).
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:432-524
#[test]
fn namespace_prefixes_and_nil_are_read_as_ncbi_reads_them() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    let prefixed = document
        .replace("<Dbtag", "<x:Dbtag")
        .replace("</Dbtag", "</x:Dbtag")
        .replace("<Object-id", "<x:Object-id")
        .replace("</Object-id", "</x:Object-id");
    let prefixed = replace(&prefixed, "<x:Dbtag>", "<x:Dbtag xmlns:x=\"u\">");
    assert_eq!(outcome(&prefixed, named::DBTAG), Ok(DBTAG.into()));
    let declared = replace(
        &document,
        "<Dbtag>",
        "<Dbtag xmlns=\"http://www.ncbi.nlm.nih.gov\" \
         xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\" xs:schemaLocation=\"a b\">",
    );
    assert_eq!(outcome(&declared, named::DBTAG), Ok(DBTAG.into()));
    let nil = |attribute: &str, root: &str| {
        let text = replace(&document, "<Dbtag_db>", &format!("<Dbtag_db {attribute}>"));
        replace(&text, "<Dbtag>", &format!("<Dbtag{root}>"))
    };
    assert_eq!(
        outcome(&nil("nil=\"true\"", ""), named::DBTAG),
        Ok(DBTAG.into())
    );
    let other = nil("o:nil=\"x\"", " xmlns:o=\"http://example.org\"");
    assert_eq!(outcome(&other, named::DBTAG), Ok(DBTAG.into()));
    let schema = " xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"";
    for bad in [nil("nil=\"x\"", ""), nil("xs:nil=\"x\"", schema)] {
        match read(bad.as_bytes(), named::DBTAG).0 {
            Err(ReadError::Exception(e)) => {
                // The report in LOSAT's words: NCBI's `(m_Pos = 0)` is `(position 0)`
                // (`crate::wording`).
                let report = String::from_utf8_lossy(&e.ncbi_what()).into_owned();
                assert!(
                    report.contains(
                        "NStr::StringToBool() - String cannot be converted to bool (position 0)"
                    ) && report.contains("Member name= db"),
                    "{report}"
                );
            }
            other => panic!("{other:?}"),
        }
    }
}

/// The standard-XML style (`SetEnforcedStdXml` on, the search strategy) is not read yet: the
/// stream refuses it at the header.
#[test]
fn the_standard_style_is_refused_until_it_is_read() {
    let document = xml(named::DBTAG, "NCBI-General", DBTAG);
    let read = ObjectStream::xml(
        SliceSource::new(document.as_bytes()),
        XmlStyle::Std,
        Names::TEST,
    )
    .read_object(named::DBTAG);
    assert!(
        matches!(read.0, Err(ReadError::Unsupported(_))),
        "{:?}",
        read.0
    );
}

/// `<Seq-inst>` of a protein of 3 residues with its DEFAULT `topology` written.
const INST: &str = "Seq-inst ::= {\n  repr raw,\n  mol aa,\n  length 3,\n  topology circular\n}\n";

/// An ENUMERATED without `value` (`ReadEnum`'s `eMissingValue`) in an OPTIONAL or DEFAULT member
/// leaves the member unset (`ReadWithSetFlagMember`): NCBI 2.17.0's blast_formatter formats
/// `<Seq-inst_strand foo="ss"/>` as the archive without it (probe x2 of review D of line RG-F3).
/// In a mandatory member it is "missing value" (probe mv-m3 of line RG-F3 agent D6: a
/// `Seq-align_type` without `value`, "Invalid input format", rc 1).
// NCBI reference (598d8ae6): c++/src/serial/member.cpp:714-756
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1414-1426
#[test]
fn missing_value_takes_an_optional_or_default_member_away() {
    let inst = xml(named::SEQ_INST, "NCBI-Sequence", INST);
    let without = |from: &str| {
        outcome(
            &replace(&inst, from, &from.replace("value=", "foo=")),
            named::SEQ_INST,
        )
    };
    assert_eq!(
        without("<Seq-inst_topology value=\"circular\"/>"),
        Ok(INST.replace(",\n  topology circular", ""))
    );
    assert_eq!(
        without("<Seq-inst_repr value=\"raw\"/>"),
        serial("missing value repr")
    );
    assert_eq!(
        without("<Seq-inst_mol value=\"aa\"/>"),
        serial("missing value mol")
    );
}

/// After a missing value the stream stands where `ReadEnum` stopped: the ENUMERATED's own tag
/// (of a named type) ended as self-closed, the tags around it not read, so the class reads on
/// from there: the member's closing tag is not read (`CloseStackTag` of a self-closed tag), and
/// the next tag is the error. NCBI 2.17.0's blast_formatter: `<Na-strand foo="plus"/>` as the
/// first or the last element of an archive's `Dense-seg_strands` (OPTIONAL) is "Invalid input
/// format", rc 1 (probes mv-m1, mv-m2 of line RG-F3 agent D6, 2026-11-23).
#[test]
fn after_a_missing_value_the_stream_reads_on_where_it_stopped() {
    let interval = xml(named::SEQ_INTERVAL, "NCBI-Seqloc", INTERVAL);
    let without = replace(&interval, "value=\"other\"", "foo=\"other\"");
    // The member `id` after `strand`: the class ends at `</Seq-interval_strand>`, `id` is missing.
    assert_eq!(
        outcome(&without, named::SEQ_INTERVAL),
        serial("member id expected")
    );
    let loc = xml(named::SEQ_LOC, "NCBI-Seqloc", LOC);
    let without = replace(&loc, "value=\"minus\"", "foo=\"minus\"");
    assert!(outcome(&without, named::SEQ_LOC).is_err());
}

/// A CHOICE element with nothing in it holds the first variant that is an anonymous list
/// (`FindEmpty`), empty, without a tag (`SetNotag`); a self-closed one, whose opening tag
/// `OpenTag` left open, is the "'>' expected" of `EndTag` (`NextTagIsClosing` ends the opening
/// tag, where `/` is no attribute), and one whose variant `FindEmpty` does not find is "choice
/// variant id expected" (NCBI 2.17.0 probes x1 and x1b of review D of line RG-F3: rc 0, rc 1).
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:308-340 (`EndTag`)
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1666-1670 (`NextTagIsClosing`)
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2411-2440
// NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:406-430
#[test]
fn an_empty_choice_element_holds_its_first_list_variant() {
    let empty = "<?xml version=\"1.0\"?>\n<Seqdesc></Seqdesc>\n";
    let (value, posts) = read(empty.as_bytes(), named::SEQDESC);
    assert!(posts.is_empty(), "{posts:?}");
    let modif = read_document(&b"Seqdesc ::= modif { }\n"[..], named::SEQDESC, Names::TEST).0;
    assert_eq!(value, modif);
    assert!(modif.is_ok());
    assert_eq!(
        outcome(
            "<?xml version=\"1.0\"?>\n<Seqdesc>\n  \n</Seqdesc>\n",
            named::SEQDESC
        ),
        outcome(empty, named::SEQDESC)
    );
    assert_eq!(
        outcome("<?xml version=\"1.0\"?>\n<Seqdesc/>\n", named::SEQDESC),
        serial("'>' expected")
    );
    // Its variant is an INTEGER: none.
    assert_eq!(
        outcome(
            "<?xml version=\"1.0\"?>\n<Object-id></Object-id>\n",
            named::OBJECT_ID
        ),
        serial("choice variant id expected")
    );
    let mut found: Vec<String> = super::table::TYPES
        .iter()
        .enumerate()
        .filter_map(|(ty, info)| match info.kind {
            Kind::Choice(variants) => find_empty(variants.members())
                .map(|at| format!("{}.{}", type_path(ty as u16), variants.members()[at].name())),
            _ => None,
        })
        .collect();
    found.sort();
    let mut expected = [
        "Seqdesc.modif",
        "Int-fuzz.alt",
        "User-field.data.strs",
        "Seq-align.segs.dendiag",
        "Auth-list.names.std",
        "Pub-set.pub",
        "SeqTable-multi-data.int",
        "SeqTable-sparse-index.indexes",
        "Seq-annot.data.ftable",
    ];
    expected.sort();
    assert_eq!(found, expected);
}

/// The name of the type `ty`, or of the member that holds an anonymous one (`Parent.member`).
fn type_path(ty: u16) -> String {
    let info = type_info(ty);
    if !info.name().is_empty() {
        return info.name().to_string();
    }
    for (parent, outer) in super::table::TYPES.iter().enumerate() {
        let (Kind::Sequence(members) | Kind::Choice(members)) = outer.kind else {
            continue;
        };
        if let Some(member) = members.members().iter().find(|member| member.ty == ty) {
            return format!("{}.{}", type_path(parent as u16), member.name());
        }
    }
    panic!("type {ty} is held by no member")
}

/// `Auth-list` with its delayed member `names` (`DELAYED`).
const AUTH: &str =
    "Auth-list ::= {\n  names std {\n    {\n      name name {\n        last \"Smith\"\n      }\n    \
     }\n  },\n  affil str \"a\"\n}\n";

/// The delayed member is skipped into its bytes (`ReadLongMember`), never parsed while the
/// object is read: a number past its storage or a control byte inside it is no error and no post
/// (NCBI 2.17.0's blast_formatter, probes x4, x5 of review D of line RG-F3: rc 0, nothing on
/// standard error). A stream of the XML format copies the bytes; another one reads them first
/// under the member's frames (`read_delayed`, `CDelayBuffer::Update`): their replaced bytes are
/// posted, and an error keeps LOSAT from writing the member.
#[test]
fn the_delayed_member_is_kept_unparsed_and_read_on_a_write_of_another_format() {
    let document = xml(named::AUTH_LIST, "NCBI-Biblio", AUTH);
    let (value, posts) = read(document.as_bytes(), named::AUTH_LIST);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.expect("read");
    let Value::Sequence(members) = &value else {
        panic!("{value:?}")
    };
    let Value::Delayed { encoding, bytes } = &members[0].1 else {
        panic!("{members:?}")
    };
    assert_eq!(*encoding, Encoding::Xml);
    let bytes = String::from_utf8(bytes.clone()).expect("ASCII");
    assert!(
        bytes.trim_start().starts_with("<Auth-list_names_std>"),
        "{bytes}"
    );
    assert!(bytes.ends_with("</Auth-list_names_std>"), "{bytes}");
    // Written in ASN.1 text: read again, written by its values.
    assert_eq!(text_of(&value, named::AUTH_LIST), AUTH);
    // Written in XML: copied, the member's closing tag right after the bytes (the last tag that
    // NCBI's writer wrote is the member's opening tag, so `CloseTagStart` puts no end of line).
    // NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1541-1555
    // NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1069-1079
    assert_eq!(
        xml_of(&value, named::AUTH_LIST, "NCBI-Biblio"),
        replace(
            &document,
            "</Auth-list_names_std>\n  </Auth-list_names>",
            "</Auth-list_names_std></Auth-list_names>"
        )
    );
    // A control byte and a number past an Int4 inside it: no error, no post while it is read.
    let ctrl = replace(&document, ">Smith<", ">S\u{1}mith<");
    let (value, posts) = read(ctrl.as_bytes(), named::AUTH_LIST);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.expect("read");
    // A text write reads it first; the read posts the byte (803.7) under its own frames, with the
    // line in the buffer and no string, and the member is written with `#` (NCBI 2.17.0's
    // psiblast -in_pssm of such an XML PSSM with -out_pssm: "Auth-list.names.std.name.str / Bad
    // char [0x1] in string at line 6", probe dl_ctrl of line RG-F3 agent D6, 2026-11-23). The
    // encoder alone reports it as a replaced byte of its own; through `Traced` it is posted.
    let mut out = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut out);
    write_document(&mut enc, named::AUTH_LIST, &value);
    assert_eq!(enc.error(), Some(SerialError::NonVisibleChar(1)));
    drop(enc);
    let mut out = Vec::new();
    let mut enc = Traced::new(AsnTextEncoder::new(&mut out), 0);
    write_document(&mut enc, named::AUTH_LIST, &value);
    assert_eq!(enc.error(), None);
    let posted = enc.take_replaced_chars();
    drop(enc);
    assert_eq!(
        String::from_utf8(out).expect("ASCII"),
        AUTH.replace("Smith", "S#mith")
    );
    let Value::Sequence(members) = &value else {
        panic!("{value:?}")
    };
    let Value::Delayed { bytes, .. } = &members[0].1 else {
        panic!("{members:?}")
    };
    let at = bytes.iter().position(|&b| b == 1).expect("the byte");
    let line = 1 + bytes[..at].iter().filter(|&&b| b == b'\n').count() as u64;
    assert_eq!(
        posted,
        vec![ReplacedChar {
            trace: "Auth-list.names.std.name.name.last".into(),
            byte: 1,
            position: StreamPosition::Line(line),
            string: Vec::new().into(),
        }]
    );
    // The same read alone.
    let (class, member) = super::delayed_member(super::members(named::AUTH_LIST)[0].ty);
    let (read_again, posts) = read_delayed(Encoding::Xml, bytes, class, member);
    assert!(read_again.is_ok(), "{read_again:?}");
    assert_eq!(posts, vec![Post::BadChar(posted[0].clone())]);
    // A number past its storage inside it: read without an error, not written in text (NCBI's
    // exception from its writer, probe dl_ovf: rc 255 after 283 bytes of the PSSM file).
    let ovf = document.replace(
        "<Person-id_name>",
        "<Person-id_dbtag><Dbtag><Dbtag_db>x</Dbtag_db><Dbtag_tag><Object-id>\
         <Object-id_id>99999999999</Object-id_id></Object-id></Dbtag_tag></Dbtag>\
         </Person-id_dbtag><Person-id_name>",
    );
    let ovf = remove_between(&ovf, "<Person-id_name>", "</Person-id_name>");
    let (value, posts) = read(ovf.as_bytes(), named::AUTH_LIST);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.expect("read");
    let mut out = Vec::new();
    let mut enc = Traced::new(AsnTextEncoder::new(&mut out), 0);
    write_document(&mut enc, named::AUTH_LIST, &value);
    assert!(
        matches!(enc.error(), Some(SerialError::Unwritten(_))),
        "{:?}",
        enc.error()
    );
}

/// `text` without the first `from` and everything up to and including the next `to`.
fn remove_between(text: &str, from: &str, to: &str) -> String {
    let start = text.find(from).expect("from");
    let end = start + text[start..].find(to).expect("to") + to.len();
    format!("{}{}", &text[..start], &text[end..])
}

/// The XML of `value` of the root type `root` (`xml`'s writer).
fn xml_of(value: &Value, root: u16, module: &'static str) -> String {
    let mut out = Vec::new();
    let mut enc = XmlEncoder::ncbi_object(&mut out, module);
    write_document(&mut enc, root, value);
    assert_eq!(enc.error(), None);
    drop(enc);
    String::from_utf8(out).expect("ASCII")
}

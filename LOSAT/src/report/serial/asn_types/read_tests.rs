//! Unit tests of the generic reader (`super::read`): one test per rule of NCBI's text reader that it
//! ports (SURVEY_D3 §2 of line RG series F; the NCBI 2.17.0 probes `c01`-`c23` and the edge
//! probes `e01`-`e07` of item D3-1 are named where a test takes their case).

use super::read::{read_document, Post, ReadError};
use super::table::named;
use super::{member_index, Value};
use crate::ncbi_exception::NCBI_BUILD_ROOT;
use crate::seqid::model::{ObjectId, ParsedSeqId};

/// The value of the document `text` of the root type `root`, or the error, with the posts.
fn read(root: u16, text: impl AsRef<[u8]>) -> (Result<Value, ReadError>, Vec<Post>) {
    read_document(text.as_ref(), root)
}

fn value(root: u16, text: &str) -> Value {
    let (value, posts) = read(root, text);
    assert!(posts.is_empty(), "{posts:?}");
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

fn serial(root: u16, text: &str) {
    match read(root, text).0 {
        Err(ReadError::Serial(_)) => {}
        other => panic!("{text}: {other:?}"),
    }
}

/// The report of the exception that stops the reader (`what()`), with the official build's tree.
fn exception(root: u16, text: &str) -> String {
    match read(root, text).0 {
        Err(ReadError::Exception(e)) => e.what(),
        other => panic!("{text}: {other:?}"),
    }
}

/// A line of NCBI's report: `T0 "<tree>/<file>", line <line>: Error: (<kind>) <text>`.
fn report_line(file: &str, line: u32, kind: &str, text: &str) -> String {
    format!("    T0 \"{NCBI_BUILD_ROOT}/{file}\", line {line}: Error: ({kind}) {text}\n")
}

fn util(file: &str, line: u32, text: &str) -> String {
    report_line(file, line, "CUtilException::eWrongData", text)
}

const ASN: &str = "src/serial/objistrasn.cpp";
const ISTR: &str = "src/serial/objistr.cpp";

fn member_frame(name: &str) -> String {
    util(
        ASN,
        1393,
        &format!("ncbi::CObjectIStreamAsn::ReadClassSequential() -  Frame type= eFrameClassMember, Member name= {name}"),
    )
}

fn class_frame(name: &str) -> String {
    util(
        ASN,
        1396,
        &format!("ncbi::CObjectIStreamAsn::ReadClassSequential() -  Frame type= eFrameClass, Object type= {name}"),
    )
}

fn named_frame(name: &str) -> String {
    util(
        ISTR,
        968,
        &format!("ncbi::CObjectIStream::Read() -  Frame type= eFrameNamed, Object type= {name}"),
    )
}

fn seq_inst(body: &str) -> String {
    format!("Seq-inst ::= {{\n  repr raw,\n  mol aa{body}\n}}\n")
}

/// The member `name` of the Seq-inst value.
fn inst_member(value: &Value, name: &str) -> Option<Value> {
    value.member(member_index(named::SEQ_INST, name)).cloned()
}

// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-700 (quoted at `read::get_number`):
// a sign right before the digits, leading zeros, the ranges of Int4 and Uint4 (`-` is not a sign
// of an unsigned), Int8 for `TGi`. Probes c02 (`pct -0`), c10 (`+7`, `007`), c17 (`to +29`, gi
// 5000000000), c18 (`length -1`).
#[test]
fn integers_follow_their_storage() {
    let length = |text: &str| inst_member(&value(named::SEQ_INST, &seq_inst(text)), "length");
    assert_eq!(length(",\n  length +29"), Some(Value::Int(29)));
    assert_eq!(length(",\n  length 007"), Some(Value::Int(7)));
    assert_eq!(
        length(",\n  length 4294967295"),
        Some(Value::Int(4294967295))
    );
    let fuzz = |text: &str| {
        inst_member(
            &value(named::SEQ_INST, &seq_inst(&format!(",\n  fuzz {text}"))),
            "fuzz",
        )
    };
    assert_eq!(
        fuzz("pct -0"),
        Some(Value::Choice(2, Box::new(Value::Int(0))))
    );
    assert_eq!(
        fuzz("pct -2147483648"),
        Some(Value::Choice(2, Box::new(Value::Int(-2147483648))))
    );
    assert!(
        exception(named::SEQ_INST, &seq_inst(",\n  fuzz pct 2147483648"))
            .contains("number overflow in line 4")
    );
    assert!(
        exception(named::SEQ_INST, &seq_inst(",\n  fuzz pct -2147483649"))
            .contains("number overflow in line 4")
    );
    assert!(
        exception(named::SEQ_INST, &seq_inst(",\n  length 4294967296"))
            .contains("number overflow in line 4")
    );
    // `-` before an unsigned, and white space after a sign: `BadNumber`.
    for text in [",\n  length -1", ",\n  length -0", ",\n  fuzz pct - 5"] {
        assert!(
            exception(named::SEQ_INST, &seq_inst(text)).contains("bad number in line 4"),
            "{text}"
        );
    }
    let gi = value(named::SEQ_ID, "Seq-id ::= gi 5000000000");
    assert_eq!(gi, Value::SeqId(ParsedSeqId::Gi(5000000000)));
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704 (quoted at `read::integer`):
// `Int2` is an Int4 that must fit, else a `CSerialException`.
#[test]
fn an_int2_out_of_range_is_a_serial_exception() {
    let table = |values: &str| format!("SeqTable-multi-data ::= int2 {{\n  {values}\n}}\n");
    assert_eq!(
        value(named::SEQTABLE_MULTI_DATA, &table("32767, -32768")),
        Value::Choice(
            member_index(named::SEQTABLE_MULTI_DATA, "int2") as u16,
            Box::new(Value::List(vec![Value::Int(32767), Value::Int(-32768)]))
        )
    );
    serial(named::SEQTABLE_MULTI_DATA, &table("32768"));
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737 (quoted at
// `read::enum_value`): an ENUMERATED by name or by a number with a name (c02 `topology 2`, c06
// `topology 7`), an INTEGER with named values by name or any number (c01 `biomol 8`,
// `completeness 099`); an Int4 overflow in either gets the line of `CEnumeratedTypeInfo::ReadEnum`
// (NCBI 2.17.0 probe e05).
#[test]
fn enumerated_and_named_integers() {
    let topology = |text: &str| {
        inst_member(
            &value(named::SEQ_INST, &seq_inst(&format!(",\n  topology {text}"))),
            "topology",
        )
    };
    assert_eq!(topology("circular"), Some(Value::Enum(2)));
    assert_eq!(topology("2"), Some(Value::Enum(2)));
    assert_eq!(topology("+002"), Some(Value::Enum(2)));
    serial(named::SEQ_INST, &seq_inst(",\n  topology 7"));
    serial(named::SEQ_INST, &seq_inst(",\n  topology round"));
    // An upper-case first letter is not a name: the number reader sees it.
    assert!(
        exception(named::SEQ_INST, &seq_inst(",\n  topology Linear"))
            .contains("bad number in line 4")
    );
    let molinfo = |text: &str| value(named::MOLINFO, &format!("MolInfo ::= {{ {text} }}"));
    let biomol = member_index(named::MOLINFO, "biomol") as u16;
    assert_eq!(
        molinfo("biomol 8"),
        Value::Sequence(vec![(biomol, Value::Enum(8))])
    );
    assert_eq!(
        molinfo("biomol 099"),
        Value::Sequence(vec![(biomol, Value::Enum(99))])
    );
    let report = exception(named::SEQ_INST, &seq_inst(",\n  topology 99999999999"));
    assert_eq!(
        report,
        [
            "NCBI C++ Exception:\n".to_string(),
            util(
                "src/util/strbuffer.cpp",
                583,
                "ncbi::CIStreamBuffer::NumberOverflow() - number overflow in line 4"
            ),
            util(
                "src/serial/enumerated.cpp",
                428,
                "ncbi::CEnumeratedTypeInfo::ReadEnum() - invalid enum value"
            ),
            member_frame("topology"),
            class_frame("Seq-inst"),
            named_frame("Seq-inst"),
        ]
        .concat()
    );
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:739-764, 444-455 (`ReadBool`,
// `ReadNull`): `TRUE`, `FALSE` and `NULL` only.
#[test]
fn booleans_and_nulls_are_words() {
    let field = |data: &str| format!("User-field ::= {{ label id 1, data {data} }}");
    let data = member_index(named::USER_FIELD, "data");
    let read_data = |text: &str| value(named::USER_FIELD, &field(text)).member(data).cloned();
    let variant =
        |name: &str| member_index(super::members(named::USER_FIELD)[data].ty, name) as u16;
    assert_eq!(
        read_data("bool TRUE"),
        Some(Value::Choice(variant("bool"), Box::new(Value::Bool(true))))
    );
    serial(named::USER_FIELD, &field("bool true"));
    serial(named::USER_FIELD, &field("bool TRUEX"));
    assert_eq!(
        read_data("int +7"),
        Some(Value::Choice(variant("int"), Box::new(Value::Int(7))))
    );
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:877-912, 938-986 (quoted at
// `read::string`): a VisibleString's byte outside 0x20..0x7E is `#` with `ERR_POST_X(7)` naming
// the stack's path, the line and the piece read (probe c09); a UTF8String keeps its bytes (c10);
// `""` is a quote and the lines are joined.
#[test]
fn visible_strings_replace_their_bad_bytes_with_ncbis_message() {
    let (title, posts) = read(named::SEQDESC, "Seqdesc ::= title \"a\x01b\"\n");
    assert_eq!(
        title.unwrap(),
        Value::Choice(
            member_index(named::SEQDESC, "title") as u16,
            Box::new(Value::Str(b"a#b".to_vec()))
        )
    );
    assert_eq!(posts.len(), 1);
    assert_eq!(
        posts[0].line("psiblast"),
        b"Error: (803.7) [psiblast] Seqdesc.title\nBad char [0x1] in string at line 1\na\x01b\n"
    );
    // Two bytes, the line of each piece, `""`, a piece of 128 bytes.
    let long = "x".repeat(130);
    let text = [
        &b"Seqdesc ::= comment \"\x7f\"\"q\nz\x80"[..],
        long.as_bytes(),
        b"\x02\"\n",
    ]
    .concat();
    let (comment, posts) = read(named::SEQDESC, text);
    let Value::Choice(_, comment) = comment.unwrap() else {
        panic!()
    };
    assert_eq!(*comment, Value::Str(format!("#\"qz#{long}#").into_bytes()));
    let messages: Vec<Vec<u8>> = posts.iter().map(|post| post.line("tblastn")).collect();
    assert_eq!(
        messages,
        [
            b"Error: (803.7) [tblastn] Seqdesc.comment\nBad char [0x7F] in string at line 1\n\x7f\n".to_vec(),
            [
                &b"Error: (803.7) [tblastn] Seqdesc.comment\nBad char [0x80] in string at line 2\nz\x80"[..],
                &long.as_bytes()[..126],
                b"\n",
            ]
            .concat(),
            b"Error: (803.7) [tblastn] Seqdesc.comment\nBad char [0x2] in string at line 2\nxxxx\x02\n".to_vec(),
        ]
    );
    // A UTF8String (`SeqTable-multi-data.string SEQUENCE OF UTF8String`) keeps its bytes.
    let utf8 = value(
        named::SEQTABLE_MULTI_DATA,
        "SeqTable-multi-data ::= string {\n  \"\u{e9}\x01\"\n}\n",
    );
    assert_eq!(
        utf8,
        Value::Choice(
            member_index(named::SEQTABLE_MULTI_DATA, "string") as u16,
            Box::new(Value::List(vec![Value::Str(
                "\u{e9}\x01".as_bytes().to_vec()
            )]))
        )
    );
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1477-1538 (quoted at
// `Reader::octets`): an odd last digit is the high nibble (c10 `'abc'H` -> `ABC0`).
#[test]
fn octets_take_an_odd_digit_as_a_high_nibble() {
    let oss = value(
        named::USER_FIELD,
        "User-field ::= { label id 1, data oss { 'abc'H } }",
    );
    let data = oss.member(member_index(named::USER_FIELD, "data")).unwrap();
    let Value::Choice(_, list) = data else {
        panic!()
    };
    assert_eq!(**list, Value::List(vec![Value::Octets(vec![0xab, 0xc0])]));
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (quoted at `Reader::real`):
// the base is `ReadUint4`, so `-10` is `BadNumber` (NCBI 2.17.0 probe e06), and a base other than
// 2 and 10 a `CSerialException`.
#[test]
fn a_real_reads_its_base_as_an_unsigned() {
    let field = |data: &str| format!("User-field ::= {{ label id 1, data {data} }}");
    let data = member_index(named::USER_FIELD, "data");
    let real = value(named::USER_FIELD, &field("real { 140, 10, -3 }"));
    let Value::Choice(_, real) = real.member(data).unwrap() else {
        panic!()
    };
    assert_eq!(**real, Value::Real(0.14));
    assert!(
        exception(named::USER_FIELD, &field("real { 1, -10, 0 }")).contains("bad number in line 1")
    );
    serial(named::USER_FIELD, &field("real { 1, 3, 0 }"));
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:343-424 (quoted at
// `read::item_index`): a member named, in brackets, or by its tag in brackets; an id that starts
// with a digit is `NStr::StringToInt`'s (a `CStringException` for `[3x]`, probe e03, whose frame
// still names the member read before); an upper-case first letter (c21), an unknown member (c04)
// and a member out of order are `CSerialException`s, and so is a missing mandatory member (c20).
#[test]
fn members_are_named_or_tagged_in_module_order() {
    let date = |body: &str| format!("Date-std ::= {{ {body} }}");
    let year = member_index(named::DATE_STD, "year") as u16;
    let month = member_index(named::DATE_STD, "month") as u16;
    let both = Value::Sequence(vec![(year, Value::Int(2024)), (month, Value::Int(5))]);
    assert_eq!(value(named::DATE_STD, &date("year 2024, month 5")), both);
    assert_eq!(value(named::DATE_STD, &date("[0] 2024, [1] 5")), both);
    assert_eq!(
        value(named::DATE_STD, &date("[year] 2024, [month] 5")),
        both
    );
    serial(named::DATE_STD, &date("Year 2024"));
    serial(named::DATE_STD, &date("year 2024, week 5"));
    serial(named::DATE_STD, &date("month 5, year 2024"));
    serial(named::DATE_STD, &date("year 2024, year 2025"));
    serial(named::DATE_STD, &date("month 5"));
    serial(named::DATE_STD, &date(""));
    let report = exception(named::DATE_STD, &date("year 2024, [1x] 5"));
    let string = |line, function: &str, text: &str| {
        report_line(
            "src/corelib/ncbistr.cpp",
            line,
            "CStringException::eConvert",
            &format!("{function}() - {text} (m_Pos = 1)"),
        )
    };
    assert_eq!(
        report,
        [
            "NCBI C++ Exception:\n".to_string(),
            string(870, "ncbi::NStr::StringToInt8", "Cannot convert string '1x' to Int8"),
            report_line(
                ASN,
                1393,
                "CStringException::eConvert",
                "ncbi::CObjectIStreamAsn::ReadClassSequential() -  Frame type= eFrameClassMember, Member name= year (m_Pos = 1)"
            ),
            report_line(
                ASN,
                1396,
                "CStringException::eConvert",
                "ncbi::CObjectIStreamAsn::ReadClassSequential() -  Frame type= eFrameClass, Object type= Date-std (m_Pos = 1)"
            ),
            report_line(
                ISTR,
                968,
                "CStringException::eConvert",
                "ncbi::CObjectIStream::Read() -  Frame type= eFrameNamed, Object type= Date-std (m_Pos = 1)"
            ),
        ]
        .concat()
    );
    // The member frame of the first member names none.
    assert!(exception(named::DATE_STD, &date("[99999999999] 2024"))
        .contains("Frame type= eFrameClassMember (m_Pos = 0)"));
}

// The frames of a number overflow in a descriptor (probe c07 of SURVEY_D3: Date-std's year in a
// `create-date`), with NCBI's lines: the class and member of each SEQUENCE, the choice and its
// variant, the array and its element (an implicit SET OF class names no type).
#[test]
fn the_report_has_a_line_per_frame() {
    let report = exception(
        named::SEQ_DESCR,
        "Seq-descr ::= {\n  title \"t\",\n  create-date std {\n    year 99999999999\n  }\n}\n",
    );
    let choice = |kind: &str, line: u32, what: &str| {
        util(
            ISTR,
            line,
            &format!("ncbi::CObjectIStream::ReadChoiceSimple() -  Frame type= {kind}, {what}"),
        )
    };
    assert_eq!(
        report,
        [
            "NCBI C++ Exception:\n".to_string(),
            util(
                "src/util/strbuffer.cpp",
                583,
                "ncbi::CIStreamBuffer::NumberOverflow() - number overflow in line 4"
            ),
            member_frame("year"),
            class_frame("Date-std"),
            choice("eFrameChoiceVariant", 1498, "Member name= std"),
            choice("eFrameChoice", 1500, "Object type= Date"),
            choice("eFrameChoiceVariant", 1498, "Member name= create-date"),
            choice("eFrameChoice", 1500, "Object type= Seqdesc"),
            util(
                ASN,
                1270,
                "ncbi::CObjectIStreamAsn::ReadContainer() -  Frame type= eFrameArrayElement"
            ),
            util(
                ASN,
                1273,
                "ncbi::CObjectIStreamAsn::ReadContainer() -  Frame type= eFrameArray, Object type= "
            ),
            named_frame("Seq-descr"),
        ]
        .concat()
    );
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1453-1475 (quoted at
// `Reader::choice`): a variant that is not one of the CHOICE's, or none, is a `CSerialException`.
#[test]
fn an_unknown_variant_is_a_serial_exception() {
    serial(named::DATE, "Date ::= when \"x\"");
    serial(named::DATE, "Date ::= }");
    serial(named::DATE, "Date ::= Str \"x\"");
    assert_eq!(
        value(named::DATE, "Date ::= [0] \"x\""),
        Value::Choice(0, Box::new(Value::Str(b"x".to_vec())))
    );
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608 (`HandleEOF`) and
// objistrasn.cpp:207-243 (`SkipWhiteSpace`): the end of the text inside a member is a
// `CSerialException`; directly in the root object it stays a `CEofException` (NCBI 2.17.0 probes
// e01, e02: psiblast ends with exit status 0 and no word).
#[test]
fn the_end_of_the_text_is_an_eof_only_in_the_root_object() {
    for text in [
        "Date-std ::=",
        "Date-std ::= -- a comment",
        "Date-std ::=\n\n",
    ] {
        assert_eq!(read(named::DATE_STD, text).0, Err(ReadError::Eof), "{text}");
    }
    for text in [
        "Date-std ::= {",
        "Date-std ::= { year",
        "Date-std ::= { year 1",
        "Date-std ::= { year 1,",
        "Date-std ::= { season \"x",
    ] {
        serial(named::DATE_STD, text);
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:245-274 (quoted at
// `Reader::skip_white_space`): comments between values are white space (c16), and a `--` inside a
// string is part of it.
#[test]
fn comments_are_white_space() {
    let descr = value(
        named::SEQ_DESCR,
        "Seq-descr ::= { -- one --\n  title \"a -- not a comment --\" -- two\n , name \"n\" }",
    );
    let Value::List(items) = descr else { panic!() };
    assert_eq!(
        items[0],
        Value::Choice(
            member_index(named::SEQDESC, "title") as u16,
            Box::new(Value::Str(b"a -- not a comment --".to_vec()))
        )
    );
    assert_eq!(items.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68 (quoted at
// `hooks::seq_descr_post_read`): an empty `Seq-descr` is NCBI's `CSerialException` (c23).
#[test]
fn an_empty_seq_descr_is_a_serial_exception() {
    serial(named::SEQ_DESCR, "Seq-descr ::= { }");
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869 and objistrasn.cpp:1024-1155
// (quoted at `read::delayed`, `read::skip_number`, `read::string`): `Auth-list.names` is skipped
// and kept as read, from its first byte that is not white space (c11): comments and spacing stay,
// a number is not range-checked (c12), a bad char is posted without the piece; the skip still
// checks the grammar.
#[test]
fn the_delayed_member_is_kept_as_read() {
    let names = "std { { name name { last \"A\" } , level 1 } } -- c --";
    let text = format!("Auth-list ::= {{\n  names   {names}\n}}");
    let list = value(named::AUTH_LIST, &text);
    assert_eq!(
        list,
        Value::Sequence(vec![(
            0,
            Value::Delayed(b"std { { name name { last \"A\" } , level 1 } }".to_vec())
        )])
    );
    let overflow =
        "Auth-list ::= { names std { { name dbtag { db \"x\", tag id 99999999999 } } } }";
    assert!(matches!(
        value(named::AUTH_LIST, overflow),
        Value::Sequence(_)
    ));
    let bad = "Auth-list ::= { names str { \"a\x01\" } }";
    let (list, posts) = read(named::AUTH_LIST, bad);
    assert!(list.is_ok());
    assert_eq!(
        posts[0].line("psiblast"),
        b"Error: (803.7) [psiblast] Auth-list.names.str.E\nBad char [0x1] in string at line 1\n"
    );
    serial(named::AUTH_LIST, "Auth-list ::= { names str { 5 } }");
    serial(
        named::AUTH_LIST,
        "Auth-list ::= { names std { { name dbtag { db \"x\", tag id - 5 } } } }",
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Variation_ref.cpp:56-151 (quoted at
// `hooks::variation_ref_post_read`): the deprecated members move into `variant-prop` (made with
// version 5) or are dropped with NCBI's messages (c13, c13b, c13c).
#[test]
fn variation_ref_moves_and_drops_its_deprecated_members() {
    let variation = |members: &str| format!("Variation-ref ::= {{ {members} data unknown NULL }}");
    let member = |name: &str| member_index(named::VARIATION_REF, name);
    let prop = |name: &str| member_index(named::VARIANTPROPERTIES, name) as u16;
    let read_variation = |members: &str| {
        let (value, posts) = read(named::VARIATION_REF, variation(members));
        (value.unwrap(), posts)
    };
    let (validated, posts) = read_variation("validated TRUE,");
    assert!(posts.is_empty());
    assert_eq!(validated.member(member("validated")), None);
    assert_eq!(
        validated.member(member("variant-prop")),
        Some(&Value::Sequence(vec![
            (prop("version"), Value::Int(5)),
            (prop("other-validation"), Value::Bool(true))
        ]))
    );
    let (state, _) = read_variation("allele-state 1,");
    assert_eq!(
        state.member(member("variant-prop")),
        Some(&Value::Sequence(vec![
            (prop("version"), Value::Int(5)),
            (prop("allele-state"), Value::Enum(1))
        ]))
    );
    let (dropped, posts) = read_variation("population-data { },");
    assert_eq!(dropped.member(member("population-data")), None);
    assert_eq!(
        posts[0].line("psiblast"),
        b"Error: [psiblast] Variation-ref.population-data is deprecated and will be ignored\n"
    );
    let (both, posts) =
        read_variation("variant-prop { version 2, other-validation FALSE }, validated TRUE,");
    assert_eq!(
        both.member(member("variant-prop")),
        Some(&Value::Sequence(vec![
            (prop("version"), Value::Int(2)),
            (prop("other-validation"), Value::Bool(false))
        ]))
    );
    assert_eq!(
        posts[0].line("tblastn"),
        b"Error: [tblastn] Both Variation-ref.validated and Variation-ref.variant-properties.other-validation are set; ignoring Variation-ref.validated\n"
    );
}

// The Seq-id read under the table goes into the Seq-id model, which keeps every member as read
// (`PDB-seq-id.chain` set to its DEFAULT stays set; SURVEY_D3 §6).
#[test]
fn a_seq_id_is_read_into_the_seq_id_model() {
    let pdb = value(
        named::SEQ_ID,
        "Seq-id ::= pdb { mol \"1abc\", chain 32, rel std { year 2000 } }",
    );
    let Value::SeqId(ParsedSeqId::Pdb(pdb)) = pdb else {
        panic!("{pdb:?}")
    };
    assert_eq!(pdb.chain, Some(32));
    assert_eq!(pdb.mol, b"1abc");
    let local = value(named::SEQ_ID, "Seq-id ::= local str \"Query_1\"");
    assert_eq!(
        local,
        Value::SeqId(ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec())))
    );
    let general = value(named::SEQ_ID, "Seq-id ::= general { db \"x\", tag id -3 }");
    let Value::SeqId(ParsedSeqId::General(dbtag)) = general else {
        panic!()
    };
    assert_eq!(dbtag.tag, ObjectId::Id(-3));
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:883-900 (quoted at `read::document`): the
// type named before `::=` must be the root's.
#[test]
fn the_file_header_names_the_root() {
    serial(named::DATE, "Seq-entry ::= str \"x\"");
    assert_eq!(
        value(named::DATE, "Date ::= str \"x\""),
        Value::Choice(0, Box::new(Value::Str(b"x".to_vec())))
    );
}

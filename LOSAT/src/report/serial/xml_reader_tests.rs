//! Tests of `xml_reader`: NCBI's XML object stream (`CObjectIStreamXml`) as rpsbproc's
//! BlastOutput input reads it (objistrxml.cpp, strbuffer.cpp, cited at the functions).

use super::asn_text_reader::{AsnSource, SliceSource, StreamSource};
use super::xml_reader::{member, optional, utf8_of, Member, XmlError, XmlReader};

const OUTER: [Member; 3] = [member("name"), optional("number"), member("inner")];
const INNER: [Member; 2] = [optional("real"), member("big")];

/// What a document of the test module reads to: `(name, number, real, big)`.
type Read = (Vec<u8>, Option<i32>, Option<f64>, i64);

fn read_one<S: AsnSource>(reader: &mut XmlReader<S>) -> Result<Read, XmlError> {
    reader.begin_document("Outer")?;
    let mut name = Vec::new();
    let (mut number, mut real, mut big) = (None, None, 0);
    reader.class("Outer", &OUTER, |reader, index| match index {
        0 => reader.string(&mut name),
        1 => {
            number = Some(reader.int4()?);
            Ok(())
        }
        _ => reader.class("Inner", &INNER, |reader, index| {
            if index == 0 {
                real = Some(reader.real()?);
            } else {
                big = reader.int8()?;
            }
            Ok(())
        }),
    })?;
    Ok((name, number, real, big))
}

fn read(text: &[u8]) -> Result<Read, XmlError> {
    read_one(&mut XmlReader::new(SliceSource::new(text)))
}

fn doc(body: &str) -> Vec<u8> {
    format!("<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<!DOCTYPE Outer PUBLIC \"-//X//Y/EN\" \"x.dtd\">\n<Outer>\n{body}\n</Outer>\n").into_bytes()
}

fn ok(name: &[u8], number: Option<i32>, real: Option<f64>, big: i64) -> Result<Read, XmlError> {
    Ok((name.to_vec(), number, real, big))
}

fn fails<T>(result: Result<T, XmlError>) -> bool {
    matches!(result, Err(XmlError::Format(_)))
}

fn refused<T>(result: Result<T, XmlError>) -> bool {
    matches!(result, Err(XmlError::Unsupported(_)))
}

// ReadClassSequential and BeginClassMember (objistr.cpp:1398-1428, objistrxml.cpp:2185-2367):
// members as `<Class_member>` in the module's order, OPTIONAL ones absent, a mandatory one
// missing, a member out of order or repeated, an unknown member, a member without the class
// prefix, text where a member should be.
#[test]
fn members_in_the_modules_order() {
    let inner = "<Outer_inner><Inner><Inner_real>1.5</Inner_real><Inner_big>7</Inner_big></Inner></Outer_inner>";
    assert_eq!(
        read(&doc(&format!(
            "<Outer_name>ab</Outer_name><Outer_number>-12</Outer_number>{inner}"
        ))),
        ok(b"ab", Some(-12), Some(1.5), 7)
    );
    assert_eq!(
        read(&doc("<Outer_name>ab</Outer_name><Outer_inner><Inner><Inner_big>7</Inner_big></Inner></Outer_inner>")),
        ok(b"ab", None, None, 7)
    );
    // The mandatory `name` missing; `big` missing.
    assert!(fails(read(&doc(&format!(
        "<Outer_number>1</Outer_number>{inner}"
    )))));
    assert!(fails(read(&doc("<Outer_name>a</Outer_name><Outer_inner><Inner><Inner_real>1</Inner_real></Inner></Outer_inner>"))));
    // Out of order, repeated, unknown, without the prefix, another class's prefix.
    assert!(fails(read(&doc(&format!(
        "<Outer_number>1</Outer_number><Outer_name>a</Outer_name>{inner}"
    )))));
    assert!(fails(read(&doc(&format!(
        "<Outer_name>a</Outer_name><Outer_name>a</Outer_name>{inner}"
    )))));
    assert!(fails(read(&doc(&format!(
        "<Outer_name>a</Outer_name><Outer_other>1</Outer_other>{inner}"
    )))));
    assert!(fails(read(&doc(&format!("<name>a</name>{inner}")))));
    assert!(fails(read(&doc(&format!(
        "<Inner_name>a</Inner_name>{inner}"
    )))));
    assert!(fails(read(&doc(&format!(
        "<Outer_name>a</Outer_name>x{inner}"
    )))));
    // A closing tag of another member.
    assert!(fails(read(&doc(&format!(
        "<Outer_name>a</Outer_number>{inner}"
    )))));
    // Another root.
    assert!(fails(read(b"<Inner><Inner_big>1</Inner_big></Inner>")));
}

// ReadFileHeader (objistrxml.cpp:582-670): declaration and DOCTYPE optional, a byte order mark,
// comments and processing instructions; an unknown `<!` tag or a DOCTYPE with an internal subset
// fails.
#[test]
fn file_header() {
    let body = "<Outer><Outer_name>a</Outer_name><Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner></Outer>";
    let expected = ok(b"a", None, None, 1);
    assert_eq!(read(body.as_bytes()), expected);
    assert_eq!(read(format!("\u{feff}{body}").as_bytes()), expected);
    assert!(fails(read(format!("<!-- c -- -->{body}").as_bytes())));
    assert_eq!(
        read(format!("<!-- c - d -->\r\n<?pi a='1'?>{body}").as_bytes()),
        expected
    );
    assert!(fails(read(format!("<!ELEMENT x>{body}").as_bytes())));
    assert!(fails(read(
        format!("<!DOCTYPE Outer [<!ENTITY a 'b'>]>{body}").as_bytes()
    )));
    assert!(fails(read(
        format!("<?xml version=\"1.0\" encoding=\"EBCDIC\"?>{body}").as_bytes()
    )));
    assert_eq!(
        read(format!("<?xml version='1.0' encoding='latin1'?>{body}").as_bytes()),
        expected
    );
    assert!(fails(read(b"<?xml version=\"1.0\"?>")));
}

// ReadTagData, ReadEscapedChar and x_ReadEncodedChar (objistrxml.cpp:722-868,1297-1354):
// entities, character references (UTF-8; 2^31 and above end the text), CDATA, NCBI's
// end-of-line handling, NUL dropped, the document's encoding for bytes from 0x80, control
// characters refused, a comment inside text fails.
#[test]
fn character_data() {
    let name = |text: &str| {
        read(&doc(&format!("<Outer_name>{text}</Outer_name><Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner>"))).map(|r| r.0)
    };
    assert_eq!(
        name("a&amp;b&lt;&gt;&quot;&apos;"),
        Ok(b"a&b<>\"'".to_vec())
    );
    assert_eq!(name("&#65;&#x42;&#xe9;&#0;c"), Ok(b"AB\xc3\xa9c".to_vec()));
    assert_eq!(name("ab&#4294967295;"), Ok(b"ab".to_vec()));
    assert_eq!(name("a<![CDATA[<&x>]]>b"), Ok(b"a<&x>b".to_vec()));
    // After a `\r`, NCBI's flag stays set until a `\n`: each later `\r` is then a `\n`.
    assert_eq!(name("a\r\nb\r\rc\rd\r"), Ok(b"a\nb\nc\nd\n".to_vec()));
    assert_eq!(name("a\rb\r"), Ok(b"ab\n".to_vec()));
    assert!(fails(name("a&nbsp;")));
    assert!(fails(name("a&#X41;")));
    assert!(fails(name("a&#;")));
    assert!(fails(name("a&0123456789012345678901234567890123;")));
    assert!(fails(name("a<!-- c -->b")));
    assert!(refused(name("a&#1;")));
    assert!(refused(name("a\u{7}")));
    let latin1 = |text: &[u8]| {
        let mut text_doc =
            b"<?xml version='1.0' encoding='ISO-8859-1'?><Outer><Outer_name>".to_vec();
        text_doc.extend_from_slice(text);
        text_doc.extend_from_slice(b"</Outer_name><Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner></Outer>");
        read(&text_doc).map(|r| r.0)
    };
    assert_eq!(latin1(b"\xe9"), Ok(b"\xc3\xa9".to_vec()));
    assert_eq!(utf8_of(0x10_FFFF), (*b"\xf4\x8f\xbf\xbf", 4));
    assert_eq!(utf8_of(0xFFFF_FFFF), ([0xff, 0xbf, 0xbf, 0xbf], 4));
    // Tabs, DEL, line ends and NULs between runs of plain bytes, entities next to each other.
    assert_eq!(name("a\tb\u{7f}c\t"), Ok(b"a\tb\x7fc\t".to_vec()));
    assert_eq!(name("\n\na\r\n\nb\rc\nd"), Ok(b"\n\na\n\nbc\nd".to_vec()));
    assert_eq!(name("a\0b\r\0\rc\0"), Ok(b"ab\nc".to_vec()));
    assert_eq!(
        name("x&amp;&amp;y&#233;&#233;z&#13;&#10;w"),
        Ok(b"x&&y\xc3\xa9\xc3\xa9z\nw".to_vec())
    );
    // Bytes from 0x80 in each encoding: as they are, or the UTF-8 bytes of their code point.
    let text: &[u8] = b"a\xe9\x80\x81\xffb";
    let cases: [(Option<&str>, &[u8]); 5] = [
        (None, text),
        (Some("UTF-8"), text),
        (Some("US-ASCII"), text),
        (Some("ISO-8859-1"), b"a\xc3\xa9\xc2\x80\xc2\x81\xc3\xbfb"),
        (Some("windows-1252"), b"a\xc3\xa9\xe2\x82\xac?\xc3\xbfb"),
    ];
    for (encoding, expected) in cases {
        assert_eq!(
            read(&encoded_doc(encoding, text)).map(|r| r.0),
            Ok(expected.to_vec()),
            "{encoding:?}"
        );
    }
}

/// A document of the test module in `encoding` (`None`: no declaration) whose `name` is `text`.
fn encoded_doc(encoding: Option<&str>, text: &[u8]) -> Vec<u8> {
    let mut doc = match encoding {
        Some(name) => format!("<?xml version='1.0' encoding='{name}'?>").into_bytes(),
        None => Vec::new(),
    };
    doc.extend_from_slice(b"<Outer><Outer_name>");
    doc.extend_from_slice(text);
    doc.extend_from_slice(
        b"</Outer_name><Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner></Outer>",
    );
    doc
}

/// The character data `text` (no `<` but CDATA sections; the entities `&amp;` and
/// `&#<decimal>;` only) as the reader takes it one character at a time: each character's UTF-8
/// bytes in the document's `encoding`, NULs dropped, and NCBI's line ends (`read_tag_data`); the
/// bytes of a CDATA section as they are, after which no `\r` is pending (`read_cd_section`).
fn one_by_one(text: &[u8], encoding: Option<&str>) -> Vec<u8> {
    let mut out = Vec::new();
    let mut cr = false;
    let mut i = 0;
    while i < text.len() {
        if let Some(section) = text[i..].strip_prefix(b"<![CDATA[") {
            let end = section.windows(3).position(|w| w == b"]]>").unwrap();
            out.extend_from_slice(&section[..end]);
            i += 9 + end + 3;
            cr = false;
            continue;
        }
        let (bytes, len) = if text[i] == b'&' {
            let end = i + text[i..].iter().position(|&c| c == b';').unwrap();
            let entity = std::str::from_utf8(&text[i + 1..end]).unwrap();
            i = end + 1;
            match entity.strip_prefix('#') {
                Some(digits) => utf8_of(digits.parse().unwrap()),
                None => utf8_of(u32::from(b'&')),
            }
        } else {
            let c = text[i];
            i += 1;
            // (Of cp1252's bytes 0x80 to 0x9F, the texts hold 0x80 and 0x81 only.)
            match (encoding, c) {
                (Some("ISO-8859-1"), 0x80..) => utf8_of(u32::from(c)),
                (Some("windows-1252"), 0x80) => utf8_of(0x20AC),
                (Some("windows-1252"), 0x81) => utf8_of(u32::from(b'?')),
                (Some("windows-1252"), 0xA0..) => utf8_of(u32::from(c)),
                _ => ([c, 0, 0, 0], 1),
            }
        };
        for &c in &bytes[..len] {
            if c == 0 {
                continue;
            }
            if cr {
                if c == b'\n' {
                    cr = false;
                }
                out.push(if c == b'\r' { b'\n' } else { c });
            } else if c == b'\r' {
                cr = true;
            } else {
                out.push(c);
            }
        }
    }
    out
}

// ReadTagData (objistrxml.cpp:1297-1354, cited at `read_tag_data`): the runs of bytes that the
// reader appends at once give what reading one character at a time gives (`one_by_one`), in
// every encoding, for random texts of plain bytes, tabs, line ends, NULs, bytes from 0x80,
// entities and CDATA sections (with `]` inside), also from a stream that keeps one byte at a
// time (a run ends at every byte).
#[test]
fn character_data_in_runs_is_read_as_one_character_at_a_time() {
    const PIECES: [&[u8]; 21] = [
        b"a",
        b"Z",
        b" ",
        b"~",
        b"\t",
        b"\x7f",
        b"\n",
        b"\r",
        b"\0",
        b"\xe9",
        b"\x80",
        b"\x81",
        b"\xc3\xa9",
        b"&amp;",
        b"&#233;",
        b"&#13;",
        b"&#10;",
        b"&#8364;",
        b"<![CDATA[]]>",
        b"<![CDATA[x]y]]z]]]>",
        b"<![CDATA[\r\0\xe9&amp;<a>\n]]>",
    ];
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let mut next = move |bound: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed % bound
    };
    for round in 0..400 {
        let mut text = Vec::new();
        for _ in 0..next(60) {
            // Long runs of one piece as well as mixtures.
            let piece = PIECES[next(PIECES.len() as u64) as usize];
            let times = if next(4) == 0 { 1 + next(20) } else { 1 };
            for _ in 0..times {
                text.extend_from_slice(piece);
            }
        }
        for encoding in [
            None,
            Some("UTF-8"),
            Some("US-ASCII"),
            Some("ISO-8859-1"),
            Some("windows-1252"),
        ] {
            let expected = one_by_one(&text, encoding);
            let doc = encoded_doc(encoding, &text);
            assert_eq!(
                read(&doc).map(|r| r.0),
                Ok(expected.clone()),
                "round {round} {encoding:?}"
            );
            let mut stream = XmlReader::new(StreamSource::new(doc.as_slice(), b""));
            stream.src.fill_one_byte_at_a_time();
            assert_eq!(
                read_one(&mut stream).map(|r| r.0),
                Ok(expected),
                "round {round} {encoding:?} (a byte at a time)"
            );
        }
    }
}

// GetInt4, GetInt8 and ReadDouble (strbuffer.cpp:617-745, objistrxml.cpp:1027-1042): white
// space, signs, overflow, the smallest values; REAL words and NCBI's conversion.
#[test]
fn numbers() {
    let number = |text: &str| {
        read(&doc(&format!("<Outer_name>a</Outer_name><Outer_number>{text}</Outer_number><Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner>"))).map(|r| r.1)
    };
    assert_eq!(number(" \x0b+2147483647 "), Ok(Some(i32::MAX)));
    assert_eq!(number("-2147483647"), Ok(Some(-i32::MAX)));
    assert_eq!(number("-2147483648"), Ok(Some(i32::MIN)));
    assert!(fails(number("2147483648")));
    assert!(fails(number("")));
    assert!(fails(number("1x")));
    assert!(fails(number("&#49;")));
    let real = |text: &str| {
        read(&doc(&format!("<Outer_name>a</Outer_name><Outer_inner><Inner><Inner_real>{text}</Inner_real><Inner_big>1</Inner_big></Inner></Outer_inner>"))).map(|r| r.2)
    };
    assert_eq!(real(" 5.04742e-43 "), Ok(Some(5.04742e-43)));
    assert_eq!(real(""), Ok(Some(0.0)));
    // `nan` and `-nan` are NCBI's `HUGE_VAL/HUGE_VAL`, sign set (printed `-nan`; review oracles
    // O1 `x_nan`, `x_mnan`).
    for text in ["nan", "-nan", "NaN"] {
        let value = real(text).unwrap().unwrap();
        assert!(value.is_nan() && value.is_sign_negative(), "{text}");
    }
    // More than 18 significant digits: NCBI's `StringToDoublePosix` adds the digits after the
    // 18th in 80-bit extended precision (line SW, SW8 B28; NCBI's corelib with
    // `fDecimalPosixFinite`, swebok-improve-20261105/sw8/nfix/probe/x87 `s2dp`: 0x3f847ae147ae147b,
    // 0x3ff0000000000000, 0x3f8eb851eb851eb8).
    for (text, bits) in [
        ("0.0100000000000000002", 0x3f84_7ae1_47ae_147b_u64),
        ("1.000000000000000001", 0x3ff0_0000_0000_0000),
        ("1.50000000000000000001e-2", 0x3f8e_b851_eb85_1eb8),
    ] {
        assert_eq!(
            real(text).map(|r| r.map(f64::to_bits)),
            Ok(Some(bits)),
            "{text}"
        );
    }
    assert_eq!(real("1e999"), Ok(Some(f64::MAX)));
    // An exponent that overflows NCBI's `int` (undefined behaviour; NCBI's wrapped exponent gives
    // DBL_MAX) in LOSAT's widened arithmetic: DBL_MIN with `fDecimalPosixFinite`
    // (PD-LOSAT-NCBI-DEFECTS exception 11).
    assert_eq!(real("1.55e-2147483647"), Ok(Some(f64::MIN_POSITIVE)));
    assert!(fails(real("1e5x")));
    assert!(fails(real("1 2")));
    let big = |text: &str| {
        read(&doc(&format!("<Outer_name>a</Outer_name><Outer_inner><Inner><Inner_big>{text}</Inner_big></Inner></Outer_inner>"))).map(|r| r.3)
    };
    assert_eq!(big("9223372036854775807"), Ok(i64::MAX));
    assert_eq!(big("-9223372036854775808"), Ok(i64::MIN));
    assert!(fails(big("9223372036854775808")));
}

// Self-closed tags (EndOpeningTagSelfClosed, objistrxml.cpp:342-378): an empty string, an empty
// list; an INTEGER or a class cannot be self-closed. Attributes and namespace prefixes are
// refused; white space inside tags is read.
#[test]
fn tags() {
    let inner = "<Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner>";
    assert_eq!(
        read(&doc(&format!("<Outer_name/>{inner}"))),
        ok(b"", None, None, 1)
    );
    assert_eq!(
        read(&doc(&format!("<Outer_name\n>a</Outer_name\t >{inner}"))),
        ok(b"a", None, None, 1)
    );
    assert!(fails(read(&doc(&format!(
        "<Outer_name>a</Outer_name><Outer_number/>{inner}"
    )))));
    assert!(fails(read(&doc(
        "<Outer_name>a</Outer_name><Outer_inner/>"
    ))));
    assert!(refused(read(&doc(&format!(
        "<Outer_name x=\"1\">a</Outer_name>{inner}"
    )))));
    assert!(refused(read(&doc(&format!(
        "<x:Outer_name>a</x:Outer_name>{inner}"
    )))));
    assert!(fails(read(&doc(&format!("<Outer_name/ >{inner}")))));
}

// Several documents from one stream (CIStreamBuffer::Close gives back what was not read):
// each read starts after the previous root's `>`; a truncated document fails.
#[test]
fn documents_one_after_another() {
    let one = "<Outer><Outer_name>a</Outer_name><Outer_inner><Inner><Inner_big>1</Inner_big></Inner></Outer_inner></Outer>";
    let two = "<?xml version=\"1.0\"?>\n<Outer><Outer_name>b</Outer_name><Outer_inner><Inner><Inner_big>2</Inner_big></Inner></Outer_inner></Outer>\n";
    let text = format!("{one}{two}<Outer><Outer_name>c");
    let mut stream = XmlReader::new(StreamSource::new(text.as_bytes(), b""));
    let mut names = Vec::new();
    loop {
        let mut name = Vec::new();
        let result = stream.begin_document("Outer").and_then(|()| {
            stream.class("Outer", &OUTER, |reader, index| match index {
                0 => reader.string(&mut name),
                1 => reader.int4().map(drop),
                _ => reader.class("Inner", &INNER, |reader, index| {
                    if index == 0 {
                        reader.real().map(drop)
                    } else {
                        reader.int8().map(drop)
                    }
                }),
            })
        });
        match result {
            Ok(()) => names.push(name),
            Err(error) => {
                assert!(matches!(error, XmlError::Format(_)));
                break;
            }
        }
    }
    assert_eq!(names, vec![b"a".to_vec(), b"b".to_vec()]);
}

// ReadTagData and ReadCDSection (objistrxml.cpp:1271-1354, cited at `read_tag_data`): the
// character data of a member read from a stream is the data read from the whole text, with its
// entities, line ends and CDATA sections, and the stream does not keep the bytes read
// (`StreamSource::release_before`): a text of megabytes is held once, in the string.
#[test]
fn long_character_data_is_not_kept_by_the_stream() {
    let mut body = String::new();
    for i in 0..60_000 {
        body.push_str(match i % 7 {
            0 => "&amp;&#65;&#x42;&#x3b1;z",
            1 => "\r\n",
            2 => "\r\r\n",
            3 => "<![CDATA[<x>]]>",
            _ => "abcdefghijklmnopqrstuvwxyz0123456789",
        });
    }
    // One CDATA section of 3 MiB.
    body.push_str(&format!("<![CDATA[{}]]>", "c".repeat(3 << 20)));
    let text = doc(&format!(
        "<Outer_name>{body}</Outer_name>\n<Outer_inner><Inner><Inner_big>7</Inner_big></Inner></Outer_inner>"
    ));
    let expected = read(&text).unwrap();
    assert!(expected.0.len() > 4_000_000);
    let mut stream = XmlReader::new(StreamSource::new(text.as_slice(), b""));
    stream.begin_document("Outer").unwrap();
    let mut name = Vec::new();
    let mut largest = 0;
    stream
        .class("Outer", &OUTER, |reader, index| match index {
            0 => {
                reader.string(&mut name)?;
                largest = reader.src.buffered();
                Ok(())
            }
            1 => reader.int4().map(drop),
            _ => reader.class("Inner", &INNER, |reader, index| {
                if index == 0 {
                    reader.real().map(drop)
                } else {
                    assert_eq!(reader.int8()?, 7);
                    Ok(())
                }
            }),
        })
        .unwrap();
    assert_eq!(name, expected.0);
    assert!(largest < 8 * 64 * 1024, "the stream kept {largest} bytes");
}

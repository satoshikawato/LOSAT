//! Tests of `xml_reader`: NCBI's XML object stream (`CObjectIStreamXml`) as rpsbproc's
//! BlastOutput input reads it (objistrxml.cpp, strbuffer.cpp, quoted at the functions).

use super::asn_text_reader::{SliceSource, StreamSource};
use super::xml_reader::{member, optional, utf8_of, Member, XmlError, XmlReader};

const OUTER: [Member; 3] = [member("name"), optional("number"), member("inner")];
const INNER: [Member; 2] = [optional("real"), member("big")];

/// What a document of the test module reads to: `(name, number, real, big)`.
type Read = (Vec<u8>, Option<i32>, Option<f64>, i64);

fn read_one(reader: &mut XmlReader<SliceSource<'_>>) -> Result<Read, XmlError> {
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
    assert_eq!(utf8_of(0x10_FFFF), b"\xf4\x8f\xbf\xbf".to_vec());
    assert_eq!(utf8_of(0xFFFF_FFFF), vec![0xff, 0xbf, 0xbf, 0xbf]);
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
    assert!(real("nan").unwrap().unwrap().is_nan());
    assert_eq!(real("1e999"), Ok(Some(f64::MAX)));
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

//! Unit tests of the objects of a stream (`super::stream`): NCBI's `EndOfData` between objects,
//! one object read and the rest left, the memory of a stream, and the encoding guess on the
//! archives and the bad inputs that NCBI BLAST+ 2.17.0's blast_formatter was given in the
//! probes of line RG-F3 (`F/f3/ncbi/C`, SURVEY_F3_NCBI §2.3: P3.001-124).

use std::io::Cursor;

use super::read::{read_document, ReadError};
use super::stream::{guess_encoding, text_stream, Encoding};
use super::table::named;
use super::Value;
use crate::report::serial::asn_input::Names;

/// Whether `error` is a `CSerialException` whose report holds `message`.
fn is_serial(error: &Option<ReadError>, message: &str) -> bool {
    matches!(error, Some(ReadError::Serial(e)) if String::from_utf8_lossy(&e.what()).contains(message))
}

/// The objects of `text` read as `root` until `EndOfData`, and the error that ended them.
fn objects(text: &[u8], root: u16) -> (Vec<Value>, Option<ReadError>) {
    let mut stream = text_stream(text, b"", Names::TEST);
    let mut out = Vec::new();
    while !stream.end_of_data() {
        match stream.read_object(root).0 {
            Ok(value) => out.push(value),
            Err(error) => return (out, Some(error)),
        }
    }
    (out, None)
}

/// White space and comments between and after the objects are skipped; the end of the input
/// inside a comment ends the data too. Anything else after an object is the next object: a word
/// without `::=` before the end is an end of the input in the root object (`CEofException`,
/// which leaves NCBI's reader), a byte that is no word a format error.
// NCBI 2.17.0 (P3.087-091): a "-- comment" line before the archive, trailing white space and a
// trailing "junk" line are read without an error message, a trailing NUL is "Invalid input
// format for BLAST Archive." after the first archive's output.
#[test]
fn end_of_data_skips_white_space_and_comments() {
    let id = "Seq-id ::= local str \"a\"";
    let text = format!("-- c\n{id}\n\n -- c --{id}\t-- unclosed");
    let (values, error) = objects(text.as_bytes(), named::SEQ_ID);
    assert_eq!((values.len(), error), (2, None));
    let (values, error) = objects(format!("{id}\njunk\n").as_bytes(), named::SEQ_ID);
    assert_eq!((values.len(), error), (1, Some(ReadError::Eof)));
    let (values, error) = objects(format!("{id}\0").as_bytes(), named::SEQ_ID);
    assert_eq!(values.len(), 1);
    assert!(is_serial(&error, "'::=' expected"), "{error:?}");
    let (values, error) = objects(b"", named::SEQ_ID);
    assert_eq!((values.len(), error), (0, None));
    let (values, error) = objects(b" \n\t\r\n", named::SEQ_ID);
    assert_eq!((values.len(), error), (0, None));
}

/// A failed read leaves the stream failed: `EndOfData` is true after it.
#[test]
fn a_failed_read_ends_the_data() {
    let mut stream = text_stream(&b"Seq-id ::= locl 5\nSeq-id ::= gi 5"[..], b"", Names::TEST);
    assert!(!stream.end_of_data());
    assert!(stream.read_object(named::SEQ_ID).0.is_err());
    assert!(stream.end_of_data());
}

/// One object read leaves the stream after its last byte: the rest goes back to the caller, a
/// second object of the text included, and the first bytes given apart from the input are read
/// first.
#[test]
fn one_object_read_leaves_the_rest() {
    let text = b"Seq-id ::= local str \"a\"\n-- c\nSeq-id ::= gi 5\n";
    let mut stream = text_stream(&text[4..], &text[..4], Names::TEST);
    let (value, posts) = stream.read_object(named::SEQ_ID);
    assert!(posts.is_empty());
    assert_eq!(
        value,
        read_document(b"Seq-id ::= local str \"a\"", named::SEQ_ID, Names::TEST).0
    );
    let (mut rest, mut input) = stream.into_rest();
    std::io::Read::read_to_end(&mut input, &mut rest).expect("a slice reads");
    assert_eq!(rest, &text[24..]);
}

/// `seekg` back to the start of a file: the input from its start again.
#[test]
fn a_seekable_input_rewinds() {
    let text = b"Blast4-request ::= { body get-matrices NULL }\n".to_vec();
    let mut stream = text_stream(Cursor::new(text.clone()), b"", Names::TEST);
    // NCBI's strategy reader tries the reply first: the header names another type.
    let (reply, _) = stream.read_object(named::BLAST4_GET_SEARCH_STRATEGY_REPLY);
    assert!(
        is_serial(
            &reply.err(),
            "incompatible type Blast4-request<>Blast4-get-search-strategy-reply"
        ),
        "the reply's read fails"
    );
    let input = stream.rewind().expect("a cursor seeks");
    let guessed = guess_encoding(input);
    assert_eq!(guessed.encoding.as_ref().ok(), Some(&Encoding::Text));
    let mut stream = guessed.into_text(Names::TEST);
    assert!(stream.read_object(named::BLAST4_REQUEST).0.is_ok());
    assert!(stream.end_of_data());
}

/// A stream keeps about one object, not the text before it: objects without blocks release
/// the bytes before them as each one starts.
#[test]
fn a_stream_keeps_one_object() {
    let one = "Seq-id ::= local str \"an identifier of some length\"\n";
    let text = one.repeat(40_000);
    let mut stream = text_stream(text.as_bytes(), b"", Names::TEST);
    let mut count = 0;
    while !stream.end_of_data() {
        assert!(stream.read_object(named::SEQ_ID).0.is_ok());
        count += 1;
    }
    assert_eq!(count, 40_000);
    assert!(
        stream.buffered() <= 10 * 64 * 1024,
        "{} bytes kept of {}",
        stream.buffered(),
        text.len()
    );
}

/// A delayed member longer than the blocks a stream releases keeps its text: the value is the
/// one read from the whole text at once.
#[test]
fn a_delayed_member_keeps_its_text_in_a_stream() {
    let names: Vec<String> = (0..60_000).map(|n| format!("\"name {n}\"")).collect();
    let text = format!(
        "Auth-list ::= {{ names str {{ {} }}, affil str \"x\" }}\n",
        names.join(",\n")
    );
    let whole = read_document(text.as_bytes(), named::AUTH_LIST, Names::TEST).0;
    let (values, error) = objects(text.as_bytes(), named::AUTH_LIST);
    assert_eq!(error, None);
    assert_eq!(values.len(), 1);
    assert_eq!(Ok(&values[0]), whole.as_ref());
    assert!(
        matches!(&values[0], Value::Sequence(members) if matches!(members[0].1, Value::Delayed { ref bytes, .. } if bytes.len() > 5 * 64 * 1024))
    );
}

/// The fixtures of this test, by name (`LOSAT/tests/fixtures/archive_read`).
fn fixture(name: &str) -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/archive_read")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The encoding of a text, binary and XML archive (`blastn -outfmt 11` with `ARCHIVE_FORMAT`
/// unset, `bin`, `xml`), and of the inputs NCBI's blast_formatter refused, as NCBI's guess
/// takes them. NCBI 2.17.0 formatted the three archives alike (P3.001-015) and printed "Invalid
/// input format for BLAST Archive." for every input below that is not text (P3.077-083, 085):
/// its guess was another format (the reader's `CBlastException`), or the binary reader failed
/// (the UTF-8 BOM: the guess is binary ASN.1, whose reader does not take the text); for the text
/// ones it is the reader's error on a type that is not `Blast4-archive` (P3.093, 121).
#[test]
fn the_encoding_is_guessed_as_ncbi_guesses_it() {
    let text = fixture("blastn_multi.asn");
    let guess = |bytes: &[u8]| -> Encoding {
        guess_encoding(bytes)
            .encoding
            .unwrap_or_else(|stop| panic!("{stop:?}"))
    };
    let mut bom = b"\xef\xbb\xbf".to_vec();
    bom.extend_from_slice(&text);
    let cases: Vec<(&str, Vec<u8>, Encoding)> = vec![
        ("text", text.clone(), Encoding::Text),
        ("binary", fixture("blastn_multi.bin"), Encoding::Binary),
        ("XML", fixture("blastn_multi.xml"), Encoding::Xml),
        ("UTF-8 BOM", bom, Encoding::Binary),
        (
            "comment",
            [&b"-- comment\n"[..], &text].concat(),
            Encoding::Text,
        ),
        (
            "CRLF",
            text.iter()
                .flat_map(|&c| if c == b'\n' { vec![b'\r', c] } else { vec![c] })
                .collect(),
            Encoding::Text,
        ),
        (
            "leading white space",
            [&b"\n\n  "[..], &text].concat(),
            Encoding::Text,
        ),
        ("first byte", text[..1].to_vec(), Encoding::Other),
        ("first 20 bytes", text[..20].to_vec(), Encoding::Text),
        ("empty", fixture("bad/empty"), Encoding::Other),
        ("newline", fixture("bad/newline"), Encoding::Other),
        ("spaces", fixture("bad/spaces"), Encoding::Other),
        ("hello world", fixture("bad/garbage"), Encoding::Other),
        ("JSON", fixture("bad/json"), Encoding::Other),
        ("bytes 0-255", fixture("bad/binjunk"), Encoding::Other),
        ("gzip", fixture("bad/gzip"), Encoding::Other),
        ("Seq-entry", fixture("bad/other_type"), Encoding::Text),
        ("strategy", fixture("bad/strategy"), Encoding::Text),
    ];
    for (name, bytes, encoding) in cases {
        assert_eq!(guess(&bytes), encoding, "{name}");
    }
    // A text read from where the guess left the stream: the guess's bytes, then the input.
    let mut stream = guess_encoding(&text[..]).into_text(Names::TEST);
    assert!(stream.read_object(named::BLAST4_ARCHIVE).0.is_ok());
    assert!(stream.end_of_data());
    let (values, error) = objects(&fixture("bad/other_type"), named::BLAST4_ARCHIVE);
    assert_eq!(values.len(), 0);
    assert!(
        is_serial(&error, "incompatible type Seq-entry<>Blast4-archive"),
        "{error:?}"
    );
}

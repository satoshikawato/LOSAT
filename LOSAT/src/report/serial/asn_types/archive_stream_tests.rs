//! Unit tests of the object streams (`super::stream`) on archives: every archive among LOSAT's
//! fixtures read and written back as read, also through the binary and XML writers and readers;
//! NCBI's own binary and XML archives read as their text; the lines the binary and XML readers
//! post; and the objects of the inputs that NCBI BLAST+ 2.17.0's blast_formatter was given in the
//! probes of line RG-F3 (`F/f3/ncbi/C`, SURVEY_F3_NCBI §2.3: P3.077-124; impl D's probes D1-D3),
//! with the objects delivered and how the stream ended.

use std::path::Path;

use super::stream::{guess_encoding, Encoding, ObjectStream, XmlStyle};
use super::table::named;
use super::write::write_document;
use super::{ReadError, Value};
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::Encoder;
use crate::report::serial::xml::XmlEncoder;

/// A fixture of this test, by name (`LOSAT/tests/fixtures/archive_read`).
fn fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/archive_read")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The archives of rpsbproc's tests that are copies of NCBI's edited by hand
/// (`fixtures/rpsbproc/README.md`: a member id in brackets, `length -0`, an `org` descriptor, no
/// `inst.length`, no query `seq-data`, a `delta` query): rpsbproc's own cases, not NCBI's writer's
/// layout (`length -0` is not even read).
const EDITED: [&str; 6] = [
    "fixtures/rpsbproc/root/bracket_member.asn",
    "fixtures/rpsbproc/root/neg0_length.asn",
    "fixtures/rpsbproc/root/org_descr.asn",
    "fixtures/rpsbproc/root/nolength.asn",
    "fixtures/rpsbproc/arch/q_nodata.asn",
    "fixtures/rpsbproc/arch/mfirst_delta.asn",
];

/// The marks of the fixtures that are templates of NCBI's output: a test puts its run's directory
/// in place of the mark (`@CWD@`, `{db_dir}`) before it compares, and NCBI broke the long strings
/// that hold the directory into lines by the real path's length, so the template's own bytes are
/// not a layout NCBI's writer gives.
const TEMPLATE_MARKS: [&[u8]; 2] = [b"@CWD@", b"{db_dir}"];

/// Every archive among LOSAT's fixtures: the files under `LOSAT/tests` that start with an archive
/// header (NCBI's `-outfmt 11` output and the files that LOSAT's tests freeze as NCBI's), with
/// their bytes. Left out: the malformed archives of rpsbproc's tests
/// (`fixtures/rpsbproc/arch/bad`), its copies edited by hand (`EDITED`), the `.merged` files of
/// `archive_replaced` (a run's stdout and stderr in one file, with NCBI's messages inside the
/// archive's text), and the templates (`TEMPLATE_MARKS`).
fn archive_fixtures() -> Vec<(std::path::PathBuf, Vec<u8>)> {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let bad = tests.join("fixtures/rpsbproc/arch/bad");
    let mut todo = vec![tests];
    let mut files = Vec::new();
    while let Some(dir) = todo.pop() {
        for entry in std::fs::read_dir(&dir).expect("a directory") {
            let path = entry.expect("a directory entry").path();
            if path.is_dir() {
                if path != bad {
                    todo.push(path);
                }
                continue;
            }
            if path
                .extension()
                .is_some_and(|ext| ext == "rs" || ext == "merged")
            {
                continue;
            }
            if EDITED.iter().any(|name| path.ends_with(name)) {
                continue;
            }
            let text = std::fs::read(&path).expect("a file");
            let template = TEMPLATE_MARKS
                .iter()
                .any(|mark| text.windows(mark.len()).any(|window| window == *mark));
            if text.starts_with(b"Blast4-archive ::=") && !template {
                files.push((path, text));
            }
        }
    }
    files.sort();
    files
}

/// The objects of a stream (any encoding) read to its end and written back as ASN.1 text, with
/// the number of objects; or why the reading stopped (an error, or a message posted).
fn as_text<R: std::io::Read>(mut stream: ObjectStream<R>) -> Result<(Vec<u8>, usize), String> {
    let mut out = Vec::new();
    let mut objects = 0;
    while !stream.end_of_data() {
        let (value, posts) = stream.read(named::BLAST4_ARCHIVE);
        if !posts.is_empty() {
            return Err(format!("messages {posts:?}"));
        }
        let value = value.map_err(|e| format!("{e:?}"))?;
        let mut enc = AsnTextEncoder::new(&mut out);
        write_document(&mut enc, named::BLAST4_ARCHIVE, &value);
        assert_eq!(enc.error(), None);
        objects += 1;
    }
    Ok((out, objects))
}

/// Where two texts first differ, with the bytes around it.
fn first_difference(expected: &[u8], got: &[u8]) -> String {
    let at = got
        .iter()
        .zip(expected)
        .position(|(a, b)| a != b)
        .unwrap_or(got.len().min(expected.len()));
    let around = |bytes: &[u8]| {
        String::from_utf8_lossy(&bytes[at.saturating_sub(40)..(at + 40).min(bytes.len())])
            .into_owned()
    };
    format!(
        "written otherwise from byte {at}: {:?} / {:?}",
        around(expected),
        around(got)
    )
}

/// Every archive of the fixtures is read whole, object by object, and its values written back
/// are its bytes.
#[test]
fn every_archive_of_the_fixtures_is_written_back_as_read() {
    let (mut files, mut objects) = (0, 0);
    let mut failures = Vec::new();
    for (path, text) in archive_fixtures() {
        files += 1;
        let error = match as_text(ObjectStream::new(&text[..], b"")) {
            Ok((out, count)) => {
                objects += count;
                (out != text).then(|| first_difference(&text, &out))
            }
            Err(error) => Some(error),
        };
        if let Some(error) = error {
            failures.push(format!("{}: {error}", path.display()));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {files} archives: {:#?}",
        failures.len(),
        &failures[..failures.len().min(10)]
    );
    assert!(
        files >= 10 && objects >= files,
        "{files} archives, {objects} objects"
    );
}

/// Every archive of the fixtures read as text, each object written by `write` (one encoding's
/// writer: `Some` of what it could not write) and the whole read back by the stream `stream` makes
/// of those bytes, written back as text: its bytes. The archives that the writer does not write
/// (a UTF8String, a delayed member: `SerialError::Unwritten`) are counted, not compared.
fn every_archive_reads_back(
    write: impl Fn(&Value, &mut Vec<u8>) -> Option<String>,
    stream: impl Fn(&[u8]) -> ObjectStream<&[u8]>,
) {
    let (mut files, mut unwritten) = (0, Vec::new());
    let mut failures = Vec::new();
    for (path, text) in archive_fixtures() {
        let mut objects = ObjectStream::new(&text[..], b"");
        let mut bytes = Vec::new();
        let mut written = true;
        while !objects.end_of_data() {
            let value = objects
                .read(named::BLAST4_ARCHIVE)
                .0
                .expect("a text fixture");
            if let Some(error) = write(&value, &mut bytes) {
                unwritten.push(format!("{}: {error}", path.display()));
                written = false;
                break;
            }
        }
        if !written {
            continue;
        }
        files += 1;
        let error = match as_text(stream(&bytes)) {
            Ok((out, _)) => (out != text).then(|| first_difference(&text, &out)),
            Err(error) => Some(error),
        };
        if let Some(error) = error {
            failures.push(format!("{}: {error}", path.display()));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {files} archives: {:#?}",
        failures.len(),
        &failures[..failures.len().min(10)]
    );
    // The writers do not write a UTF8String (the binary writer's tag depends on NCBI's
    // parameter `SERIAL_WRITE_UTF8STRING_TAG`; the XML documents LOSAT writes hold none): the
    // archives that hold one are left out here, and NCBI's own binary and XML archives of them
    // are read below (`ncbis_binary_archives_read_as_their_text`, `..._xml_...`).
    assert!(
        files >= 100,
        "{files} archives read back, {} not written: {:#?}",
        unwritten.len(),
        &unwritten[..unwritten.len().min(10)]
    );
}

/// The values of every archive of the fixtures, written by NCBI's binary writer
/// (`AsnBinaryEncoder`) and read by the binary reader, are the values of the text.
#[test]
fn every_archive_of_the_fixtures_reads_back_from_binary() {
    every_archive_reads_back(
        |value, out| {
            let mut enc = AsnBinaryEncoder::new(out);
            write_document(&mut enc, named::BLAST4_ARCHIVE, value);
            enc.error().map(|error| format!("{error:?}"))
        },
        |bytes| ObjectStream::binary(bytes, b""),
    );
}

/// The values of every archive of the fixtures, written by NCBI's XML writer of the archive
/// (`XmlEncoder::archive`, `ARCHIVE_FORMAT=xml`) and read by the XML reader in the ASN.1 style,
/// are the values of the text. A value whose writing posts a replaced character is not written
/// as it is (counted as not written).
#[test]
fn every_archive_of_the_fixtures_reads_back_from_xml() {
    every_archive_reads_back(
        |value, out| {
            let mut enc = XmlEncoder::archive(out);
            write_document(&mut enc, named::BLAST4_ARCHIVE, value);
            match enc.error() {
                Some(error) => Some(format!("{error:?}")),
                None => enc.has_posts().then(|| "a character replaced".to_string()),
            }
        },
        |bytes| ObjectStream::xml(bytes, b"", XmlStyle::Asn),
    );
}

/// The archives that NCBI BLAST+ 2.17.0 wrote of one search in ASN.1 text and in another
/// encoding (`ARCHIVE_FORMAT=bin` / `=xml`; the probes of line RG-F3 and of impl D: one object, a
/// batch, a psiblast PSSM, `-subject`, an error, a query list of Seq-locs) are read through the
/// encoding's guess as the same values: written back as text, the bytes of the text archive.
fn ncbis_archives_read_as_their_text(encoding: Encoding, pairs: &[(&str, &str)]) {
    for &(text_name, other_name) in pairs {
        let text = fixture(text_name);
        let other = fixture(other_name);
        let guessed = guess_encoding(&other[..]);
        assert_eq!(guessed.encoding.as_ref().ok(), Some(&encoding));
        let (out, objects) = as_text(guessed.into_stream().expect("a stream"))
            .unwrap_or_else(|error| panic!("{other_name}: {error}"));
        assert!(objects >= 1);
        assert!(
            out == text,
            "{other_name}: {}",
            first_difference(&text, &out)
        );
    }
}

#[test]
fn ncbis_binary_archives_read_as_their_text() {
    ncbis_archives_read_as_their_text(
        Encoding::Binary,
        &[
            ("blastn_multi.asn", "blastn_multi.bin"),
            ("blastn_batch.asn", "blastn_batch.bin"),
            ("psiblast_in_pssm.asn", "psiblast_in_pssm.bin"),
            ("psiblast_it3.asn", "psiblast_it3.bin"),
            ("blastp_subject.asn", "blastp_subject.bin"),
            ("error_blastp.asn", "error_blastp.bin"),
            ("blastp_query_loc.asn", "blastp_query_loc.bin"),
        ],
    );
}

/// The XML archives (the psiblast ones are 200 kB and not copied: the round trip through
/// `XmlEncoder::archive` covers the PSSM).
#[test]
fn ncbis_xml_archives_read_as_their_text() {
    ncbis_archives_read_as_their_text(
        Encoding::Xml,
        &[
            ("blastn_multi.asn", "blastn_multi.xml"),
            ("blastn_batch.asn", "blastn_batch.xml"),
            ("blastp_subject.asn", "blastp_subject.xml"),
            ("error_blastp.asn", "error_blastp.xml"),
            ("blastp_query_loc.asn", "blastp_query_loc.xml"),
        ],
    );
}

/// The lines that the reader of a stream posts, in NCBI's form (blast_formatter: no diagnostic
/// prefix), and the objects it reads before the end of the data or an error.
fn posted(bytes: &[u8]) -> (Vec<u8>, usize) {
    let mut stream = guess_encoding(bytes).into_stream().expect("a stream");
    let (mut lines, mut objects) = (Vec::new(), 0);
    while !stream.end_of_data() {
        let (value, posts) = stream.read(named::BLAST4_ARCHIVE);
        for post in posts {
            lines.extend(post.line(None));
        }
        if value.is_err() {
            break;
        }
        objects += 1;
    }
    (lines, objects)
}

/// A binary archive (`blastp -subject`, impl D's probes of 2026-11-20 in `F3/tmp/d/probe/`)
/// whose VisibleStrings NCBI's reader posts: a UTF8String's tag for a VisibleString once in the
/// process (`ERR_POST_X_ONCE`, the program and the archive read twice: one line), and each byte
/// outside 0x20..0x7E of the query's title with the position after the string's bytes in the
/// stream (the second archive's from the start of the stream). The lines are NCBI's standard
/// error of `blast_formatter -archive`.
#[test]
fn binary_posts_are_ncbis() {
    let archive = fixture("blastp_subject.bin");
    let replace = |bytes: &[u8], from: &[u8], to: &[u8]| -> Vec<u8> {
        let at = bytes
            .windows(from.len())
            .position(|window| window == from)
            .expect("the bytes to replace");
        [&bytes[..at], to, &bytes[at + from.len()..]].concat()
    };
    let utf8 = replace(
        &replace(&archive, b"\x1a\x06blastp", b"\x0c\x06blastp"),
        b"\x1a\x05plain",
        b"\x0c\x05plain",
    );
    let warning = b"Warning: (803.10) CObjectIStreamAsnBinary: UTF8String data for VisibleString \
                    member Blast4-archive.request.body.queue-search.program, ASN.1 specification \
                    may need an update\n";
    assert_eq!(posted(&utf8), (warning.to_vec(), 1));
    assert_eq!(posted(&[&utf8[..], &utf8].concat()), (warning.to_vec(), 2));
    let title = b"Blast4-archive.request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.\
                  E.title";
    let bad = replace(
        &archive,
        b"\x1a\x0eq1 small query",
        b"\x1a\x0eq1 sm\x01ll\x7fquery",
    );
    let lines = |at: usize| {
        [
            &b"Error: (803.7) "[..],
            title,
            format!("\nBad char [0x1] in string at byte {at}\n").as_bytes(),
            b"q1 sm\x01ll\x7fquery\n",
            b"Error: (803.7) ",
            title,
            format!("\nBad char [0x7F] in string at byte {at}\n").as_bytes(),
            b"q1 sm#ll\x7fquery\n",
        ]
        .concat()
    };
    assert_eq!(posted(&bad), (lines(193), 1));
    assert_eq!(
        posted(&[&bad[..], &bad].concat()),
        ([lines(193), lines(193 + 3914)].concat(), 2)
    );
}

/// An XML archive (`blastp -subject`, `ARCHIVE_FORMAT=xml`, the fixture `blastp_subject.xml`)
/// changed as impl D's probe D3 of 2026-11-20 changed it, formatted by NCBI's blast_formatter
/// 2.17.0 (`F3/ncbi/D3/out3*`): a control character in the query's title is posted once for 0x01
/// and not for 0x7F, with the line it is on (the second archive's from the start of the stream),
/// unknown attributes, comments, a processing instruction of the forms `<?pi?>` and
/// `<?pi a="1"?>` and a self-closed title read; a processing instruction with a word that is no
/// attribute (`<?pi x?>`, `<?pi hello world?>`) ends the input as a serial error before any
/// archive.
#[test]
fn xml_variants_read_as_ncbi_read_them() {
    let archive = fixture("blastp_subject.xml");
    let title: &[u8] = b"<Seqdesc_title>q1 small query</Seqdesc_title>";
    let with = |bytes: &[u8], to: &[u8]| -> Vec<u8> {
        let at = bytes
            .windows(title.len())
            .position(|window| window == title)
            .expect("the title");
        [&bytes[..at], to, &bytes[at + title.len()..]].concat()
    };
    let line = archive[..archive
        .windows(title.len())
        .position(|window| window == title)
        .expect("the title")]
        .iter()
        .filter(|&&c| c == b'\n')
        .count()
        + 1;
    assert_eq!(line, 56);
    let lines = archive.iter().filter(|&&c| c == b'\n').count();
    let bad = with(
        &archive,
        b"<Seqdesc_title>q1 sm\x01ll\x7fquery</Seqdesc_title>",
    );
    let post = |line: usize| {
        format!(
            "Error: (803.7) Blast4-archive.request.body.queue-search.queries.bioseq-set.seq-set.\
             seq.descr.title\nBad char [0x1] in string at line {line}\n"
        )
        .into_bytes()
    };
    assert_eq!(posted(&bad), (post(line), 1));
    assert_eq!(
        posted(&[&bad[..], &bad].concat()),
        ([post(line), post(line + lines)].concat(), 2)
    );
    use End::{Data, Serial};
    let cases: Vec<(&str, Vec<u8>, End)> = vec![
        (
            "attributes",
            with(
                &archive,
                b"<Seqdesc_title a=\"1\" b='2'>q1 small query</Seqdesc_title>",
            ),
            Data,
        ),
        (
            "comment",
            with(
                &archive,
                b"<!-- c -->\n<Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Data,
        ),
        (
            "comment after",
            with(
                &archive,
                b"<Seqdesc_title>q1 small query</Seqdesc_title><!-- c -->",
            ),
            Data,
        ),
        (
            "two comments with a hyphen",
            with(
                &archive,
                b"<!-- a --><!-- a-b --><Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Data,
        ),
        (
            "PI without words",
            with(
                &archive,
                b"<?pi?><Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Data,
        ),
        (
            "PI with an attribute",
            with(
                &archive,
                b"<?pi a=\"1\"?><Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Data,
        ),
        (
            "PI after",
            with(
                &archive,
                b"<Seqdesc_title>q1 small query</Seqdesc_title><?pi a=\"1\"?>",
            ),
            Data,
        ),
        (
            "comment and PI",
            with(
                &archive,
                b"<!-- c -->\n<?pi a=\"1\"?><Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Data,
        ),
        (
            "self-closed title",
            with(&archive, b"<Seqdesc_title/>"),
            Data,
        ),
        (
            "PI with a word",
            with(
                &archive,
                b"<?pi x?><Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Serial,
        ),
        (
            "PI with words",
            with(
                &archive,
                b"<?pi hello world?><Seqdesc_title>q1 small query</Seqdesc_title>",
            ),
            Serial,
        ),
    ];
    for (name, bytes, end) in cases {
        assert_eq!(
            outcome(&bytes),
            Ok((usize::from(end == Data), end)),
            "{name}"
        );
    }
}

/// How the objects of an archive text end.
#[derive(Debug, PartialEq)]
enum End {
    /// `EndOfData`: blast_formatter ends with 0.
    Data,
    /// The reader's `CEofException`, which blast_formatter's `CATCH_ALL` ends without a message
    /// and with 0.
    Eof,
    /// The reader's `CSerialException`: "Invalid input format for BLAST Archive.", 1.
    Serial,
}

/// The archive objects read of `text` and how they end; the encoding guessed for an input that
/// LOSAT reads in no encoding.
fn outcome(text: &[u8]) -> Result<(usize, End), Encoding> {
    let guessed = guess_encoding(text);
    let encoding = match guessed.encoding.as_ref() {
        Ok(encoding) => *encoding,
        Err(stop) => panic!("{stop:?}"),
    };
    let Some(mut stream) = guessed.into_stream() else {
        return Err(encoding);
    };
    let mut count = 0;
    while !stream.end_of_data() {
        match stream.read(named::BLAST4_ARCHIVE).0 {
            Ok(_) => count += 1,
            Err(ReadError::Eof) => return Ok((count, End::Eof)),
            Err(ReadError::Serial(_)) => return Ok((count, End::Serial)),
            Err(other) => panic!("{other:?}"),
        }
    }
    Ok((count, End::Data))
}

/// The inputs of NCBI's probes, made from one archive (`blastn -outfmt 11` of three queries,
/// `blastn_multi.asn`), the same archive as binary ASN.1 and XML, and a batch of three
/// archives (`blastn_batch.asn`, `.bin`): the archives LOSAT reads before the end of the data or
/// an error, and the error, are the archives NCBI formatted and the code it ended with.
#[test]
fn the_objects_of_a_stream_end_as_ncbis_did() {
    let text = fixture("blastn_multi.asn");
    let batch = fixture("blastn_batch.asn");
    let binary = fixture("blastn_multi.bin");
    let binary_batch = fixture("blastn_batch.bin");
    let binary_with = |after: &[u8]| [&binary[..], after].concat();
    let binary_cut = |count: usize| binary[..binary.len() - count].to_vec();
    let subject = fixture("blastp_subject.bin");
    let subject_with = |after: &[u8]| [&subject[..], after].concat();
    // The root `30 80 ... 00 00` as `30 82 HH LL ...` with the length of its content and `more`.
    let definite = |more: usize| {
        let content = &subject[2..subject.len() - 2];
        let length = u16::try_from(content.len() + more).expect("a 2-byte length");
        [&b"\x30\x82"[..], &length.to_be_bytes(), content].concat()
    };
    let xml = fixture("blastn_multi.xml");
    let with = |before: &[u8], after: &[u8]| [before, &text[..], after].concat();
    let crlf: Vec<u8> = text
        .iter()
        .flat_map(|&c| if c == b'\n' { vec![b'\r', c] } else { vec![c] })
        .collect();
    use End::{Data, Eof, Serial};
    let cases: Vec<(&str, Vec<u8>, Result<(usize, End), Encoding>)> = vec![
        // 0, the archives formatted (P3.087, 089-092, 096-097, 115-116).
        ("text", text.clone(), Ok((1, Data))),
        ("CRLF", crlf, Ok((1, Data))),
        ("leading white space", with(b"\n\n  ", b""), Ok((1, Data))),
        (
            "a comment before",
            with(b"-- comment\n", b""),
            Ok((1, Data)),
        ),
        (
            "trailing white space",
            with(b"", b"\n\n  \n"),
            Ok((1, Data)),
        ),
        ("twice", with(b"", &text[..]), Ok((2, Data))),
        ("no last newline", text[..12628].to_vec(), Ok((1, Data))),
        ("batch cut after 1", batch[..9366].to_vec(), Ok((1, Data))),
        ("batch cut after 2", batch[..11323].to_vec(), Ok((2, Data))),
        // 0 after the archive, no message: a line "junk" is a header cut short (P3.088).
        ("junk after", with(b"", b"junk\n"), Ok((1, Eof))),
        // 1 after the archives read before the error (P3.090, 117-118).
        ("NUL after", with(b"", b"\0"), Ok((1, Serial))),
        ("batch cut in 2", batch[..10344].to_vec(), Ok((1, Serial))),
        ("batch cut in 3", batch[..13187].to_vec(), Ok((2, Serial))),
        // 1, nothing formatted (P3.093, 098-104, 121).
        ("first 20 bytes", text[..20].to_vec(), Ok((0, Serial))),
        ("cut in the request", text[..1800].to_vec(), Ok((0, Serial))),
        (
            "cut before the results",
            text[..3601].to_vec(),
            Ok((0, Serial)),
        ),
        ("cut in the results", text[..8115].to_vec(), Ok((0, Serial))),
        ("no last 2 bytes", text[..12627].to_vec(), Ok((0, Serial))),
        ("no last 10 bytes", text[..12619].to_vec(), Ok((0, Serial))),
        ("Seq-entry", fixture("bad/other_type"), Ok((0, Serial))),
        ("search strategy", fixture("bad/strategy"), Ok((0, Serial))),
        // 1, nothing formatted: the guess is no text (P3.077-086, 099).
        ("first byte", text[..1].to_vec(), Err(Encoding::Other)),
        ("empty", fixture("bad/empty"), Err(Encoding::Other)),
        ("newline", fixture("bad/newline"), Err(Encoding::Other)),
        ("spaces", fixture("bad/spaces"), Err(Encoding::Other)),
        ("hello world", fixture("bad/garbage"), Err(Encoding::Other)),
        ("JSON", fixture("bad/json"), Err(Encoding::Other)),
        ("bytes 0-255", fixture("bad/binjunk"), Err(Encoding::Other)),
        ("gzip", fixture("bad/gzip"), Err(Encoding::Other)),
        // The guess is binary ASN.1, whose reader stops at the first byte (P3.085).
        ("UTF-8 BOM", with(b"\xef\xbb\xbf", b""), Ok((0, Serial))),
        // Binary ASN.1 (P3.006-010, 104; P7.001-032): formatted; 1 cut at 50 % and 90 % or by 2
        // to 16 bytes (the end of the data inside the member frame); 0 and nothing formatted
        // without its last byte (in the root class's `EndClass`, the one frame left after it is
        // the root's); 1 after the archive with `00 00` after it (no SEQUENCE's tag); 0 after it
        // with a stray `30` (the end of the data in the root class's length).
        ("binary", binary.clone(), Ok((1, Data))),
        (
            "binary cut at 50 %",
            binary[..2461].to_vec(),
            Ok((0, Serial)),
        ),
        (
            "binary cut at 90 %",
            binary[..4429].to_vec(),
            Ok((0, Serial)),
        ),
        ("binary no last byte", binary_cut(1), Ok((0, Eof))),
        ("binary no last 2 bytes", binary_cut(2), Ok((0, Serial))),
        ("binary no last 3 bytes", binary_cut(3), Ok((0, Serial))),
        ("binary no last 4 bytes", binary_cut(4), Ok((0, Serial))),
        ("binary no last 5 bytes", binary_cut(5), Ok((0, Serial))),
        ("binary no last 6 bytes", binary_cut(6), Ok((0, Serial))),
        ("binary no last 8 bytes", binary_cut(8), Ok((0, Serial))),
        ("binary no last 10 bytes", binary_cut(10), Ok((0, Serial))),
        ("binary no last 11 bytes", binary_cut(11), Ok((0, Serial))),
        ("binary no last 12 bytes", binary_cut(12), Ok((0, Serial))),
        ("binary no last 16 bytes", binary_cut(16), Ok((0, Serial))),
        ("binary 00 00 after", binary_with(b"\0\0"), Ok((1, Serial))),
        ("binary 30 after", binary_with(b"\x30"), Ok((1, Eof))),
        ("binary twice", binary_with(&binary), Ok((2, Data))),
        (
            "binary twice, no last byte",
            binary_with(&binary_cut(1)),
            Ok((1, Eof)),
        ),
        ("binary batch", binary_batch.clone(), Ok((3, Data))),
        // Impl D's probes on `blastp_subject.bin` (2026-11-20): one byte `00` after it, or `30 80`
        // (the end of the data in the next root's member frame), or `30 80 00 00` (no member
        // `request`): 1, then invalid input. The root with a definite length `82 HH LL`: NCBI
        // counts it from the first length byte, so a length of the content alone ends the
        // class 2 bytes early (invalid input) and one 2 longer is the archive.
        ("binary 00 after", subject_with(b"\0"), Ok((1, Serial))),
        (
            "binary 30 80 after",
            subject_with(b"\x30\x80"),
            Ok((1, Serial)),
        ),
        (
            "binary 30 80 00 00 after",
            subject_with(b"\x30\x80\0\0"),
            Ok((1, Serial)),
        ),
        ("binary definite root", definite(0), Ok((0, Serial))),
        ("binary definite root + 2", definite(2), Ok((1, Data))),
        (
            "binary definite root twice",
            [definite(0), definite(0)].concat(),
            Ok((0, Serial)),
        ),
        (
            "binary batch no last byte",
            binary_batch[..binary_batch.len() - 1].to_vec(),
            Ok((2, Eof)),
        ),
        // XML (P3.011-015, 092, 106-107): formatted; 1 cut at 50 % or without the last 20 bytes;
        // formatted without its last byte (the end of line after the root's closing tag); 1
        // after the archive with an element after it (the next object's root `junk`).
        ("XML", xml.clone(), Ok((1, Data))),
        ("XML cut at 50 %", xml[..25305].to_vec(), Ok((0, Serial))),
        ("XML no last byte", xml[..50610].to_vec(), Ok((1, Data))),
        (
            "XML no last 20 bytes",
            xml[..50591].to_vec(),
            Ok((0, Serial)),
        ),
        (
            "XML element after",
            [&xml[..], b"<junk/>"].concat(),
            Ok((1, Serial)),
        ),
    ];
    for (name, bytes, expected) in cases {
        assert_eq!(outcome(&bytes), expected, "{name}");
    }
}

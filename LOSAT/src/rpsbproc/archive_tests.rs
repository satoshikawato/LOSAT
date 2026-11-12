//! The Blast4-archive reader on NCBI rpsblast 2.17.0 archives (the fixtures of
//! `LOSAT/tests/fixtures/rpsbproc/arch/`) and on the edited archives that NCBI's rpsbproc read
//! as the end of its input (PS16a probe 3: `arch/bad/`).

use super::*;
use crate::report::serial::asn_text_reader::AsnTextErrorKind;

const ONE: &[u8] = include_bytes!("../../tests/fixtures/rpsbproc/arch/one.asn");

/// The objects of `text`, read as rpsbproc reads them: the outcome of each `next_archive`
/// until the end (`Ok(None)`) or a refusal.
fn read_all(text: &[u8]) -> Vec<std::result::Result<Archive, AsnTextError>> {
    let start = text
        .iter()
        .position(|&b| (b as i8) > b' ' as i8)
        .expect("a byte that is not white space");
    let mut stream = ArchiveStream::new(&text[start + 1..], &text[start..=start]);
    let mut out = Vec::new();
    loop {
        match stream.next_archive() {
            Ok(Some(archive)) => out.push(Ok(archive)),
            Ok(None) => return out,
            Err(error) => {
                out.push(Err(error));
                assert!(matches!(stream.next_archive(), Ok(None)));
                return out;
            }
        }
    }
}

#[test]
fn ncbi_archive_is_read() {
    let archives = read_all(ONE);
    assert_eq!(archives.len(), 1);
    let archive = archives[0].as_ref().unwrap();
    assert_eq!(archive.ident.as_deref(), Some(&b"2.17.0+"[..]));
    let RequestBody::QueueSearch(search) = &archive.body else {
        panic!("queue-search");
    };
    assert_eq!(search.program, b"blastp");
    assert!(matches!(&search.subject, Subject::Database(db) if db == b"Cdd"));
    assert!(search.has_algorithm_options);
    assert!(matches!(&search.matrix_name, Some(ParamValue::String(name)) if name == b"BLOSUM62"));
    assert!(matches!(
        search.evalue_threshold,
        Some(ParamValue::Cutoff(Cutoff::EValue(value))) if value == 0.01
    ));
    let Queries::BioseqSet(entries) = &search.queries else {
        panic!("bioseq-set");
    };
    assert_eq!(entries.len(), 1);
    let SeqEntry::Seq(query) = &entries[0] else {
        panic!("seq");
    };
    assert_eq!(
        query.ids,
        vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))]
    );
    assert_eq!(query.mol, 3);
    assert_eq!(query.length, Some(79));
    assert!(query
        .title
        .as_ref()
        .unwrap()
        .starts_with(b"NP_414783.1 protein YkfF"));
    assert!(query.unread_descriptor.is_none());
    let aligns = archive.results.alignments.as_ref().unwrap();
    assert_eq!(aligns.len(), 1);
    let Segs::Denseg(seg) = &aligns[0].segs else {
        panic!("denseg");
    };
    assert_eq!(seg.starts, vec![4, 0]);
    assert_eq!(seg.lens, vec![70]);
    assert!(matches!(&seg.ids[1], ParsedSeqId::General(tag) if tag.db == b"CDD"));
    assert_eq!(aligns[0].score, Some(337));
    assert_eq!(aligns[0].num_ident, Some(60));
    assert_eq!(aligns[0].e_value, Some(504742318164089.0 * 1e-57));
    assert!(aligns[0].score_of_other_type.is_none());
    assert!(archive.results.masks.as_ref().unwrap().is_empty());
}

// PS16a probe 3: NCBI's rpsbproc reads two objects without a line between them, and stops
// without an error at an object that its typed reader refuses.
#[test]
fn objects_follow_each_other_and_a_bad_object_ends_the_input() {
    let twice = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/twice_nonl.asn");
    assert_eq!(read_all(twice).len(), 2);
    for bad in [
        &include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/dupmember.asn")[..],
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/unkmember.asn"),
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/trailcomma.asn"),
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/trunc.asn"),
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/noresults.asn"),
    ] {
        assert!(read_all(bad).is_empty());
    }
    let second_bad = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/second_bad.asn");
    assert_eq!(read_all(second_bad).len(), 1);
    let trail = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/trailgarb.asn");
    assert_eq!(read_all(trail).len(), 1);
    // Comments and CRLF line ends are white space.
    for good in [
        &include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/comment.asn")[..],
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/crlf.asn"),
    ] {
        assert_eq!(read_all(good).len(), 1);
    }
}

// A form that LOSAT does not read is refused once the object reads whole; when NCBI's reader
// fails later in the same object, the object ends the input as NCBI's does.
#[test]
fn refusals_wait_for_the_whole_object() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let long_real = text.replace(
        "{ 504742318164089, 10, -57 }",
        "{ 5047423181640891234, 10, -61 }",
    );
    assert_ne!(long_real, text);
    let read = read_all(long_real.as_bytes());
    assert_eq!(read.len(), 1);
    let error = read[0].as_ref().unwrap_err();
    assert_eq!(error.kind, AsnTextErrorKind::Unsupported);
    assert!(error.message.contains("18 significant digits"), "{error}");
    let then_bad = long_real.replace("gapped FALSE", "gapped MAYBE");
    assert_ne!(then_bad, long_real);
    assert!(read_all(then_bad.as_bytes()).is_empty());
    // A byte outside the visible range in a VisibleString, here the title (NCBI replaces it by
    // '#' and posts an error).
    let invisible = text.replacen(
        "title \"NP_414783.1 protein YkfF",
        "title \"NP_414783.1 protein Yk\u{e9}F",
        1,
    );
    assert_ne!(invisible, text);
    let read = read_all(invisible.as_bytes());
    assert!(read[0].as_ref().unwrap_err().message.contains("byte 0xc3"));
}

// The `CFastaReader` DefLine is a UTF8String, which NCBI reads with any byte (PS16a review
// S2-4, oracles O3-4 and O3-5): a TAB or UTF-8 there is read, and so are the other values of a
// User-field by their types.
#[test]
fn a_user_field_string_takes_any_byte() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let defline = text.replacen(
        "data str \">NP_414783.1 protein YkfF",
        "data str \">NP_414783.1 protein\t\u{3b1}-YkfF",
        1,
    );
    assert_ne!(defline, text);
    let read = read_all(defline.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok());
    let fields = text.replacen(
        "label str \"DefLine\",",
        "label str \"Other\",\n data strs { \"a\tb\", \"\u{e9}\" } },\n {\n label id 5,\n num 2,\n data fields { { label str \"x\", data object { type id 1, data { { label id 2, data reals { { 5, 10, 0 } } } } } } } },\n {\n label str \"DefLine\",",
        1,
    );
    assert_ne!(fields, text);
    let read = read_all(fields.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok(), "{:?}", read[0]);
    // A User-field value that does not read as its type is NCBI's error.
    let bad_int = fields.replacen("label id 5,", "label id 5, data int x }, { label id 6,", 1);
    assert!(read_all(bad_int.as_bytes()).is_empty());
}

/// A reader that gives one byte per `read`.
struct OneByte<'a>(&'a [u8]);

impl std::io::Read for OneByte<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let Some((&first, rest)) = self.0.split_first() else {
            return Ok(0);
        };
        buf[0] = first;
        self.0 = rest;
        Ok(1)
    }
}

// The stream source refills at every byte, and a text of many blocks drops the bytes it has read
// (`StreamSource::release_before`); the objects read as from one slice.
#[test]
fn a_stream_refills_and_releases() {
    let three = include_bytes!("../../tests/fixtures/rpsbproc/arch/three.asn");
    let mut text = Vec::new();
    while text.len() <= 8 * 64 * 1024 {
        text.extend_from_slice(three);
    }
    let copies = text.len() / three.len();
    let whole = read_all(&text);
    assert_eq!(whole.len(), copies);
    let mut stream = ArchiveStream::new(OneByte(&text[1..]), &text[..1]);
    let mut largest = 0;
    for expected in &whole {
        let archive = stream.next_archive().unwrap().unwrap();
        assert_eq!(
            format!("{archive:?}"),
            format!("{:?}", expected.as_ref().unwrap())
        );
        largest = largest.max(stream.buffered());
    }
    assert!(matches!(stream.next_archive(), Ok(None)));
    // The bytes kept: at most the four blocks read before a release and what follows them.
    assert!(largest < 5 * 64 * 1024, "{largest} of {}", text.len());
}

// The type name of the root object: another type is NCBI's error, no name is read as the
// archive (`SkipFileHeader`, objistr.cpp:892-898).
#[test]
fn the_root_type_name_is_checked() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let other = text.replacen("Blast4-archive ::=", "Blast4-request ::=", 1);
    assert!(read_all(other.as_bytes()).is_empty());
    let unnamed = text.replacen("Blast4-archive ::=", "::=", 1);
    assert_eq!(read_all(unnamed.as_bytes()).len(), 1);
}

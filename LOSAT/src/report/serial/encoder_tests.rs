use super::{Encoder, NullElement, Stream, StreamPosition, Traced};
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::xml::XmlEncoder;

/// The frozen standard error of NCBI BLAST+ 2.17.0 psiblast for the rows of a NULL element of
/// `Seq-hist.replaces.ids` in `-outfmt 11` (line RG D3c-3c): the ASN.1 text and binary archive of
/// `tests/fixtures/pssm/d3b_hist_ids.asn.gz`, the XML archive of `d3c_r_hist_ids.asn.gz`.
const TEXT_ERR: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/db_search/db.ptr.ps.o11.hist_ids.err"
));
const BINARY_ERR: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/db_search/db.ptr.ps.o11.hist_ids.bin.err"
));
const XML_ERR: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/db_search/db.ptr.ps.o11.r_hist_ids.xml.err"
));

/// The trace of the element and the frames up to its list (`Traced::stack_trace`,
/// `Traced::frame_names`), as the archive's writer leaves them at that element.
const TRACE: &str =
    "Blast4-archive.request.body.queue-search.queries.pssm.pssm.query.seq.inst.hist.replaces.ids";
const FRAMES: [&str; 13] = [
    "Blast4-archive",
    "request",
    "body",
    "queue-search",
    "queries",
    "pssm",
    "pssm",
    "query",
    "seq",
    "inst",
    "hist",
    "replaces",
    "ids",
];

fn null_element(stream: Stream, position: StreamPosition, extra: &[&str]) -> NullElement {
    let element_frames = match stream {
        // The ASN.1 streams' frames of the container and its element; NCBI's XML stream has
        // neither for an element of a named type (`Seq-id`).
        Stream::AsnText | Stream::AsnBinary => 2,
        Stream::Xml => 0,
    };
    NullElement {
        stream,
        element_frame: false,
        trace: format!("{TRACE}{}", if element_frames > 0 { ".E" } else { "" }),
        position,
        frames: FRAMES
            .iter()
            .chain(extra)
            .map(|frame| frame.to_string())
            .collect(),
    }
}

/// NCBI's text of what `cli::NativeError::null_element` reports (`post`, then `CATCH_ALL`'s
/// `Error: ` and `what()` and an end of line): the frozen bytes of NCBI, before LOSAT's words.
fn ncbi_report(thrown: &NullElement) -> Vec<u8> {
    let mut report = thrown.post(Some("psiblast"));
    report.extend_from_slice(b"Error: ");
    report.extend_from_slice(&thrown.exception().ncbi_what());
    report.push(b'\n');
    report
}

// D3c-3c review R4: the place of the throw (`NullElement::site`: file, line and function of each
// stream), which LOSAT's words of the report leave out and the rows therefore do not compare.
#[test]
fn the_report_of_each_stream_is_ncbis_frozen_text() {
    let text = null_element(Stream::AsnText, StreamPosition::Line(28), &["[]", "[]"]);
    assert_eq!(
        String::from_utf8_lossy(&ncbi_report(&text)),
        String::from_utf8_lossy(TEXT_ERR)
    );
    let binary = null_element(Stream::AsnBinary, StreamPosition::Byte(216), &["[]", "[]"]);
    assert_eq!(
        String::from_utf8_lossy(&ncbi_report(&binary)),
        String::from_utf8_lossy(BINARY_ERR)
    );
    let xml = null_element(Stream::Xml, StreamPosition::Line(59), &[]);
    assert_eq!(
        String::from_utf8_lossy(&ncbi_report(&xml)),
        String::from_utf8_lossy(XML_ERR)
    );
}

// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1255-1300
// The branch of an element that has its own frame (an anonymous element type) throws at the
// second site of the XML stream; the ASN.1 streams have one site whatever the element.
#[test]
fn an_element_frame_moves_only_the_xml_site() {
    let site = |stream: Stream, element_frame: bool| {
        let thrown = NullElement {
            element_frame,
            ..null_element(stream, StreamPosition::Line(1), &[])
        };
        let what = String::from_utf8(thrown.exception().ncbi_what()).unwrap();
        let place = what.split_once(".cpp\", ").unwrap().1;
        place.split_once(": Error").unwrap().0.to_string()
    };
    assert_eq!(site(Stream::Xml, false), "line 1272");
    assert_eq!(site(Stream::Xml, true), "line 1296");
    for element_frame in [false, true] {
        assert_eq!(site(Stream::AsnText, element_frame), "line 597");
        assert_eq!(site(Stream::AsnBinary, element_frame), "line 1244");
    }
}

/// `Traced` over `enc` writes a SEQUENCE `Doc` whose member `list` holds one element stated as
/// NULL of the type `element_type` (empty: an anonymous type), inside the named type `named` when
/// there is one, and gives NCBI's throw there.
fn thrown<E: Encoder>(
    mut enc: Traced<E>,
    named: Option<&'static str>,
    element_type: &str,
) -> NullElement {
    enc.begin_document("Doc");
    enc.begin_sequence("Doc");
    enc.member("list", 0);
    if let Some(named) = named {
        enc.begin_named(named);
    }
    enc.begin_list(false);
    enc.null_element(element_type);
    *enc.take_thrown().expect("the stream throws at the element")
}

// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1186-1209
// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1366-1376
// The frames that the exception names: NCBI's XML stream has a frame for an element of an
// anonymous type and for a named type, none for its containers; the ASN.1 streams the reverse.
#[test]
fn the_xml_frames_of_an_anonymous_element() {
    let mut buf = Vec::new();
    let xml = thrown(Traced::xml(XmlEncoder::archive(&mut buf), 0), None, "");
    assert_eq!(xml.stream, Stream::Xml);
    assert!(xml.element_frame);
    assert_eq!(xml.trace, "Doc.list.E");
    assert_eq!(xml.frames, ["Doc", "list", "[]"]);
    let what = String::from_utf8(xml.exception().ncbi_what()).unwrap();
    assert!(what.contains("objostrxml.cpp\", line 1296: "), "{what}");
    assert!(what.ends_with(" ( at Doc.list.[])\n"), "{what}");

    let mut buf = Vec::new();
    let named = thrown(
        Traced::xml(XmlEncoder::archive(&mut buf), 0),
        Some("Seq-descr"),
        "Seqdesc",
    );
    assert!(!named.element_frame);
    assert_eq!(named.trace, "Doc.list");
    assert_eq!(named.frames, ["Doc", "list", "Seq-descr"]);

    let mut buf = Vec::new();
    let text = thrown(
        Traced::new(AsnTextEncoder::new(&mut buf), 0),
        Some("Seq-descr"),
        "",
    );
    assert_eq!(text.stream, Stream::AsnText);
    assert!(text.element_frame);
    assert_eq!(text.trace, "Doc.list.E");
    assert_eq!(text.frames, ["Doc", "list", "[]", "[]"]);
}

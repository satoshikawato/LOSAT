use super::*;
use crate::seqid::BlastDefline;

fn defline(ids: Vec<ParsedSeqId>, title: Option<&str>) -> BlastDefline {
    BlastDefline {
        title: title.map(|t| t.as_bytes().to_vec()),
        ids,
        ..Default::default()
    }
}

fn text_ids(spec: &str) -> Vec<ParsedSeqId> {
    parse_ids(spec.as_bytes(), 0, &mut |_| {})
        .unwrap()
        .into_iter()
        .flatten()
        .collect()
}

#[test]
fn the_title_joins_the_later_deflines_with_their_ids() {
    let set = DeflineSet {
        deflines: vec![
            defline(
                text_ids("ref|NM_000001.1|"),
                Some("alpha mRNA [Homo sapiens]"),
            ),
            defline(text_ids("gb|AB000001.1|"), Some("alpha cDNA")),
            defline(text_ids("emb|X00001.2|"), None),
            defline(vec![], Some("no ids")),
        ],
        utf8: false,
    };
    assert_eq!(
        set.bioseq_title(),
        b"alpha mRNA [Homo sapiens] >gb|AB000001.1| alpha cDNA >emb|X00001.2| no ids"
    );
    let configured = configure_defline_title(&set.bioseq_title(), false).unwrap();
    assert_eq!(
        configured,
        b"alpha mRNA [Homo sapiens] >AB000001.1 alpha cDNA >X00001.2 no ids"
    );
    let with_ctrl_a = configure_defline_title(&set.bioseq_title(), true).unwrap();
    assert_eq!(
        with_ctrl_a,
        b"alpha mRNA [Homo sapiens]\x01AB000001.1 alpha cDNA\x01X00001.2 no ids"
    );
    // A part that is not a Seq-id keeps its ` >`.
    assert_eq!(
        configure_defline_title(b"a >not an id >gi|5 x", true).unwrap(),
        b"a >not an id\x01gi|5 x"
    );
    assert_eq!(
        rejoin_defline_title(b"a >not an id >gi|5 x", true).unwrap(),
        b"a >not an id\x01gi|5 x"
    );
}

#[test]
fn the_short_header_has_three_forms() {
    let header = |ids: Vec<ParsedSeqId>, range, minus| {
        let set = DeflineSet {
            deflines: vec![defline(ids, Some("title"))],
            utf8: false,
        };
        let config = Config {
            strand_minus: minus,
            ..Config::default()
        };
        let mut out = Vec::new();
        short_header(&set, range, &config, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    };
    assert_eq!(header(vec![ordinal_id(3)], None, false), ">title\n");
    assert_eq!(
        header(vec![ordinal_id(3)], Some((2, 8)), true),
        ">:c9-3 title\n"
    );
    assert_eq!(
        header(text_ids("lcl|mynuc6"), None, false),
        ">mynuc6 title\n"
    );
    assert_eq!(
        header(text_ids("lcl|mynuc6"), Some((4, 19)), false),
        ">mynuc6:5-20 title\n"
    );
    assert_eq!(
        header(text_ids("gi|654321|gb|AB000002.3|"), None, false),
        ">AB000002.3 title\n"
    );
    assert_eq!(
        header(text_ids("gi|654321|gb|AB000002.3|"), Some((4, 19)), true),
        ">AB000002.3:c20-5 title\n"
    );
}

/// A front line without Seq-ids (an edited header; NCBI's `ProcessFastaDeflines` finds no first
/// id) gives no header line at all; a later one only adds its title (line RG D1c-1c-b).
#[test]
fn a_front_line_without_ids_writes_no_short_header() {
    let set = DeflineSet {
        deflines: vec![
            defline(vec![], Some("delta1")),
            defline(text_ids("ref|NP_000006.1|"), Some("delta2")),
        ],
        utf8: false,
    };
    let mut out = Vec::new();
    short_header(&set, None, &Config::default(), &mut out).unwrap();
    assert!(out.is_empty());
    let set = DeflineSet {
        deflines: vec![
            defline(text_ids("ref|NP_000002.1|"), Some("beta1")),
            defline(vec![], Some("beta2")),
        ],
        utf8: false,
    };
    short_header(&set, None, &Config::default(), &mut out).unwrap();
    assert_eq!(out, b">NP_000002.1 beta1 beta2\n");
}

#[test]
fn split_keeps_empty_pieces() {
    assert_eq!(
        split_by_pattern(b"a >b > >c", b" >"),
        [&b"a"[..], b"b", b"", b"c"]
    );
    assert_eq!(split_by_pattern(b" >x", b" >"), [&b""[..], b"x"]);
    assert_eq!(split_by_pattern(b"", b" >"), [&b""[..]]);
}

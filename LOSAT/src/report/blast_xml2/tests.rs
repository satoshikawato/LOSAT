//! Unit tests of the XML2 / JSON documents: the pieces of NCBI BLAST+ 2.17.0's outfmt 15 and 16
//! (probes `blastp_sub_invalidonly`, `blastp_sub_nohfirst`; `outfmt-20261019/S/probe5/out`) written
//! from records.

use super::*;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::{
    QueryMask, StructuredXml, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats, REFERENCE_GAPPED_BLAST,
};
use crate::report::serial::encoder::Encoder;
use crate::report::serial::json::JsonEncoder;
use crate::report::serial::xml::XmlEncoder;

const KA: &str = "Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options";

/// `blastp_sub_invalidonly_fmt16.out`: one query of an unsearched batch, two `-subject` records.
const INVALID_ONLY_16: &str = r####"<?xml version="1.0" encoding="US-ASCII"?>
<BlastXML2
    xmlns="http://www.ncbi.nlm.nih.gov"
    xmlns:xs="http://www.w3.org/2001/XMLSchema-instance"
    xs:schemaLocation="http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd"
>
<BlastOutput2>
  <report>
    <Report>
      <program>blastp</program>
      <version>BLASTP 2.17.0+</version>
      <reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</reference>
      <search-target>
        <Target>
          <subjects>Subject_1</subjects>
          <subjects>Subject_2</subjects>
        </Target>
      </search-target>
      <params>
        <Parameters>
          <matrix>BLOSUM62</matrix>
          <expect>10</expect>
          <gap-open>11</gap-open>
          <gap-extend>1</gap-extend>
          <filter>F</filter>
          <cbs>2</cbs>
        </Parameters>
      </params>
      <results>
        <Results>
          <bl2seq>
            <Search>
              <query-id>Query_1</query-id>
              <query-title>inv all X</query-title>
              <query-len>30</query-len>
              <stat>
                <Statistics>
                  <hsp-len>0</hsp-len>
                  <eff-space>0</eff-space>
                  <kappa>-1</kappa>
                  <lambda>-1</lambda>
                  <entropy>-1</entropy>
                </Statistics>
              </stat>
              <message>Query_1 inv all X:  Query_1 inv all X: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options  No hits found &#xa;No hits found</message>
            </Search>
          </bl2seq>
        </Results>
      </results>
    </Report>
  </report>
</BlastOutput2>
</BlastXML2>
"####;

/// `blastp_sub_invalidonly_fmt15.out`.
const INVALID_ONLY_15: &str = r####"{
"BlastOutput2": [
{
    "report": {
      "program": "blastp",
      "version": "BLASTP 2.17.0+",
      "reference": "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs\", Nucleic Acids Res. 25:3389-3402.",
      "search_target": {
        "subjects": [
          "Subject_1",
          "Subject_2"
        ]
      },
      "params": {
        "matrix": "BLOSUM62",
        "expect": 10,
        "gap_open": 11,
        "gap_extend": 1,
        "filter": "F",
        "cbs": 2
      },
      "results": {
        "bl2seq": [
          {
            "query_id": "Query_1",
            "query_title": "inv all X",
            "query_len": 30,
            "stat": {
              "hsp_len": 0,
              "eff_space": 0,
              "kappa": -1,
              "lambda": -1,
              "entropy": -1
            },
            "message": "Query_1 inv all X:  Query_1 inv all X: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options  No hits found \u000ANo hits found"
          }
        ]
      }
    }
  }
]
}
"####;

fn run<'a>(matrix: Option<&'a FormattingMatrix>, db_name: &'static [u8]) -> XmlRun<'a> {
    XmlRun {
        program: "blastp",
        version: "BLASTP 2.17.0+",
        reference: REFERENCE_GAPPED_BLAST,
        db_name,
        matrix: Some("BLOSUM62"),
        expect: 10.0,
        gap_open: 11,
        gap_extend: 1,
        reward_penalty: None,
        filter: b"F",
        db_num: 4,
        db_len: 5000,
        positives: matrix,
        mask_char: b'X',
    }
}

fn params(target: Xml2Target) -> Xml2Run {
    Xml2Run {
        target,
        gapped: true,
        cbs: 2,
        query_gencode: None,
        db_gencode: None,
    }
}

fn invalid_query(message: &[u8]) -> XmlQuery<'_> {
    XmlQuery {
        id: b"Query_1",
        title: b"inv all X",
        len: 30,
        stats: XmlStats::NONE,
        message: Some(message),
        searched: false,
        masks: &[],
    }
}

fn invalid_message() -> Vec<u8> {
    format!("Query_1 inv all X:  Query_1 inv all X: {KA}  No hits found").into_bytes()
}

fn document(kind: Xml2Kind, target: Xml2Target, finish: bool) -> String {
    let message = invalid_message();
    let run = run(None, b"");
    let mut writer = Xml2Writer::new(kind, params(target));
    let mut out = Vec::new();
    writer
        .query(&run, &invalid_query(&message), &mut out, |_| {
            panic!("a query of an unsearched batch has no hits")
        })
        .unwrap();
    if finish {
        writer.finish(&mut out).unwrap();
    }
    String::from_utf8(out).unwrap()
}

#[test]
fn an_unsearched_query_is_the_sample_in_16_and_15() {
    // No `hits` member, ` \nNo hits found` after the message (`&#xa;` in XML, `\u000A` in JSON),
    // the `Subject_N` list without a wrapper (XML) or as an array (JSON), `bl2seq` for -subject.
    assert_eq!(
        document(Xml2Kind::Single16, Xml2Target::Subjects(2), true),
        INVALID_ONLY_16
    );
    assert_eq!(
        document(Xml2Kind::Single15, Xml2Target::Subjects(2), true),
        INVALID_ONLY_15
    );
}

#[test]
fn an_unsearched_query_lists_no_masks() {
    // NCBI makes the result of a batch it did not search without the search, so the masks the
    // query carries (an N-only query with a lower-case run) are not written; the document is
    // the one of a query without masks, in 16 and 15.
    let message = invalid_message();
    let masks = [QueryMask {
        from: 10,
        to: 19,
        frame: 0,
    }];
    for (kind, expected) in [
        (Xml2Kind::Single16, INVALID_ONLY_16),
        (Xml2Kind::Single15, INVALID_ONLY_15),
    ] {
        let run = run(None, b"");
        let mut writer = Xml2Writer::new(kind, params(Xml2Target::Subjects(2)));
        let mut out = Vec::new();
        let query = XmlQuery {
            masks: &masks,
            ..invalid_query(&message)
        };
        writer
            .query(&run, &query, &mut out, |_| {
                panic!("a query of an unsearched batch has no hits")
            })
            .unwrap();
        writer.finish(&mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }
}

#[test]
fn a_subject_report_is_held_until_the_epilog() {
    // `-subject` keeps the last report back (blast_format.cpp:2363-2404): without the epilog
    // (a failed batch) the document ends before it.
    assert_eq!(
        document(Xml2Kind::Single16, Xml2Target::Subjects(2), false),
        INVALID_ONLY_16
            .split("\n<BlastOutput2>")
            .next()
            .unwrap()
            .to_string()
    );
    assert_eq!(
        document(Xml2Kind::Single15, Xml2Target::Subjects(2), false),
        "{\n\"BlastOutput2\": [\n"
    );
    // `-db` writes each report at once
    let db = document(Xml2Kind::Single16, Xml2Target::Database, false);
    assert!(db.ends_with("\n</BlastOutput2>\n"), "{db}");
}

#[test]
fn the_framing_joins_reports_and_ends_the_document() {
    let report = |text: &str| text.as_bytes().to_vec();
    // 16: the declaration and root start tag, the reports (each `\n<BlastOutput2>` ...
    // `</BlastOutput2>\n`, an empty line between two), the end tag.
    let mut doc = Xml2Document::new(Xml2Kind::Single16, false);
    let mut out = Vec::new();
    doc.report(&mut report("\n<BlastOutput2>a</BlastOutput2>\n"), &mut out)
        .unwrap();
    doc.report(&mut report("\n<BlastOutput2>b</BlastOutput2>\n"), &mut out)
        .unwrap();
    doc.epilog(&mut out).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!(
            "{}\n<BlastOutput2>a</BlastOutput2>\n\n<BlastOutput2>b</BlastOutput2>\n</BlastXML2>\n",
            framing::header16()
        )
    );
    // 15: `,\n` before the second report and later ones
    let mut doc = Xml2Document::new(Xml2Kind::Single15, false);
    let mut out = Vec::new();
    doc.report(&mut report("{\n  a\n}\n"), &mut out).unwrap();
    doc.report(&mut report("{\n  b\n}\n"), &mut out).unwrap();
    doc.epilog(&mut out).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "{\n\"BlastOutput2\": [\n{\n  a\n}\n,\n{\n  b\n}\n]\n}\n"
    );
    // a run that formatted no query: the prolog and the epilog (an empty pipe)
    for (kind, expected) in [
        (
            Xml2Kind::Single16,
            format!("{}</BlastXML2>\n", framing::header16()),
        ),
        (
            Xml2Kind::Single15,
            "{\n\"BlastOutput2\": [\n]\n}\n".to_string(),
        ),
    ] {
        let mut out = Vec::new();
        Xml2Document::new(kind, false).epilog(&mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }
    // the 16 header (xmlns attributes each on a line, `>` on its own line)
    assert!(framing::header16().starts_with(
        "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastXML2\n    xmlns=\"http://www.ncbi.nlm.nih.gov\"\n"
    ));
    assert!(framing::header16().ends_with("NCBI_BlastOutput2.xsd\"\n>"));
}

#[test]
fn a_document_started_early_writes_its_start_once() {
    // BLASTX writes the document's start with its prolog, before the first batch: the first
    // report and the epilog add no second one, and a run that stops at its first batch leaves
    // the start alone.
    for kind in [Xml2Kind::Single16, Xml2Kind::Single15] {
        for delayed in [false, true] {
            let lazy = |early: bool| {
                let mut doc = Xml2Document::new(kind, delayed);
                let mut out = Vec::new();
                if early {
                    doc.start(&mut out);
                }
                let start = out.len();
                doc.report(
                    &mut b"\n<BlastOutput2>a</BlastOutput2>\n".to_vec(),
                    &mut out,
                )
                .unwrap();
                doc.report(
                    &mut b"\n<BlastOutput2>b</BlastOutput2>\n".to_vec(),
                    &mut out,
                )
                .unwrap();
                doc.epilog(&mut out).unwrap();
                (start, out)
            };
            let (start, early) = lazy(true);
            assert_eq!(start, framing::prolog_bytes(kind).len(), "{kind:?}");
            assert_eq!(early, lazy(false).1, "{kind:?} delayed {delayed}");
            // a second call changes nothing
            let mut doc = Xml2Document::new(kind, delayed);
            let mut out = Vec::new();
            doc.start(&mut out);
            doc.start(&mut out);
            assert_eq!(out, framing::prolog_bytes(kind), "{kind:?}");
        }
    }
}

#[test]
fn a_delayed_document_writes_one_report_late_and_drops_the_last_on_a_failure() {
    let mut doc = Xml2Document::new(Xml2Kind::Single15, true);
    let mut out = Vec::new();
    doc.report(&mut b"{1}\n".to_vec(), &mut out).unwrap();
    assert_eq!(out, framing::HEAD15);
    doc.report(&mut b"{2}\n".to_vec(), &mut out).unwrap();
    assert_eq!(out, b"{\n\"BlastOutput2\": [\n{1}\n");
    // a failed batch: no epilog, `{2}` is never written
    let mut doc = Xml2Document::new(Xml2Kind::Single15, true);
    let mut out = Vec::new();
    doc.report(&mut b"{1}\n".to_vec(), &mut out).unwrap();
    doc.report(&mut b"{2}\n".to_vec(), &mut out).unwrap();
    doc.report(&mut b"{3}\n".to_vec(), &mut out).unwrap();
    assert_eq!(out, b"{\n\"BlastOutput2\": [\n{1}\n,\n{2}\n");
    // the epilog writes the held report
    doc.epilog(&mut out).unwrap();
    assert_eq!(out, b"{\n\"BlastOutput2\": [\n{1}\n,\n{2}\n,\n{3}\n]\n}\n");
}

#[test]
fn a_run_written_by_several_documents_in_turn_has_the_bytes_of_one() {
    // The BLASTP block search writes a run's reports with a document per block, handing the
    // framing state (the count, the start, the held report of a `-subject` run) on.
    for kind in [Xml2Kind::Single16, Xml2Kind::Single15] {
        for delayed in [false, true] {
            let reports: Vec<Vec<u8>> = (1..=5)
                .map(|n| format!("\n<BlastOutput2>{n}</BlastOutput2>\n").into_bytes())
                .collect();
            let mut whole = Xml2Document::new(kind, delayed);
            let mut expected = Vec::new();
            for report in &reports {
                whole.report(&mut report.clone(), &mut expected).unwrap();
            }
            whole.epilog(&mut expected).unwrap();
            for blocks in [vec![1, 4], vec![2, 2, 1], vec![1, 1, 1, 1, 1]] {
                let mut out = Vec::new();
                let mut carry = Xml2Carry::default();
                let mut next = 0;
                for (index, size) in blocks.iter().enumerate() {
                    let mut doc = Xml2Document::new(kind, delayed);
                    doc.resume(carry);
                    for report in &reports[next..next + size] {
                        doc.report(&mut report.clone(), &mut out).unwrap();
                    }
                    next += size;
                    if index + 1 == blocks.len() {
                        doc.epilog(&mut out).unwrap();
                    }
                    carry = doc.take_carry();
                }
                assert_eq!(
                    out, expected,
                    "{kind:?} delayed {delayed} blocks {blocks:?}"
                );
            }
        }
    }
}

fn one_hsp(hits_title: &'static [u8], target: Xml2Target, taxid: i32, kind: Xml2Kind) -> String {
    one_hsp_named(
        hits_title,
        target,
        taxid,
        kind,
        crate::report::blast_xml::records::UNNAMED,
    )
}

fn one_hsp_named(
    hits_title: &'static [u8],
    target: Xml2Target,
    taxid: i32,
    kind: Xml2Kind,
    taxonomy: &crate::report::blast_xml::HitTaxonomy,
) -> String {
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let run = run(
        Some(&matrix),
        if target == Xml2Target::Database {
            b"db/sprot"
        } else {
            b""
        },
    );
    let mut writer = Xml2Writer::new(kind, params(target));
    let masks = [
        QueryMask {
            from: 1,
            to: 2,
            frame: 0,
        },
        QueryMask {
            from: 5,
            to: 6,
            frame: -1,
        },
    ];
    let query = XmlQuery {
        id: b"Query_1",
        title: b"",
        len: 12,
        stats: XmlStats::NONE,
        message: None,
        searched: true,
        masks: &masks,
    };
    let mut out = Vec::new();
    writer
        .query(&run, &query, &mut out, |sink| {
            sink.hit(
                &XmlHit {
                    id: b"gnl|BL_ORD_ID|3",
                    title: hits_title,
                    accession: b"3",
                    len: 20,
                    taxid,
                    taxonomy,
                },
                |hsps| {
                    hsps.hsp(&XmlHsp {
                        score: 21,
                        bit_score: 12.697412345,
                        evalue: 1.515236,
                        query_from: 3,
                        query_to: 8,
                        hit_from: 5,
                        hit_to: 10,
                        query_frame: 0,
                        hit_frame: 0,
                        qrow: b"WCNYFW",
                        srow: b"WENYCW",
                        ops: None,
                        query_first: 2,
                        query_step: 1,
                    })
                },
            )
        })
        .unwrap();
    writer.finish(&mut out).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn hit_descriptions_follow_the_database_and_subject_rules() {
    // -db: the accession (the ordinal id), the stored title (empty: `<title></title>`), the taxid
    // only when non-zero
    let db = one_hsp(b"", Xml2Target::Database, 0, Xml2Kind::Single16);
    assert!(db.contains(
        "<HitDescr>\n                      <id>gnl|BL_ORD_ID|3</id>\n                      <accession>3</accession>\n                      <title></title>\n                    </HitDescr>"
    ), "{db}");
    assert!(db.contains("<db>db/sprot</db>"), "{db}");
    assert!(
        db.contains("<db-num>4</db-num>\n                  <db-len>5000</db-len>"),
        "{db}"
    );
    assert!(
        !db.contains("<query-title>"),
        "an empty title has no member"
    );
    let tax = one_hsp(b"t", Xml2Target::Database, 9606, Xml2Kind::Single16);
    assert!(
        tax.contains("<title>t</title>\n                      <taxid>9606</taxid>"),
        "{tax}"
    );
    // -subject: a local id has no accession and no taxid; the Subject_N list; no db-num
    let sub = one_hsp(
        b"s1 title",
        Xml2Target::Subjects(1),
        9606,
        Xml2Kind::Single16,
    );
    assert!(sub.contains(
        "<HitDescr>\n                      <id>gnl|BL_ORD_ID|3</id>\n                      <title>s1 title</title>\n                    </HitDescr>"
    ), "{sub}");
    assert!(
        sub.contains("<Target>\n          <subjects>Subject_1</subjects>\n        </Target>"),
        "{sub}"
    );
    assert!(!sub.contains("db-num"), "{sub}");
    // the same in JSON
    let json = one_hsp(b"t", Xml2Target::Database, 9606, Xml2Kind::Single15);
    assert!(json.contains("\"accession\": \"3\",\n                  \"title\": \"t\",\n                  \"taxid\": 9606"), "{json}");
}

#[test]
fn query_masking_lists_the_plus_regions_and_a_nucleotide_row_is_masked_in_lower_case() {
    let out = one_hsp(b"t", Xml2Target::Database, 0, Xml2Kind::Single16);
    // the region of a minus frame is not listed
    assert!(out.contains(
        "<query-masking>\n                <Range>\n                  <from>1</from>\n                  <to>2</to>\n                </Range>\n              </query-masking>"
    ), "{out}");
    assert!(!out.contains("<from>5</from>"), "{out}");
    // the Hsp has no strand and no frame for blastp; the query row has no mask in its range
    assert!(!out.contains("-strand") && !out.contains("-frame"), "{out}");
    assert!(out.contains("<qseq>WCNYFW</qseq>"), "{out}");
}

/// BLASTN's `Hsp` (OF5a part 2): `query-strand` Plus and a `hit-strand` from the order of the
/// subject's coordinates (a Minus HSP has `hit-from` above `hit-to`), no `positive`, no frame; the
/// regions of the query (frame 0) are applied to the nucleotide row in lower case, and the
/// scores come from `reward_penalty`.
fn blastn_document(kind: Xml2Kind, minus: bool) -> String {
    let mut run = run(None, b"");
    run.program = "blastn";
    run.version = "BLASTN 2.17.0+";
    run.matrix = None;
    run.reward_penalty = Some((1, -2));
    run.gap_open = 0;
    run.gap_extend = 0;
    run.filter = b"L;m;";
    let mut writer = Xml2Writer::new(kind, params(Xml2Target::Subjects(2)));
    let masks = [QueryMask {
        from: 3,
        to: 4,
        frame: 0,
    }];
    let query = XmlQuery {
        id: b"Query_1",
        title: b"",
        len: 12,
        stats: XmlStats::NONE,
        message: None,
        searched: true,
        masks: &masks,
    };
    let (hit_from, hit_to) = if minus { (10, 5) } else { (5, 10) };
    let mut out = Vec::new();
    writer
        .query(&run, &query, &mut out, |sink| {
            sink.hit(
                &XmlHit {
                    id: b"Subject_1",
                    title: b"",
                    accession: b"",
                    len: 20,
                    taxid: 0,
                    taxonomy: crate::report::blast_xml::records::UNNAMED,
                },
                |hsps| {
                    hsps.hsp(&XmlHsp {
                        score: 21,
                        bit_score: 12.7,
                        evalue: 1.5,
                        query_from: 3,
                        query_to: 8,
                        hit_from,
                        hit_to,
                        query_frame: 1,
                        hit_frame: if minus { -1 } else { 1 },
                        qrow: b"ACGTAC",
                        srow: b"ACGTAC",
                        ops: None,
                        query_first: 2,
                        query_step: 1,
                    })
                },
            )
        })
        .unwrap();
    writer.finish(&mut out).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn a_blastn_hsp_has_strands_no_positive_and_a_masked_row() {
    for minus in [false, true] {
        let out = blastn_document(Xml2Kind::Single16, minus);
        assert!(out.contains("<query-strand>Plus</query-strand>"), "{out}");
        let strand = if minus { "Minus" } else { "Plus" };
        assert!(
            out.contains(&format!("<hit-strand>{strand}</hit-strand>")),
            "{out}"
        );
        assert!(!out.contains("<positive>"), "{out}");
        assert!(!out.contains("-frame"), "{out}");
        // the row starts at the 0-based query position 2; the region 3 to 4 (inclusive) masks its
        // second and third letters
        assert!(out.contains("<qseq>AcgTAC</qseq>"), "{out}");
        assert!(out.contains("<sc-match>1</sc-match>"), "{out}");
        assert!(out.contains("<sc-mismatch>-2</sc-mismatch>"), "{out}");
        assert!(out.contains("<gap-open>0</gap-open>"), "{out}");
        assert!(out.contains("<filter>L;m;</filter>"), "{out}");
        // an empty subject title and no accession for a `-subject` hit
        assert!(out.contains("<title></title>"), "{out}");
        assert!(!out.contains("<accession>"), "{out}");
        // no description: no query-title
        assert!(!out.contains("<query-title>"), "{out}");
    }
}

#[test]
fn the_json_report_is_an_array_element_and_the_xml_std_mode_has_its_own_tags() {
    // JsonEncoder::array_element: the object at one level, no wrapper, a final end of line
    let mut out = Vec::new();
    let mut enc = JsonEncoder::array_element(&mut out);
    enc.begin_document("BlastOutput2");
    enc.begin_sequence("BlastOutput2");
    enc.member("report", 0);
    enc.variant("Target", "db", 0);
    enc.visible_str(b"d\"b\n");
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "{\n    \"report\": {\n      \"db\": \"d\\\"b\\u000A\"\n    }\n  }\n"
    );
    // the 14 side file: declaration, the root with the namespace attributes, member tags by
    // name only, a CHOICE as its type tag then the variant, an LF as `&#xa;`
    let mut out = Vec::new();
    let mut enc = XmlEncoder::xml2_document(&mut out);
    enc.begin_document("BlastOutput2");
    enc.begin_sequence("BlastOutput2");
    enc.member("report", 0);
    enc.begin_sequence("Report");
    enc.member("search-target", 3);
    enc.variant("Target", "db", 0);
    enc.visible_str(b"a\nb");
    enc.end_variant();
    enc.end_member();
    enc.member("hits", 4);
    enc.begin_list(false);
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!(
            "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastOutput2{}>\n  <report>\n    <Report>\n      <search-target>\n        <Target>\n          <db>a&#xa;b</db>\n        </Target>\n      </search-target>\n      <hits></hits>\n    </Report>\n  </report>\n</BlastOutput2>\n",
            crate::report::serial::xml::XML2_ROOT_ATTRIBUTES
        )
    );
}

#[test]
fn the_scientific_name_follows_the_taxid_and_an_unknown_taxid_leaves_the_id_alone() {
    use crate::report::blast_xml::HitTaxonomy;
    let named = HitTaxonomy::Named(b"Quote \"x\", y; z".to_vec());
    let xml = one_hsp_named(b"t", Xml2Target::Database, 9606, Xml2Kind::Single16, &named);
    assert!(
        xml.contains("<taxid>9606</taxid>\n                      <sciname>Quote &quot;x&quot;, y; z</sciname>\n                    </HitDescr>"),
        "{xml}"
    );
    let json = one_hsp_named(b"t", Xml2Target::Database, 9606, Xml2Kind::Single15, &named);
    assert!(
        json.contains("\"taxid\": 9606,\n                  \"sciname\": \"Quote \\\"x\\\", y; z\""),
        "{json}"
    );
    // an empty name is written
    let empty = HitTaxonomy::Named(Vec::new());
    let xml = one_hsp_named(b"t", Xml2Target::Database, 9606, Xml2Kind::Single16, &empty);
    assert!(xml.contains("<sciname></sciname>"), "{xml}");
    // not found: the description is the id alone, the hit keeps its length
    let missing = HitTaxonomy::NotFound(Vec::new());
    let xml = one_hsp_named(
        b"t",
        Xml2Target::Database,
        9606,
        Xml2Kind::Single16,
        &missing,
    );
    assert!(
        xml.contains("<HitDescr>\n                      <id>gnl|BL_ORD_ID|3</id>\n                    </HitDescr>"),
        "{xml}"
    );
    assert!(xml.contains("<len>20</len>"), "{xml}");
    // a name without a taxid is not written
    let xml = one_hsp_named(b"t", Xml2Target::Database, 0, Xml2Kind::Single16, &named);
    assert!(!xml.contains("sciname"), "{xml}");
}

// ---- Outfmt 13 and 14: the names of the files and the master file (OF5c) ----

/// Side files in memory: the files of a test, and `fail_at` the n-th one that cannot be opened.
#[derive(Default)]
struct MemoryFiles {
    names: Vec<u8>,
    files: std::rc::Rc<std::cell::RefCell<Vec<(u32, Vec<u8>)>>>,
    fail_at: Option<u32>,
}

struct Collector(u32, std::rc::Rc<std::cell::RefCell<Vec<(u32, Vec<u8>)>>>);

impl std::io::Write for Collector {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        let mut files = self.1.borrow_mut();
        let at = files.iter().position(|(n, _)| *n == self.0).unwrap();
        files[at].1.extend_from_slice(data);
        Ok(data.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl SideFiles for MemoryFiles {
    fn create(&mut self, n: u32) -> std::io::Result<Box<dyn std::io::Write + '_>> {
        if self.fail_at == Some(n) {
            return Err(std::io::Error::other("is a directory"));
        }
        self.files.borrow_mut().push((n, Vec::new()));
        Ok(Box::new(Collector(n, self.files.clone())))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        let mut name = self.names.clone();
        name.extend_from_slice(format!("_{n}.x").as_bytes());
        name
    }
}

/// The master files of NCBI BLAST+ 2.17.0 for two reports (probes `naming/upper14`, `naming/xml13`,
/// `blastp_sub_emptypipe_fmt1[34]`: `S/probe5/out`).
const MASTER14_HEAD_TEXT: &str = "<?xml version=\"1.0\"?>\n<BlastXML2\nxmlns=\"http://www.ncbi.nlm.nih.gov\"\nxmlns:xi=\"http://www.w3.org/2003/XInclude\"\nxmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\nxs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\">\n";

#[test]
fn the_master_files_list_the_side_files_and_the_json_master_has_no_final_end_of_line() {
    for (kind, reports, expected) in [
        (
            Xml2Kind::Files14,
            2,
            format!("{MASTER14_HEAD_TEXT}\t<xi:include href=\"b_1.x\"/>\n\t<xi:include href=\"b_2.x\"/>\n</BlastXML2>\n"),
        ),
        (
            Xml2Kind::Files14,
            0,
            format!("{MASTER14_HEAD_TEXT}</BlastXML2>\n"),
        ),
        (
            Xml2Kind::Files13,
            2,
            "{\n\t\"BlastJSON\": [\n\t\t{\"File\": \"b_1.x\" },\n\t\t{\"File\": \"b_2.x\" }\n\t]\n}".to_string(),
        ),
        (
            Xml2Kind::Files13,
            1,
            "{\n\t\"BlastJSON\": [\n\t\t{\"File\": \"b_1.x\" }\n\t]\n}".to_string(),
        ),
        (
            Xml2Kind::Files13,
            0,
            "{\n\t\"BlastJSON\": [\n\t]\n}".to_string(),
        ),
    ] {
        let files = MemoryFiles {
            names: b"b".to_vec(),
            ..MemoryFiles::default()
        };
        let stored = files.files.clone();
        let mut document = Xml2Document::with_files(kind, false, files);
        let mut master = Vec::new();
        for n in 1..=reports {
            document
                .report(&mut format!("report {n}\n").into_bytes(), &mut master)
                .unwrap();
        }
        // the reports go to the files, nothing to the master until the epilog
        assert!(master.is_empty(), "{kind:?}");
        document.epilog(&mut master).unwrap();
        assert_eq!(String::from_utf8(master).unwrap(), expected, "{kind:?} {reports}");
        let stored = stored.borrow();
        assert_eq!(stored.len(), reports);
        for (index, (n, bytes)) in stored.iter().enumerate() {
            assert_eq!(*n as usize, index + 1);
            assert_eq!(bytes, format!("report {n}\n").as_bytes());
        }
    }
}

#[test]
fn the_files_of_a_subject_run_are_written_one_report_late_and_a_lost_report_has_no_file() {
    // `-subject` holds the last report back: a failed batch (no epilog) leaves one file and no
    // master (probes `blastp_sub_badbatch_fmt1[34]`); `-db` writes at once (`blastp_db_badbatch_*`).
    for (delayed, files_written) in [(true, 1), (false, 2)] {
        let files = MemoryFiles::default();
        let stored = files.files.clone();
        let mut document = Xml2Document::with_files(Xml2Kind::Files14, delayed, files);
        let mut master = Vec::new();
        document.report(&mut b"q1\n".to_vec(), &mut master).unwrap();
        document.report(&mut b"q2\n".to_vec(), &mut master).unwrap();
        assert!(master.is_empty());
        assert_eq!(stored.borrow().len(), files_written, "delayed {delayed}");
    }
}

#[test]
fn a_side_file_that_cannot_be_opened_is_ncbis_argument_error_and_earlier_files_stay() {
    // probes `naming/unwritable13`, `unwritable14`: `Command line argument error: Cannot open
    // output file`, exit status 1; the first file is written, the master stays empty.
    let files = MemoryFiles {
        fail_at: Some(2),
        ..MemoryFiles::default()
    };
    let stored = files.files.clone();
    let mut document = Xml2Document::with_files(Xml2Kind::Files13, false, files);
    let mut master = Vec::new();
    document.report(&mut b"q1\n".to_vec(), &mut master).unwrap();
    let error = document
        .report(&mut b"q2\n".to_vec(), &mut master)
        .unwrap_err();
    let native = error.downcast::<crate::cli::NativeError>().unwrap();
    assert_eq!(native.exit, 1);
    assert_eq!(
        native.message,
        "Command line argument error: Cannot open output file\n"
    );
    assert_eq!(stored.borrow().len(), 1);
    assert!(master.is_empty());
}

#[test]
fn no_files_means_no_side_file_for_a_caller_that_has_no_output_file() {
    let mut document = Xml2Document::for_run(Xml2Kind::Files14, false, None);
    let error = document
        .report(&mut b"q1\n".to_vec(), &mut Vec::new())
        .unwrap_err();
    assert!(error.downcast_ref::<crate::cli::NativeError>().is_some());
    // outfmt 15 and 16 have no side files, with or without a name
    for kind in [Xml2Kind::Single16, Xml2Kind::Single15] {
        let mut out = Vec::new();
        let mut document = Xml2Document::for_run(kind, false, Some(std::path::Path::new("o.out")));
        document.report(&mut b"x\n".to_vec(), &mut out).unwrap();
        document.epilog(&mut out).unwrap();
        assert!(String::from_utf8(out).unwrap().contains("x\n"));
    }
}

/// The eight naming probes of `S/probe5/out/naming` (DESIGN_OF5 1.7): the `-out` argument, the
/// format, the side file `FsSideFiles` writes and the name the master lists. The extension is
/// dropped only if it is `.xml` (14) or `.json` (13) in any case; the master always uses the XML
/// rule and drops the directory.
#[test]
fn file_names_follow_ncbis_base_name_rules() {
    let cases = [
        ("d/X.XML", Xml2Kind::Files14, "d/X_1.xml", "X_1.xml"),
        ("d/Y.Json", Xml2Kind::Files13, "d/Y_1.json", "Y.Json_1.json"),
        ("d/k.xml", Xml2Kind::Files13, "d/k.xml_1.json", "k_1.json"),
        (
            "d/m.a.json",
            Xml2Kind::Files13,
            "d/m.a_1.json",
            "m.a.json_1.json",
        ),
        (
            "d/j.json",
            Xml2Kind::Files14,
            "d/j.json_1.xml",
            "j.json_1.xml",
        ),
        ("d/z", Xml2Kind::Files14, "d/z_1.xml", "z_1.xml"),
        ("d/.xml", Xml2Kind::Files14, "d/_1.xml", "_1.xml"),
        ("rel.xml", Xml2Kind::Files14, "rel_1.xml", "rel_1.xml"),
        // more spellings: the file name's last dot, a dot in the directory, a trailing dot
        ("a.b/c", Xml2Kind::Files14, "a.b/c_1.xml", "c_1.xml"),
        ("a.b/c.XmL", Xml2Kind::Files14, "a.b/c_1.xml", "c_1.xml"),
        ("c.", Xml2Kind::Files13, "c._1.json", "c._1.json"),
        (
            "/abs/o.JSON",
            Xml2Kind::Files13,
            "/abs/o_1.json",
            "o.JSON_1.json",
        ),
    ];
    for (out, kind, side, listed) in cases {
        let files = FsSideFiles::new(std::path::Path::new(out), kind);
        assert_eq!(
            files.path(1),
            std::path::PathBuf::from(side),
            "{out} {kind:?}"
        );
        assert_eq!(
            String::from_utf8(files.listed(1)).unwrap(),
            listed,
            "{out} {kind:?}"
        );
    }
    assert_eq!(
        base_name(std::path::Path::new("d/m.a.json"), false, true),
        b"d/m.a"
    );
    assert_eq!(
        base_name(std::path::Path::new("d/m.a.json"), true, false),
        b"m.a.json"
    );
    assert_eq!(base_name(std::path::Path::new("x.xml"), true, false), b"x");
}

#[test]
fn side_files_on_the_file_system_are_truncated_and_numbered() {
    let dir = std::env::temp_dir().join(format!("losat_xml2_files_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // a stale file of an earlier run: a longer file is cut to the new report (`IOS_BASE::out`)
    std::fs::write(dir.join("o.out_1.xml"), "stale stale stale stale\n").unwrap();
    std::fs::write(dir.join("o.out_3.xml"), "stale\n").unwrap();
    let out = dir.join("o.out");
    let mut document = Xml2Document::for_run(Xml2Kind::Files14, false, Some(&out));
    let mut master = Vec::new();
    document
        .report(&mut b"one\n".to_vec(), &mut master)
        .unwrap();
    document
        .report(&mut b"two\n".to_vec(), &mut master)
        .unwrap();
    document.epilog(&mut master).unwrap();
    assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one\n");
    assert_eq!(std::fs::read(dir.join("o.out_2.xml")).unwrap(), b"two\n");
    // the file of a number the run did not reach stays as it was (probe `stale14`)
    assert_eq!(std::fs::read(dir.join("o.out_3.xml")).unwrap(), b"stale\n");
    let master = String::from_utf8(master).unwrap();
    assert!(master.contains("href=\"o.out_1.xml\""), "{master}");
    assert!(master.contains("href=\"o.out_2.xml\""), "{master}");
    // a directory of the name of a side file: NCBI's error
    let blocked = dir.join("p.out");
    std::fs::create_dir_all(dir.join("p.out_2.xml")).unwrap();
    let mut document = Xml2Document::for_run(Xml2Kind::Files14, false, Some(&blocked));
    document
        .report(&mut b"one\n".to_vec(), &mut Vec::new())
        .unwrap();
    let error = document
        .report(&mut b"two\n".to_vec(), &mut Vec::new())
        .unwrap_err();
    assert!(error.downcast_ref::<crate::cli::NativeError>().is_some());
    assert_eq!(std::fs::read(dir.join("p.out_1.xml")).unwrap(), b"one\n");
    std::fs::remove_dir_all(&dir).unwrap();
}

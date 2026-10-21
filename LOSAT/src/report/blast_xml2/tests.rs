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

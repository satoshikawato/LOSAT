//! Unit tests of the outfmt 5 writer: the pieces of NCBI BLAST+ 2.17.0's documents
//! (`ncbi-samples/out/blastp_subject_fmt5.out`, probe `blastp_sub_nohfirst`) written from
//! records.

use super::matrix::FormattingMatrix;
use super::*;
use crate::report::serial::encoder::Encoder;
use crate::report::serial::xml::XmlEncoder;

fn run<'a>(matrix: &'a FormattingMatrix) -> XmlRun<'a> {
    XmlRun {
        program: "blastp",
        app: "blastp",
        version: "BLASTP 2.17.0+",
        reference: REFERENCE_GAPPED_BLAST,
        db_name: b"",
        matrix: Some("BLOSUM62"),
        expect: 10.0,
        gap_open: 11,
        gap_extend: 1,
        reward_penalty: None,
        filter: b"F",
        pattern: None,
        db_num: 0,
        db_len: 0,
        db_len_overflow: false,
        positives: Some(matrix),
        mask_char: b'X',
    }
}

const STATS: XmlStats = XmlStats {
    hsp_len: 11,
    eff_space: 8738,
    kappa: 0.041,
    lambda: 0.267,
    entropy: 0.14,
};

const HEADER: &str = concat!(
    "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n",
    "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n",
    "<BlastOutput>\n",
    "  <BlastOutput_program>blastp</BlastOutput_program>\n",
    "  <BlastOutput_version>BLASTP 2.17.0+</BlastOutput_version>\n",
    "  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>\n",
    "  <BlastOutput_db></BlastOutput_db>\n",
    "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>\n",
    "  <BlastOutput_query-def>noh first query</BlastOutput_query-def>\n",
    "  <BlastOutput_query-len>32</BlastOutput_query-len>\n",
    "  <BlastOutput_param>\n",
    "    <Parameters>\n",
    "      <Parameters_matrix>BLOSUM62</Parameters_matrix>\n",
    "      <Parameters_expect>10</Parameters_expect>\n",
    "      <Parameters_gap-open>11</Parameters_gap-open>\n",
    "      <Parameters_gap-extend>1</Parameters_gap-extend>\n",
    "      <Parameters_filter>F</Parameters_filter>\n",
    "    </Parameters>\n",
    "  </BlastOutput_param>\n",
    "<BlastOutput_iterations>"
);

/// Iteration 1 of probe `blastp_sub_nohfirst`: no hits.
const NO_HITS: &str = concat!(
    "\n<Iteration>\n",
    "  <Iteration_iter-num>1</Iteration_iter-num>\n",
    "  <Iteration_query-ID>Query_1</Iteration_query-ID>\n",
    "  <Iteration_query-def>noh first query</Iteration_query-def>\n",
    "  <Iteration_query-len>32</Iteration_query-len>\n",
    "<Iteration_hits>\n",
    "</Iteration_hits>\n",
    "  <Iteration_stat>\n",
    "    <Statistics>\n",
    "      <Statistics_db-num>0</Statistics_db-num>\n",
    "      <Statistics_db-len>0</Statistics_db-len>\n",
    "      <Statistics_hsp-len>11</Statistics_hsp-len>\n",
    "      <Statistics_eff-space>8738</Statistics_eff-space>\n",
    "      <Statistics_kappa>0.041</Statistics_kappa>\n",
    "      <Statistics_lambda>0.267</Statistics_lambda>\n",
    "      <Statistics_entropy>0.14</Statistics_entropy>\n",
    "    </Statistics>\n",
    "  </Iteration_stat>\n",
    "  <Iteration_message>No hits found</Iteration_message>\n",
    "</Iteration>"
);

/// Iteration 2 of `blastp_subject_fmt5.out` (query 2) with its Hit of two HSPs.
const WITH_HITS: &str = concat!(
    "\n<Iteration>\n",
    "  <Iteration_iter-num>2</Iteration_iter-num>\n",
    "  <Iteration_query-ID>Query_2</Iteration_query-ID>\n",
    "  <Iteration_query-def>qp2 unrelated protein (no hits expected)</Iteration_query-def>\n",
    "  <Iteration_query-len>45</Iteration_query-len>\n",
    "<Iteration_hits>\n",
    "<Hit>\n",
    "  <Hit_num>1</Hit_num>\n",
    "  <Hit_id>sp1</Hit_id>\n",
    "  <Hit_def>sp1 subject protein alpha</Hit_def>\n",
    "  <Hit_accession>Subject_1</Hit_accession>\n",
    "  <Hit_len>137</Hit_len>\n",
    "  <Hit_hsps>\n",
    "    <Hsp>\n",
    "      <Hsp_num>1</Hsp_num>\n",
    "      <Hsp_bit-score>12.6974</Hsp_bit-score>\n",
    "      <Hsp_score>21</Hsp_score>\n",
    "      <Hsp_evalue>1.51524</Hsp_evalue>\n",
    "      <Hsp_query-from>33</Hsp_query-from>\n",
    "      <Hsp_query-to>38</Hsp_query-to>\n",
    "      <Hsp_hit-from>35</Hsp_hit-from>\n",
    "      <Hsp_hit-to>40</Hsp_hit-to>\n",
    "      <Hsp_query-frame>0</Hsp_query-frame>\n",
    "      <Hsp_hit-frame>0</Hsp_hit-frame>\n",
    "      <Hsp_identity>4</Hsp_identity>\n",
    "      <Hsp_positive>4</Hsp_positive>\n",
    "      <Hsp_gaps>0</Hsp_gaps>\n",
    "      <Hsp_align-len>6</Hsp_align-len>\n",
    "      <Hsp_qseq>WCNYFW</Hsp_qseq>\n",
    "      <Hsp_hseq>WENYCW</Hsp_hseq>\n",
    "      <Hsp_midline>W NY W</Hsp_midline>\n",
    "    </Hsp>\n",
    "    <Hsp>\n",
    "      <Hsp_num>2</Hsp_num>\n",
    "      <Hsp_bit-score>12.3122</Hsp_bit-score>\n",
    "      <Hsp_score>20</Hsp_score>\n",
    "      <Hsp_evalue>1.88522</Hsp_evalue>\n",
    "      <Hsp_query-from>24</Hsp_query-from>\n",
    "      <Hsp_query-to>32</Hsp_query-to>\n",
    "      <Hsp_hit-from>1</Hsp_hit-from>\n",
    "      <Hsp_hit-to>9</Hsp_hit-to>\n",
    "      <Hsp_query-frame>0</Hsp_query-frame>\n",
    "      <Hsp_hit-frame>0</Hsp_hit-frame>\n",
    "      <Hsp_identity>3</Hsp_identity>\n",
    "      <Hsp_positive>5</Hsp_positive>\n",
    "      <Hsp_gaps>0</Hsp_gaps>\n",
    "      <Hsp_align-len>9</Hsp_align-len>\n",
    "      <Hsp_qseq>LWKEHLTFQ</Hsp_qseq>\n",
    "      <Hsp_hseq>MWMEHSVLK</Hsp_hseq>\n",
    "      <Hsp_midline>+W EH   +</Hsp_midline>\n",
    "    </Hsp>\n",
    "  </Hit_hsps>\n",
    "</Hit>\n",
    "</Iteration_hits>\n",
    "  <Iteration_stat>\n",
    "    <Statistics>\n",
    "      <Statistics_db-num>0</Statistics_db-num>\n",
    "      <Statistics_db-len>0</Statistics_db-len>\n",
    "      <Statistics_hsp-len>11</Statistics_hsp-len>\n",
    "      <Statistics_eff-space>8738</Statistics_eff-space>\n",
    "      <Statistics_kappa>0.041</Statistics_kappa>\n",
    "      <Statistics_lambda>0.267</Statistics_lambda>\n",
    "      <Statistics_entropy>0.14</Statistics_entropy>\n",
    "    </Statistics>\n",
    "  </Iteration_stat>\n",
    "</Iteration>"
);

const TRAILER: &str = "\n</BlastOutput_iterations>\n</BlastOutput>\n\n";

fn write_hits(writer: &mut BlastXmlWriter, matrix: &FormattingMatrix, out: &mut Vec<u8>) {
    let run = run(matrix);
    let first = XmlQuery {
        id: b"Query_1",
        title: b"noh first query",
        len: 32,
        stats: STATS,
        message: Some(b"No hits found"),
        searched: true,
        masks: &[],
        joins_held: false,
    };
    writer.query(&run, &first, out, |_| Ok(())).unwrap();
    let second = XmlQuery {
        id: b"Query_2",
        title: b"qp2 unrelated protein (no hits expected)",
        len: 45,
        stats: STATS,
        message: None,
        searched: true,
        masks: &[],
        joins_held: false,
    };
    writer
        .query(&run, &second, out, |sink| {
            sink.hit(
                &XmlHit {
                    id: b"sp1",
                    title: b"sp1 subject protein alpha",
                    accession: b"Subject_1",
                    len: 137,
                    taxid: 0,
                    taxonomy: crate::report::blast_xml::records::UNNAMED,
                    descriptions: &[],
                    defline_thrown: false,
                    unloaded: false,
                    lookup_throws: None,
                },
                |hsps| {
                    hsps.hsp(&XmlHsp {
                        score: 21,
                        bit_score: 12.697412345,
                        evalue: 1.515236,
                        query_from: 33,
                        query_to: 38,
                        hit_from: 35,
                        hit_to: 40,
                        query_frame: 0,
                        hit_frame: 0,
                        qrow: b"WCNYFW",
                        srow: b"WENYCW",
                        ops: None,
                        query_first: 32,
                        query_step: 1,
                        alignment_identities: 0,
                    })?;
                    hsps.hsp(&XmlHsp {
                        score: 20,
                        bit_score: 12.3122,
                        evalue: 1.88522,
                        query_from: 24,
                        query_to: 32,
                        hit_from: 1,
                        hit_to: 9,
                        query_frame: 0,
                        hit_frame: 0,
                        qrow: b"LWKEHLTFQ",
                        srow: b"MWMEHSVLK",
                        ops: None,
                        query_first: 23,
                        query_step: 1,
                        alignment_identities: 0,
                    })
                },
            )
        })
        .unwrap();
}

#[test]
fn a_document_is_ncbis_pieces_joined() {
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let mut out = Vec::new();
    let mut writer = BlastXmlWriter::new(0);
    write_hits(&mut writer, &matrix, &mut out);
    writer.finish(&mut out).unwrap();
    let expected = format!("{HEADER}{NO_HITS}{WITH_HITS}{TRAILER}");
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn a_later_call_continues_the_document() {
    // The first call formatted one query; the second writes the iteration number 2 with no
    // header, and `finish` writes the trailer once.
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let mut out = Vec::new();
    let mut writer = BlastXmlWriter::new(1);
    let run = run(&matrix);
    let query = XmlQuery {
        id: b"Query_2",
        title: b"qp2 unrelated protein (no hits expected)",
        len: 45,
        stats: STATS,
        message: None,
        searched: true,
        masks: &[],
        joins_held: false,
    };
    writer.query(&run, &query, &mut out, |_| Ok(())).unwrap();
    let text = String::from_utf8(out.clone()).unwrap();
    assert!(text.starts_with("\n<Iteration>\n  <Iteration_iter-num>2</Iteration_iter-num>"));
    assert!(!text.contains("<BlastOutput>"));
    writer.finish(&mut out).unwrap();
    assert!(out.ends_with(TRAILER.as_bytes()));
    // an epilog call for a run with no query writes the epilog's end of line alone
    // (blast_format.cpp:2242-2243)
    let mut none = Vec::new();
    BlastXmlWriter::new(0).finish(&mut none).unwrap();
    assert_eq!(none, b"\n");
}

#[test]
fn blastn_parameters_and_the_statistics_without_a_valid_context() {
    let mut run = run_blastn();
    run.db_num = 3;
    run.db_len = 10_000_000_000;
    let mut out = Vec::new();
    let mut writer = BlastXmlWriter::new(0);
    let query = XmlQuery {
        id: b"Query_1",
        // an empty generated title is written as `No definition line`
        title: b"",
        len: 0,
        stats: XmlStats::NONE,
        message: Some(b"x"),
        searched: true,
        masks: &[],
        joins_held: false,
    };
    writer.query(&run, &query, &mut out, |_| Ok(())).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("      <Parameters_expect>0.001</Parameters_expect>\n      <Parameters_sc-match>1</Parameters_sc-match>\n      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>\n      <Parameters_gap-open>0</Parameters_gap-open>"));
    assert!(!text.contains("Parameters_matrix"));
    assert!(text.contains("<Statistics_db-num>3</Statistics_db-num>\n      <Statistics_db-len>10000000000</Statistics_db-len>\n      <Statistics_hsp-len>0</Statistics_hsp-len>\n      <Statistics_eff-space>0</Statistics_eff-space>\n      <Statistics_kappa>-1</Statistics_kappa>\n      <Statistics_lambda>-1</Statistics_lambda>\n      <Statistics_entropy>-1</Statistics_entropy>"));
}

// PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits (review C of line RG
// D1a-4, SF1): a summed
// `LENGTH` of exactly -2^63 (an alias file's `LENGTH` 2^63), which NCBI's serial XML writes as
// `-0` (probe q13 of review C), is written as its digits.
#[test]
fn a_db_len_of_minus_two_to_the_63_is_written_as_its_digits() {
    let mut run = run_blastn();
    run.db_num = 3;
    run.db_len = i64::MIN;
    let mut out = Vec::new();
    let query = XmlQuery {
        id: b"Query_1",
        title: b"",
        len: 0,
        stats: XmlStats::NONE,
        message: Some(b"x"),
        searched: true,
        masks: &[],
        joins_held: false,
    };
    BlastXmlWriter::new(0)
        .query(&run, &query, &mut out, |_| Ok(()))
        .unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(
        text.contains("<Statistics_db-len>-9223372036854775808</Statistics_db-len>"),
        "{text}"
    );
}

fn run_blastn() -> XmlRun<'static> {
    XmlRun {
        program: "blastn",
        app: "blastn",
        version: "BLASTN 2.17.0+",
        reference: REFERENCE_MEGABLAST,
        db_name: b"db/snucl",
        matrix: None,
        expect: 0.001,
        gap_open: 0,
        gap_extend: 0,
        reward_penalty: Some((1, -2)),
        filter: b"L;m;",
        pattern: None,
        db_num: 0,
        db_len: 0,
        db_len_overflow: false,
        positives: None,
        mask_char: b'N',
    }
}

#[test]
fn strings_are_written_as_ncbis_us_ascii_stream_writes_them() {
    use crate::ncbi_str::utf8::Utf8Error;
    use crate::report::serial::encoder::{SerialError, StreamPosition};
    let def = |value: &[u8]| {
        let mut out = Vec::new();
        let mut enc = XmlEncoder::object(&mut out);
        enc.begin_sequence("Hit");
        enc.member("def", 2);
        enc.visible_str(value);
        enc.end_member();
        enc.end_sequence();
        let error = enc.error();
        let posts = enc.take_posts();
        let thrown = enc.take_thrown();
        (out, error, posts, thrown)
    };
    let (out, error, posts, _) = def(b"q2 & <x> \"quoted\" 'apos' a&b\t\x7f");
    assert_eq!((error, posts.len()), (None, 0));
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "\n<Hit>\n  <Hit_def>q2 &amp; &lt;x&gt; &quot;quoted&quot; &apos;apos&apos; a&amp;b&#x9;\x7f</Hit_def>\n</Hit>"
    );
    // Overlong forms of ASCII characters are written as those characters (K4h probe q.mq
    // classes ovl3, ovl3lt, ovl3c1 and taxonomy names ovl2, c080); a NUL ends the string.
    let (out, error, posts, _) = def(b"a\xe0\x81\x81b\xe0\x80\xbc\xc0\x80c\xc1\x81\x00d");
    assert_eq!(error, None);
    assert!(posts.is_empty());
    assert!(
        out.ends_with(b"<Hit_def>aAb&lt;cA</Hit_def>\n</Hit>"),
        "{out:?}"
    );
    // A control character, raw or decoded, is `#` with a post at the stream's line.
    let (out, error, posts, _) = def(b"\x01s1 a\xe0\x80\x81b");
    assert_eq!(error, None);
    assert!(out.ends_with(b"<Hit_def>#s1 a#b</Hit_def>\n</Hit>"));
    assert_eq!(posts.len(), 2);
    assert_eq!(posts[0].replaced.trace, "Hit.def");
    assert_eq!(posts[0].replaced.position, StreamPosition::Line(3));
    assert_eq!(posts[1].replaced.byte, 1);
    assert_eq!(
        posts[0].replaced.post(Some("blastp")),
        b"Error: (803.7) [blastp] Hit.def\nBad char [0x1] in string at line 3\n".to_vec()
    );
    // The three exceptions: the document is cut before the failing character.
    for (value, expected, kept) in [
        (
            &b"q1 a\xc3\xa9b z"[..],
            Utf8Error::SymbolToChar,
            &b"q1 a"[..],
        ),
        (b"S\xa9s", Utf8Error::DecodeFirst, b"S"),
        (b"q1 one\xc3", Utf8Error::DecodeNext, b"q1 one"),
    ] {
        let (out, error, posts, thrown) = def(value);
        assert_eq!(error, Some(SerialError::Utf8(expected)));
        assert!(posts.is_empty());
        let thrown = thrown.expect("thrown");
        let cut = [b"\n<Hit>\n  <Hit_def>".as_slice(), kept].concat();
        assert_eq!(&out[..thrown.at], cut.as_slice());
    }
    // The trace of an object written with `WriteObject` has no named frame.
    let (_, _, _, thrown) = def(b"a\xe9");
    let what = String::from_utf8(thrown.expect("thrown").exception().what()).unwrap();
    // NCBI's report (its header, then the three exceptions of the chain) in LOSAT's words
    // (`crate::wording::losat_text`).
    assert_eq!(what.lines().count(), 3, "{what}");
    assert_eq!(
        what,
        crate::wording::losat_text(
            "NCBI C++ Exception:\n\
             \x20   T0 \"/src/corelib/ncbistr.cpp\", line 7202: Error: (CStringException::eBadArgs) ncbi::CUtf8::DecodeNext() - Source string is not in UTF8 format (m_Pos = 0)\n\
             \x20   T0 \"/src/serial/objostrxml.cpp\", line 1538: Error: (CStringException::eBadArgs) ncbi::CObjectOStreamXml::WriteClassMember() -  Frame type= eFrameClassMember, Member name= def (m_Pos = 0)\n\
             \x20   T0 \"/src/serial/objostrxml.cpp\", line 1519: Error: (CStringException::eBadArgs) ncbi::CObjectOStreamXml::WriteClass() -  Frame type= eFrameClass, Object type= Hit (m_Pos = 0)\n"
        )
    );
}

#[test]
fn a_standard_xml_trace_names_the_choice_frames() {
    use crate::report::serial::encoder::StreamPosition;
    let mut out = Vec::new();
    let mut enc = XmlEncoder::xml2_report(&mut out);
    enc.begin_document("BlastOutput2");
    enc.begin_sequence("BlastOutput2");
    enc.member("report", 0);
    enc.begin_sequence("Report");
    enc.member("search-target", 3);
    enc.variant("Target", "db", 0);
    enc.visible_str(b"d\x01b");
    enc.end_variant();
    enc.end_member();
    enc.member("results", 5);
    enc.variant("Results", "search", 1);
    enc.begin_sequence("Search");
    enc.member("query-title", 1);
    enc.visible_str(b"q1 a\xc3\xa9");
    let posts = enc.take_posts();
    assert_eq!(
        posts[0].replaced.trace,
        "BlastOutput2.report.search-target.db"
    );
    assert_eq!(posts[0].replaced.position, StreamPosition::Line(7));
    let exception = enc.take_thrown().expect("thrown").exception();
    assert_eq!(
        exception.get_msg(),
        b" Frame type= eFrameNamed, Object type= BlastOutput2"
    );
    // The frames of NCBI's report, in LOSAT's words (`crate::wording::losat_text`): the lines
    // after the first exception's, each with the indentation of a further exception.
    let what = String::from_utf8(exception.what()).unwrap();
    let frames: Vec<&str> = what.lines().skip(1).map(str::trim_start).collect();
    assert_eq!(
        frames,
        [
            "Frame type= eFrameClassMember, Member name= query-title (m_Pos = 0)",
            "Frame type= eFrameClass, Object type= Search (m_Pos = 0)",
            "Frame type= eFrameChoiceVariant, Member name= search (m_Pos = 0)",
            "Frame type= eFrameChoice, Object type= Results (m_Pos = 0)",
            "Frame type= eFrameClassMember, Member name= results (m_Pos = 0)",
            "Frame type= eFrameClass, Object type= Report (m_Pos = 0)",
            "Frame type= eFrameClassMember, Member name= report (m_Pos = 0)",
            "Frame type= eFrameClass, Object type= BlastOutput2 (m_Pos = 0)",
            "Frame type= eFrameNamed, Object type= BlastOutput2 (m_Pos = 0)",
        ]
        .map(|frame| crate::wording::losat_text(frame).into_owned())
    );
}

#[test]
fn reals_are_written_with_fifteen_digits() {
    let mut out = Vec::new();
    let mut enc = XmlEncoder::object(&mut out);
    enc.begin_sequence("Statistics");
    for (name, value) in [
        ("kappa", 0.710602795216363),
        ("lambda", 1.37406312246009),
        ("entropy", 1.30724660390929),
        ("eff-space", 1e15),
        ("a", 1e-5),
        ("b", 0.001),
        ("c", 28680.0),
        ("d", 2.60916e-178),
        ("e", 0.133956144488482),
    ] {
        enc.member(name, 0);
        enc.real(value);
        enc.end_member();
    }
    enc.end_sequence();
    let text = String::from_utf8(out).unwrap();
    for line in [
        "<Statistics_kappa>0.710602795216363</Statistics_kappa>",
        "<Statistics_lambda>1.37406312246009</Statistics_lambda>",
        "<Statistics_entropy>1.30724660390929</Statistics_entropy>",
        "<Statistics_eff-space>1e+15</Statistics_eff-space>",
        "<Statistics_a>1e-05</Statistics_a>",
        "<Statistics_b>0.001</Statistics_b>",
        "<Statistics_c>28680</Statistics_c>",
        "<Statistics_d>2.60916e-178</Statistics_d>",
        "<Statistics_e>0.133956144488482</Statistics_e>",
    ] {
        assert!(text.contains(line), "{line} in {text}");
    }
}

#[test]
fn the_translation_frame_is_the_records_not_the_searched_intervals() {
    use super::hsp::translation_frame;
    // plus: (start - 1) % 3 + 1
    assert_eq!(translation_frame(true, 1, 30, 100), 1);
    assert_eq!(translation_frame(true, 5, 31, 100), 2);
    assert_eq!(translation_frame(true, 3, 29, 100), 3);
    // minus: -((length - end) % 3 + 1)
    assert_eq!(translation_frame(false, 10, 100, 100), -1);
    assert_eq!(translation_frame(false, 10, 99, 100), -2);
    assert_eq!(translation_frame(false, 10, 98, 100), -3);
    assert_eq!(translation_frame(false, 10, 97, 100), -1);
}

/// The header of the first query (BLAST XML) that throws (NCBI 2.17.0, K4f probe
/// `blastp.db.f5.none.b3`: the posts, then `Error: NCBI C++ Exception:`, exit 255) writes nothing
/// and makes no `Hit`; its 803.7 post goes to the diagnostics before the exception, which the
/// writer of the posts hands back alone (a library caller gets the post on its sink, review 1 of
/// K4h L5), and an `Iteration` that throws later cuts its query's loads.
#[test]
fn a_header_exception_comes_back_after_its_posts() {
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let run = run(&matrix);
    let query = |title: &'static [u8]| XmlQuery {
        id: b"Query_1",
        title,
        len: 32,
        stats: STATS,
        message: None,
        searched: true,
        masks: &[],
        joins_held: false,
    };
    let mut writer = BlastXmlWriter::new(0);
    let mut out = Vec::new();
    let mut hits_made = false;
    writer
        .query(&run, &query(b"q\x01 caf\xe9"), &mut out, |_| {
            hits_made = true;
            Ok(())
        })
        .unwrap();
    assert!(out.is_empty() && !hits_made);
    let mut posts = crate::report::query_warnings::HitPosts::default();
    posts.document(writer.take_posts());
    let mut diagnostics = Vec::new();
    let error = crate::report::query_warnings::write_posts_to(
        &out,
        &posts,
        0,
        &mut Vec::new(),
        &mut diagnostics,
        None,
    )
    .unwrap_err();
    let native = error.downcast_ref::<crate::cli::NativeError>().unwrap();
    assert_eq!(native.exit, 255);
    // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
    assert!(
        native.message.starts_with(crate::wording::losat_text(
            "Error: NCBI C++ Exception:\n    T0 \"/src/corelib/ncbistr.cpp\", line 7190: Error: (CStringException::eBadArgs) ncbi::CUtf8::DecodeFirst() - Source string is not in UTF8 format (m_Pos = 0)\n"
        ).as_bytes()),
        "{}",
        String::from_utf8_lossy(&native.message)
    );
    let diagnostics = String::from_utf8_lossy(&diagnostics);
    assert!(
        diagnostics.starts_with("Error: (803.7) [blastp] BlastOutput.query-def\n"),
        "{diagnostics}"
    );
    // A later query's `Iteration` that throws: the run goes on, its loads are cut.
    let mut writer = BlastXmlWriter::new(0);
    let mut out = Vec::new();
    writer
        .query(&run, &query(b"plain"), &mut out, |_| Ok(()))
        .unwrap();
    assert!(!writer.take_posts().cut);
    writer
        .query(&run, &query(b"caf\xe9"), &mut out, |_| {
            panic!("no `Hit` of a thrown `Iteration`")
        })
        .unwrap();
    let posts = writer.take_posts();
    assert!(posts.cut && posts.ends_run.is_none() && posts.posts.len() == 1);
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:126
// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:100
/// A run of another program's search strategy writes the program of the strategy into its
/// document (`BlastOutput_program`) and the running application's name into its posts: the
/// 803.7 post of a header string and the post of an exception that a report raised
/// (`caught_post`), both written by the application's process.
#[test]
fn the_posts_of_a_cross_program_run_name_the_running_application() {
    use crate::ncbi_exception::{NcbiException, ThrowSite};
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let run = XmlRun {
        program: "tblastx",
        app: "blastn",
        ..run(&matrix)
    };
    let query = XmlQuery {
        id: b"Query_1",
        title: b"q\x01 caf\xe9",
        len: 32,
        stats: STATS,
        message: None,
        searched: true,
        masks: &[],
        joins_held: false,
    };
    let mut writer = BlastXmlWriter::new(0);
    let mut out = Vec::new();
    writer.query(&run, &query, &mut out, |_| Ok(())).unwrap();
    let mut posts = crate::report::query_warnings::HitPosts::default();
    posts.document(writer.take_posts());
    let mut diagnostics = Vec::new();
    crate::report::query_warnings::write_posts_to(
        &out,
        &posts,
        0,
        &mut Vec::new(),
        &mut diagnostics,
        None,
    )
    .unwrap_err();
    let diagnostics = String::from_utf8_lossy(&diagnostics);
    assert!(
        diagnostics.starts_with("Error: (803.7) [blastn] BlastOutput.query-def\n"),
        "{diagnostics}"
    );
    let exception = NcbiException::new(
        ThrowSite {
            file: "src/objects/seq/../seqalign/Seq_align.cpp",
            line: 227,
            kind: "CSeqalignException::eInvalidAlignment",
            function: None,
        },
        "CSeq_align::GetSeqRange(): Row is empty",
    );
    let post = super::writer::caught_post(Err(CaughtException(exception).into()), &writer, run.app)
        .unwrap()
        .expect("the post of a caught exception");
    assert!(post.starts_with(b"Error: [blastn] "), "{post:?}");
    // Without an exception that the report function catches there is no post, and another error
    // goes back to the caller.
    assert!(super::writer::caught_post(Ok(()), &writer, run.app)
        .unwrap()
        .is_none());
    assert!(super::writer::caught_post(Err(anyhow::anyhow!("other")), &writer, run.app).is_err());
}

// Line DBID pfoff-4 (DECISIONS 62 (15); NCBI 2.17.0, probe `R/pfoff/p4`, `bp.U2.b1.dflt.5`): the
// `Hit` of a database sequence whose lookup throws where its batch's prefetch is off is NCBI's
// `catch`: the Seq-align id's FASTA form, `Unknown`, `Hit_len` 0, and each HSP with its number, the
// Seq-align's scores, the four coordinates 0, its `num_ident` (-1 without one) and empty rows; the
// writer reads no alignment of its HSPs, and the next `Hit` is written as usual. The XML2 and JSON
// documents raise the exception there instead (`CaughtException`), before anything of the `Hit`.
#[test]
fn a_hit_whose_lookup_throws_is_ncbis_catch() {
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let run = run(&matrix);
    let exception =
        crate::report::db_subject::LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO)
            .exception()
            .unwrap();
    let thrown = XmlHit {
        id: b"ref|NP_100001.1|",
        title: b"Unknown",
        accession: b"Unknown",
        len: 140,
        taxid: 0,
        taxonomy: crate::report::blast_xml::records::UNNAMED,
        descriptions: &[],
        defline_thrown: false,
        unloaded: false,
        lookup_throws: Some(&exception),
    };
    let query = XmlQuery {
        id: b"Query_2",
        title: b"qp2",
        len: 150,
        stats: STATS,
        message: None,
        searched: true,
        masks: &[],
        joins_held: false,
    };
    let mut writer = BlastXmlWriter::new(1);
    let mut out = Vec::new();
    writer
        .query(&run, &query, &mut out, |sink| {
            sink.hit(&thrown, |hsps| {
                assert!(!hsps.reads_alignment());
                hsps.hsp(&XmlHsp::unaligned(351, 139.8134567, 1.145634e-47, 66))?;
                hsps.hsp(&XmlHsp::unaligned(26, 14.62344, 2.880961, -1))
            })
        })
        .unwrap();
    let text = String::from_utf8(out).unwrap();
    let hit = concat!(
        "<Hit>\n",
        "  <Hit_num>1</Hit_num>\n",
        "  <Hit_id>ref|NP_100001.1|</Hit_id>\n",
        "  <Hit_def>Unknown</Hit_def>\n",
        "  <Hit_accession>Unknown</Hit_accession>\n",
        "  <Hit_len>0</Hit_len>\n",
        "  <Hit_hsps>\n",
        "    <Hsp>\n",
        "      <Hsp_num>1</Hsp_num>\n",
        "      <Hsp_bit-score>139.813</Hsp_bit-score>\n",
        "      <Hsp_score>351</Hsp_score>\n",
        "      <Hsp_evalue>1.14563e-47</Hsp_evalue>\n",
        "      <Hsp_query-from>0</Hsp_query-from>\n",
        "      <Hsp_query-to>0</Hsp_query-to>\n",
        "      <Hsp_hit-from>0</Hsp_hit-from>\n",
        "      <Hsp_hit-to>0</Hsp_hit-to>\n",
        "      <Hsp_identity>66</Hsp_identity>\n",
        "      <Hsp_qseq></Hsp_qseq>\n",
        "      <Hsp_hseq></Hsp_hseq>\n",
        "    </Hsp>\n",
        "    <Hsp>\n",
        "      <Hsp_num>2</Hsp_num>\n",
        "      <Hsp_bit-score>14.6234</Hsp_bit-score>\n",
        "      <Hsp_score>26</Hsp_score>\n",
        "      <Hsp_evalue>2.88096</Hsp_evalue>\n",
        "      <Hsp_query-from>0</Hsp_query-from>\n",
        "      <Hsp_query-to>0</Hsp_query-to>\n",
        "      <Hsp_hit-from>0</Hsp_hit-from>\n",
        "      <Hsp_hit-to>0</Hsp_hit-to>\n",
        "      <Hsp_identity>-1</Hsp_identity>\n",
        "      <Hsp_qseq></Hsp_qseq>\n",
        "      <Hsp_hseq></Hsp_hseq>\n",
        "    </Hsp>\n",
        "  </Hit_hsps>\n",
        "</Hit>\n",
        "</Iteration_hits>",
    );
    assert!(text.contains(hit), "{text}");
    // The run goes on: no cut of the query's loads, nothing posted.
    let posts = writer.take_posts();
    assert!(!posts.cut && posts.ends_run.is_none() && posts.posts.is_empty());
    // A `Hit` that loads its subject reads its HSPs' alignments.
    let mut out = Vec::new();
    let loaded = XmlHit {
        lookup_throws: None,
        ..thrown
    };
    writer
        .query(&run, &query, &mut out, |sink| {
            sink.hit(&loaded, |hsps| {
                assert!(hsps.reads_alignment());
                Ok(())
            })
        })
        .unwrap();
    // XML2 and JSON raise the exception at the `Hit`.
    for kind in [
        crate::report::blast_xml2::Xml2Kind::Single16,
        crate::report::blast_xml2::Xml2Kind::Single15,
    ] {
        let mut document = crate::report::blast_xml2::Xml2Writer::new(
            kind,
            crate::report::blast_xml2::Xml2Run {
                target: crate::report::blast_xml2::Xml2Target::Database,
                subject_names: None,
                gapped: true,
                cbs: 2,
                query_gencode: None,
                db_gencode: None,
                program_type: None,
                query_dropped: false,
            },
        );
        let error = document
            .query(&run, &query, &mut Vec::new(), |sink| {
                sink.hit(&thrown, |_| panic!("no HSP of a `Hit` that raises"))
            })
            .unwrap_err();
        let caught = error.downcast::<CaughtException>().unwrap();
        assert_eq!(caught.0, exception);
    }
}

// Line PS, PS15b-4: PHI-BLAST's pattern as read is `Parameters_pattern`, after `Parameters_filter`
// (blastxml.asn's order; blastxml_format.cpp:617-618).
#[test]
fn a_phi_pattern_follows_the_filter_in_the_parameters() {
    let matrix = FormattingMatrix::new("BLOSUM62").unwrap();
    let run = XmlRun {
        pattern: Some(b"[AG]-x(4)-G-K-[ST]"),
        ..run(&matrix)
    };
    let query = XmlQuery {
        id: b"Query_1",
        title: b"noh first query",
        len: 32,
        stats: STATS,
        message: Some(b"No hits found"),
        searched: true,
        masks: &[],
        joins_held: false,
    };
    let mut out = Vec::new();
    let mut writer = BlastXmlWriter::new(0);
    writer.query(&run, &query, &mut out, |_| Ok(())).unwrap();
    let filter = "      <Parameters_filter>F</Parameters_filter>\n";
    let header = HEADER.replacen(
        filter,
        &format!("{filter}      <Parameters_pattern>[AG]-x(4)-G-K-[ST]</Parameters_pattern>\n"),
        1,
    );
    assert_ne!(header, HEADER);
    assert!(
        String::from_utf8(out).unwrap().starts_with(&header),
        "the header with the pattern"
    );
}

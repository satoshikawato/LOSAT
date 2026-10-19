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
        version: "BLASTP 2.17.0+",
        reference: REFERENCE_GAPPED_BLAST,
        db_name: b"",
        matrix: Some("BLOSUM62"),
        expect: 10.0,
        gap_open: 11,
        gap_extend: 1,
        reward_penalty: None,
        filter: b"F",
        db_num: 0,
        db_len: 0,
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
        def: b"noh first query",
        len: 32,
        stats: STATS,
        message: Some(b"No hits found"),
    };
    writer.query(&run, &first, out, |_| Ok(())).unwrap();
    let second = XmlQuery {
        id: b"Query_2",
        def: b"qp2 unrelated protein (no hits expected)",
        len: 45,
        stats: STATS,
        message: None,
    };
    writer
        .query(&run, &second, out, |sink| {
            sink.hit(
                &XmlHit {
                    id: b"sp1",
                    def: b"sp1 subject protein alpha",
                    accession: b"Subject_1",
                    len: 137,
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
                        masks: &[],
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
                        masks: &[],
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
    writer.finish(&mut out);
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
        def: b"qp2 unrelated protein (no hits expected)",
        len: 45,
        stats: STATS,
        message: None,
    };
    writer.query(&run, &query, &mut out, |_| Ok(())).unwrap();
    let text = String::from_utf8(out.clone()).unwrap();
    assert!(text.starts_with("\n<Iteration>\n  <Iteration_iter-num>2</Iteration_iter-num>"));
    assert!(!text.contains("<BlastOutput>"));
    writer.finish(&mut out);
    assert!(out.ends_with(TRAILER.as_bytes()));
    // an epilog call for a run with no query writes nothing
    let mut none = Vec::new();
    BlastXmlWriter::new(0).finish(&mut none);
    assert!(none.is_empty());
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
        def: b"No definition line",
        len: 0,
        stats: XmlStats::NONE,
        message: Some(b"x"),
    };
    writer.query(&run, &query, &mut out, |_| Ok(())).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("      <Parameters_expect>0.001</Parameters_expect>\n      <Parameters_sc-match>1</Parameters_sc-match>\n      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>\n      <Parameters_gap-open>0</Parameters_gap-open>"));
    assert!(!text.contains("Parameters_matrix"));
    assert!(text.contains("<Statistics_db-num>3</Statistics_db-num>\n      <Statistics_db-len>10000000000</Statistics_db-len>\n      <Statistics_hsp-len>0</Statistics_hsp-len>\n      <Statistics_eff-space>0</Statistics_eff-space>\n      <Statistics_kappa>-1</Statistics_kappa>\n      <Statistics_lambda>-1</Statistics_lambda>\n      <Statistics_entropy>-1</Statistics_entropy>"));
}

fn run_blastn() -> XmlRun<'static> {
    XmlRun {
        program: "blastn",
        version: "BLASTN 2.17.0+",
        reference: REFERENCE_MEGABLAST,
        db_name: b"db/snucl",
        matrix: None,
        expect: 0.001,
        gap_open: 0,
        gap_extend: 0,
        reward_penalty: Some((1, -2)),
        filter: b"L;m;",
        db_num: 0,
        db_len: 0,
        positives: None,
        mask_char: b'N',
    }
}

#[test]
fn strings_are_escaped_and_non_ascii_is_an_error() {
    let mut out = Vec::new();
    let mut enc = XmlEncoder::object(&mut out);
    enc.begin_sequence("Hit");
    enc.member("def", 2);
    enc.visible_str(b"q2 & <x> \"quoted\" 'apos' a&b\t");
    enc.end_member();
    enc.end_sequence();
    assert!(enc.error().is_none());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "\n<Hit>\n  <Hit_def>q2 &amp; &lt;x&gt; &quot;quoted&quot; &apos;apos&apos; a&amp;b&#x9;</Hit_def>\n</Hit>"
    );
    let mut out = Vec::new();
    let mut enc = XmlEncoder::object(&mut out);
    enc.visible_str("caf\u{e9}".as_bytes());
    assert!(enc.error().is_some());
    let mut enc = XmlEncoder::object(&mut out);
    enc.visible_str(b"a\x01");
    assert!(enc.error().is_some());
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

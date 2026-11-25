use super::*;

#[test]
fn thresholds_are_written_as_a_cpp_stream_writes_a_double() {
    for (value, text) in [
        (11.0, "11"),
        (999999.0, "999999"),
        (1000000.0, "1e+06"),
        (1234567.0, "1.23457e+06"),
        (2147483647.0, "2.14748e+09"),
        (12.5, "12.5"),
        (0.0001, "0.0001"),
        (0.00001, "1e-05"),
    ] {
        assert_eq!(cpp_default_double(value), text, "{value}");
    }
}

/// The rows of a protein alignment `query` x `subject` (`-` gaps) from (`q_start`, `s_start`),
/// with the edit script of their columns (`gap_info`: a gap in the subject is an insertion,
/// one in the query a deletion), from which the report makes the HSP's Dense-seg.
fn rows_hit(query: &str, subject: &str, q_start: usize, s_start: usize) -> Hit {
    use crate::common::GapEditOp;
    let residues = |row: &str| row.bytes().filter(|&b| b != b'-').count();
    let mut ops: Vec<GapEditOp> = Vec::new();
    for (q, s) in query.bytes().zip(subject.bytes()) {
        let op = match (q == b'-', s == b'-') {
            (false, true) => GapEditOp::Ins(1),
            (true, false) => GapEditOp::Del(1),
            _ => GapEditOp::Sub(1),
        };
        match (ops.last_mut(), op) {
            (Some(GapEditOp::Sub(num)), GapEditOp::Sub(_))
            | (Some(GapEditOp::Ins(num)), GapEditOp::Ins(_))
            | (Some(GapEditOp::Del(num)), GapEditOp::Del(_)) => *num += 1,
            _ => ops.push(op),
        }
    }
    Hit {
        identity: 0.0,
        length: query.len(),
        mismatch: 0,
        gapopen: 1,
        q_start,
        q_end: q_start + residues(query) - 1,
        s_start,
        s_end: s_start + residues(subject) - 1,
        e_value: 1e-10,
        bit_score: 50.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
        raw_score: 100,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: Some(ops),
        num_positives: 0,
    }
}

fn rows_text(query: &str, subject: &str, q_start: usize, s_start: usize, line: usize) -> String {
    let hit = rows_hit(query, subject, q_start, s_start);
    let config = PairwiseConfig {
        line_length: line,
        program: "rpsblast".to_string(),
        ..PairwiseConfig::default()
    };
    let mut out = Vec::new();
    write_alignment_with_sequences(
        &mut out,
        &hit,
        query,
        subject,
        [&[], &[]],
        Default::default(),
        None,
        &config,
    )
    .unwrap();
    String::from_utf8(out).unwrap()
}

// showalign.cpp:1209-1216 and 1497-1502 (cited at `end_gap_span` and
// `write_alignment_with_sequences`): a row's end gaps are spaces, a chunk without the row's
// letters leaves the row out (the Query row with its middle line), and the next chunk's row
// shows its first residue's number.
#[test]
fn end_gaps_are_spaces_and_rows_without_letters_are_left_out() {
    let middle = |letters: &str| format!("{}{letters}\n", " ".repeat(10));
    // The subject starts with six end gaps: the first chunk has no Sbjct line, the second
    // starts it at its residue 1 after two spaces.
    let text = rows_text("MKVLWAALLV", "------AALL", 1, 1, 4);
    let expected = [
        "Query  1  MKVL  4\n".to_string(),
        middle("    "),
        "\n".to_string(),
        "Query  5  WAAL  8\n".to_string(),
        middle("  A "),
        "Sbjct  1    AA  2\n\n".to_string(),
        "Query  9  LV  10\n".to_string(),
        middle("L+"),
        "Sbjct  3  LL  4\n\n".to_string(),
    ]
    .concat();
    assert_eq!(text, expected);
    // The query starts with three end gaps.
    let text = rows_text("---KVLW", "MQLKVLW", 7, 1, 60);
    let expected = [
        "Query  7     KVLW  10\n".to_string(),
        middle("   KVLW"),
        "Sbjct  1  MQLKVLW  7\n\n".to_string(),
    ]
    .concat();
    assert_eq!(text, expected);
    // A query chunk of end gaps only: no Query line and no middle line.
    let text = rows_text("----KVLW", "MQLWKVLW", 1, 1, 4);
    let expected = [
        "Sbjct  1  MQLW  4\n\n".to_string(),
        "Query  1  KVLW  4\n".to_string(),
        middle("KVLW"),
        "Sbjct  5  KVLW  8\n\n".to_string(),
    ]
    .concat();
    assert_eq!(text, expected);
}

fn tblastx_hit(sum_n: i32) -> PairwiseHit {
    let query: String = "ACDEFGHIKL".repeat(6) + "M";
    let subject: String = "ACDEFGHIKL".repeat(6) + "N";
    PairwiseHit {
        hit: Hit {
            identity: 0.0,
            length: 61,
            mismatch: 1,
            gapopen: 0,
            q_start: 200,
            q_end: 18,
            s_start: 3,
            s_end: 185,
            e_value: 2e-30,
            bit_score: 120.5,
            num_ident: 60,
            query_frame: -2,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: 255,
            sort_query_offset: 0,
            sort_query_end: 61,
            sort_subject_offset: 0,
            sort_subject_end: 61,
            has_sort_offsets: true,
            gap_info: None,
            num_positives: 60,
        },
        query_seq: Some(query),
        subject_seq: Some(subject),
        query_frame: Some(-2),
        subject_frame: Some(3),
        positives: Some(60),
        gaps: Some(0),
        subject_length: Some(400),
        subject_title: None,
        comp_adjust_method: None,
        sum_n: Some(sum_n),
        subject_stored_gaps: Vec::new(),
        query_stored_gaps: Vec::new(),
        query_case: Default::default(),
        seqalign_num_ident: 60,
    }
}

// Line DBID item ID4t: TBLASTX's outfmt 0 makes its rows a window of queries at a time and
// writes each window's report (`algorithm::tblastx::blast_engine::run_impl`), as BLASTN's
// (ID4f, above): the reports of consecutive query ranges, each given only its own hits, are
// the report of all queries, on one thread and on a pool.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn tblastx_reports_of_consecutive_query_ranges_are_the_whole_report() {
    let config = PairwiseConfig {
        program: "tblastx".to_string(),
        ..PairwiseConfig::default()
    };
    let karlin = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.0,
        beta: 0.0,
    };
    let queries: Vec<TblastxPairwiseQuery> = (0..4)
        .map(|q_idx| TblastxPairwiseQuery {
            query_name: format!("q{q_idx} query {q_idx}").into_bytes(),
            query_length: 300,
            karlin: Some(karlin),
            effective_search_space: 10_000,
        })
        .collect();
    let report = TblastxPairwiseReport {
        version: "2.17.0+".to_string(),
        database_name: b"subjects".to_vec(),
        database_num_sequences: 3,
        database_total_letters: 900,
        matrix_name: "BLOSUM62".to_string(),
        word_threshold: 13.0,
        window_size: 40,
        num_descriptions: 500,
        num_alignments: 250,
        // Query 1 is not searched (a batch without a valid context).
        unsearched: vec![false, true, false, false],
        epilog: true,
        prolog: true,
        parsed_subjects: Vec::new(),
        subject_sources: Vec::new(),
        sum_stats: true,
    };
    // Query 0: subjects 1 (2 HSPs) and 2; query 1: no hits; queries 2 and 3: subject 1.
    let hits: Vec<PairwiseHit> = [(0, 1, 4), (0, 1, 4), (0, 2, 1), (2, 1, 1), (3, 1, 1)]
        .into_iter()
        .map(|(q_idx, s_idx, sum_n)| {
            let mut hit = tblastx_hit(sum_n);
            hit.hit.q_idx = q_idx;
            hit.hit.s_idx = s_idx;
            hit
        })
        .collect();
    let titles: Vec<Arc<[u8]>> = (0..3)
        .map(|s_idx| Arc::from(format!("s{s_idx} subject {s_idx}").as_bytes()))
        .collect();
    // The warnings written before queries 0 and 2 (`QueryWarnings`, by query index).
    let before: Vec<crate::report::diag_log::Posts> = (0..queries.len())
        .map(|q_idx| match q_idx {
            0 | 2 => crate::report::diag_log::Posts::of_lines(
                format!(
                    "Warning: [tblastx] Query_{} q{q_idx}: a warning\n",
                    q_idx + 1
                )
                .into_bytes(),
            ),
            _ => crate::report::diag_log::Posts::new(),
        })
        .collect();
    crate::utils::threading::with_search_pool(4, "test", |pool| {
        for x_pool in [None, Some(pool)] {
            // The report and the warnings of the queries `range`.
            let write = |range: std::ops::Range<usize>| {
                let range_hits: Vec<PairwiseHit> = hits
                    .iter()
                    .filter(|hit| range.contains(&(hit.hit.q_idx as usize)))
                    .cloned()
                    .collect();
                let mut written = Vec::new();
                let mut diagnostics = Vec::new();
                write_tblastx_pairwise_report(
                    &range_hits,
                    &mut written,
                    &config,
                    &queries,
                    range,
                    &titles,
                    &report,
                    None,
                    None,
                    Some(&mut crate::report::query_warnings::QueryWarnings::of_posts(
                        &before,
                        None,
                        &mut diagnostics,
                    )),
                    x_pool,
                    None,
                    None,
                    &crate::report::db_subject::HitLoads::default(),
                )
                .unwrap();
                (written, diagnostics)
            };
            let (whole, whole_warnings) = write(0..queries.len());
            let text = String::from_utf8(whole.clone()).unwrap();
            assert!(text.starts_with("TBLASTX 2.17.0+"), "{text}");
            assert_eq!(text.matches("Query= ").count(), queries.len(), "{text}");
            assert!(text.contains("No hits found"), "{text}");
            assert!(text.contains("Expect(4)"), "{text}");
            assert!(text.ends_with("Window for multiple hits: 40\n"), "{text}");
            assert_eq!(text.matches("Window for multiple hits").count(), 1, "{text}");
            assert_eq!(
                String::from_utf8(whole_warnings.clone()).unwrap(),
                "Warning: [tblastx] Query_1 q0: a warning\nWarning: [tblastx] Query_3 q2: a warning\n"
            );
            for cuts in [vec![1], vec![2], vec![3], vec![1, 2, 3], vec![1, 3]] {
                let mut joined = Vec::new();
                let mut joined_warnings = Vec::new();
                let mut start = 0;
                for end in cuts.into_iter().chain([queries.len()]) {
                    let (written, warnings) = write(start..end);
                    joined.extend(written);
                    joined_warnings.extend(warnings);
                    start = end;
                }
                assert_eq!(String::from_utf8(joined).unwrap(), text);
                assert_eq!(joined_warnings, whole_warnings);
            }
        }
        Ok(())
    })
    .unwrap();
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3593-3606 and 310-332: Expect(n)
// only for a linked set of two or more, both frames, no Strand line.
#[test]
fn tblastx_hsp_info_prints_the_linked_set_and_both_frames() {
    let mut out = Vec::new();
    write_tblastx_hsp_info(&mut out, &tblastx_hit(4)).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        " Score = 120 bits (255),  Expect(4) = 2e-30\n Identities = 60/61 (98%), Positives = 60/61 (98%), Gaps = 0/61 (0%)\n Frame = -2/+3\n\n"
    );
    let mut out = Vec::new();
    write_tblastx_hsp_info(&mut out, &tblastx_hit(1)).unwrap();
    assert!(String::from_utf8(out)
        .unwrap()
        .starts_with(" Score = 120 bits (255),  Expect = 2e-30\n"));
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:230,305,358 and c++/src/objtools/align_format/showalign.cpp:1598-1626: 60
// residues span 180 nucleotides, each row in the direction of its frame, and the width
// is that of the largest 0-based coordinate (199: 3 digits).
#[test]
fn tblastx_rows_step_three_nucleotides_in_each_frame_direction() {
    let mut out = Vec::new();
    // The query's record has 300 nucleotides (frame -2 at 200).
    write_tblastx_alignment(
        &mut out,
        &tblastx_hit(1),
        None,
        300,
        &PairwiseConfig::default(),
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let row = "ACDEFGHIKL".repeat(6);
    assert_eq!(lines[0], format!("Query  200  {row}  21"));
    assert_eq!(lines[1], format!("            {row}"));
    assert_eq!(lines[2], format!("Sbjct  3    {row}  182"));
    assert_eq!(lines[3], "");
    assert_eq!(lines[4], "Query  20   M  18");
    assert_eq!(lines[5], "             ");
    assert_eq!(lines[6], "Sbjct  183  N  185");
}

// D5: NCBI BLAST+ 2.17.0's description rows of `-subject` records without parsed ids
// (`blastp -subject ps_spec.faa -outfmt 0 -html`, probes `bp.s.f0.h` and, with
// `T/t4/rc/localid/.ncbirc`, `bp.s.f0.rc_localid.h`): the title encoded after the cut on raw
// bytes, the padding of raw bytes, the score link of `Subject_N`, and the local-id link.
#[test]
fn html_description_rows_are_ncbis() {
    let mut first = tblastx_hit(1);
    first.hit.bit_score = 313.0;
    first.hit.e_value = 2e-109;
    let mut third = tblastx_hit(1);
    third.hit.s_idx = 2;
    third.hit.bit_score = 247.0;
    third.hit.e_value = 2e-87;
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&first]), (2, vec![&third])].into_iter().collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> = [
        (
            0,
            b"s1 a<b>c & d \"e\" &amp; &lt;f&gt; &#62; &auml; short".to_vec(),
        ),
        (
            2,
            [
                &b"s3 "[..],
                &[b'x'; 70],
                b" & tail after a word of seventy letters",
            ]
            .concat(),
        ),
    ]
    .into_iter()
    .collect();
    let rows = |file: Option<&[u8]>| {
        let registry = crate::report::formatter_registry::InitConfig::of(
            "blastp",
            file.map(|bytes| (".ncbirc".into(), bytes.to_vec())),
            false,
        )
        .unwrap();
        let links = Arc::new(crate::report::html_ids::HtmlLinks {
            registry,
            db_name: Vec::new(),
            nucleotide: false,
            ncbirc: Default::default(),
        });
        let report = crate::report::html_ids::HtmlReport::new(&links, []);
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0, 2],
            true,
            &hits,
            &titles,
            &[],
            &[],
            false,
            true,
            None,
            &mut crate::report::query_warnings::no_places,
            Some((
                &report,
                crate::report::db_subject::AfterThrow::SkipsAlignment { html: true },
            )),
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .skip(3)
            .take(2)
            .map(String::from)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        rows(None),
        [
            "s1 a&lt;b&gt;c &amp; d &quot;e&quot; &amp; &lt;f&gt; &gt; ä short                                     <a href=#Subject_1>313</a>     2e-109",
            "s3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx...  <a href=#Subject_3>247</a>     2e-87 ",
        ]
    );
    let file: &[u8] = b"[LOCAL_ID]\nTOOL_URL=http://local.example/<@seq_id@>?db=<@db_name@>&taxid=<@taxid@>&p=<@protocol@>&t=<@title@>\n";
    assert_eq!(
        rows(Some(file)),
        [
            "http://local.example/Subject_1?db=&taxid=0&p=https:&t=Subject_1</a>s1 a&lt;b&gt;c &amp; d &quot;e&quot; &amp; &lt;f&gt; &gt; ä short                                     <a href=#Subject_1>313</a>     2e-109",
            "http://local.example/Subject_3?db=&taxid=0&p=https:&t=Subject_3</a>s3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx...  <a href=#Subject_3>247</a>     2e-87 ",
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:829-836,961-965: the N column of
// an ungapped sum-statistics search, from the first HSP of each subject.
#[test]
fn description_table_shows_the_linked_set_column() {
    let first = tblastx_hit(12);
    let mut second = tblastx_hit(1);
    second.hit.s_idx = 1;
    second.hit.bit_score = 30.2;
    second.hit.e_value = 0.5;
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&first]), (1, vec![&second])]
            .into_iter()
            .collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> =
        [(0, b"s0".to_vec()), (1, b"s1".to_vec())]
            .into_iter()
            .collect();
    let mut out = Vec::new();
    write_blastn_description_table(
        &mut out,
        &[0, 1],
        true,
        &hits,
        &titles,
        &[],
        &[],
        true,
        false,
        None,
        &mut crate::report::query_warnings::no_places,
        None,
        &crate::report::db_subject::HitLoads::default(),
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert!(lines[0].ends_with("Score     E"), "{text}");
    assert!(lines[1].ends_with("(Bits)  Value  N"), "{text}");
    assert!(lines[3].ends_with("  120     2e-30  12"), "{text}");
    assert!(lines[4].ends_with("  30.2    0.50   1 "), "{text}");
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:771-840,930-975 (`-sorthits`): the
// Total score, Query cover and Max ident columns, as NCBI BLAST+ 2.17.0 `blastp -sorthits 0`
// writes them (probe t1 `blastp.s.h0kx`): the total in the HSPs' order, the cover of the
// union of the query ranges, the identity of the best HSP (its Seq-align's `num_ident`,
// `PairwiseHit::seqalign_num_ident`), and the trailing spaces of the nine letters of
// "Max ident".
#[test]
fn sorthits_adds_the_total_cover_and_identity_columns() {
    let mut a = tblastx_hit(1);
    (
        a.hit.bit_score,
        a.hit.e_value,
        a.seqalign_num_ident,
        a.hit.length,
    ) = (28.1, 0.18, 27, 99);
    (a.hit.q_start, a.hit.q_end) = (1, 99);
    let mut b1 = a.clone();
    b1.hit.s_idx = 1;
    (
        b1.hit.bit_score,
        b1.hit.e_value,
        b1.seqalign_num_ident,
        b1.hit.length,
    ) = (24.6, 2.3, 16, 25);
    (b1.hit.q_start, b1.hit.q_end) = (10, 34);
    let mut b2 = b1.clone();
    (
        b2.hit.bit_score,
        b2.hit.e_value,
        b2.seqalign_num_ident,
        b2.hit.length,
    ) = (24.6, 2.4, 10, 20);
    (b2.hit.q_start, b2.hit.q_end) = (49, 30);
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&a]), (1, vec![&b1, &b2])].into_iter().collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> =
        [(0, b"s0".to_vec()), (1, b"s1".to_vec())]
            .into_iter()
            .collect();
    let columns = Some(SortColumns {
        query_length: 897,
        tblastx: false,
    });
    let mut out = Vec::new();
    write_blastn_description_table(
        &mut out,
        &[0, 1],
        true,
        &hits,
        &titles,
        &[],
        &[],
        false,
        true,
        columns,
        &mut crate::report::query_warnings::no_places,
        None,
        &crate::report::db_subject::HitLoads::default(),
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        lines[0],
        format!("{}Score   Total  Query    E     Max", " ".repeat(70))
    );
    assert_eq!(
        lines[1],
        "Sequences producing significant alignments:                          (Bits)   Score  cover   Value  Ident"
    );
    assert_eq!(lines[2], "");
    let label = |name: &str| format!("{name:<68}");
    // 99 of 897 letters (11 %), 27 of 99 identities (27 %).
    assert_eq!(
        lines[3],
        format!("{}  28.1     28.1   11%    0.18   27%       ", label("s0"))
    );
    // 24.6 + 24.6; positions 10-49 (40 of 897: 4 %); the first of the equal best HSPs.
    assert_eq!(
        lines[4],
        format!("{}  24.6     49.2   4%     2.3    64%       ", label("s1"))
    );
    // TBLASTX: the identity over three nucleotides per amino acid (60 of 183: 33 %), after
    // the N column ("NIdent").
    let t = tblastx_hit(2);
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&t])].into_iter().collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> =
        [(0, b"t".to_vec())].into_iter().collect();
    let columns = Some(SortColumns {
        query_length: 300,
        tblastx: true,
    });
    let mut out = Vec::new();
    write_blastn_description_table(
        &mut out,
        &[0],
        true,
        &hits,
        &titles,
        &[],
        &[],
        true,
        false,
        columns,
        &mut crate::report::query_warnings::no_places,
        None,
        &crate::report::db_subject::HitLoads::default(),
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert!(
        lines[1].ends_with("(Bits)   Score  cover   Value  NIdent"),
        "{text}"
    );
    // Positions 18-200: 183 of 300 (61 %).
    assert!(
        lines[3].ends_with("  120      120    61%    2e-30  2  33%       "),
        "{text}"
    );
    // The identity is the Seq-align's, not that of the rows: NCBI BLAST+ 2.17.0 `tblastn
    // -query` (the first record of t1 `tn_q.faa`) `-subject tn_s.fna -sorthits 0` shows
    // "94%" for an HSP of "Identities = 284/286" whose engine counted the query's SEG-masked
    // letters as mismatches, and "99%" with `-seg no`; here 270 engine identities of 286.
    let mut masked = tblastx_hit(1);
    (
        masked.hit.num_ident,
        masked.seqalign_num_ident,
        masked.hit.length,
    ) = (284, 270, 286);
    (masked.hit.q_start, masked.hit.q_end) = (1, 286);
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&masked])].into_iter().collect();
    let columns = Some(SortColumns {
        query_length: 286,
        tblastx: false,
    });
    let mut out = Vec::new();
    write_blastn_description_table(
        &mut out,
        &[0],
        true,
        &hits,
        &titles,
        &[],
        &[],
        false,
        true,
        columns,
        &mut crate::report::query_warnings::no_places,
        None,
        &crate::report::db_subject::HitLoads::default(),
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(
        text.lines()
            .nth(3)
            .unwrap()
            .ends_with("  100%   2e-30  94%       "),
        "{text}"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1133-1138,940-944 and
// align_format_util.cpp:932-938: a last row whose total is wider than every earlier total
// pads the total column with a negative `size_t`, where NCBI writes spaces without end; the
// table without `-sorthits` only widens the score column.
#[test]
fn sorthits_refuses_a_last_total_wider_than_its_column() {
    let mut a = tblastx_hit(1);
    a.hit.bit_score = 60000.0;
    let mut b = a.clone();
    b.hit.s_idx = 1;
    let mut c = a.clone();
    c.hit.s_idx = 1;
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&a]), (1, vec![&b, &c])].into_iter().collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> =
        [(0, b"s0".to_vec()), (1, b"s1".to_vec())]
            .into_iter()
            .collect();
    let columns = Some(SortColumns {
        query_length: 300,
        tblastx: false,
    });
    let table = |order: &[u32], last_row_counted: bool, columns: Option<SortColumns>| {
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            order,
            last_row_counted,
            &hits,
            &titles,
            &[],
            &[],
            false,
            true,
            columns,
            &mut crate::report::query_warnings::no_places,
            None,
            &crate::report::db_subject::HitLoads::default(),
        )
        .map(|()| out)
    };
    let error = table(&[0, 1], true, columns).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::Unsupported);
    assert_eq!(error.to_string(), TotalScoreRefusal.to_string());
    assert!(is_total_score_refusal(&error.into()));
    assert!(!is_total_score_refusal(
        &io::Error::new(io::ErrorKind::Unsupported, "other").into()
    ));
    assert!(!is_total_score_refusal(&anyhow::anyhow!("other")));
    assert!(table(&[0, 1], true, None).is_ok());
    // Not the last subject of the list, or an earlier total as wide: written.
    assert!(table(&[0, 1], false, columns).is_ok());
    assert!(table(&[1, 1], true, columns).is_ok());
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4259,4299-4302: BLASTX's row
// shows the subject's highest bit score and that HSP's E-value (the first of equal
// maxima), not its first HSP by E-value, and the N column the first HSP's linked-set
// size. NCBI BLAST+ 2.17.0 `blastx -query tests/fasta/outfmt0/bx_best_hsp_query.fna
// -subject tests/fasta/outfmt0/bx_best_hsp_subject.faa` (fixture `best_hsp.blastx`):
// HSPs of 119 bits (Expect(2) = 8e-38), 45.4 and 131 (5e-35); the row is `131     5e-35`.
#[test]
fn blastx_description_row_is_the_highest_bit_score() {
    let mut first = tblastx_hit(2);
    first.hit.bit_score = 119.2;
    first.hit.e_value = 8e-38;
    let mut best = tblastx_hit(1);
    best.hit.bit_score = 131.3;
    best.hit.e_value = 5e-35;
    let mut tie = best.clone();
    tie.hit.e_value = 0.5;
    for hit in [&mut first, &mut best, &mut tie] {
        hit.subject_title = Some("pep YOL075C".to_string());
    }
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&first, &best, &tie])].into_iter().collect();
    let ids = [DescriptionLabel::title_only(
        b"YOL075C pep YOL075C".to_vec(),
    )];
    for (show_sum_n, end) in [(false, ""), (true, "  2")] {
        let mut out = Vec::new();
        write_subject_summary_table_with_sum_n(
            &mut out,
            &[0],
            true,
            &hits,
            &ids,
            show_sum_n,
            None,
            &mut crate::report::query_warnings::no_places,
            None,
            None,
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[3],
            format!("{:<68}  131     5e-35{end}", "YOL075C pep YOL075C"),
            "{text}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452 (cited at
// `QueryWarnings`): a query's warning goes to `cerr`, whose tie flushes the report before it.
// Line RG T3d: BLASTX's report goes through the shared per-query writer
// (`x_pairwise_par::write`): a query's warning comes after the reports written before it,
// also from a writer that holds its bytes until it is flushed (LOSAT's buffer of the command
// line, `cli::NcbiOstream::buffered`); a query whose description table reads an HSP with empty
// rows stops after its `Length=` (`EmptyRowsStop::RowRange`): its error comes after the
// report's last bytes, and nothing is written after it, no text of a later read either
// (`QueryWarnings::end_run`); the parallel writer (`x_pairwise_par::write` given a pool, which
// `LOSAT_X_PAIRWISEPAR` does not turn off here) writes the same bytes.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn blastx_query_warnings_follow_the_reports_before_them() {
    use std::sync::{Arc as SyncArc, Mutex};
    #[derive(Clone, Default)]
    struct Log(SyncArc<Mutex<Vec<u8>>>);
    impl Write for Log {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    /// Holds what is written until it is flushed into the log.
    struct Held(Vec<u8>, Log);
    impl Write for Held {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.0.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            self.1.write_all(&std::mem::take(&mut self.0))
        }
    }
    impl crate::report::ostream::ReportWrite for Held {}
    let first = tblastx_hit(1);
    let mut empty = tblastx_hit(1);
    empty.hit.q_idx = 2;
    empty.hit.sort_query_end = empty.hit.sort_query_offset;
    let query = |name: &[u8]| BlastpPairwiseQuery {
        query_name: name.to_vec(),
        query_length: 300,
        ungapped_karlin: Default::default(),
        effective_search_space: 1,
        valid: true,
        batch_skipped: false,
    };
    let queries = [query(b"q0"), query(b"q1"), query(b"q2")];
    let report = BlastpPairwiseReport {
        version: String::new(),
        database_name: Vec::new(),
        database_num_sequences: 1,
        database_total_letters: 400,
        matrix_name: Some("BLOSUM62".into()),
        gap_open: 11,
        gap_extend: 1,
        word_threshold: 12.0,
        window_size: 40,
        gapped_karlin: Default::default(),
        gumbel: None,
        gapped: true,
        show_linked_set_size: false,
        num_descriptions: usize::MAX,
        num_alignments: usize::MAX,
    };
    let options = BlastxPairwiseOptions {
        word_threshold: 12.0,
        gapped: true,
        sum_stats: false,
        num_descriptions: 500,
        num_alignments: 0,
    };
    let posts: Vec<crate::report::diag_log::Posts> = [&b""[..], b"W1\n", b"W2\n"]
        .map(|lines| crate::report::diag_log::Posts::of_lines(lines.to_vec()))
        .into();
    // Query 2's table would read subject 0 (its text `T0`): the stop comes before it.
    let reads = vec![
        Vec::new(),
        Vec::new(),
        vec![crate::report::query_warnings::SubjectRead::of_text(
            0, b"T0\n",
        )],
    ];
    let config = PairwiseConfig {
        program: "blastx".into(),
        ..PairwiseConfig::default()
    };
    let run = |pool: Option<&crate::utils::threading::SearchPool<'_>>| -> (String, String) {
        let log = Log::default();
        let mut sink = log.clone();
        let mut held = Held(Vec::new(), log.clone());
        let mut warnings =
            crate::report::query_warnings::QueryWarnings::of_posts(&posts, None, &mut sink)
                .with_subject_reads(&reads);
        let error = write_blastx_pairwise_report(
            &[first.clone(), empty.clone()],
            &mut held,
            Some(&mut warnings),
            &config,
            &queries,
            &[DescriptionLabel::title_only(b"s0 title".to_vec())],
            &report,
            &options,
            None,
            &[],
            &[],
            None,
            None,
            0,
            None,
            pool,
            None,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap_err();
        // The merged stream when the error comes back: the report's last bytes are there
        // (written to the stream, or held by its writer until the flush that comes before the
        // error's post, as `cerr`'s tie flushes `cout`).
        let at_error = [log.0.lock().unwrap().as_slice(), held.0.as_slice()].concat();
        drop(warnings);
        held.flush().unwrap();
        let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
        assert_eq!(
            text.as_bytes(),
            at_error.as_slice(),
            "bytes written after the error (pool {})",
            pool.is_some()
        );
        (text, error.to_string())
    };
    let (serial, error) = run(None);
    assert_eq!(error, "RowRange(Empty)");
    let warning = serial.find("W1\n").expect("query 1's warning");
    assert!(
        serial[..warning].contains("Query= q0")
            && serial[..warning].ends_with("Effective search space used: 1\n"),
        "{serial}"
    );
    assert!(
        serial[warning..].starts_with("W1\n\n\nQuery= q1"),
        "{serial}"
    );
    assert!(
        serial.ends_with("W2\n\n\nQuery= q2\n\nLength=300\n"),
        "{serial}"
    );
    crate::utils::threading::x_with_search_pool_local(4, "test", |pool| {
        assert_eq!(run(Some(pool)), (serial.clone(), error.clone()));
        Ok(())
    })
    .unwrap();
}

/// The title that the FASTA reader gives a record with the defline `defline`.
fn reader_title(defline: &[u8]) -> Vec<u8> {
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
    let bytes = [&b">"[..], defline, b"\nACGT\n"].concat();
    let mut source =
        FastaInputSource::from_bytes(&bytes, ReaderConfig::query("BLASTN", false, false));
    read_all(&mut source, &mut |_| Ok(())).unwrap()[0]
        .title
        .clone()
}

// NCBI BLAST+ 2.17.0, `blastn -task blastn` outfmt 0 (inventory range RP,
// `scratch_RP/out/<case>/blastn_f0.out`): the `Query= ` line is the raw title wrapped
// by `NStr::Wrap` at 68 bytes (a line can end inside a UTF-8 sequence), and the
// heading is the title wrapped by `s_WrapOutputLine` (60 bytes, then the next white
// space). `w_<kind>_<k>`: `id `, k letters `a`, one character, `zz tail words`.
#[test]
fn query_lines_and_headings_wrap_bytes_as_ncbi() {
    let cjk = &b"\xe6\xbc\xa2"[..];
    let eacute = &b"\xc3\xa9"[..];
    let a = |count: usize| vec![b'a'; count];
    for (defline, query_block) in [
        (b"id desc".to_vec(), b"id desc\n".to_vec()),
        (b"id\tdesc".to_vec(), b"id\n".to_vec()),
        (b"".to_vec(), b"".to_vec()),
        (b"id\xe9x desc".to_vec(), b"id\xe9x desc\n".to_vec()),
        (
            b"id &amp; &lt;b&gt; \"q\" desc".to_vec(),
            b"id &amp; &lt;b&gt; \"q\" desc\n".to_vec(),
        ),
        (b"\x01id desc".to_vec(), b"\x01id desc\n".to_vec()),
        (b"id desc   ".to_vec(), b"id desc\n".to_vec()),
        (b"TPA: id desc".to_vec(), b"TPA: id desc\n".to_vec()),
        (
            [&b"id "[..], &a(59), cjk, b"zz tail words"].concat(),
            [&b"id "[..], &a(59), cjk, b"zz\ntail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(60), cjk, b"zz tail words"].concat(),
            [&b"id "[..], &a(60), cjk, b"zz\ntail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(63), cjk, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(63), cjk, b"zz\ntail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(64), cjk, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(64), cjk, b"z\nz tail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(65), cjk, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(65), cjk, b"\nzz tail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(66), cjk, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(66), b"\xe6\xbc\n\xa2zz tail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(67), cjk, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(67), b"\xe6\n\xbc\xa2zz tail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(66), eacute, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(66), eacute, b"\nzz tail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(67), eacute, b"zz tail words"].concat(),
            [&b"id\n"[..], &a(67), b"\xc3\n\xa9zz tail words\n"].concat(),
        ),
        (
            [&b"id "[..], &a(66), b"Zzz tail words"].concat(),
            [&b"id\n"[..], &a(66), b"Zz\nz tail words\n"].concat(),
        ),
    ] {
        let title = reader_title(&defline);
        let mut out = Vec::new();
        write_blastp_query_header(&mut out, &title, 330, false).unwrap();
        assert_eq!(
            out,
            [&b"Query= "[..], &query_block, b"\nLength=330\n"].concat(),
            "{defline:?}"
        );
    }
    for (middle, heading) in [
        (
            [&a(59), cjk].concat(),
            [&a(59), cjk, b"zz \ntail words"].concat(),
        ),
        (
            [&a(66), cjk].concat(),
            [&a(66), cjk, b"zz \ntail words"].concat(),
        ),
        (
            [&a(67), eacute].concat(),
            [&a(67), eacute, b"zz \ntail words"].concat(),
        ),
        (
            [&a(66)[..], b"Z"].concat(),
            [&a(66)[..], b"Zzz \ntail words"].concat(),
        ),
    ] {
        let title = [&b"id "[..], &middle, b"zz tail words"].concat();
        let mut out = Vec::new();
        write_subject_header(&mut out, &title, None, Some(450)).unwrap();
        assert_eq!(
            out,
            [&b"> id "[..], &heading, b"\nLength=450\n\n"].concat(),
            "{title:?}"
        );
    }
}

// NCBI BLAST+ 2.17.0 outfmt 0 for a subject whose title NCBI cannot decode (inventory
// range RP, `scratch_RP/out/u06_undef81/tblastn_f0.out`: `id \x81 x`, three HSPs):
// the description is `Unknown`, and each HSP is replaced by one line, without a heading.
#[test]
fn a_title_of_unknown_encoding_is_unknown_and_its_alignments_are_skipped() {
    let first = tblastx_hit(1);
    let mut second = tblastx_hit(1);
    second.hit.bit_score = 13.9;
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&first, &second])].into_iter().collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> =
        [(0, b"id \x81 x".to_vec())].into_iter().collect();
    let mut out = Vec::new();
    write_blastn_description_table(
        &mut out,
        &[0],
        true,
        &hits,
        &titles,
        &[],
        &[],
        false,
        false,
        None,
        &mut crate::report::query_warnings::no_places,
        None,
        &crate::report::db_subject::HitLoads::default(),
    )
    .unwrap();
    let row = out.split(|&byte| byte == b'\n').nth(3).unwrap().to_vec();
    assert!(row.starts_with(&[&b"Unknown"[..], &[b' '; 61], b"  120"].concat()));
    let mut out = Vec::new();
    write_subject_alignments_part(
        &mut out,
        None,
        0,
        &[],
        &titles[&0],
        false,
        &hits[&0],
        &[0, 1],
        0..hits[&0].len(),
        None,
        &mut crate::report::query_warnings::no_places,
        |_, _, _| panic!("no HSP is written"),
    )
    .unwrap();
    assert_eq!(
        out,
        b"Sequence with id Subject_1 no longer exists in database...alignment skipped\n\
              Sequence with id Subject_1 no longer exists in database...alignment skipped\n"
    );
}

// DBID pfoff-3 (DECISIONS 62 (13); NCBI BLAST+ 2.17.0, probe `R/pfoff/p3`, `tx.U2.b1.dflt.0`):
// a subject whose load throws where its batch's prefetch is off has the description line of
// its Seq-align id and `Unknown`, and in outfmt 0 each of its HSPs one line that names the
// id's accession without its version, with no heading (`DatabaseHeading::Thrown`).
#[test]
fn a_subject_whose_load_throws_is_unknown_and_its_alignments_are_skipped() {
    use crate::report::db_subject::{
        AfterThrow, FoundInstead, HitLoads, LookupFailure, ParsedSubject, RunThrow,
    };
    let subject = |accession: &[u8], throws: bool| {
        let subject = ParsedSubject::new(crate::seqid::DeflineSet {
            deflines: vec![crate::seqid::BlastDefline {
                title: Some(b"ntitle".to_vec()),
                ids: vec![crate::seqid::ParsedSeqId::Other(
                    crate::seqid::model::TextseqId {
                        accession: Some(accession.to_vec()),
                        version: Some(1),
                        ..Default::default()
                    },
                )],
                ..Default::default()
            }],
            utf8: false,
        })
        .unwrap();
        Some(Arc::new(if throws {
            subject.found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            ))))
        } else {
            subject
        }))
    };
    // 0 loads, 1 throws.
    let parsed = vec![subject(b"NM_100000", false), subject(b"NM_100001", true)];
    // The report's reads as its writer makes them (`PairwiseReads`): both description lines,
    // both alignments.
    let reads =
        crate::report::db_subject::ReportKeys::Pairwise(crate::report::db_subject::PairwiseReads {
            hits: vec![vec![0], vec![1]],
            described: vec![0, 1],
            aligned: 2,
        });
    let loads = HitLoads::of(
        crate::report::db_subject::LoadsBefore::NONE,
        Some(RunThrow::every_row(AfterThrow::SkipsAlignment {
            html: false,
        })),
        [(true, reads)],
        |key| parsed[key].as_deref(),
        |_, _| 0,
    );
    assert!(loads.loaded(0));
    assert!(loads.thrown(1).is_some());
    let mut second = tblastx_hit(1);
    second.hit.s_idx = 1;
    let mut third = second.clone();
    third.hit.bit_score = 13.9;
    let first = tblastx_hit(1);
    let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
        [(0, vec![&first]), (1, vec![&second, &third])]
            .into_iter()
            .collect();
    let titles: std::collections::HashMap<u32, Vec<u8>> =
        [(0, b"ntitle".to_vec()), (1, b"ntitle".to_vec())]
            .into_iter()
            .collect();
    let mut out = Vec::new();
    write_blastn_description_table(
        &mut out,
        &[0, 1],
        true,
        &hits,
        &titles,
        &parsed,
        &[],
        false,
        false,
        None,
        &mut crate::report::query_warnings::no_places,
        None,
        &loads,
    )
    .unwrap();
    let rows: Vec<&[u8]> = out.split(|&byte| byte == b'\n').skip(3).take(2).collect();
    assert!(
        rows[0].starts_with(b"NM_100000.1 ntitle "),
        "{:?}",
        String::from_utf8_lossy(rows[0])
    );
    assert!(
        rows[1].starts_with(b"NM_100001.1 Unknown "),
        "{:?}",
        String::from_utf8_lossy(rows[1])
    );
    // With `-html` the row's link loads the subject again in the `catch` (`HtmlRows::end`):
    // the run ends at the second row with the load's exception, after the first row.
    let registry =
        crate::report::formatter_registry::InitConfig::of("blastn", None, false).unwrap();
    let links = Arc::new(crate::report::html_ids::HtmlLinks {
        registry,
        db_name: Vec::new(),
        nucleotide: true,
        ncbirc: Default::default(),
    });
    let report = crate::report::html_ids::HtmlReport::new(&links, []);
    let mut html = Vec::new();
    let error = write_blastn_description_table(
        &mut html,
        &[0, 1],
        true,
        &hits,
        &titles,
        &parsed,
        &[],
        false,
        false,
        None,
        &mut crate::report::query_warnings::no_places,
        Some((&report, AfterThrow::SkipsAlignment { html: true })),
        &loads,
    )
    .unwrap_err();
    let exception = crate::ncbi_exception::NcbiException::of_io_error(&error).unwrap();
    assert_eq!(exception.what(), loads.thrown(1).unwrap().what());
    // The two header lines, the blank line and the first row; nothing of the second.
    let text = String::from_utf8(html).unwrap();
    let lines: Vec<&str> = text.split_terminator('\n').collect();
    assert_eq!(lines.len(), 4, "{text:?}");
    assert!(lines[3].contains("NM_100000.1"), "{text:?}");
    assert!(!text.contains("NM_100001"), "{text:?}");
    let mut out = Vec::new();
    write_database_subject_alignments_part(
        &mut out,
        None,
        DatabaseHeading::loaded(parsed[1].as_deref().unwrap(), false, 2, true),
        false,
        &hits[&1],
        &[1, 2],
        0..2,
        None,
        &mut crate::report::query_warnings::no_places,
        |_, _, _| panic!("no HSP is written"),
    )
    .unwrap();
    assert_eq!(
        out,
        b"Sequence with id NM_100001 no longer exists in database...alignment skipped\n\
          Sequence with id NM_100001 no longer exists in database...alignment skipped\n"
    );
}

/// A BLASTP pairwise report written a window of queries at a time (the rows of a caller that
/// makes them window by window) is the report written at once: the prolog with the first
/// window, the epilog with the last, the query indices global. So is a psiblast round's report
/// (`psiblast_round`, the prolog with the first window, no epilog), also in two parts.
#[test]
fn a_blastp_report_written_in_windows_is_the_report_written_at_once() {
    let karlin = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        ..KarlinParams::default()
    };
    let query = |length: usize| BlastpPairwiseQuery {
        query_name: Vec::new(),
        query_length: length,
        ungapped_karlin: karlin,
        effective_search_space: 1_000_000,
        valid: true,
        batch_skipped: false,
    };
    let queries = [query(30), query(40), query(50), query(60)];
    let query_titles: Vec<Arc<[u8]>> = (0..4)
        .map(|i| Arc::from(format!("Query_{i}").into_bytes().into_boxed_slice()))
        .collect();
    let subject_titles: Vec<Arc<[u8]>> = vec![Arc::from(&b"gnl|CDD|1 profile"[..])];
    let parsed: Vec<Option<Arc<crate::report::db_subject::ParsedSubject>>> = vec![None];
    let hit = |q_idx: u32| {
        let mut pairwise = PairwiseHit::from(Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value: 1e-5,
            bit_score: 25.0,
            num_ident: 10,
            query_frame: 0,
            query_length: 30,
            q_idx,
            s_idx: 0,
            raw_score: 50,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 10,
        });
        pairwise.query_seq = Some("MKVLAAGIVG".to_string());
        pairwise.subject_seq = Some("MKVLAAGIVG".to_string());
        pairwise.subject_length = Some(180);
        pairwise.subject_title = Some("profile".to_string());
        pairwise
    };
    // Queries 0 and 2 have no hit.
    let hits = [hit(1), hit(3)];
    let report = BlastpPairwiseReport {
        version: "BLASTP 2.17.0+".to_string(),
        database_name: b"db".to_vec(),
        database_num_sequences: 1,
        database_total_letters: 180,
        matrix_name: Some("BLOSUM62".to_string()),
        gap_open: 11,
        gap_extend: 1,
        word_threshold: 0.0,
        window_size: 40,
        gapped_karlin: karlin,
        gumbel: None,
        gapped: true,
        show_linked_set_size: false,
        num_descriptions: 500,
        num_alignments: 250,
    };
    let config = PairwiseConfig {
        line_length: 60,
        show_frame: false,
        program: "rpsblast".to_string(),
        protein_matrix: ScoringMatrix::Blosum62,
        sort: Default::default(),
        view: PairwiseView::Pairwise,
        unreadable_query: None,
        html: None,
    };
    let gapped_k = [0.041, 0.047, 0.041, 0.041];
    let write = |range: std::ops::Range<usize>, rows: &[PairwiseHit], out: &mut Vec<u8>| {
        write_blastp_pairwise_window(
            rows,
            range,
            out,
            &config,
            &queries,
            &query_titles,
            &subject_titles,
            &parsed,
            &[],
            &report,
            Some(&gapped_k),
            None,
            None,
            None,
            true,
            true,
            None,
            None,
            None,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
    };
    let mut at_once = Vec::new();
    write(0..4, &hits, &mut at_once);
    assert!(String::from_utf8_lossy(&at_once).contains("Query_3"));
    for windows in [
        vec![0..1, 1..4],
        vec![0..2, 2..3, 3..4],
        vec![0..1, 1..2, 2..3, 3..4],
    ] {
        let mut windowed = Vec::new();
        for window in windows {
            let rows: Vec<PairwiseHit> = hits
                .iter()
                .filter(|hit| window.contains(&(hit.hit.q_idx as usize)))
                .cloned()
                .collect();
            write(window, &rows, &mut windowed);
        }
        assert_eq!(
            String::from_utf8_lossy(&windowed),
            String::from_utf8_lossy(&at_once)
        );
    }
    let round_config = PairwiseConfig {
        program: "psiblast".to_string(),
        ..config.clone()
    };
    for found_before in [None, Some(&[true][..])] {
        let round = |range: std::ops::Range<usize>, rows: &[PairwiseHit], out: &mut Vec<u8>| {
            psiblast_round::write_psiblast_round_report(
                rows,
                range,
                out,
                &round_config,
                &queries,
                &query_titles,
                &subject_titles,
                &parsed,
                &[],
                &crate::report::db_subject::HitLoads::default(),
                &report,
                None,
                None,
                None,
                true,
                psiblast_round::ReportRound {
                    number: 2,
                    found_before,
                },
            )
            .unwrap();
        };
        let mut at_once = Vec::new();
        round(0..4, &hits, &mut at_once);
        let text = String::from_utf8_lossy(&at_once).into_owned();
        assert!(text.contains("Results from round 2") && text.contains("Query_3"));
        for windows in [vec![0..1, 1..4], vec![0..2, 2..3, 3..4]] {
            let mut windowed = Vec::new();
            for window in windows {
                let rows: Vec<PairwiseHit> = hits
                    .iter()
                    .filter(|hit| window.contains(&(hit.hit.q_idx as usize)))
                    .cloned()
                    .collect();
                round(window, &rows, &mut windowed);
            }
            assert_eq!(String::from_utf8_lossy(&windowed), text);
        }
    }
}

// NCBI BLAST+ 2.17.0 psiblast -outfmt 0 with a PSSM query of the molecule `other` (rows
// db.ps.inmsa.f.fa.mol.other.0 and db.ps.inmsa.f.fa.mol.other.id.0): the description table is
// as ever, and each HSP is one line `Sequence with id <id> no longer exists in
// database...alignment skipped` with no heading and no score lines; the id is the Seq-align
// id's text: `BL_ORD_ID:<oid>` of a database sequence without ids, `Subject_N` of a `-subject`
// record, in every round.
#[test]
fn an_unreadable_query_writes_one_line_for_each_hsp_with_the_subjects_id() {
    let karlin = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        ..KarlinParams::default()
    };
    let queries = [BlastpPairwiseQuery {
        query_name: Vec::new(),
        query_length: 61,
        ungapped_karlin: karlin,
        effective_search_space: 1_000_000,
        valid: true,
        batch_skipped: false,
    }];
    let query_titles: Vec<Arc<[u8]>> = vec![Arc::from(&b"q"[..])];
    let subject_titles: Vec<Arc<[u8]>> = vec![Arc::from(&b"first"[..]), Arc::from(&b"second"[..])];
    let parsed: Vec<Option<Arc<crate::report::db_subject::ParsedSubject>>> = vec![None, None];
    let hit = |s_idx: u32, bit_score: f64| {
        let mut pairwise = PairwiseHit::from(Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value: 1e-5,
            bit_score,
            num_ident: 10,
            query_frame: 0,
            query_length: 61,
            q_idx: 0,
            s_idx,
            raw_score: 50,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 10,
        });
        pairwise.query_seq = Some("MKVLAAGIVG".to_string());
        pairwise.subject_seq = Some("MKVLAAGIVG".to_string());
        pairwise.subject_length = Some(185);
        pairwise
    };
    let hits = [hit(0, 120.5), hit(0, 80.5), hit(1, 100.5)];
    let report = BlastpPairwiseReport {
        version: "PSIBLAST (test)".to_string(),
        database_name: b"db".to_vec(),
        database_num_sequences: 2,
        database_total_letters: 370,
        matrix_name: Some("BLOSUM62".to_string()),
        gap_open: 11,
        gap_extend: 1,
        word_threshold: 0.0,
        window_size: 40,
        gapped_karlin: karlin,
        gumbel: None,
        gapped: true,
        show_linked_set_size: false,
        num_descriptions: 500,
        num_alignments: 250,
    };
    let write = |unreadable_query: Option<UnreadableQuery>, number: u32| {
        let config = PairwiseConfig {
            program: "psiblast".to_string(),
            unreadable_query,
            ..PairwiseConfig::default()
        };
        let mut out = Vec::new();
        psiblast_round::write_psiblast_round_report(
            &hits,
            0..1,
            &mut out,
            &config,
            &queries,
            &query_titles,
            &subject_titles,
            &parsed,
            &[],
            &crate::report::db_subject::HitLoads::default(),
            &report,
            None,
            None,
            None,
            true,
            psiblast_round::ReportRound {
                number,
                found_before: None,
            },
        )
        .unwrap();
        String::from_utf8(out).unwrap()
    };
    let lines = |text: &str| -> Vec<String> {
        text.lines()
            .filter(|line| line.starts_with("Sequence with id "))
            .map(str::to_string)
            .collect()
    };
    let skipped = |ids: [&str; 3]| -> Vec<String> {
        ids.map(|id| {
            format!("Sequence with id {id} no longer exists in database...alignment skipped")
        })
        .to_vec()
    };
    // A `-subject` set: the record's position.
    let records = write(Some(UnreadableQuery::default()), 1);
    assert_eq!(
        lines(&records),
        skipped(["Subject_1", "Subject_1", "Subject_2"])
    );
    // A database's sequences without Seq-ids: their OIDs, which are not their ranks.
    let database = write(
        Some(UnreadableQuery {
            oids: Some(Arc::from([45_u32, 481])),
        }),
        2,
    );
    assert_eq!(
        lines(&database),
        skipped(["BL_ORD_ID:45", "BL_ORD_ID:45", "BL_ORD_ID:481"])
    );
    // Nothing else of the alignments is written: no heading, no score line, no alignment.
    for text in [&records, &database] {
        assert!(text.contains("Results from round "), "{text}");
        assert!(
            !text.contains("Length=185") && !text.contains(" Score = "),
            "{text}"
        );
        assert!(!text.contains("Query  "), "{text}");
    }
    // A readable query writes the alignments as before.
    let readable = write(None, 1);
    assert!(
        lines(&readable).is_empty() && readable.contains(" Score = "),
        "{readable}"
    );
}

// NCBI BLAST+ 2.17.0 psiblast two-part description table (`psiblast_round::write_part`) uses the
// same shown-HSP rule as the main description table (`shown_hsp` / `table_scores`): start from 0.0,
// keep the first HSP with a strictly higher bit score, and show 0.0 / 0.0 when every HSP of a
// subject has bit score <= 0 or NaN (PR #1638 review item 1).
#[test]
fn psiblast_round_two_part_table_shows_zero_when_all_bit_scores_are_negative() {
    let karlin = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        ..KarlinParams::default()
    };
    let queries = [BlastpPairwiseQuery {
        query_name: Vec::new(),
        query_length: 61,
        ungapped_karlin: karlin,
        effective_search_space: 1_000_000,
        valid: true,
        batch_skipped: false,
    }];
    let query_titles: Vec<Arc<[u8]>> = vec![Arc::from(&b"q"[..])];
    let subject_titles: Vec<Arc<[u8]>> = vec![
        Arc::from(&b"neg_repeat"[..]),
        Arc::from(&b"neg_new"[..]),
        Arc::from(&b"pos_new"[..]),
    ];
    let parsed: Vec<Option<Arc<crate::report::db_subject::ParsedSubject>>> = vec![None, None, None];
    let hit = |s_idx: u32, bit_score: f64, e_value: f64| {
        let mut pairwise = PairwiseHit::from(Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value,
            bit_score,
            num_ident: 10,
            query_frame: 0,
            query_length: 61,
            q_idx: 0,
            s_idx,
            raw_score: 50,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 10,
        });
        pairwise.query_seq = Some("MKVLAAGIVG".to_string());
        pairwise.subject_seq = Some("MKVLAAGIVG".to_string());
        pairwise.subject_length = Some(185);
        pairwise
    };
    let hits = [
        hit(0, -1.0, 1e-5),
        hit(0, -1.0, 2e-3),
        hit(1, -1.0, 3e-4),
        hit(2, 42.5, 1e-9),
    ];
    let report = BlastpPairwiseReport {
        version: "PSIBLAST (test)".to_string(),
        database_name: b"db".to_vec(),
        database_num_sequences: 3,
        database_total_letters: 555,
        matrix_name: Some("BLOSUM62".to_string()),
        gap_open: 11,
        gap_extend: 1,
        word_threshold: 0.0,
        window_size: 40,
        gapped_karlin: karlin,
        gumbel: None,
        gapped: true,
        show_linked_set_size: false,
        num_descriptions: 500,
        num_alignments: 250,
    };
    let config = PairwiseConfig {
        program: "psiblast".to_string(),
        ..PairwiseConfig::default()
    };
    let mut out = Vec::new();
    psiblast_round::write_psiblast_round_report(
        &hits,
        0..1,
        &mut out,
        &config,
        &queries,
        &query_titles,
        &subject_titles,
        &parsed,
        &[],
        &crate::report::db_subject::HitLoads::default(),
        &report,
        None,
        None,
        None,
        true,
        psiblast_round::ReportRound {
            number: 2,
            found_before: Some(&[true, false, false]),
        },
    )
    .unwrap();
    let text = String::from_utf8(out).unwrap();
    let row = |name: &str| -> String {
        text.lines()
            .find(|line| line.starts_with(name))
            .unwrap_or_else(|| panic!("missing row for {name} in:\n{text}"))
            .to_string()
    };
    fn tokens(line: &str) -> Vec<&str> {
        line.split_whitespace().collect()
    }
    assert_eq!(tokens(&row("neg_repeat")), vec!["neg_repeat", "0.0", "0.0"]);
    assert_eq!(tokens(&row("neg_new")), vec!["neg_new", "0.0", "0.0"]);
    assert_eq!(tokens(&row("pos_new")), vec!["pos_new", "42.5", "1e-09"]);
}

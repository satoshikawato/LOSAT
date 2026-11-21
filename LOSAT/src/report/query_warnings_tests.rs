use crate::blastinput::fasta_reader::FastaRecord;

#[test]
fn messages_of_a_query_are_in_ncbis_sorted_order() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1043-1054 sorts a query's messages (same error id and
    // severity, so by text) and removes repeats.
    let query = FastaRecord::new("Query_1", b"q1", b"OOO");
    let o = replaced_o_message(b"OOO").unwrap();
    let line = query_warning(
        "blastp",
        0,
        &query,
        &[o.clone(), INVALID_QUERY_MESSAGE.to_string(), o.clone()],
    );
    let text = String::from_utf8(line).unwrap();
    assert!(
        text.starts_with("Warning: [blastp] Query_1 q1: Could not calculate"),
        "{text}"
    );
    assert_eq!(text.matches("One or more O").count(), 1, "{text}");
}

use super::*;

// NCBI 2.17.0 (M4a-1b probes s19, s20; M4a-1c E): the list form of `GetSequenceProtein`'s errors.
#[test]
fn the_invalid_residues_error_lists_the_first_twenty_positions() {
    assert_eq!(invalid_residues_error(&[]), None);
    assert_eq!(
        invalid_residues_error(&[26]).unwrap(),
        "Warning: Invalid residues found at positions 26 "
    );
    let twenty: Vec<usize> = (0..20).collect();
    let listed = (0..20)
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    assert_eq!(
        invalid_residues_error(&twenty).unwrap(),
        format!("Warning: Invalid residues found at positions {listed} ")
    );
    let more: Vec<usize> = (0..22).collect();
    assert_eq!(
        invalid_residues_error(&more).unwrap(),
        format!("Warning: Invalid residues found at positions {listed},... (only first 20 shown) ")
    );
}

// Probes C31 and D01 of line DI-8: NCBI adds the empty-database message to the sorted
// messages of the set-up, though "Sequence ..." sorts after "Filtering ...".
#[test]
fn the_empty_database_message_comes_after_the_set_up_messages() {
    let messages = [
        EMPTY_DATABASE_MESSAGE.to_string(),
        "Sequence contains no data".to_string(),
        EMPTY_DATABASE_MESSAGE.to_string(),
    ];
    let query = FastaRecord::new("Query_2", b"empty", b"");
    assert_eq!(
        query_warning("blastp", 1, &query, &messages),
        b"Warning: [blastp] Query_2 empty: Sequence contains no data Filtering resulted in an empty database. \n"
    );
}

// NCBI 2.17.0 psiblast `-in_pssm` (the review of PS4a, S2 and N3; `PS/tmp/fix4a` probes
// s2.idO, n3.localid_O, n3.localneg_O, n3.localid_notitle_O): the label is the PSSM's Seq-id,
// whatever the query's place.
#[test]
fn a_pssm_querys_warning_starts_with_its_seq_id() {
    let o = replaced_o_message(b"MKVOAAOK").unwrap();
    let title = b"YP_009220601.1 TK-TMK, chimeric thymidine and thymidylate kinase";
    for (seq_id, title, shown) in [
        (
            &b"myq"[..],
            &title[..],
            &b"myq YP_009220601.1 TK-TMK.. "[..],
        ),
        (b"5", title, b"5 YP_009220601.1 TK-TMK, .. "),
        (b"-3", title, b"-3 YP_009220601.1 TK-TMK,.. "),
        (b"5", b"", b"5"),
    ] {
        let mut expected = b"Warning: [psiblast] ".to_vec();
        expected.extend_from_slice(shown);
        expected.extend_from_slice(b": ");
        expected.extend_from_slice(o.as_bytes());
        expected.extend_from_slice(b" \n");
        assert_eq!(
            query_warning_with_seq_id("psiblast", seq_id, title, std::slice::from_ref(&o)),
            expected
        );
    }
    assert!(query_warning_with_seq_id("psiblast", b"5", b"", &[]).is_empty());
}

// Query_1 + FASTA title is shortened after 35 bytes, then the invalid-Karlin
// warning keeps NCBI's trailing space and newline.
#[test]
fn invalid_query_warning_matches_ncbi_bytes() {
    let short = FastaRecord::new("Query_1", b"nohit_query", b"W");
    assert_eq!(
        invalid_query_warning("tblastn", 0, &short),
        b"Warning: [tblastn] Query_1 nohit_query: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n"
    );
    let long = FastaRecord::new("Query_2", b"long_header abcdefghijklmnopqrstuvwxyz", b"W");
    assert!(invalid_query_warning("blastn", 1, &long)
        .starts_with(b"Warning: [blastn] Query_2 long_header abcde.. : "));
}

// NCBI BLAST+ 2.17.0, the invalid-query warning of BLASTN for one query with each
// defline (inventory range RP, `scratch_RP/warn/iq_*_nuc_q.fa` and
// `warn/out/iq_*_blastn_6.err`): the query ID between `[blastn] ` and `: Could`.
#[test]
fn warning_query_ids_follow_the_ncbi_oracle() {
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
    let w = |count: usize| vec![b'w'; count];
    let long = [&b"id "[..], &w(100)].concat();
    let cjk = [&[b'a'; 23][..], b"\xe6\xbc\xa2", &[b'z'; 14]].concat();
    let z27 = vec![b'z'; 27];
    let z28 = vec![b'z'; 28];
    for (defline, shown) in [
        (&b"sid sdesc"[..], &b"Query_1 sid sdesc"[..]),
        (b"", b"Query_1"),
        (b"id \xc3\xa9 x", b"Query_1 id \xc3\xa9 x"),
        (b"id \xe9 x", b"Query_1 id \xe9 x"),
        (b"id\tx", b"Query_1 id"),
        (&long, b"Query_1 id wwwwwwwwwwwwww.. "),
        (&cjk, b"Query_1 aaaaaaaaaaaaaaaaa.. "),
        (b"id x\xff", b"Query_1 id x\xff"),
        (&z27, b"Query_1 zzzzzzzzzzzzzzzzzzzzzzzzzzz"),
        (&z28, b"Query_1 zzzzzzzzzzzzzzzzz.. "),
    ] {
        let bytes = [&b">"[..], defline, b"\nACGTACGT\n"].concat();
        let mut source =
            FastaInputSource::from_bytes(&bytes, ReaderConfig::query("BLASTN", false, false));
        let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
        assert_eq!(
            warning_query_id(records[0].local_id.as_bytes(), &records[0].title),
            shown,
            "{defline:?}"
        );
    }
}

/// The report (`R:`) and the warnings' texts in the order written.
#[derive(Clone, Default)]
struct Log(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

impl Write for Log {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl crate::report::ostream::ReportWrite for Log {}

#[test]
fn the_texts_of_the_reads_follow_the_report_up_to_each_read() {
    let reads: Vec<QueryReads> = vec![
        vec![
            SubjectRead::of_text(5, b"A\n"),
            SubjectRead::of_text(7, b"B\n"),
            SubjectRead::of_text(8, b"C\n"),
        ],
        vec![SubjectRead::of_text(3, b"D\n")],
        vec![SubjectRead::of_text(4, b"E\n")],
    ];
    let before: Vec<Posts> = [b"W0\n", b"W1\n", b"W2\n"]
        .map(|lines| Posts::of_lines(lines.to_vec()))
        .into();
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut warnings = QueryWarnings::of_posts(&before, None, &mut sink).with_subject_reads(&reads);
    assert!(warnings.reads_with_texts(0) && !warnings.reads_with_texts(3));
    warnings.before_query(0, &mut report).unwrap();
    // A subject without a text writes nothing; a read writes the texts up to its own.
    warnings.subject_read(0, 9, &mut report).unwrap();
    report.write_all(b"R1\n").unwrap();
    warnings.subject_read(0, 7, &mut report).unwrap();
    report.write_all(b"R2\n").unwrap();
    // Read before: nothing again.
    warnings.subject_read(0, 5, &mut report).unwrap();
    warnings.subject_read(0, 7, &mut report).unwrap();
    // C is never read: it comes after query 0's report, before query 1's warnings.
    warnings.before_query(1, &mut report).unwrap();
    report.write_all(b"R3\n").unwrap();
    warnings.subject_read(1, 3, &mut report).unwrap();
    // Query 2's text, never read, comes after the last report.
    warnings.before_query(2, &mut report).unwrap();
    warnings.after_last_query(3, &mut report).unwrap();
    drop(warnings);
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    assert_eq!(text, "W0\nR1\nA\nB\nR2\nC\nW1\nR3\nD\nW2\nE\n");
}

#[test]
fn the_posts_of_a_hit_follow_the_report_up_to_the_place_of_its_load() {
    // XML2 and JSON: every `Hit` loads its subject at the report's start, after outfmt 15's
    // separator; the report written before is flushed first (a buffered output stream).
    let log = Log::default();
    let mut report = io::BufWriter::with_capacity(1 << 16, log.clone());
    let mut sink = log.clone();
    report.write_all(b"{0}\n").unwrap();
    let mut posts = HitPosts::new(false);
    posts.hit(2, 1, |texts| texts.write_all(b"T1\n")).unwrap();
    posts.hit(2, 2, |_| Ok(())).unwrap();
    posts.hit(2, 3, |texts| texts.write_all(b"T3\n")).unwrap();
    posts
        .write(b",\n{1}\n", &mut report, &mut sink, |_, _, _| Ok(()))
        .unwrap();
    report.flush().unwrap();
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    assert_eq!(text, "{0}\n,\nT1\nT3\n{1}\n");
    // BLAST XML: each `Hit` loads its subject at its own start; a query whose model's reads
    // write texts marks every `Hit`, and a `Hit`'s read comes before the texts made there.
    let reads: Vec<QueryReads> = vec![vec![SubjectRead::of_text(5, b"A\n")]];
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut warnings = QueryWarnings::of_posts(&[], None, &mut sink).with_subject_reads(&reads);
    let mut posts = HitPosts::new(warnings.reads_with_texts(0));
    posts.hit(0, 5, |texts| texts.write_all(b"P5\n")).unwrap();
    posts.hit(4, 6, |texts| texts.write_all(b"P6\n")).unwrap();
    write_hit_posts(b"<H5><H6>", &posts, 0, &mut report, Some(&mut warnings)).unwrap();
    drop(warnings);
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    assert_eq!(text, "A\nP5\n<H5>P6\n<H6>");
    // Without warnings, texts that would be lost are LOSAT's internal error; marks alone
    // write the report.
    let mut alone = Vec::new();
    let lost = write_hit_posts(b"<H5><H6>", &posts, 0, &mut alone, None).unwrap_err();
    assert!(
        lost.to_string().starts_with("LOSAT internal error"),
        "{lost}"
    );
    let mut marks = HitPosts::new(true);
    marks.hit(0, 5, |_| Ok(())).unwrap();
    let mut alone = Vec::new();
    write_hit_posts(b"<H5><H6>", &marks, 0, &mut alone, None).unwrap();
    assert_eq!(alone, b"<H5><H6>");
}

/// A JSON report that raised an exception its report function catches (line DBID pfoff-1: the
/// `Hit` of a subject whose lookup throws): the separator before it stays, the texts of the
/// loads of the `Hit`s before it come after it, the read of the subject that throws and those
/// after it are not written (`cut`: NCBI loads none of them), and the exception's post comes
/// last (blastp 2.17.0, outfmt 15 of `pP pU`, `R/mixdb2/review1/out/ncbi/bp.lim1.U2.15`).
#[test]
fn the_post_of_a_caught_exception_follows_the_loads_before_it() {
    let reads: Vec<QueryReads> = vec![vec![
        SubjectRead::of_text(1, b"R1\n"),
        SubjectRead::of_text(2, b"R2\n"),
    ]];
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut warnings = QueryWarnings::of_posts(&[], None, &mut sink).with_subject_reads(&reads);
    let mut posts = HitPosts::new(warnings.reads_with_texts(0));
    posts.hit(2, 1, |texts| texts.write_all(b"T1\n")).unwrap();
    posts.document(crate::report::blast_xml::DocumentPosts {
        cut: true,
        ..Default::default()
    });
    posts.caught(b"Error: [blastp] E\n".to_vec());
    write_hit_posts(b",\n", &posts, 0, &mut report, Some(&mut warnings)).unwrap();
    drop(warnings);
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    assert_eq!(text, ",\nR1\nT1\nError: [blastp] E\n");
    // Without warnings the post would be lost: LOSAT's internal error.
    let mut alone = HitPosts::new(false);
    alone.caught(b"Error: [blastp] E\n".to_vec());
    assert!(write_hit_posts(b"", &alone, 0, &mut Vec::new(), None).is_err());
}

/// A read with what its load looks up (`ReadLoad`), the run's texts `texts`.
fn loaded(
    texts: &std::sync::Arc<LoadTexts>,
    subject: u32,
    text: &[u8],
    (volume, offset, opening): (u32, bool, usize),
) -> SubjectRead {
    SubjectRead {
        subject,
        text: text.to_vec(),
        load: Some(ReadLoad {
            texts: std::sync::Arc::clone(texts),
            volume,
            looks_up: true,
            offset,
            opening,
        }),
    }
}

/// A BLAST XML report that made no `Hit` after its last mark (`DocumentPosts::cut`: its
/// `Iteration` or a `Hit` threw) loads none of its later subjects: the texts of their loads go
/// to NCBI's next load that makes the same ask, of whichever subject (NCBI 2.17.0, review 1 of
/// K4h M2, `q3w.faa` with `PRE_FETCH_SEQS_LIMIT=1`: the offset error at query 3's `Hit`, not
/// after query 2's post; agent3 probes m2.*): the texts of opening at the next load that looks
/// a taxid up, a volume's offset error at its next load that looks that taxid up. A text that
/// no later load asks for is never written, and a run that NCBI's XML stream ends writes no
/// text after its posts.
#[test]
fn a_cut_report_leaves_its_texts_to_the_next_load_that_asks() {
    let texts = std::sync::Arc::new(LoadTexts::default());
    // Query 0 reads 5 (the run's first ask: the texts of opening, then volume 0's offset
    // error) and 6 (volume 1's offset error); query 1 reads 9 (volume 1, the offset taxid)
    // and 6 again; query 2 reads 7 (volume 0, another taxid).
    let reads: Vec<QueryReads> = vec![
        vec![
            loaded(&texts, 5, b"O\nE0\n", (0, true, 2)),
            loaded(&texts, 6, b"E1\n", (1, true, 0)),
        ],
        vec![
            loaded(&texts, 9, b"", (1, true, 0)),
            loaded(&texts, 6, b"", (1, true, 0)),
        ],
        vec![loaded(&texts, 7, b"", (0, false, 0))],
    ];
    // Query 0 reads 5, then a `Hit` throws: 6's offset error goes to 9's load.
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut cursor = ReadCursor::new(&reads).unwrap();
    cursor.subject_read(0, 5, &mut report, &mut sink).unwrap();
    cursor.cut(0, &mut report, &mut sink).unwrap();
    cursor.before_query(1, &mut report, &mut sink).unwrap();
    assert!(cursor.has_texts(1));
    report.write_all(b"<Q1>").unwrap();
    cursor.subject_read(1, 9, &mut report, &mut sink).unwrap();
    cursor.subject_read(1, 6, &mut report, &mut sink).unwrap();
    report.write_all(b"<H6>").unwrap();
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    assert_eq!(text, "O\nE0\n<Q1>E1\n<H6>");
    // Query 0's `Iteration` throws: its texts of opening and volume 1's offset error go to 9's
    // load, volume 0's offset error waits for a load that looks its taxid up (7 does not).
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut cursor = ReadCursor::new(&reads).unwrap();
    cursor.cut(0, &mut report, &mut sink).unwrap();
    report.write_all(b"<Q1>").unwrap();
    cursor.subject_read(1, 9, &mut report, &mut sink).unwrap();
    cursor.before_query(2, &mut report, &mut sink).unwrap();
    cursor.subject_read(2, 7, &mut report, &mut sink).unwrap();
    assert!(cursor.unwritten().is_empty());
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    assert_eq!(text, "<Q1>O\nE1\n");
    // A batch modelled after that report: its first read of volume 0 that looks the taxid up
    // makes volume 0's ask, whose text then counts as its own (cut again, it moves again).
    let mut load = reads[0][0].load.clone().unwrap();
    load.opening = 0;
    assert_eq!(load.place(b"X\n".to_vec()), b"X\nE0\n");
    assert_eq!(load.opening, 0);
    load.left_out(b"X\nE0\n");
    let next = loaded(&texts, 8, b"", (0, true, 0));
    assert_eq!(next.placed(), b"X\nE0\n");
    texts.clear();
    // The end of the run: nothing after it.
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut cursor = ReadCursor::new(&reads).unwrap();
    cursor.cut(0, &mut report, &mut sink).unwrap();
    cursor.end_run();
    cursor.subject_read(1, 9, &mut report, &mut sink).unwrap();
    cursor.before_query(2, &mut report, &mut sink).unwrap();
    assert!(cursor.unwritten().is_empty() && !cursor.has_texts(1));
    assert!(log.0.lock().unwrap().is_empty());
}

/// NCBI's posts of the bytes replaced in a batch's search strategy request follow the reader's
/// messages of the batch in its first query and come before the rest of its lines (NCBI
/// 2.17.0, probe RG K/k4b/probe2 `blastp_b3_f6`: the O warning of query 3 after the posts of
/// its batch); a batch without posts keeps its lines, and the posts are whole posts.
/// A batch whose prefetch is on, modelled before the report of a batch before it that the
/// stream cut (a program that models every batch before it writes): its loads' texts are made
/// where the posts are written (`Posts::loads`), after that report, so the texts that the cut
/// left go to the prefetch's loads that make their asks, in the order of the loads, after the
/// batch's reader messages and request posts (NCBI 2.17.0, agent4 m2 probe `*.qm.*.mix.f5`:
/// the offset error of query 2's prefetch after query 1's `Failed s_BlastXMLAddIteration`).
#[test]
fn a_prefetch_takes_the_texts_that_a_cut_report_before_it_left() {
    let texts = std::sync::Arc::new(LoadTexts::default());
    // Query 0's batch reads 5 (the texts of opening, then volume 0's offset error) one by one;
    // query 1's batch prefetches 6 (volume 1, its offset error) and 7 (volume 0, the offset
    // taxid), whose model gave 7 no text: 5 asked before it.
    let reads: Vec<QueryReads> = vec![vec![loaded(&texts, 5, b"O\nE0\n", (0, true, 2))]];
    let mut prefetch = Posts::new();
    prefetch.loads(vec![
        loaded(&texts, 6, b"E1\n", (1, true, 0)),
        loaded(&texts, 7, b"", (0, true, 0)),
    ]);
    let mut request = Posts::new();
    request.post(crate::report::diag_log::ERROR, b"P1\n".to_vec());
    let posts = query_posts(
        vec![Vec::new(), b"R1\nW1\n".to_vec()],
        vec![(1, request)],
        vec![(1, prefetch)],
        |query| [0, 3][query],
    );
    let log = Log::default();
    let mut report = log.clone();
    let mut sink = log.clone();
    let mut warnings = QueryWarnings::of_posts(&posts, None, &mut sink).with_subject_reads(&reads);
    warnings.before_query(0, &mut report).unwrap();
    // Query 0's report is cut before its read of 5.
    let (cursor, sink) = warnings.reads_and_sink();
    cursor.unwrap().cut(0, &mut report, sink).unwrap();
    warnings.before_query(1, &mut report).unwrap();
    drop(warnings);
    let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
    // The opening at 6, the first load that looks a taxid up; volume 0's offset error at 7.
    assert_eq!(text, "R1\nP1\nO\nE1\nE0\nW1\n");
    // Without a cut, the model's texts.
    let texts = std::sync::Arc::new(LoadTexts::default());
    let mut prefetch = Posts::new();
    prefetch.loads(vec![
        loaded(&texts, 6, b"E1\n", (1, true, 0)),
        loaded(&texts, 7, b"", (0, true, 0)),
    ]);
    assert_eq!(prefetch.text(), b"E1\n");
    // A run whose lookups write no text keeps the prefetch's lines as they are.
    let mut plain = Posts::new();
    plain.loads(vec![SubjectRead::of_text(6, b"T\n")]);
    assert_eq!(plain, Posts::of_lines(b"T\n".to_vec()));
}

/// The ledger that `LoadTexts::left_before` models for a stopping batch's prefetch is the one
/// that the writers leave where they have written the reports and posts before it (review 5
/// of line DBID loadord, M1): a cut report's texts, less what the posts and reads after it take
/// (the opening at the first load that looks a taxid up, a volume's offset error at its next
/// load of that taxid); a load of the stop placed from it writes what one placed after the
/// writers writes.
#[test]
fn the_texts_left_before_a_batch_are_those_the_writers_leave() {
    // Query 0's report (cut) reads 5 (the texts of opening, volume 0's offset error) and 6
    // (volume 1's); query 1's batch prefetches 7 (volume 1, the offset taxid: it takes the
    // opening and volume 1's error); query 2's report reads 8 (volume 0, another taxid).
    let shapes = |texts: &std::sync::Arc<LoadTexts>| {
        let reads: Vec<QueryReads> = vec![
            vec![
                loaded(texts, 5, b"O\nE0\n", (0, true, 2)),
                loaded(texts, 6, b"E1\n", (1, true, 0)),
            ],
            Vec::new(),
            vec![loaded(texts, 8, b"", (0, false, 0))],
        ];
        let mut prefetch = Posts::new();
        prefetch.loads(vec![loaded(texts, 7, b"", (1, true, 0))]);
        (reads, prefetch)
    };
    let pending = |texts: &LoadTexts| {
        let pending = texts.pending.lock().unwrap();
        (pending.opening.clone(), pending.offsets.clone())
    };
    // The writers: query 0 cut, then query 1's posts, then query 2's report.
    let texts = std::sync::Arc::new(LoadTexts::default());
    let (reads, prefetch) = shapes(&texts);
    let posts = query_posts(
        vec![Vec::new(), Vec::new(), Vec::new()],
        Vec::new(),
        vec![(1, prefetch.clone())],
        |_| 0,
    );
    let log = Log::default();
    let (mut report, mut sink) = (log.clone(), log.clone());
    let mut warnings = QueryWarnings::of_posts(&posts, None, &mut sink).with_subject_reads(&reads);
    warnings.before_query(0, &mut report).unwrap();
    let (cursor, sink) = warnings.reads_and_sink();
    cursor.unwrap().cut(0, &mut report, sink).unwrap();
    warnings.before_query(1, &mut report).unwrap();
    warnings.subject_read(2, 8, &mut report).unwrap();
    warnings.before_query(3, &mut report).unwrap();
    drop(warnings);
    assert_eq!(log.0.lock().unwrap().as_slice(), b"O\nE1\n");
    let written = pending(&texts);
    assert_eq!(written, (Vec::new(), vec![(0, b"E0\n".to_vec())]));
    // The model of the same, before any of it is written, in a run's own ledger.
    let model = std::sync::Arc::new(LoadTexts::default());
    let (reads, prefetch) = shapes(&model);
    let prefetched = |query: usize| (query == 1).then_some(&prefetch);
    let left = LoadTexts::left_before(&reads, &|query| query == 0, &prefetched, 3);
    assert_eq!(pending(&left), written);
    assert!(pending(&model) == (Vec::new(), Vec::new()) && !left.is_empty());
    // The stop's load of volume 0 that looks the taxid up writes volume 0's error.
    let mut stop = loaded(&model, 9, b"", (0, true, 0));
    stop.place_left(&left);
    assert_eq!(stop.text, b"E0\n");
    assert!(left.is_empty());
    // Up to query 1 only (a stop in query 1's batch): its prefetch is that batch's own.
    let left = LoadTexts::left_before(&reads, &|query| query == 0, &prefetched, 1);
    assert_eq!(
        pending(&left),
        (
            b"O\n".to_vec(),
            vec![(0, b"E0\n".to_vec()), (1, b"E1\n".to_vec())]
        )
    );
    // No cut: nothing left.
    assert!(LoadTexts::left_before(&reads, &|_| false, &prefetched, 3).is_empty());
}

#[test]
fn request_posts_follow_the_reader_messages_of_their_batch() {
    use crate::report::diag_log::ERROR;
    let lines = vec![
        b"R0\nW0\n".to_vec(),
        b"W1\n".to_vec(),
        b"R2\nT2\nW2\n".to_vec(),
    ];
    let post = |text: &[u8]| {
        let mut posts = Posts::new();
        posts.post(ERROR, text);
        posts
    };
    let reader = [3, 0, 3];
    let posts = query_posts(
        lines,
        vec![(0, post(b"P0\n")), (2, post(b"P2a\n")), (2, post(b"P2b\n"))],
        Vec::new(),
        |query| reader[query],
    );
    let texts: Vec<Vec<u8>> = posts.iter().map(Posts::text).collect();
    assert_eq!(
        texts,
        [
            b"R0\nP0\nW0\n".to_vec(),
            b"W1\n".to_vec(),
            b"R2\nP2a\nP2b\nT2\nW2\n".to_vec()
        ]
    );
    let mut whole = Posts::of_lines(b"R2\n".to_vec());
    whole.append(&post(b"P2a\n"));
    whole.append(&post(b"P2b\n"));
    whole.lines(b"T2\nW2\n");
    assert_eq!(posts[2], whole);
    // A request after the last query's lines keeps its posts after them.
    let posts = query_posts(
        vec![b"W0\n".to_vec()],
        vec![(1, post(b"P1\n"))],
        Vec::new(),
        |_| 0,
    );
    assert_eq!(posts[0].text(), b"W0\nP1\n");
}

#[test]
fn a_query_without_texts_reads_nothing() {
    let reads: Vec<QueryReads> = vec![Vec::new(), Vec::new()];
    let mut sink = Vec::new();
    let warnings = QueryWarnings::of_posts(&[], None, &mut sink).with_subject_reads(&reads);
    assert!(!warnings.reads_with_texts(0));
    assert!(ReadCursor::new(&reads).is_none());
}

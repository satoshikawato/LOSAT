use super::*;
use crate::common::Hit;
use std::io::Write;

fn hit(q_idx: u32, s_idx: u32, length: usize) -> PairwiseHit {
    PairwiseHit::from(Hit {
        identity: 100.0,
        length,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: length,
        s_start: 1,
        s_end: length,
        e_value: 1e-10,
        bit_score: 50.0,
        num_ident: length,
        query_frame: 1,
        query_length: 0,
        q_idx,
        s_idx,
        raw_score: 25,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: length,
    })
}

/// Marks every part: a query's head and footer, each subject's heading once (with its first
/// HSP) and each HSP; fails at the HSP of length `fail_at` after writing part of it.
struct Marks {
    num_alignments: usize,
    fail_at: Option<usize>,
    /// The subjects `(previous, s_idx)` whose hits are shown as one (`joins`).
    joined: Option<(u32, u32)>,
    /// The alignments are one section (`merges_alignments`, the query-anchored reports).
    merged: bool,
}

impl QuerySections for Marks {
    fn num_alignments(&self) -> usize {
        self.num_alignments
    }

    fn line_length(&self) -> usize {
        crate::report::pairwise::DEFAULT_LINE_LENGTH
    }

    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        self.joined == Some((previous, s_idx))
    }

    fn sort(&self) -> (crate::report::hit_sort::ReportSort, bool) {
        (crate::report::hit_sort::ReportSort::default(), false)
    }

    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&QueryGroup<'_>>,
        place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        let subjects = group.map_or(0, |group| group.subject_order.len());
        writeln!(writer, "head {q_idx} subjects {subjects}")?;
        // One line per subject, each read before it (as a description line).
        for &s_idx in group.map_or(&[][..], |group| group.subject_order.as_slice()) {
            place(ReportPlace::Load(s_idx), writer)?;
            writeln!(writer, "line {s_idx}")?;
        }
        Ok(())
    }

    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        _probe: Option<&mut FormatProbe<'_>>,
        group: &QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        _place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        if part.start == 0 {
            writeln!(writer, "subject {s_idx}")?;
        }
        for hsp in part {
            let length = group.subject_hits[&s_idx][hsp].hit.length;
            write!(writer, "  hsp {s_idx}.{hsp} length {length}")?;
            if self.fail_at == Some(length) {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "marked"));
            }
            writeln!(writer)?;
        }
        Ok(())
    }

    fn view(&self) -> crate::report::pairwise::PairwiseView {
        if self.merged {
            crate::report::pairwise::PairwiseView::QueryAnchored {
                identity: true,
                anchored: true,
            }
        } else {
            crate::report::pairwise::PairwiseView::Pairwise
        }
    }

    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &QueryGroup<'_>,
    ) -> io::Result<()> {
        writeln!(writer, "section {q_idx}")?;
        for &s_idx in group.aligned(self.num_alignments) {
            for hit in &group.subject_hits[&s_idx] {
                writeln!(writer, "  row {s_idx} length {}", hit.hit.length)?;
            }
        }
        Ok(())
    }

    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        writeln!(writer, "foot {q_idx}")
    }
}

/// Query 0: subjects 5 (3 HSPs), 2 (1), 7 (4); query 1: no hits; query 2: subject 1 (6 HSPs),
/// subject 5 (2). Lengths 10-410 (estimates 440-2040 bytes).
fn report_hits() -> Vec<PairwiseHit> {
    let mut hits = Vec::new();
    let mut length = 10;
    for (q_idx, s_idx, count) in [(0, 5, 3), (0, 2, 1), (0, 7, 4), (2, 1, 6), (2, 5, 2)] {
        for _ in 0..count {
            hits.push(hit(q_idx, s_idx, length));
            length += 25;
        }
    }
    hits
}

fn by_query(hits: &[PairwiseHit]) -> Vec<Vec<(HspIndex, &PairwiseHit)>> {
    let mut by_query = vec![Vec::new(); 3];
    for (index, hit) in hits.iter().enumerate() {
        by_query[hit.hit.q_idx as usize].push((index, hit));
    }
    by_query
}

fn titles() -> Vec<Arc<[u8]>> {
    (0..8)
        .map(|s| Arc::from(format!("s{s}").as_bytes()))
        .collect()
}

#[test]
fn pieces_of_any_size_concatenate_to_the_serial_report() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    for num_alignments in [0, 1, 2, 10] {
        let marks = Marks {
            num_alignments,
            fail_at: None,
            joined: None,
            merged: false,
        };
        let mut serial = Vec::new();
        write_serial(&mut serial, &marks, 0, &by_query, &titles, None, None).unwrap();
        let groups: Vec<_> = by_query
            .iter()
            .map(|query_hits| QueryGroup::new(query_hits, &titles, |_, _| false))
            .collect();
        for piece_bytes in [1, 500, 900, 2_000, 5_000, usize::MAX] {
            let planned = pieces(&marks, 0, &groups, piece_bytes);
            let mut joined = Vec::new();
            for &piece in &planned {
                let formatted = format_piece(&marks, groups[piece.query].as_ref(), piece, false);
                formatted.result.unwrap();
                assert!(formatted.places.is_empty());
                joined.extend_from_slice(&formatted.buffer.bytes);
            }
            assert_eq!(
                String::from_utf8(joined).unwrap(),
                String::from_utf8(serial.clone()).unwrap(),
                "num_alignments {num_alignments}, piece_bytes {piece_bytes}"
            );
            let aligned_hsps: usize = groups
                .iter()
                .flatten()
                .map(|group| {
                    group
                        .aligned(num_alignments)
                        .iter()
                        .map(|s_idx| group.subject_hits[s_idx].len())
                        .sum::<usize>()
                })
                .sum();
            let alignment_pieces = planned
                .iter()
                .filter(|piece| matches!(piece.part, Part::Alignments(..)))
                .count();
            if piece_bytes == 1 {
                // Every HSP its own piece, a subject's HSPs in several.
                assert_eq!(alignment_pieces, aligned_hsps);
            }
            if piece_bytes == usize::MAX {
                // One piece per query with aligned subjects.
                assert!(alignment_pieces <= 2);
            }
            // Three heads, two footers.
            assert_eq!(planned.len() - alignment_pieces, 5);
        }
    }
}

#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn the_parallel_writer_writes_the_serial_report_in_any_chunks() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    // Also with subject 2 shown in the hit of subject 5 before it (line DBID N3).
    for joined in [None, Some((5, 2))] {
        let marks = Marks {
            num_alignments: 10,
            fail_at: None,
            joined,
            merged: false,
        };
        let mut serial = Vec::new();
        write_serial(&mut serial, &marks, 0, &by_query, &titles, None, None).unwrap();
        let heads = String::from_utf8(serial.clone()).unwrap();
        let expected_head = if joined.is_some() {
            "head 0 subjects 2"
        } else {
            "head 0 subjects 3"
        };
        assert!(heads.contains(expected_head), "{heads}");
        parallel_reports_equal(&marks, &by_query, &titles, &serial);
    }
}

/// The parallel writer's reports of `marks` in several chunks and pieces are `serial`.
/// One log of the report and the warnings' texts in the order written.
#[derive(Clone, Default)]
struct Log(Arc<std::sync::Mutex<Vec<u8>>>);

impl ReportWrite for Log {}

impl Write for Log {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn the_texts_of_the_reads_come_where_the_serial_writer_reads_the_subjects() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    let marks = Marks {
        num_alignments: 10,
        fail_at: None,
        joined: None,
        merged: false,
    };
    // Query 0 shows subjects 5, 2 and 7, query 2 subjects 1 and 5; subject 6's text is not
    // placed (no line reads it) and comes after query 2's report.
    let reads: Vec<crate::report::query_warnings::QueryReads> = vec![
        vec![
            crate::report::query_warnings::SubjectRead::of_text(2, b"T2\n"),
            crate::report::query_warnings::SubjectRead::of_text(7, b"T7\n"),
        ],
        Vec::new(),
        vec![
            crate::report::query_warnings::SubjectRead::of_text(5, b"T5\n"),
            crate::report::query_warnings::SubjectRead::of_text(6, b"T6\n"),
        ],
    ];
    let before: Vec<crate::report::diag_log::Posts> = [b"W0\n", b"W1\n", b"W2\n"]
        .map(|lines| crate::report::diag_log::Posts::of_lines(lines.to_vec()))
        .into();
    let run = |parallel: Option<(usize, usize)>| -> String {
        let log = Log::default();
        let (mut report, mut sink) = (log.clone(), log.clone());
        let mut warnings =
            QueryWarnings::of_posts(&before, None, &mut sink).with_subject_reads(&reads);
        match parallel {
            None => write_serial(
                &mut report,
                &marks,
                0,
                &by_query,
                &titles,
                None,
                Some(&mut warnings),
            )
            .unwrap(),
            Some((budget, piece_bytes)) => {
                // The warnings' sink need not be `Send` (BLASTX's diagnostics are not).
                crate::utils::threading::x_with_search_pool_local(4, "test", |pool| {
                    write_parallel_with(
                        &mut report,
                        &marks,
                        0,
                        &by_query,
                        &titles,
                        Some(&mut warnings),
                        pool,
                        budget,
                        piece_bytes,
                    )?;
                    Ok(())
                })
                .unwrap();
            }
        }
        warnings.after_last_query(3, &mut report).unwrap();
        drop(warnings);
        let text = log.0.lock().unwrap().clone();
        String::from_utf8(text).unwrap()
    };
    let serial = run(None);
    assert!(
        serial.starts_with("W0\nhead 0 subjects 3\nline 5\nT2\nline 2\nT7\nline 7\n"),
        "{serial}"
    );
    assert!(
        serial.contains("W2\nhead 2 subjects 2\nline 1\nT5\nline 5\n"),
        "{serial}"
    );
    assert!(serial.ends_with("T6\n"), "{serial}");
    for chunks in [(1, 1), (1_000, 500), (usize::MAX, 2_000)] {
        assert_eq!(run(Some(chunks)), serial, "{chunks:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1998-2003 (cited at
// `write_serial`): the query-anchored reports load every aligned subject before the section.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn a_merged_section_is_one_piece_after_the_reads_of_its_subjects() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    for num_alignments in [0, 2, 10] {
        let marks = Marks {
            num_alignments,
            fail_at: None,
            joined: None,
            merged: true,
        };
        let mut serial = Vec::new();
        write_serial(&mut serial, &marks, 0, &by_query, &titles, None, None).unwrap();
        let text = String::from_utf8(serial.clone()).unwrap();
        assert!(!text.contains("subject "), "{text}");
        assert_eq!(text.contains("section 0"), num_alignments > 0, "{text}");
        assert!(text.contains("foot 0") && text.contains("foot 2"), "{text}");
        parallel_reports_equal(&marks, &by_query, &titles, &serial);
    }
}

fn parallel_reports_equal(
    marks: &Marks,
    by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    titles: &[Arc<[u8]>],
    serial: &[u8],
) {
    crate::utils::threading::with_search_pool(4, "test", |pool| {
        for (budget, piece_bytes) in [(1, 1), (1_000, 500), (3_000, 900), (usize::MAX, 2_000)] {
            let mut parallel = Vec::new();
            let (chunks, pieces) = write_parallel_with(
                &mut parallel,
                marks,
                0,
                by_query,
                titles,
                None,
                pool,
                budget,
                piece_bytes,
            )
            .unwrap();
            assert_eq!(
                parallel, serial,
                "budget {budget}, piece_bytes {piece_bytes}"
            );
            if budget == 1 {
                assert_eq!(chunks, pieces);
            }
        }
        Ok(())
    })
    .unwrap();
}

// Line DBID N3: the hits of a subject whose first hit follows a hit of another subject that
// shows them (a copy of one record, `QuerySections::joins`) are the other's, after its own;
// the counts of descriptions and alignments count the two as one, and the subject's own HSPs
// are still told apart (`own_hsps`, the heading).
#[test]
fn a_joined_subject_s_hits_are_those_of_the_subject_before_it() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    let group = QueryGroup::new(&by_query[0], &titles, |previous, s_idx| {
        (previous, s_idx) == (5, 2)
    })
    .unwrap();
    assert_eq!(group.subject_order, [5, 7]);
    assert_eq!(group.subject_hit_indices[&5], [0, 1, 2, 3]);
    assert_eq!((group.subject_hits[&5].len(), group.own_hsps(5)), (4, 3));
    assert_eq!(group.aligned(1), [5]);
    assert_eq!(group.described(2), [5, 7]);
    // Only a subject that follows is asked: 7 follows 2, not 5.
    let group = QueryGroup::new(&by_query[0], &titles, |previous, s_idx| {
        (previous, s_idx) == (5, 7)
    })
    .unwrap();
    assert_eq!(group.subject_order, [5, 2, 7]);
    // A chain of three.
    let group = QueryGroup::new(&by_query[0], &titles, |previous, s_idx| {
        matches!((previous, s_idx), (5, 2) | (2, 7))
    })
    .unwrap();
    assert_eq!(group.subject_order, [5]);
    assert_eq!((group.subject_hits[&5].len(), group.own_hsps(5)), (8, 3));
}

// Review 3 of line RG T3a, D1; review 2 of line DBID pfoff-3, M1, and review 1's H1: the
// subjects whose description lines a report shows, which NCBI looks up (by their display ids,
// and where the batch's prefetch is off by their Seq-align ids): the first
// `-num_descriptions` subjects of the group, or in psiblast's later round the subjects of the
// lines of its two-part table, made of its first `-num_descriptions` HSPs (at least one) put
// in the part of the subjects found before and the part of the new ones, that order
// (`two_part_lines`, blast_format.cpp:541-554); a subject shown in another's hit is that
// one's.
#[test]
fn the_described_subjects_are_the_first_subjects_or_those_of_the_first_hsps() {
    let hits = report_hits();
    let described = |q_idx, found_before: Option<&[bool]>, count, joins: fn(u32, u32) -> bool| {
        crate::report::pairwise::blastp_described_subjects(
            Default::default(),
            found_before,
            count,
            q_idx,
            &hits,
            joins,
        )
    };
    let none: fn(u32, u32) -> bool = |_, _| false;
    // Query 0 shows subject 5 with 3 HSPs, 2 with 1 and 7 with 4.
    assert!(described(0, None, 0, none).is_empty());
    assert_eq!(described(0, None, 2, none), [5, 2]);
    assert_eq!(described(0, None, 500, none), [5, 2, 7]);
    // The two-part table counts HSPs: 2 HSPs are both subject 5's, where 2 hits are 5 and 2.
    let new = [false; 8];
    assert_eq!(described(0, Some(&new), 0, none), [5]);
    assert_eq!(described(0, Some(&new), 2, none), [5]);
    assert_eq!(described(0, Some(&new), 3, none), [5]);
    assert_eq!(described(0, Some(&new), 4, none), [5, 2]);
    assert_eq!(described(0, Some(&new), 5, none), [5, 2, 7]);
    // A subject found before comes first, its part before that of the new ones.
    let mut found = new;
    found[2] = true;
    assert_eq!(described(0, Some(&found), 4, none), [2, 5]);
    assert_eq!(described(0, Some(&found), 3, none), [5]);
    found[7] = true;
    assert_eq!(described(0, Some(&found), 500, none), [2, 7, 5]);
    assert_eq!(
        described(0, None, 2, |previous, s_idx| (previous, s_idx) == (5, 2)),
        [5, 7]
    );
    assert!(described(1, None, 500, none).is_empty());
    assert_eq!(described(2, Some(&new), 6, none), [1]);
    assert_eq!(described(2, Some(&new), 7, none), [1, 5]);
}

// Review 2 of line DBID pfoff-3, L5: the model of what a pairwise report reads, made before
// the rows exist (`QueryGroup::read_order`, which each program's `pairwise_read_order` feeds
// with its HSPs), is the writer's own answer on its rows (`writer_reads`: the group, each hit's
// subjects, the description lines, the alignments), whatever the sorts, the counts, the
// copies joined in a hit and psiblast's previous set; then a description line or an alignment
// that the writer reads is never one that the loads leave unread, nor the reverse.
#[test]
fn the_model_of_the_reads_is_the_writers_group() {
    // Subjects 5 (3 HSPs), 2, 7 (2), 3 and its copy 4 (joined after 3): values that each sort
    // orders differently from the rows.
    let values = [
        (5, 40.0, 1e-8, 30, 9),
        (5, 70.0, 1e-12, 10, 1),
        (5, 30.0, 1e-6, 50, 30),
        (2, 90.0, 1e-20, 12, 40),
        (7, 20.0, 1e-3, 80, 3),
        (7, 60.0, 1e-10, 5, 22),
        (3, 50.0, 1e-9, 25, 15),
        (4, 55.0, 1e-9, 25, 2),
    ];
    let rows: Vec<PairwiseHit> = values
        .iter()
        .map(|&(s_idx, bits, e_value, ident, from)| {
            let mut row = hit(0, s_idx, 100);
            row.hit.bit_score = bits;
            row.hit.e_value = e_value;
            row.seqalign_num_ident = ident;
            (row.hit.q_start, row.hit.s_start) = (from, from + 3);
            row
        })
        .collect();
    let joins = |previous: u32, s_idx: u32| (previous, s_idx) == (3, 4);
    let mut found = [false; 8];
    found[7] = true;
    found[4] = true;
    for hits in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
        for hsps in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
            for (descriptions, alignments) in [(0, 0), (1, 0), (2, 1), (3, 5), (500, 2)] {
                for found_before in [None, Some(&found[..])] {
                    for tblastx in [false, true] {
                        let reading = crate::report::pairwise::PairwiseReading {
                            sort: super::super::hit_sort::ReportSort { hits, hsps },
                            shown: crate::report::db_subject::Shown {
                                descriptions,
                                alignments,
                            },
                            html: false,
                        };
                        let model = QueryGroup::read_order(
                            rows.iter().map(|row| (row.hit.s_idx, row)),
                            |row| super::super::hit_sort::HspScores::of(row),
                            joins,
                            reading,
                            tblastx,
                            found_before,
                        );
                        let writer = writer_reads(&rows, joins, reading, tblastx, found_before);
                        assert_eq!(model, writer, "{reading:?} {found_before:?} {tblastx}");
                        // The subjects that the reads load are those the writer shows.
                        let mut shown: Vec<usize> = writer
                            .described
                            .iter()
                            .copied()
                            .chain(0..writer.aligned)
                            .flat_map(|hit| writer.hits[hit].clone())
                            .collect();
                        shown.sort_unstable();
                        shown.dedup();
                        let mut read = model.subjects();
                        read.sort_unstable();
                        assert_eq!(read, shown);
                    }
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:541-554
// Each HSP goes to the part of its subject (found before or new) until `-num_descriptions` HSPs
// (at least one) are put; an HSP of the subject of the part's last line joins that line.
#[test]
fn the_two_part_lines_split_the_first_hsps_by_their_subjects() {
    let hsps = [(5, 5), (5, 6), (2, 2), (5, 5), (7, 7)];
    let found = |own: &u32| *own == 6 || *own == 7;
    let [repeated, new] = two_part_lines(hsps, 4, found);
    assert_eq!(repeated, [(5, vec![6])]);
    assert_eq!(new, [(5, vec![5]), (2, vec![2]), (5, vec![5])]);
    let [repeated, new] = two_part_lines(hsps, 0, found);
    assert!(repeated.is_empty());
    assert_eq!(new, [(5, vec![5])]);
    let [repeated, _] = two_part_lines(hsps, 500, found);
    assert_eq!(repeated, [(5, vec![6]), (7, vec![7])]);
}

#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn windows_of_queries_written_in_turn_write_the_whole_report() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    let marks = Marks {
        num_alignments: 10,
        fail_at: None,
        joined: None,
        merged: false,
    };
    let mut serial = Vec::new();
    write_serial(&mut serial, &marks, 0, &by_query, &titles, None, None).unwrap();
    let estimates: Vec<usize> = by_query
        .iter()
        .map(|query_hits| query_estimate(query_hits, marks.line_length()))
        .collect();
    assert_eq!(
        estimates,
        [0, 1, 2].map(|q_idx| query_estimate_of(
            hits.iter()
                .filter(|hit| hit.hit.q_idx == q_idx)
                .map(|hit| hit.hit.length),
            60
        ))
    );
    crate::utils::threading::with_search_pool(4, "test", |pool| {
        for budget in [1, 3_000, 10_000, usize::MAX] {
            let windows = query_windows(&estimates, budget);
            assert_eq!(windows.first().map(|window| window.start), Some(0));
            assert_eq!(
                windows.last().map(|window| window.end),
                Some(by_query.len())
            );
            assert!(windows.iter().all(|window| !window.is_empty()));
            assert!(windows.windows(2).all(|pair| pair[0].end == pair[1].start));
            if budget == 1 {
                assert_eq!(windows.len(), by_query.len());
            }
            if budget == usize::MAX {
                assert_eq!(windows.len(), 1);
            }
            let (mut serial_windows, mut parallel_windows) = (Vec::new(), Vec::new());
            for window in windows {
                // Only the window's hits, numbered from 0 (a caller that makes the rows of
                // one window at a time).
                let window_hits: Vec<PairwiseHit> = hits
                    .iter()
                    .filter(|hit| window.contains(&(hit.hit.q_idx as usize)))
                    .cloned()
                    .collect();
                let mut window_by_query = vec![Vec::new(); window.len()];
                for (index, hit) in window_hits.iter().enumerate() {
                    window_by_query[hit.hit.q_idx as usize - window.start].push((index, hit));
                }
                write_serial(
                    &mut serial_windows,
                    &marks,
                    window.start,
                    &window_by_query,
                    &titles,
                    None,
                    None,
                )
                .unwrap();
                write_parallel_with(
                    &mut parallel_windows,
                    &marks,
                    window.start,
                    &window_by_query,
                    &titles,
                    None,
                    pool,
                    budget,
                    900,
                )
                .unwrap();
            }
            assert_eq!(serial_windows, serial, "budget {budget}");
            assert_eq!(parallel_windows, serial, "budget {budget}");
        }
        Ok(())
    })
    .unwrap();
}

#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn a_failed_part_writes_its_bytes_before_the_error_as_the_serial_writer() {
    let hits = report_hits();
    let by_query = by_query(&hits);
    let titles = titles();
    // The HSP of length 185 is subject 7's second (query 0).
    let marks = Marks {
        num_alignments: 10,
        fail_at: Some(185),
        joined: None,
        merged: false,
    };
    let mut serial = Vec::new();
    assert!(write_serial(&mut serial, &marks, 0, &by_query, &titles, None, None).is_err());
    crate::utils::threading::with_search_pool(4, "test", |pool| {
        for piece_bytes in [1, 900, usize::MAX] {
            let mut parallel = Vec::new();
            let result = write_parallel_with(
                &mut parallel,
                &marks,
                0,
                &by_query,
                &titles,
                None,
                pool,
                usize::MAX,
                piece_bytes,
            );
            assert_eq!(
                result.unwrap_err().kind(),
                io::ErrorKind::InvalidData,
                "piece_bytes {piece_bytes}"
            );
            assert_eq!(parallel, serial, "piece_bytes {piece_bytes}");
        }
        Ok(())
    })
    .unwrap();
}

// NCBI's default line length (`kDfltLineLength`, 60): the estimate is the one before the line
// length was counted, `400 + 4 * length`, so the default reports' windows, chunks and pieces
// do not move; other line lengths count their lines.
#[test]
fn the_estimate_at_the_default_line_length_is_400_plus_4_per_column() {
    let default = crate::report::pairwise::DEFAULT_LINE_LENGTH;
    for length in (0..20_000).chain([1 << 20, (1 << 28) + 7]) {
        assert_eq!(
            hit_estimate(length, default),
            400 + 4 * length,
            "length {length}"
        );
    }
    // One line of 3 rows per column at `-line_length 1`; the letters alone (and the 400) at
    // the largest line length.
    assert_eq!(hit_estimate(1000, 1), 400 + 63 * 1000);
    assert_eq!(hit_estimate(1000, 7), 400 + 3000 + 60_000 / 7);
    assert_eq!(hit_estimate(1000, 2147483647), 400 + 3000);
    assert_eq!(hit_estimate(1000, 0), hit_estimate(1000, 1));
    assert_eq!(
        query_estimate_of([10, 20], 1),
        QUERY_BYTES + 2 * DESCRIPTION_BYTES + hit_estimate(10, 1) + hit_estimate(20, 1)
    );
    // A query's estimate at 60 is its header's and each hit's, as before.
    let lengths = [7, 59, 60, 61, 1000, 12_345];
    assert_eq!(
        query_estimate_of(lengths, default),
        lengths.iter().fold(QUERY_BYTES, |bytes, &length| bytes
            + DESCRIPTION_BYTES
            + 400
            + 4 * length)
    );
    // Elsewhere the lines are counted over the query's columns in all.
    assert_eq!(
        query_estimate_of([10, 11], 7),
        QUERY_BYTES + 2 * DESCRIPTION_BYTES + 800 + 3 * 21 + 60 * 21 / 7
    );
}

// DBID pfoff-3 (DECISIONS 62 (13)): the writers end the run before the load of the aligned
// subject whose lookup throws where the report's writer ends there (`aligned_end`): outfmt 0
// at its first alignment only, outfmt 1 to 4 at any; the error carries NCBI's exception
// (`NativeError::of_writer`: the application's report of the database error, exit status 2).
#[test]
fn the_writers_end_the_run_at_an_aligned_subject_whose_load_throws() {
    use crate::report::db_subject::AfterThrow;
    let exception =
        crate::report::db_subject::LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO)
            .exception()
            .unwrap();
    // 7 throws; 9 throws with a local Seq-align id.
    let thrown = |s_idx: u32| match s_idx {
        7 => Some((&exception, false)),
        9 => Some((&exception, true)),
        _ => None,
    };
    let end = |after, aligned: &[u32]| aligned_end(after, aligned, thrown).map(|(at, _)| at);
    assert_eq!(end(AfterThrow::SkipsAlignment, &[7, 1]), Some(0));
    assert_eq!(end(AfterThrow::SkipsAlignment, &[1, 7]), None);
    assert_eq!(end(AfterThrow::SkipsAlignment, &[9, 7]), None);
    assert_eq!(end(AfterThrow::EndsAtAlignment, &[1, 9, 7]), Some(1));
    let at = aligned_end(AfterThrow::EndsAtAlignment, &[1, 7], thrown);
    assert!(end_here(at, 0).is_ok());
    let error = end_here(at, 1).unwrap_err();
    let carried = crate::ncbi_exception::NcbiException::of_io_error(&error).unwrap();
    assert_eq!(
        carried.get_msg(),
        crate::seqdb::tax_files::NO_ACCESSION_INFO.as_bytes()
    );
}

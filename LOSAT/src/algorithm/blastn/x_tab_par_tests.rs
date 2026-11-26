use super::super::hsp::{BlastnHsp, BlastnHspList};
use super::*;
use std::io::Write;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

fn hsp(rng: &mut Rng, q_idx: u32, s_idx: u32) -> BlastnHsp {
    let length = 1 + rng.below(3000) as usize;
    let num_ident = rng.below(length as u64 + 1) as usize;
    let e_value = match rng.below(5) {
        0 => 0.0,
        1 => 1e-180 * (1 + rng.below(9)) as f64,
        2 => (rng.below(1000) as f64) * 1e-5,
        3 => 1.0 + rng.below(9) as f64 / 7.0,
        _ => 10f64.powi(-(rng.below(300) as i32)),
    };
    let q_start = 1 + rng.below(5000) as usize;
    let s_start = 1 + rng.below(3_000_000) as usize;
    let minus = rng.below(2) == 1;
    BlastnHsp {
        identity: 100.0 * num_ident as f64 / length as f64,
        length,
        mismatch: length - num_ident,
        gapopen: rng.below(9) as usize,
        q_start,
        q_end: q_start + length - 1,
        s_start,
        s_end: if minus {
            s_start.saturating_sub(length - 1).max(1)
        } else {
            s_start + length - 1
        },
        e_value,
        bit_score: rng.below(6000) as f64 / 3.0 + 0.5,
        num_ident,
        query_frame: if minus { -1 } else { 1 },
        query_length: q_start + length,
        q_idx,
        s_idx,
        raw_score: rng.below(3000) as i32,
        internal_q_offset_0: (q_start - 1) as i32,
        internal_q_end_0: (q_start + length - 1) as i32,
        internal_s_offset_0: s_start - 1,
        internal_s_end_0: s_start + length - 1,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: num_ident,
    }
}

/// Hit lists of `queries` queries: some without a list, some with empty HSP lists, the others
/// with up to `max_rows` rows over up to three subjects.
fn hit_lists(rng: &mut Rng, queries: usize, max_rows: u64) -> Vec<Option<BlastnHitList>> {
    (0..queries as u32)
        .map(|q_idx| {
            if rng.below(4) == 0 {
                return None;
            }
            let subjects = rng.below(4) as u32;
            let mut list = BlastnHitList::new(subjects.max(1) as usize);
            for s_idx in 0..subjects {
                let rows = rng.below(max_rows + 1);
                let hsps = (0..rows).map(|_| hsp(rng, q_idx, s_idx)).collect();
                list.hsplist_array.push(BlastnHspList {
                    oid: s_idx,
                    query_index: q_idx,
                    hsps,
                    best_evalue: 0.0,
                });
            }
            list.hsplist_count = list.hsplist_array.len();
            Some(list)
        })
        .collect()
}

fn ids(prefix: &str, n: usize) -> Vec<Arc<[u8]>> {
    (0..n)
        .map(|i| Arc::<[u8]>::from(format!("{prefix}{i}").into_bytes()))
        .collect()
}

struct Case {
    hit_lists: Vec<Option<BlastnHitList>>,
    query_ids: Vec<Arc<[u8]>>,
    subject_ids: Vec<crate::algorithm::common::tabular::SubjectId>,
    query_titles: Vec<Arc<[u8]>>,
    unsearched: Vec<bool>,
    warnings: Vec<crate::report::diag_log::Posts>,
}

fn case(rng: &mut Rng, queries: usize, max_rows: u64) -> Case {
    Case {
        hit_lists: hit_lists(rng, queries, max_rows),
        query_ids: ids("query_", queries),
        subject_ids: ids("lcl|subject_", 3).into_iter().map(Some).collect(),
        query_titles: ids("query title ", queries),
        unsearched: (0..queries).map(|_| rng.below(9) == 0).collect(),
        warnings: (0..queries)
            .map(|q| {
                if rng.below(5) == 0 {
                    crate::report::diag_log::Posts::of_lines(
                        format!("Warning: [blastn] Query_{q} warning\n").into_bytes(),
                    )
                } else {
                    crate::report::diag_log::Posts::new()
                }
            })
            .collect(),
    }
}

/// What the report's writer and the warnings' sink saw, in order: the report's bytes (writes
/// one after another joined), the report's flushes and the sink's bytes.
#[derive(Debug, PartialEq)]
enum Event {
    Report(Vec<u8>),
    Flush,
    Sink(Vec<u8>),
}

/// The report's writer (`sink` false) or the warnings' sink, both adding to one log.
struct Side {
    log: Arc<std::sync::Mutex<Vec<Event>>>,
    sink: bool,
}

impl crate::report::ostream::ReportWrite for Side {}

impl Write for Side {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        let mut events = self.log.lock().unwrap();
        match (self.sink, events.last_mut()) {
            (false, Some(Event::Report(bytes))) | (true, Some(Event::Sink(bytes))) => {
                bytes.extend_from_slice(buf)
            }
            (false, _) => events.push(Event::Report(buf.to_vec())),
            (true, _) => events.push(Event::Sink(buf.to_vec())),
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        if !self.sink {
            self.log.lock().unwrap().push(Event::Flush);
        }
        Ok(())
    }
}

/// The log of the serial writer (`pool` `None`) or of the pool's: the report, its flushes and
/// the warnings in the order they were written.
fn run(
    case: &Case,
    format: BlastnOutputFormat,
    epilog: Option<usize>,
    pool: Option<&SearchPool<'_>>,
) -> Vec<Event> {
    let log = Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut report = Side {
        log: log.clone(),
        sink: false,
    };
    let mut sink = Side {
        log: log.clone(),
        sink: true,
    };
    let mut warnings = QueryWarnings::of_posts(&case.warnings, None, &mut sink);
    let args = (
        &case.hit_lists[..],
        &case.query_ids[..],
        &case.subject_ids[..],
        &case.query_titles[..],
        &case.unsearched[..],
    );
    match pool {
        Some(pool) => write_on_pool(
            pool,
            args.0,
            &mut report,
            args.1,
            args.2,
            format,
            args.3,
            "User specified sequence set (Input: subject.fa)",
            args.4,
            epilog,
            Some(&mut warnings),
        ),
        None => write_output_blastn_hitlists_to_writer(
            args.0,
            &mut report,
            args.1,
            args.2,
            &crate::report::db_subject::HitLoads::default(),
            format,
            args.3,
            "User specified sequence set (Input: subject.fa)",
            args.4,
            epilog,
            None,
            Some(&mut warnings),
            None,
        )
        .map(|stop| assert_eq!(stop, None)),
    }
    .unwrap();
    // The warnings end (`QueryWarnings`'s drop writes any read's text not placed) before the
    // sink is read.
    drop(warnings);
    drop((report, sink));
    Arc::try_unwrap(log).unwrap().into_inner().unwrap()
}

/// The report's bytes of a log.
fn report_bytes(events: &[Event]) -> Vec<u8> {
    events
        .iter()
        .filter_map(|event| match event {
            Event::Report(bytes) => Some(&bytes[..]),
            _ => None,
        })
        .flatten()
        .copied()
        .collect()
}

#[test]
fn rows_on_the_pool_are_the_serial_writers_bytes() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let cases: Vec<Case> = [
        (0usize, 0u64),
        (1, 0),
        (1, 5),
        (7, 3),
        (60, 40),
        (400, 9),
        (25, 3000),
    ]
    .into_iter()
    .map(|(queries, max_rows)| case(&mut rng, queries, max_rows))
    .collect();
    // Warnings between the queries' rows (the order of report, flush and sink is compared).
    let mut sinks = 0;
    for threads in [2usize, 3, 8] {
        crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
            for case in &cases {
                let queries = case.hit_lists.len();
                for format in [
                    BlastnOutputFormat::Tabular,
                    BlastnOutputFormat::TabularWithComments,
                ] {
                    for epilog in [None, Some(queries + 3)] {
                        let serial = run(case, format, epilog, None);
                        let pooled = run(case, format, epilog, Some(pool));
                        assert_eq!(
                            serial, pooled,
                            "threads {threads} queries {queries} {format:?} {epilog:?}"
                        );
                        sinks += serial
                            .iter()
                            .filter(|event| matches!(event, Event::Sink(_)))
                            .count();
                    }
                }
            }
            Ok(())
        })
        .unwrap();
    }
    assert!(sinks > 100, "{sinks} warnings");
}

#[test]
fn a_query_larger_than_a_window_and_many_small_ones() {
    // The first query has more rows than a window of two threads, the others a few each.
    let mut rng = Rng(7);
    let mut case = case(&mut rng, 50, 4);
    let mut big = BlastnHitList::new(1);
    big.hsplist_array.push(BlastnHspList {
        oid: 1,
        query_index: 0,
        hsps: (0..(ROWS_PER_THREAD * 2 + 17))
            .map(|_| hsp(&mut rng, 0, 1))
            .collect(),
        best_evalue: 0.0,
    });
    big.hsplist_count = 1;
    case.hit_lists[0] = Some(big);
    for format in [
        BlastnOutputFormat::Tabular,
        BlastnOutputFormat::TabularWithComments,
    ] {
        let serial = run(&case, format, None, None);
        crate::utils::threading::with_search_pool(2, "blastn", |pool| {
            assert_eq!(serial, run(&case, format, None, Some(pool)), "{format:?}");
            Ok(())
        })
        .unwrap();
        assert!(report_bytes(&serial).len() > ROWS_PER_THREAD * 2 * 40);
    }
    let first = pieces(&case, 0)
        .take_while(|piece| piece.q_idx == 0)
        .count();
    assert_eq!(first, (ROWS_PER_THREAD * 2 + 17).div_ceil(PIECE_ROWS));
}

fn pieces(case: &Case, header_bytes: usize) -> Pieces<'_> {
    Pieces {
        hit_lists: &case.hit_lists,
        query_ids: &case.query_ids,
        subject_ids: &case.subject_ids,
        header_bytes,
        q_idx: 0,
        list: 0,
        hsp: 0,
        header: true,
    }
}

#[test]
fn pieces_take_every_row_once_in_order() {
    let mut rng = Rng(0x5eed);
    for (queries, max_rows) in [(0usize, 0u64), (30, 0), (40, 7), (12, 1500), (3, 5000)] {
        let case = case(&mut rng, queries, max_rows);
        let mut seen = vec![0usize; queries];
        let mut last = None;
        for piece in pieces(&case, 0) {
            assert!(piece.rows <= PIECE_ROWS);
            // A query's first piece has its header; the queries come in order, none skipped.
            assert_eq!(
                piece.header,
                last != Some(piece.q_idx),
                "{queries} {max_rows}"
            );
            if piece.header {
                assert_eq!(piece.q_idx, last.map_or(0, |q| q + 1));
                assert_eq!((piece.list, piece.hsp), (0, 0));
            }
            // The piece starts where the query's rows written so far end.
            let lists = case.hit_lists[piece.q_idx]
                .as_ref()
                .map_or(&[][..], |list| &list.hsplist_array[..]);
            let before: usize = lists[..piece.list]
                .iter()
                .map(|list| list.hsps.len())
                .sum::<usize>()
                + piece.hsp;
            assert_eq!(before, seen[piece.q_idx]);
            seen[piece.q_idx] += piece.rows;
            last = Some(piece.q_idx);
        }
        assert_eq!(last.map_or(0, |q| q + 1), queries);
        for (q_idx, rows) in seen.into_iter().enumerate() {
            let all: usize = case.hit_lists[q_idx].as_ref().map_or(0, |list| {
                list.hsplist_array.iter().map(|hsps| hsps.hsps.len()).sum()
            });
            assert_eq!(rows, all, "query {q_idx}");
        }
    }
}

/// Counts the rows a format observer is told of.
struct Rows(usize);

impl crate::api::local_blast::FormatObserver for Rows {
    fn hsp_begin(&mut self, _format: usize, _hsp: crate::api::local_blast::HspIndex) {
        self.0 += 1;
    }
    fn hsp_end(&mut self, _format: usize, _hsp: crate::api::local_blast::HspIndex) {}
}

#[test]
fn a_format_observer_keeps_the_serial_writer() {
    let mut rng = Rng(11);
    let case = case(&mut rng, 30, 20);
    let rows: usize = case
        .hit_lists
        .iter()
        .flatten()
        .flat_map(|list| &list.hsplist_array)
        .map(|hsps| hsps.hsps.len())
        .sum();
    assert!(rows > 0);
    let report = |pool: Option<&SearchPool<'_>>, observer: Option<&mut Rows>| {
        let mut report = Vec::new();
        let mut probe = observer.map(|rows| FormatProbe::new(rows, 0));
        write_output(
            pool,
            &case.hit_lists,
            &mut report,
            &case.query_ids,
            &case.subject_ids,
            &crate::report::db_subject::HitLoads::default(),
            BlastnOutputFormat::Tabular,
            &case.query_titles,
            "subject",
            &case.unsearched,
            None,
            probe.as_mut(),
            None,
            None,
        )
        .unwrap();
        report
    };
    let serial = report(None, None);
    crate::utils::threading::with_search_pool(4, "blastn", |pool| {
        let mut observer = Rows(0);
        assert_eq!(serial, report(Some(pool), Some(&mut observer)));
        // The observer saw every row: the serial writer wrote them, not the pool's path.
        assert_eq!(observer.0, rows);
        Ok(())
    })
    .unwrap();
}

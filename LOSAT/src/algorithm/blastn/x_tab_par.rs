//! EXPERIMENT (LOSAT_X_BNTABPAR): BLASTN's outfmt 6 / 7 rows of a report formatted on the search's
//! pool.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:828-832
//! ```c
//! ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
//!         const CSeq_align& s = **itr;
//!         tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
//!         tabinfo.Print();
//! }
//! ```
//! NCBI prints the rows of a query's alignments one after another, and the queries of a batch one
//! after another, on the thread of the application's loop. `write_output_blastn_hitlists_to_writer`
//! (`hsp.rs`) does the same on the thread of `run_in_pool`'s batch loop, while the other threads of
//! the search's pool wait: 31,168 cDNAs against 3 Mb with `-task blastn` (1.02 M rows, 90 query
//! batches written one at a time by LOSAT_X_RPTSTREAM, LOSAT_X_QBATCHPAR=0) spent 0.88 s in the
//! report calls at 8, 16 and 32 threads alike (ws, 2026-10-24), 87 % of it in the rows
//! (`write_hit_fields_bytes`, local perf).
//!
//! What the switch changes (read once per process, `mode`), when the caller gives a pool of two
//! threads or more (`pool`: the search's pool while no LOSAT_X_QBATCHPAR slot exists), there is no
//! format observer (`FormatProbe`, which flushes around every row; the Web adapter) and the
//! BLASTN trace is off (`tracing::enabled`; it logs rows in the serial order): the rows are cut
//! into pieces (`Piece`) of at most `PIECE_ROWS` rows of one query, the first piece of each query
//! with its outfmt 7 header (a query without rows is one piece of no row); the pieces are taken in
//! windows of at least `ROWS_PER_THREAD` rows per thread of the pool (a piece of no row counts as
//! one), and the pool formats each piece of a window into a buffer of its own with
//! `write_query_header` and `write_hsp_rows`, the two parts of the serial function of one query
//! (`write_query_rows`); the calling thread then writes the window's buffers in order, each
//! query's warnings (`QueryWarnings`) before its first piece, before the next window is formatted.
//! The last line (`write_tabular_epilog`) stays on the calling thread.
//!
//! Why the bytes are those of the serial writer: the pieces of a query, in order, are its header
//! and then its HSP lists' rows in order, which is what `write_query_rows` writes for it, into
//! buffers instead of the writer; the serial writer is that call for each query in order with the
//! query's warnings before it (`QueryWarnings::before_query`, which flushes the report written so
//! far and then writes to its own sink, as it does here). A row reads only its HSP, the ids and
//! the fixed configuration, and the header only its query's title, hit list and the subject
//! title; no value is carried from one row or query to the next (the running HSP index is read
//! only by the format observer, which keeps the serial writer). Writing to a buffer cannot fail,
//! so a failed write is one of the writer's, which stops the report as before.
//!
//! Off (the serial writer): the switch off (`LOSAT_X_BNTABPAR=0`); no pool or one of one thread
//! (`-num_threads 1`, builds without the `parallel` feature, Wasm without threads); while
//! LOSAT_X_QBATCHPAR slots exist (their pools would share the CPUs; their reports already overlap
//! the other slots' searches); a format observer; the BLASTN trace.
//!
//! Memory: the formatted rows of one window (`ROWS_PER_THREAD` x threads rows plus one piece,
//! about 6 MiB at 32 threads with the ~100-byte rows of the cDNA input), also for one query of
//! many rows or many queries without rows.
//!
//! With LOSAT_X_REACH the decision is recorded once per report call of at least one query
//! (`x_reach!`; the call of the epilog alone formats no row).
//!
//! No NCBI counterpart for the windows and buffers: they choose where the same bytes are made and
//! change no value NCBI computes.

use std::io::{self, Write};
use std::sync::{Arc, OnceLock};

use crate::api::local_blast::FormatProbe;
use crate::report::query_warnings::QueryWarnings;
use crate::report::OutputConfig;
use crate::utils::threading::SearchPool;

use super::hsp::{
    write_hsp_rows, write_output_blastn_hitlists_to_writer, write_query_header,
    write_tabular_epilog, BlastnHitList, BlastnOutputFormat,
};

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// Rows per thread of the pool in one window of pieces.
const ROWS_PER_THREAD: usize = 2048;

/// Rows of one piece at most (four pieces per thread in a window of one large query).
const PIECE_ROWS: usize = ROWS_PER_THREAD / 4;

/// Bytes of a row besides its two ids (eleven tabs, the newline and ten numbers; a tabular row of
/// the cDNA input is ~100 bytes).
const ROW_BYTES: usize = 72;

/// Bytes of outfmt 7's header lines besides the two titles.
const HEADER_BYTES: usize = 320;

// No NCBI counterpart: reads the switch once (LOSAT_X_BNTABPAR is a strict-set switch: on unless
// it is `0`). It chooses between two LOSAT code paths with the same result.
fn mode() -> bool {
    static MODE: OnceLock<bool> = OnceLock::new();
    *MODE.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNTABPAR"))
}

/// The pool the tabular rows are formatted on (`Some` with the switch on and a pool of two
/// threads or more), or `None` for the serial writer.
pub(crate) fn pool<'p, 'q>(pool: &'p SearchPool<'q>) -> Option<&'p SearchPool<'q>> {
    const PARALLEL: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));
    (PARALLEL && mode() && pool.enabled() && pool.threads() > 1).then_some(pool)
}

/// `write_output_blastn_hitlists_to_writer`'s report: its rows formatted on `pool` (from `pool`)
/// where the switch's rule allows, else by the serial writer.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_output<W: Write>(
    pool: Option<&SearchPool<'_>>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: Option<usize>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    match pool {
        Some(pool) if probe.is_none() && !super::tracing::enabled() => {
            if !hit_lists.is_empty() {
                crate::utils::xreach::x_reach!("LOSAT_X_BNTABPAR", true);
            }
            write_on_pool(
                pool,
                hit_lists,
                writer,
                query_ids,
                subject_ids,
                output_format,
                query_titles,
                subject_title,
                unsearched,
                epilog,
                warnings,
            )
        }
        _ => {
            if mode() && !hit_lists.is_empty() {
                crate::utils::xreach::x_reach!("LOSAT_X_BNTABPAR", false);
            }
            write_output_blastn_hitlists_to_writer(
                hit_lists,
                writer,
                query_ids,
                subject_ids,
                output_format,
                query_titles,
                subject_title,
                unsearched,
                epilog,
                probe,
                warnings,
            )
        }
    }
}

/// Rows of one query, from HSP `hsp` of HSP list `list` of its hit list on: the first piece of the
/// query (`header`) begins with its outfmt 7 header.
struct Piece {
    q_idx: usize,
    header: bool,
    list: usize,
    hsp: usize,
    rows: usize,
    /// The estimated bytes of the piece (the capacity of its buffer).
    bytes: usize,
}

/// The pieces of a report's queries in order (`Piece`), from query 0 on.
struct Pieces<'a> {
    hit_lists: &'a [Option<BlastnHitList>],
    query_ids: &'a [Arc<[u8]>],
    subject_ids: &'a [Arc<[u8]>],
    /// Header bytes of a query's first piece (0 for outfmt 6).
    header_bytes: usize,
    q_idx: usize,
    list: usize,
    hsp: usize,
    header: bool,
}

impl Iterator for Pieces<'_> {
    type Item = Piece;

    fn next(&mut self) -> Option<Piece> {
        let q_idx = self.q_idx;
        let hit_list = self.hit_lists.get(q_idx)?;
        let lists = hit_list
            .as_ref()
            .map_or(&[][..], |list| &list.hsplist_array[..]);
        let id_bytes = |s_idx: u32| {
            self.subject_ids
                .get(s_idx as usize)
                .map_or(7, |id| id.len())
        };
        let query_bytes = self.query_ids.get(q_idx).map_or(7, |id| id.len()) + ROW_BYTES;
        let (list, hsp, header) = (self.list, self.hsp, self.header);
        let mut rows = 0;
        let mut bytes = if header { self.header_bytes } else { 0 };
        while self.list < lists.len() && rows < PIECE_ROWS {
            let hsps = &lists[self.list].hsps;
            let take = (hsps.len() - self.hsp).min(PIECE_ROWS - rows);
            if take > 0 {
                bytes += take * (query_bytes + id_bytes(hsps[self.hsp].s_idx));
            }
            rows += take;
            self.hsp += take;
            if self.hsp == hsps.len() {
                self.list += 1;
                self.hsp = 0;
            }
        }
        while self.list < lists.len() && lists[self.list].hsps.is_empty() {
            self.list += 1;
        }
        if self.list == lists.len() {
            (self.q_idx, self.list, self.hsp, self.header) = (q_idx + 1, 0, 0, true);
        } else {
            self.header = false;
        }
        Some(Piece {
            q_idx,
            header,
            list,
            hsp,
            rows,
            bytes,
        })
    }
}

/// The serial writer's bytes, with each window's pieces formatted on `pool` (module docs).
#[allow(clippy::too_many_arguments)]
fn write_on_pool<W: Write>(
    pool: &SearchPool<'_>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: Option<usize>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = OutputConfig::ncbi_compat();
    let window_rows = ROWS_PER_THREAD.saturating_mul(pool.threads().max(1));
    let mut pieces = Pieces {
        hit_lists,
        query_ids,
        subject_ids,
        header_bytes: match output_format {
            BlastnOutputFormat::TabularWithComments => HEADER_BYTES + subject_title.len(),
            _ => 0,
        },
        q_idx: 0,
        list: 0,
        hsp: 0,
        header: true,
    };
    let format = |piece: &Piece| -> io::Result<Vec<u8>> {
        let hit_list = hit_lists[piece.q_idx].as_ref();
        let title = query_titles.get(piece.q_idx).map_or(0, |title| title.len());
        let mut buffer = Vec::with_capacity(piece.bytes + if piece.header { title } else { 0 });
        if piece.header {
            write_query_header(
                &mut buffer,
                piece.q_idx,
                hit_list,
                output_format,
                query_titles,
                subject_title,
                unsearched,
            )?;
        }
        // The running HSP index is the format observer's, which this path never has.
        let mut hsp_index = 0;
        let lists = hit_list.map_or(&[][..], |list| &list.hsplist_array[..]);
        let (mut list, mut hsp, mut rows) = (piece.list, piece.hsp, piece.rows);
        while rows > 0 {
            let hsps = &lists[list].hsps;
            let take = (hsps.len() - hsp).min(rows);
            write_hsp_rows(
                &mut buffer,
                &hsps[hsp..hsp + take],
                query_ids,
                subject_ids,
                &config,
                &mut hsp_index,
                None,
            )?;
            rows -= take;
            list += 1;
            hsp = 0;
        }
        Ok(buffer)
    };
    loop {
        let mut window = Vec::new();
        let mut rows = 0;
        while rows < window_rows {
            match pieces.next() {
                Some(piece) => {
                    rows += piece.rows.max(1);
                    window.push(piece);
                }
                None => break,
            }
        }
        if window.is_empty() {
            break;
        }
        let buffers = map(pool, &window, &format);
        for (piece, buffer) in window.iter().zip(buffers) {
            if piece.header {
                if let Some(warnings) = warnings.as_deref_mut() {
                    warnings.before_query(piece.q_idx, &mut *writer)?;
                }
            }
            writer.write_all(&buffer?)?;
        }
    }
    write_tabular_epilog(writer, output_format, epilog)
}

/// `f` of every piece of `pieces`, in order, on `pool`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn map<R: Send>(
    pool: &SearchPool<'_>,
    pieces: &[Piece],
    f: &(impl Fn(&Piece) -> R + Sync),
) -> Vec<R> {
    pool.install(|| pieces.par_iter().map(f).collect())
}

/// `f` of every piece of `pieces`, in order (builds without threads, where `pool` gives `None`
/// and this is never called).
#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn map<R: Send>(
    _pool: &SearchPool<'_>,
    pieces: &[Piece],
    f: &(impl Fn(&Piece) -> R + Sync),
) -> Vec<R> {
    pieces.iter().map(f).collect()
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod tests {
    use super::super::hsp::{BlastnHsp, BlastnHspList};
    use super::*;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
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
            internal_q_offset_0: q_start - 1,
            internal_q_end_0: q_start + length - 1,
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
        subject_ids: Vec<Arc<[u8]>>,
        query_titles: Vec<Arc<[u8]>>,
        unsearched: Vec<bool>,
        warnings: Vec<Vec<u8>>,
    }

    fn case(rng: &mut Rng, queries: usize, max_rows: u64) -> Case {
        Case {
            hit_lists: hit_lists(rng, queries, max_rows),
            query_ids: ids("query_", queries),
            subject_ids: ids("lcl|subject_", 3),
            query_titles: ids("query title ", queries),
            unsearched: (0..queries).map(|_| rng.below(9) == 0).collect(),
            warnings: (0..queries)
                .map(|q| {
                    if rng.below(5) == 0 {
                        format!("Warning: [blastn] Query_{q} warning\n").into_bytes()
                    } else {
                        Vec::new()
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
        let mut warnings = QueryWarnings {
            before: &case.warnings,
            sink: &mut sink,
        };
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
                format,
                args.3,
                "User specified sequence set (Input: subject.fa)",
                args.4,
                epilog,
                None,
                Some(&mut warnings),
            ),
        }
        .unwrap();
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
                BlastnOutputFormat::Tabular,
                &case.query_titles,
                "subject",
                &case.unsearched,
                None,
                probe.as_mut(),
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
}

//! EXPERIMENT (LOSAT_X_BNTABPAR): BLASTN's outfmt 6 / 7 rows of a report formatted on the search's
//! pool.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:828-832
//! NCBI prints the rows of a query's alignments one after another, and the queries of a batch one
//! after another, on the thread of the application's loop. `write_output_blastn_hitlists_to_writer`
//! (`hsp.rs`) does the same on the thread of `run_in_pool`'s batch loop, while the other threads of
//! the search's pool wait: 31,168 cDNAs against 3 Mb with `-task blastn` (1.02 M rows, 90 query
//! batches written one at a time by LOSAT_X_RPTSTREAM, no query-batch slots) spent 0.88 s in the
//! report calls at 8, 16 and 32 threads alike (ws, 2026-10-24), 87 % of it in the rows
//! (`write_hit_fields_bytes`, local perf).
//!
//! What the switch changes (read once per process, `mode`), when the caller gives a pool of two
//! threads or more (`pool`: the search's pool), there is no format observer (`FormatProbe`, which
//! flushes around every row; the Web adapter), the BLASTN trace is off (`tracing::enabled`; it logs rows in the serial order) and no subject's
//! rows are written as a failed lookup leaves them or stop (a load that throws, an id that fails:
//! `write_query_rows`, which carries the result set's rows from one to the next): the rows are cut
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
//! (`-num_threads 1`, builds without the `parallel` feature, Wasm without threads); a format
//! observer; the BLASTN trace; a subject whose load throws or whose id fails.
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

use std::io;
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
pub(crate) fn write_output<W: crate::report::ostream::ReportWrite>(
    pool: Option<&SearchPool<'_>>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    loads: &crate::report::db_subject::HitLoads,
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: Option<usize>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut QueryWarnings<'_>>,
    left_of_record: Option<(&crate::blastinput::seq_range::Placements, bool)>,
) -> io::Result<Option<crate::report::tabular::RowStop>> {
    match pool {
        // The serial writer writes the rows of subjects whose load throws or whose id fails, the
        // rows that stop (`write_output_blastn_hitlists_to_writer`), and reports with a `-db` HSP
        // starting left of its query record.
        Some(pool)
            if probe.is_none()
                && !super::tracing::enabled()
                && !loads.any_throws()
                && !subject_ids.iter().any(Option::is_none)
                && left_of_record.is_none() =>
        {
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
            .map(|()| None)
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
                loads,
                output_format,
                query_titles,
                subject_title,
                unsearched,
                epilog,
                probe,
                warnings,
                left_of_record,
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
    subject_ids: &'a [crate::algorithm::common::tabular::SubjectId],
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
                .map_or(7, |id| id.as_deref().map_or(0, <[u8]>::len))
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
fn write_on_pool<W: crate::report::ostream::ReportWrite>(
    pool: &SearchPool<'_>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
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
    // The queries whose reports read subjects that write texts: their pieces mark where the
    // serial writer reads each subject (its first row), and the texts are written there
    // (`QueryWarnings::subject_read`).
    let marked: Vec<bool> = (0..hit_lists.len())
        .map(|q_idx| {
            warnings
                .as_deref()
                .is_some_and(|warnings| warnings.reads_with_texts(q_idx))
        })
        .collect();
    let format = |piece: &Piece| -> io::Result<(Vec<u8>, Vec<(usize, u32)>)> {
        let mut reads = Vec::new();
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
            if hsp == 0 && take > 0 && marked[piece.q_idx] {
                reads.push((buffer.len(), lists[list].oid));
            }
            write_hsp_rows(
                &mut buffer,
                &hsps[hsp..hsp + take],
                query_ids,
                subject_ids,
                &config,
                &mut hsp_index,
                None,
                None,
            )?;
            rows -= take;
            list += 1;
            hsp = 0;
        }
        Ok((buffer, reads))
    };
    // Whether a query's report was written (its `CBlastTabularInfo` lives).
    let mut started = false;
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
                // `~CBlastTabularInfo` of the query before flushes the stream.
                if std::mem::replace(&mut started, true) {
                    writer.ncbi_flush()?;
                }
                if let Some(warnings) = warnings.as_deref_mut() {
                    warnings.before_query(piece.q_idx, &mut *writer)?;
                }
            }
            let (buffer, reads) = buffer?;
            crate::report::query_warnings::write_marked(
                &buffer,
                &reads,
                piece.q_idx,
                &mut *writer,
                warnings.as_deref_mut(),
            )?;
        }
    }
    if started {
        writer.ncbi_flush()?;
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
#[path = "x_tab_par_tests.rs"]
mod tests;

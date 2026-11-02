//! EXPERIMENT (LOSAT_X_PAIRWISEPAR): TBLASTX's report stage
//! (`write_tblastx_outputs` in blast_engine/run_impl.rs) on the search's pool (until 2026-11-02
//! under a switch of its own, LOSAT_X_TBXRPTPAR, folded into LOSAT_X_PAIRWISEPAR,
//! `report/x_pairwise_par.rs`).
//!
//! The reports are written once, after the last query batch, by the thread of the batch loop.
//! Before this switch every part of that stage ran on that thread while the other threads of
//! the pool waited: on the genome 3 x 3 input at 16 threads (ws, 2026-10-19) the stage took
//! 1.31 s of 8.6 s with outfmt 6 and 3.37 s of 10.6 s with outfmt 0, of which the displayed
//! identities took 0.99 s (each HSP's codons translated again from the nucleotides), the
//! outfmt 0 rows 1.80 s (translated a second time, with the masks lowercased), the tabular rows
//! 0.17 s (outfmt 0's text 0.41 s) and the final order (then `report_order` over every batch's
//! HSPs at once; now `hit_collector.rs` inside each batch's search) 0.15 s.
//!
//! What the switch changes (`x_pairwise_par::pool`: on, and a pool of two threads or more,
//! `pool`):
//! 1. The displayed identities (`report::set_displayed_identities`) run over chunks of the final
//!    list on the pool, each chunk with the serial function.
//! 2. The outfmt 0 rows (`report::pairwise_hits`, also the Web adapter's hit records) are an
//!    indexed parallel map of `PairwiseRows::hit` over the final list, or over the HSPs of one
//!    window of queries (`run_impl::write_tblastx_outputs`), collected in its order.
//! 3. The outfmt 6 / 7 rows (`report::write_tabular`) take BLASTN's pieces and windows
//!    (`algorithm/blastn/x_tab_par.rs`): pieces of at most `PIECE_ROWS` rows of one query, the
//!    first piece of each query with its outfmt 7 header (a query without rows is one piece of
//!    no row), in windows of at least `ROWS_PER_THREAD` rows per thread (a piece of no row counts
//!    as one); the pool formats each piece of a window into a buffer of its own, and the loop's
//!    thread writes the window's buffers in order, each query's warnings before its first piece,
//!    before the next window is formatted (`x_pairwise_par::write_in_order`), so at most one
//!    window of rows is held unwritten. The last line stays on the loop's thread. With the Web
//!    adapter's format observer (`FormatProbe`, which flushes around every row) the rows are
//!    written by the serial function. Until 2026-11-02 only a query of more than 1,024 rows had
//!    its rows formatted on the pool, in windows of its own (4 chunks of 1,024 rows per thread).
//! With a Seq-align format requested (outfmt 8, 12), the engine counts of the Seq-aligns
//! (`report::engine_counts`) run over chunks of the final list on the pool as item 1 does.
//! The outfmt 0 text (`write_tblastx_pairwise_report`) and the writes stay on the loop's thread
//! (the final order is made in each batch's search, `hit_collector.rs`).
//!
//! Why the bytes are those of the serial stage: every value is computed by the serial code's
//! function of one HSP (or of one row) on the same inputs (the HSP, the query and subject
//! letters, the genetic codes, the placements, the masks), which reads nothing another HSP
//! writes; the order of the final list is not changed (chunks are disjoint slices; the map is
//! collected by index); a query's pieces in order are its outfmt 7 header and its rows in order,
//! which is what the serial function writes for it after the query's warnings.
//! No floating-point value is summed across HSPs. The same bytes go to the writer in the same
//! order, in larger writes; a failed write stops the run as before.
//!
//! Off (the serial stage): the switch off (`LOSAT_X_PAIRWISEPAR=0`);
//! a pool of one thread or none (`-num_threads 1`, builds without the `parallel`
//! feature, Wasm without threads).
//!
//! Memory: the outfmt 0 rows are one list of the final list's length in both forms; the tabular
//! windows hold at most `ROWS_PER_THREAD` x threads rows formatted plus one piece (~3.3 MB at 16
//! threads for rows of ~100 bytes).

use std::io::{self, Write};
use std::sync::Arc;

use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::pairwise::PairwiseHit;
use crate::report::x_pairwise_par::{self, PIECE_ROWS};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::threading::SearchPool;

use super::blast_engine::TblastxHsp;
use super::report::{self, PairwiseRows, SubjectLetters};
use super::scoring::TbxScore;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// The HSPs of one task (the displayed identities, the engine counts, the outfmt 0 rows).
const CHUNK: usize = 1024;

// No NCBI counterpart: chooses which threads compute the report's values; it does not change any
// value NCBI computes.
/// The pool the report stage runs on (`x_pairwise_par::pool`: `Some` with the switch on and a
/// pool of two threads or more), or `None` for the serial stage.
pub(crate) fn pool<'p, 'q>(pool: &'p SearchPool<'q>) -> Option<&'p SearchPool<'q>> {
    let pooled = x_pairwise_par::pool(pool);
    if x_pairwise_par::on() {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRWISEPAR", pooled.is_some());
    }
    pooled
}

/// Rows `rows` of query `q_idx`'s rows; `header`: the query's first piece (with its outfmt 7
/// header).
struct Piece {
    q_idx: usize,
    header: bool,
    rows: std::ops::Range<usize>,
}

/// `report::write_tabular` without a format observer, with the rows formatted on `pool` (item
/// 3).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tabular<W: Write>(
    pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    writer: &mut W,
    comments: bool,
    query_titles: &[&[u8]],
    database: &str,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    unsearched: &[bool],
    epilog: Option<usize>,
    warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = crate::report::outfmt6::OutputConfig::ncbi_compat();
    let mut by_query: Vec<Vec<&TblastxHsp>> = vec![Vec::new(); query_titles.len()];
    for hsp in hits {
        by_query[hsp.hit.q_idx as usize].push(hsp);
    }
    // Each query's pieces of at most `PIECE_ROWS` rows; a query without rows is one piece.
    let pieces = by_query.iter().enumerate().flat_map(|(q_idx, query_hits)| {
        let rows = query_hits.len();
        (0..rows.div_ceil(PIECE_ROWS).max(1)).map(move |piece| Piece {
            q_idx,
            header: piece == 0,
            rows: piece * PIECE_ROWS..((piece + 1) * PIECE_ROWS).min(rows),
        })
    });
    let format = |piece: &Piece| {
        let query_hits = &by_query[piece.q_idx];
        let mut buffer = Vec::with_capacity(piece.rows.len() * 96);
        let result = (|| {
            if comments && piece.header {
                report::write_outfmt7_query_header(
                    &mut buffer,
                    query_titles[piece.q_idx],
                    database,
                    (!unsearched[piece.q_idx]).then_some(query_hits.len()),
                )?;
            }
            for hsp in &query_hits[piece.rows.clone()] {
                report::write_tabular_row(&mut buffer, hsp, query_ids, subject_ids, &config)?;
            }
            Ok(())
        })();
        (buffer, result)
    };
    // The query's warnings come before its lines (`QueryWarnings`), before its first piece.
    x_pairwise_par::write_in_order(
        writer,
        pool,
        pieces,
        x_pairwise_par::row_budget(pool),
        |piece| piece.rows.len().max(1),
        |piece| piece.header.then_some(piece.q_idx),
        format,
        warnings,
    )?;
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2233-2238
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
    //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
    //         tabinfo.PrintNumProcessed(m_QueriesFormatted);
    // ```
    // The same last line as `report::write_tabular`.
    if let (true, Some(queries_before)) = (comments, epilog) {
        writeln!(
            writer,
            "# BLAST processed {} queries",
            queries_before + query_titles.len()
        )?;
    }
    Ok(())
}

// The parallel forms of items 1 and 2; without threads (the second forms) they are the serial
// loops.
/// `report::set_displayed_identities` over chunks of `hits` on `pool` (item 1).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn set_displayed_identities<S: TbxScore>(
    pool: &SearchPool<'_>,
    hits: &mut [TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) {
    pool.install(|| {
        hits.par_chunks_mut(CHUNK).for_each(|chunk| {
            report::set_displayed_identities(chunk, queries, subjects, query_code, db_code, sc)
        })
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
pub(crate) fn set_displayed_identities<S: TbxScore>(
    _pool: &SearchPool<'_>,
    hits: &mut [TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) {
    report::set_displayed_identities(hits, queries, subjects, query_code, db_code, sc)
}

/// `report::engine_counts` (the Seq-aligns' engine counts, outfmt 8 and 12 only) over chunks of
/// `hits` on `pool`, as item 1: each chunk counted by the serial function, the chunks' counts
/// joined in the list's order.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn engine_counts<S: TbxScore>(
    pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) -> Vec<(u32, u32)> {
    let chunks: Vec<Vec<(u32, u32)>> = pool.install(|| {
        hits.par_chunks(CHUNK)
            .map(|chunk| report::engine_counts(chunk, queries, subjects, query_code, db_code, sc))
            .collect()
    });
    chunks.concat()
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
pub(crate) fn engine_counts<S: TbxScore>(
    _pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
    sc: S,
) -> Vec<(u32, u32)> {
    report::engine_counts(hits, queries, subjects, query_code, db_code, sc)
}

/// The outfmt 0 rows of `hits` (`report::pairwise_hits`) as an indexed map on `pool` (item 2).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn pairwise_hits(
    pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    rows: &PairwiseRows<'_>,
) -> Vec<PairwiseHit> {
    pool.install(|| {
        hits.par_iter()
            .with_min_len(CHUNK / 4)
            .map_init(Vec::new, |window, hsp| rows.hit(hsp, window))
            .collect()
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
pub(crate) fn pairwise_hits(
    _pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    rows: &PairwiseRows<'_>,
) -> Vec<PairwiseHit> {
    let mut window = Vec::new();
    hits.iter().map(|hsp| rows.hit(hsp, &mut window)).collect()
}

//! EXPERIMENT (LOSAT_X_TBXRPTPAR, shadow LOSAT_X_TBXRPTPARSHADOW): TBLASTX's report stage
//! (`write_tblastx_outputs` in blast_engine/run_impl.rs) on the search's pool.
//!
//! The reports are written once, after the last query batch, by the thread of the batch loop.
//! Before this switch every part of that stage ran on that thread while the other threads of
//! the pool waited: on the genome 3 x 3 input at 16 threads (ws, 2026-10-19) the stage took
//! 1.31 s of 8.6 s with outfmt 6 and 3.37 s of 10.6 s with outfmt 0, of which the displayed
//! identities took 0.99 s (each HSP's codons translated again from the nucleotides), the
//! outfmt 0 rows 1.80 s (translated a second time, with the masks lowercased), the tabular rows
//! 0.17 s (outfmt 0's text 0.41 s) and the final order (`report_order`) 0.15 s.
//!
//! What the switch changes (read once per process, `mode`), when the pool has two threads or
//! more (`pool`):
//! 1. The displayed identities (`report::set_displayed_identities`) run over chunks of the final
//!    list on the pool, each chunk with the serial function.
//! 2. The outfmt 0 rows (`report::pairwise_hits`, also the Web adapter's hit records) are an
//!    indexed parallel map of `PairwiseRows::hit` over the final list, collected in its order.
//! 3. The outfmt 6 / 7 rows of a query (`report::write_tabular`) with more than `CHUNK` rows are
//!    formatted into one buffer per chunk of `CHUNK` rows by the pool's threads, in windows of
//!    `WINDOW_CHUNKS` chunks per thread; the loop's thread writes a window's buffers in order
//!    before the next window is formatted, so at most one window of rows is held unwritten. The
//!    warnings before each query, the outfmt 7 headers and the last line stay on the loop's
//!    thread, in the serial order. With the Web adapter's format observer (`FormatProbe`, which
//!    flushes around every row) the rows are written by the serial function.
//! The final order (`report_order`), the outfmt 0 text (`write_tblastx_pairwise_report`) and the
//! writes stay on the loop's thread.
//!
//! Why the bytes are those of the serial stage: every value is computed by the serial code's
//! function of one HSP (or of one row) on the same inputs (the HSP, the query and subject
//! letters, the genetic codes, the placements, the masks), which reads nothing another HSP
//! writes; the order of the final list is not changed (chunks are disjoint slices; the map is
//! collected by index); the rows of a query are the serial rows, concatenated in chunk order.
//! No floating-point value is summed across HSPs. The same bytes go to the writer in the same
//! order, in larger writes; a failed write stops the run as before.
//!
//! Off (the serial stage): the switch off (`LOSAT_X_TBXRPTPAR=0`, or `LOSAT_X_DEFAULT=0` and it
//! is not set); a pool of one thread or none (`-num_threads 1`, builds without the `parallel`
//! feature, Wasm without threads).
//!
//! Memory: the outfmt 0 rows are one list of the final list's length in both forms; the tabular
//! windows hold at most `CHUNK` x `WINDOW_CHUNKS` x threads rows formatted (~13 MB at 16
//! threads for rows of ~100 bytes).
//!
//! The shadow (`LOSAT_X_TBXRPTPARSHADOW`, implies the switch) runs the serial functions too, on
//! a copy: it compares every HSP's displayed identities, positives, mismatches and identity
//! (bits), every outfmt 0 row (`Debug`), and the bytes of every tabular report (both written
//! into memory, without warnings); on a difference it prints `[X_SHADOW] TBXRPTPAR ...
//! DIFFERENT` and panics; after the reports it prints `[X_SHADOW] TBXRPTPAR hsps= pairwise=
//! bytes= identical` (the HSPs, outfmt 0 rows and tabular bytes compared).

use std::io::{self, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::pairwise::PairwiseHit;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::threading::SearchPool;

use super::blast_engine::TblastxHsp;
use super::report::{self, PairwiseRows, SubjectLetters};

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// The rows (or HSPs) of one task.
const CHUNK: usize = 1024;

/// The chunks of tabular rows formatted per thread before the loop's thread writes them.
const WINDOW_CHUNKS: usize = 4;

// No NCBI counterpart: what the shadow compared (LOSAT_X_TBXRPTPARSHADOW); it does not change any
// value NCBI computes.
static SHADOW_HSPS: AtomicU64 = AtomicU64::new(0);
static SHADOW_PAIRWISE: AtomicU64 = AtomicU64::new(0);
static SHADOW_BYTES: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: chooses which threads compute the report's values; it does not change any
// value NCBI computes. Read once per process.
/// 0 = off, 1 = LOSAT_X_TBXRPTPAR (a strict-set switch: on unless `LOSAT_X_TBXRPTPAR` is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set; `utils/xswitch.rs`), 2 = LOSAT_X_TBXRPTPARSHADOW.
pub(crate) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_TBXRPTPARSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_TBXRPTPAR") {
            1
        } else {
            0
        }
    })
}

/// The pool the report stage runs on (`Some` with the switch on and a pool of two threads or
/// more), or `None` for the serial stage.
pub(crate) fn pool<'p, 'q>(pool: &'p SearchPool<'q>) -> Option<&'p SearchPool<'q>> {
    const PARALLEL: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));
    (PARALLEL && mode() != 0 && pool.enabled() && pool.threads() > 1).then_some(pool)
}

/// `report::set_displayed_identities` over chunks of `hits` on `pool` (item 1).
pub(crate) fn set_displayed_identities(
    pool: &SearchPool<'_>,
    hits: &mut [TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) {
    let reference = (mode() == 2).then(|| {
        let mut copy = hits.to_vec();
        report::set_displayed_identities(&mut copy, queries, subjects, query_code, db_code);
        copy
    });
    identities_on(pool, hits, queries, subjects, query_code, db_code);
    if let Some(reference) = reference {
        for (index, (hsp, serial)) in hits.iter().zip(&reference).enumerate() {
            let (a, b) = (&hsp.hit, &serial.hit);
            if (
                a.num_ident,
                a.num_positives,
                a.mismatch,
                a.identity.to_bits(),
            ) != (
                b.num_ident,
                b.num_positives,
                b.mismatch,
                b.identity.to_bits(),
            ) {
                different(&format!(
                    "HSP {index}: identities {} / {}, positives {} / {}, mismatches {} / {}, identity {} / {}",
                    a.num_ident, b.num_ident, a.num_positives, b.num_positives, a.mismatch,
                    b.mismatch, a.identity, b.identity
                ));
            }
        }
        SHADOW_HSPS.fetch_add(hits.len() as u64, Ordering::Relaxed);
    }
}

/// The outfmt 0 rows of `hits` (`report::pairwise_hits`) as an indexed map on `pool` (item 2).
pub(crate) fn pairwise_hits(
    pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    rows: &PairwiseRows<'_>,
) -> Vec<PairwiseHit> {
    let pairwise = pairwise_on(pool, hits, rows);
    if mode() == 2 {
        let mut window = Vec::new();
        for (index, (hit, hsp)) in pairwise.iter().zip(hits).enumerate() {
            let (a, b) = (
                format!("{hit:?}"),
                format!("{:?}", rows.hit(hsp, &mut window)),
            );
            if a != b {
                different(&format!("outfmt 0 row {index}: {a} / {b}"));
            }
        }
        SHADOW_PAIRWISE.fetch_add(pairwise.len() as u64, Ordering::Relaxed);
    }
    pairwise
}

/// `report::write_tabular` without a format observer, with the rows of each query formatted on
/// `pool` (item 3).
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
    epilog: bool,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    if mode() == 2 {
        let (mut serial, mut parallel) = (Vec::new(), Vec::new());
        report::write_tabular(
            hits,
            &mut serial,
            comments,
            query_titles,
            database,
            query_ids,
            subject_ids,
            unsearched,
            epilog,
            None,
            None,
        )?;
        write_tabular_rows(
            pool,
            hits,
            &mut parallel,
            comments,
            query_titles,
            database,
            query_ids,
            subject_ids,
            unsearched,
            epilog,
            None,
        )?;
        if serial != parallel {
            let at = serial
                .iter()
                .zip(&parallel)
                .position(|(a, b)| a != b)
                .unwrap_or(serial.len().min(parallel.len()));
            different(&format!(
                "tabular bytes {} / {}, first difference at byte {at}",
                parallel.len(),
                serial.len()
            ));
        }
        SHADOW_BYTES.fetch_add(serial.len() as u64, Ordering::Relaxed);
    }
    write_tabular_rows(
        pool,
        hits,
        writer,
        comments,
        query_titles,
        database,
        query_ids,
        subject_ids,
        unsearched,
        epilog,
        warnings.as_deref_mut(),
    )
}

#[allow(clippy::too_many_arguments)]
fn write_tabular_rows<W: Write>(
    pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    writer: &mut W,
    comments: bool,
    query_titles: &[&[u8]],
    database: &str,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    unsearched: &[bool],
    epilog: bool,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = crate::report::outfmt6::OutputConfig::ncbi_compat();
    let mut by_query: Vec<Vec<&TblastxHsp>> = vec![Vec::new(); query_titles.len()];
    for hsp in hits {
        by_query[hsp.hit.q_idx as usize].push(hsp);
    }
    for (q_idx, query_hits) in by_query.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        if comments {
            report::write_outfmt7_query_header(
                writer,
                query_titles[q_idx],
                database,
                (!unsearched[q_idx]).then_some(query_hits.len()),
            )?;
        }
        if query_hits.len() <= CHUNK {
            for hsp in query_hits {
                report::write_tabular_row(writer, hsp, query_ids, subject_ids, &config)?;
            }
            continue;
        }
        for window in query_hits.chunks(CHUNK * WINDOW_CHUNKS * pool.threads()) {
            for buffer in format_rows(pool, window, query_ids, subject_ids, &config)? {
                writer.write_all(&buffer)?;
            }
        }
    }
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2233-2238
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
    //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
    //         tabinfo.PrintNumProcessed(m_QueriesFormatted);
    // ```
    // The same last line as `report::write_tabular`.
    if comments && epilog {
        writeln!(writer, "# BLAST processed {} queries", query_titles.len())?;
    }
    Ok(())
}

/// The rows of `rows` into one buffer (`report::write_tabular_row`).
fn format_chunk(
    rows: &[&TblastxHsp],
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    config: &crate::report::outfmt6::OutputConfig,
) -> io::Result<Vec<u8>> {
    let mut buffer = Vec::with_capacity(rows.len() * 96);
    for hsp in rows {
        report::write_tabular_row(&mut buffer, hsp, query_ids, subject_ids, config)?;
    }
    Ok(buffer)
}

// The parallel forms of items 1-3; without threads (the second forms) they are the serial loops.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn identities_on(
    pool: &SearchPool<'_>,
    hits: &mut [TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) {
    pool.install(|| {
        hits.par_chunks_mut(CHUNK).for_each(|chunk| {
            report::set_displayed_identities(chunk, queries, subjects, query_code, db_code)
        })
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn identities_on(
    _pool: &SearchPool<'_>,
    hits: &mut [TblastxHsp],
    queries: &[FastaRecord],
    subjects: SubjectLetters<'_>,
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) {
    report::set_displayed_identities(hits, queries, subjects, query_code, db_code)
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn pairwise_on(
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
fn pairwise_on(
    _pool: &SearchPool<'_>,
    hits: &[TblastxHsp],
    rows: &PairwiseRows<'_>,
) -> Vec<PairwiseHit> {
    let mut window = Vec::new();
    hits.iter().map(|hsp| rows.hit(hsp, &mut window)).collect()
}

/// The rows of `window`, one buffer per chunk of `CHUNK` rows, in order.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn format_rows(
    pool: &SearchPool<'_>,
    window: &[&TblastxHsp],
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    config: &crate::report::outfmt6::OutputConfig,
) -> io::Result<Vec<Vec<u8>>> {
    pool.install(|| {
        window
            .par_chunks(CHUNK)
            .map(|rows| format_chunk(rows, query_ids, subject_ids, config))
            .collect()
    })
}

/// The rows of `window`, one buffer per chunk of `CHUNK` rows, in order.
#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn format_rows(
    _pool: &SearchPool<'_>,
    window: &[&TblastxHsp],
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    config: &crate::report::outfmt6::OutputConfig,
) -> io::Result<Vec<Vec<u8>>> {
    window
        .chunks(CHUNK)
        .map(|rows| format_chunk(rows, query_ids, subject_ids, config))
        .collect()
}

// No NCBI counterpart: the shadow's line after the reports (`LOSAT_X_TBXRPTPARSHADOW`); it does not
// change any value NCBI computes.
/// The shadow's line after the reports (`LOSAT_X_TBXRPTPARSHADOW`).
pub(crate) fn print_shadow_summary() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] TBXRPTPAR hsps={} pairwise={} bytes={} identical",
            SHADOW_HSPS.load(Ordering::Relaxed),
            SHADOW_PAIRWISE.load(Ordering::Relaxed),
            SHADOW_BYTES.load(Ordering::Relaxed)
        );
    }
}

fn different(what: &str) -> ! {
    eprintln!("[X_SHADOW] TBXRPTPAR DIFFERENT ({what})");
    panic!("LOSAT_X_TBXRPTPARSHADOW: the report stage on the pool differs from the serial stage");
}

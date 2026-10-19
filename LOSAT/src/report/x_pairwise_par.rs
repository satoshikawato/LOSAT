//! EXPERIMENT (LOSAT_X_PAIRWISEPAR, shadow LOSAT_X_PAIRWISEPARSHADOW): the per-query parts of the
//! pairwise reports (outfmt 0) of BLASTN, TBLASTX and BLASTP formatted on the search's pool.
//!
//! The writers (`write_blastn_pairwise_report`, `write_tblastx_pairwise_report`,
//! `write_blastp_pairwise_report` in `pairwise.rs`) run on the thread of the search's loop while
//! the other threads of the pool wait. megablast q100 x 16S (`-subject`, 12 NCBI query batches,
//! LOSAT_X_RPTSTREAM writes each batch's reports after its search) spent 1.4-2.0 s of 7-10 s on
//! that thread at 8 and 16 threads (local, 2026-10-19); TBLASTX's outfmt 0 text 0.41 s of the
//! genome 3 x 3 input at 16 threads, after LOSAT_X_TBXRPTPAR.
//!
//! The three writers have one shape (`write`): after the prolog, for each query its warnings
//! (`QueryWarnings`), then its head (`QuerySections::head`: `\n\n`, the query's header, and
//! either "No hits found" with its footer or the description table and a blank line), the
//! alignments of each of the first `num_alignments` subjects in the order of the table
//! (`QuerySections::subject`) and its footer (`QuerySections::foot`); the epilog after the last.
//! Each part reads only its query's hits (grouped by subject, `QueryGroup`), the subject titles
//! and the report's fixed data.
//!
//! What the switch changes (read once per process, `mode`), when the caller gives a pool of two
//! threads or more (`pool`) and there is no format observer (`FormatProbe`, which flushes around
//! every subject and HSP): the queries are taken in windows whose estimated bytes (`estimate`)
//! reach `WINDOW_BYTES_PER_THREAD` per thread (at least one query). The pool groups each query's
//! hits, then formats the window's pieces (a query's head, runs of its alignments of about
//! `PIECE_BYTES` estimated bytes, cut between HSPs, also within one subject's, and its footer)
//! into one buffer each; the calling thread then writes the buffers in order, each query's
//! warnings before its head, before the next window.
//!
//! Why the bytes are those of the serial writer: every piece is the serial writer's calls on the
//! same arguments (the same functions, `write_serial` calls them on the writer), into a buffer
//! instead of the writer. A subject's alignments cut between HSPs are its heading with the first
//! HSP, then each HSP's text, which reads only that HSP and the subject's title
//! (`write_subject_alignments_part`). The pieces are collected by index and written in the serial
//! order; the warnings come where the serial writer writes them (before the query's head, after a
//! flush). No value is summed across queries. The same bytes go to the writer in larger writes; a failed
//! write stops the report as before.
//!
//! Off (the serial writer, `write_serial`): the switch off (`LOSAT_X_PAIRWISEPAR=0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set); no pool or one of one thread (`-num_threads 1`,
//! builds without the `parallel` feature, Wasm without threads); a format observer (the Web
//! adapter); the callers' own conditions (BLASTN: while LOSAT_X_QBATCHPAR slots run).
//!
//! Memory: one window of formatted text, ~`WINDOW_BYTES_PER_THREAD` x threads by the estimate
//! (16 MiB at 16 threads), or one query's text when one query is larger.
//!
//! The shadow (`LOSAT_X_PAIRWISEPARSHADOW`, implies the switch) also writes the report by the
//! serial writer into memory and compares the bytes (without warnings); on a difference it
//! prints `[X_SHADOW] PAIRWISEPAR ... DIFFERENT` and panics; after each report it prints
//! `[X_SHADOW] PAIRWISEPAR queries= windows= pieces= bytes= identical`.

use std::collections::HashMap;
use std::io::{self, Write};
use std::sync::{Arc, OnceLock};

use super::pairwise::PairwiseHit;
use super::query_warnings::QueryWarnings;
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// The estimated bytes of a run of subjects' alignments formatted as one piece.
const PIECE_BYTES: usize = 256 << 10;

/// The estimated bytes of a window per thread of the pool.
const WINDOW_BYTES_PER_THREAD: usize = 1 << 20;

// No NCBI counterpart: chooses which threads format the report; it does not change any value NCBI
// computes. Read once per process.
/// 0 = off, 1 = LOSAT_X_PAIRWISEPAR (a strict-set switch: on unless `LOSAT_X_PAIRWISEPAR` is `0`,
/// or `LOSAT_X_DEFAULT=0` and it is not set; `utils/xswitch.rs`), 2 = LOSAT_X_PAIRWISEPARSHADOW.
pub(crate) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_PAIRWISEPARSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_PAIRWISEPAR") {
            1
        } else {
            0
        }
    })
}

/// The pool the pairwise writers format on (`Some` with the switch on and a pool of two threads
/// or more), or `None` for the serial writer.
pub(crate) fn pool<'p, 'q>(pool: &'p SearchPool<'q>) -> Option<&'p SearchPool<'q>> {
    const PARALLEL: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));
    (PARALLEL && mode() != 0 && pool.enabled() && pool.threads() > 1).then_some(pool)
}

/// One query's hits by subject, as the serial writers group them: the subjects in the order of
/// their first hit, each subject's hits and their indices in the report's hit list, and the
/// subjects' title bytes (by subject index).
pub(crate) struct QueryGroup<'h> {
    pub(crate) subject_order: Vec<u32>,
    pub(crate) subject_hits: HashMap<u32, Vec<&'h PairwiseHit>>,
    pub(crate) subject_hit_indices: HashMap<u32, Vec<HspIndex>>,
    pub(crate) subject_titles: HashMap<u32, Vec<u8>>,
}

impl<'h> QueryGroup<'h> {
    /// The group of `query_hits`; `None` for a query without hits.
    pub(crate) fn new(
        query_hits: &[(HspIndex, &'h PairwiseHit)],
        subject_titles: &[Arc<[u8]>],
    ) -> Option<Self> {
        if query_hits.is_empty() {
            return None;
        }
        let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
        let mut subject_hit_indices: HashMap<u32, Vec<HspIndex>> = HashMap::new();
        let mut subject_order: Vec<u32> = Vec::new();
        for &(hsp_index, hit) in query_hits {
            let s_idx = hit.hit.s_idx;
            if !subject_hits.contains_key(&s_idx) {
                subject_order.push(s_idx);
            }
            subject_hits.entry(s_idx).or_default().push(hit);
            subject_hit_indices
                .entry(s_idx)
                .or_default()
                .push(hsp_index);
        }
        let subject_titles: HashMap<u32, Vec<u8>> = subject_order
            .iter()
            .map(|&s_idx| {
                let title = subject_titles
                    .get(s_idx as usize)
                    .map_or(&b""[..], |title| title.as_ref());
                (s_idx, title.to_vec())
            })
            .collect();
        Some(QueryGroup {
            subject_order,
            subject_hits,
            subject_hit_indices,
            subject_titles,
        })
    }

    /// The subjects of the description table (the first `num_descriptions`).
    pub(crate) fn described(&self, num_descriptions: usize) -> &[u32] {
        &self.subject_order[..self.subject_order.len().min(num_descriptions)]
    }

    /// The subjects with alignments (the first `num_alignments`).
    ///
    /// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    /// (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
    pub(crate) fn aligned(&self, num_alignments: usize) -> &[u32] {
        &self.subject_order[..self.subject_order.len().min(num_alignments)]
    }
}

/// The parts of one program's per-query report (`write`).
pub(crate) trait QuerySections: Sync {
    /// The subjects with alignments per query (`QueryGroup::aligned`).
    fn num_alignments(&self) -> usize;

    /// `\n\n`, the query's header, then "No hits found" and the query's footer (`group` is
    /// `None`), or the description table and the blank line after it.
    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&QueryGroup<'_>>,
    ) -> io::Result<()>;

    /// The alignments of the HSPs `part` of subject `s_idx` of the query (`group`), the
    /// subject's heading with its first HSP: the parts of `0..` its HSPs in order write the
    /// subject's alignments.
    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
    ) -> io::Result<()>;

    /// The footer of a query with hits.
    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()>;
}

/// Writes the per-query reports of `hits_by_query` (each query's hits with their indices in
/// the report's hit list) between a writer's prolog and epilog: on `pool` (`pool`) without a
/// `probe`, else by the serial writer.
pub(crate) fn write<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut QueryWarnings<'_>>,
    pool: Option<&SearchPool<'_>>,
) -> io::Result<()> {
    let pool = match pool {
        Some(pool) if probe.is_none() => pool,
        _ => {
            return write_serial(
                writer,
                sections,
                hits_by_query,
                subject_titles,
                probe,
                warnings,
            )
        }
    };
    if mode() == 2 {
        let (mut serial, mut parallel) = (Vec::new(), Vec::new());
        write_serial(
            &mut serial,
            sections,
            hits_by_query,
            subject_titles,
            None,
            None,
        )?;
        let (windows, pieces) = write_parallel(
            &mut parallel,
            sections,
            hits_by_query,
            subject_titles,
            None,
            pool,
        )?;
        if serial != parallel {
            let at = serial
                .iter()
                .zip(&parallel)
                .position(|(a, b)| a != b)
                .unwrap_or(serial.len().min(parallel.len()));
            different(&format!(
                "bytes {} / {}, first difference at byte {at}",
                parallel.len(),
                serial.len()
            ));
        }
        eprintln!(
            "[X_SHADOW] PAIRWISEPAR queries={} windows={windows} pieces={pieces} bytes={} identical",
            hits_by_query.len(),
            serial.len()
        );
    }
    write_parallel(
        writer,
        sections,
        hits_by_query,
        subject_titles,
        warnings,
        pool,
    )?;
    Ok(())
}

/// The serial writer: each query's warnings and parts on `writer`, in order.
fn write_serial<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    for (q_idx, query_hits) in hits_by_query.iter().enumerate() {
        // The query's warnings (and those of its batch's reading) come before its preamble
        // (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let group = QueryGroup::new(query_hits, subject_titles);
        sections.head(writer, q_idx, group.as_ref())?;
        if let Some(group) = &group {
            for &s_idx in group.aligned(sections.num_alignments()) {
                let hsps = group.subject_hits[&s_idx].len();
                sections.subject(writer, probe.as_deref_mut(), group, s_idx, 0..hsps)?;
            }
            sections.foot(writer, q_idx)?;
        }
    }
    Ok(())
}

/// The estimated bytes of a hit's alignment: its header lines and three rows of ~80 bytes per
/// 60 columns. Only the windows' and pieces' sizes depend on it.
fn estimate(hit: &PairwiseHit) -> usize {
    400 + 4 * hit.hit.length
}

/// What a piece of a query's report holds.
#[derive(Clone, Copy)]
enum Part {
    Head,
    /// The alignments from (position in the aligned subjects, HSP of that subject) up to the
    /// second, exclusive: a run of whole subjects, or part of one subject's HSPs.
    Alignments((usize, usize), (usize, usize)),
    Foot,
}

#[derive(Clone, Copy)]
struct Piece {
    query: usize,
    part: Part,
}

/// The pieces of the queries `first..` (`groups`, in order), each query's head first.
fn pieces<S: QuerySections>(
    sections: &S,
    first: usize,
    groups: &[Option<QueryGroup<'_>>],
) -> Vec<Piece> {
    let mut pieces = Vec::new();
    for (offset, group) in groups.iter().enumerate() {
        let query = first + offset;
        pieces.push(Piece {
            query,
            part: Part::Head,
        });
        let Some(group) = group else {
            continue;
        };
        let aligned = group.aligned(sections.num_alignments());
        let (mut start, mut bytes, mut end) = ((0, 0), 0, (0, 0));
        for (position, s_idx) in aligned.iter().enumerate() {
            for (hsp, hit) in group.subject_hits[s_idx].iter().enumerate() {
                bytes += estimate(hit);
                end = (position, hsp + 1);
                if bytes >= PIECE_BYTES {
                    pieces.push(Piece {
                        query,
                        part: Part::Alignments(start, end),
                    });
                    (start, bytes) = (end, 0);
                }
            }
        }
        if start != end {
            pieces.push(Piece {
                query,
                part: Part::Alignments(start, end),
            });
        }
        pieces.push(Piece {
            query,
            part: Part::Foot,
        });
    }
    pieces
}

/// One piece into its own buffer, by the serial writer's calls.
fn format_piece<S: QuerySections>(
    sections: &S,
    group: Option<&QueryGroup<'_>>,
    piece: Piece,
) -> io::Result<Vec<u8>> {
    let mut buffer = Vec::new();
    match piece.part {
        Part::Head => sections.head(&mut buffer, piece.query, group)?,
        Part::Alignments((first, first_hsp), (last, end_hsp)) => {
            let group = group.expect("a query with alignment pieces has hits");
            let aligned = group.aligned(sections.num_alignments());
            for (position, &s_idx) in aligned.iter().enumerate().take(last + 1).skip(first) {
                let from = if position == first { first_hsp } else { 0 };
                let to = if position == last {
                    end_hsp
                } else {
                    group.subject_hits[&s_idx].len()
                };
                if from < to {
                    sections.subject(&mut buffer, None, group, s_idx, from..to)?;
                }
            }
        }
        Part::Foot => sections.foot(&mut buffer, piece.query)?,
    }
    Ok(buffer)
}

/// The parallel writer (module header): windows of queries, their pieces formatted on `pool`,
/// written in order by the calling thread. Returns the numbers of windows and pieces.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn write_parallel<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    mut warnings: Option<&mut QueryWarnings<'_>>,
    pool: &SearchPool<'_>,
) -> io::Result<(usize, usize)> {
    let budget = WINDOW_BYTES_PER_THREAD.saturating_mul(pool.threads());
    let (mut windows, mut pieces_written) = (0, 0);
    let mut first = 0;
    while first < hits_by_query.len() {
        let mut end = first;
        let mut bytes = 0usize;
        while end < hits_by_query.len() && (end == first || bytes < budget) {
            bytes += hits_by_query[end]
                .iter()
                .map(|&(_, hit)| estimate(hit))
                .sum::<usize>();
            end += 1;
        }
        let window = &hits_by_query[first..end];
        let (planned, buffers) = pool.install(|| {
            let groups: Vec<Option<QueryGroup<'_>>> = window
                .par_iter()
                .map(|query_hits| QueryGroup::new(query_hits, subject_titles))
                .collect();
            let planned = pieces(sections, first, &groups);
            let buffers: Vec<io::Result<Vec<u8>>> = planned
                .par_iter()
                .map(|&piece| format_piece(sections, groups[piece.query - first].as_ref(), piece))
                .collect();
            (planned, buffers)
        });
        for (piece, buffer) in planned.iter().zip(buffers) {
            let buffer = buffer?;
            if matches!(piece.part, Part::Head) {
                if let Some(warnings) = warnings.as_deref_mut() {
                    warnings.before_query(piece.query, &mut *writer)?;
                }
            }
            writer.write_all(&buffer)?;
        }
        windows += 1;
        pieces_written += planned.len();
        first = end;
    }
    Ok((windows, pieces_written))
}

/// Without threads `pool` returns `None`; the pieces are formatted in order here (the shadow's
/// arm only).
#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn write_parallel<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    mut warnings: Option<&mut QueryWarnings<'_>>,
    _pool: &SearchPool<'_>,
) -> io::Result<(usize, usize)> {
    let groups: Vec<Option<QueryGroup<'_>>> = hits_by_query
        .iter()
        .map(|query_hits| QueryGroup::new(query_hits, subject_titles))
        .collect();
    let planned = pieces(sections, 0, &groups);
    for &piece in &planned {
        let buffer = format_piece(sections, groups[piece.query].as_ref(), piece)?;
        if matches!(piece.part, Part::Head) {
            if let Some(warnings) = warnings.as_deref_mut() {
                warnings.before_query(piece.query, &mut *writer)?;
            }
        }
        writer.write_all(&buffer)?;
    }
    Ok((1, planned.len()))
}

fn different(what: &str) -> ! {
    eprintln!("[X_SHADOW] PAIRWISEPAR DIFFERENT ({what})");
    panic!("LOSAT_X_PAIRWISEPARSHADOW: the pairwise report formatted on the pool differs from the serial writer");
}

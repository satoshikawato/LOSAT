//! EXPERIMENT (LOSAT_X_PAIRWISEPAR): the per-query parts of the pairwise reports (outfmt 0, and the
//! query-anchored 1 to 4) of BLASTN, TBLASTX, BLASTP, TBLASTN and BLASTX (RPSTBLASTN with BLASTX's
//! writer) formatted on the search's pool, and BLASTN's and BLASTP's outfmt 0
//! rows (`algorithm/blastn/pairwise.rs`, `pairwise_hits`; `algorithm/blastp/blast_engine.rs`,
//! `build_pairwise_hits`) made on it.
//!
//! The writers (`write_blastn_pairwise_report`, `write_tblastx_pairwise_report`,
//! `write_blastp_pairwise_report` in `pairwise.rs`) run on the thread of the search's loop while
//! the other threads of the pool wait. megablast q100 x 16S (`-subject`, 12 NCBI query batches,
//! LOSAT_X_RPTSTREAM writes each batch's reports after its search) spent 1.4-2.0 s of 7-10 s on
//! that thread at 8 and 16 threads (local, 2026-10-19), 0.26 s of it in the rows; TBLASTX's
//! outfmt 0 text 0.41 s of the genome 3 x 3 input at 16 threads, after LOSAT_X_TBXRPTPAR.
//!
//! The writers of the five programs have one shape (`write`; TBLASTN gives it no pool, BLASTX
//! the pool of a query batch whose searches are done, `write_blastx_pairwise_report`): after the
//! prolog, for each query its warnings
//! (`QueryWarnings`), then its head (`QuerySections::head`: `\n\n`, the query's header, and
//! either "No hits found" with its footer or the description table and a blank line), the
//! alignments of each of the first `num_alignments` subjects in the order of the table
//! (`QuerySections::subject`) and its footer (`QuerySections::foot`); the epilog after the last.
//! Each part reads only its query's hits (grouped by subject, `QueryGroup`), the subject titles
//! and the report's fixed data.
//!
//! What the switch changes (read once per process, `mode`), when the caller gives a pool of two
//! threads or more (`pool`) and there is no format observer (`FormatProbe`, which flushes around
//! every subject and HSP; the Web adapter): the queries are taken in windows whose estimated
//! bytes (`query_estimate`) reach `WINDOW_BYTES_PER_THREAD` per thread (at least one query), and
//! the pool groups each window's hits. The window's pieces (a query's head, runs of its
//! alignments of about `PIECE_BYTES` estimated bytes, cut between HSPs, also within one
//! subject's, or the one section of a query-anchored report, `QuerySections::alignments`, and its
//! footer) are formatted by the pool into one buffer each, in chunks of at
//! most the same estimated bytes (at least one piece); the calling thread writes a chunk's
//! buffers in order, each query's warnings before its head, before the next chunk is formatted.
//! BLASTN's rows: each subject's HSP list of each query on the pool, collected in the list order
//! (also with a format observer: the rows are values, not writes). BLASTP's rows: each hit of the
//! final list on the pool, collected in the list order, the first error in that order.
//!
//! Why the bytes are those of the serial writer: every piece is the serial writer's calls on the
//! same arguments (the same functions, `write_serial` calls them on the writer), into a buffer
//! instead of the writer. A subject's alignments cut between HSPs are its heading with the first
//! HSP, then each HSP's text, which reads only that HSP and the subject's title
//! (`write_subject_alignments_part`). The pieces are collected by index and written in the serial
//! order; the warnings come where the serial writer writes them (before the query's head, after a
//! flush). A piece whose call fails has its bytes up to the failure written, then the error, as
//! the serial writer's buffer would. No value is summed across queries. The same bytes go to the
//! writer in larger writes; a failed write stops the report as before. Each row is the serial
//! function of its HSP list.
//!
//! Off (the serial writer, `write_serial`): the switch off (`LOSAT_X_PAIRWISEPAR=0`);
//! no pool or one of one thread (`-num_threads 1`,
//! builds without the `parallel` feature, Wasm without threads); a format observer (the text
//! only); a caller that gives no pool.
//!
//! Memory: formatted text of at most `WINDOW_BYTES_PER_THREAD` x threads estimated bytes (16 MiB
//! at 16 threads) plus one piece (`PIECE_BYTES` plus one HSP) at a time, also for one large
//! query; the groups of one window of queries. The estimates count the report's lines
//! (`hit_estimate` at its `-line_length`), so the budget holds at every line length.
//!
//! BLASTN's and TBLASTX's outfmt 0 make their rows a window of queries at a time (line DBID, ID4f
//! and ID4t): the windows of `query_windows` over the queries' estimates (`query_estimate_of` of
//! their HSPs' lengths, before any row is made) with `window_bytes` of the pool (one thread's
//! budget without a pool, also with the switch off), each written (`write` of the window,
//! `first_query` its first query) before the next window's rows are made, as NCBI formats one
//! query's result set at a time (`PrintOneResultSet`). With the same budget a row window is one
//! window of the parallel writer. Not where the caller's hit records or a format observer read
//! every row (the Web adapter; `algorithm::blastn::blast_engine::run`,
//! `algorithm::tblastx::blast_engine::run_impl::rows_up_front`).

use std::collections::HashMap;
use std::io::{self, Write};
use std::sync::{Arc, OnceLock};

use super::pairwise::PairwiseHit;
use super::query_warnings::{QueryWarnings, SubjectReadHook};
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
/// LOSAT_X_PAIRWISEPAR, a strict-set switch: on unless it is `0` (`utils/xswitch.rs`).
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRWISEPAR"))
}

/// The pool the pairwise writers format on (`Some` with the switch on and a pool of two threads
/// or more), or `None` for the serial writer.
pub(crate) fn pool<'p, 'q>(pool: &'p SearchPool<'q>) -> Option<&'p SearchPool<'q>> {
    const PARALLEL: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));
    (PARALLEL && on() && pool.enabled() && pool.threads() > 1).then_some(pool)
}

/// One query's hits by subject, as the serial writers group them: the subjects in the order of
/// their first hit, each subject's hits and their indices in the report's hit list, and the
/// subjects' title bytes (by subject index). The hits of a subject that the report shows in the
/// hit of the subject before it (`QuerySections::joins`) are that subject's.
pub(crate) struct QueryGroup<'h> {
    pub(crate) subject_order: Vec<u32>,
    pub(crate) subject_hits: HashMap<u32, Vec<&'h PairwiseHit>>,
    pub(crate) subject_hit_indices: HashMap<u32, Vec<HspIndex>>,
    pub(crate) subject_titles: HashMap<u32, Vec<u8>>,
    /// The rank of each subject's alignment heading (`heading_rank`), made at the first `-html`
    /// heading of the query.
    heading_ranks: std::sync::OnceLock<HashMap<u32, usize>>,
}

impl<'h> QueryGroup<'h> {
    /// The group of `query_hits`; `None` for a query without hits. `joins(previous, s_idx)`:
    /// whether the hits of subject `s_idx`, whose first hit follows a hit of subject
    /// `previous`, are shown in `previous`'s hit (`QuerySections::joins`).
    pub(crate) fn new(
        query_hits: &[(HspIndex, &'h PairwiseHit)],
        subject_titles: &[Arc<[u8]>],
        joins: impl Fn(u32, u32) -> bool,
    ) -> Option<Self> {
        if query_hits.is_empty() {
            return None;
        }
        let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
        let mut subject_hit_indices: HashMap<u32, Vec<HspIndex>> = HashMap::new();
        let mut subject_order: Vec<u32> = Vec::new();
        let shown = shown_subjects(query_hits.iter().map(|(_, hit)| hit.hit.s_idx), joins);
        for (&(hsp_index, hit), shown) in query_hits.iter().zip(shown) {
            if !subject_hits.contains_key(&shown) {
                subject_order.push(shown);
            }
            subject_hits.entry(shown).or_default().push(hit);
            subject_hit_indices
                .entry(shown)
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
            heading_ranks: std::sync::OnceLock::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1929, 2307
    /// The rank of subject `s_idx`'s alignment heading (NCBI's `m_cur_align`, the `blast_rank` of
    /// its `-html` links, from 1 for each query): `x_PrintDefLine` counts each heading it starts,
    /// which is once for a subject and once for each alignment of a subject whose heading throws
    /// (`throws`: its title cannot be made; NCBI tries again at each alignment, as
    /// `write_hsps_or_skips` notes). The ranks are made once, at the first call (only `-html`
    /// headings call it); a heading that the registry fails throws in every subject, so it shows
    /// no rank.
    pub(crate) fn heading_rank(&self, s_idx: u32, throws: impl Fn(u32) -> bool) -> usize {
        let ranks = self.heading_ranks.get_or_init(|| {
            let mut started = 0;
            self.subject_order
                .iter()
                .map(|&subject| {
                    let rank = started + 1;
                    started += if throws(subject) {
                        self.subject_hits[&subject].len()
                    } else {
                        1
                    };
                    (subject, rank)
                })
                .collect()
        });
        ranks[&s_idx]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1523-1530
    /// The group of `query_hits` as `new` makes it, in the order of `sections`' `-sorthits` and
    /// `-sorthsps` (`sort`): NCBI sorts one query's Seq-align set before its description table
    /// and its alignments.
    pub(crate) fn for_sections<S: QuerySections>(
        query_hits: &[(HspIndex, &'h PairwiseHit)],
        subject_titles: &[Arc<[u8]>],
        sections: &S,
    ) -> Option<Self> {
        let (sort, tblastx) = sections.sort();
        Self::sorted(
            query_hits,
            subject_titles,
            |previous, s_idx| sections.joins(previous, s_idx),
            sort,
            tblastx,
        )
    }

    /// The group of `query_hits` grouped by `joins` (`new`) and sorted by `sort` (`sort`): the
    /// one place where the writer's grouping and order of one query's subjects are built, for
    /// the report itself (`for_sections`) and for what is decided before it is written (the
    /// described subjects of BLASTP's refusals, `pairwise::blastp_report_group`), so that the
    /// refusal model cannot drift from the report.
    pub(crate) fn sorted(
        query_hits: &[(HspIndex, &'h PairwiseHit)],
        subject_titles: &[Arc<[u8]>],
        joins: impl Fn(u32, u32) -> bool,
        sort: super::hit_sort::ReportSort,
        tblastx: bool,
    ) -> Option<Self> {
        let mut group = Self::new(query_hits, subject_titles, joins)?;
        group.sort(sort, tblastx);
        Some(group)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1743-1781
    /// Puts the subjects and each subject's HSPs (with their indices) in the order of NCBI's
    /// `SortSeqalignForSortableFormat` (`hit_sort::sort_query_hits`): a hit of NCBI is one
    /// subject of the group (the HSPs of the subjects shown in its hit too, whose Seq-ids
    /// `Match`). The description table and the alignments then read the sorted group, as NCBI's
    /// `x_DisplayDeflines` and `PruneSeqalign` read the sorted Seq-align set. Nothing is done
    /// (or allocated) when `sort` sorts nothing.
    pub(crate) fn sort(&mut self, sort: super::hit_sort::ReportSort, tblastx: bool) {
        let Some(sorted) = sort_group(
            &self.subject_order,
            &self.subject_hits,
            |hit| super::hit_sort::HspScores::of(hit),
            sort,
            tblastx,
        ) else {
            return;
        };
        for (s_idx, hsp_order) in self.subject_order.iter().zip(&sorted.hsps) {
            if let Some(hits) = self.subject_hits.get_mut(s_idx) {
                *hits = hsp_order.iter().map(|&k| hits[k]).collect();
            }
            if let Some(indices) = self.subject_hit_indices.get_mut(s_idx) {
                *indices = hsp_order.iter().map(|&k| indices[k]).collect();
            }
        }
        self.subject_order = sorted
            .hits
            .iter()
            .map(|&hit| self.subject_order[hit])
            .collect();
    }

    /// The number of the hits of subject `s_idx` that are its own, not those of the subjects
    /// shown in its hit.
    pub(crate) fn own_hsps(&self, s_idx: u32) -> usize {
        self.subject_hits[&s_idx]
            .iter()
            .filter(|hit| hit.hit.s_idx == s_idx)
            .count()
    }

    /// The subjects of the description table (the first `num_descriptions`).
    pub(crate) fn described(&self, num_descriptions: usize) -> &[u32] {
        &self.subject_order[..self.subject_order.len().min(num_descriptions)]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:541-554
    /// The HSPs of psiblast's two-part description table (`x_SplitSeqAlign`), each with the
    /// subject whose hit shows it: the first `num_descriptions` in the group's order, and the first
    /// alone for 0 (the count is checked after an HSP is put in its part).
    pub(crate) fn first_hsps(
        &self,
        num_descriptions: usize,
    ) -> impl Iterator<Item = (u32, &'h PairwiseHit)> + '_ {
        self.subject_order
            .iter()
            .flat_map(move |&shown| {
                self.subject_hits[&shown]
                    .iter()
                    .map(move |&hit| (shown, hit))
            })
            .take(num_descriptions.max(1))
    }

    /// The subjects with alignments (the first `num_alignments`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    /// (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
    pub(crate) fn aligned(&self, num_alignments: usize) -> &[u32] {
        &self.subject_order[..self.subject_order.len().min(num_alignments)]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1523-1530
    /// The subjects of one query's pairwise report in the order its description lines and
    /// alignment headings read them, for the models of the reads that make the texts of those
    /// reads before the rows exist (`common::database::subject_read_texts`): `hits`, each HSP of
    /// the query's hit list with its subject as the report's rows name it, in the order of the
    /// rows, are grouped as [`Self::new`] groups the rows (`joins`) and the groups put in the
    /// order of `-sorthits` (`sort`, TBLASTX's `tblastx`) as [`Self::sort`] puts them, by the
    /// same `sort_query_hits` of the HSPs' values (`scores`, read only where `sort` sorts: the
    /// row's `HspScores::of`). Each subject comes once, those shown in a group's hit after the
    /// group's own subject, in the order of their first HSP.
    pub(crate) fn read_order<T>(
        hits: impl IntoIterator<Item = (u32, T)>,
        scores: impl Fn(&T) -> super::hit_sort::HspScores,
        joins: impl Fn(u32, u32) -> bool,
        sort: super::hit_sort::ReportSort,
        tblastx: bool,
    ) -> Vec<u32> {
        let hits: Vec<(u32, T)> = hits.into_iter().collect();
        let shown = shown_subjects(hits.iter().map(|&(s_idx, _)| s_idx), joins);
        let mut order: Vec<u32> = Vec::new();
        let mut grouped: HashMap<u32, Vec<&T>> = HashMap::new();
        let mut members: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for ((s_idx, hit), shown) in hits.iter().zip(shown) {
            grouped
                .entry(shown)
                .or_insert_with(|| {
                    order.push(shown);
                    Vec::new()
                })
                .push(hit);
            if seen.insert(*s_idx) {
                members.entry(shown).or_default().push(*s_idx);
            }
        }
        if let Some(sorted) = sort_group(&order, &grouped, |hit| scores(hit), sort, tblastx) {
            order = sorted.hits.iter().map(|&group| order[group]).collect();
        }
        order
            .iter()
            .flat_map(|shown| members.remove(shown).unwrap_or_default())
            .collect()
    }
}

/// The subject whose hit shows each of one query's hits, whose subjects are `s_idxs` in the
/// report's order: its own, or, for a subject whose first hit follows a hit of subject `previous`
/// and whose hits `joins(previous, s_idx)` shows in `previous`'s hit, the subject that shows
/// `previous`'s (asked where the subject changes, once per subject). [`QueryGroup::new`] and
/// [`QueryGroup::read_order`] group by it.
fn shown_subjects(
    s_idxs: impl IntoIterator<Item = u32>,
    joins: impl Fn(u32, u32) -> bool,
) -> Vec<u32> {
    let mut shown_in: HashMap<u32, u32> = HashMap::new();
    let mut previous: Option<(u32, u32)> = None;
    s_idxs
        .into_iter()
        .map(|s_idx| {
            let shown = match previous {
                Some((same, shown)) if same == s_idx => shown,
                _ => match shown_in.get(&s_idx) {
                    Some(&shown) => shown,
                    None => {
                        let shown = match previous {
                            Some((before, shown)) if joins(before, s_idx) => shown,
                            _ => s_idx,
                        };
                        shown_in.insert(s_idx, shown);
                        shown
                    }
                },
            };
            previous = Some((s_idx, shown));
            shown
        })
        .collect()
}

/// `hit_sort::sort_query_hits` of a group's subjects `order`, each with its hits `hits` (in the
/// order of the rows) and the values of a hit `scores`: the subjects (indices of `order`) and
/// each one's hits in the order of `-sorthits` and `-sorthsps`, or `None` (nothing read) where
/// `sort` sorts nothing. [`QueryGroup::sort`] and [`QueryGroup::read_order`] sort by it.
fn sort_group<T>(
    order: &[u32],
    hits: &HashMap<u32, Vec<T>>,
    scores: impl Fn(&T) -> super::hit_sort::HspScores,
    sort: super::hit_sort::ReportSort,
    tblastx: bool,
) -> Option<super::hit_sort::SortedHits> {
    if !sort.sorts() {
        return None;
    }
    let scores: Vec<Vec<super::hit_sort::HspScores>> = order
        .iter()
        .map(|s_idx| hits[s_idx].iter().map(&scores).collect())
        .collect();
    super::hit_sort::sort_query_hits(sort, tblastx, &scores)
}

/// The parts of one program's per-query report (`write`).
pub(crate) trait QuerySections: Sync {
    /// The subjects with alignments per query (`QueryGroup::aligned`).
    fn num_alignments(&self) -> usize;

    /// The columns of a line of the report's alignments (`-line_length`, `m_LineLen`), which the
    /// estimates of its windows, chunks and pieces count (`hit_estimate`).
    fn line_length(&self) -> usize;

    /// Whether the hits of subject `s_idx` that follow those of subject `previous` are shown in
    /// `previous`'s hit (`report::db_subject::joins_as_copy`; `QueryGroup::new`).
    fn joins(&self, _previous: u32, _s_idx: u32) -> bool {
        false
    }

    /// `-sorthits` and `-sorthsps` of the report and whether it is TBLASTX's
    /// (`QueryGroup::sort`); none by default.
    fn sort(&self) -> (super::hit_sort::ReportSort, bool) {
        (super::hit_sort::ReportSort::default(), false)
    }

    /// `\n\n`, the query's header, then "No hits found" and the query's footer (`group` is
    /// `None`), or the description table and the blank line after it; `read` before each
    /// description line, where NCBI loads its subject (`QueryWarnings::subject_read`).
    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&QueryGroup<'_>>,
        read: SubjectReadHook<'_, W>,
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

    /// Whether the query's alignments are one section (`alignments`) in place of each subject's
    /// (`subject`): the query-anchored reports, outfmt 1 to 4 (`report::query_anchored`).
    fn merges_alignments(&self) -> bool {
        false
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2002-2009
    /// Where `merges_alignments`: the merge's work on aligned subject `s_idx` of the query
    /// (`group`) after it reads the subject and before it reads the next one, all before the
    /// section's first byte. NCBI makes the Dense-segs of the subject's Std-segs there, which can
    /// throw (BLASTX's `EmptyRowsStop::StdSegLengths`); nothing by default.
    fn merged_subject<W: Write>(
        &self,
        _writer: &mut W,
        _group: &QueryGroup<'_>,
        _s_idx: u32,
    ) -> io::Result<()> {
        Ok(())
    }

    /// The one section of the alignments of the query's aligned subjects (`group`), where
    /// `merges_alignments`.
    fn alignments<W: Write>(
        &self,
        _writer: &mut W,
        _q_idx: usize,
        _group: &QueryGroup<'_>,
    ) -> io::Result<()> {
        Err(io::Error::other(
            "LOSAT internal error: this report has no query-anchored section",
        ))
    }

    /// The footer of a query with hits.
    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()>;
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:286-304 (each app's `Run`)
/// Writes the per-query reports of `hits_by_query` (each query's hits with their indices in
/// the report's hit list; the first is query `first_query`) between a writer's prolog and
/// epilog: on `pool` (`pool`) without a `probe`, else by the serial writer. An error of the
/// writer (NCBI's exception out of `PrintOneResultSet`, a failed write of its stream) ends the
/// run: NCBI leaves its formatting loop for `CATCH_ALL` and loads no subject after it, so no text
/// of a read is written after it (`QueryWarnings::end_run`), also not when the warnings end.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    first_query: usize,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    pool: Option<&SearchPool<'_>>,
) -> io::Result<()> {
    let written = match pool {
        Some(pool) if probe.is_none() => {
            crate::utils::xreach::x_reach!("LOSAT_X_PAIRWISEPAR", true);
            write_parallel(
                writer,
                sections,
                first_query,
                hits_by_query,
                subject_titles,
                warnings.as_deref_mut(),
                pool,
            )
            .map(|_| ())
        }
        _ => {
            if on() {
                crate::utils::xreach::x_reach!("LOSAT_X_PAIRWISEPAR", false);
            }
            write_serial(
                writer,
                sections,
                first_query,
                hits_by_query,
                subject_titles,
                probe,
                warnings.as_deref_mut(),
            )
        }
    };
    if written.is_err() {
        if let Some(warnings) = warnings {
            warnings.end_run();
        }
    }
    written
}

/// The serial writer: each query's warnings and parts on `writer`, in order.
fn write_serial<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    first_query: usize,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    for (q_idx, query_hits) in (first_query..).zip(hits_by_query) {
        // The query's warnings (and those of its batch's reading) come before its preamble
        // (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let group = QueryGroup::for_sections(query_hits, subject_titles, sections);
        // NCBI loads each subject that the query's batch did not prefetch at its description
        // line, or one the table does not show at its alignment's heading
        // (`QueryWarnings::subject_read`).
        let mut read = |s_idx: u32, writer: &mut W| match warnings.as_deref_mut() {
            Some(warnings) => warnings.subject_read(q_idx, s_idx, writer),
            None => Ok(()),
        };
        sections.head(writer, q_idx, group.as_ref(), &mut read)?;
        if let Some(group) = &group {
            if sections.merges_alignments() {
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1997-2003
                // The merge loads every aligned subject before the section's first byte
                // (`merged_subject` after each load); without one (`-num_alignments 0`) there is
                // no section (`DisplaySeqalign` returns at once).
                let aligned = group.aligned(sections.num_alignments());
                if !aligned.is_empty() {
                    for &s_idx in aligned {
                        read(s_idx, &mut *writer)?;
                        sections.merged_subject(writer, group, s_idx)?;
                    }
                    sections.alignments(writer, q_idx, group)?;
                }
            } else {
                for &s_idx in group.aligned(sections.num_alignments()) {
                    read(s_idx, &mut *writer)?;
                    let hsps = group.subject_hits[&s_idx].len();
                    sections.subject(writer, probe.as_deref_mut(), group, s_idx, 0..hsps)?;
                }
            }
            sections.foot(writer, q_idx)?;
        }
    }
    Ok(())
}

/// The estimated bytes of a hit's alignment in a report of `line_length` columns a line
/// (`hit_estimate`). Only the windows', chunks' and pieces' sizes depend on the estimates.
fn estimate(hit: &PairwiseHit, line_length: usize) -> usize {
    hit_estimate(hit.hit.length, line_length)
}

/// The estimated bytes of the alignment of a hit of alignment length `length` in a report of
/// `line_length` columns a line (`alignments_estimate` of one hit).
#[inline]
pub(crate) fn hit_estimate(length: usize, line_length: usize) -> usize {
    alignments_estimate(1, length as u64, line_length)
}

/// The estimated bytes of the alignments of `hits` hits of `columns` alignment columns in all, in
/// a report of `line_length` columns a line (`-line_length`, `QuerySections::line_length`): 400
/// per hit for its header lines, three letters per column (the query row, the middle line and the
/// subject row) and about 60 bytes of labels, coordinates and margins per line of three rows,
/// `columns / line_length` lines (a line's rows of coordinates of width `w` take `35 + 5 * w`
/// bytes besides their letters, 60 at `w` = 5). It holds at every line length: about 63 bytes
/// per column at `-line_length 1`, the 3 letters per column and one line of margins (in the 400)
/// at the largest. At the default 60 it is `400 * hits + 4 * columns`, the estimate before the
/// line length was counted, so the default reports' windows, chunks and pieces are as they were.
/// The query-anchored reports (outfmt 1 to 4) write fewer rows per hit and are estimated by it
/// too. Counted in `u64`, which no alignment length overflows (a 32-bit target's estimate above
/// `usize::MAX` saturates); a line length of 0, which `-line_length` refuses, counts as 1.
#[inline]
fn alignments_estimate(hits: u64, columns: u64, line_length: usize) -> usize {
    let bytes = 400 * hits + 3 * columns + 60 * columns / line_length.max(1) as u64;
    usize::try_from(bytes).unwrap_or(usize::MAX)
}

/// The estimated bytes of a query's header with "No hits found" or the start of its table, and
/// of its footer (the statistics lines).
const QUERY_BYTES: usize = 512;

/// The estimated bytes of one row of a query's description table.
const DESCRIPTION_BYTES: usize = 100;

/// The estimated bytes of a query's report of `line_length` columns a line (windows of queries).
fn query_estimate(query_hits: &[(HspIndex, &PairwiseHit)], line_length: usize) -> usize {
    query_estimate_of(
        query_hits.iter().map(|&(_, hit)| hit.hit.length),
        line_length,
    )
}

/// `query_estimate` of a query whose hits have the alignment lengths `lengths`, in a report of
/// `line_length` columns a line (the report's `-line_length`, as its sections give it,
/// `QuerySections::line_length`): the one estimate of the writer's windows and of the windows of
/// rows that its callers make before it (BLASTN's, BLASTP's and TBLASTX's). The lines are counted
/// over the query's columns in all (`alignments_estimate`), one division per query.
pub(crate) fn query_estimate_of(
    lengths: impl IntoIterator<Item = usize>,
    line_length: usize,
) -> usize {
    let (hits, columns) = lengths
        .into_iter()
        .fold((0u64, 0u64), |(hits, columns), length| {
            (hits + 1, columns + length as u64)
        });
    let heads = QUERY_BYTES as u64 + DESCRIPTION_BYTES as u64 * hits;
    usize::try_from(heads)
        .unwrap_or(usize::MAX)
        .saturating_add(alignments_estimate(hits, columns, line_length))
}

/// The estimated bytes of a window of queries: `WINDOW_BYTES_PER_THREAD` per thread of `pool`
/// (one without a pool).
pub(crate) fn window_bytes(pool: Option<&SearchPool<'_>>) -> usize {
    WINDOW_BYTES_PER_THREAD.saturating_mul(pool.map_or(1, |pool| pool.threads()))
}

/// The windows of the queries of `estimates` (each query's estimated bytes, in order): from each
/// window's first query, the queries up to the one that brings the window to `budget` (at least
/// one query).
pub(crate) fn query_windows(estimates: &[usize], budget: usize) -> Vec<std::ops::Range<usize>> {
    let mut windows = Vec::new();
    let mut first = 0;
    while first < estimates.len() {
        let mut end = first;
        let mut bytes = 0usize;
        while end < estimates.len() && (end == first || bytes < budget) {
            bytes = bytes.saturating_add(estimates[end]);
            end += 1;
        }
        windows.push(first..end);
        first = end;
    }
    windows
}

/// What a piece of a query's report holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Part {
    Head,
    /// The alignments from (position in the aligned subjects, HSP of that subject) up to the
    /// second, exclusive: a run of whole subjects, or part of one subject's HSPs.
    Alignments((usize, usize), (usize, usize)),
    /// The one section of a query-anchored report (`QuerySections::alignments`).
    Section,
    Foot,
}

#[derive(Clone, Copy, Debug)]
struct Piece {
    query: usize,
    part: Part,
    /// Its estimated bytes (chunks of pieces).
    estimate: usize,
}

/// The pieces of the queries `first..` (`groups`, in order): each query's head, its alignments
/// cut after the HSP that brings a run to `piece_bytes` estimated bytes, and its footer.
fn pieces<S: QuerySections>(
    sections: &S,
    first: usize,
    groups: &[Option<QueryGroup<'_>>],
    piece_bytes: usize,
) -> Vec<Piece> {
    let line_length = sections.line_length();
    let mut pieces = Vec::new();
    for (offset, group) in groups.iter().enumerate() {
        let query = first + offset;
        let Some(group) = group else {
            pieces.push(Piece {
                query,
                part: Part::Head,
                estimate: QUERY_BYTES,
            });
            continue;
        };
        pieces.push(Piece {
            query,
            part: Part::Head,
            estimate: QUERY_BYTES + DESCRIPTION_BYTES * group.subject_order.len(),
        });
        let aligned = group.aligned(sections.num_alignments());
        if sections.merges_alignments() {
            if !aligned.is_empty() {
                pieces.push(Piece {
                    query,
                    part: Part::Section,
                    estimate: aligned
                        .iter()
                        .flat_map(|s_idx| &group.subject_hits[s_idx])
                        .fold(0usize, |bytes, hit| {
                            bytes.saturating_add(estimate(hit, line_length))
                        }),
                });
            }
            pieces.push(Piece {
                query,
                part: Part::Foot,
                estimate: QUERY_BYTES,
            });
            continue;
        }
        let (mut start, mut bytes, mut end) = ((0, 0), 0usize, (0, 0));
        for (position, s_idx) in aligned.iter().enumerate() {
            for (hsp, hit) in group.subject_hits[s_idx].iter().enumerate() {
                bytes = bytes.saturating_add(estimate(hit, line_length));
                end = (position, hsp + 1);
                if bytes >= piece_bytes {
                    pieces.push(Piece {
                        query,
                        part: Part::Alignments(start, end),
                        estimate: bytes,
                    });
                    (start, bytes) = (end, 0);
                }
            }
        }
        if start != end {
            pieces.push(Piece {
                query,
                part: Part::Alignments(start, end),
                estimate: bytes,
            });
        }
        pieces.push(Piece {
            query,
            part: Part::Foot,
            estimate: QUERY_BYTES,
        });
    }
    pieces
}

/// One piece into its own buffer, by the serial writer's calls; with the places in the buffer
/// where the serial writer reads a subject (`QueryWarnings::subject_read`; only when `marks`),
/// what a failed call left in the buffer (the serial writer has written it before the error)
/// and the error.
fn format_piece<S: QuerySections>(
    sections: &S,
    group: Option<&QueryGroup<'_>>,
    piece: Piece,
    marks: bool,
) -> Formatted {
    let mut buffer = Vec::new();
    let mut reads: Vec<(usize, u32)> = Vec::new();
    let mut mark = |s_idx: u32, buffer: &mut Vec<u8>| {
        if marks {
            reads.push((buffer.len(), s_idx));
        }
        Ok(())
    };
    let result = (|| match piece.part {
        Part::Head => sections.head(&mut buffer, piece.query, group, &mut mark),
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
                    if from == 0 {
                        mark(s_idx, &mut buffer)?;
                    }
                    sections.subject(&mut buffer, None, group, s_idx, from..to)?;
                }
            }
            Ok(())
        }
        Part::Section => {
            let group = group.expect("a query with a section has hits");
            for &s_idx in group.aligned(sections.num_alignments()) {
                mark(s_idx, &mut buffer)?;
                sections.merged_subject(&mut buffer, group, s_idx)?;
            }
            sections.alignments(&mut buffer, piece.query, group)
        }
        Part::Foot => sections.foot(&mut buffer, piece.query),
    })();
    Formatted {
        buffer,
        reads,
        result,
    }
}

/// A piece formatted into its own buffer (`format_piece`).
struct Formatted {
    buffer: Vec<u8>,
    /// Where the serial writer reads a subject: the buffer's length then and the subject.
    reads: Vec<(usize, u32)>,
    result: io::Result<()>,
}

/// `f` of every item, in order: on `pool`, or in order on this thread without threads (where
/// `pool` gives `None`: the fallback only lets the writer compile there).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn map<T: Sync, R: Send>(
    pool: &SearchPool<'_>,
    items: &[T],
    f: impl Fn(&T) -> R + Sync + Send,
) -> Vec<R> {
    pool.install(|| items.par_iter().map(f).collect())
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn map<T: Sync, R: Send>(
    _pool: &SearchPool<'_>,
    items: &[T],
    f: impl Fn(&T) -> R + Sync + Send,
) -> Vec<R> {
    items.iter().map(f).collect()
}

/// The parallel writer (module header): windows of queries grouped on `pool`, their pieces
/// formatted on `pool` in chunks of at most the budget's estimated bytes (at least one piece),
/// each chunk written in order by the calling thread before the next. Returns the numbers of
/// chunks and pieces.
fn write_parallel<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    first_query: usize,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    warnings: Option<&mut QueryWarnings<'_>>,
    pool: &SearchPool<'_>,
) -> io::Result<(usize, usize)> {
    write_parallel_with(
        writer,
        sections,
        first_query,
        hits_by_query,
        subject_titles,
        warnings,
        pool,
        window_bytes(Some(pool)),
        PIECE_BYTES,
    )
}

/// `write_parallel` with windows and chunks of `budget` and pieces of `piece_bytes` estimated
/// bytes.
#[allow(clippy::too_many_arguments)]
fn write_parallel_with<S: QuerySections, W: Write>(
    writer: &mut W,
    sections: &S,
    first_query: usize,
    hits_by_query: &[Vec<(HspIndex, &PairwiseHit)>],
    subject_titles: &[Arc<[u8]>],
    mut warnings: Option<&mut QueryWarnings<'_>>,
    pool: &SearchPool<'_>,
    budget: usize,
    piece_bytes: usize,
) -> io::Result<(usize, usize)> {
    let (mut chunks, mut pieces_written) = (0, 0);
    let estimates: Vec<usize> = hits_by_query
        .iter()
        .map(|query_hits| query_estimate(query_hits, sections.line_length()))
        .collect();
    for window in query_windows(&estimates, budget) {
        let groups = map(pool, &hits_by_query[window.clone()], |query_hits| {
            QueryGroup::for_sections(query_hits, subject_titles, sections)
        });
        let first = first_query + window.start;
        let planned = pieces(sections, first, &groups, piece_bytes);
        let mut at = 0;
        while at < planned.len() {
            let mut to = at;
            let mut bytes = 0usize;
            while to < planned.len() && (to == at || bytes < budget) {
                bytes = bytes.saturating_add(planned[to].estimate);
                to += 1;
            }
            let chunk = &planned[at..to];
            // The pieces of a query whose reads write texts mark where the serial writer reads
            // its subjects; the texts are written there (`QueryWarnings::subject_read`).
            let marked: Vec<bool> = chunk
                .iter()
                .map(|piece| {
                    warnings
                        .as_deref()
                        .is_some_and(|warnings| warnings.reads_with_texts(piece.query))
                })
                .collect();
            let indexed: Vec<(Piece, bool)> = chunk.iter().copied().zip(marked).collect();
            let buffers = map(pool, &indexed, |&(piece, marks)| {
                format_piece(sections, groups[piece.query - first].as_ref(), piece, marks)
            });
            for (piece, formatted) in chunk.iter().zip(buffers) {
                if piece.part == Part::Head {
                    if let Some(warnings) = warnings.as_deref_mut() {
                        warnings.before_query(piece.query, &mut *writer)?;
                    }
                }
                super::query_warnings::write_marked(
                    &formatted.buffer,
                    &formatted.reads,
                    piece.query,
                    &mut *writer,
                    warnings.as_deref_mut(),
                )?;
                formatted.result?;
            }
            chunks += 1;
            at = to;
        }
        pieces_written += planned.len();
    }
    Ok((chunks, pieces_written))
}

#[cfg(test)]
#[path = "x_pairwise_par_tests.rs"]
mod tests;

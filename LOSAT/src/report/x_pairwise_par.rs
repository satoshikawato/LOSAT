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
        })
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
mod tests {
    use super::*;
    use crate::common::Hit;

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

        fn head<W: Write>(
            &self,
            writer: &mut W,
            q_idx: usize,
            group: Option<&QueryGroup<'_>>,
            read: SubjectReadHook<'_, W>,
        ) -> io::Result<()> {
            let subjects = group.map_or(0, |group| group.subject_order.len());
            writeln!(writer, "head {q_idx} subjects {subjects}")?;
            // One line per subject, each read before it (as a description line).
            for &s_idx in group.map_or(&[][..], |group| group.subject_order.as_slice()) {
                read(s_idx, writer)?;
                writeln!(writer, "line {s_idx}")?;
            }
            Ok(())
        }

        fn subject<W: Write>(
            &self,
            writer: &mut W,
            _probe: Option<&mut FormatProbe<'_>>,
            group: &QueryGroup<'_>,
            s_idx: u32,
            part: std::ops::Range<usize>,
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

        fn merges_alignments(&self) -> bool {
            self.merged
        }

        fn alignments<W: Write>(
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

        fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
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
                    let formatted =
                        format_piece(&marks, groups[piece.query].as_ref(), piece, false);
                    formatted.result.unwrap();
                    assert!(formatted.reads.is_empty());
                    joined.extend_from_slice(&formatted.buffer);
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

    // Review 3 of line RG T3a, D1: the subjects whose description lines a report shows, each of
    // which NCBI looks up again by its display id (`pairwise::blastp_described_subjects`): the
    // first `-num_descriptions` subjects of the group, or in psiblast's two-part table the
    // subjects of its first `-num_descriptions` HSPs, at least the first (`first_hsps`,
    // blast_format.cpp:541-554); a subject shown in another's hit is that one's.
    #[test]
    fn the_described_subjects_are_the_first_subjects_or_those_of_the_first_hsps() {
        let hits = report_hits();
        let by_query = by_query(&hits);
        let titles = titles();
        let group = QueryGroup::new(&by_query[0], &titles, |_, _| false).unwrap();
        let first = |count| -> Vec<u32> {
            group
                .first_hsps(count)
                .map(|(shown, hit)| {
                    assert_eq!(shown, hit.hit.s_idx);
                    shown
                })
                .collect()
        };
        assert_eq!(first(0), [5]);
        assert_eq!(first(1), [5]);
        assert_eq!(first(4), [5, 5, 5, 2]);
        assert_eq!(first(100).len(), 8);
        let described = |q_idx, two_part, count, joins: fn(u32, u32) -> bool| {
            crate::report::pairwise::blastp_described_subjects(
                Default::default(),
                two_part,
                count,
                q_idx,
                &hits,
                joins,
            )
        };
        let none: fn(u32, u32) -> bool = |_, _| false;
        assert!(described(0, false, 0, none).is_empty());
        assert_eq!(described(0, false, 2, none), [5, 2]);
        assert_eq!(described(0, false, 500, none), [5, 2, 7]);
        assert_eq!(described(0, true, 0, none), [5]);
        assert_eq!(described(0, true, 4, none), [5, 2]);
        assert_eq!(described(0, true, 5, none), [5, 2, 7]);
        assert_eq!(
            described(0, false, 2, |previous, s_idx| (previous, s_idx) == (5, 2)),
            [5, 7]
        );
        assert!(described(1, false, 500, none).is_empty());
        assert_eq!(described(2, true, 6, none), [1]);
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
}

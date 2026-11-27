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
use std::io;
use std::sync::{Arc, OnceLock};

use super::pairwise::PairwiseHit;
use super::query_warnings::{QueryWarnings, ReportHook, ReportPlace};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::report::ostream::{DisplaySeqalign, ReportWrite};
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

/// How many alignment headings NCBI's display starts for one subject of a query
/// (`QueryGroup::heading_rank`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum HeadingStarts {
    /// One heading, for the subject's first alignment.
    Once,
    /// One for each alignment: the heading's title cannot be made, so each alignment tries it.
    EachAlignment,
    /// None: the subject's load throws before its heading is started.
    Never,
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1929, 1955-1981, 2307
    /// The rank of subject `s_idx`'s alignment heading (NCBI's `m_cur_align`, the `blast_rank` of
    /// its `-html` links, from 1 for each query): `x_PrintDefLine` counts each heading it starts
    /// (`starts`): once for a subject, once for each alignment of a subject whose heading throws
    /// (its title cannot be made; NCBI tries again at each alignment, as `write_hsps_or_skips`
    /// notes), and never for a subject whose load throws (the display's load of the alignment's
    /// subject fails before the heading is started; line RG T4a2 review N1). The ranks are made
    /// once, at the first call (only `-html` headings call it); a heading that the registry fails
    /// throws in every subject, so it shows no rank.
    pub(crate) fn heading_rank(&self, s_idx: u32, starts: impl Fn(u32) -> HeadingStarts) -> usize {
        let ranks = self.heading_ranks.get_or_init(|| {
            let mut started = 0;
            self.subject_order
                .iter()
                .map(|&subject| {
                    let rank = started + 1;
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1950-1982,3965-3972
                    // An HSP whose query row starts left of its record throws in `x_PrepareRowData`
                    // before `x_PrintDefLine`, so it neither starts a heading nor increments `m_cur_align`.
                    let shown_hsps = self.subject_hits[&subject]
                        .iter()
                        .filter(|hit| !hit.query_left_of_record())
                        .count();
                    started += if shown_hsps == 0 {
                        0
                    } else {
                        match starts(subject) {
                            HeadingStarts::Once => 1,
                            HeadingStarts::EachAlignment => shown_hsps,
                            HeadingStarts::Never => 0,
                        }
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

    /// The subjects of the description table (the first `num_descriptions`; psiblast's two-part
    /// table after the first round shows those of `described_subjects`).
    pub(crate) fn described(&self, num_descriptions: usize) -> &[u32] {
        first(&self.subject_order, num_descriptions)
    }

    /// Each HSP of the group in the group's order, with the subject whose hit shows it.
    pub(crate) fn hsps(&self) -> impl Iterator<Item = (u32, &'h PairwiseHit)> + '_ {
        self.subject_order.iter().flat_map(move |&shown| {
            self.subject_hits[&shown]
                .iter()
                .map(move |&hit| (shown, hit))
        })
    }

    /// The subjects whose description lines the group's report shows, in the order of the lines
    /// (`described_subjects` of the group's subjects and their HSPs' own subjects).
    pub(crate) fn described_subjects(
        &self,
        num_descriptions: usize,
        found_before: Option<&[bool]>,
    ) -> Vec<u32> {
        described_subjects(
            &self.subject_order,
            |shown| self.subject_hits[&shown].iter().map(|hit| hit.hit.s_idx),
            num_descriptions,
            found_before,
        )
    }

    /// The subjects with alignments (the first `num_alignments`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    /// (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
    pub(crate) fn aligned(&self, num_alignments: usize) -> &[u32] {
        first(&self.subject_order, num_alignments)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1523-1530
    /// What one query's pairwise report reads, for the models of the reads that are made before
    /// the rows exist (the loads, `report::db_subject::HitLoads::of`; the texts of the reads,
    /// `common::database::subject_read_texts`): `hits`, each HSP of the query's hit list with its
    /// subject as the report's rows name it, in the order of the rows, are grouped as
    /// [`Self::new`] groups the rows (`joins`), and the groups and each group's HSPs put in the
    /// order of `-sorthits` and `-sorthsps` (`reading.sort`, TBLASTX's `tblastx`) as
    /// [`Self::sort`] puts them, by the same `sort_query_hits` of the HSPs' values (`scores`, read
    /// only where the sort sorts: the row's `HspScores::of`); of those hits, the description
    /// lines show those of `described_subjects` (the counts of `reading.shown`; `found_before`:
    /// psiblast's previous round's set, its two-part table), the alignments the first ones
    /// (`aligned`): the writer's own answers on the same groups.
    pub(crate) fn read_order<T>(
        hits: impl IntoIterator<Item = (u32, T)>,
        scores: impl Fn(&T) -> super::hit_sort::HspScores,
        joins: impl Fn(u32, u32) -> bool,
        reading: super::pairwise::PairwiseReading,
        tblastx: bool,
        found_before: Option<&[bool]>,
    ) -> super::db_subject::PairwiseReads {
        let hits: Vec<(u32, T)> = hits.into_iter().collect();
        let shown = shown_subjects(hits.iter().map(|&(s_idx, _)| s_idx), joins);
        let mut order: Vec<u32> = Vec::new();
        // Each group's HSPs, with their own subjects, and its subjects, each once.
        let mut grouped: HashMap<u32, Vec<(u32, &T)>> = HashMap::new();
        let mut members: HashMap<u32, Vec<usize>> = HashMap::new();
        let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for ((s_idx, hit), shown) in hits.iter().zip(shown) {
            grouped
                .entry(shown)
                .or_insert_with(|| {
                    order.push(shown);
                    Vec::new()
                })
                .push((*s_idx, hit));
            if seen.insert(*s_idx) {
                members.entry(shown).or_default().push(*s_idx as usize);
            }
        }
        if let Some(sorted) = sort_group(
            &order,
            &grouped,
            |(_, hit)| scores(hit),
            reading.sort,
            tblastx,
        ) {
            for (shown, hsp_order) in order.iter().zip(&sorted.hsps) {
                if let Some(hsps) = grouped.get_mut(shown) {
                    *hsps = hsp_order.iter().map(|&k| hsps[k]).collect();
                }
            }
            order = sorted.hits.iter().map(|&group| order[group]).collect();
        }
        let described = described_subjects(
            &order,
            |shown| grouped[&shown].iter().map(|&(own, _)| own),
            reading.shown.descriptions,
            found_before,
        );
        let position: HashMap<u32, usize> = order
            .iter()
            .enumerate()
            .map(|(hit, &shown)| (shown, hit))
            .collect();
        super::db_subject::PairwiseReads {
            described: described.iter().map(|shown| position[shown]).collect(),
            aligned: first(&order, reading.shown.alignments).len(),
            hits: order
                .iter()
                .map(|shown| members.remove(shown).unwrap_or_default())
                .collect(),
        }
    }
}

/// The writer's reads of the rows of one query (`rows`, all of the query `q_idx`, in the order of
/// the rows; `joins` and `reading` as the writer's sections give them, `tblastx` TBLASTX's sort):
/// its group (`QueryGroup::sorted`), the subjects shown in each hit (its own first, then the
/// copies, each once, in the order of their first row: the writer reads a hit's copies with the
/// hit and has no order of them, so that order is the model's rule, `PairwiseReads::hits`), the
/// hits of its description lines
/// (`QueryGroup::described_subjects`) and of its alignments (`QueryGroup::aligned`). The tests of
/// each program's `pairwise_read_order` compare their model with it (review 2 of line DBID
/// pfoff-3, L5): the groups, their order and what the description lines and alignments show. The
/// order of a hit's copies here is the model's rule copied, not the writer's, so those tests do
/// not check it (review 3, L-4).
#[cfg(test)]
pub(crate) fn writer_reads(
    rows: &[PairwiseHit],
    joins: impl Fn(u32, u32) -> bool,
    reading: super::pairwise::PairwiseReading,
    tblastx: bool,
    found_before: Option<&[bool]>,
) -> super::db_subject::PairwiseReads {
    let query_hits: Vec<(HspIndex, &PairwiseHit)> = rows.iter().enumerate().collect();
    let titles: Vec<Arc<[u8]>> = Vec::new();
    let Some(group) = QueryGroup::sorted(&query_hits, &titles, joins, reading.sort, tblastx) else {
        return Default::default();
    };
    let position = |shown: &u32| group.subject_order.iter().position(|s| s == shown).unwrap();
    let shown_in: HashMap<u32, u32> = group
        .hsps()
        .map(|(shown, hit)| (hit.hit.s_idx, shown))
        .collect();
    let mut members: HashMap<u32, Vec<usize>> = HashMap::new();
    for row in rows {
        let hit = members.entry(shown_in[&row.hit.s_idx]).or_default();
        if !hit.contains(&(row.hit.s_idx as usize)) {
            hit.push(row.hit.s_idx as usize);
        }
    }
    super::db_subject::PairwiseReads {
        hits: group
            .subject_order
            .iter()
            .map(|shown| members.remove(shown).unwrap_or_default())
            .collect(),
        described: group
            .described_subjects(reading.shown.descriptions, found_before)
            .iter()
            .map(position)
            .collect(),
        aligned: group.aligned(reading.shown.alignments).len(),
    }
}

/// The first `count` subjects of `order` (`-num_descriptions`, `-num_alignments`).
fn first(order: &[u32], count: usize) -> &[u32] {
    &order[..order.len().min(count)]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:541-554
// `x_SplitSeqAlign` walks the query's Seq-aligns (one per HSP) in order, puts each in the part of
// the subjects of the previous round's set or in the new part by its subject's id, and stops once
// it has put `-num_descriptions` of them, counting after it puts one.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:578-602
// `x_DisplayDeflines` then shows the part of the subjects found before, a blank line, and the new
// part, each a description table of its own.
/// psiblast's two-part description table of one query after the first round (`x_SplitSeqAlign`):
/// of the query's HSPs `hsps` (Seq-aligns, in the report's order, each with the subject whose hit
/// shows it), the first `num_descriptions`, counted after an HSP is put in its part (so the count
/// is of HSPs, not hits, and `-num_descriptions 0` keeps the first), each in the part of the
/// subjects found before where `found_before` holds for it, else in the new part. The lines of a
/// part are its HSPs of one shown subject that follow each other in it, each line the subject and
/// those HSPs. Returns the part of the subjects found before, then the new part: the table shows
/// them in that order (`psiblast_round`).
pub(crate) fn two_part_lines<T>(
    hsps: impl IntoIterator<Item = (u32, T)>,
    num_descriptions: usize,
    found_before: impl Fn(&T) -> bool,
) -> [Vec<(u32, Vec<T>)>; 2] {
    let mut parts: [Vec<(u32, Vec<T>)>; 2] = [Vec::new(), Vec::new()];
    for (shown, hsp) in hsps.into_iter().take(num_descriptions.max(1)) {
        let part = &mut parts[usize::from(!found_before(&hsp))];
        match part.last_mut() {
            Some((last, hsps)) if *last == shown => hsps.push(hsp),
            _ => part.push((shown, vec![hsp])),
        }
    }
    parts
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:686-700
// The description table makes a line of each subject of the Seq-aligns in their order, up to
// `-num_descriptions` lines.
/// The subjects whose description lines one query's pairwise report shows, each once, in the order
/// of the lines: of the report's subjects `order` (in its order, each the subject whose hit shows
/// its HSPs, `hsps(shown)` the own subjects of those HSPs in the hit's order), the first
/// `num_descriptions` (`QueryGroup::described`), or after psiblast's first round (`found_before`:
/// whether each subject, by subject index, is in the previous round's inclusion set) those of the
/// lines of its two-part table (`two_part_lines`, which count HSPs). The one answer of which
/// subjects a report's description lines show: for its writer (`psiblast_round`), for the
/// refusals decided before it (`pairwise::blastp_described_subjects`) and for the models of its
/// loads (`QueryGroup::read_order`).
pub(crate) fn described_subjects<I: IntoIterator<Item = u32>>(
    order: &[u32],
    hsps: impl Fn(u32) -> I,
    num_descriptions: usize,
    found_before: Option<&[bool]>,
) -> Vec<u32> {
    let Some(found_before) = found_before else {
        return first(order, num_descriptions).to_vec();
    };
    let parts = two_part_lines(
        order
            .iter()
            .flat_map(|&shown| hsps(shown).into_iter().map(move |own| (shown, own))),
        num_descriptions,
        |&own| found_before.get(own as usize).copied().unwrap_or(false),
    );
    let mut described: Vec<u32> = Vec::new();
    for (shown, _) in parts.iter().flatten() {
        if !described.contains(shown) {
            described.push(*shown);
        }
    }
    described
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
    /// (`QueryGroup::sort`). Every program gives its own: the order of the hits, which the reads
    /// of a batch whose prefetch is off follow (`read_order`), is the report's (review 3 of line
    /// DBID pfoff-3, L-1).
    fn sort(&self) -> (super::hit_sort::ReportSort, bool);

    /// `\n\n`, the query's header, then "No hits found" and the query's footer (`group` is
    /// `None`), or the description table and the blank line after it; `place` at each
    /// `ReportPlace` of it: before each description line, where NCBI loads its subject
    /// (`ReportPlace::Load`), and with `-html` the reads of NCBI's formatter registries
    /// (`ReportPlace::Html`).
    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&QueryGroup<'_>>,
        place: ReportHook<'_, W>,
    ) -> io::Result<()>;

    /// The alignments of the HSPs `part` of subject `s_idx` of the query (`group`), the
    /// subject's heading with its first HSP: the parts of `0..` its HSPs in order write the
    /// subject's alignments; `place` at the `-html` reads of NCBI's formatter registries that the
    /// heading's links make (`ReportPlace::Html`, before the heading's bytes).
    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: ReportHook<'_, W>,
    ) -> io::Result<()>;

    /// The run's `-html` links (`PairwiseConfig::html_links`), whose reads of NCBI's formatter
    /// registries write their posts (`HtmlLinks::posts_of`); `None` (the default) for the plain
    /// report.
    fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        None
    }

    /// Whether subject `s_idx`'s Seq-align id is a local id (`ParsedSubject::alignment_id_is_local`),
    /// which `GetDbType` does not load before the `-html` alignments read `./.ncbirc`
    /// (`alignments_start`); `false` by default.
    fn alignment_id_is_local(&self, _s_idx: u32) -> bool {
        false
    }

    /// The report's view (`PairwiseView`), which says what the report does after a subject whose
    /// lookup throws (`alignment_end`): every program gives its own (review 3 of line DBID
    /// pfoff-3, L-1).
    fn view(&self) -> super::pairwise::PairwiseView;

    /// Whether the query's alignments are one section (`alignments`) in place of each subject's
    /// (`subject`): the query-anchored reports, outfmt 1 to 4 (`report::query_anchored`).
    fn merges_alignments(&self) -> bool {
        self.view() != super::pairwise::PairwiseView::Pairwise
    }

    /// The load of subject `s_idx` where it throws (its batch's prefetch is off and its lookup
    /// fails, `report::db_subject::HitLoads::thrown`): its exception, and whether its Seq-align id
    /// is a local id (`pairwise::thrown_load`); `None` by default (every subject loads). Its
    /// description line and alignments are written as NCBI's `catch`es write them, and the
    /// query's alignments may end the run there (`alignment_end`).
    fn thrown(&self, _s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        None
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2002-2009
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1817-1824
    /// NCBI's Dense-segs of the Std-segs of aligned subject `s_idx`'s HSPs (`group`), made before
    /// any byte of the query's alignments, which can throw (BLASTX's
    /// `EmptyRowsStop::StdSegLengths`): where `merges_alignments`, the merge's work on the subject
    /// after it reads it and before it reads the next one; in the `-html` pairwise report, for
    /// the aligned subjects that `html_links_densegs` does not pass over, before the first
    /// alignment. Nothing by default.
    fn subject_densegs<W: crate::report::ostream::ReportWrite>(
        &self,
        _writer: &mut W,
        _group: &QueryGroup<'_>,
        _s_idx: u32,
    ) -> io::Result<()> {
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4590-4610
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2145-2154
    /// The Seq-id text (`CSeq_id::GetSeqIdString`) of aligned subject `s_idx`'s Seq-align id
    /// without and with its version, which the `-html` pairwise report compares with the last
    /// subject whose Dense-segs it made (`html_links_densegs`). `None` (the default) where the
    /// report's Dense-segs never throw (`subject_densegs` makes nothing there): then that pass
    /// makes none for the subject.
    fn alignment_id_texts(&self, _s_idx: u32) -> Option<[Vec<u8>; 2]> {
        None
    }

    /// The one section of the alignments of the query's aligned subjects (`group`), where
    /// `merges_alignments`.
    fn alignments<W: crate::report::ostream::ReportWrite>(
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
    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()>;
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:286-304 (each app's `Run`)
/// Writes the per-query reports of `hits_by_query` (each query's hits with their indices in
/// the report's hit list; the first is query `first_query`) between a writer's prolog and
/// epilog: on `pool` (`pool`) without a `probe`, else by the serial writer. An error of the
/// writer (NCBI's exception out of `PrintOneResultSet`, a failed write of its stream) ends the
/// run: NCBI leaves its formatting loop for `CATCH_ALL` and loads no subject after it, so no text
/// of a read is written after it (`QueryWarnings::end_run`), also not when the warnings end.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write<S: QuerySections, W: ReportWrite>(
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

/// Where the alignments of the query of `group` end the run (`aligned_end` of what the report's
/// view does after a throw, with `-html` where the report has its links,
/// `PairwiseView::after_throw`, the answer that the model of the loads takes too,
/// `FormatThrow::of`, and of the loads of its aligned subjects, `QuerySections::thrown`).
pub(crate) fn alignment_end<'s, S: QuerySections>(
    sections: &'s S,
    group: &QueryGroup<'_>,
) -> Option<(usize, &'s crate::ncbi_exception::NcbiException)> {
    aligned_end(
        sections.view().after_throw(sections.html_links().is_some()),
        group.aligned(sections.num_alignments()),
        |s_idx| sections.thrown(s_idx),
    )
}

/// Where the alignments of the subjects `aligned` (in their order) of a report whose writer does
/// `after` after a throw end the run (`AfterThrow::alignment_end` of the loads `thrown`: the
/// exception of a subject whose load throws, and whether its Seq-align id is local): the position
/// among the aligned subjects and the exception, which the writers return there (`end_here`)
/// after the loads before it, as NCBI's application ends with it (`CATCH_ALL`).
pub(super) fn aligned_end<'t>(
    after: crate::report::db_subject::AfterThrow,
    aligned: &[u32],
    thrown: impl Fn(u32) -> Option<(&'t crate::ncbi_exception::NcbiException, bool)>,
) -> Option<(usize, &'t crate::ncbi_exception::NcbiException)> {
    let end = after.alignment_end(
        aligned.iter().copied(),
        |s_idx| thrown(s_idx).is_some(),
        |s_idx| thrown(s_idx).is_some_and(|(_, local)| local),
    )?;
    Some((end, thrown(aligned[end])?.0))
}

/// The exception of `end` (`alignment_end`) as the writer's error where the alignments reach its
/// position `position` (before the load of that aligned subject).
pub(super) fn end_here(
    end: Option<(usize, &crate::ncbi_exception::NcbiException)>,
    position: usize,
) -> io::Result<()> {
    match end {
        Some((at, exception)) if at == position => Err(exception.clone().into_io_error()),
        _ => Ok(()),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1881-1885
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1907-1914
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1921-1934
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4577-4613
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4188-4212
/// The start of the alignment of aligned subject `s_idx` at position `position` of the query's
/// alignments (`group`), at its `place`s, before its heading: the end of the run where `end`
/// (`alignment_end`) is at that position, then its load (`ReportPlace::Load`). Before the first
/// alignment, with `-html`, `x_InitAlignParams` reads `./.ncbirc` after `GetDbType` loads the
/// first alignment's subject, so before that load (and before the end of the run there) where
/// its Seq-align id is local (`GetDbType` does not load it; the merge of a query-anchored
/// section, outfmt 1 to 4, or the pairwise heading does); then the pairwise report (outfmt 0)
/// makes the Dense-segs of the HSPs of the aligned subjects that `html_links_densegs` does not
/// pass over, for its links (`x_PreProcessSeqAlign`, `QuerySections::subject_densegs`), which
/// can throw, before any alignment. One order for both writers (line RG T4b review a B1, N3).
pub(crate) fn alignment_start<S: QuerySections, W: ReportWrite>(
    sections: &S,
    group: &QueryGroup<'_>,
    end: Option<(usize, &crate::ncbi_exception::NcbiException)>,
    position: usize,
    s_idx: u32,
    writer: &mut W,
    place: ReportHook<'_, W>,
) -> io::Result<()> {
    if position > 0 || sections.html_links().is_none() {
        end_here(end, position)?;
        return place(ReportPlace::Load(s_idx), writer);
    }
    let ncbirc = ReportPlace::Html(super::html_ids::HtmlPost::Ncbirc);
    if sections.alignment_id_is_local(s_idx) {
        place(ncbirc, writer)?;
        html_links_densegs(sections, group, writer)?;
        end_here(end, position)?;
        place(ReportPlace::Load(s_idx), writer)
    } else {
        end_here(end, position)?;
        place(ReportPlace::Load(s_idx), writer)?;
        place(ncbirc, writer)?;
        html_links_densegs(sections, group, writer)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4590-4610
/// The Dense-segs that the `-html` pairwise report (outfmt 0) makes of the HSPs of the aligned
/// subjects of the query (`group`), in their order, before its first alignment, for the
/// alignments' links (`x_PreProcessSeqAlign`, `x_CalcUrlLinksParams`; `alignment_start`). A
/// subject is passed over where both texts of the last subject whose Dense-segs were made are
/// not empty, its text without the version is that subject's and its text with the version is
/// not (`QuerySections::alignment_id_texts`: another version of one accession); the passed
/// subject does not become the last one, and its Dense-segs are made at its own alignment, as in
/// the plain report (line RG T4b review b F1). The plain report makes each at its alignment, the
/// query-anchored section at its merge.
fn html_links_densegs<S: QuerySections, W: ReportWrite>(
    sections: &S,
    group: &QueryGroup<'_>,
    writer: &mut W,
) -> io::Result<()> {
    if sections.merges_alignments() {
        return Ok(());
    }
    let mut last: [Vec<u8>; 2] = Default::default();
    for &s_idx in group.aligned(sections.num_alignments()) {
        let Some(texts) = sections.alignment_id_texts(s_idx) else {
            continue;
        };
        let [id, versioned] = &last;
        if !id.is_empty() && !versioned.is_empty() && *id == texts[0] && *versioned != texts[1] {
            continue;
        }
        sections.subject_densegs(writer, group, s_idx)?;
        last = texts;
    }
    Ok(())
}

/// What the report of query `query` writes to the run's diagnostics at `place`, after the report
/// written so far (`writer`, flushed first): a load's text (`QueryWarnings::subject_read`) or what
/// a `-html` read of NCBI's formatter registries posts (`html_at_place`). Inlined at each hook of
/// the writers: the plain report's places are the loads that it reads alone (line RG T4a2 fix
/// round 2, review c C5: an out-of-line call at each load of a plain run).
#[inline]
pub(crate) fn at_place<W: io::Write>(
    warnings: Option<&mut QueryWarnings<'_>>,
    links: Option<&super::html_ids::HtmlLinks>,
    query: usize,
    place: ReportPlace,
    writer: &mut W,
) -> io::Result<()> {
    match place {
        ReportPlace::Load(s_idx) => match warnings {
            Some(warnings) => warnings.subject_read(query, s_idx, writer),
            None => Ok(()),
        },
        ReportPlace::Html(post) => html_at_place(warnings, links, post, writer),
    }
}

/// `at_place` of a `-html` read `post`: what it writes (`HtmlLinks::posts_of` of the run's `links`
/// for the report of `warnings`, `QueryWarnings::report`), the posts after the report written so
/// far and the `GETURL_DEBUG` line. One that posts without warnings is LOSAT's internal error
/// (`query_warnings::HTML_POSTS_WITHOUT_WARNINGS`).
fn html_at_place<W: io::Write>(
    warnings: Option<&mut QueryWarnings<'_>>,
    links: Option<&super::html_ids::HtmlLinks>,
    post: super::html_ids::HtmlPost,
    writer: &mut W,
) -> io::Result<()> {
    let Some(links) = links else {
        return Ok(());
    };
    let report = warnings.as_deref().map_or(0, QueryWarnings::report);
    let (posts, line) = links.posts_of(post, report);
    if posts.is_empty() && line.is_empty() {
        return Ok(());
    }
    let Some(warnings) = warnings else {
        return Err(io::Error::other(
            super::query_warnings::HTML_POSTS_WITHOUT_WARNINGS,
        ));
    };
    if !posts.is_empty() {
        warnings.html_post(&posts, writer)?;
    }
    if !line.is_empty() {
        warnings.html_stderr_line(&line, writer)?;
    }
    Ok(())
}

/// The serial writer: each query's warnings and parts on `writer`, in order.
fn write_serial<S: QuerySections, W: ReportWrite>(
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
        // (`QueryWarnings::subject_read`); with `-html` its formatter reads its registries
        // (`ReportPlace::Html`).
        let links = sections.html_links();
        let mut place = |place: ReportPlace, writer: &mut W| {
            at_place(warnings.as_deref_mut(), links, q_idx, place, writer)
        };
        sections.head(writer, q_idx, group.as_ref(), &mut place)?;
        if let Some(group) = &group {
            let mut display = DisplaySeqalign::default();
            if sections.merges_alignments() {
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1997-2003
                // The merge loads every aligned subject before the section's first byte
                // (`subject_densegs` after each load); without one (`-num_alignments 0`) there is
                // no section (`DisplaySeqalign` returns at once). A load that throws ends the run
                // there (`alignment_end`). With `-html` the section reads `./.ncbirc` where the
                // pairwise alignments do, at the first load (`alignment_start`).
                let aligned = group.aligned(sections.num_alignments());
                let end = alignment_end(sections, group);
                if !aligned.is_empty() {
                    for (position, &s_idx) in aligned.iter().enumerate() {
                        alignment_start(
                            sections,
                            group,
                            end,
                            position,
                            s_idx,
                            &mut *writer,
                            &mut place,
                        )?;
                        sections.subject_densegs(writer, group, s_idx)?;
                    }
                    display.enter(writer)?;
                    sections.alignments(writer, q_idx, group)?;
                }
            } else {
                // `GetDbType`'s load of the first alignment's subject may end the run before any
                // alignment (`alignment_end`).
                let end = alignment_end(sections, group);
                for (position, &s_idx) in
                    group.aligned(sections.num_alignments()).iter().enumerate()
                {
                    alignment_start(
                        sections,
                        group,
                        end,
                        position,
                        s_idx,
                        &mut *writer,
                        &mut place,
                    )?;
                    display.enter(writer)?;
                    let hsps = group.subject_hits[&s_idx].len();
                    sections.subject(
                        writer,
                        probe.as_deref_mut(),
                        group,
                        s_idx,
                        0..hsps,
                        &mut place,
                    )?;
                }
            }
            display.leave(writer)?;
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
/// where the serial writer writes to the diagnostics (`ReportPlace`, `at_place`; only when
/// `marks`), what a failed call left in the buffer (the serial writer has written it before the
/// error) and the error.
fn format_piece<S: QuerySections>(
    sections: &S,
    group: Option<&QueryGroup<'_>>,
    piece: Piece,
    marks: bool,
) -> Formatted {
    let mut buffer = crate::report::ostream::CutBuffer::default();
    let mut places: Vec<(usize, ReportPlace)> = Vec::new();
    let mut mark = |place: ReportPlace, buffer: &mut crate::report::ostream::CutBuffer| {
        if marks {
            places.push((buffer.bytes.len(), place));
        }
        Ok(())
    };
    let result = (|| match piece.part {
        Part::Head => sections.head(&mut buffer, piece.query, group, &mut mark),
        Part::Alignments((first, first_hsp), (last, end_hsp)) => {
            let group = group.expect("a query with alignment pieces has hits");
            let aligned = group.aligned(sections.num_alignments());
            let end = alignment_end(sections, group);
            for (position, &s_idx) in aligned.iter().enumerate().take(last + 1).skip(first) {
                let from = if position == first { first_hsp } else { 0 };
                let to = if position == last {
                    end_hsp
                } else {
                    group.subject_hits[&s_idx].len()
                };
                if from < to {
                    if from == 0 {
                        alignment_start(
                            sections,
                            group,
                            end,
                            position,
                            s_idx,
                            &mut buffer,
                            &mut mark,
                        )?;
                    }
                    sections.subject(&mut buffer, None, group, s_idx, from..to, &mut mark)?;
                }
            }
            Ok(())
        }
        Part::Section => {
            let group = group.expect("a query with a section has hits");
            let end = alignment_end(sections, group);
            for (position, &s_idx) in group.aligned(sections.num_alignments()).iter().enumerate() {
                alignment_start(
                    sections,
                    group,
                    end,
                    position,
                    s_idx,
                    &mut buffer,
                    &mut mark,
                )?;
                sections.subject_densegs(&mut buffer, group, s_idx)?;
            }
            sections.alignments(&mut buffer, piece.query, group)
        }
        Part::Foot => sections.foot(&mut buffer, piece.query),
    })();
    Formatted {
        buffer,
        places,
        result,
    }
}

/// A piece formatted into its own buffer (`format_piece`), with NCBI's flushes made in it (a
/// query footer's).
struct Formatted {
    buffer: crate::report::ostream::CutBuffer,
    /// Where the serial writer writes to the diagnostics: the buffer's length then and the place.
    places: Vec<(usize, ReportPlace)>,
    result: io::Result<()>,
}

impl Formatted {
    /// Writes the piece of query `query` to `writer`, with NCBI's insertions, regions and flushes
    /// of its buffer (`CutCursor`), and what the serial writer writes to the diagnostics at its
    /// places, after the bytes before each (`at_place`).
    fn write<W: ReportWrite>(
        &self,
        query: usize,
        writer: &mut W,
        mut warnings: Option<&mut QueryWarnings<'_>>,
        links: Option<&super::html_ids::HtmlLinks>,
    ) -> io::Result<()> {
        let report = &self.buffer.bytes;
        let mut cuts = self.buffer.cuts.cursor();
        let mut written = 0;
        for &(at, place) in &self.places {
            cuts.write(writer, report, written, at)?;
            written = at;
            at_place(warnings.as_deref_mut(), links, query, place, writer)?;
        }
        cuts.write(writer, report, written, report.len())
    }
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
fn write_parallel<S: QuerySections, W: ReportWrite>(
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
fn write_parallel_with<S: QuerySections, W: ReportWrite>(
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
    let mut display = DisplaySeqalign::default();
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
            // The pieces of a query whose reads write texts, and of a `-html` report, mark where
            // the serial writer writes to the diagnostics; the texts are written there
            // (`at_place`).
            let links = sections.html_links();
            let marked: Vec<bool> = chunk
                .iter()
                .map(|piece| {
                    links.is_some()
                        || warnings
                            .as_deref()
                            .is_some_and(|warnings| warnings.reads_with_texts(piece.query))
                })
                .collect();
            let indexed: Vec<(Piece, bool)> = chunk.iter().copied().zip(marked).collect();
            let buffers = map(pool, &indexed, |&(piece, marks)| {
                format_piece(sections, groups[piece.query - first].as_ref(), piece, marks)
            });
            for (piece, formatted) in chunk.iter().zip(buffers) {
                match piece.part {
                    Part::Alignments(..) | Part::Section => display.enter(&mut *writer)?,
                    Part::Head | Part::Foot => display.leave(&mut *writer)?,
                }
                if piece.part == Part::Head {
                    if let Some(warnings) = warnings.as_deref_mut() {
                        warnings.before_query(piece.query, &mut *writer)?;
                    }
                }
                formatted.write(piece.query, &mut *writer, warnings.as_deref_mut(), links)?;
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

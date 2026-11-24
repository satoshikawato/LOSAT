//! The subject passes of one blastn query batch (`SubjectPasses`): the search of one subject by stage
//! (`search_subject`, which was the closure `process_subject` of `search_query_batch`), the pass over a
//! subject set's subjects in subject order (`search_subjects`, serial, on the pool, or with
//! `LOSAT_X_BNSTREAM`), the preliminary pass (`prelim_pass`: a subject set or a database, into the HSP
//! collector or, with `-culling_limit` or `-best_hit_*`, NCBI's culling or best hit writer, with the `-ungapped` finish of each
//! subject's list) and the traceback pass (`traceback_pass`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427 (the loop of
//! `BLAST_PreliminarySearchEngine` on the subjects)
//! and c++/src/algo/blast/core/blast_traceback.c:1521-1523 (the loop of `BLAST_ComputeTraceback_MT`
//! on the subjects' batches)
//!
//! This module is a child of `run` for the types it keeps private (`SubjectSource`,
//! `PreparedSubjects`, `PrelimHit`, `SubjectHits`, `SubjectScratch`, `SubjectCtx`, the hit-list
//! helpers, `database`). Its methods are the statements that `search_query_batch` held between the
//! batch's set-up and its result, in their order; what they read of the batch is `SubjectPasses`,
//! `PrelimBatch` and `TraceBatch`, made once per batch (per subject: `special_subject`'s copies).

use super::batch_setup::BatchCutoffs;
use super::calculate_blastn_context_statistics;
use super::prelim_subject::{self, PrelimBatch, ScanBatch};
use super::traceback_subject::{self, TraceBatch};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::X_BNSTREAM_SUBJECTS_PER_THREAD;
use super::{
    besthit, culling, database, listed_prelim_hits, prelim_hits_by_subject,
    ungapped_reap_query_coverage, ungapped_subject_stage, update_hitlists_with_subject_hits,
    x_bnstream, BlastnArgs, BlastnHitList, BlastnSubjectRef, FastaRecord, GapAlignScratch,
    LookupTables, PrelimHit, PrelimHitCollector, PrelimHitLists, SubjectCtx, SubjectHits,
    SubjectScratch, SubjectSource,
};
use crate::algorithm::blastn::coordination::subject_trace_id;
use anyhow::Result;

// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1298
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// The scratch of one worker: the gapped aligner's and the subject search's.
type SubjectState = (GapAlignScratch, SubjectScratch);

/// The statistics of one subject (`search_query_batch`'s `subject_statistics`): each context's
/// cutoffs and effective search space for the subject's length with the given option spaces.
pub(super) type SubjectStatisticsFn<'a> =
    dyn Fn(usize, &[i64]) -> Result<(BatchCutoffs, Vec<i64>)> + Sync + 'a;

/// What a database's subjects read with the statistics of each subject (`DbStatistics::per_subject`):
/// their option spaces and the function of a subject's statistics. Built once per batch by
/// `search_query_batch`. The preliminary stage's own spaces and cutoffs where they are not the
/// setup's (a statistics' total of 0 or less) are the batch's (`PrelimBatch`), not these.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443 and
/// c++/src/algo/blast/core/blast_traceback.c:1598-1628 (cited at `stats::setup::subject_eff_lengths`):
/// once per subject, the preliminary stage reads the subject's own cutoffs and spaces, the
/// traceback the subject's hit-saving cutoffs and the setup's spaces.
pub(super) struct SpecialStatistics<'a> {
    /// The preliminary stage's option spaces.
    pub(super) option_spaces: &'a [i64],
    /// The traceback's own option spaces (`-searchsp`, which the memento gives it back).
    pub(super) traceback_option_spaces: &'a [i64],
    pub(super) subject_statistics: &'a SubjectStatisticsFn<'a>,
    /// The first overflow of a subject's search space, which the batch returns after its passes.
    pub(super) subject_space_overflow: &'a crate::stats::setup::SubjectSpaceOverflow,
}

/// What the finish of a database sequence's HSP list in an `-ungapped` search reads besides the
/// batch (`SubjectPasses::x_ungapped_database_stage`), made once per preliminary pass by
/// `prelim_into`: the database, the fetch's exclusion of the sequences that the negative seqid list
/// keeps without a line, and the first error of the pass's fetches.
struct UngappedStage<'a> {
    database: &'a database::BlastnDatabase,
    excluded: Option<crate::algorithm::common::database::FetchExcluded<'a>>,
    /// The first error of reading a database sequence again (`BlastnDatabase::with_blastna`) or of
    /// its fetch (`database::fetch_refusal`): the one of the smallest ordinal stays.
    error: std::sync::Mutex<Option<(u32, anyhow::Error)>>,
}

/// What the passes over a batch's subjects read of the batch, built once by `search_query_batch`
/// (shared references and `Copy` values, so the struct is `Sync`: the passes run on the pool).
///
/// NCBI's preliminary stage keeps, for each query, the preliminary HSPs of at most
/// `prelim_hitlist_size` subjects (a split batch merges those of its query chunks), and
/// its traceback reads only those.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1675-1681,1782-1786
/// A database has its own drivers of the two passes (`database::search_blocks`,
/// `database::traceback_listed`), with a scratch per worker made as the subject set's
/// parallel pass makes them.
pub(super) struct SubjectPasses<'a> {
    pub(super) args: &'a BlastnArgs,
    /// The subject set or the database.
    pub(super) subjects: SubjectSource<'a>,
    /// The values the preliminary and the traceback stage read (`search_subject` dispatches): the
    /// queries, their contexts, effective search spaces and cutoffs, and the LOSAT_TIMING counters
    /// among them.
    pub(super) prelim_batch: &'a PrelimBatch<'a>,
    pub(super) trace_batch: &'a TraceBatch<'a>,
    pub(super) any_valid_context: bool,
    /// The statistics of each subject of a database (`None` on the default path, whose subjects
    /// read the batches' cutoffs and spaces); `special_subject` searches with copies of the
    /// stage's batch that hold the subject's own.
    pub(super) special: Option<&'a SpecialStatistics<'a>>,
    /// The sizes of the preliminary and of the final hit lists.
    pub(super) prelim_hitlist_size: usize,
    pub(super) hitlist_size: usize,
    /// The algorithm of `-db_hard_mask`, whose masked letters the traceback and the ungapped
    /// re-evaluation of a database sequence read as 'B' (`database::mask_hard_letters`; decision H2).
    pub(super) db_hard_mask: Option<i32>,
    /// The batch's first fetch with `check_oid_exclusion` that read a header with a UTF8String
    /// where a VisibleString is expected (`QueryBatch::fetched_utf8`; the batch's one note, which
    /// `special_subject`'s copies share).
    pub(super) fetched_utf8: &'a crate::seqdb::Utf8Reads,
}

impl SubjectPasses<'_> {
    /// The records of a subject set (none for a database). Their count is the number of subjects that
    /// `search_subjects` searches and the number of lists into which the subject set's traceback pass
    /// groups the kept HSPs (`prelim_hits_by_subject`); their titles name the subjects in LOSAT's debug
    /// and trace lines (a database's sequences are "unknown" there).
    ///
    /// No NCBI counterpart: an accessor of LOSAT's subject set; NCBI reads each subject through a
    /// `BlastSeqSrc` (the loop of `BLAST_PreliminarySearchEngine` cited at the head of this module).
    fn subject_records(&self) -> &[FastaRecord] {
        self.subjects
            .prepared()
            .map_or(&[], |prepared| prepared.records)
    }

    /// One subject for one stage: the preliminary search (`prelim_source` is `None`: the subject's
    /// preliminary HSPs, `scan` the pass's values and tables) or the traceback of the HSPs the
    /// preliminary stage kept (`Some`: the subject's final HSPs). `state` is the worker's scratch;
    /// the cutoffs and spaces are the stage's batch's (`scan.batch`, `trace_batch`). Nothing is
    /// searched without a valid context.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1469-1470 (the subject's
    /// preliminary search) and c++/src/algo/blast/core/blast_traceback.c:1521,1632-1633 (its traceback)
    fn search_subject<const UNGAPPED: bool>(
        &self,
        s_idx: usize,
        subject: BlastnSubjectRef<'_>,
        state: &mut SubjectState,
        scan: Option<&ScanBatch<'_>>,
        prelim_source: Option<&Vec<PrelimHit>>,
    ) -> (Option<SubjectHits>, Vec<PrelimHit>) {
        let Self {
            any_valid_context,
            prelim_batch,
            trace_batch,
            ..
        } = *self;
        let (blastn_debug, blastn_trace_enabled) =
            (prelim_batch.blastn_debug, prelim_batch.blastn_trace_enabled);
        let (gap_scratch, subject_scratch) = state;
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
        // With no valid context, NCBI does not search.
        if !any_valid_context {
            return (None, Vec::new());
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
        let s_len_full = subject.length;
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
        // Made only for the debug and trace lines, which alone read it.
        let s_id_trace = if blastn_debug || blastn_trace_enabled {
            self.subject_records().get(s_idx).map_or_else(
                || "unknown".to_string(),
                |record| subject_trace_id(record.title()),
            )
        } else {
            String::new()
        };
        let s_id = s_id_trace.as_str();

        // What both stages read of this call (`SubjectCtx`; the batch's values are the stage's
        // batch's: the batch's own, or `special_subject`'s copy).
        let subject_ctx = SubjectCtx {
            s_idx,
            s_id,
            s_len_full,
        };
        if let Some(kept) = prelim_source {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1521,1632-1633
            // The traceback of the HSPs that the preliminary stage kept
            // (`traceback_subject::traceback_subject`).
            let mut subject_hits: Option<SubjectHits> = None;
            traceback_subject::traceback_subject::<UNGAPPED>(
                trace_batch,
                &subject_ctx,
                subject,
                kept,
                gap_scratch,
                subject_scratch,
                &mut subject_hits,
            );
            (subject_hits, Vec::new())
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1469-1470
            // The preliminary search of the subject (`prelim_subject::prelim_subject`).
            // The preliminary stage gives the subject's HSPs to the collector (`search_subjects`).
            (
                None,
                prelim_subject::prelim_subject::<UNGAPPED>(
                    scan.expect("the preliminary stage's calls pass the pass's values"),
                    &subject_ctx,
                    subject,
                    gap_scratch,
                    subject_scratch,
                ),
            )
        }
    }

    /// `search_subject` with the statistics of each subject (`special`): the subject's own cutoffs
    /// and spaces of the preliminary stage's option spaces, or, in the traceback, the subject's own
    /// cutoffs of the traceback's option spaces with the setup's spaces (the traceback batch's) for its
    /// e-values. The stage's batch is copied with them (`PrelimBatch` and the pass's `ScanBatch`, or
    /// `TraceBatch`), so the default path's batches keep the batch's values. The preliminary stage of
    /// an `-ungapped` search (`ungapped`) finishes the subject's list with the same copy
    /// (`x_ungapped_database_stage`), so the subject's own statistics are made once, here. An
    /// overflow of the subject's search space is noted (`subject_space_overflow`, LOSAT's refusal
    /// after the pass) and the subject has no HSPs. The callers call it apart
    /// (`stats::setup::apart`), so that the default path keeps its code.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443 and
    /// c++/src/algo/blast/core/blast_traceback.c:1598-1628 (cited at `stats::setup::subject_eff_lengths`)
    fn special_subject<const UNGAPPED: bool>(
        &self,
        s_idx: usize,
        subject: BlastnSubjectRef<'_>,
        state: &mut SubjectState,
        scan: Option<&ScanBatch<'_>>,
        prelim_source: Option<&Vec<PrelimHit>>,
        special: &SpecialStatistics<'_>,
        ungapped: Option<&UngappedStage<'_>>,
    ) -> (Option<SubjectHits>, Vec<PrelimHit>) {
        if !self.any_valid_context {
            return (None, Vec::new());
        }
        // An ungapped search has no traceback: the values of its preliminary stage stay
        // (blast_traceback.c:1526, `perform_traceback`).
        let traceback = prelim_source.is_some() && !UNGAPPED;
        let spaces = if traceback {
            special.traceback_option_spaces
        } else {
            special.option_spaces
        };
        let (cutoffs, subject_spaces) = match (special.subject_statistics)(subject.length, spaces) {
            Ok(statistics) => statistics,
            Err(error) => {
                if let Ok(overflow) = error.downcast::<crate::stats::setup::SearchSpaceOverflow>() {
                    special.subject_space_overflow.note(overflow);
                }
                return (None, Vec::new());
            }
        };
        if prelim_source.is_some() {
            let trace_batch = TraceBatch {
                query_eff_searchsp: if traceback {
                    self.trace_batch.query_eff_searchsp
                } else {
                    &subject_spaces
                },
                cutoff_scores: &cutoffs.0,
                hit_saving_cutoff_scores: &cutoffs.1,
                ..*self.trace_batch
            };
            SubjectPasses {
                trace_batch: &trace_batch,
                ..*self
            }
            .search_subject::<UNGAPPED>(s_idx, subject, state, scan, prelim_source)
        } else {
            let prelim_batch = PrelimBatch {
                query_eff_searchsp: &subject_spaces,
                cutoff_scores: &cutoffs.0,
                hit_saving_cutoff_scores: &cutoffs.1,
                x_dropoff_scores: &cutoffs.2,
                reduced_cutoff_scores: &cutoffs.3,
                ..*self.prelim_batch
            };
            let scan = scan.map(|scan| ScanBatch {
                batch: &prelim_batch,
                cutoff_scores: &cutoffs.0,
                x_dropoff_scores: &cutoffs.2,
                reduced_cutoff_scores: &cutoffs.3,
                ..*scan
            });
            let passes = SubjectPasses {
                prelim_batch: &prelim_batch,
                ..*self
            };
            let mut found =
                passes.search_subject::<UNGAPPED>(s_idx, subject, state, scan.as_ref(), None);
            if let Some(stage) = ungapped {
                passes.x_ungapped_database_stage(stage, s_idx as u32, &mut found.1);
            }
            found
        }
    }

    /// One pass over the subjects in subject order: the preliminary stage gives each subject's
    /// preliminary HSPs (`kept` is `None`), the traceback saves the final HSPs of the kept
    /// preliminary HSPs in the final hit lists of `hitlist_size` subjects
    /// (`Blast_HSPResultsInsertHSPList`).
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1700-1707
    /// Each subject's results go to `sink` in OID order, which fills the hit lists as NCBI's
    /// `Blast_HitListUpdate` does, one HSP list at a time. One thread hands each subject over as
    /// soon as it is searched, so only the HSPs the hit lists keep stay in memory. Two or more
    /// threads search all subjects of the batch and then hand them over, or, with
    /// `LOSAT_X_BNSTREAM`, hand each subject over in OID order while they search the next ones.
    fn search_subjects<const UNGAPPED: bool>(
        &self,
        kept: Option<&[Vec<PrelimHit>]>,
        scan: Option<&ScanBatch<'_>>,
        sink: &mut (dyn FnMut(usize, Option<SubjectHits>, Vec<PrelimHit>) + Send),
    ) {
        let use_parallel = self.prelim_batch.use_parallel;
        let new_state = || self.new_scratch();
        let search_one = |state: &mut SubjectState, s_idx: usize| {
            self.search_subject::<UNGAPPED>(
                s_idx,
                self.subjects
                    .prepared()
                    .expect("the subjects of a subject set")
                    .subject(s_idx),
                state,
                scan,
                kept.map(|kept| &kept[s_idx]),
            )
        };
        let num_subjects = self.subject_records().len();
        if use_parallel {
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            {
                if x_bnstream() {
                    crate::utils::xreach::x_reach!("LOSAT_X_BNSTREAM", true);
                    let threads = self.prelim_batch.parallel_pool.threads().max(1);
                    let window = threads * X_BNSTREAM_SUBJECTS_PER_THREAD;
                    self.prelim_batch.parallel_pool.install(|| {
                        x_bnstream::search_in_order(
                            num_subjects,
                            window,
                            false,
                            &new_state,
                            &search_one,
                            &mut |s_idx, (hits, prelim)| sink(s_idx, hits, prelim),
                        )
                    });
                } else {
                    let results: Vec<(Option<SubjectHits>, Vec<PrelimHit>)> =
                        self.prelim_batch.parallel_pool.install(|| {
                            (0..num_subjects)
                                .into_par_iter()
                                .map_init(new_state, search_one)
                                .collect()
                        });
                    for (s_idx, (hits, prelim)) in results.into_iter().enumerate() {
                        sink(s_idx, hits, prelim);
                    }
                }
            }
            #[cfg(any(
                not(feature = "parallel"),
                all(target_arch = "wasm32", not(feature = "wasm-threads"))
            ))]
            {
                unreachable!("use_parallel is false when parallel threads are disabled");
            }
        } else {
            if x_bnstream() {
                crate::utils::xreach::x_reach!("LOSAT_X_BNSTREAM", false);
            }
            let mut state = new_state();
            for s_idx in 0..num_subjects {
                let (hits, prelim) = search_one(&mut state, s_idx);
                sink(s_idx, hits, prelim);
            }
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:991-1041
    fn new_scratch(&self) -> SubjectState {
        (
            GapAlignScratch::new(),
            SubjectScratch::new(
                self.prelim_batch.offset_array_size,
                self.prelim_batch.window_size,
            ),
        )
    }

    /// An ungapped (`-ungapped`) search only (the other instance never calls it): the reap of a subject's
    /// HSP list by the query coverage, made as the HSPs are given to the collector (it takes nothing of
    /// the subject), the last step of the finish that NCBI's preliminary stage makes of an ungapped list.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1538
    fn x_ungapped_query_coverage(&self, hits: &mut Vec<PrelimHit>) {
        ungapped_reap_query_coverage(
            hits,
            self.args.qcov_hsp_perc.unwrap_or(0.0),
            |context_idx| self.prelim_batch.query_contexts[context_idx].len,
        );
    }

    /// The finish of a database sequence's HSP list in an `-ungapped` search (the other instance never
    /// calls it), when its HSP list is not empty only: a sequence that the fetch excludes
    /// (`stage.excluded`) loses its HSPs unread; any other is read again, unpacked (BLASTNA), and
    /// re-evaluated there with the word cutoffs and spaces it was searched with (`self.prelim_batch`:
    /// the batch's, or `special_subject`'s copy with the subject's own). A read error, or a header
    /// that the fetch cannot read (`database::fetch_refusal`), goes to `stage.error` (the one of the
    /// smallest ordinal stays) and the sequence has no HSPs.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1480,1490-1493
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2651-2665
    fn x_ungapped_database_stage(
        &self,
        stage: &UngappedStage<'_>,
        ordinal: u32,
        hits: &mut Vec<PrelimHit>,
    ) {
        // The preliminary stage's word cutoffs and spaces (`Blast_HSPListReevaluateUngapped` and the
        // e-values of an ungapped search, blast_engine.c:1492-1497).
        let &PrelimBatch {
            encoded_query_concat_blastna,
            encoded_query_concat_nomask,
            score_matrix,
            search_karlin_ref,
            round_down_evalue_score,
            prelim_reap,
            percent_identity,
            query_eff_searchsp,
            cutoff_scores,
            ..
        } = self.prelim_batch;
        if hits.is_empty() {
            return;
        }
        let excluded = match stage.excluded.as_ref() {
            Some(excluded) => excluded
                .excludes(ordinal)
                .map_err(|error| database::fetch_refusal((ordinal, error))),
            None => Ok(false),
        };
        let finished = excluded.and_then(|excluded| {
            if excluded {
                hits.clear();
                return Ok(());
            }
            stage
                .database
                .with_blastna(ordinal, self.db_hard_mask, |s_blastna| {
                    ungapped_subject_stage(
                        hits,
                        s_blastna,
                        encoded_query_concat_blastna,
                        encoded_query_concat_nomask,
                        score_matrix,
                        cutoff_scores,
                        percent_identity,
                        |prelim| {
                            let eff_searchsp = query_eff_searchsp
                                .get(prelim.context_idx as usize)
                                .copied()
                                .unwrap_or(0);
                            calculate_blastn_context_statistics(
                                prelim.prelim_score,
                                &search_karlin_ref[prelim.context_idx as usize].gapped,
                                eff_searchsp,
                                round_down_evalue_score,
                            )
                            .1
                        },
                        prelim_reap,
                    );
                    self.x_ungapped_query_coverage(hits);
                })
        });
        if let Err(error) = finished {
            hits.clear();
            let mut first = stage.error.lock().expect("the error cell is not poisoned");
            if first
                .as_ref()
                .is_none_or(|(smallest, _)| ordinal < *smallest)
            {
                *first = Some((ordinal, error));
            }
        }
    }

    /// The preliminary pass of the batch: the lists of `split_prelim_lists` (a split batch's query chunks
    /// searched before), else the HSPs of every subject, in the order of their ordinal ids, into the
    /// preliminary hit lists of `prelim_hitlist_size` subjects per query. With `-culling_limit` NCBI's
    /// culling writer takes the place of the HSP collector, with `-best_hit_*` its best hit writer. A database's read error is returned as it is
    /// (the caller reports it with `database_search_error`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1674-1679 (the preliminary search
    /// of a batch in `Blast_RunPreliminarySearchWithInterrupt`, into the batch's HSP stream)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:351-357 (the stream's writer
    /// without a filtering option: the HSP collector, `PrelimHitCollector` here; `-culling_limit` chooses
    /// the culling writer in the body below, and `-best_hit_overhang` / `-best_hit_score_edge` the best
    /// hit writer)
    pub(super) fn prelim_pass<const UNGAPPED: bool>(
        &self,
        split_prelim_lists: Option<PrelimHitLists>,
        lookup_tables: &LookupTables,
    ) -> Result<PrelimHitLists> {
        if let Some(lists) = split_prelim_lists {
            return Ok(lists);
        }
        // `-culling_limit`: NCBI's culling writer takes the place of the HSP collector
        // (`BlastHSPWriter` of `CreateHspWriter`), for a batch and for a query chunk alike
        // (a chunk's queries are its parts, so the lengths are the chunk's). Chosen once;
        // the default `else` arm below is unchanged.
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:318-342
        if self.args.culling_limit_value() > 0 {
            let mut writer = culling::PrelimCulling::new(
                self.prelim_batch.query_lengths,
                self.prelim_hitlist_size,
                self.args.culling_limit_value(),
            );
            self.prelim_into::<UNGAPPED, _>(lookup_tables, &mut |oid, prelim| {
                writer.add(oid, prelim)
            })?;
            Ok(writer.finish())
        } else if let Some(params) = self.args.best_hit() {
            // `-best_hit_overhang` / `-best_hit_score_edge`: NCBI's best hit writer, in the same
            // place, with the overhang of the preliminary stage (0).
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:305-329
            let mut writer = besthit::PrelimBestHit::new(
                self.prelim_batch.query_lengths,
                self.prelim_hitlist_size,
                params,
            );
            self.prelim_into::<UNGAPPED, _>(lookup_tables, &mut |oid, prelim| {
                writer.add(oid, prelim)
            })?;
            Ok(writer.finish())
        } else {
            let mut collector =
                PrelimHitCollector::new(self.prelim_batch.queries.len(), self.prelim_hitlist_size);
            self.prelim_into::<UNGAPPED, _>(lookup_tables, &mut |oid, prelim| {
                collector.add(oid, prelim)
            })?;
            Ok(collector.hit_lists)
        }
    }

    /// The preliminary search of every subject, each subject's HSPs (with `-ungapped`, finished) given to
    /// `sink` in the order of the subjects' ids: a subject set's by `search_subjects`, a database's by
    /// `database::search_blocks` (one instance of it for each sink, each `UNGAPPED` and each kernel: the
    /// batch's statistics, or the statistics of each subject apart, `special_subject`, chosen once per
    /// call). The subject set's and the database's drivers are the same for the HSP collector and the
    /// culling writer.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413 (the loop on the subjects,
    /// which `search_subjects` and `database::search_blocks` replace)
    fn prelim_into<const UNGAPPED: bool, M: FnMut(usize, Vec<PrelimHit>) + Send>(
        &self,
        lookup_tables: &LookupTables,
        sink: &mut M,
    ) -> Result<()> {
        // What the preliminary search of every subject and chunk of the batch reads besides
        // `PrelimBatch`: its lookup tables (`ScanBatch`), built once for the pass.
        let scan = ScanBatch::new::<UNGAPPED>(self.prelim_batch, lookup_tables);
        match self.subjects {
            SubjectSource::Prepared(_) => {
                self.search_subjects::<UNGAPPED>(None, Some(&scan), &mut |oid, _, mut prelim| {
                    if UNGAPPED {
                        self.x_ungapped_query_coverage(&mut prelim);
                    }
                    sink(oid, prelim)
                });
                Ok(())
            }
            SubjectSource::Database(database) => {
                // The subjects come one at a time in the order of their ordinal ids
                // (`database::search_blocks`), as the subject set's come in OID order.
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
                // With no valid context NCBI does not search, and reads no sequence.
                if !self.any_valid_context {
                    return Ok(());
                }
                let mut add =
                    |ordinal: u32, prelim_hits: Vec<PrelimHit>| sink(ordinal as usize, prelim_hits);
                let (parallel_pool, use_parallel) = (
                    self.prelim_batch.parallel_pool,
                    self.prelim_batch.use_parallel,
                );
                let new_scratch = || self.new_scratch();
                // An ungapped search finishes each sequence's HSP list after its search, on the
                // passes it was searched with (`x_ungapped_database_stage`, the default kernel's
                // or `special_subject`'s).
                let ungapped = UNGAPPED.then(|| UngappedStage {
                    database,
                    excluded: database::fetch_excluded(database),
                    error: std::sync::Mutex::new(None),
                });
                // The database's blocks with the kernel `$search` (a subject's search, given
                // its ordinal, its sequence and the worker's scratch), written into the blocks'
                // callback, so that the callback captures what the kernel reads and no closure
                // of it (one load less per subject).
                macro_rules! database_prelim_blocks {
                    (|$ordinal:ident, $subject:ident, $state:ident| $search:expr) => {
                        database::search_blocks(
                            database,
                            parallel_pool,
                            use_parallel,
                            &new_scratch,
                            &|$state: &mut SubjectState,
                              $ordinal: u32,
                              $subject: BlastnSubjectRef<'_>| {
                                let (_, prelim_hits) = $search;
                                (!prelim_hits.is_empty()).then_some(prelim_hits)
                            },
                            &mut add,
                        )
                    };
                }
                // The kernel, chosen once per call: the batch's statistics, the same with the
                // finish of each list of an ungapped search, or the statistics of each subject
                // apart (`special_subject`, which finishes on its copy of the batch). With a
                // MegaBLAST database index, one kernel that searches only the sequences with
                // seeds and chooses among those for each.
                let searched = match (self.prelim_batch.mb_index, self.special, ungapped.as_ref()) {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1420-1423
                    // NCBI neither reads nor searches a sequence that has no seeds in the index.
                    (Some(index), special, stage) => {
                        database_prelim_blocks!(|ordinal, subject, state| {
                            if !index.has_seeds(ordinal) {
                                return None;
                            }
                            match (special, stage) {
                                (None, None) => self.search_subject::<UNGAPPED>(
                                    ordinal as usize,
                                    subject,
                                    state,
                                    Some(&scan),
                                    None,
                                ),
                                (None, Some(stage)) => {
                                    let mut found = self.search_subject::<UNGAPPED>(
                                        ordinal as usize,
                                        subject,
                                        state,
                                        Some(&scan),
                                        None,
                                    );
                                    self.x_ungapped_database_stage(stage, ordinal, &mut found.1);
                                    found
                                }
                                (Some(special), stage) => crate::stats::setup::apart(|| {
                                    self.special_subject::<UNGAPPED>(
                                        ordinal as usize,
                                        subject,
                                        state,
                                        Some(&scan),
                                        None,
                                        special,
                                        stage,
                                    )
                                }),
                            }
                        })
                    }
                    (None, None, None) => {
                        database_prelim_blocks!(|ordinal, subject, state| self
                            .search_subject::<UNGAPPED>(
                                ordinal as usize,
                                subject,
                                state,
                                Some(&scan),
                                None
                            ))
                    }
                    (None, None, Some(stage)) => {
                        database_prelim_blocks!(|ordinal, subject, state| {
                            let mut found = self.search_subject::<UNGAPPED>(
                                ordinal as usize,
                                subject,
                                state,
                                Some(&scan),
                                None,
                            );
                            self.x_ungapped_database_stage(stage, ordinal, &mut found.1);
                            found
                        })
                    }
                    (None, Some(special), stage) => {
                        database_prelim_blocks!(|ordinal, subject, state| {
                            crate::stats::setup::apart(|| {
                                self.special_subject::<UNGAPPED>(
                                    ordinal as usize,
                                    subject,
                                    state,
                                    Some(&scan),
                                    None,
                                    special,
                                    stage,
                                )
                            })
                        })
                    }
                };
                // The re-evaluation's fetches read every sequence with HSPs, in the order of
                // the subject loop: the first that took a UTF8String for a VisibleString
                // (`FetchExcluded::excludes`) is the batch's first read of one.
                if let Some(at) = ungapped
                    .as_ref()
                    .and_then(|stage| stage.excluded.as_ref())
                    .and_then(crate::algorithm::common::database::FetchExcluded::take_utf8)
                {
                    self.fetched_utf8.note(0, at);
                }
                // An ungapped search's first error of reading a sequence again, after the
                // blocks' own.
                searched.and_then(|()| {
                    match ungapped.as_ref().and_then(|stage| {
                        stage
                            .error
                            .lock()
                            .expect("the error cell is not poisoned")
                            .take()
                    }) {
                        Some((_, error)) => Err(error),
                        None => Ok(()),
                    }
                })
            }
        }
    }

    /// The traceback pass of the batch: the final HSPs of the HSPs the preliminary pass kept
    /// (`prelim_lists`), in the final hit lists of `hitlist_size` subjects per query. A database's read
    /// error is returned as it is (the caller reports it with `database_search_error`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1700-1707
    pub(super) fn traceback_pass<const UNGAPPED: bool>(
        &self,
        prelim_lists: PrelimHitLists,
    ) -> Result<Vec<Option<BlastnHitList>>> {
        let (queries_ref, hitlist_size) = (self.prelim_batch.queries, self.hitlist_size);
        match self.subjects {
            SubjectSource::Prepared(_) => {
                let kept = prelim_hits_by_subject(prelim_lists, self.subject_records().len());
                let mut hit_lists: Vec<Option<BlastnHitList>> =
                    Vec::with_capacity(queries_ref.len());
                hit_lists.resize_with(queries_ref.len(), || None);
                self.search_subjects::<UNGAPPED>(Some(&kept), None, &mut |_, hits, _| {
                    if let Some(hits) = hits {
                        update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
                    }
                });
                Ok(hit_lists)
            }
            SubjectSource::Database(database) => {
                let (parallel_pool, use_parallel) = (
                    self.prelim_batch.parallel_pool,
                    self.prelim_batch.use_parallel,
                );
                let new_scratch = || self.new_scratch();
                let mut listed = listed_prelim_hits(prelim_lists);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1567-1577
                // The gapped traceback's fetch excludes a listed sequence that the negative
                // seqid list keeps without a line (`database::fetch_excluded`): the lists of
                // every query on it leave before the final hit lists (an `-ungapped` search
                // fetches nothing here; its re-evaluation excluded them, `prelim_pass`).
                if !UNGAPPED {
                    if let Some(excluded) = database::fetch_excluded(database) {
                        let fetched = excluded
                            .retain(&mut listed, |(ordinal, _)| Some(*ordinal))
                            .map_err(database::fetch_refusal)?;
                        // The traceback fetches every listed sequence (its first read of a
                        // UTF8String header: `FetchedUtf8::traceback_first`).
                        if let Some(ordinal) =
                            crate::algorithm::common::database::FetchedUtf8::traceback_first(
                                &fetched,
                            )
                        {
                            excluded.note_fetched(self.fetched_utf8, ordinal);
                        }
                    }
                }
                let found = database::traceback_listed(
                    database,
                    self.db_hard_mask,
                    parallel_pool,
                    use_parallel,
                    &listed,
                    &new_scratch,
                    &|state: &mut SubjectState,
                      ordinal: u32,
                      subject: BlastnSubjectRef<'_>,
                      kept: &Vec<PrelimHit>| {
                        // One test per listed subject (a subject with HSPs): the batch's statistics
                        // (`TraceBatch`), or the statistics of each subject apart (`special_subject`).
                        let (mut subject_hits, _) = match self.special {
                            None => self.search_subject::<UNGAPPED>(
                                ordinal as usize,
                                subject,
                                state,
                                None,
                                Some(kept),
                            ),
                            Some(special) => crate::stats::setup::apart(|| {
                                self.special_subject::<UNGAPPED>(
                                    ordinal as usize,
                                    subject,
                                    state,
                                    None,
                                    Some(kept),
                                    special,
                                    None,
                                )
                            }),
                        };
                        // Every listed sequence's lists are kept until the last one is traced
                        // (`found` below), so each keeps only its length (`SubjectHits::fit`). A
                        // subject set's lists go to the hit lists at once (`search_subjects` above),
                        // whose copies and `update_hitlists_with_subject_hits` hold their counts.
                        if let Some(hits) = subject_hits.as_mut() {
                            hits.fit();
                        }
                        subject_hits
                    },
                    &|kept: &Vec<PrelimHit>| {
                        kept.iter()
                            .map(|hit| (hit.prelim_ss, hit.prelim_se))
                            .collect()
                    },
                )?;
                // The final HSPs of the listed sequences in the order of their ordinal ids, as
                // the subject set's traceback saves those of its subjects
                // (`Blast_HSPResultsInsertHSPList`, `search_subjects`).
                let mut hit_lists: Vec<Option<BlastnHitList>> =
                    Vec::with_capacity(queries_ref.len());
                hit_lists.resize_with(queries_ref.len(), || None);
                for hits in found {
                    update_hitlists_with_subject_hits(&mut hit_lists, hits, hitlist_size);
                }
                Ok(hit_lists)
            }
        }
    }
}

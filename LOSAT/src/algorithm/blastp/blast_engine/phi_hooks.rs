//! psiblast `-phi_pattern` (item PS8a): round 1 as a PHI-BLAST search in blastp's call
//! (`run_call`), with the stages of `algorithm::phi_search` where blastp's differ. The pattern
//! search takes the place of the word lookup and the ungapped stage: the queries' contexts without
//! a word table (`phi_lookup`), each query's setup (`phi_start`: the Karlin block of NCBI's table,
//! the occurrences, the search space and the cutoff), the preliminary stage over the subjects
//! (`phi_prelim`: every occurrence against every subject hit, kept by blastp's preliminary hit
//! lists), and the traceback of the kept lists once every subject is scanned (`phi_traceback`: the
//! e-values read the subjects' count of pattern hits). The final lists are blastp's; the reports
//! split them by query occurrence (`BlastpHsp::phi_occurrence`, `occurrence_lists`). Every
//! function is out of line and cold: blastp's default path only tests `ResolvedBlastpArgs::phi`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1775-1778

use anyhow::Result;
use std::cell::RefCell;
use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::Mutex;

use super::{
    database, BlastpRuntimeLookup, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch,
    PrelimCollector, PrelimFilter,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{update_best_evalue, BlastpHitList, BlastpHspList};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::phi_search::engine::{self, PhiHsp, PhiSearch, PhiSetup, PhiSetupError};
use crate::algorithm::phi_search::pattern::PhiError;
use crate::algorithm::phi_search::stats;
use crate::algorithm::pssm_search::gapalign::PssmGapScratch;
use crate::algorithm::pssm_search::lookup::build_pssm_lookup_in;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::tblastx::lookup::{
    build_lookup_segments_for_context, prepare_matrix_lookup_query_for_word_size, QueryContext,
};
use crate::algorithm::tblastx::translation::QueryFrame;

thread_local! {
    static PHI_SCRATCH: RefCell<PssmGapScratch> = RefCell::new(PssmGapScratch::new());
}

/// The queries' contexts of a PHI-BLAST search (their validity and the letters that the search
/// reads, as blastp's query preparation makes them) and a word table without words: the pattern
/// scan finds the seeds (`phi_prelim`), and the table is what the call's other readers of the
/// scan's shape take.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:121-130
#[cold]
#[inline(never)]
pub(super) fn phi_lookup(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let (_, _, contexts) = prepare_matrix_lookup_query_for_word_size(
        query_frames,
        args.word_size,
        false,
        args.scoring.matrix,
        None,
    );
    let table = build_pssm_lookup_in(&QueryPssm { rows: Vec::new() }, 3, 0, &[]);
    (BlastpRuntimeLookup::Aa(table), contexts)
}

/// The Karlin block of NCBI's PHI-BLAST table for the search's matrix and gap costs
/// (`stats::phi_karlin`), or NCBI's message where the table has none: NCBI then writes the message
/// in the report of each query and searches nothing (`PhiCall::no_statistics`). The table's names
/// are compared with the matrix name as given (`strcmp`): `-matrix blosum62` has no entry (probe
/// `k.lpam.0`: `Matrix pam250 not allowed in PHI-BLAST`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:165-308
#[cold]
#[inline(never)]
pub(super) fn phi_karlin(args: &ResolvedBlastpArgs) -> Result<stats::PhiKarlin, String> {
    stats::phi_karlin(
        &args.matrix_name,
        args.scoring.gap_open,
        args.scoring.gap_extend,
    )
}

/// The searches of the call's queries and what the preliminary stage gathers for the traceback.
pub(super) struct PhiCall<'p> {
    /// Each query's search, by query.
    pub(super) searches: Vec<PhiSearch<'p>>,
    /// The pattern hits of the scanned subjects (NCBI's `Int8` `lookup_hits`; the traceback reads
    /// them as an `Int4`, `num_patterns_db`).
    lookup_hits: AtomicI64,
    /// LOSAT's refusal of a subject's scan (`phi_search::scan`), the first one met.
    refusal: Mutex<Option<String>>,
    /// NCBI's message where its table has no Karlin block for the matrix and gap costs
    /// (`phi_karlin`): no query is searched (`searches` is empty), and each query's report is
    /// the message, no hits and the footer of -1 values (`s_PHIScoreBlkFill` leaves the score
    /// block without Karlin blocks; NCBI's setup goes on, so the pattern's and the occurrences'
    /// errors still stop the run).
    pub(super) no_statistics: Option<String>,
}

/// How the call's PHI-BLAST setup ends.
pub(super) enum PhiStart<'p> {
    /// Every query is searched.
    Search(PhiCall<'p>),
    /// NCBI's error of a query's setup: the pattern's compilation (`query` 0, the first query that
    /// NCBI sets up) or the query's occurrences; NCBI's engine stops with it after the reports of
    /// the queries before (`BLAST engine error: Warning: <message> `).
    EngineError { query: usize, text: String },
}

/// The setup of each query of the call: NCBI's errors of the pattern (`PhiError::Ncbi`) and of
/// the query's occurrences (`PhiSetupError`), LOSAT's refusals of the query's scan, or the
/// searches. `db_length` and `db_num_seqs` are the database's (`-dbsize` replaces the length).
/// The Karlin block of NCBI's table, or its message where there is none, is the caller's
/// (`phi_karlin`); without it the queries' occurrences are found for their errors alone.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1065-1118
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:170-185
#[cold]
#[inline(never)]
pub(super) fn phi_start<'p>(
    args: &'p ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    karlin: Result<stats::PhiKarlin, String>,
    db_length: i64,
    db_num_seqs: i64,
) -> Result<PhiStart<'p>> {
    let phi = args
        .phi
        .as_deref()
        .expect("the setup of a PHI-BLAST search");
    let pattern = match phi {
        Ok(pattern) => pattern,
        Err(PhiError::Ncbi(message)) => {
            return Ok(PhiStart::EngineError {
                query: 0,
                text: format!("Warning: {message} "),
            })
        }
        // `validate_requested_blastp_support` has refused it.
        Err(PhiError::Refusal(message)) => anyhow::bail!("{message}"),
    };
    let rows = crate::utils::matrix::protein_score_rows(args.scoring.matrix);
    let score = |q: u8, s: u8| rows.score(q, s);
    let mut searches = Vec::with_capacity(query_frames.len());
    for (query, frames) in query_frames.iter().enumerate() {
        let frame = &frames[0];
        let letters = &frame.aa_seq[1..=frame.aa_len];
        let segments: Vec<(usize, usize)> =
            build_lookup_segments_for_context(0, frame.aa_len, &frame.seg_masks)
                .into_iter()
                .map(|(left, right)| (left as usize, right as usize))
                .collect();
        // Without NCBI's Karlin block the query is set up for its errors alone.
        let set_up = match &karlin {
            Ok(karlin) => engine::setup(&PhiSetup {
                pattern,
                query: letters,
                segments: &segments,
                score: &score,
                karlin: *karlin,
                db_length,
                db_num_seqs,
                evalue: args.evalue,
                gap_open: args.scoring.gap_open,
                gap_extend: args.scoring.gap_extend,
                x_drop_bits: args.x_drop_gapped_bits,
                x_drop_final_bits: args.x_drop_final_bits,
            })
            .map(Some),
            Err(_) => engine::occurrences(pattern, letters, &segments).map(|_| None),
        };
        match set_up {
            Ok(search) => searches.extend(search),
            Err(PhiSetupError::Refusal(message)) => anyhow::bail!("{message}"),
            Err(error) => {
                return Ok(PhiStart::EngineError {
                    query,
                    text: format!("Warning: {} ", String::from_utf8_lossy(&error.message())),
                })
            }
        }
    }
    Ok(PhiStart::Search(PhiCall {
        searches,
        lookup_hits: AtomicI64::new(0),
        refusal: Mutex::new(None),
        no_statistics: karlin.err(),
    }))
}

/// The preliminary HSP that the hit lists keep for a PHI HSP of query `q_idx` and subject
/// `s_idx`: its extent and score (the e-value of the preliminary stage is 0, NCBI computes none).
fn preliminary_of(hsp: &PhiHsp, q_idx: usize, s_idx: usize, query_length: usize) -> XPrelimHsp {
    XPrelimHsp {
        hsp: BlastpPreliminaryHsp {
            query_start: hsp.query_start,
            query_end: hsp.query_end,
            subject_start: hsp.subject_start,
            subject_end: hsp.subject_end,
            gapped_query_start: hsp.q_pattern_offset,
            gapped_subject_start: hsp.s_pattern_offset,
            raw_score: hsp.score,
            query_context: q_idx as i32,
            query_frame: 0,
            subject_frame: 0,
            query_length,
            q_idx: q_idx as u32,
            s_idx: s_idx as u32,
        },
        e_value: 0.0,
    }
}

/// The preliminary stage of a PHI-BLAST search (`PHIBlastWordFinder` and `PHIGetGappedScore`):
/// each subject's pattern hits against every query occurrence (`PhiSearch::prelim`), the HSPs at
/// the cutoff or above kept by the preliminary hit lists (blastp's collector, or the writer of
/// `-culling_limit` / `-best_hit_*`), and the subject's hits counted for the traceback. A
/// database is read in blocks as blastp's (`database::search_blocks`), a subject set in turn.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1016-1046
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn phi_prelim(
    phi: &PhiCall<'_>,
    contexts: &[QueryContext],
    window: i32,
    new_scratch: &(dyn Fn() -> BlastpSubjectScratch + Sync),
    database: Option<&database::BlastpDatabase>,
    subjects: &[EncodedProtein],
    pool: &crate::utils::threading::SearchPool<'_>,
    prelim_filter: Option<&mut PrelimFilter>,
    prelim_hitlist_size: usize,
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> Result<()> {
    // Without NCBI's Karlin block nothing is searched (`PhiCall::no_statistics`).
    if phi.searches.is_empty() {
        return Ok(());
    }
    let process_subject = |_: &mut BlastpSubjectScratch,
                           s_idx: usize,
                           subject: BlastpSubjectRef<'_>,
                           subject_result: &mut BlastpSubjectResult| {
        let letters = &subject.aa_seq[1..];
        let mut counted = false;
        for (q_idx, search) in phi.searches.iter().enumerate() {
            if !contexts[q_idx].is_valid {
                continue;
            }
            let found = PHI_SCRATCH
                .with(|scratch| search.prelim(letters, subject.aa_len, &mut scratch.borrow_mut()));
            let (hsps, hits) = match found {
                Ok(found) => found,
                Err(PhiError::Refusal(message)) | Err(PhiError::Ncbi(message)) => {
                    let mut refusal = phi.refusal.lock().unwrap_or_else(|e| e.into_inner());
                    refusal.get_or_insert(message);
                    return;
                }
            };
            // The hits are the subject's, whatever the query: counted once per subject.
            if !counted {
                phi.lookup_hits.fetch_add(hits as i64, Ordering::Relaxed);
                counted = true;
            }
            if hsps.is_empty() {
                continue;
            }
            let query_length = contexts[q_idx].aa_len;
            let list = BlastpHspList {
                oid: s_idx as u32,
                query_index: q_idx as u32,
                hsps: hsps
                    .iter()
                    .map(|hsp| preliminary_of(hsp, q_idx, s_idx, query_length))
                    .collect(),
                best_evalue: 0.0,
            };
            subject_result.x_compact_hsp_lists.push((q_idx, list));
        }
    };
    let searched = match (database, prelim_filter) {
        (Some(database), None) => database::search_blocks(
            database,
            pool,
            window,
            &new_scratch,
            &process_subject,
            &mut PrelimCollector {
                x_compact_hit_lists,
                prelim_hitlist_size,
            },
        ),
        (Some(database), Some(filter)) => database::search_blocks(
            database,
            pool,
            window,
            &new_scratch,
            &process_subject,
            filter,
        ),
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1498
        (None, mut prelim_filter) => {
            let mut scratch = new_scratch();
            let mut subject_result = BlastpSubjectResult::new(0);
            for (s_idx, subject) in subjects.iter().enumerate() {
                process_subject(
                    &mut scratch,
                    s_idx,
                    BlastpSubjectRef::of(subject),
                    &mut subject_result,
                );
                if let Some(filter) = prelim_filter.as_deref_mut() {
                    filter.add(&mut subject_result);
                    continue;
                }
                for (q_idx, list) in subject_result.x_compact_hsp_lists.drain(..) {
                    x_compact_hit_lists[q_idx]
                        .get_or_insert_with(|| BlastpHitList::new(prelim_hitlist_size))
                        .update(list);
                }
            }
            Ok(())
        }
    };
    searched?;
    if let Some(message) = phi.refusal.lock().unwrap_or_else(|e| e.into_inner()).take() {
        anyhow::bail!("{message}");
    }
    Ok(())
}

/// The traceback of a PHI-BLAST search (`s_PHITracebackFromHSPList` of each kept list): the
/// subject's preliminary HSPs found again (`PhiSearch::prelim` gives the same list), each traced
/// with the final X-drop, its e-value from the count of pattern hits of every scanned subject
/// (`num_patterns_db`, NCBI's `(Int4)` of its `Int8` sum), reaped by `-evalue`; the final HSPs
/// with their counts (`pssm_search::traceback0::final_hsp`) and their occurrence, into the final
/// hit lists as blastp's traceback inserts them (`traceback::hit_lists_of_traced`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1636-1641
#[cold]
#[inline(never)]
pub(super) fn phi_traceback(
    phi: &PhiCall<'_>,
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    subjects: XRedoSubjects<'_>,
    compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    sorts_results: bool,
) -> Result<Vec<Option<BlastpHitList>>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1775-1778
    let num_patterns_db = phi.lookup_hits.load(Ordering::Relaxed) as i32;
    let mut traced: Vec<(usize, Option<BlastpHspList>)> = Vec::new();
    let mut scratch = PssmGapScratch::new();
    for (q_idx, hit_list) in compact_hit_lists.iter_mut().enumerate() {
        let Some(hit_list) = hit_list.take() else {
            continue;
        };
        let search = &phi.searches[q_idx];
        let ctx = &contexts[q_idx];
        let query_nomask = super::query_nomask_sequence(ctx);
        let count = hit_list.hsplist_count;
        for list in hit_list.hsplist_array.into_iter().take(count) {
            let read = subjects.read(list.oid)?;
            let residues = read.residues();
            let mut subject = Vec::with_capacity(residues.len() + 1);
            subject.extend_from_slice(residues);
            subject.push(0);
            let (hsps, _) = search
                .prelim(&subject, residues.len(), &mut scratch)
                .map_err(|error| match error {
                    PhiError::Refusal(message) | PhiError::Ncbi(message) => {
                        anyhow::anyhow!("{message}")
                    }
                })?;
            let kept = search
                .traceback(
                    &subject,
                    residues.len(),
                    hsps,
                    num_patterns_db,
                    args.evalue,
                    &mut scratch,
                )
                .map_err(|error| match error {
                    PhiError::Refusal(message) | PhiError::Ncbi(message) => {
                        anyhow::anyhow!("{message}")
                    }
                })?;
            if kept.is_empty() {
                traced.push((q_idx, None));
                continue;
            }
            let hsps = kept
                .into_iter()
                .map(|hsp| {
                    let mut preliminary =
                        preliminary_of(&hsp, q_idx, list.oid as usize, ctx.aa_len).hsp;
                    preliminary.query_end = hsp.query_end;
                    preliminary.subject_end = hsp.subject_end;
                    let mut final_hsp = crate::algorithm::pssm_search::traceback0::final_hsp(
                        &preliminary,
                        hsp.edit_script,
                        query_nomask,
                        residues,
                        args.scoring.matrix,
                    );
                    final_hsp.e_value = hsp.evalue;
                    final_hsp.bit_score = hsp.bit_score;
                    final_hsp.phi_occurrence = hsp.occurrence as u16;
                    final_hsp
                })
                .collect();
            let mut traced_list = BlastpHspList {
                oid: list.oid,
                query_index: list.query_index,
                hsps,
                best_evalue: 0.0,
            };
            update_best_evalue(&mut traced_list);
            traced.push((q_idx, Some(traced_list)));
        }
    }
    Ok(super::traceback::hit_lists_of_traced(
        traced,
        compact_hit_lists.len(),
        args.traceback_hitlist_size(),
        sorts_results,
    ))
}

/// One occurrence's part of a query's final hit list (`PHIBlast_HSPResultsSplit`): the subjects
/// with an HSP of the occurrence, each with those HSPs in their order, sorted by
/// `s_EvalueCompareHSPLists` (best e-value, then top score, then the higher subject first; the
/// lists' best e-values are all 0).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3572-3627
#[cold]
#[inline(never)]
pub(super) fn occurrence_list(
    hit_list: Option<&BlastpHitList>,
    occurrence: usize,
) -> Option<BlastpHitList> {
    let hit_list = hit_list?;
    let lists: Vec<BlastpHspList> = hit_list
        .hsplist_array
        .iter()
        .take(hit_list.hsplist_count)
        .filter_map(|list| {
            let hsps: Vec<_> = list
                .hsps
                .iter()
                .filter(|hsp| usize::from(hsp.phi_occurrence) == occurrence)
                .cloned()
                .collect();
            // NCBI's lists of the split keep the best e-value of a new list, 0
            // (`Blast_HSPListNew`), which the sort below reads: they sort by their first HSP's
            // score, then by subject.
            (!hsps.is_empty()).then(|| BlastpHspList {
                oid: list.oid,
                query_index: list.query_index,
                hsps,
                best_evalue: 0.0,
            })
        })
        .collect();
    if lists.is_empty() {
        return None;
    }
    let mut part = BlastpHitList::new(lists.len());
    part.hsplist_count = lists.len();
    part.hsplist_array = lists;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3616-3619
    part.sort_by_evalue();
    Some(part)
}

/// The footers of PHI-BLAST's outfmt 0 report: each query's ungapped and gapped Karlin blocks are
/// the pattern's (`s_PHIScoreBlkFill`: Lambda, K and H of the table, in both), its search space
/// the pattern's (`stats::search_space`), and the Gumbel columns, where the report prints them
/// (no `-searchsp`), zeros (PHI-BLAST's Gumbel block is not filled). Without the table's Karlin
/// block (`PhiCall::no_statistics`) each footer is that of a query not searched: -1 in both
/// blocks, no Gumbel columns, a search space of 0, also with `-searchsp`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:504-506
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
pub(super) fn phi_footers(
    phi: &PhiCall<'_>,
    queries: &mut [crate::report::BlastpPairwiseQuery],
    report: &mut crate::report::BlastpPairwiseReport,
) {
    if phi.no_statistics.is_some() {
        for query in queries.iter_mut() {
            query.batch_skipped = true;
        }
        return;
    }
    let karlin =
        |search: &PhiSearch<'_>, base: crate::stats::KarlinParams| crate::stats::KarlinParams {
            lambda: search.karlin.lambda,
            k: search.karlin.k,
            h: search.karlin.h,
            ..base
        };
    for (query, search) in queries.iter_mut().zip(&phi.searches) {
        query.ungapped_karlin = karlin(search, query.ungapped_karlin);
        query.effective_search_space = search.search_space;
    }
    if let Some(search) = phi.searches.first() {
        report.gapped_karlin = karlin(search, report.gapped_karlin);
    }
    report.gumbel = report.gumbel.map(|gumbel| crate::stats::BlastGumbelBlk {
        a: 0.0,
        alpha: 0.0,
        sigma: 0.0,
        a_un: 0.0,
        alpha_un: 0.0,
        ..gumbel
    });
}

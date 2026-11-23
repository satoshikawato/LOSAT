//! TBLASTX search coordination and execution
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:438-1497
//!
//! This module contains the main `run()` function
//! that coordinate the TBLASTX search process:
//! - Query and subject sequence reading
//! - Lookup table construction
//! - Parallel subject scanning with two-hit extension
//! - HSP reevaluation, linking, and culling
//! - Output generation
//!
//! # Module Structure
//!
//! - `mod.rs` - Shared structures and helpers (WorkerState, atomic ops, reevaluation)
//! - `run_impl.rs` - Main run() function for backbone mode

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1427 (main search loop)
mod run_impl;

// Re-export the main run function
pub use run_impl::check_options;
pub use run_impl::run;
pub use run_impl::run_local;
pub use run_impl::run_local_database;
// EXPERIMENT (PAIRPAR): the summary line `main` prints at exit.
pub use run_impl::x_pair_par_print_summary;
// The seed release's `LOSAT_X_STATS` line `main` prints at exit. No NCBI
// counterpart: it prints counters and changes no value NCBI computes.
pub use run_impl::x_seed_release_print_summary;
// What every report of a run needs (`blast_xml`).
pub(crate) use run_impl::ReportDatabase;
pub(crate) use run_impl::TblastxReportRun;
// The error of a query batch's search that comes after the reports of the batches before (the
// refusals of `extension::query_batch`).
pub(crate) use run_impl::ErrorAfterBatchReports;

// Imports used by submodules
pub(crate) use anyhow::Result;
pub(crate) use indicatif::ProgressBar;
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1297
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use rayon::prelude::*;
pub(crate) use std::collections::HashSet;
pub(crate) use std::sync::atomic::{AtomicI32, AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
pub(crate) use std::sync::mpsc::channel;
pub(crate) use std::time::Instant;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
pub(crate) use crate::api::local_blast::{FormatProbe, OutputSink, ReportOutputs};
pub(crate) use crate::common::Hit;
pub(crate) use crate::config::ScoringMatrix;
pub(crate) use crate::stats::{lookup_protein_params_ungapped, KarlinParams};
pub(crate) use crate::utils::genetic_code::GeneticCode;
pub(crate) use crate::utils::seg::SegMasker;

pub(crate) use super::args::{TblastxArgs, TblastxMatrix};
pub(crate) use super::chaining::UngappedHit;
pub(crate) use super::constants::GAP_TRIGGER_BIT_SCORE;
pub(crate) use super::diagnostics::{
    diagnostics_enabled, print_summary as print_diagnostics_summary, DiagnosticCounters,
};
pub(crate) use super::extension::{
    convert_coords, extend_hit_two_hit_fitting, extend_hit_two_hit_with, QueryBatch,
};
// EXPERIMENT (LOSAT_X_TBXUNGAPSENT): extension/x_two_hit_fast.rs.
pub(crate) use super::extension::{x_tbx_sent_contexts, x_tbx_sent_subject};
#[cfg(test)]
pub(crate) use super::lookup::build_ncbi_lookup;
pub(crate) use super::lookup::QueryContext;
pub(crate) use super::ncbi_cutoffs::{
    cutoff_score_for_update_tblastx, cutoff_score_max_for_tblastx, gap_trigger_raw_score,
    x_drop_raw_score, BLAST_GAP_DECAY_RATE,
};
pub(crate) use super::reevaluate::{
    get_num_identities_and_positives_ungapped, reevaluate_ungapped_hit_ncbi_translated_with,
};
pub(crate) use super::report;
pub(crate) use super::scoring::{Blosum62Score, MatrixScore, TbxScore};
pub(crate) use super::stage_dump;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:553-558
pub(crate) use super::sum_stats_linking::{
    apply_sum_stats_even_gap_linking_with_parallel, calculate_link_hsp_cutoffs_ncbi,
    compute_avg_query_length_ncbi, find_smallest_lambda_params, link_one_hit_with_int4_ends,
    LinkingParams,
};
pub(crate) use super::tracing::{
    trace_final_hit_if_match, trace_hsp_target, trace_match_target, trace_ungapped_hit_if_match,
};
pub(crate) use super::translation::{generate_frames, QueryFrame};
pub(crate) use crate::stats::karlin::bit_score as calc_bit_score;

// Import OffsetPair from scan submodule
pub(crate) use super::scan::OffsetPair;

// Import DiagStruct from blast_extend module (NCBI blast_extend.c equivalent)
pub(crate) use super::blast_extend::DiagStruct;

// Import InitHSP and related functions from blast_gapalign module (NCBI blast_gapalign.c equivalent)
pub(crate) use super::blast_gapalign::{
    get_ungapped_hsp_list, get_ungapped_hsp_list_int4, sort_init_hsps_by_score_ncbi,
    sort_ungapped_hits_by_score_int4, sort_ungapped_hits_by_score_ncbi, trace_init_hsp_if_match,
    ungapped_hits_is_sorted_by_score_ncbi, InitHSP,
};

// Import subject scanning functions from blast_aascan module (NCBI blast_aascan.c equivalent)
#[cfg(test)]
pub(crate) use super::blast_aascan::s_blast_aa_scan_subject;

// ---------------------------------------------------------------------------
// Shared structures
// ---------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-148
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-99
/// One HSP of the final TBLASTX result: the tabular record, and the parts of NCBI's
/// `BlastHSP` that the pairwise report adds (the Seq-align of `blast_seqalign.cpp`).
#[derive(Debug, Clone)]
pub(crate) struct TblastxHsp {
    pub hit: Hit,
    /// `BlastHSP::subject.frame`.
    pub subject_frame: i8,
    /// `BlastHSP::num` (the Seq-align score "sum_n" when it is above 1).
    pub num: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:72-100
/// What the reports show of one query apart from its HSPs.
#[derive(Debug, Clone, Default)]
pub(crate) struct TblastxQueryStats {
    /// The Karlin block of the query's first valid context; `None` when no context of
    /// the query is valid.
    pub karlin: Option<KarlinParams>,
    /// That context's effective search space (0 without a valid context).
    pub eff_searchsp: i64,
    /// That context's length in residues (0 without a valid context): the BLAST XML (outfmt 5)
    /// computes the context's length adjustment again from it at report time
    /// (`blast_xml::query_stats`), so that no search path keeps the adjustment.
    pub aa_len: usize,
    /// The effective search space of the query's context 0 (frame +1), 0 when that context is
    /// not valid: the `EffectiveSearchSpace` option of a BLAST archive (outfmt 11) takes the
    /// first context of the batch's first query, not its first valid one.
    pub context0_search_space: i64,
    /// That context's length in residues, 0 when it is not valid (`archive::RequestLength::space`).
    pub context0_aa_len: i64,
    /// The SEG intervals (residues of the frame, end exclusive) of each translated frame.
    pub seg_masks: Vec<(i8, Vec<(usize, usize)>)>,
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-195
/// The result of one query batch. The query indices of `hits` and `queries` are the
/// queries' positions in the batch. `hits` are in report order (`hit_collector.rs`: each
/// query's `-max_target_seqs` subjects in NCBI's order), so the batches' hits follow each other
/// in the run's report.
#[derive(Debug, Default)]
pub(crate) struct TblastxBatch {
    pub hits: Vec<TblastxHsp>,
    /// False when the batch has no valid context, so NCBI did not search it
    /// (`CheckInternalData`, local_blast.cpp:177-208).
    pub searched: bool,
    pub queries: Vec<TblastxQueryStats>,
    /// EXPERIMENT (LOSAT_X_TBXPIPE): the batch's `LOSAT_TIMING` block, which a search of the
    /// batch pipe keeps here for the batch loop to print when it takes the batch; empty
    /// otherwise (the block is printed at the end of the search). No NCBI counterpart.
    pub x_timing: Vec<String>,
    /// The first of the batch's re-evaluation fetches (`check_oid_exclusion`) that read a header
    /// with a UTF8String where a VisibleString is expected (`FetchExcluded::excludes`, in
    /// increasing ordinal): NCBI's warning comes there, before its Seq-aligns' reads.
    pub fetched_utf8: Option<crate::seqdb::SeqLocation>,
}

pub(crate) struct WorkerState<S> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1475
    /// The channel to the writer thread (`hit_collector.rs`), which takes each subject's final
    /// HSPs with the subject's index; `None` when the subject loop keeps the collector itself.
    pub tx: Option<std::sync::mpsc::Sender<(usize, Vec<TblastxHsp>)>>,
    /// The final HSPs of the subject being searched, sent at the end of the subject.
    pub hits: Vec<TblastxHsp>,
    pub offset_pairs: Vec<OffsetPair>,
    pub diag_array: Vec<DiagStruct>,
    pub diag_offset: i32,
    /// The scores of the search's matrix (`scoring.rs`), which the subject's closures take as an
    /// argument.
    pub sc: S,
}

// ---------------------------------------------------------------------------
// Atomic helpers
// ---------------------------------------------------------------------------

#[inline]
pub(crate) fn atomic_max_usize(dst: &AtomicUsize, val: usize) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
pub(crate) fn atomic_min_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val < cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
pub(crate) fn atomic_max_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

// ---------------------------------------------------------------------------
// HSP reevaluation
// ---------------------------------------------------------------------------

/// NCBI Blast_HSPListReevaluateUngapped equivalent
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2609-2737, c++/src/algo/blast/core/blast_engine.c:1492-1497
///
/// NCBI code flow:
/// 1. For each HSP in the list
/// 2. Get context and compute query_start (context start position)
/// 3. Call Blast_HSPReevaluateWithAmbiguitiesUngapped
/// 4. Call Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
/// 5. Delete HSPs that fail tests
/// 6. Sort by score (scores may have changed)
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2694-2706
///
/// This function performs batch reevaluation on all HSPs in the list,
/// updating scores and trimming boundaries. HSPs that fail reevaluation
/// or subsequent tests are removed from the list.
pub(crate) fn reevaluate_ungapped_hsp_list<S: TbxScore>(
    sc: S,
    ungapped_hits: Vec<UngappedHit>,
    contexts: &[QueryContext],
    s_frames_scoring: &[QueryFrame],
    s_frames_report: &[QueryFrame],
    cutoff_scores: &[i32],
    timing_enabled: bool,
    reeval_ns: &AtomicU64,
    reeval_calls: &AtomicU64,
    is_long_sequence: bool,
    stats_hsp_filtered_by_reeval: &mut usize,
) -> Vec<UngappedHit> {
    let mut kept_hits = Vec::new();

    for mut hit in ungapped_hits {
        // No HSP of length 0 or less arrives (its end, an `Int4` at or below its start, would read
        // here as a `u32`): the one-hit word finder's path drops those HSPs after its first linking
        // call, where NCBI's re-evaluation (this function) deletes them
        // (`link_one_hit_with_int4_ends`).
        debug_assert!(
            !hit.is_int4_inverted(),
            "an HSP of length 0 or less reached the re-evaluation"
        );
        let ctx = &contexts[hit.ctx_idx as usize];
        let s_frame_scoring = &s_frames_scoring[usize::from(hit.s_f_idx)];
        let s_frame_report = &s_frames_report[usize::from(hit.s_f_idx)];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:686-688
        let cutoff = cutoff_scores[hit.ctx_idx as usize];

        // Start of this context's residues in the query buffer (`query_blk->sequence` offset by `query_info->contexts[context].query_offset`).
        // query->sequence points past the leading NULLB sentinel.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:311-315, c++/src/algo/blast/core/blast_util.c:112-116.
        let query_full = &ctx.aa_seq;
        let query = &query_full[1..query_full.len() - 1];
        let subject_full = &s_frame_scoring.aa_seq;
        let subject = &subject_full[1..subject_full.len() - 1];

        // Context-relative 0-based query and subject start offsets (`hsp->query.offset`).
        let qs = hit.q_aa_start as usize;
        let ss = hit.s_aa_start as usize;
        let len_u = (hit.q_aa_end as usize).saturating_sub(hit.q_aa_start as usize);

        if len_u == 0 {
            continue;
        }

        // NCBI: Blast_HSPReevaluateWithAmbiguitiesUngapped (c++/src/algo/blast/core/blast_hits.c:2702-2705)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:675-733
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:688
        // If score < cutoff_score, the HSP is deleted (s_UpdateReevaluatedHSPUngapped returns FALSE)
        let t0 = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let (new_qs, new_ss, new_len, new_score) = if let Some(result) =
            reevaluate_ungapped_hit_ncbi_translated_with(sc, query, subject, qs, ss, len_u, cutoff)
        {
            result
        } else {
            // Deleted by NCBI reevaluation (score < cutoff_score)
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:730-732 (s_UpdateReevaluatedHSPUngapped returns FALSE)
            if is_long_sequence {
                *stats_hsp_filtered_by_reeval += 1;
            }
            if let Some(t) = t0 {
                let elapsed = t.elapsed();
                reeval_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                reeval_calls.fetch_add(1, AtomicOrdering::Relaxed);
            }
            continue;
        };
        if let Some(t) = t0 {
            let elapsed = t.elapsed();
            reeval_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
            reeval_calls.fetch_add(1, AtomicOrdering::Relaxed);
        }

        // NCBI: Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2720
        // NCBI uses sequence_nomask that is also offset past the leading NULLB.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381-1382, c++/src/algo/blast/core/blast_util.c:112-116.
        let query_nomask_full = ctx.aa_seq_nomask.as_ref().unwrap_or(&ctx.aa_seq);
        let query_nomask = &query_nomask_full[1..query_nomask_full.len() - 1];
        // NCBI uses the HSP's scoring subject for reevaluation above, then
        // recomputes identities from the reporting sequence buffers.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2717
        //
        // For local `-subject` TBLASTX, NCBI BLAST+ scores the translated
        // subject with the search-time subject buffer, while qseq/sseq and
        // identity reporting follow the requested DB genetic code. Keep those
        // buffers separate so code-4 TGA->W display parity does not change the
        // raw-score path.
        let subject_report_full = &s_frame_report.aa_seq;
        let subject_report = &subject_report_full[1..subject_report_full.len() - 1];
        let num_ident = if let Some((num_ident, _align_length)) =
            get_num_identities_and_positives_ungapped(
                query_nomask,
                subject_report,
                new_qs,
                new_ss,
                new_len,
            ) {
            num_ident
        } else {
            continue;
        };

        // NCBI: Blast_HSPTest (c++/src/algo/blast/core/blast_hits.c:2719)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1444-1449 (defaults are zero via calloc in BlastHitSavingOptionsNew)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-999 (Blast_HSPTest uses percent_identity/min_hit_length)
        let percent_identity = 0.0;
        let min_hit_length = 0;
        if crate::core::blast_hits::hsp_test(
            num_ident as i64,
            new_len as i64,
            percent_identity,
            min_hit_length,
        ) {
            continue;
        }

        // Update hit with reevaluated coordinates and score (0-based offsets)
        let new_qe = new_qs + new_len;
        let new_se = new_ss + new_len;
        hit.q_aa_start = u32::try_from(new_qs)
            .expect("q_aa_start: a query offset fits u32 (NCBI's BlastSeg offset is Int4)");
        hit.q_aa_end = u32::try_from(new_qe)
            .expect("q_aa_end: a query offset fits u32 (NCBI's BlastSeg end is Int4)");
        hit.s_aa_start = u32::try_from(new_ss)
            .expect("s_aa_start: a subject offset fits u32 (NCBI's BlastSeg offset is Int4)");
        hit.s_aa_end = u32::try_from(new_se)
            .expect("s_aa_end: a subject offset fits u32 (NCBI's BlastSeg end is Int4)");
        hit.raw_score = new_score;
        hit.num_ident = u32::try_from(num_ident)
            .expect("num_ident: an identity count fits u32 (NCBI's num_ident is Int4)");

        trace_ungapped_hit_if_match("after_reevaluate", &hit);
        kept_hits.push(hit);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2733-2734
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353
    if !ungapped_hits_is_sorted_by_score_ncbi(&kept_hits) {
        sort_ungapped_hits_by_score_ncbi(&mut kept_hits);
    }

    kept_hits
}

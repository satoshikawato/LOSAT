//! The kernel of a BLASTP subject's preliminary search (`BlastpPrelimKernel::search_subject`, which
//! `run_call`'s `process_subject` and `process_own_subject` call once per subject): the scan of the
//! subject for the words of the call's lookup table with the two-hit test (the one-hit test with
//! `-window_size 0`) and the ungapped extension of each hit that passes it, with BLOSUM62 or the
//! search's matrix, on the standard table, the compressed one of blastp-fast or the ranges of
//! `-lcase_masking`; the diagonal offset moved past the subject; then its gapped half
//! (`gapped_subject`): the subject's ungapped hits sorted by score and, with `-chaining` and
//! BLOSUM62, chained; the queries whose gapped extensions start restricted; the gapped extension of
//! each hit outside the interval tree of the HSPs kept so far, with the exact retry of a query whose
//! restricted score is inconclusive; the purge of the HSPs with common endpoints; and the subject's
//! HSP lists per query with their E-values, reaped by the preliminary E-value. The two methods are
//! the kernel's statements at their place, in their order.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200,214-215,223-224
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:440,713
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-526
//!
//! This module is a child of `blast_engine` for the private items of that file it reads (the
//! subject's scratch, reference and result, the preliminary values, the initial hits, the lookup
//! table, the clocks and the helpers of the scan and of the gapped stage) and its sibling modules
//! `x_hits_flat`, `x_ih_pool` and `x_soft_subject`.

use std::sync::Arc;

use super::{
    blastp_adjust_hsp_offsets_and_get_query_data, blastp_advance_preliminary_hit_index,
    blastp_context_idx_for_query_offset_with_bounds, blastp_for_each_offset_pair,
    blastp_for_each_offset_pair_of, blastp_get_ungapped_hsp_context_with_bounds,
    blastp_restricted_align_enabled, blastp_save_init_hsp, blastp_timing_start,
    blastp_trace_init_hsp_overlaps_target, build_restricted_align_array_into,
    chain_blastp_init_hsps, extend_preliminary_blastp_hit, ncbi_qsort_init_hsps_by_score,
    ncbi_qsort_preliminary_hits, preliminary_hits_is_sorted_by_score_ncbi,
    preliminary_score_compare, preliminary_tree_hsp, prepare_preliminary_hits_for_kappa_in_place,
    rebuild_preliminary_interval_tree_from_list, remove_preliminary_hits_for_query,
    should_retry_exact_preliminary_gapped_alignment, split_preliminary_hits_by_query_into,
    x_bp_tree_max_cells, x_build_compact_preliminary_hsp_list, x_hits_flat, x_ih_pool, x_outlined,
    x_soft_subject, BlastpContextLookupBounds, BlastpOffsetPairScanMode, BlastpPrelimValues,
    BlastpRuntimeLookup, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, BlastpTiming,
    InitHSP, BLASTP_RESTRICTED_MULT,
};
use crate::algorithm::blastn::interval_tree::TreeHsp;
use crate::algorithm::blastp;
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::extension::{
    extend_one_hit_blosum62, extend_one_hit_for, extend_one_hit_of, extend_two_hit_blosum62,
    extend_two_hit_table,
};
use crate::algorithm::blastp::gapalign::{BlastpGappedAlignmentMode, BlastpPreliminaryHsp};
use crate::algorithm::blastp::hsp::{reap_hsplist_by_evalue, update_best_evalue};
use crate::algorithm::blastp::kappa::{blastp_trace_hsp_target, blastp_trace_matches_pair};
use crate::algorithm::blastp::x_seed_fused::XContextTable;
use crate::algorithm::blastp::x_ungapped_sent::XSentQuery;
use crate::algorithm::common::diag_entry::{DiagRing, DiagTable};
use crate::algorithm::common::query_batch_tables::QueryBatchTables;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{X_BLASTP, X_BLASTP_MATRIX};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::config::ScoringMatrix;
use crate::core::blast_hits::reap_by_query_coverage;
use crate::stats::{BlastGumbelBlk, KarlinParams};

/// What the kernel reads of `run_call` (`search_subject` and its gapped half `gapped_subject`): each
/// field is `run_call`'s value of the same name, which the kernel closure captured before (references
/// of the captures' types, `Copy` values by value), so that the moved statements read them as they
/// did in the closure. `run_call` makes it once, before the kernel's callers `process_subject` and
/// `process_own_subject`.
pub(super) struct BlastpPrelimKernel<'a> {
    // The call's options and its queries.
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) query_ids: &'a Vec<Arc<[u8]>>,
    pub(super) contexts: &'a Vec<QueryContext>,
    pub(super) context_lookup_bounds: BlastpContextLookupBounds,
    pub(super) concatenated_query: &'a Vec<u8>,
    // The scan's lookup table, the context of each offset of the concatenated query, the ungapped
    // X-drop of each context, the window and the word size, the query block of LOSAT_X_BPUNGAPSENT's
    // sentinel walks, the diagonals' mask, the span of the queries (the interval tree's bounds) and
    // the scan ranges of `-lcase_masking`'s subjects. `x_ungap_sent` is made from
    // `x_dropoff_per_context` (`call_setup::scan_shape`): the SAFETY of the sentinel walks reads the
    // one through the other.
    pub(super) lookup: &'a BlastpRuntimeLookup,
    pub(super) x_context_table: &'a XContextTable,
    pub(super) x_dropoff_per_context: &'a Vec<i32>,
    pub(super) window: i32,
    pub(super) word_size: i32,
    pub(super) x_ungap_sent: &'a Option<XSentQuery>,
    pub(super) diag_mask: i32,
    pub(super) total_query_span: i32,
    pub(super) x_soft_subjects: Option<&'a [Vec<(i32, i32)>]>,
    // The gapped extension's X-drop and the statistics of the preliminary E-values.
    pub(super) x_drop_gapped: i32,
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) gapped_gumbel: &'a Option<BlastGumbelBlk>,
    // The trace's first query (LOSAT_DEBUG_BLASTP_TRACE_HSP) and the clocks (LOSAT_TIMING).
    pub(super) x_trace_query_base: u32,
    pub(super) timing_enabled: bool,
    pub(super) timing: &'a Option<Arc<BlastpTiming>>,
    /// Set where the interval tree insertion of an HSP of a whole subject never ends (its HSPs lie
    /// inside the tree's bounds and are in order, so never): the call then stops with an error, as
    /// BLASTX's search does (`blastx::search::subject_hits_with`).
    pub(super) tree_without_end: &'a std::sync::atomic::AtomicBool,
}

/// One subject with letters as `search_subject` hands it to its gapped half (`gapped_subject`): the
/// subject and its preliminary statistics, with the values that `search_subject` made of them,
/// the subject's letters after the sentinel byte, its length (NCBI's `subject->length`) and the
/// word cutoffs (one producer each). A hard chunk shorter than the word (`one_word_gapped`) has the
/// stored letters of its one word and the chunk's length.
#[derive(Clone, Copy)]
struct ScannedSubject<'a> {
    subject: BlastpSubjectRef<'a>,
    subject_sequence: &'a [u8],
    subject_length: i32,
    values: &'a BlastpPrelimValues<'a>,
    word_cutoff_scores: &'a [i32],
}

impl BlastpPrelimKernel<'_> {
    /// The kernel for one subject (`run_call`'s `process_subject` and `process_own_subject` call it
    /// with the subject's preliminary statistics `values`): `subject_result` reset for the subject,
    /// and for a subject with letters the scan and the ungapped extensions into the scratch's list of
    /// initial hits, the diagonal offset moved past the subject, the interval tree reset for it, and
    /// the gapped half (`gapped_subject`), which gives the list back to the scratch and writes the
    /// subject's HSP lists into `subject_result` (the body is [`Self::search_subject_in`]).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-498
    #[inline(always)]
    pub(super) fn search_subject(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
    ) {
        // `-lcase_masking`: this subject's scan ranges, when it has lower-case letters.
        let x_soft = self
            .x_soft_subjects
            .and_then(|all| all.get(s_idx))
            .filter(|ranges| !ranges.is_empty())
            .map(Vec::as_slice);
        self.search_subject_in_any(
            scratch,
            s_idx,
            subject,
            subject_result,
            values,
            x_soft,
            false,
        );
    }

    /// [`Self::search_subject_in`] for a search with a window of 0 or more, and for a negative
    /// `-window_size` (a search strategy's; the command line takes none) [`Self::search_subject_in_neg`]:
    /// the choice is made here, once for a subject, and the scan's loops of `search_subject_in` read
    /// no window's sign. The masked kernels (`db_masked.rs`) call this one.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn search_subject_in_any(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        scan_ranges: Option<&[(i32, i32)]>,
        chunk: bool,
    ) {
        if self.window < 0 {
            self.search_subject_in_neg(
                scratch,
                s_idx,
                subject,
                subject_result,
                values,
                scan_ranges,
                chunk,
            );
        } else {
            self.search_subject_in(
                scratch,
                s_idx,
                subject,
                subject_result,
                values,
                scan_ranges,
                chunk,
                None,
            );
        }
    }

    /// [`Self::search_subject`] with the scan ranges `scan_ranges` (`-lcase_masking`'s of a
    /// subject, or with `-db_soft_mask` a database sequence's; `None`: every word), and, with
    /// `chunk`, for one hard chunk of a sequence with `-db_hard_mask` (`db_masked.rs`): the chunk's
    /// HSPs are left in the scratch's HSP list (`gapped_subject`), sorted by score, for the merge
    /// of the sequence's chunks. With `tables` (a `-db_hard_mask` sequence that `db_masked.rs`
    /// keeps apart, `common::hard_chunk_walk::DiagonalWalk::apart`) the scan reads those tables,
    /// NCBI's of each query batch, at their offset (the reference scan over `scan_ranges`, or the
    /// whole subject), and their offset moves past the subject in place of the scratch's, whose
    /// table it does not touch. One body for every caller (no code placement in NCBI): the
    /// masked kernels call it as the default path does, so that the functions it calls keep one
    /// call site each and stay inlined into it. The four macros name the method's locals (their
    /// hygiene), so they are defined after them.
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn search_subject_in(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        scan_ranges: Option<&[(i32, i32)]>,
        chunk: bool,
        mut tables: Option<&mut QueryBatchTables>,
    ) {
        // `run_call`'s values as its former kernel closure read them: a reference or a `Copy` value
        // each.
        let Self {
            args,
            contexts,
            context_lookup_bounds,
            concatenated_query,
            lookup,
            x_context_table,
            x_dropoff_per_context,
            window,
            word_size,
            x_ungap_sent,
            diag_mask,
            total_query_span,
            timing_enabled,
            timing,
            ..
        } = *self;
        // A negative window goes to `search_subject_in_neg` (`search_subject_in_any`).
        debug_assert!(
            window >= 0,
            "a negative window is searched by search_subject_in_neg"
        );
        subject_result.reset(s_idx);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1354,1429-1432
        // A subject without letters (its warning was written when the subjects were set
        // up) is not searched, and the diagonal offset does not move for it.
        if subject.aa_len == 0 {
            return;
        }
        // The word cutoffs that the scans' per-hit calls read, as a local.
        let word_cutoff_scores = values.word_cutoff_scores;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:114-116,127
        let subject_sequence = &subject.aa_seq[1..];
        let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1439-1442
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1016-1017,1021-1022
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:939-941,945-946,369,372-374
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:223-224
        let mut init_hsps = std::mem::take(&mut scratch.init_hsps);
        init_hsps.clear();
        let subject_diag_offset = match &tables {
            Some(tables) => tables.offset,
            None => scratch.diag_offset,
        };

        let scan_ungapped_start = blastp_timing_start(timing_enabled);
        // The scan ranges (`-lcase_masking`, `-db_soft_mask`), and those of a search on `tables`
        // (the whole subject where the caller gives none).
        let x_soft = scan_ranges;
        let x_whole = [(0, subject.aa_len as i32)];
        let x_tables_ranges = x_soft.unwrap_or(&x_whole[..]);
        // The subject scan and the two-hit test in one loop (a search with a window).
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-505,509
        // Dispatch point: the two-hit test reads the lookup cells in place
        // (`x_seed_fused::x_scan_two_hit`, and `x_seed_fused_compressed::x_scan_two_hit_compressed`
        // for the compressed table of blastp-fast) and this closure takes the hits that need the
        // context check and the extension (`x_candidate_of!` below; `x_candidate_sent!` with the
        // sentinel walks).
        // EXPERIMENT (LOSAT_X_BPDIAG): inside it the hits of a cell are tested in batches.
        // The two-hit test of a hit that needs the context check and the extension, with the
        // extension `$extend` (and its leading arguments): `extend_two_hit_blosum62` for BLOSUM62,
        // `extend_two_hit_table` of the search's matrix for another one (item SC2b). Each
        // expansion is a closure type of its own, so BLOSUM62's loops keep their own instances.
        macro_rules! x_candidate_of {
            ($extend:path $(, $lead:expr)?) => {
                |query_offset: i32,
                                   subject_offset: i32,
                                   last_hit: i32,
                                   diff: i32,
                                   diag_entry: &mut DiagStruct| {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,566-573
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-510
                    // The scan's query offset is unsigned, so the context test compares the
                    // offset less the distance in unsigned arithmetic: a distance past the
                    // query's start (an entry of an earlier hard chunk of the sequence, series NP,
                    // M4b-3b) passes it.
                    let (ctx_idx, frame_base) = x_context_table
                        .get_with_start(query_offset)
                        .unwrap_or_else(|| {
                            let idx = blastp_context_idx_for_query_offset_with_bounds(
                                contexts.as_slice(),
                                context_lookup_bounds,
                                query_offset,
                            );
                            (idx, contexts[idx].frame_base)
                        });
                    if (query_offset.wrapping_sub(diff) as u32) < frame_base as u32 {
                        diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                        return;
                    }
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1109-1112
                    // The first hit's end is signed: negative where that hit lies before the hard
                    // chunk's start (series NP, M4b-3b).
                    let Some(ungapped_hit) = $extend(
                        $($lead,)?
                        concatenated_query.as_slice(),
                        subject_raw,
                        last_hit + word_size,
                        usize::try_from(subject_offset)
                            .expect("NCBI BLAST two-hit subject offset must be non-negative"),
                        usize::try_from(query_offset)
                            .expect("NCBI BLAST two-hit query offset must be non-negative"),
                        x_drop,
                        args.word_size,
                    ) else {
                        return;
                    };
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-591,596-606
                    if ungapped_hit.ungapped_data.score >= cutoff {
                        blastp_save_init_hsp(
                            &mut init_hsps,
                            ungapped_hit.ungapped_data.q_start,
                            ungapped_hit.ungapped_data.s_start,
                            query_offset,
                            subject_offset,
                            ungapped_hit.ungapped_data.length,
                            ungapped_hit.ungapped_data.score,
                        );
                    }
                    if ungapped_hit.right_extend {
                        diag_entry.set_flag(1);
                        diag_entry.set_last_hit(
                            ungapped_hit.s_last_off - (word_size - 1) + subject_diag_offset,
                        );
                    } else {
                        diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                    }
                }
            };
        }
        // The one-hit test of a hit (`-window_size 0`) on the diagonal table `$table`
        // (`common::diag_entry::DiagTable`: the scratch's ring, or `tables`) with the extension
        // `$extend` (and its leading arguments): `extend_one_hit_blosum62` for BLOSUM62, the
        // reference `extend_one_hit_of` of the search's matrix for another one (item SC2b).
        macro_rules! x_one_hit_of {
            ($table:ident, $extend:path $(, $lead:expr)?) => {
                |pair: OffsetPair| {
                    let query_offset =
                        i32::try_from(pair.q_off).expect("NCBI BLAST query offset must fit in Int4");
                    let subject_offset =
                        i32::try_from(pair.s_off).expect("NCBI BLAST subject offset must fit in Int4");
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:776-778
                    let diag_entry = $table.one_hit_entry(query_offset, subject_offset);
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:777-782,787-800,802
                    let diff = subject_offset - (diag_entry.last_hit() - subject_diag_offset);
                    if diff < 0 {
                        return;
                    }

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:783-786
                    let ctx_idx = blastp_context_idx_for_query_offset_with_bounds(
                        contexts.as_slice(),
                        context_lookup_bounds,
                        query_offset,
                    );
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1039-1042
                    let Some(ungapped_hit) = $extend(
                        $($lead,)?
                        concatenated_query.as_slice(),
                        subject_raw,
                        usize::try_from(query_offset)
                            .expect("NCBI BLAST one-hit query offset must be non-negative"),
                        usize::try_from(subject_offset)
                            .expect("NCBI BLAST one-hit subject offset must be non-negative"),
                        x_drop,
                        args.word_size,
                    ) else {
                        return;
                    };

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:794-800
                    if ungapped_hit.ungapped_data.score >= cutoff {
                        blastp_save_init_hsp(
                            &mut init_hsps,
                            ungapped_hit.ungapped_data.q_start,
                            ungapped_hit.ungapped_data.s_start,
                            query_offset,
                            subject_offset,
                            ungapped_hit.ungapped_data.length,
                            ungapped_hit.ungapped_data.score,
                        );
                    }
                    diag_entry
                        .set_last_hit(ungapped_hit.s_last_off - (word_size - 1) + subject_diag_offset);
                }
            };
        }
        // The one-hit scan of the search's table with the test `$one_hit` (`x_one_hit_of!`); the
        // reference scan of `x_tables_ranges` with `tables`.
        macro_rules! x_one_hit_ranges {
            ($one_hit:ident, $ranges:expr) => {
                x_soft_subject::x_for_each_hit(
                    &lookup,
                    subject_sequence,
                    subject.aa_len,
                    $ranges,
                    &mut scratch.offset_pairs,
                    false,
                    &mut |pair, _| $one_hit(pair),
                )
            };
        }
        macro_rules! x_one_hit_scan {
            ($one_hit:ident, $for_each_offset_pair:ident) => {
                if let Some(ranges) = x_soft {
                    x_soft_subject::x_for_each_hit(
                        &lookup,
                        subject_sequence,
                        subject.aa_len,
                        ranges,
                        &mut scratch.offset_pairs,
                        false,
                        &mut |pair, _| $one_hit(pair),
                    );
                } else {
                    $for_each_offset_pair(
                        &lookup,
                        subject_sequence,
                        subject.aa_len,
                        &mut scratch.offset_pairs,
                        BlastpOffsetPairScanMode::OneHit,
                        $one_hit,
                    );
                }
            };
        }
        // EXPERIMENT (LOSAT_X_BPUNGAPSENT): the closure of `x_candidate_of!` with the two-hit
        // extension of `x_ungapped_sent.rs` inlined (word size 3, walks stopped by the sentinels),
        // for the table of single letters, on a block and a subject that pass the path's checks,
        // with the block's byte table (BLOSUM62's, or another matrix's: item SC2c). The macro's
        // closure is a type of its own, so the loops of the closure above keep their own
        // instances.
        macro_rules! x_candidate_sent {
            ($sent_query:expr, $sent_subject:expr) => {
                |query_offset: i32,
                 subject_offset: i32,
                 last_hit: i32,
                 diff: i32,
                 diag_entry: &mut DiagStruct| {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,566-573
                    // (as in the closure above, the unsigned test included)
                    let (ctx_idx, frame_base) = x_context_table
                        .get_with_start(query_offset)
                        .unwrap_or_else(|| {
                            let idx = blastp_context_idx_for_query_offset_with_bounds(
                                contexts.as_slice(),
                                context_lookup_bounds,
                                query_offset,
                            );
                            (idx, contexts[idx].frame_base)
                        });
                    if (query_offset.wrapping_sub(diff) as u32) < frame_base as u32 {
                        diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                        return;
                    }
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];
                    let s_left_off = last_hit + word_size;
                    let s_right_off = usize::try_from(subject_offset)
                        .expect("NCBI BLAST two-hit subject offset must be non-negative");
                    let q_right_off = usize::try_from(query_offset)
                        .expect("NCBI BLAST two-hit query offset must be non-negative");
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                    // (`x_ungapped_sent::two_hit`)
                    // SAFETY: x_drop is a context's X-drop, which `XSentQuery::new` (or `new_of`)
                    // checked.
                    let extended = unsafe {
                        blastp::x_ungapped_sent::two_hit(
                            $sent_query,
                            $sent_subject,
                            s_left_off,
                            s_right_off,
                            q_right_off,
                            x_drop,
                        )
                    };
                    let Some(ungapped_hit) = extended else {
                        return;
                    };
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-591,596-606
                    // (as in the closure above)
                    if ungapped_hit.ungapped_data.score >= cutoff {
                        blastp_save_init_hsp(
                            &mut init_hsps,
                            ungapped_hit.ungapped_data.q_start,
                            ungapped_hit.ungapped_data.s_start,
                            query_offset,
                            subject_offset,
                            ungapped_hit.ungapped_data.length,
                            ungapped_hit.ungapped_data.score,
                        );
                    }
                    if ungapped_hit.right_extend {
                        diag_entry.set_flag(1);
                        diag_entry.set_last_hit(
                            ungapped_hit.s_last_off - (word_size - 1) + subject_diag_offset,
                        );
                    } else {
                        diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                    }
                }
            };
        }
        // Another matrix than BLOSUM62 (items SC2b, SC2c): the same two-hit tests on the same
        // scans, with the extension of the matrix's byte-table walks (the sentinel walks of
        // `x_candidate_sent!` where the block and the subject pass their checks, else
        // `extend_two_hit_table`; the reference where they do not cover a call), in a function of
        // their own (`x_outlined`), so that this kernel only gains its call. Its one-hit search is
        // the branch below.
        if args.window_size != 0 && args.scoring.matrix != ScoringMatrix::Blosum62 {
            x_outlined(|| {
                let x_of = (
                    args.scoring.matrix,
                    blastp::x_ungapped_fast::matrix_table(args.scoring.matrix),
                );
                let x_candidate = x_candidate_of!(extend_two_hit_table, x_of);
                if let Some(tables) = tables.as_deref_mut() {
                    x_soft_subject::x_scan_two_hit_ranges(
                        lookup,
                        subject_sequence,
                        subject.aa_len,
                        x_tables_ranges,
                        &mut scratch.offset_pairs,
                        tables,
                        subject_diag_offset,
                        window,
                        word_size,
                        x_candidate,
                    );
                } else if let Some(ranges) = x_soft {
                    x_soft_subject::x_scan_two_hit_ranges(
                        lookup,
                        subject_sequence,
                        subject.aa_len,
                        ranges,
                        &mut scratch.offset_pairs,
                        &mut DiagRing {
                            entries: &mut scratch.diag_array,
                            mask: diag_mask,
                        },
                        subject_diag_offset,
                        window,
                        word_size,
                        x_candidate,
                    );
                } else {
                    // EXPERIMENT (LOSAT_X_BPUNGAPSENT): the standard table's scans (either width)
                    // with the sentinel walks of the matrix's byte table, on a block and a subject
                    // that pass the path's checks (item SC2c).
                    let x_sent = match (&x_ungap_sent, &lookup) {
                        (
                            Some(sent_query),
                            BlastpRuntimeLookup::Aa(_) | BlastpRuntimeLookup::Aa16(_),
                        ) => {
                            let sent_subject =
                                blastp::x_ungapped_sent::XSentSubject::new(subject.aa_seq);
                            crate::utils::xreach::x_reach!(
                                "LOSAT_X_BPUNGAPSENT",
                                sent_subject.is_some()
                            );
                            sent_subject.map(|sent_subject| (sent_query, sent_subject))
                        }
                        (Some(_), _) => {
                            crate::utils::xreach::x_reach!("LOSAT_X_BPUNGAPSENT", false);
                            None
                        }
                        (None, _) => None,
                    };
                    // The tables of another matrix: the standard table (either width) or the
                    // compressed one (`blastp_matrix_compressed_lookup`: the flat table, the
                    // reference one where the flat build declines or LOSAT_X_BPCLUT is off).
                    match &lookup {
                        BlastpRuntimeLookup::Aa(aa_lookup) if x_sent.is_some() => {
                            let (sent_query, sent_subject) = x_sent.expect("checked");
                            blastp::x_seed_fused::x_scan_two_hit(
                                aa_lookup,
                                subject_sequence,
                                subject.aa_len,
                                &mut scratch.diag_array,
                                diag_mask,
                                subject_diag_offset,
                                window,
                                word_size,
                                x_candidate_sent!(sent_query, sent_subject),
                            )
                        }
                        BlastpRuntimeLookup::Aa16(aa_lookup) if x_sent.is_some() => {
                            let (sent_query, sent_subject) = x_sent.expect("checked");
                            blastp::x_seed_fused::x_scan_two_hit(
                                aa_lookup,
                                subject_sequence,
                                subject.aa_len,
                                &mut scratch.diag_array,
                                diag_mask,
                                subject_diag_offset,
                                window,
                                word_size,
                                x_candidate_sent!(sent_query, sent_subject),
                            )
                        }
                        BlastpRuntimeLookup::Aa(aa_lookup) => blastp::x_seed_fused::x_scan_two_hit(
                            aa_lookup,
                            subject_sequence,
                            subject.aa_len,
                            &mut scratch.diag_array,
                            diag_mask,
                            subject_diag_offset,
                            window,
                            word_size,
                            x_candidate,
                        ),
                        BlastpRuntimeLookup::Aa16(aa_lookup) => {
                            blastp::x_seed_fused::x_scan_two_hit(
                                aa_lookup,
                                subject_sequence,
                                subject.aa_len,
                                &mut scratch.diag_array,
                                diag_mask,
                                subject_diag_offset,
                                window,
                                word_size,
                                x_candidate,
                            )
                        }
                        BlastpRuntimeLookup::Compressed(compressed_lookup) => {
                            blastp::x_seed_fused_compressed::x_scan_two_hit_compressed(
                                compressed_lookup,
                                subject_sequence,
                                subject.aa_len,
                                &mut scratch.diag_array,
                                diag_mask,
                                subject_diag_offset,
                                window,
                                word_size,
                                x_candidate,
                            )
                        }
                        // A copy of its own of the flat scan (`X_BLASTP_MATRIX`), so that
                        // BLOSUM62's keeps its one caller.
                        BlastpRuntimeLookup::CompressedFlat(flat_lookup) => {
                            blastp::x_seed_fused_clut::x_scan_two_hit_compressed_flat::<
                                X_BLASTP_MATRIX,
                                _,
                            >(
                                flat_lookup,
                                subject_sequence,
                                subject.aa_len,
                                &mut scratch.diag_array,
                                diag_mask,
                                subject_diag_offset,
                                window,
                                word_size,
                                x_candidate,
                            )
                        }
                    }
                }
            });
        } else if args.window_size != 0 {
            let x_candidate = x_candidate_of!(extend_two_hit_blosum62);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:555-568
            let x_sent = match (&x_ungap_sent, &lookup) {
                (Some(sent_query), BlastpRuntimeLookup::Aa(aa_lookup)) => {
                    let sent_subject = blastp::x_ungapped_sent::XSentSubject::new(subject.aa_seq);
                    crate::utils::xreach::x_reach!("LOSAT_X_BPUNGAPSENT", sent_subject.is_some());
                    sent_subject.map(|sent_subject| (sent_query, sent_subject, aa_lookup))
                }
                (Some(_), BlastpRuntimeLookup::Aa16(_)) => None,
                (Some(_), _) => {
                    crate::utils::xreach::x_reach!("LOSAT_X_BPUNGAPSENT", false);
                    None
                }
                (None, _) => None,
            };
            // EXPERIMENT (LOSAT_X_BPLK16): the same choice for the table of 16-bit hits.
            let x_sent16 = match (&x_ungap_sent, &lookup) {
                (Some(sent_query), BlastpRuntimeLookup::Aa16(aa_lookup)) => {
                    let sent_subject = blastp::x_ungapped_sent::XSentSubject::new(subject.aa_seq);
                    crate::utils::xreach::x_reach!("LOSAT_X_BPUNGAPSENT", sent_subject.is_some());
                    sent_subject.map(|sent_subject| (sent_query, sent_subject, aa_lookup))
                }
                _ => None,
            };
            if let Some(tables) = tables.as_deref_mut() {
                x_soft_subject::x_scan_two_hit_ranges(
                    lookup,
                    subject_sequence,
                    subject.aa_len,
                    x_tables_ranges,
                    &mut scratch.offset_pairs,
                    tables,
                    subject_diag_offset,
                    window,
                    word_size,
                    x_candidate,
                );
            } else if let Some(ranges) = x_soft {
                x_soft_subject::x_scan_two_hit_ranges(
                    lookup,
                    subject_sequence,
                    subject.aa_len,
                    ranges,
                    &mut scratch.offset_pairs,
                    &mut DiagRing {
                        entries: &mut scratch.diag_array,
                        mask: diag_mask,
                    },
                    subject_diag_offset,
                    window,
                    word_size,
                    x_candidate,
                );
            } else if let Some((sent_query, sent_subject, aa_lookup)) = x_sent16 {
                blastp::x_seed_fused::x_scan_two_hit(
                    aa_lookup,
                    subject_sequence,
                    subject.aa_len,
                    &mut scratch.diag_array,
                    diag_mask,
                    subject_diag_offset,
                    window,
                    word_size,
                    x_candidate_sent!(sent_query, sent_subject),
                );
            } else if let Some((sent_query, sent_subject, aa_lookup)) = x_sent {
                blastp::x_seed_fused::x_scan_two_hit(
                    aa_lookup,
                    subject_sequence,
                    subject.aa_len,
                    &mut scratch.diag_array,
                    diag_mask,
                    subject_diag_offset,
                    window,
                    word_size,
                    x_candidate_sent!(sent_query, sent_subject),
                );
            } else {
                match &lookup {
                    BlastpRuntimeLookup::Aa(aa_lookup) => blastp::x_seed_fused::x_scan_two_hit(
                        aa_lookup,
                        subject_sequence,
                        subject.aa_len,
                        &mut scratch.diag_array,
                        diag_mask,
                        subject_diag_offset,
                        window,
                        word_size,
                        x_candidate,
                    ),
                    // EXPERIMENT (LOSAT_X_BPLK16): the same loop on the table of 16-bit hits.
                    BlastpRuntimeLookup::Aa16(aa_lookup) => blastp::x_seed_fused::x_scan_two_hit(
                        aa_lookup,
                        subject_sequence,
                        subject.aa_len,
                        &mut scratch.diag_array,
                        diag_mask,
                        subject_diag_offset,
                        window,
                        word_size,
                        x_candidate,
                    ),
                    BlastpRuntimeLookup::Compressed(compressed_lookup) => {
                        blastp::x_seed_fused_compressed::x_scan_two_hit_compressed(
                            compressed_lookup,
                            subject_sequence,
                            subject.aa_len,
                            &mut scratch.diag_array,
                            diag_mask,
                            subject_diag_offset,
                            window,
                            word_size,
                            x_candidate,
                        )
                    }
                    // EXPERIMENT (LOSAT_X_BPCLUT): the same loop on the flat table, each cell read as one
                    // slice (`x_seed_fused_clut.rs`).
                    BlastpRuntimeLookup::CompressedFlat(flat_lookup) => {
                        blastp::x_seed_fused_clut::x_scan_two_hit_compressed_flat::<X_BLASTP, _>(
                            flat_lookup,
                            subject_sequence,
                            subject.aa_len,
                            &mut scratch.diag_array,
                            diag_mask,
                            subject_diag_offset,
                            window,
                            word_size,
                            x_candidate,
                        )
                    }
                }
            }
        } else if args.scoring.matrix != ScoringMatrix::Blosum62 {
            x_outlined(|| match tables.as_deref_mut() {
                Some(x_tables) => {
                    let mut x_one_hit =
                        x_one_hit_of!(x_tables, extend_one_hit_of, args.scoring.matrix);
                    x_one_hit_ranges!(x_one_hit, x_tables_ranges);
                }
                None => {
                    let mut x_ring = DiagRing {
                        entries: &mut scratch.diag_array,
                        mask: diag_mask,
                    };
                    let mut x_one_hit =
                        x_one_hit_of!(x_ring, extend_one_hit_of, args.scoring.matrix);
                    x_one_hit_scan!(x_one_hit, blastp_for_each_offset_pair_of);
                }
            });
        } else if let Some(x_tables) = tables.as_deref_mut() {
            // A sequence kept apart on `tables`: the reference scan of its ranges.
            let mut x_one_hit = x_one_hit_of!(x_tables, extend_one_hit_blosum62);
            x_one_hit_ranges!(x_one_hit, x_tables_ranges);
        } else {
            // One-hit search (`-window_size 0`): the buffer of offset pairs, as NCBI's
            // `s_BlastAaWordFinder_OneHit`.
            let mut x_ring = DiagRing {
                entries: &mut scratch.diag_array,
                mask: diag_mask,
            };
            let mut x_one_hit = x_one_hit_of!(x_ring, extend_one_hit_blosum62);
            x_one_hit_scan!(x_one_hit, blastp_for_each_offset_pair);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:168-173
        match tables {
            Some(tables) => tables.finish_chunk(subject.aa_len as i32, window),
            None => scratch.finish_subject(subject.aa_len, window),
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
        let interval_tree = &mut scratch.preliminary_interval_tree;
        interval_tree.reset_with_bounds_cells(
            0,
            total_query_span,
            0,
            subject.aa_len as i32 + 1,
            x_bp_tree_max_cells(),
        );
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.scan_ungapped_ns,
                &timing.scan_ungapped_calls,
                scan_ungapped_start,
            );
        }
        // The gapped half (`gapped_subject`): from the order of the ungapped hits to the subject's HSP
        // lists in `subject_result`, with the subject's letters and word cutoffs made above. An
        // ungapped search (`-ungapped`) takes `ungapped_subject` instead: its HSPs are the ungapped
        // hits (decision D2 of task O4b: the one place where the mode enters the kernel).
        let scanned = ScannedSubject {
            subject,
            subject_sequence,
            subject_length: subject.aa_len as i32,
            values,
            word_cutoff_scores,
        };
        if args.ungapped {
            self.ungapped_subject(scratch, s_idx, scanned, subject_result, init_hsps, chunk);
        } else {
            let ends =
                self.gapped_subject(scratch, s_idx, scanned, subject_result, init_hsps, chunk);
            // The HSPs of a subject's letters lie inside the tree's bounds and are in order.
            if !ends {
                self.tree_without_end
                    .store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }
    }

    /// [`Self::search_subject_in`] for a negative `-window_size` (a hand-edited search strategy: the
    /// command line takes none): the one-hit scan of NCBI's word finder, which runs for a window of 0
    /// or less (`multiple_hits` is a window above 0), on NCBI's diagonal table of a window below 0,
    /// which keeps its state from one subject to the next. The table starts at offset `window` with
    /// every entry 0 (it can be shorter than the query: `compute_diag_table_shape`, the diagonals
    /// then share entries), an entry the scan sets holds the end of its extension plus the offset
    /// in the 31 bits of `last_hit`, a hit is extended only where it lies at or past the entry's
    /// `last_hit` less the offset (`neg_window_diff`), and after each subject the offset moves on by
    /// the subject's length plus the window, or goes back to the window with every entry set to
    /// minus the window where it has reached a quarter of `INT4_MAX` (`neg_window_exit`). A fresh entry reads
    /// as `-window` letters ahead, so a word of the first `-window` letters of a subject is not
    /// extended on a diagonal that no earlier subject has written. The caller runs the subjects on
    /// one thread in NCBI's order (`run_call`'s pool), which is what makes the table's state NCBI's.
    /// The offset and the differences are 64-bit (`BlastpSubjectScratch::neg_offset`): where NCBI's
    /// `Int4` arithmetic is defined they are its values, where it would overflow they are the widened
    /// ones. The scan, the extension and the gapped half are `search_subject_in`'s for a window of 0.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:713-808
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:36-66,87-106,128-157,162-172
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn search_subject_in_neg(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        scan_ranges: Option<&[(i32, i32)]>,
        chunk: bool,
    ) {
        let Self {
            args,
            contexts,
            context_lookup_bounds,
            concatenated_query,
            lookup,
            x_dropoff_per_context,
            word_size,
            diag_mask,
            total_query_span,
            timing_enabled,
            timing,
            ..
        } = *self;
        subject_result.reset(s_idx);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1354,1429-1432
        // A subject without letters is not searched, and the diagonal offset does not move for it.
        if subject.aa_len == 0 {
            return;
        }
        let word_cutoff_scores = values.word_cutoff_scores;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:114-116,127
        let subject_sequence = &subject.aa_seq[1..];
        let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
        let mut init_hsps = std::mem::take(&mut scratch.init_hsps);
        init_hsps.clear();
        let diag_offset = scratch.neg_offset;
        let matrix = args.scoring.matrix;
        let mut ring = DiagRing {
            entries: &mut scratch.diag_array,
            mask: diag_mask,
        };

        let scan_ungapped_start = blastp_timing_start(timing_enabled);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:772-802
        let mut one_hit = |pair: OffsetPair| {
            let query_offset =
                i32::try_from(pair.q_off).expect("NCBI BLAST query offset must fit in Int4");
            let subject_offset =
                i32::try_from(pair.s_off).expect("NCBI BLAST subject offset must fit in Int4");
            let diag_entry = ring.one_hit_entry(query_offset, subject_offset);
            if neg_window_diff(subject_offset, diag_entry.last_hit(), diag_offset) < 0 {
                return;
            }
            let ctx_idx = blastp_context_idx_for_query_offset_with_bounds(
                contexts.as_slice(),
                context_lookup_bounds,
                query_offset,
            );
            let x_drop = x_dropoff_per_context[ctx_idx];
            let cutoff = word_cutoff_scores[ctx_idx];
            let Some(ungapped_hit) = extend_one_hit_for(
                matrix,
                concatenated_query.as_slice(),
                subject_raw,
                usize::try_from(query_offset)
                    .expect("NCBI BLAST one-hit query offset must be non-negative"),
                usize::try_from(subject_offset)
                    .expect("NCBI BLAST one-hit subject offset must be non-negative"),
                x_drop,
                args.word_size,
            ) else {
                return;
            };
            if ungapped_hit.ungapped_data.score >= cutoff {
                blastp_save_init_hsp(
                    &mut init_hsps,
                    ungapped_hit.ungapped_data.q_start,
                    ungapped_hit.ungapped_data.s_start,
                    query_offset,
                    subject_offset,
                    ungapped_hit.ungapped_data.length,
                    ungapped_hit.ungapped_data.score,
                );
            }
            diag_entry.set_last_hit(neg_window_last_hit(
                ungapped_hit.s_last_off,
                word_size,
                diag_offset,
            ));
        };
        if let Some(ranges) = scan_ranges {
            x_soft_subject::x_for_each_hit(
                lookup,
                subject_sequence,
                subject.aa_len,
                ranges,
                &mut scratch.offset_pairs,
                false,
                &mut |pair, _| one_hit(pair),
            );
        } else if matrix == ScoringMatrix::Blosum62 {
            blastp_for_each_offset_pair(
                lookup,
                subject_sequence,
                subject.aa_len,
                &mut scratch.offset_pairs,
                BlastpOffsetPairScanMode::OneHit,
                one_hit,
            );
        } else {
            blastp_for_each_offset_pair_of(
                lookup,
                subject_sequence,
                subject.aa_len,
                &mut scratch.offset_pairs,
                BlastpOffsetPairScanMode::OneHit,
                one_hit,
            );
        }

        self.neg_window_exit_of(scratch, subject.aa_len as i64);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
        let interval_tree = &mut scratch.preliminary_interval_tree;
        interval_tree.reset_with_bounds_cells(
            0,
            total_query_span,
            0,
            subject.aa_len as i32 + 1,
            x_bp_tree_max_cells(),
        );
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.scan_ungapped_ns,
                &timing.scan_ungapped_calls,
                scan_ungapped_start,
            );
        }
        let scanned = ScannedSubject {
            subject,
            subject_sequence,
            subject_length: subject.aa_len as i32,
            values,
            word_cutoff_scores,
        };
        if args.ungapped {
            self.ungapped_subject(scratch, s_idx, scanned, subject_result, init_hsps, chunk);
        } else {
            let ends =
                self.gapped_subject(scratch, s_idx, scanned, subject_result, init_hsps, chunk);
            // The HSPs of a subject's letters lie inside the tree's bounds and are in order.
            if !ends {
                self.tree_without_end
                    .store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500, 560-562
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4745-4776
    // A hard chunk shorter than the word goes through the word finder and its ungapped hits
    // become its HSP list, as for any chunk.
    /// The HSPs of a hard chunk shorter than the word in an ungapped search (`-ungapped`,
    /// `db_masked.rs`): the seeds `init_hsps` that its one word saved
    /// (`common::hard_chunk_scan::OneWord`) become the chunk's HSPs as in `ungapped_subject`, left
    /// in the
    /// scratch's HSP list, sorted by score, for the merge of the sequence's chunks (`values` and
    /// `subject_result` as `search_subject_in`'s; the list of hits goes back to the scratch).
    #[cold]
    #[inline(never)]
    pub(super) fn one_word_hsps(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        init_hsps: Vec<InitHSP>,
    ) {
        // The subject of a chunk's HSPs is not read: they stay in the scratch for the merge.
        let scanned = ScannedSubject {
            subject: BlastpSubjectRef {
                aa_seq: &[0, 0],
                aa_len: 0,
            },
            subject_sequence: &[0],
            subject_length: 0,
            values,
            word_cutoff_scores: values.word_cutoff_scores,
        };
        self.ungapped_subject(scratch, s_idx, scanned, subject_result, init_hsps, true);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-545
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3820-3831
    // A hard chunk shorter than the word goes through the gapped stage as any chunk: its interval
    // tree ends one past the chunk's length (`Int4`, negative between overlapping masks) and its
    // extensions take that length, over the stored letters of its one word.
    /// The HSPs of a hard chunk shorter than the word in a gapped search (`db_masked.rs`): the
    /// seeds `init_hsps` that its one word saved (`common::hard_chunk_scan::OneWord`) through the
    /// gapped
    /// stage (`gapped_subject`) on the subject of the stored letters `letters` (the word's and the
    /// sentinel after them) and of the chunk's length `length`, left in the scratch's HSP list,
    /// sorted by score, for the merge of the sequence's chunks (`values` and `subject_result` as
    /// `search_subject_in`'s). `false` where NCBI's insertion of an HSP into the tree never ends
    /// (`BlastIntervalTree::try_add_hsp`; LOSAT's refusal,
    /// `common::hard_chunk_scan::short_chunk_tree_refusal`).
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    pub(super) fn one_word_gapped(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        init_hsps: Vec<InitHSP>,
        letters: &[u8],
        length: i32,
    ) -> bool {
        scratch.preliminary_interval_tree.reset_with_bounds_cells(
            0,
            self.total_query_span,
            0,
            length.saturating_add(1),
            x_bp_tree_max_cells(),
        );
        // The subject's own letters are not read: the chunk's HSPs stay in the scratch.
        let scanned = ScannedSubject {
            subject: BlastpSubjectRef {
                aa_seq: &[0, 0],
                aa_len: 0,
            },
            subject_sequence: letters,
            subject_length: length,
            values,
            word_cutoff_scores: values.word_cutoff_scores,
        };
        self.gapped_subject(scratch, s_idx, scanned, subject_result, init_hsps, true)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-172
    /// `Blast_ExtendWordExit` of a negative window after a subject or hard chunk of `length` letters
    /// (NCBI's `Int4`: negative for a chunk between overlapping masks): the scratch's diagonal
    /// offset moved on, or the table cleared (`neg_window_exit`).
    pub(super) fn neg_window_exit_of(&self, scratch: &mut BlastpSubjectScratch, length: i64) {
        let window = self.window;
        let (next_offset, cleared) = neg_window_exit(scratch.neg_offset, window, length);
        if cleared {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-106
            for diag in &mut scratch.diag_array {
                *diag = DiagStruct::default();
                diag.set_last_hit(neg_window_cleared_last_hit(window));
            }
        }
        scratch.neg_offset = next_offset;
    }

    /// The ungapped half of the kernel for one subject with letters (`-ungapped`, in place of
    /// `gapped_subject`): the subject's ungapped hits become its HSPs, in the order of the sorted
    /// list of hits, sorted by score, without the purge of common endpoints and without a
    /// containment test; the HSPs of a negative length are deleted (`Blast_HSPTest`); the list of
    /// each query then takes its E-values from the context's own Karlin block (`kbp_std`, no sum
    /// statistics, no Gumbel block), is reaped by the preliminary E-value and by
    /// `-qcov_hsp_perc`, and goes into `subject_result` (`ungapped_lists`). Arguments and the
    /// scratch's arrays as `gapped_subject`'s: the list of hits goes back to the scratch before the
    /// method returns; with `chunk` (one hard chunk of a sequence with `-db_hard_mask`,
    /// `db_masked.rs`) the chunk's HSPs, sorted by score, are left in the scratch's HSP list for the
    /// merge of the sequence's chunks, which ends with `ungapped_lists` on the joined list.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1482-1541
    /// (BLASTP has no `link_hsp_params` without psiblast's `-sum_stats`, which has no `-ungapped`;
    /// the bit scores are made by the final stage, `UngappedFinish::run`.)
    #[cold]
    #[inline(never)]
    fn ungapped_subject(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        scanned: ScannedSubject<'_>,
        subject_result: &mut BlastpSubjectResult,
        mut init_hsps: Vec<InitHSP>,
        chunk: bool,
    ) {
        let Self {
            contexts,
            context_lookup_bounds,
            concatenated_query,
            ..
        } = *self;
        let ScannedSubject {
            subject, values, ..
        } = scanned;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234-235
        ncbi_qsort_init_hsps_by_score(&mut init_hsps);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4752-4770
        let mut preliminary_hsp_list = std::mem::take(&mut scratch.preliminary_hsp_list);
        preliminary_hsp_list.clear();
        for init_hsp in &init_hsps {
            let mut init_hsp = *init_hsp;
            let (context_idx, ctx, _) = blastp_adjust_hsp_offsets_and_get_query_data(
                concatenated_query,
                contexts,
                context_lookup_bounds,
                &mut init_hsp,
            );
            let ungapped = init_hsp.ungapped_data;
            preliminary_hsp_list.push(BlastpPreliminaryHsp {
                query_start: ungapped.q_start,
                query_end: ungapped.length + ungapped.q_start,
                subject_start: ungapped.s_start,
                subject_end: ungapped.length + ungapped.s_start,
                gapped_query_start: init_hsp.q_seed_absolute,
                gapped_subject_start: init_hsp.s_seed,
                raw_score: ungapped.score,
                query_context: i32::try_from(context_idx)
                    .expect("NCBI BLAST preliminary blastp context index must fit in Int4"),
                query_frame: 0,
                subject_frame: 0,
                query_length: ctx.aa_len,
                q_idx: ctx.q_idx,
                s_idx: s_idx as u32,
            });
        }
        scratch.init_hsps = init_hsps;
        x_ih_pool::at_subject_end(&mut scratch.init_hsps);
        scratch.init_hsps.clear();
        if !preliminary_hits_is_sorted_by_score_ncbi(&preliminary_hsp_list) {
            ncbi_qsort_preliminary_hits(&mut preliminary_hsp_list, preliminary_score_compare);
        }
        // A hard chunk's HSPs stay in the list for the merge of the sequence's chunks
        // (`db_masked.rs`).
        if chunk {
            scratch.preliminary_hsp_list = preliminary_hsp_list;
            return;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2713-2733
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-999
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:586-591
        // `-subject_besthit` on the subject's list of every query, sorted by score, before the test
        // below deletes the HSPs of a negative length: the step sees them, as NCBI's does (it runs
        // in the engine; the reevaluation that deletes them is the caller's).
        self.subject_besthit(&mut preliminary_hsp_list);
        self.ungapped_lists(
            scratch,
            s_idx,
            subject.aa_len,
            values,
            &mut preliminary_hsp_list,
            subject_result,
        );
        scratch.preliminary_hsp_list = preliminary_hsp_list;
    }

    /// The HSP lists per query of an ungapped search's subject of `subject_len` letters from its
    /// HSPs `preliminary_hsp_list` (sorted by score, after `-subject_besthit`; drained): the HSPs of
    /// a negative length deleted, then each query's list with its E-values from the context's own
    /// Karlin block, reaped by the preliminary E-value and by `-qcov_hsp_perc`, into
    /// `subject_result`. The end of `ungapped_subject`, and of the merge of a hard-masked
    /// sequence's chunks (`db_masked.rs`), whose length is the whole sequence's.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1482-1541
    #[cold]
    #[inline(never)]
    pub(super) fn ungapped_lists(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject_len: usize,
        values: &BlastpPrelimValues<'_>,
        preliminary_hsp_list: &mut Vec<BlastpPreliminaryHsp>,
        subject_result: &mut BlastpSubjectResult,
    ) {
        let Self {
            args,
            query_ids,
            contexts,
            ..
        } = *self;
        // BLASTP has no `-perc_identity` and its minimum hit length is 0: the test deletes an HSP
        // of a negative length (`align_length` of an ungapped HSP is the length of its query range,
        // blast_hits.c:773-780) and keeps the others, whose scores did not change, in their order.
        delete_failed_hsp_tests(preliminary_hsp_list);
        // The lists of the queries (`Blast_HSPListSplitByQuery` of the HSP stream's writer), then
        // each list's E-values from the context's own block (`context_block`), the reap by the
        // preliminary E-value and by `-qcov_hsp_perc`, as NCBI's lines cited above do on the
        // subject's list (a HSP's evalue depends on its own HSP only).
        let mut preliminary_hits_by_query = std::mem::take(&mut scratch.preliminary_hits_by_query);
        split_preliminary_hits_by_query_into(
            preliminary_hsp_list,
            contexts,
            query_ids.len(),
            &mut preliminary_hits_by_query,
        );
        for (q_idx, preliminary_hits) in preliminary_hits_by_query
            .iter_mut()
            .take(query_ids.len())
            .enumerate()
        {
            if preliminary_hits.is_empty() {
                continue;
            }
            let ctx = &contexts[q_idx];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1869-1873
            // `kbp` is `sbp->kbp` (`kbp_std`) for an ungapped search and there is no Gumbel block.
            let mut compact_hsp_list = x_build_compact_preliminary_hsp_list(
                preliminary_hits.iter(),
                s_idx as u32,
                q_idx as u32,
                ctx.aa_len,
                subject_len,
                &ctx.karlin_params,
                None,
                &values.search_spaces[q_idx],
                None,
            );
            reap_hsplist_by_evalue(&mut compact_hsp_list, values.prelim_evalue);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2010-2045 (cited at
            // `core::blast_hits::reap_by_query_coverage`)
            if reap_by_query_coverage(&mut compact_hsp_list.hsps, args.query_cov_hsp_perc, |hsp| {
                (
                    i64::from(hsp.hsp.query_start),
                    i64::from(hsp.hsp.query_end),
                    hsp.hsp.query_length as i64,
                )
            }) {
                update_best_evalue(&mut compact_hsp_list);
            }
            if !compact_hsp_list.hsps.is_empty() {
                compact_hsp_list.hsps.shrink_to_fit();
                subject_result
                    .x_compact_hsp_lists
                    .push((q_idx, compact_hsp_list));
            }
            preliminary_hits.clear();
        }
        scratch.preliminary_hits_by_query = preliminary_hits_by_query;
    }

    /// The gapped half of the kernel for one subject with letters (`search_subject`'s last statement,
    /// after the scan): `scanned` is the subject with its statistics and what `search_subject` made
    /// of them, `init_hsps` are the subject's ungapped hits (the scratch's list, which goes back to
    /// the scratch here), and the subject's HSP lists go into `subject_result`. The arrays that the method
    /// takes of the scratch (the restricted alignment's flags, the HSP list and its counts per query,
    /// the lists per query) go back to it before the method returns; no statement leaves it early.
    /// Returns `false` where NCBI's insertion of an HSP into the interval tree never ends
    /// (`BlastIntervalTree::try_add_hsp`: an HSP of a hard chunk shorter than the word out of the
    /// tree's bounds, `one_word_gapped`), after which the stage stops; a whole subject's HSPs lie
    /// inside its bounds.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3666
    fn gapped_subject(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        scanned: ScannedSubject<'_>,
        subject_result: &mut BlastpSubjectResult,
        mut init_hsps: Vec<InitHSP>,
        chunk: bool,
    ) -> bool {
        // `run_call`'s values as its former kernel closure read them: a reference or a `Copy` value
        // each.
        let Self {
            args,
            query_ids,
            contexts,
            context_lookup_bounds,
            concatenated_query,
            x_drop_gapped,
            x_trace_query_base,
            timing_enabled,
            timing,
            ..
        } = *self;
        // The subject and its statistics with the letters and word cutoffs that `search_subject` made
        // of them (`scanned`), the hit-saving cutoffs (only this half reads them), and the interval
        // tree that the scan reset for the subject.
        let ScannedSubject {
            subject,
            subject_sequence,
            subject_length,
            values,
            word_cutoff_scores,
        } = scanned;
        let hit_cutoff_scores = values.hit_cutoff_scores;
        let mut ends = true;
        let interval_tree = &mut scratch.preliminary_interval_tree;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234-235
        ncbi_qsort_init_hsps_by_score(&mut init_hsps);
        let init_hsps = if args.chaining && args.matrix_name == "BLOSUM62" {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3718-3722,3726
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:569-576
            // The name is the typed one (`-matrix blosum62` does not chain; the C setter
            // `BlastScoringOptionsSetMatrix` upper-cases it, the C++ one does not).
            chain_blastp_init_hsps(
                init_hsps,
                contexts,
                word_cutoff_scores,
                hit_cutoff_scores,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                &mut scratch.chaining_keep,
                &mut scratch.chaining_nodes,
            )
        } else {
            init_hsps
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3740-3744,3748,3753,3788
        let mut restricted_align_array = std::mem::take(&mut scratch.restricted_align_array);
        if blastp_restricted_align_enabled(args.app, args.evalue) {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3736-3744,3748,3753,3788-3792
            build_restricted_align_array_into(
                &init_hsps,
                contexts,
                context_lookup_bounds,
                hit_cutoff_scores,
                query_ids.len(),
                &mut scratch.restricted_align_found,
                &mut restricted_align_array,
            );
        } else {
            restricted_align_array.clear();
            restricted_align_array.resize(query_ids.len(), false);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3971,4002-4003
        //
        // Keep the accepted preliminary HSPs in one ordered list for the
        // current subject so exact-retry rebuilds preserve NCBI's HSPList
        // order instead of HashMap iteration order.
        let mut preliminary_hsp_list = std::mem::take(&mut scratch.preliminary_hsp_list);
        preliminary_hsp_list.clear();
        let mut preliminary_hsp_counts_by_query =
            std::mem::take(&mut scratch.preliminary_hsp_counts_by_query);
        preliminary_hsp_counts_by_query.clear();
        preliminary_hsp_counts_by_query.resize(query_ids.len(), 0);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1482,1490-1497,1542
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3700-3701
        let mut redo_index: Option<usize> = None;
        let mut redo_query: Option<usize> = None;
        let mut first_init_hsp_index_by_query: Vec<Option<usize>> = Vec::new();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3878-3880
        let mut next_same_query_index: Vec<Option<usize>> = Vec::new();
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:79-80
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4211
        let preliminary_gapped_start = blastp_timing_start(timing_enabled);
        let mut hit_index = 0usize;
        while hit_index < init_hsps.len() {
            let mut init_hsp = init_hsps[hit_index];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3860,3867-3870
            let (context_idx, ctx, query_sequence) = blastp_adjust_hsp_offsets_and_get_query_data(
                concatenated_query,
                contexts,
                context_lookup_bounds,
                &mut init_hsp,
            );
            let query_index = usize::try_from(ctx.q_idx)
                .expect("NCBI BLAST requires blastp query indices to fit in usize");
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3872-3880
            if let (Some(redo_index), Some(redo_query)) = (redo_index, redo_query) {
                if hit_index < redo_index && query_index != redo_query {
                    hit_index += 1;
                    continue;
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3910-3919
            let ungapped_tree_hsp = TreeHsp {
                query_offset: init_hsp.ungapped_data.q_start,
                query_end: init_hsp.query_end_absolute(),
                subject_offset: init_hsp.ungapped_data.s_start,
                subject_end: init_hsp.subject_end(),
                score: init_hsp.ungapped_data.score,
                query_frame: 1,
                query_length: ctx.aa_len as i32,
                query_context_offset: ctx.frame_base,
                subject_frame_sign: 0,
            };
            let tree_contains = interval_tree.contains_hsp(&ungapped_tree_hsp, ctx.frame_base, 0);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3910-3919
            if let Some(target) = blastp_trace_hsp_target() {
                if blastp_trace_matches_pair(target, ctx.q_idx + x_trace_query_base, s_idx as u32)
                    && blastp_trace_init_hsp_overlaps_target(target, &init_hsp)
                {
                    eprintln!(
                        "[BLASTP_TRACE_HSP] init hit_index={} q_idx={} s_idx={} ungapped=q{}..{} s{}..{} seed=({}, {}) score={} cutoff={} mode={:?} tree_contains={}",
                        hit_index,
                        ctx.q_idx + x_trace_query_base,
                        s_idx,
                        init_hsp.ungapped_data.q_start + 1,
                        init_hsp.query_end_absolute(),
                        init_hsp.ungapped_data.s_start + 1,
                        init_hsp.subject_end(),
                        init_hsp.q_seed_absolute,
                        init_hsp.s_seed,
                        init_hsp.ungapped_data.score,
                        hit_cutoff_scores[context_idx],
                        if restricted_align_array
                            .get(query_index)
                            .copied()
                            .unwrap_or(false)
                        {
                            BlastpGappedAlignmentMode::Restricted
                        } else {
                            BlastpGappedAlignmentMode::Exact
                        },
                        tree_contains
                    );
                }
            }
            if tree_contains {
                hit_index = blastp_advance_preliminary_hit_index(
                    hit_index,
                    redo_index,
                    &next_same_query_index,
                );
                continue;
            }
            let restricted_alignment = restricted_align_array
                .get(query_index)
                .copied()
                .unwrap_or(false);
            let cutoff = hit_cutoff_scores[context_idx];
            let restricted_cutoff = (BLASTP_RESTRICTED_MULT * cutoff as f64) as i32;
            let mode = if restricted_alignment {
                BlastpGappedAlignmentMode::Restricted
            } else {
                BlastpGappedAlignmentMode::Exact
            };
            let preliminary_hsp = extend_preliminary_blastp_hit(
                &mut init_hsp,
                s_idx as u32,
                context_idx,
                ctx,
                query_sequence,
                subject_sequence,
                subject_length,
                args.scoring.matrix,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                x_drop_gapped,
                mode,
                &mut scratch.preliminary_gap_scratch,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3936-3945,3949-3953
            if let Some(target) = blastp_trace_hsp_target() {
                if blastp_trace_matches_pair(target, ctx.q_idx + x_trace_query_base, s_idx as u32)
                    && usize::try_from(preliminary_hsp.query_start + 1).ok() == Some(target.q_start)
                    && usize::try_from(preliminary_hsp.query_end).ok() == Some(target.q_end)
                    && usize::try_from(preliminary_hsp.subject_start + 1).ok()
                        == Some(target.s_start)
                    && usize::try_from(preliminary_hsp.subject_end).ok() == Some(target.s_end)
                {
                    eprintln!(
                        "[BLASTP_TRACE_HSP] preliminary_gapped hit_index={} q={}..{} s={}..{} gapped_start=({}, {}) raw_score={} cutoff={} mode={:?}",
                        hit_index,
                        preliminary_hsp.query_start + 1,
                        preliminary_hsp.query_end,
                        preliminary_hsp.subject_start + 1,
                        preliminary_hsp.subject_end,
                        preliminary_hsp.gapped_query_start,
                        preliminary_hsp.gapped_subject_start,
                        preliminary_hsp.raw_score,
                        cutoff,
                        mode
                    );
                }
            }
            if should_retry_exact_preliminary_gapped_alignment(
                restricted_alignment,
                preliminary_hsp.raw_score,
                cutoff,
                restricted_cutoff,
            ) {
                if let Some(restricted_slot) = restricted_align_array.get_mut(query_index) {
                    *restricted_slot = false;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3978,4001
                let accepted_for_query = preliminary_hsp_counts_by_query
                    .get(query_index)
                    .copied()
                    .unwrap_or(0);
                let (removed_hsps, rebuild_nodes) = if accepted_for_query == 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3978,4001
                    //
                    // With no accepted HSP for this query, NCBI's rebuild
                    // would copy every existing non-query HSP back into the
                    // tree. The current tree already contains exactly that set.
                    (0, Some(0))
                } else {
                    let removed_hsps = remove_preliminary_hits_for_query(
                        &mut preliminary_hsp_list,
                        query_index,
                        contexts,
                    );
                    debug_assert_eq!(
                        removed_hsps, accepted_for_query,
                        "NCBI redo removal should drop all accepted HSPs for the retried query"
                    );
                    if let Some(count) = preliminary_hsp_counts_by_query.get_mut(query_index) {
                        *count = 0;
                    }
                    let rebuild_nodes = rebuild_preliminary_interval_tree_from_list(
                        interval_tree,
                        &preliminary_hsp_list,
                        contexts,
                    );
                    (removed_hsps, rebuild_nodes)
                };
                let Some(rebuild_nodes) = rebuild_nodes else {
                    ends = false;
                    break;
                };
                if let Some(timing) = timing.as_ref() {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3977,4007-4008
                    BlastpTiming::record_count(&timing.preliminary_exact_retry_count, 1);
                    BlastpTiming::record_count(
                        &timing.preliminary_retry_removed_hsps,
                        removed_hsps,
                    );
                    if accepted_for_query != 0 {
                        BlastpTiming::record_count(&timing.preliminary_interval_tree_rebuilds, 1);
                        BlastpTiming::record_count(
                            &timing.preliminary_interval_tree_rebuild_nodes,
                            rebuild_nodes,
                        );
                    }
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4007-4010
                redo_index = Some(hit_index);
                redo_query = Some(query_index);
                if first_init_hsp_index_by_query.is_empty() {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3878-3880
                    //
                    // NCBI restarts at the first initial hit and skips until
                    // it finds the retried query. Rust precomputes that first
                    // per-query hit lazily and jumps to the same first
                    // processable index, preserving the subsequent hit order.
                    first_init_hsp_index_by_query.resize(query_ids.len(), None);
                    next_same_query_index.clear();
                    next_same_query_index.resize(init_hsps.len(), None);
                    let mut last_init_hsp_index_by_query = vec![None; query_ids.len()];
                    for (index, init_hsp) in init_hsps.iter().enumerate() {
                        let context_idx = blastp_get_ungapped_hsp_context_with_bounds(
                            contexts,
                            context_lookup_bounds,
                            init_hsp,
                        );
                        let init_query_index = usize::try_from(contexts[context_idx].q_idx)
                            .expect("NCBI BLAST requires blastp query indices to fit in usize");
                        if let Some(slot) = first_init_hsp_index_by_query.get_mut(init_query_index)
                        {
                            if slot.is_none() {
                                *slot = Some(index);
                            }
                        }
                        if let Some(last_index) =
                            last_init_hsp_index_by_query.get_mut(init_query_index)
                        {
                            if let Some(previous_index) = *last_index {
                                next_same_query_index[previous_index] = Some(index);
                            }
                            *last_index = Some(index);
                        }
                    }
                }
                hit_index = first_init_hsp_index_by_query
                    .get(query_index)
                    .and_then(|index| *index)
                    .unwrap_or(0);
                continue;
            }
            if preliminary_hsp.raw_score < cutoff {
                hit_index = blastp_advance_preliminary_hit_index(
                    hit_index,
                    redo_index,
                    &next_same_query_index,
                );
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:510-795 (cited at
            // `BlastIntervalTree::try_add_hsp`)
            if !interval_tree.try_add_hsp(
                preliminary_tree_hsp(&preliminary_hsp, ctx.frame_base),
                ctx.frame_base,
            ) {
                ends = false;
                break;
            }
            preliminary_hsp_list.push(preliminary_hsp);
            if let Some(count) = preliminary_hsp_counts_by_query.get_mut(query_index) {
                *count += 1;
            }
            hit_index =
                blastp_advance_preliminary_hit_index(hit_index, redo_index, &next_same_query_index);
        }
        if let Some(timing) = timing.as_ref() {
            BlastpTiming::record_call(
                &timing.preliminary_gapped_ns,
                &timing.preliminary_gapped_calls,
                preliminary_gapped_start,
            );
        }
        scratch.init_hsps = init_hsps;
        x_ih_pool::at_subject_end(&mut scratch.init_hsps);
        scratch.init_hsps.clear();

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:545,555
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478,2504
        let endpoint_duplicates_purged =
            prepare_preliminary_hits_for_kappa_in_place(&mut preliminary_hsp_list);
        if let Some(timing) = timing.as_ref() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478,2504
            BlastpTiming::record_count(
                &timing.preliminary_endpoint_duplicates_purged,
                endpoint_duplicates_purged,
            );
        }
        // The subject's HSP lists per query (`subject_lists`) after `-subject_besthit`; a hard
        // chunk's HSPs stay in the list for the merge of the sequence's chunks, which filters the
        // joined list after each chunk (`db_masked.rs`).
        if !chunk {
            self.subject_besthit(&mut preliminary_hsp_list);
            self.subject_lists(
                scratch,
                s_idx,
                subject.aa_len,
                values,
                &mut preliminary_hsp_list,
                subject_result,
            );
        }
        scratch.restricted_align_array = restricted_align_array;
        scratch.restricted_align_array.clear();
        scratch.preliminary_hsp_counts_by_query = preliminary_hsp_counts_by_query;
        scratch.preliminary_hsp_counts_by_query.clear();
        scratch.preliminary_hsp_list = preliminary_hsp_list;
        ends
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:586-591
    // `-subject_besthit` on the subject's list of every query, sorted by score, after each
    // chunk's merge: before the e-values and the reap by the preliminary e-value.
    /// `-subject_besthit`, when the search asks for it, over a subject's HSPs sorted by score
    /// (`core::blast_hits::subject_besthit`): at the end of `gapped_subject`, and after each hard
    /// chunk's merge into a masked sequence's list (`db_masked.rs`).
    #[inline(always)]
    pub(super) fn subject_besthit(&self, preliminary_hsp_list: &mut Vec<BlastpPreliminaryHsp>) {
        if self.args.subject_besthit {
            crate::core::blast_hits::subject_besthit(preliminary_hsp_list, false, |hsp| {
                crate::core::blast_hits::BestHitExtent {
                    context: i64::from(hsp.query_context),
                    query_frame: 0,
                    offset: i64::from(hsp.query_start),
                    end: i64::from(hsp.query_end),
                    query_length: 0,
                }
            });
        }
    }

    /// The HSP lists per query of a subject of `subject_len` letters from its HSPs
    /// `preliminary_hsp_list` (purged and sorted by score; drained), with the E-values of the
    /// statistics `values` and the reap by the preliminary E-value, into `subject_result`: the end
    /// of `gapped_subject`, and of the merge of a hard-masked sequence's chunks (`db_masked.rs`),
    /// whose statistics are those of the whole sequence.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:869-899
    /// After the chunks of a subject, the E-values of its HSPs with the subject's length and the
    /// reap by the preliminary E-value.
    #[inline(always)]
    pub(super) fn subject_lists(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject_len: usize,
        values: &BlastpPrelimValues<'_>,
        preliminary_hsp_list: &mut Vec<BlastpPreliminaryHsp>,
        subject_result: &mut BlastpSubjectResult,
    ) {
        let Self {
            query_ids,
            contexts,
            gapped_params,
            gapped_gumbel,
            ..
        } = *self;
        if x_hits_flat::x_bphitsflat() {
            x_hits_flat::x_lists_by_query(
                preliminary_hsp_list,
                contexts,
                query_ids.len(),
                &mut scratch.x_hits_order,
                &mut scratch.x_hits_starts,
                s_idx as u32,
                subject_len,
                gapped_params,
                gapped_gumbel.as_ref(),
                values.search_spaces,
                values.prelim_evalue,
                values.sum_stats.map(|sum_stats| &sum_stats.preliminary),
                subject_result,
            );
        } else {
            let mut preliminary_hits_by_query =
                std::mem::take(&mut scratch.preliminary_hits_by_query);
            split_preliminary_hits_by_query_into(
                preliminary_hsp_list,
                contexts,
                query_ids.len(),
                &mut preliminary_hits_by_query,
            );
            for (q_idx, preliminary_hits) in preliminary_hits_by_query
                .iter_mut()
                .take(query_ids.len())
                .enumerate()
            {
                if preliminary_hits.is_empty() {
                    continue;
                }
                let ctx = &contexts[q_idx];
                // The compact list, with the E-values of `Blast_HSPListGetEvalues` and the reap by the
                // preliminary E-value (`x_prelim_compact.rs`); a list that keeps HSPs gives back the
                // room of the reaped ones (`x_prelim_free.rs`).
                let mut compact_hsp_list = x_build_compact_preliminary_hsp_list(
                    preliminary_hits.iter(),
                    s_idx as u32,
                    q_idx as u32,
                    ctx.aa_len,
                    subject_len,
                    gapped_params,
                    gapped_gumbel.as_ref(),
                    &values.search_spaces[q_idx],
                    values.sum_stats.map(|sum_stats| &sum_stats.preliminary),
                );
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:885-888,899
                //
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:654,662-668
                reap_hsplist_by_evalue(&mut compact_hsp_list, values.prelim_evalue);
                if !compact_hsp_list.hsps.is_empty() {
                    compact_hsp_list.hsps.shrink_to_fit();
                    subject_result
                        .x_compact_hsp_lists
                        .push((q_idx, compact_hsp_list));
                }
                preliminary_hits.clear();
            }
            scratch.preliminary_hits_by_query = preliminary_hits_by_query;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2713-2733 (cited at
// `BlastpPrelimKernel::ungapped_subject`) and c++/src/algo/blast/core/blast_hits.c:993-999
/// `Blast_HSPTest` of the HSPs of an ungapped search's subject (`Blast_HSPListReevaluateUngapped`):
/// BLASTP has no `-perc_identity` and its minimum hit length is 0, so the test deletes an HSP of a
/// negative length (`align_length` of an ungapped HSP is the length of its query range,
/// blast_hits.c:773-780) and keeps the others, an HSP of length 0 among them, in their order (their
/// scores did not change).
fn delete_failed_hsp_tests(hsps: &mut Vec<BlastpPreliminaryHsp>) {
    hsps.retain(|hsp| hsp.query_end - hsp.query_start >= 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:776-782
/// The `diff` of the one-hit test for a negative window (`search_subject_in_neg`): the subject
/// offset (an `Uint4` in NCBI) less the entry's `last_hit` less the diagonal offset. Where NCBI's
/// `Int4` subtraction of the offset is defined the difference wraps in 32 bits as its unsigned
/// arithmetic does; where that subtraction would overflow it is the 64-bit difference.
fn neg_window_diff(subject_offset: i32, last_hit: i32, diag_offset: i64) -> i64 {
    let behind = i64::from(last_hit) - diag_offset;
    match i32::try_from(behind) {
        Ok(behind) => i64::from((subject_offset as u32).wrapping_sub(behind as u32) as i32),
        Err(_) => i64::from(subject_offset) - behind,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:799-800
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
/// What a one-hit extension stores in its diagonal entry for a negative window: the last subject
/// offset it reached less the word's overlap plus the diagonal offset, as the 31 bits that the
/// entry's `last_hit` holds (`DiagStruct::set_last_hit` keeps them), from the 64-bit sum.
fn neg_window_last_hit(s_last_off: i32, word_size: i32, diag_offset: i64) -> i32 {
    (i64::from(s_last_off) - i64::from(word_size - 1) + diag_offset) as i32
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-172
/// The diagonal offset after a subject or chunk of `length` letters (`Blast_ExtendWordExit`) and
/// whether the table is cleared: from an offset of a quarter of `INT4_MAX` or more it goes back to `window`
/// with every entry cleared, else it moves on by the length plus `window` (64-bit, so that a very
/// negative window does not wrap).
fn neg_window_exit(diag_offset: i64, window: i32, length: i64) -> (i64, bool) {
    if diag_offset >= i64::from(i32::MAX / 4) {
        (i64::from(window), true)
    } else {
        (diag_offset + length + i64::from(window), false)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-106
/// The `last_hit` of a cleared entry: minus the window, in the entry's 31 bits.
fn neg_window_cleared_last_hit(window: i32) -> i32 {
    (-i64::from(window)) as i32
}

#[cfg(test)]
#[path = "prelim_kernel_tests.rs"]
mod tests;

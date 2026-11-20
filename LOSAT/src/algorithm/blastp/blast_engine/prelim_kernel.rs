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
//! ```c
//! Int2 BlastAaWordFinder(BLAST_SequenceBlk * subject,
//! ...
//!     if (ewp->diag_table->multiple_hits) {
//!         status = s_BlastAaWordFinder_TwoHit(subject, query,
//! ...
//!     } else {
//!         status = s_BlastAaWordFinder_OneHit(subject, query,
//! ```
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:440,713
//! ```c
//! s_BlastAaWordFinder_TwoHit(const BLAST_SequenceBlk * subject,
//! ...
//! s_BlastAaWordFinder_OneHit(const BLAST_SequenceBlk * subject,
//! ```
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:522-526
//! ```c
//!         if (aux_struct->GetGappedScore) {
//!             status = aux_struct->GetGappedScore(program_number, query,
//!                     query_info,
//!                     subject, gap_align, score_params, ext_params, hit_params,
//!                     word_params, init_hitlist, &hsp_list, gapped_stats, NULL);
//! ```
//!
//! This module is a child of `blast_engine` for the private items of that file it reads (the
//! subject's scratch, reference and result, the preliminary values, the initial hits, the lookup
//! table, the clocks and the helpers of the scan and of the gapped stage) and its sibling modules
//! `x_hits_flat`, `x_ih_pool` and `x_soft_subject`.

use std::sync::Arc;

use super::{
    blastp_adjust_hsp_offsets_and_get_query_data, blastp_advance_preliminary_hit_index,
    blastp_context_idx_for_query_offset_with_bounds, blastp_diag_coord,
    blastp_for_each_offset_pair, blastp_for_each_offset_pair_of,
    blastp_get_ungapped_hsp_context_with_bounds, blastp_restricted_align_enabled,
    blastp_save_init_hsp, blastp_timing_start, blastp_trace_init_hsp_overlaps_target,
    build_restricted_align_array_into, chain_blastp_init_hsps, extend_preliminary_blastp_hit,
    ncbi_qsort_init_hsps_by_score, preliminary_tree_hsp,
    prepare_preliminary_hits_for_kappa_in_place, rebuild_preliminary_interval_tree_from_list,
    remove_preliminary_hits_for_query, should_retry_exact_preliminary_gapped_alignment,
    split_preliminary_hits_by_query_into, x_bp_tree_max_cells,
    x_build_compact_preliminary_hsp_list, x_hits_flat, x_ih_pool, x_outlined, x_soft_subject,
    BlastpContextLookupBounds, BlastpOffsetPairScanMode, BlastpPrelimValues, BlastpRuntimeLookup,
    BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch, BlastpTiming, InitHSP,
    BLASTP_RESTRICTED_MULT,
};
use crate::algorithm::blastn::interval_tree::{IndexMethod, TreeHsp};
use crate::algorithm::blastp;
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::extension::{
    extend_one_hit_blosum62, extend_one_hit_of, extend_two_hit_blosum62, extend_two_hit_table,
};
use crate::algorithm::blastp::gapalign::{BlastpGappedAlignmentMode, BlastpPreliminaryHsp};
use crate::algorithm::blastp::hsp::reap_hsplist_by_evalue;
use crate::algorithm::blastp::kappa::{blastp_trace_hsp_target, blastp_trace_matches_pair};
use crate::algorithm::blastp::x_seed_fused::XContextTable;
use crate::algorithm::blastp::x_ungapped_sent::XSentQuery;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{X_BLASTP, X_BLASTP_MATRIX};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::config::ScoringMatrix;
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
}

/// One subject with letters as `search_subject` hands it to its gapped half (`gapped_subject`): the
/// subject and its preliminary statistics, with the two values that `search_subject` made of them,
/// the subject's letters after the sentinel byte and the word cutoffs (one producer each).
#[derive(Clone, Copy)]
struct ScannedSubject<'a> {
    subject: BlastpSubjectRef<'a>,
    subject_sequence: &'a [u8],
    values: &'a BlastpPrelimValues<'a>,
    word_cutoff_scores: &'a [i32],
}

impl BlastpPrelimKernel<'_> {
    /// The kernel for one subject (`run_call`'s `process_subject` and `process_own_subject` call it
    /// with the subject's preliminary statistics `values`): `subject_result` reset for the subject,
    /// and for a subject with letters the scan and the ungapped extensions into the scratch's list of
    /// initial hits, the diagonal offset moved past the subject, the interval tree reset for it, and
    /// the gapped half (`gapped_subject`), which gives the list back to the scratch and writes the
    /// subject's HSP lists into `subject_result`. The four macros name the method's locals (their
    /// hygiene), so they are defined after them.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-498
    /// ```c
    ///         if (aux_struct->WordFinder) {
    ///             aux_struct->WordFinder(subject, query, query_info, lookup, matrix,
    ///                                    word_params, aux_struct->ewp,
    ///                                    aux_struct->offset_pairs,
    ///                                    kScanSubjectOffsetArraySize,
    ///                                    init_hitlist, ungapped_stats);
    /// ```
    pub(super) fn search_subject(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
    ) {
        self.search_subject_in::<false>(scratch, s_idx, subject, subject_result, values, None);
    }

    /// [`Self::search_subject`] with the scan ranges `db_soft` of a database sequence with
    /// `-db_soft_mask` (`None`: those of `-lcase_masking`'s subjects, if any), and, with `CHUNK`,
    /// for one hard chunk of a sequence with `-db_hard_mask` (`db_masked.rs`): the chunk's HSPs
    /// are left in the scratch's HSP list (`gapped_subject`), sorted by score, for the merge of
    /// the sequence's chunks. Inlined into each caller, so that `search_subject` keeps its code.
    #[inline(always)]
    pub(super) fn search_subject_in<const CHUNK: bool>(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
        db_soft: Option<&[(i32, i32)]>,
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
            x_soft_subjects,
            timing_enabled,
            timing,
            ..
        } = *self;
        subject_result.reset(s_idx);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1354,1429-1432
        // ```c
        //     Int4 min_subj_seq_length = 1;
        // ...
        //        if (seq_arg.seq->length < min_subj_seq_length) {
        //            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
        //            continue;
        //        }
        // ```
        // A subject without letters (its warning was written when the subjects were set
        // up) is not searched, and the diagonal offset does not move for it.
        if subject.aa_len == 0 {
            return;
        }
        // The word cutoffs that the scans' per-hit calls read, as a local.
        let word_cutoff_scores = values.word_cutoff_scores;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:114-116,127
        // ```c
        // (*seq_blk)->sequence_start = (Uint1 *) buffer;
        // /* The first byte is a sentinel byte. */
        // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
        // ...
        // (*seq_blk)->length = length;
        // ```
        let subject_sequence = &subject.aa_seq[1..];
        let subject_raw = &subject.aa_seq[1..subject.aa_seq.len() - 1];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1439-1442
        // ```c
        // if ((status = BLAST_OneSubjectUpdateParameters(program_number,
        //                seq_arg.seq->length, score_options, query_info,
        //                sbp, hit_params, word_params,
        //                eff_len_params)) != 0)
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1016-1017,1021-1022
        // ```c
        // BlastHitSavingParametersUpdate(program_number, sbp, query_info, subject_length, 0, /* FIXME */
        //                                hit_params);
        // ...
        // BlastInitialWordParametersUpdate(program_number, hit_params, sbp, query_info,
        //      subject_length, word_params);
        // ```
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:939-941,945-946,369,372-374
        // ```c
        // new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
        //             query_info->contexts[context].query_length,
        //             avg_subject_length);
        // ...
        // params->cutoffs[context].cutoff_score = new_cutoff;
        // params->cutoffs[context].cutoff_score_max = new_cutoff;
        // ...
        // new_cutoff = gap_trigger;
        // ...
        // new_cutoff = MIN(new_cutoff,
        //                  hit_params->cutoffs[context].cutoff_score_max);
        // curr_cutoffs->cutoff_score = new_cutoff;
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:223-224
        // ```c
        // init_hitlist->init_hsp_array = (BlastInitHSP *)
        //     malloc(MIN_INIT_HITLIST_SIZE * sizeof(BlastInitHSP));
        // ```
        let mut init_hsps = std::mem::take(&mut scratch.init_hsps);
        init_hsps.clear();
        let subject_diag_offset = scratch.diag_offset;

        let scan_ungapped_start = blastp_timing_start(timing_enabled);
        // `-lcase_masking`: this subject's scan ranges, when it has lower-case letters.
        // `-db_soft_mask`: a database sequence's scan ranges (`db_masked.rs`).
        let x_soft = match db_soft {
            Some(ranges) => Some(ranges),
            None => x_soft_subjects
                .and_then(|all| all.get(s_idx))
                .filter(|ranges| !ranges.is_empty())
                .map(Vec::as_slice),
        };
        // The subject scan and the two-hit test in one loop (a search with a window).
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-505,509
        // ```c
        // while (scan_range[1] <= scan_range[2]) {
        //     /* scan the subject sequence for hits */
        //     hits = scansub(lookup_wrap, subject,
        //                               offset_pairs, array_size, scan_range);
        //     ...
        //     for (i = 0; i < hits; ++i) {
        // ```
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
                    // ```c
                    // curr_context = BSearchContextInfo(query_offset, query_info);
                    // ...
                    // if (query_offset - diff <
                    //     query_info->contexts[curr_context].query_offset) {
                    //     diag_array[diag_coord].last_hit =
                    //         subject_offset + diag_offset;
                    //     continue;
                    // }
                    // ```
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
                    if query_offset - diff < frame_base {
                        diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                        return;
                    }
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1109-1112
                    // ```c
                    // if (use_pssm)
                    //     score += matrix[q_right_off + i][s[s_right_off + i]];
                    // else
                    //     score += matrix[q[q_right_off + i]][s[s_right_off + i]];
                    // ```
                    let Some(ungapped_hit) = $extend(
                        $($lead,)?
                        concatenated_query.as_slice(),
                        subject_raw,
                        usize::try_from(last_hit + word_size)
                            .expect("NCBI BLAST two-hit seed boundary must be non-negative"),
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
                    // ```c
                    // if (score >= cutoffs->cutoff_score)
                    //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                    //                      query_offset, subject_offset, hsp_len,
                    //                      score);
                    // ...
                    // if (right_extend) {
                    //     diag_array[diag_coord].flag = 1;
                    //     diag_array[diag_coord].last_hit =
                    //         s_last_off - (wordsize - 1) + diag_offset;
                    // }
                    // else {
                    //     diag_array[diag_coord].last_hit =
                    //         subject_offset + diag_offset;
                    // }
                    // ```
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
        // The one-hit test of a hit (`-window_size 0`) with the extension `$extend` (and its
        // leading arguments): `extend_one_hit_blosum62` for BLOSUM62, the reference
        // `extend_one_hit_of` of the search's matrix for another one (item SC2b).
        macro_rules! x_one_hit_of {
            ($extend:path $(, $lead:expr)?) => {
                |pair: OffsetPair| {
                    let query_offset =
                        i32::try_from(pair.q_off).expect("NCBI BLAST query offset must fit in Int4");
                    let subject_offset =
                        i32::try_from(pair.s_off).expect("NCBI BLAST subject offset must fit in Int4");
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:776-778
                    // ```c
                    // diag_coord = (subject_offset - query_offset) & diag_mask;
                    // diff = subject_offset -
                    //     (diag_array[diag_coord].last_hit - diag_offset);
                    // ```
                    let diag_coord = blastp_diag_coord(subject_offset - query_offset, diag_mask);
                    let diag_entry = &mut scratch.diag_array[diag_coord];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:777-782,787-800,802
                    // ```c
                    // diff = subject_offset -
                    //     (diag_array[diag_coord].last_hit - diag_offset);
                    // if (diff >= 0) {
                    // ...
                    //     score = s_BlastAaExtendOneHit(matrix, subject, query,
                    //                                   subject_offset, query_offset,
                    //                                   cutoffs->x_dropoff,
                    //                                   &hsp_q, &hsp_s, &hsp_len,
                    //                                   wordsize, use_pssm, &s_last_off);
                    //     if (score >= cutoffs->cutoff_score) {
                    //         BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                    //                          query_offset, subject_offset, hsp_len,
                    //                          score);
                    //     }
                    //     diag_array[diag_coord].last_hit =
                    //         s_last_off - (wordsize - 1) + diag_offset;
                    // ...
                    // }
                    // ```
                    let diff = subject_offset - (diag_entry.last_hit() - subject_diag_offset);
                    if diff < 0 {
                        return;
                    }

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:783-786
                    // ```c
                    // Int4 curr_context = BSearchContextInfo(query_offset,
                    //                                        query_info);
                    // BlastUngappedCutoffs *cutoffs = word_params->cutoffs +
                    //                                         curr_context;
                    // ```
                    let ctx_idx = blastp_context_idx_for_query_offset_with_bounds(
                        contexts.as_slice(),
                        context_lookup_bounds,
                        query_offset,
                    );
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];

                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1039-1042
                    // ```c
                    // if (use_pssm)
                    //     sum += matrix[q_off + i][s[s_off + i]];
                    // else
                    //     sum += matrix[q[q_off + i]][s[s_off + i]];
                    // ```
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
                    // ```c
                    // if (score >= cutoffs->cutoff_score) {
                    //     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                    //                      query_offset, subject_offset, hsp_len,
                    //                      score);
                    // }
                    // diag_array[diag_coord].last_hit =
                    //     s_last_off - (wordsize - 1) + diag_offset;
                    // ```
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
        // The one-hit scan of the search's table with the test `$one_hit` (`x_one_hit_of!`).
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
                    // (as in the closure above)
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
                    if query_offset - diff < frame_base {
                        diag_entry.set_last_hit(subject_offset + subject_diag_offset);
                        return;
                    }
                    let x_drop = x_dropoff_per_context[ctx_idx];
                    let cutoff = word_cutoff_scores[ctx_idx];
                    let s_left_off = usize::try_from(last_hit + word_size)
                        .expect("NCBI BLAST two-hit seed boundary must be non-negative");
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
                if let Some(ranges) = x_soft {
                    x_soft_subject::x_scan_two_hit_ranges(
                        lookup,
                        subject_sequence,
                        subject.aa_len,
                        ranges,
                        &mut scratch.offset_pairs,
                        &mut scratch.diag_array,
                        diag_mask,
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
            // ```c
            //     if (lookup_wrap->lut_type == eAaLookupTable) {
            //         ...
            //            lut->scansub_callback = (void *)s_BlastAaScanSubject;
            //         ...
            //     }
            //     else if (lookup_wrap->lut_type == eCompressedAaLookupTable) {
            //         ...
            //         lut->scansub_callback = (void *)s_BlastCompressedAaScanSubject;
            //     }
            // ```
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
            if let Some(ranges) = x_soft {
                x_soft_subject::x_scan_two_hit_ranges(
                    lookup,
                    subject_sequence,
                    subject.aa_len,
                    ranges,
                    &mut scratch.offset_pairs,
                    &mut scratch.diag_array,
                    diag_mask,
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
            x_outlined(|| {
                let mut x_one_hit = x_one_hit_of!(extend_one_hit_of, args.scoring.matrix);
                x_one_hit_scan!(x_one_hit, blastp_for_each_offset_pair_of);
            });
        } else {
            // One-hit search (`-window_size 0`): the buffer of offset pairs, as NCBI's
            // `s_BlastAaWordFinder_OneHit`.
            let mut x_one_hit = x_one_hit_of!(extend_one_hit_blosum62);
            x_one_hit_scan!(x_one_hit, blastp_for_each_offset_pair);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:168-173
        // ```c
        // if (ewp->diag_table->offset >= INT4_MAX / 4) {
        //     ewp->diag_table->offset = ewp->diag_table->window;
        //     s_BlastDiagClear(ewp->diag_table);
        // } else {
        //     ewp->diag_table->offset += subject_length + ewp->diag_table->window;
        // }
        // ```
        scratch.finish_subject(subject.aa_len, window);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
        // ```c
        // tree = Blast_IntervalTreeInit(0, query->length+1,
        //                               0, subject->length+1);
        // ```
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
        // lists in `subject_result`, with the subject's letters and word cutoffs made above.
        let scanned = ScannedSubject {
            subject,
            subject_sequence,
            values,
            word_cutoff_scores,
        };
        self.gapped_subject::<CHUNK>(scratch, s_idx, scanned, subject_result, init_hsps);
    }

    /// The gapped half of the kernel for one subject with letters (`search_subject`'s last statement,
    /// after the scan): `scanned` is the subject with its statistics and what `search_subject` made
    /// of them, `init_hsps` are the subject's ungapped hits (the scratch's list, which goes back to
    /// the scratch here), and the subject's HSP lists go into `subject_result`. The arrays that the method
    /// takes of the scratch (the restricted alignment's flags, the HSP list and its counts per query,
    /// the lists per query) go back to it before the method returns; no statement leaves it early.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3664-3666
    /// ```c
    /// Int2 BLAST_GetGappedScore (EBlastProgramType program_number,
    ///         BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
    ///         BLAST_SequenceBlk* subject,
    /// ```
    fn gapped_subject<const CHUNK: bool>(
        &self,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        scanned: ScannedSubject<'_>,
        subject_result: &mut BlastpSubjectResult,
        mut init_hsps: Vec<InitHSP>,
    ) {
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
            values,
            word_cutoff_scores,
        } = scanned;
        let hit_cutoff_scores = values.hit_cutoff_scores;
        let interval_tree = &mut scratch.preliminary_interval_tree;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234-235
        // ```c
        // Blast_InitHitListSortByScore(init_hitlist);
        // return status;
        // ```
        ncbi_qsort_init_hsps_by_score(&mut init_hsps);
        let init_hsps = if args.chaining && args.matrix_name == "BLOSUM62" {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3718-3722,3726
            // ```c
            // if (ext_params->options->chaining &&
            //     is_prot && strcmp(score_params->options->matrix, "BLOSUM62") == 0) {
            //     status = s_ChainingAlignment(query_info, gap_align, score_params,
            //                                  hit_params, word_params, init_hitlist);
            //     ...
            // }
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:569-576
            // ```c
            // CBlastOptionsLocal::SetMatrixName(const char* matrix)
            // {
            //     if (!matrix)
            //         return;
            //
            //     sfree(m_ScoringOpts->matrix);
            //     m_ScoringOpts->matrix = strdup(matrix);
            // }
            // ```
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
        // ```c
        // restricted_align_array =
        //     (Boolean*)calloc(Blast_GetQueryIndexFromContext(
        //                                         query_info->last_context,
        //                                         program_number) + 1,
        //                      sizeof(Boolean));
        // for (i = 0;i < init_hitlist->total;i++) {
        //     ...
        //     init_hsp = &init_hitlist->init_hsp_array[i];
        //     ...
        // }
        // ```
        let mut restricted_align_array = std::mem::take(&mut scratch.restricted_align_array);
        if blastp_restricted_align_enabled(args.app, args.evalue) {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3736-3744,3748,3753,3788-3792
            // ```c
            // if (hit_params->restricted_align && !score_params->options->is_ooframe) {
            //     int i;
            //     Boolean* found = calloc(query_info->last_context + 1, sizeof(Boolean));
            //
            //     restricted_align_array =
            //         (Boolean*)calloc(Blast_GetQueryIndexFromContext(
            //                                             query_info->last_context,
            //                                             program_number) + 1,
            //                          sizeof(Boolean));
            //     for (i = 0;i < init_hitlist->total;i++) {
            //         ...
            //         init_hsp = &init_hitlist->init_hsp_array[i];
            //         ...
            //     }
            //     if (found) {
            //         free(found);
            //     }
            // }
            // ```
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
        // ```c
        // BlastHSPList* new_hsp_list = Blast_HSPListNew(kHspNumMax);
        // ...
        // Blast_HSPListFree(hsp_list);
        // hsp_list = new_hsp_list;
        // ```
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
        // ```c
        // if (!gapped_calculation) {
        //     ...
        //     /* The following must be performed for any ungapped
        //        search with a nucleotide database. */
        //     status =
        //        Blast_HSPListReevaluateUngapped(
        //                  program_number, hsp_list, query,
        //                  seq_arg.seq, word_params, hit_params,
        //                  query_info, sbp, score_params, seq_src,
        //                  seq_arg.seq->gen_code_string);
        // ...
        // }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3700-3701
        // ```c
        // Int4 redo_index = -1; /* these are used for recomputing alignmnets if the */
        // Int4 redo_query = -1; /* approximate alignment score is inconclusive */
        // ```
        let mut redo_index: Option<usize> = None;
        let mut redo_query: Option<usize> = None;
        let mut first_init_hsp_index_by_query: Vec<Option<usize>> = Vec::new();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3878-3880
        // ```c
        // if (index < redo_index && query_index != redo_query) {
        //     continue;
        // }
        // ```
        let mut next_same_query_index: Vec<Option<usize>> = Vec::new();
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:79-80
        // ```c
        // BlastGapDP* dp_mem; /**< scratch structures for dynamic programming */
        // Int4 dp_mem_alloc;  /**< current number of structures allocated */
        // ```
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4211
        // ```c
        // s_BlastProtGappedAlignment(EBlastProgramType program,
        //    BLAST_SequenceBlk* query_blk, BLAST_SequenceBlk* subject_blk,
        //    BlastGapAlignStruct* gap_align,
        // ```
        let preliminary_gapped_start = blastp_timing_start(timing_enabled);
        let mut hit_index = 0usize;
        while hit_index < init_hsps.len() {
            let mut init_hsp = init_hsps[hit_index];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3860,3867-3870
            // ```c
            // tmp_init_hsp = init_hsp_array[index];
            // ...
            // s_AdjustHspOffsetsAndGetQueryData(query, query_info, init_hsp,
            //                                   &query_tmp, &context);
            // query_index = Blast_GetQueryIndexFromContext(context, program_number);
            // ```
            let (context_idx, ctx, query_sequence) = blastp_adjust_hsp_offsets_and_get_query_data(
                concatenated_query,
                contexts,
                context_lookup_bounds,
                &mut init_hsp,
            );
            let query_index = usize::try_from(ctx.q_idx)
                .expect("NCBI BLAST requires blastp query indices to fit in usize");
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3872-3880
            // ```c
            // /* If redo_index > -1 and redo_query > -1, the main loop is recomputing
            //    gappaed alignments for a single query, becuase the approximate
            //    alignment score was inconclusive. This recomputing was triggered when
            //    index was equal to redo_index. Until index reaches redo_index again,
            //    skip all concatenated queries with index different from redo_query */
            // if (index < redo_index && query_index != redo_query) {
            //     continue;
            // }
            // ```
            if let (Some(redo_index), Some(redo_query)) = (redo_index, redo_query) {
                if hit_index < redo_index && query_index != redo_query {
                    hit_index += 1;
                    continue;
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3910-3919
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.query.frame = query_info->contexts[context].frame;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // tmp_hsp.subject.frame = subject->frame;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
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
            // ```c
            // tmp_hsp.query.offset = q_start;
            // tmp_hsp.query.end = q_end;
            // tmp_hsp.query.frame = query_info->contexts[context].frame;
            // tmp_hsp.subject.offset = s_start;
            // tmp_hsp.subject.end = s_end;
            // tmp_hsp.subject.frame = subject->frame;
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation))
            // ```
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
                args.scoring.matrix,
                args.scoring.gap_open,
                args.scoring.gap_extend,
                x_drop_gapped,
                mode,
                &mut scratch.preliminary_gap_scratch,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3936-3945,3949-3953
            // ```c
            // max_offset =
            //    BlastGetStartForGappedAlignment(query_tmp.sequence,
            //       subject->sequence, gap_align->sbp,
            //       init_hsp->ungapped_data->q_start,
            //       init_hsp->ungapped_data->length,
            //       init_hsp->ungapped_data->s_start,
            //       init_hsp->ungapped_data->length);
            // init_hsp->offsets.qs_offsets.s_off +=
            //     max_offset - init_hsp->offsets.qs_offsets.q_off;
            // init_hsp->offsets.qs_offsets.q_off = max_offset;
            // ...
            // status = s_BlastProtGappedAlignment(program_number, &query_tmp,
            //                                     subject, gap_align,
            //                                     score_params, init_hsp,
            //                                     restricted_alignment,
            //                                     fence_hit);
            // ```
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
                // ```c
                // Blast_IntervalTreeReset(tree);
                // /* remove all HSPs computed for the current query */
                // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
                // ...
                // }
                // ```
                let accepted_for_query = preliminary_hsp_counts_by_query
                    .get(query_index)
                    .copied()
                    .unwrap_or(0);
                let (removed_hsps, rebuild_nodes) = if accepted_for_query == 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3978,4001
                    // ```c
                    // Blast_IntervalTreeReset(tree);
                    // /* remove all HSPs computed for the current query */
                    // for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
                    //     ...
                    // }
                    // ```
                    //
                    // With no accepted HSP for this query, NCBI's rebuild
                    // would copy every existing non-query HSP back into the
                    // tree. The current tree already contains exactly that set.
                    (0, 0)
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
                if let Some(timing) = timing.as_ref() {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3977,4007-4008
                    // ```c
                    // Blast_IntervalTreeReset(tree);
                    // /* remove all HSPs computed for the current query */
                    // ...
                    // redo_index = index;
                    // redo_query = query_index;
                    // ```
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
                // ```c
                // redo_index = index;
                // redo_query = query_index;
                // index = -1;
                // continue;
                // ```
                redo_index = Some(hit_index);
                redo_query = Some(query_index);
                if first_init_hsp_index_by_query.is_empty() {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3878-3880
                    // ```c
                    // if (index < redo_index && query_index != redo_query) {
                    //     continue;
                    // }
                    // ```
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
            interval_tree.add_hsp(
                preliminary_tree_hsp(&preliminary_hsp, ctx.frame_base),
                ctx.frame_base,
                IndexMethod::QueryAndSubject,
            );
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
        // ```c
        // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
        // ...
        // Blast_HSPListSortByScore(hsp_list);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478,2504
        // ```c
        // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
        // ...
        // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
        // ```
        let endpoint_duplicates_purged =
            prepare_preliminary_hits_for_kappa_in_place(&mut preliminary_hsp_list);
        if let Some(timing) = timing.as_ref() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478,2504
            // ```c
            // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
            // ...
            // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
            // ```
            BlastpTiming::record_count(
                &timing.preliminary_endpoint_duplicates_purged,
                endpoint_duplicates_purged,
            );
        }
        // The subject's HSP lists per query (`subject_lists`); a hard chunk's HSPs stay in the
        // list for the merge of the sequence's chunks (`db_masked.rs`).
        if !CHUNK {
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
                // ```c
                // Blast_HSPListGetEvalues(program_number, query_info,
                //                                  stat_length, hsp_list_out,
                //                                  score_options->gapped_calculation,
                //                                  isRPS, gap_align->sbp, 0, scale_factor);
                // ...
                // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
                // ```
                //
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:654,662-668
                // ```c
                // cutoff = hit_params->prelim_evalue;
                // ...
                // if (hsp->evalue > cutoff) {
                //     hsp_array[index] = Blast_HSPFree(hsp_array[index]);
                // } else {
                //     if (index > hsp_cnt)
                //         hsp_array[hsp_cnt] = hsp_array[index];
                //     hsp_cnt++;
                // }
                // ```
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

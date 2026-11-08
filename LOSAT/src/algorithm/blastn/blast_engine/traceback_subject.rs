//! The subject's side of the blastn traceback of one subject: the subject's letters in blastna,
//! which the traceback aligns on (`traceback_subject_blastna`), and the traceback of an ungapped
//! (`-ungapped`) search, which aligns nothing (`ungapped_traceback_hits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1526-1528,1556-1573
//! (`BLAST_ComputeTraceback_MT` fetches the subject of each batch for the traceback) and
//! :1691-1696 (the search without a traceback keeps its preliminary HSPs).
//!
//! This module is a child of `run` for the types it keeps private (`BlastnSubjectRef`,
//! `PrelimHit`, the e-value helpers). The two functions are steps of `search_query_batch`'s
//! `process_subject`, which calls them where their code was; the values they read are the fields
//! of `traceback_finish::TracebackCtx`.

use super::traceback_finish::TracebackCtx;
use super::{
    adjust_blastn_offsets, calculate_blastn_context_statistics, encode_iupac_to_blastna, ungapped,
    x_blastna_into_pool, BlastnHsp, BlastnSubjectRef, PrelimHit,
};
use crate::utils::threading::SearchPool;

/// The subject's letters in blastna for its traceback: a database sequence's own, or the letters
/// made into `s_seq_blastna_vec`, which the caller keeps while it reads them.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1528, 1567-1573
/// ```c
///             const EBlastEncoding encoding = Blast_TracebackGetEncoding(program_number);
/// ...
///                 seq_arg.oid = batch->hsplist_array[0]->oid;
///                 seq_arg.encoding = encoding;
/// ...
///                 if (BlastSeqSrcGetSequence(seqsrc, &seq_arg) < 0) {
/// ```
/// The traceback fetches the uncompressed (blastna) subject of each subject with
/// preliminary HSPs; the preliminary stage reads only the packed one. The blastna letters
/// are made here, for this subject only, instead of being kept for every subject of the
/// search. A database sequence comes with its blastna letters, read with its
/// ambiguity codes (`database::traceback_listed`).
pub(super) fn traceback_subject_blastna<'v>(
    subject: BlastnSubjectRef<'v>,
    parallel_pool: &SearchPool<'_>,
    s_seq_blastna_vec: &'v mut Vec<u8>,
) -> &'v [u8] {
    match (subject.blastna, subject.store) {
        (Some(blastna), _) => blastna,
        // EXPERIMENT (LOSAT_X_BNSTORE): made from the packed subject, on the threads of the
        // pool (one thread decoding a long subject left the others waiting at every batch).
        (None, Some(store)) => {
            *s_seq_blastna_vec = Vec::new();
            x_blastna_into_pool(store, s_seq_blastna_vec, parallel_pool);
            s_seq_blastna_vec.as_slice()
        }
        (None, None) => {
            *s_seq_blastna_vec = encode_iupac_to_blastna(
                subject
                    .letters
                    .expect("a subject reaches the traceback with its letters"),
            );
            s_seq_blastna_vec.as_slice()
        }
    }
}

/// The hits of the traceback of an ungapped (`-ungapped`) search on one subject: it aligns
/// nothing and purges nothing. The HSPs the preliminary stage kept (re-evaluated, with their
/// e-values) become the subject's hits in the kept order, with the identities counted again on
/// the subject's BLASTNA letters (the ranges the traceback decoded include the HSPs').
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1526,1691-1696
/// ```c
///             const Boolean perform_traceback = score_params->options->gapped_calculation;
/// ...
///                 } else {
///                     /* traceback skipped; compute bit scores for searches
///                        where the traceback phase is seperated from the
///                        preliminary search. */
///                     Blast_HSPListGetBitScores(hsp_list, FALSE, sbp);
///                 }
/// ```
pub(super) fn ungapped_traceback_hits(
    ctx: &TracebackCtx<'_>,
    combined_prelim_hits: &[PrelimHit],
) -> Vec<BlastnHsp> {
    let &TracebackCtx {
        s_idx,
        s_seq_blastna,
        encoded_queries_nomask,
        query_eff_searchsp,
        search_karlin_ref,
        round_down_evalue_score,
        query_lengths,
        ..
    } = ctx;

    let mut final_hits: Vec<BlastnHsp> = Vec::with_capacity(combined_prelim_hits.len());
    for prelim in combined_prelim_hits.iter() {
        let context_idx = prelim.context_idx as usize;
        let length = prelim.prelim_qe - prelim.prelim_qs;
        let matches = ungapped::count_identities_ungapped(
            encoded_queries_nomask[context_idx],
            s_seq_blastna,
            prelim.prelim_qs,
            prelim.prelim_ss,
            length,
        );
        let eff_searchsp = query_eff_searchsp[context_idx];
        let (bit_score, eval) = calculate_blastn_context_statistics(
            prelim.prelim_score,
            &search_karlin_ref[context_idx].gapped,
            eff_searchsp,
            round_down_evalue_score,
        );
        let identity = if length > 0 {
            ((matches as f64 / length as f64) * 100.0).min(100.0)
        } else {
            0.0
        };
        let query_length = query_lengths[prelim.query_idx as usize];
        let (hit_q_start, hit_q_end, hit_s_start, hit_s_end) = adjust_blastn_offsets(
            prelim.prelim_qs,
            prelim.prelim_qe,
            prelim.prelim_ss,
            prelim.prelim_se,
            query_length,
            prelim.query_frame,
        );
        final_hits.push(BlastnHsp {
            identity,
            length,
            mismatch: length - matches,
            gapopen: 0,
            q_start: hit_q_start,
            q_end: hit_q_end,
            s_start: hit_s_start,
            s_end: hit_s_end,
            e_value: eval,
            bit_score,
            num_ident: matches,
            query_frame: prelim.query_frame,
            query_length,
            q_idx: prelim.query_idx,
            s_idx: s_idx as u32,
            raw_score: prelim.prelim_score,
            internal_q_offset_0: prelim.prelim_qs,
            internal_q_end_0: prelim.prelim_qe,
            internal_s_offset_0: prelim.prelim_ss,
            internal_s_end_0: prelim.prelim_se,
            internal_query_context_offset: prelim.query_context_offset,
            gap_info: None,
            num_positives: matches,
        });
    }
    final_hits
}

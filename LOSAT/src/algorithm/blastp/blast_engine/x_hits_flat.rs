//! EXPERIMENT (LOSAT_X_BPHITSFLAT): a subject's preliminary HSPs are grouped by query through one
//! array of their positions, not copied into one vector per query that each thread keeps.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:112-132
//! ```c
//!       hsp_list_array = calloc(results->num_queries, sizeof(BlastHSPList*));
//! ...
//!       for (index = 0; index < hsp_list->hspcnt; index++) {
//!          Int4 query_index;
//!          hsp = hsp_list->hsp_array[index];
//!          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
//!
//!          if (!(tmp_hsp_list = hsp_list_array[query_index])) {
//!             hsp_list_array[query_index] = tmp_hsp_list =
//!                Blast_HSPListNew(params->hsp_num_max);
//! ...
//!          Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
//! ```
//! NCBI makes the per-query lists for each subject and hands them to the hit lists, so a thread
//! keeps nothing of them between subjects. The subject scratch kept one vector per query
//! (`BlastpSubjectScratch::preliminary_hits_by_query`, cleared, never shrunk), each as long as
//! the most HSPs that query had in any one subject of the thread: the search held, per thread,
//! the sum over the queries of each query's largest list. With 100 UniRef50 queries against
//! UniRef50 at 32 threads that was 14 MiB of capacity (v4-D12 `DECISIONS` 71).
//!
//! With the switch on, a stable counting sort puts the positions of the subject's HSPs in query
//! order ([`order_by_query`]); the HSPs of each query are read in place from the subject's list,
//! in the order of that list, which is the order the per-query vectors held them in. The lists
//! made from them, their order and every value are unchanged; a thread keeps an array of 4 bytes
//! per HSP of its largest subject and one of 4 bytes per query.
use std::sync::OnceLock;

use super::{
    blastp_query_index_from_preliminary_context, reap_hsplist_by_evalue,
    x_build_compact_preliminary_hsp_list, BlastGumbelBlk, BlastpPreliminaryHsp,
    BlastpSubjectResult, KarlinParams, QueryContext,
};

// No NCBI counterpart: reads the LOSAT_X_BPHITSFLAT switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPHITSFLAT): strict set (on unless `0`).
pub(super) fn x_bphitsflat() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPHITSFLAT"))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:116-132
/// ```c
///       for (index = 0; index < hsp_list->hspcnt; index++) {
///          Int4 query_index;
///          hsp = hsp_list->hsp_array[index];
///          query_index = Blast_GetQueryIndexFromContext(hsp->context, program);
/// ...
///          Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
/// ```
/// The positions in `hits` of the HSPs of each query, in the order of `hits` (a stable counting
/// sort): query `q` has `order[starts[q]..starts[q + 1]]`. A query index at or past
/// `num_queries` is left out, as `split_preliminary_hits_by_query_into` leaves it out.
pub(super) fn order_by_query(
    hits: &[BlastpPreliminaryHsp],
    contexts: &[QueryContext],
    num_queries: usize,
    order: &mut Vec<u32>,
    starts: &mut Vec<u32>,
) {
    starts.clear();
    starts.resize(num_queries + 1, 0);
    for hsp in hits {
        let query_index = blastp_query_index_from_preliminary_context(contexts, hsp.query_context);
        if query_index < num_queries {
            starts[query_index] += 1;
        }
    }
    // The end of each query's range; filling from the last HSP back moves each to its start.
    let mut end = 0u32;
    for start in starts.iter_mut().take(num_queries) {
        end += *start;
        *start = end;
    }
    starts[num_queries] = end;
    order.clear();
    order.resize(end as usize, 0);
    for (position, hsp) in hits.iter().enumerate().rev() {
        let query_index = blastp_query_index_from_preliminary_context(contexts, hsp.query_context);
        if query_index < num_queries {
            starts[query_index] -= 1;
            order[starts[query_index] as usize] =
                u32::try_from(position).expect("a subject's HSP count fits in u32");
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:143-151
/// ```c
///       for (index = 0; index < results->num_queries; index++) {
///          if (hsp_list_array[index]) {
///             if (!results->hitlist_array[index]) {
///                results->hitlist_array[index] =
///                   Blast_HitListNew(params->prelim_hitlist_size);
///             }
///             Blast_HitListUpdate(results->hitlist_array[index],
///                                 hsp_list_array[index]);
/// ```
/// The subject's compact lists, query by query, from the HSPs of `preliminary_hits` (emptied):
/// the loop over `preliminary_hits_by_query` in `run_resolved_in_pool` with the per-query vectors
/// replaced by [`order_by_query`].
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn x_lists_by_query(
    preliminary_hits: &mut Vec<BlastpPreliminaryHsp>,
    contexts: &[QueryContext],
    num_queries: usize,
    order: &mut Vec<u32>,
    starts: &mut Vec<u32>,
    oid: u32,
    subject_length: usize,
    karlin_params: &KarlinParams,
    gumbel_params: &BlastGumbelBlk,
    prelim_evalue: f64,
    subject_result: &mut BlastpSubjectResult,
) {
    order_by_query(preliminary_hits, contexts, num_queries, order, starts);
    for q_idx in 0..num_queries {
        let run = &order[starts[q_idx] as usize..starts[q_idx + 1] as usize];
        if run.is_empty() {
            continue;
        }
        let mut compact_hsp_list = x_build_compact_preliminary_hsp_list(
            run.iter()
                .map(|&position| &preliminary_hits[position as usize]),
            oid,
            q_idx as u32,
            contexts[q_idx].aa_len,
            subject_length,
            karlin_params,
            gumbel_params,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:654-663
        // ```c
        //    cutoff = hit_params->prelim_evalue;
        // ...
        //          if (hsp->evalue > cutoff) {
        //             hsp_array[index] = Blast_HSPFree(hsp_array[index]);
        // ```
        // The reap by the preliminary E-value, as in the per-query loop.
        reap_hsplist_by_evalue(&mut compact_hsp_list, prelim_evalue);
        if !compact_hsp_list.hsps.is_empty() {
            compact_hsp_list.hsps.shrink_to_fit();
            subject_result
                .x_compact_hsp_lists
                .push((q_idx, compact_hsp_list));
        }
    }
    preliminary_hits.clear();
}

//! EXPERIMENT (LOSAT_X_BNHITCUT, a switch of the strict set: on unless it is `0`): at two threads or
//! more, a thread that has searched a `-db` sequence drops the HSPs that the query's hit list would
//! refuse, before they wait for the owner of the window (`x_db_stream.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3266-3292
//! ```c
//!       Blast_HSPListSortByEvalue(hsp_list);
//!       hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
//!
//!       ASSERT(hit_list->heapified);
//!       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
//!       if (evalue_order < 0) {
//!          /* This hit list is less significant than any of those already saved;
//!             discard it. */
//!          Blast_HSPListFree(hsp_list);
//! ```
//!
//! Why the reports do not change: the owner publishes the best e-value of the list at the top of a
//! query's heap after each update, once the hit list is a heap. A thread drops a query's HSPs of a
//! sequence only if each of them compares worse than a published value by the e-value alone
//! (`evalue_comp`, so that a tie, which the score and the ordinal id decide, stays). The list at the
//! top of the heap is the worst the hit list holds, and an update puts in a list no worse than that
//! one, so its e-value never gets worse. The sequence's list reaches the hit list after the value
//! was published, and there it compares worse than the list at the top, which `Blast_HitListUpdate`
//! frees without changing the hit list. The value a thread reads may be older than the owner's: an
//! older value drops less. The other queries' lists do not read this query's, and the HSPs that stay
//! keep their order.
//!
//! The count of good initial extensions (the next query batch's size) is added up in the search,
//! before this step.
//!
//! Why the memory: a query's hit list keeps few sequences of many (20 with
//! `-max_target_seqs 10`), so where most sequences hit most queries (amplicons against SILVA) almost
//! all the HSPs that wait for the owner are going to be freed by it, in the arena of the thread that
//! found them, which keeps the space: about 4 MiB per thread (`v4-D12/STATE.md`, s29).
//!
//! No NCBI counterpart for the step itself: NCBI's threads write each list to the HSP stream one at a
//! time under a lock (blast_engine.c:1553-1554), so that nothing waits.

use std::cmp::Ordering;
use std::sync::atomic::{AtomicU64, Ordering as MemoryOrdering};
use std::sync::OnceLock;

use super::{evalue_comp, PrelimHit};

/// Whether LOSAT_X_BNHITCUT is on (a strict-set switch: on unless it is `0`).
pub(super) fn x_bn_hit_cut() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNHITCUT"))
}

/// For each query of a batch, the best e-value of the list at the top of its hit list's heap, as
/// the bits of an `f64`; infinity while the hit list is no heap.
pub(super) struct HitCuts {
    worst: Vec<AtomicU64>,
}

impl HitCuts {
    pub(super) fn new(num_queries: usize) -> Self {
        Self {
            worst: (0..num_queries)
                .map(|_| AtomicU64::new(f64::INFINITY.to_bits()))
                .collect(),
        }
    }

    /// The owner: `best_evalue` of the list at the top of `query`'s heap.
    pub(super) fn publish(&self, query: usize, best_evalue: f64) {
        self.worst[query].store(best_evalue.to_bits(), MemoryOrdering::Relaxed);
    }

    /// A search thread: drops the HSPs of the queries all of whose HSPs in `hits` compare worse than
    /// the published value. `seen` is the thread's bit per query (all 0 between calls).
    pub(super) fn cut(&self, hits: &mut Vec<PrelimHit>, seen: &mut Vec<u64>) {
        if seen.is_empty() {
            seen.resize(self.worst.len().div_ceil(64), 0);
        }
        for hit in hits.iter() {
            let query = hit.query_idx as usize;
            let worst = f64::from_bits(self.worst[query].load(MemoryOrdering::Relaxed));
            if evalue_comp(worst, hit.prelim_evalue) != Ordering::Less {
                seen[query / 64] |= 1u64 << (query % 64);
            }
        }
        hits.retain(|hit| {
            let query = hit.query_idx as usize;
            seen[query / 64] & (1u64 << (query % 64)) != 0
        });
        for hit in hits.iter() {
            seen[hit.query_idx as usize / 64] = 0;
        }
    }
}

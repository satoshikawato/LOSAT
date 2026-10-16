//! EXPERIMENT (LOSAT_X_BPREDOFIT): blastp's redone HSP lists hold no more memory than their HSPs
//! need.
//!
//! The composition-based redo builds a new HSP list for every subject of a query that it
//! redoes, from the alignments it found (`kappa.rs` `hsplist_from_distinct_alignments`, NCBI's
//! `s_HSPListFromDistinctAlignments` below), and the lists are kept in the queries' redone heaps
//! until the reports of the search are written. The list's vector grows by `push`, whose first
//! allocation has room for four HSPs (`BlastpHsp`, 168 bytes each), and most lists hold one. So
//! a block of queries with many matches holds about four times the room its redone HSPs need, at
//! the end of the redo, where such a block's heap peaks: human proteome block 1 (566 proteins,
//! 50,404 letters) against the human proteome kept 78,337 redone HSPs in 52.6 MB of vectors
//! (massif, one thread, `LOSAT_X_BPQBLOCK=50000`).
//!
//! With the switch on, a list is built with room for exactly the alignments it is made from, and
//! a list from which HSPs were removed afterwards (`reap_contained_hits`, `evaluate_and_purge`)
//! gives back their room (`Vec::shrink_to_fit`) before it is kept. Neither changes a value: the
//! list holds the same HSPs in the same order, and only the vector's capacity differs.
//!
//! `LOSAT_X_BPREDOFIT` is a switch of the strict set (`utils/xswitch.rs`): on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. With `LOSAT_TIMING=1` each search
//! of a query set prints `[TIMING] x_bpredofit: lists=<n> hsps=<m>` (the redone lists kept with
//! fitted room, and their HSPs).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:321-356
//! ```c
//! for (align = *alignments; NULL != align; align = align->next) {
//!     GapEditScript * editScript = align->context;
//!     align->context = NULL;
//!     status = Blast_HSPListSaveHSP(hsp_list, new_hsp);
//! }
//! BlastCompo_AlignmentsFree(alignments, s_FreeEditScript);
//! ```

use std::sync::atomic::{AtomicUsize, Ordering};

use crate::core::composition_adjustment::redo_alignment::BlastCompoAlignment;

// No NCBI counterpart: reads the LOSAT_X_BPREDOFIT switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPREDOFIT): a strict-set switch, on unless it is `0`, or `LOSAT_X_DEFAULT=0`
/// and it is not set.
pub(crate) fn x_bp_redo_fit() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPREDOFIT"))
}

/// The redone lists kept with fitted room since the last timing line, and their HSPs.
static FITTED_LISTS: AtomicUsize = AtomicUsize::new(0);
static FITTED_HSPS: AtomicUsize = AtomicUsize::new(0);

// No NCBI counterpart: counts the alignments of the list `s_HSPListFromDistinctAlignments` walks,
// to size the vector; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPREDOFIT): the number of alignments in the linked list.
pub(crate) fn x_alignment_count(mut align: Option<&BlastCompoAlignment>) -> usize {
    let mut count = 0;
    while let Some(current) = align {
        count += 1;
        align = current.next.as_deref();
    }
    count
}

// No NCBI counterpart: gives back unused room of a vector; the same elements in the same order.
/// EXPERIMENT (LOSAT_X_BPREDOFIT): the kept list's vector without spare room, and counted.
pub(crate) fn x_fit_redone_hsps<T>(hsps: &mut Vec<T>) {
    if hsps.len() < hsps.capacity() {
        hsps.shrink_to_fit();
    }
    FITTED_LISTS.fetch_add(1, Ordering::Relaxed);
    FITTED_HSPS.fetch_add(hsps.len(), Ordering::Relaxed);
}

// No NCBI counterpart: a LOSAT timing line.
/// EXPERIMENT (LOSAT_X_BPREDOFIT): with the switch on, the timing line of the lists kept since
/// the last one.
pub(crate) fn x_bp_redo_fit_print_timing() {
    if x_bp_redo_fit() {
        eprintln!(
            "[TIMING] x_bpredofit: lists={} hsps={}",
            FITTED_LISTS.swap(0, Ordering::Relaxed),
            FITTED_HSPS.swap(0, Ordering::Relaxed)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_fitted_vector_keeps_its_elements_in_order() {
        let mut hsps = Vec::with_capacity(4);
        hsps.extend([3u32, 1, 2]);
        hsps.remove(1);
        x_fit_redone_hsps(&mut hsps);
        assert_eq!(hsps, [3, 2]);
        assert_eq!(hsps.capacity(), hsps.len());
    }
}

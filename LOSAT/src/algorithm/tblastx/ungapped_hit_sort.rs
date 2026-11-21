//! Stable index-replay sorting for TBLASTX `UngappedHit` records.

use super::chaining::UngappedHit;
use super::sum_stats_linking::x_hsp_inplace;
use std::cmp::Ordering;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-486
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:990-994
// LOSAT owns complete UngappedHit values rather than NCBI pointer arrays.
// Stable-sort their original indices with the unchanged NCBI comparator, then
// replay complete original records. Comparator-equal records retain incoming
// order without adding a LOSAT-only tie-break field.
pub(crate) fn sort_ungapped_hits_by_index_replay<C>(hits: &mut [UngappedHit], compare: C)
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync,
{
    if hits.len() <= 1 {
        return;
    }

    // The index sort, then the records moved into that order inside `hits` instead of copied out
    // and cloned back (tblastx/x_hsp_inplace.rs, item 1 of its header).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
    // NCBI's qsort also moves the entries inside the one array.
    x_hsp_inplace::sort_dispatch(hits, compare);
}

#[cfg(test)]
#[path = "ungapped_hit_sort_tests.rs"]
mod tests;

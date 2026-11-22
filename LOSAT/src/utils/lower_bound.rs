//! The position that NCBI's calls of the C++ `std::lower_bound` return, found by asking about the
//! same elements in the same order. Where the elements are not in order for the value searched
//! for (an unsorted `CRangeCollection` of `GetUniqSeqCoverage`, the units of a binary
//! windowmasker file, the `indexes` of a Seq-table's sparse index), the position depends on which
//! elements the search compares and in which order, so every caller of LOSAT uses this one search.
//!
//! NCBI reference (598d8ae6): the calls of `lower_bound` that use it:
//! c++/include/util/range_coll.hpp:331-334, 355, 369 (`x_CombineWith`, `x_Subtract`),
//! c++/src/algo/winmask/seq_masker_uset_array.cpp:77 (`CSeqMaskerUsetArray::get_info`),
//! c++/src/objects/seqtable/SeqTable_sparse_index.cpp:539-540 (`GetIndexAt`).

/// The position `lower_bound` returns for the indices `start..end`, where `less(i)` answers whether
/// the element at index `i` is less than the value searched for: the first index for which `less`
/// answers false when it answers true for every index before that one and false from it on, and
/// `end` when it answers true for all. For other answers the result is where this search ends,
/// which is where NCBI's `lower_bound` ends: `less` is asked about the same indices, in the same
/// order, as the C++ `std::lower_bound` that NCBI's code calls asked a comparator that recorded
/// its calls, built with g++ 13 (`LOSAT/tests/unit/lower_bound_observed/README.md`).
///
/// The search keeps a window of `n` indices from `lo`, at first `start..end`. While the window is
/// not empty, it asks about the index `n / 2` places into the window. When the answer is true, the
/// window becomes the indices after that one; otherwise, the indices before it. It returns `lo` of
/// the empty window.
pub(crate) fn lower_bound_by(
    start: usize,
    end: usize,
    mut less: impl FnMut(usize) -> bool,
) -> usize {
    let mut lo = start;
    let mut n = end - start;
    while n > 0 {
        let half = n / 2;
        let probe = lo + half;
        if less(probe) {
            lo = probe + 1;
            n -= half + 1;
        } else {
            n = half;
        }
    }
    lo
}

#[cfg(test)]
#[path = "lower_bound_tests.rs"]
mod tests;

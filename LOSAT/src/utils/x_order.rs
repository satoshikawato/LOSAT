//! No NCBI counterpart: a stable sort that moves no element into a scratch copy, and the in-place
//! permutation it applies (also the one of TBLASTX's in-place linking, `tblastx/x_hsp_inplace.rs`).

use std::cmp::Ordering;

/// No NCBI counterpart (memory only): sorts `v` by `compare` into the order of a stable sort
/// (comparator-equal elements keep their input order) without a scratch copy of the elements.
/// NCBI sorts arrays of pointers to its HSPs (`qsort` of `hsp_array`, `hspcnt` elements of the size
/// of a `BlastHSP*`, blast_hits.c:1379-1381), so its sorts move pointers. `slice::sort_by` on the
/// elements themselves allocates a scratch for 21 or more elements: a copy of the whole slice up to
/// 8 MB (43,478 HSPs of 184 bytes), half of the slice above that (Sakai x EDL933 `-subject` blastn:
/// 38,579 HSPs, 6.8 MiB, at the end of the traceback, while the traceback's buffers are still
/// held). Here the positions are sorted instead (ties broken by the position, which is the stable
/// order and lets the in-place unstable sort run), and the sorted order is then applied with
/// `apply_permutation_in_place`: one `usize` per element. The order of a list of up to
/// `ORDER_ON_STACK` elements lives on the stack, so sorting the short lists that most subjects give
/// makes no heap allocation, as `slice::sort_by` makes none for them (16S x 100 megablast sorts
/// many lists of a few HSPs; a heap order for each was about 0.1 % of its time at 1 thread).
pub(crate) fn x_sort_stable_by_order<T, F>(v: &mut [T], mut compare: F)
where
    F: FnMut(&T, &T) -> Ordering,
{
    let n = v.len();
    if n <= 1 {
        return;
    }
    let mut on_stack = [0usize; ORDER_ON_STACK];
    let mut on_heap: Vec<usize> = Vec::new();
    // order[k]: the input position of the element that goes to position k
    let order: &mut [usize] = if n <= ORDER_ON_STACK {
        &mut on_stack[..n]
    } else {
        on_heap.resize(n, 0);
        &mut on_heap
    };
    for (k, slot) in order.iter_mut().enumerate() {
        *slot = k;
    }
    order.sort_unstable_by(|&a, &b| compare(&v[a], &v[b]).then(a.cmp(&b)));
    apply_permutation_in_place(v, order);
}

/// No NCBI counterpart: the longest list whose `x_sort_stable_by_order` order is kept on the stack
/// (256 bytes).
const ORDER_ON_STACK: usize = 32;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// NCBI's qsort moves the pointers inside the one array. `source[dst]` names the record that
// belongs at `dst`; each cycle of that permutation is walked once with swaps
// (`tblastx/x_hsp_inplace.rs`, module doc, item 1).
/// Rearranges `items` so that `items[dst]` is the element that was at `source[dst]`. `source` must
/// be a permutation of `0..items.len()`; it is left as the identity.
pub(crate) fn apply_permutation_in_place<T>(items: &mut [T], source: &mut [usize]) {
    assert_eq!(items.len(), source.len(), "permutation length");
    for start in 0..source.len() {
        let mut dst = start;
        loop {
            let src = source[dst];
            // `dst` is settled by the swap below (or needs none): mark it.
            source[dst] = dst;
            if src == start {
                // The cycle is closed: the record that started at `start` has been carried here.
                break;
            }
            items.swap(dst, src);
            dst = src;
        }
    }
}

#[cfg(test)]
#[path = "x_order_tests.rs"]
mod tests;

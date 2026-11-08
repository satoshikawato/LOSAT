//! No NCBI counterpart: a stable sort that moves no element into a scratch copy, and the in-place
//! permutation it applies (also the one of TBLASTX's in-place linking, `tblastx/x_hsp_inplace.rs`).

use std::cmp::Ordering;

/// No NCBI counterpart (memory only): sorts `v` by `compare` into the order of a stable sort
/// (comparator-equal elements keep their input order) without a scratch copy of the elements.
/// NCBI sorts arrays of pointers to its HSPs (`qsort(hsp_array, hspcnt, sizeof(BlastHSP*), ...)`,
/// blast_hits.c:1379-1381), so its sorts move pointers. `slice::sort_by` on the elements
/// themselves allocates a scratch for 21 or more elements: a copy of the whole slice up to 8 MB
/// (43,478 HSPs of 184 bytes), half of the slice above that (Sakai x EDL933 `-subject` blastn:
/// 38,579 HSPs, 6.8 MiB, at the end of the traceback, while the traceback's buffers are still
/// held). Here the positions are sorted instead (ties broken by the position, which is the stable
/// order and lets the in-place unstable sort run), and the sorted order is then applied with
/// `apply_permutation_in_place`: one `usize` per element.
pub(crate) fn x_sort_stable_by_order<T, F>(v: &mut [T], mut compare: F)
where
    F: FnMut(&T, &T) -> Ordering,
{
    let n = v.len();
    if n <= 1 {
        return;
    }
    // order[k]: the input position of the element that goes to position k
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&a, &b| compare(&v[a], &v[b]).then(a.cmp(&b)));
    apply_permutation_in_place(v, &mut order);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1380-1381
// ```c
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
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
mod tests {
    use super::*;

    // No NCBI counterpart: a fixed xorshift so the cases repeat.
    fn next(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn x_sort_stable_by_order_matches_stable_sort() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        for n in (0..70).chain([127, 128, 1000, 4099]) {
            for keys in [1u64, 2, 5, 1000, u64::MAX] {
                // (key, input position): ties on the key must keep the input order
                let input: Vec<(u64, usize)> = (0..n)
                    .map(|i| (next(&mut state) % keys.max(1), i))
                    .collect();
                let mut want = input.clone();
                want.sort_by(|a, b| a.0.cmp(&b.0));
                let mut got = input.clone();
                x_sort_stable_by_order(&mut got, |a, b| a.0.cmp(&b.0));
                assert_eq!(got, want, "n {n} keys {keys}");
                let mut got_rev = input.clone();
                x_sort_stable_by_order(&mut got_rev, |a, b| b.0.cmp(&a.0));
                let mut want_rev = input;
                want_rev.sort_by(|a, b| b.0.cmp(&a.0));
                assert_eq!(got_rev, want_rev, "n {n} keys {keys} descending");
            }
        }
    }

    #[test]
    fn x_sort_stable_by_order_keeps_sorted_and_reverses() {
        let mut sorted: Vec<u32> = (0..300).collect();
        x_sort_stable_by_order(&mut sorted, |a, b| a.cmp(b));
        assert_eq!(sorted, (0..300).collect::<Vec<_>>());
        let mut reversed: Vec<u32> = (0..300).rev().collect();
        x_sort_stable_by_order(&mut reversed, |a, b| a.cmp(b));
        assert_eq!(reversed, (0..300).collect::<Vec<_>>());
        let mut equal = vec![(7u8, 'a'), (7, 'b'), (7, 'c')];
        x_sort_stable_by_order(&mut equal, |_, _| Ordering::Equal);
        assert_eq!(equal, vec![(7, 'a'), (7, 'b'), (7, 'c')]);
    }
}
